/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package readers;
import htsjdk.samtools.*;

import java.io.IOException;
import java.util.*;

import java.io.BufferedWriter;
import java.io.FileWriter;

import channel.*;
import htsjdk.samtools.cram.build.ContainerParser;
import htsjdk.samtools.cram.build.Cram2SamRecordFactory;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CramCompressionRecord;
import htsjdk.samtools.cram.structure.CramHeader;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;

/**
 *
 * @author jiadonglin
 */
public class SignalReader {
    
    private long startTime;
    private long endTime;
    
    private BufferedWriter writer = null;
    

    public String[] chromNameMap = new String[24];
    public int[] chromLengthMap = new int[24];    
                   
    private final int isizeUpper;
    private final int isizeLower;
    private final int readLen;
    private final int fragMean;
    
    private final int minMapQ;
    private int clusteringDist; 
    
//    private int superitemMaxSpanRange;
    
    
    private Map<String, List<SAMRecord>> rpTracker;
    private ChannelParser breakChannel;
    private ChannelParser isizeLargeChannel;
    private ChannelParser isizeSmallChannel;
    private ChannelParser oemChannel;
    private ChannelParser oriChannel;
    
    private int superitemCount;
    // keep normal read per base every 1Mb length.
    private int readDepthContainerBuffer = 1000000;
    private int[] readDepthContainer = new int[readDepthContainerBuffer];
    
    private int numOfARPs;
    private int numOfRPs;
    
    private Map<String, List<int[]>> maskedRegion;
    
    public SignalReader(int fragMean, int fragStd, int cutStd, int readLen, int maxDist, int minMapQ){
        isizeUpper = fragMean + cutStd * fragStd;
        isizeLower = fragMean - cutStd * fragStd;
        
        this.readLen = readLen;
        this.fragMean = fragMean;
        this.minMapQ = minMapQ;
        this.clusteringDist = maxDist;
        if (maxDist == -1){
            this.clusteringDist = fragMean - 2 * readLen;
        }
//        superitemMaxSpanRange = 2 * fragMean;
                
    }
    /**
     * 
     * @param bamFile the alignment file
     * @param fastaIndexFile the reference index file
     * @param chrom a user specified chromosome (optional)
     * @param chromStart start of a genome region (optional)
     * @param chromEnd end of a genome region (optional)
     * @param superitemOutPath output path of the created superitems
     * @param abnormalSigOut output path of the abnormal signals (optional)
     * @throws IOException 
     */
    public void doWork(String bamFile, String fastaIndexFile, String chrom, int chromStart, int chromEnd, String superitemOutPath, String abnormalSigOut) throws IOException{
        startTime = System.currentTimeMillis();
        
        // the channel constructor need a name, it can be whatever you like, just used for naming some output file.
        breakChannel = new ChannelParser(clusteringDist, "break", abnormalSigOut);
        isizeLargeChannel = new ChannelParser(clusteringDist, "isize_large", abnormalSigOut);
        isizeSmallChannel = new ChannelParser(clusteringDist, "isize_small", abnormalSigOut);
        oemChannel = new ChannelParser(clusteringDist, "oem", abnormalSigOut);
        oriChannel = new ChannelParser(clusteringDist, "ori", abnormalSigOut);        
        
        if (!superitemOutPath.isEmpty()){
            createSuperItemWriter(superitemOutPath);            
        }
        extractSignalsFromBAM(bamFile, fastaIndexFile, chrom, chromStart, chromEnd);
                
        endTime = System.currentTimeMillis();
        writer.close();
        printSuperItemGeneratorStats();
            
    }
    
    /**
     * Processing the BAM file and extract abnormal alignments
     * @param bamFile
     * @param fastaIndexFile
     * @param chrom
     * @param regionStart
     * @param regionEnd
     * @throws IOException 
     */
    
    private void extractSignalsFromBAM(String bamFile, String fastaIndexFile, String chrom, int regionStart, int regionEnd) throws IOException{        
               
        rpTracker = new HashMap<>();
        int windowStart = 0;
        int windowEnd = 0;
        
        fileReader myFileReader = new fileReader();        
        final SamReader samReader = myFileReader.openBamFile(bamFile, ValidationStringency.SILENT, false);
        
//        CigarOps cigarOper = new CigarOps();
        myFileReader.readFaIdxFile(fastaIndexFile, chromNameMap, chromLengthMap);       
        // access user specified region
        if (chrom != null){  
            
            SAMFileHeader samFileHeader = samReader.getFileHeader();
            
            SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();
            SAMSequenceRecord refSequenceRecord = sequenceDictionary.getSequence(chrom);
            int refSequenceLength = refSequenceRecord.getSequenceLength();           
            int nWindows = refSequenceLength / readDepthContainerBuffer;
            
                  
            
            if (regionStart != 0 && regionEnd != 0){
                refSequenceLength = regionEnd - regionStart + 1;
                if (refSequenceLength <= readDepthContainerBuffer){
                    nWindows = 1;
                    readDepthContainerBuffer = refSequenceLength;
                }else{
                    nWindows = refSequenceLength/readDepthContainerBuffer;
                }
                windowStart = regionStart;
                
                
            }
            
            int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
            
            for (int i = 0; i < nWindows; i++){  
                windowEnd = windowStart + readDepthContainerBuffer;                 
                SAMRecordIterator iterator = samReader.query(chrom, windowStart, windowEnd, false);               
                
                analysisAlignment(iterator, windowStart, chrom);                
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                windowStart = windowEnd;
                
//                readDepthPreStepBuffer = copyFromReadDepthBuffer();
                readDepthPreStepBuffer = readDepthContainer;
                
                readDepthContainer = new int[readDepthContainerBuffer];
               
                writeAllSuperItems();
                    
            }
            // process remaining alignment in BAM
            SAMRecordIterator iterator = samReader.query(chrom, windowStart, refSequenceLength, false);
            analysisAlignment(iterator, windowStart, chrom);
            processRemainingSignals();
            int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
            
            System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
            
            
            writeAllSuperItems();
           
           
        } 
        // read whole genome
        else{            
            int length = chromLengthMap.length;
            SAMRecordIterator iterator;
            for (int i = 0;i < length; i ++){
                int refSequenceLength = chromLengthMap[i];
                String curChrom = chromNameMap[i];
                System.out.println("Start processing chrom: " + curChrom + ", chrom length: " + refSequenceLength);
                
                int nWindows = refSequenceLength / readDepthContainerBuffer;

                int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
                for (int k = 0; k < nWindows; k++){
                    windowEnd = windowStart + readDepthContainerBuffer;                 
                    iterator = samReader.query(curChrom, windowStart, windowEnd, false);
                    
                    analysisAlignment(iterator, windowStart,curChrom);                
                    processRemainingSignals();
                    int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                    System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                    windowStart = windowEnd;

                    readDepthPreStepBuffer = copyFromReadDepthBuffer();
                    readDepthContainer = new int[readDepthContainerBuffer];
                    
                    writeAllSuperItems();                                                                
                    
                }
                iterator = samReader.query(curChrom, windowStart, refSequenceLength, false);
                analysisAlignment(iterator, windowStart, curChrom);
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
                
                writeAllSuperItems();
                windowStart = 0;
                windowEnd = 0;
            
                rpTracker.clear();
            }
        }                                               
    }      
    
    private void extractSignalsFromCRAM(String cramFile, String refFile, String chrom, 
            int chromS, int chromE) throws IOException, IllegalAccessException{
        rpTracker = new HashMap<>();
        InputStream is = new BufferedInputStream(new FileInputStream(new File(cramFile)));
        CramHeader cramHeader = CramIO.readCramHeader(is);
        Container c = null;
        ContainerParser parser = new ContainerParser(cramHeader.getSamFileHeader());
        
        // Access specific genome region
        if (chrom != null){            
            if (cramHeader.getSamFileHeader().getSequenceIndex(chrom) < 0){
                System.err.println("Reference sequence not found for name: " + chrom);
                return;
            }
            ArrayList<CramCompressionRecord> cramRecords = new ArrayList<>(10000);
            while(true){
                parser.getRecords(c, cramRecords, ValidationStringency.SILENT); 
                Cram2SamRecordFactory c2sFactory = new Cram2SamRecordFactory(cramHeader.getSamFileHeader());

            }
        }
        
        
    }
    
    /**
     * Analysis each BAM record through different channels
     * @param iterator
     * @param windowStart
     * @param cigarOper 
     */
    private void analysisAlignment(SAMRecordIterator iterator, int windowStart, String curRefName){
//        CigarOps corasenCigar = new CigarOps();
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();                        
//            int recordChrIdx = record.getReferenceIndex();
//            String recordChrName = record.getReferenceName();
                       
           
            List<CigarElement> cigarElements = record.getCigar().getCigarElements();
            
            if (badReads(cigarElements, record.getMappingQuality()) || record.getDuplicateReadFlag()){
                continue;
            }
            // count the number of normal read per base
            int isGoodAlign = exactAlignment(cigarElements);
            if (isGoodAlign != -1){
                updateReadDepthArray(record.getAlignmentStart(), isGoodAlign, windowStart);
            }                                      
           
            SEClippedParser(record, cigarElements, curRefName);   
            RPUnmappedParser(record, curRefName);
            RPisizeParser(record, cigarElements, curRefName);
            
            if (!rpTracker.containsKey(record.getReadName())){
                List<SAMRecord> records = new ArrayList<>();
                records.add(record);
                rpTracker.put(record.getReadName(), records);
            }else{
                List<SAMRecord> records = rpTracker.get(record.getReadName());
                records.add(record);
                numOfRPs += 1;
                RPoriParser(records, cigarElements, curRefName);
                rpTracker.remove(record.getReadName());
            }                
            
            
        } 
        iterator.close();
    }
    
    /**
     * Copy the read depth buffer in current window and save it for next window usage
     * @return 
     */
    
    private int[] copyFromReadDepthBuffer(){
        int[] newBuffer = new int[readDepthContainerBuffer];
//        int startPosToCopy = readDepthContainerBuffer - readLen;
        for (int i = 0; i < readDepthContainerBuffer; i ++){
            int val = readDepthContainer[i];
//            newBuffer[i - startPosToCopy] = val;
            newBuffer[i] = val;
        }
        return newBuffer;
    }
    /**
     * Discard reads of clipped length longer than 70% of the read length and reads with low mapping quality
     * @param cigarElements
     * @return 
     */
    private boolean badReads(List<CigarElement> cigarElements, int quality){
        int clippedLength = 0;
        boolean isBad = false;
        for (CigarElement element : cigarElements){
            String operation = element.getOperator().toString();
            int optLength = element.getLength();
            if (operation.equals("S") || operation.equals("H")){
                clippedLength += optLength;
            }
        }
        if (clippedLength > 0.7 * readLen){
            isBad = true;
        }
        if (quality <= minMapQ){
            isBad = true;
        }
        return isBad;
    }
    /**
     * Get exact matched read
     * @param cigarElements
     * @return 
     */
    private int exactAlignment(List<CigarElement> cigarElements){
        if (cigarElements.size() == 1){
            String cigarOperation = cigarElements.get(0).getOperator().toString();
            int opLength = cigarElements.get(0).getLength();
            return cigarOperation.equals("M") ? opLength : -1;
        }
        else return -1;
    }
    
    /**
     * Keep track of read depth at each site`
     * @param pos
     * @param length
     * @param windowStart 
     */
    private void updateReadDepthArray(int pos, int length, int windowStart){
        
        for (int i = 0; i < length ; i ++){
            if ( (pos + i - windowStart) >= readDepthContainerBuffer){
                continue;
            }
            else if (pos + i < windowStart) {
                continue;
            }
            else{
                readDepthContainer[pos + i - windowStart] += 1;            
            }
            
        }
        
    }
  
    /**
     * Process soft clipped reads
     * @param record
     * @param cigarElements
     * @param cigarOper 
     */
    
    private void SEClippedParser(SAMRecord record, List<CigarElement> cigarElements, String curRefName) {
        // For a mapped read and read of relatively high mapQ
        if (!record.getReadUnmappedFlag()){
            CigarOps cigarOper = new CigarOps(record);
            String firstOperation = cigarElements.get(0).getOperator().toString();
                        
            
            cigarOper.calQueryPosFromCigar(cigarElements, 1, record.getReadNegativeStrandFlag(),readLen);
            int qsPos = cigarOper.getqsPos();
            
            String cigarStr = cigarOper.getCigarStr();
            int mutCoord = record.getAlignmentStart();           
                                    
            if (!cigarStr.equals("M") && !cigarStr.isEmpty()){
                if (firstOperation.equals("M")){                    
                    mutCoord += qsPos;                
                }
                if (cigarOper.isCoIDread() && firstOperation.equals("S")){
                    mutCoord += qsPos;
                }
                
                String ori = record.getReadNegativeStrandFlag() ? "-" : "+";                               
                MutSignal mutSignal = new MutSignal(record, cigarStr, mutCoord, ori);
                mutSignal.setIsizeNormal(isizeUpper, isizeLower);

                breakChannel.addSignals(mutSignal, fragMean, readLen, curRefName);
            }
        }
        
    }
    
     /**
     * One end unmapped read
     * @param record
     * @return 
     */
    private void RPUnmappedParser(SAMRecord record, String curRefName){
        
        // read unmapped
        if (record.getReadUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
           
            oemChannel.addSignals(mutSignal, fragMean, readLen, curRefName);
            numOfARPs += 1;

        }else if (record.getMateUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
                        
            oemChannel.addSignals(mutSignal, fragMean, readLen, curRefName);
            numOfARPs += 1;
        }

    }
    /**
     * Process PE of abnormal insert size
     * @param record
     * @param cigarElements 
     */
    private void RPisizeParser(SAMRecord record, List<CigarElement> cigarElements, String curRefName){
        // only process read-pair mapped on the same chrom.
       
        CigarElement leftMostCigarElement = cigarElements.get(0);
        String leftMostCigarOperator = leftMostCigarElement.getOperator().toString();
      
        int mutCoord = record.getAlignmentStart();
        if (leftMostCigarOperator.equals("M") && !record.getReadNegativeStrandFlag()){
            mutCoord += leftMostCigarElement.getLength();
        }

        int insertSize = record.getInferredInsertSize();

        String ori = record.getReadNegativeStrandFlag() ? "-" : "+";
        if (Math.abs(insertSize) >= isizeUpper){    
//                System.out.println(record.getReadName() + " isize: " +Math.abs(insertSize));

            MutSignal mutSignal = new MutSignal(record, "ARP_LARGE_INSERT", mutCoord, ori);               
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), "ARP_LARGE_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
            isizeLargeChannel.addSignals(mutSignal, fragMean, readLen, curRefName);            
            numOfARPs += 1;
        }
        else if (Math.abs(insertSize) <= isizeLower && insertSize != 0){

            MutSignal mutSignal = new MutSignal(record, "ARP_SMALL_INSERT", mutCoord, ori);       
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                record.getInferredInsertSize(), "ARP_SMALL_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);

            isizeSmallChannel.addSignals(mutSignal, fragMean, readLen, curRefName); 
            numOfARPs += 1;
            
        }
                
    }
    /**
     * Process read pairs with abnormal orientation
     * @param records
     * @param cigarElements 
     */
    private void RPoriParser(List<SAMRecord> records, List<CigarElement> cigarElements, String curRefName){
        
        SAMRecord leftMostRecord = records.get(0);
        SAMRecord rightMostRecord = records.get(records.size() - 1);
        // For read-pair, it should be proper paired. Its read and mate are all mapped.
        if (leftMostRecord.getReadPairedFlag() && !leftMostRecord.getReadUnmappedFlag() && !leftMostRecord.getMateUnmappedFlag()){
            int mutCoord = leftMostRecord.getAlignmentStart();
            if (leftMostRecord.getReadNegativeStrandFlag() == leftMostRecord.getMateNegativeStrandFlag()){
                String mutType;
                String ori = leftMostRecord.getReadNegativeStrandFlag() ? "-":"+";
                if (leftMostRecord.getReadNegativeStrandFlag()){
                    mutType = "ARP_RR";
                }else{
                    CigarElement leftMostCigarElement = cigarElements.get(0);
                    String leftMostCigarOperation = leftMostCigarElement.getOperator().toString();

                    if (leftMostCigarOperation.equals("M")){
                        mutCoord += leftMostCigarElement.getLength();
                    }
                    mutType = "ARP_FF";
                }
                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, ori);
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), ori);
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);   

                oriChannel.addSignals(readMutSignal, fragMean, readLen, curRefName);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen, leftMostRecord.getMateReferenceName());
                numOfARPs += 1;
            }
            else if (leftMostRecord.getReadNegativeStrandFlag() && !leftMostRecord.getMateNegativeStrandFlag()){
                String mutType = "ARP_RF";

                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, "-");
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);
                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), "+");
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                oriChannel.addSignals(readMutSignal, fragMean, readLen, curRefName);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen, leftMostRecord.getMateReferenceName());

                numOfARPs += 1;
            }
        }       
    }

           
    private void processRemainingSignals() {
        breakChannel.processFinalSignals(fragMean, readLen);
        isizeLargeChannel.processFinalSignals(fragMean, readLen);
        isizeSmallChannel.processFinalSignals(fragMean, readLen);
        oemChannel.processFinalSignals(fragMean, readLen);
        oriChannel.processFinalSignals(fragMean, readLen);
                             
    }
    
    private void writeAllSuperItems() throws IOException{
        if (writer != null){
            breakChannel.writeSuperItemsInChannel(writer);
            isizeLargeChannel.writeSuperItemsInChannel(writer);
            isizeSmallChannel.writeSuperItemsInChannel(writer);
            oemChannel.writeSuperItemsInChannel(writer);
            oriChannel.writeSuperItemsInChannel(writer);

        }
        
    }
    /**
     * calculate normal read aligned at a specific position and the number of superitems that generated within this window.
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer
     * @return 
     */
    private int assignReadDepthAndCountSuperItem(int windowStart, int windowSize, int[] preReadDepthBuffer){
        breakChannel.setSuperitemWeightRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeLargeChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeSmallChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oriChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oemChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);

        
        int curWindowSuperItem = 0;
        curWindowSuperItem += breakChannel.getSuperitemCount();
        curWindowSuperItem += isizeLargeChannel.getSuperitemCount();
        curWindowSuperItem += isizeSmallChannel.getSuperitemCount();
        curWindowSuperItem += oemChannel.getSuperitemCount();
        curWindowSuperItem += oriChannel.getSuperitemCount();

        
        
        return curWindowSuperItem;
    }
    
    private void createSuperItemWriter(String superitemOutPath) throws IOException{
        writer = new BufferedWriter(new FileWriter(superitemOutPath));
        writer.write("type\tchromIdx\tnread\tpos\tsaPos\tori\tweight\tratio\tsumMapQ\tplusRead\tminusRead\tsplitRead\tsplitMapQ\titxRead\tregion\tmateRegion\tqnames\tmConsensus\tcConsensus\n");

    }
    
    private boolean goodClippedRead(CigarOps cigarOps){
        int[] clippedLength = cigarOps.getClippedStatus();
        return clippedLength[0] >= 0.05 * readLen || clippedLength[1] >= 0.05 * readLen; 
    }
    
    public int getWGARPNum(){
        return numOfARPs;
    }    
        
    public void printSuperItemGeneratorStats(){
        StringBuilder sb = new StringBuilder();
        sb.append("\n==============  SuperItem Generation =============\n");
        sb.append("Time: " + (endTime - startTime) + "ms");  
        sb.append("\nTotal superitems: " + superitemCount);   
        sb.append("\nTotal number of read-pairs: " + numOfRPs);
        sb.append("\nTotal number of doiscordant read-pairs: " + numOfARPs);
        sb.append("\n======================================================");
        
        System.out.println(sb.toString());
    }
    
}
