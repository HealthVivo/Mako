/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.io.FileWriter;
import java.util.Map.Entry;

import dataStructure.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import readers.fileReader;

import utils.MemoryLogger;
import utils.superItemLink;
import strMatch.stringMatcher;


/**
 *
 * @author jiadonglin
 */
public class CFSPM {
    
    private long startTime;
    private long endTime;
    
    final private int minsuppAbsolute;

    BufferedWriter rawPatternWriter = null;
    BufferedWriter idWriter;
    BufferedWriter intMeaWriter;
    stringMatcher strMatcher;
    
    boolean showPatternIdentifiers = false;
    boolean hasMaskRegion = false;
    
    private int patternCount;
    private final double confAF; // Minimum AF for a Super-Item to estimate breakpoints
    private SequentialPatterns patterns = null;
    
    // Save all generate patterns during pattern growth
    final private List<List<pseudoSequentialPattern>> patternCandidates = new ArrayList<>();
     
    
    String[] chrIdxNameMap;
    final private int patternSpanMaxRegion;
    private SequenceDatabase database;
    private Map<String, List<ItemSeqIdentifier>> itemAppearMap;
    private Map<String, List<int[]>> maskedRegion;
      
    
    public CFSPM(int minSup, int maxRegionSpan, double af){
        minsuppAbsolute = minSup;
        patternSpanMaxRegion = maxRegionSpan;
        confAF = af;
    }
    
    
    /**
     * Main function
     * @param database  Superitem sequence database
     * @param rawPatternOut    discovered frequent abnormal patterns
     * @param mergedPatternOut  merged patterns
     * @param susRegionWriter   patterns indicate suspect regions waiting for further processing
     * @param svRegionWriter    predicted SV breakpoint regions with corresponding pattern information
     * @param faFilePath    reference genome file, used for fast local re-aligment
     * @return
     * @throws IOException 
     */
    public SequentialPatterns runAlgorithm(SequenceDatabase database, String rawPatternOut, String mergedPatternOut, 
            BufferedWriter susRegionWriter, BufferedWriter svRegionWriter, String faFilePath) throws IOException{
        MemoryLogger.getInstance().reset();
        
        this.database = database;
        startTime = System.currentTimeMillis();
                        
        System.out.println("Loading reference genome frome file ....");
        fileReader myReader = new fileReader();
        ReferenceSequenceFile refSeqFile = myReader.readFastaFile(faFilePath);
              
        // find frequent patterns
        prefixSpan(rawPatternOut);

        // doing post-processing of FSPMs
        merge merger = new merge(confAF, patternCount, chrIdxNameMap);
        stringMatcher matcher = new stringMatcher();
        merger.wgsPatternMerge(mergedPatternOut, matcher, database, patternCandidates, refSeqFile, svRegionWriter, susRegionWriter);
                        
        if (rawPatternWriter != null){
            rawPatternWriter.close();
        }
        endTime = System.currentTimeMillis();
        return patterns;
    }
    /**
     * Main function
     * @param outputPath output file of discovered frequent patterns
     * @throws IOException 
     */
    private void prefixSpan(String rawPatternOut) throws IOException{
        if (rawPatternOut == null){
            rawPatternWriter = null;            
        }
        // save results in file
        else{
            patterns = null;
            rawPatternWriter = new BufferedWriter(new FileWriter(rawPatternOut));
        }
        // Infomation of single superitem
        System.out.println("Collect information of single superitem .....");
        itemAppearMap = findSequencesContainItems(database);
        
        /**
         * Start creating initial pseudo-projected database.
         */
        System.out.println("Start creating initial pseudo sequence database ....");
        List<PseudoSequence> initialContext = new ArrayList<>();
        for (Sequence sequence : database.getSequences()){
            
            if (sequence.size() != 0){
                initialContext.add(new PseudoSequence(sequence, 0, 0, 0));
            }
        }
        for (Entry<String, List<ItemSeqIdentifier>> entry : itemAppearMap.entrySet()){
            String item = entry.getKey();
            Map<Integer, pseudoSequentialPattern> prefixTracker = new HashMap<>();
            if (entry.getValue().size() >= minsuppAbsolute){
                                
                List<ItemSeqIdentifier> itemAppearIdx = entry.getValue();
                List<PseudoSequence> projectedDatabase = buildProjectedContext(prefixTracker, item, initialContext, false);
                
                // Create a prefix with initial sequence ID 0
                SequentialPattern prefix = new SequentialPattern(0);
                prefix.addItemset(new Itemset(item));
                prefix.setItemAppear(itemAppearIdx);
                
                depthFirstRecursion(prefixTracker, prefix, projectedDatabase);
            }
        }
        
    }
    
    
    private void depthFirstRecursion(Map<Integer, pseudoSequentialPattern> prefixTracker, SequentialPattern prefix, List<PseudoSequence> psSeqDatabase) throws IOException{        
        Set<Pair> pairs = itemCountsInProjectedDB(prefix, psSeqDatabase);
        for (Pair pair : pairs){
            if (pair.getCount() >= minsuppAbsolute){
                SequentialPattern newPrefix;
                // If the frequent item is of form (_A), append it to the last itemset of the current prefix. 
                if(pair.isPostfix()){
                    newPrefix = appendItemToPrefixLastItemset(prefix, pair.getItem());
                }else{
                    newPrefix = appendItemToPrefixSequence(prefix, pair.getItem());
                }
                newPrefix.setItemAppear(pair.getItemAppear());
                savePatternCandidate(prefixTracker, newPrefix); 
//                if (newPrefix.size() > 1){
//                    List<SuperItem> superitems = newPrefix.getSuperItemsInPattern(database);                    
//                    
//                    for (int i = 0; i < superitems.size(); i++){   
//                        
//                        if (superitems.get(i).getPos() == 50227212){
//                            System.out.println(newPrefix.toString());
//                            System.out.println(superitems.size());
//                            System.out.println(superitems.get(i).toConciseString() + "\t" + i);     
//                        }   
//                        
//                    }                    
//                }  
               
                // Build pseudo-projected database of appended item.
                List<PseudoSequence> projectedDB = buildProjectedContext(prefixTracker, pair.getItem(), psSeqDatabase, pair.isPostfix());
                         
                
                // save patterns in file, not superitems
                savePatternToFile(newPrefix);
                prefixTracker = new HashMap<>();
                
                depthFirstRecursion(prefixTracker, newPrefix, projectedDB);
            }
        }
        MemoryLogger.getInstance().checkMemory();
    }
    
     
    /**
     * Pair is used to record and separate two conditions:
     * 1) (_A) -> (true, A)
     * 2) (A) -> (false, A)
     * Only count the item at start of the sequence, since the pattern growth in a continuous fashion.
     * @param prefixPattern
     * @param sequences
     * @return a set of pairs
     */
    private Set<Pair> itemCountsInProjectedDB(SequentialPattern prefixPattern, List<PseudoSequence> sequences){       
        Map<Pair, Pair> mapPairs = new HashMap<>();
        for (PseudoSequence sequence : sequences){
            Sequence oriSequence = database.getSequenceByID(sequence.getId());
            for (int j = 0; j < sequence.getSizeOfItemsetAt(0, oriSequence); j++){
                SuperItem si = sequence.getItemAtItemsetAt(j, 0, oriSequence);
                String item = si.getType();
                
                Pair paire = new Pair(sequence.isPostfix(0), item);
                Pair oldPaire = mapPairs.get(paire);
                if (oldPaire == null){
                    mapPairs.put(paire, paire);
                }
                /** 
                * same item found, use the previous one. 
                * previous pair object record the item appearance index.
                */
                else{
                    paire = oldPaire;
                }
                // Update item index for each item.
                paire.addItemAppearIdx(sequence.getId(), sequence.getFirstItemsetIdx(), 0, j);
            }
        }
        return mapPairs.keySet();
    }
    
    
    private Map<String, List<ItemSeqIdentifier>> findSequencesContainItems(SequenceDatabase database) {
        Map<String, List<ItemSeqIdentifier>> itemAppearMap = new HashMap<>();
        for(Sequence sequence : database.getSequences()){
            int sequenceID = sequence.getId();
            for(int i = 0; i < sequence.getItemsets().size();i++){
                List<SuperItem> itemset = sequence.getItemsets().get(i);
                for(int j = 0; j < itemset.size(); j ++){
                    SuperItem superitem = itemset.get(j);
                    ItemSeqIdentifier itemIdentity = new ItemSeqIdentifier(sequenceID, sequenceID, i, j);
                    String sitype = superitem.getType();
                    List<ItemSeqIdentifier> itemAppearIdx = itemAppearMap.get(sitype);
                    if(itemAppearIdx == null){
                        itemAppearIdx = new ArrayList<>();
                        itemAppearIdx.add(itemIdentity);
                        itemAppearMap.put(sitype, itemAppearIdx);
                    }else{
                        itemAppearMap.get(sitype).add(itemIdentity);
                    }
                }
            }
        }
        return itemAppearMap;
    }
    

    /**
     * Create a pseudo projected database. 
     * For the tree, root is empty. Level one is the initial projection, where you have to keep all possible suffix string of a give prefix.
     * 
     * @param item superitem type, but has to locate exact object while doing projection
     * @param database a pseudo database
     * @param isSuffix if the item is a suffix or not
     * @return 
     */
 
    private List<PseudoSequence> buildProjectedContext(Map<Integer, pseudoSequentialPattern> prefixTracker, String item, List<PseudoSequence> psSeqDatabase, boolean inSuffix){
        
        List<PseudoSequence> newPseudoProjectedDatabase = new ArrayList<>();
        
        for (PseudoSequence psSequence : psSeqDatabase) {
            // simple check if this is level one projection. Single item projection.
            if (psSequence.getOriSeqSize() == psSequence.getSize()){
                for (int i = 0; i < psSequence.getSize();i++){
                    int seqId = psSequence.getId();
                    int itemsetIdx = i;
                    Sequence oriSequence = database.getSequenceByID(seqId);
                    
                    int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                    
                    if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                        SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx, idxOfItemInItemset);                                               
                        pseudoSequentialPattern curPsPattern = prefixTracker.get(curSuperItem.getPos());
                        
                        if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx, idxOfItemInItemset + 1);
                            newSequence.setGenomeStartPos(curSuperItem.getPos());
//                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            
                            boolean nextSuperItemInRange = ableToBuildProjection(curPsPattern, i + 1, oriSequence, newSequence);
                            
                            if (curSuperItem.isARPsuperitem()){
//                                searchReadPairLinkedSuperItem(curSuperItem, itemsetIdx, seqId, psSequence, this.database);
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthReadPairLink(curSuperItem, si);
                                        if (isLinkable) {
                                            nextSuperItemInRange = isLinkable;
                                            break;
                                        }
                                    }                                                                        
                                }    
                            }
                                                        
                            if (newSequence.getSize() > 0 && nextSuperItemInRange){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                            
                        }
                        else if (itemsetIdx != psSequence.getSize() - 1){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx + 1, 0);
                            newSequence.setGenomeStartPos(curSuperItem.getPos());
//                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            
//                            boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, newSequence);
                            boolean nextSuperItemInRange = ableToBuildProjection(curPsPattern, i + 1, oriSequence, newSequence);
                            if (curSuperItem.isARPsuperitem()){
//                                searchReadPairLinkedSuperItem(curSuperItem, itemsetIdx, seqId, psSequence, this.database);
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthReadPairLink(curSuperItem, si);
                                        if (isLinkable){
                                            nextSuperItemInRange = isLinkable;
                                            break;
                                        }
                                    }
                                    
                                }                                                                
                            }
                            if (newSequence.getSize() > 0 && nextSuperItemInRange){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                            
                        }
                        
                    }
                }         
            }
            // In the deeper level of the tree, only build pseudo-sequence of the start item. Otherwise, the pattern is not consecutive.
            else{
                int seqId = psSequence.getId();
                int itemsetIdx = 0;
                Sequence oriSequence = database.getSequenceByID(seqId);

                int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                
                if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                    SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx + psSequence.getFirstItemsetIdx(), idxOfItemInItemset); 
                    
                    pseudoSequentialPattern curPsPattern = prefixTracker.get(curSuperItem.getPos());
//                    if (curPsPattern.getPatternRightMostPos() == 119703535){
//                        System.out.println(curPsPattern.toString(database));
//                    }
                    // If there exist more than one SuperItems at one position, go this condition.
                    if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
//                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);
                        
//                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);
                        boolean nextSuperItemInRange = ableToBuildProjection(curPsPattern, nextSuperItemIdx, oriSequence, psSequence);
                        
                        if (curSuperItem.isARPsuperitem()){
//                            searchReadPairLinkedSuperItem(curSuperItem, itemsetIdx, seqId, psSequence, this.database);
                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
                                int deeperSearchRange = nextSuperItemIdx + 3;
                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
                                    boolean isLinkable = patternGrowthReadPairLink(curSuperItem, si);
                                    if (isLinkable){
                                        nextSuperItemInRange = isLinkable;
                                        break;
                                    }
                                }

                            }                                                                
                        }
                        if (nextSuperItemInRange){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx, idxOfItemInItemset + 1);
                            if (newSequence.getSize() > 0){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                        }                 
                    }
                    // There is only one SuperItem at each genome position.
                    else if (itemsetIdx != psSequence.getSize() - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
//                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);
//                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);
                        boolean nextSuperItemInRange = ableToBuildProjection(curPsPattern, nextSuperItemIdx, oriSequence, psSequence);
//                        if (curSuperItem.isARPsuperitem()){
////                            searchReadPairLinkedSuperItem(curSuperItem, itemsetIdx, seqId, psSequence, this.database);
//                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
//                                int deeperSearchRange = nextSuperItemIdx + 3;
//                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
//                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
//                                    boolean isLinkable = patternGrowthReadPairLink(curSuperItem, si);
//                                    if (isLinkable){
//                                        nextSuperItemInRange = isLinkable;
//                                        break;
//                                    }
//                                }
//
//                            }                                                                
//                        }
                        
                        if (nextSuperItemInRange){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx + 1, 0);
                            if (newSequence.getSize() > 0){
                                newPseudoProjectedDatabase.add(newSequence);
                            }
                        }                 
                    }                   
                }                     
            }                        
        }
        return newPseudoProjectedDatabase;
    }
    
    
    /**
     * Genome start position of current sequence, limit the pattern span region.
     * The idea is that, I will check it while building the pseudo-projection.
     * 1) Once the genome cord of the left most superitem of the new projection minus the genome start cord of 
     * its corresponding initial projection is beyond the max region threshold, then this new pseudo-projection will be discarded.
     * 2) Though the distance calculate in 1) is beyond the threshold, if these two superitems are connected by read-pairs, I wont
     * discard the new appended superitem.
     * @return 
     */
    private boolean ableToBuildProjection(pseudoSequentialPattern psPattern, int nextItemIdx, Sequence oriSequence, PseudoSequence psSequence){
        boolean isNextItemInRange = false;
        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextItemIdx, 0);
        int genomeStartPos = psSequence.getGenomeStartPos();
        int spannedRegion = nextSuperItem.getPos() - genomeStartPos;
        
        if (spannedRegion <= patternSpanMaxRegion){
            isNextItemInRange = true;
        }
        boolean tmp = searchReadPairLinkedSuperItem(psPattern, nextItemIdx, oriSequence);
        return isNextItemInRange || tmp;
    }
    
    private boolean searchReadPairLinkedSuperItem(pseudoSequentialPattern psPattern, int nextItemIdx, Sequence oriSequence){
        
        if (psPattern == null){
            return false;
        }
        int oriSeqSize = oriSequence.size();
        List<SuperItem> superitems = psPattern.getSuperItemsOfPattern(database);                                        
        for (SuperItem curSi : superitems){

            if (curSi.isARPsuperitem()){
                int farestSearchRange = curSi.getSuperitemMateRegion().end;                                
                for (int k = nextItemIdx; k < oriSeqSize; k++){                    
                    SuperItem newSi = oriSequence.superItemAtPos(k, 0);
                    if (newSi.getPos() > farestSearchRange){
                        break;
                    }
                    boolean isLinkable = patternGrowthReadPairLink(curSi, newSi);
                    if (isLinkable){                        
                        return true;
                    }
                }
            }
        }
        return false;
    }
    /**
     * Aims at using read-pair info for pattern growth.
     * @param curSuperItem
     * @param nextSuperItem
     * @return 
     */
    private boolean patternGrowthReadPairLink(SuperItem curSuperItem, SuperItem nextSuperItem){
        superItemLink linker = new superItemLink();
        return linker.twoSuperItemLinkCheck(curSuperItem, nextSuperItem);
    }
    
    
    

    /**
     * Append to <(A)>
     * @param prefix
     * @param item
     * @return 
     */
    private SequentialPattern appendItemToPrefixSequence(SequentialPattern prefix, String item){
        SequentialPattern newPrefix = prefix.cloneSequence(); 
        newPrefix.addItemset(new Itemset(item)); 
        return newPrefix;
    }
    
    /**
     * Append to <(A_)>
     * @param prefix
     * @param item
     * @return 
     */
    private SequentialPattern appendItemToPrefixLastItemset(SequentialPattern prefix, String item){
        SequentialPattern newPrefix = prefix.cloneSequence();
        Itemset itemset = newPrefix.get(newPrefix.size() - 1);
        itemset.addItem(item);
        return newPrefix;
    }
    /**
     * Write frequent patterns to file if needed.
     * @param prefix
     * @throws IOException 
     */
    private void savePatternToFile(SequentialPattern prefix) throws IOException{
        // increase the pattern count
//        patternCount ++;
        if (rawPatternWriter != null){
            String outStr = prefix.patternDetailsString(database);
            rawPatternWriter.write(outStr);
            rawPatternWriter.newLine();
        }
    }
    
    /**
     * Save all patterns during pattern growth for further pattern merge process     
     * @param prefix     
     */
    private void savePatternCandidate(Map<Integer, pseudoSequentialPattern> prefixTracker, SequentialPattern prefix){
        int numOfTypes = prefix.getNumOfTypes();
        int patternLength = prefix.length();
                
        List<ItemSeqIdentifier> itemSeqIdentifiers = prefix.getItemAppear();
        for (ItemSeqIdentifier itemIdentity : itemSeqIdentifiers){
            List<pseudoSuperItem> curPattern = new ArrayList<>();
            int seqId = itemIdentity.getSeqID();
            int superitemSetStartIdx = itemIdentity.getSubSeqID() - patternLength + 1;     
            
            for (int i = 0; i < patternLength; i ++){
                int superitemSetIdx = superitemSetStartIdx + i;
                int length = database.getSequenceByID(seqId).getItemsets().get(superitemSetIdx).size();
                for(int j = 0; j < length;j ++){                       
                    pseudoSuperItem psItem = new pseudoSuperItem(seqId, superitemSetIdx, j);                        
                    psItem.setPsSuperitemLeftPos(database);
                    curPattern.add(psItem);
                   
                }                
            }
            
            pseudoSequentialPattern pattern = new pseudoSequentialPattern(curPattern, database);
            prefixTracker.put(pattern.getPatternRightMostPos(), pattern);
            
            if (patternLength >= 2 && numOfTypes == 1){
                int supportArps = numOfSupportARPs(pattern);
                if (supportArps > 1){
                    int chromIdx = pattern.ChromId;
                    while (patternCandidates.size() < chromIdx + 1){
                        patternCandidates.add(new ArrayList<>());
                    }
                    patternCandidates.get(chromIdx).add(pattern);
                    patternCount ++;
                }
            }else{

                int chromIdx = pattern.ChromId;
                while (patternCandidates.size() < chromIdx + 1){
                    patternCandidates.add(new ArrayList<>());
                }
                patternCandidates.get(chromIdx).add(pattern);
                patternCount ++;
            }                                

        }                       
    }
                               
   /**
     * For arp pattern of length larger than 2 with same items, check connections
     * @param pattern
     * @return 
     */
    public int numOfSupportARPs(pseudoSequentialPattern pattern){
        List<SuperItem> superitemList = pattern.getSuperItemsOfPattern(database);
        int length = superitemList.size();
        int maximuSup = 0;
        for (int i = 0; i < length; i++){
            int supportedARPs = 0;
            for (int j = 0; j < length; j ++){
                if (i != j){
                    SuperItem superitemOne = pattern.getSuperItemFromOriginal(database, i);
                    SuperItem superitemTwo = pattern.getSuperItemFromOriginal(database, j);

                    String[] qnameOneByteList = superitemOne.getQNames();
                    String[] qnameTwoByteList = superitemTwo.getQNames();
                    Set<String> uniqueQName = new HashSet<>();

                    
                    for (String qname : qnameOneByteList){

                        uniqueQName.add(qname);
                    }       
                    for (String qname : qnameTwoByteList){
                        if (uniqueQName.contains(qname)){
                            supportedARPs += 1;                
                        }
                        uniqueQName.add(qname);
                    }
                }
                
            }
            if (supportedARPs > maximuSup){
                maximuSup = supportedARPs;
            }
            
        }
        
        return maximuSup;
    }
    
   
    public void setParams(String[] idxChrMap, String maskFile) throws IOException{
        System.out.println("Loading default excludable regions ....");
        chrIdxNameMap = idxChrMap;   
        if (maskFile != null){
            setMaskRegion(maskFile);
            hasMaskRegion = true;
        }                
    }
    
    public void setMaskRegion(String maskFile) throws IOException{
        maskedRegion = new HashMap<>();
        FileInputStream fin;
        BufferedReader myInput;
                        
        fin = new FileInputStream(new File(maskFile));
        myInput = new BufferedReader(new InputStreamReader(fin));
        String thisLine;
        while((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split("\t");
            String chrName = tokens[0];
            int[] region = new int[]{Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2])};
            List<int[]> regions = maskedRegion.get(chrName);
            if (regions == null){
                regions = new ArrayList<>();
                regions.add(region);
                maskedRegion.put(chrName, regions);
            }else{
                regions.add(region);
            }
        }                
    }
    
    private boolean patternInMaskedRegion(pseudoSequentialPattern pattern){
        String chrom = chrIdxNameMap[pattern.ChromId];
//        if (!chrom.contains("chr")){
//            chrom = "chr" + chrom;
//        }
        List<int[]> regions = maskedRegion.get(chrom);
        for (int[] region : regions){
            if (pattern.patternLeftMostIntervalStart >= region[0] && pattern.patternLeftMostIntervalStart <= region[1]){
                return true;
            }
            else if(pattern.patternRightMostPos >= region[0] && pattern.patternRightMostPos <= region[1]){
                return true;
            }
        }
        return false;
    }
    public void printAlgoStatistics(){
        StringBuilder r = new StringBuilder(200);
        r.append("\n=============  C2Pspan =============\n Total time ~ ");
        r.append(endTime - startTime);
        r.append(" ms\n");     
        r.append(" Max memory (mb) : " );
        r.append(MemoryLogger.getInstance().getMaxMemory());        
        r.append('\n');
        r.append("===================================================\n");
        System.out.println(r.toString());
    }
    
        
}