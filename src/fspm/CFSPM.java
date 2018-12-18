/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
//import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
//import htsjdk.samtools.util.Interval;

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
import utils.svOutInfo;
import utils.MemoryLogger;
import utils.Linkage;
import utils.stringMatcher;


/**
 *
 * @author jiadonglin
 */
public class CFSPM {
    
    private long startTime;
    private long endTime;
    
    final private int minsuppAbsolute;

    BufferedWriter freqPatternWriter = null;
    BufferedWriter idWriter;
    BufferedWriter intMeaWriter;
    stringMatcher strMatcher;
    
    boolean showPatternIdentifiers = false;
    boolean hasMaskRegion = false;
    
    /** original sequence count **/
//    private int sequenceCount;

    private int patternCount;
//    private int nonARPpatternCount = 0;
    private SequentialPatterns patterns = null;
    
    // Save all generate patterns during pattern growth
    final private List<List<pseudoSequentialPattern>> patternCandidates = new ArrayList<>();
       
    private ReferenceSequenceFile refSeqFile;
    String[] chrIdxNameMap;
    final private int patternSpanMaxRegion;
    private SequenceDatabase database;
    private Map<String, List<ItemSeqIdentifier>> itemAppearMap;
    private Map<String, List<int[]>> maskedRegion;
      
    
    public CFSPM(int minSup, int maxRegionSpan){
        this.minsuppAbsolute = minSup;
        this.patternSpanMaxRegion = maxRegionSpan;
       
    }
    
    
    /**
     * Main function
     * @param database  Superitem sequence database
     * @param freqPatternOut    discovered frequent abnormal patterns
     * @param mergedPatternOut  merged patterns
     * @param susRegionWriter   patterns indicate suspect regions waiting for further processing
     * @param svRegionWriter    predicted SV breakpoint regions with corresponding pattern information
     * @param faFilePath    reference genome file, used for fast local re-aligment
     * @return
     * @throws IOException 
     */
    public SequentialPatterns runAlgorithm(SequenceDatabase database, String freqPatternOut, String mergedPatternOut, 
            BufferedWriter susRegionWriter, BufferedWriter svRegionWriter, String faFilePath) throws IOException{
        MemoryLogger.getInstance().reset();
        
        this.database = database;
        startTime = System.currentTimeMillis();
                        
        System.out.println("Loading reference genome frome file ....");
        fileReader myReader = new fileReader();
        refSeqFile = myReader.readFastaFile(faFilePath);
              
        // find frequent patterns
        prefixSpan(freqPatternOut);

        // doing post-processing of FSPMs
        wgsPatternMerge(mergedPatternOut, svRegionWriter, susRegionWriter);
                        
        if (freqPatternWriter != null){
            freqPatternWriter.close();
        }
        endTime = System.currentTimeMillis();
        return patterns;
    }
    /**
     * Main function
     * @param outputPath output file of discovered frequent patterns
     * @throws IOException 
     */
    private void prefixSpan(String freqPatternOut) throws IOException{
        if (freqPatternOut == null){
            freqPatternWriter = null;            
        }
        // save results in file
        else{
            patterns = null;
            freqPatternWriter = new BufferedWriter(new FileWriter(freqPatternOut));
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
            
            if (entry.getValue().size() >= minsuppAbsolute){
                                
                List<ItemSeqIdentifier> itemAppearIdx = entry.getValue();
                List<PseudoSequence> projectedDatabase = buildProjectedContext(new SequentialPattern(0), item, initialContext, false);
                
                // Create a prefix with initial sequence ID 0
                SequentialPattern prefix = new SequentialPattern(0);
                prefix.addItemset(new Itemset(item));
                prefix.setItemAppear(itemAppearIdx);
                
                depthFirstRecursion(prefix, projectedDatabase);
            }
        }
        
    }
    
    
    private void depthFirstRecursion(SequentialPattern prefix, List<PseudoSequence> database) throws IOException{
        Set<Pair> pairs = itemCountsInProjectedDB(prefix, database);
        for (Pair pair : pairs){
            if (pair.getCount() >= minsuppAbsolute){
                SequentialPattern newPrefix;
                // If the frequent item is of form (_A), append it to the last itemset of the current prefix. 
                if(pair.isPostfix()){
                    newPrefix = appendItemToPrefixLastItemset(prefix, pair.getItem());
                }else{
                    newPrefix = appendItemToPrefixSequence(prefix, pair.getItem());
                }
                // Build pseudo-projected database of appended item.
                List<PseudoSequence> projectedDB = buildProjectedContext(newPrefix, pair.getItem(), database, pair.isPostfix());
                newPrefix.setItemAppear(pair.getItemAppear());
                                                       
                savePatternCandidate(newPrefix);                      
                // save patterns in file, not superitems
                savePatternToFile(newPrefix);
                depthFirstRecursion(newPrefix, projectedDB);
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
            Sequence oriSequence = this.database.getSequenceByID(sequence.getId());
            for (int j = 0; j < sequence.getSizeOfItemsetAt(0, oriSequence); j++){
                String item = sequence.getItemAtItemsetAt(j, 0, oriSequence).getType();
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
 
    private List<PseudoSequence> buildProjectedContext(SequentialPattern prefix, String item, List<PseudoSequence> database, boolean inSuffix){
        List<PseudoSequence> newPseudoProjectedDatabase = new ArrayList<>();
        
        for (PseudoSequence psSequence : database) {
            // simple check if this is level one projection. Single item projection.
            if (psSequence.getOriSeqSize() == psSequence.getSize()){
                for (int i = 0; i < psSequence.getSize();i++){
                    int seqId = psSequence.getId();
                    int itemsetIdx = i;
                    Sequence oriSequence = this.database.getSequenceByID(seqId);
                    
                    int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                    
                    if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                        SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx, idxOfItemInItemset);                                               
                                        
                        if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                            PseudoSequence newSequence = new PseudoSequence(psSequence, itemsetIdx, idxOfItemInItemset + 1);
                            newSequence.setGenomeStartPos(curSuperItem.getPos());
                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, newSequence);
                            
//                            if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                                boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                                nextSuperItemInRange = isLinkable;
//                            }
                            if (curSuperItem.isARPsuperitem()){
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
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
                            SuperItem nextSuperItem = oriSequence.superItemAtPos(i + 1, 0);
                            
                            boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, newSequence);
                            
//                            if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                                boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                                nextSuperItemInRange = isLinkable;
//                            }
                            if (curSuperItem.isARPsuperitem()){
                                if (i + 3 < psSequence.getSize()){
                                    int deeperSearchRange = i + 3;
                                    for (int k = i + 1; k < deeperSearchRange; k++){
                                        SuperItem si = oriSequence.superItemAtPos(k, 0);
                                        boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
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
                Sequence oriSequence = this.database.getSequenceByID(seqId);

                int idxOfItemInItemset = psSequence.indexOf(itemsetIdx, item, oriSequence);
                
                if (idxOfItemInItemset != -1 && psSequence.isPostfix(itemsetIdx) == inSuffix){
                    SuperItem curSuperItem = oriSequence.superItemAtPos(itemsetIdx, idxOfItemInItemset); 
                    
                    if (idxOfItemInItemset != psSequence.getSizeOfItemsetAt(itemsetIdx, oriSequence) - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);

                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);

                        if (curSuperItem.isARPsuperitem()){
                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
                                int deeperSearchRange = nextSuperItemIdx + 3;
                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
                                    boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
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
                    else if (itemsetIdx != psSequence.getSize() - 1){
                        int nextSuperItemIdx = psSequence.getFirstItemsetIdx() + 1;
                        SuperItem nextSuperItem = oriSequence.superItemAtPos(nextSuperItemIdx, 0);
                        boolean nextSuperItemInRange = ableToBuildProjection(nextSuperItem, psSequence);

//                        if (curSuperItem.isARPsuperitem() && nextSuperItem.isARPsuperitem()){
//                            boolean isLinkable = patternGrowthLinkage(curSuperItem, nextSuperItem);
//                            nextSuperItemInRange = isLinkable;
//                        }                                                
                        if (curSuperItem.isARPsuperitem()){
                            if (nextSuperItemIdx + 3 < psSequence.getSize()){
                                int deeperSearchRange = nextSuperItemIdx + 3;
                                for (int k = nextSuperItemIdx + 1; k < deeperSearchRange; k++){
                                    SuperItem si = oriSequence.superItemAtPos(k, 0);
                                    boolean isLinkable = patternGrowthLinkage(curSuperItem, si);
                                    if (isLinkable){
                                        nextSuperItemInRange = isLinkable;
                                        break;
                                    }
                                }

                            }                                                                
                        }
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
    private boolean ableToBuildProjection(SuperItem nextSuperItem, PseudoSequence psSequence){
        boolean isNextItemInRange = false;
        int genomeStartPos = psSequence.getGenomeStartPos();
        int spannedRegion = nextSuperItem.getPos() - genomeStartPos;
        
        if (spannedRegion <= patternSpanMaxRegion){
            isNextItemInRange = true;
        }
        
        return isNextItemInRange;
    }
    /**
     * Aims at using read-pair info for pattern growth.
     * @param curSuperItem
     * @param nextSuperItem
     * @return 
     */
    private boolean patternGrowthLinkage(SuperItem curSuperItem, SuperItem nextSuperItem){
        Linkage linker = new Linkage();
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
        if (freqPatternWriter != null){
            String outStr = prefix.patternDetailsString(database);
            freqPatternWriter.write(outStr);
            freqPatternWriter.newLine();
        }
    }
    
    /**
     * Save all patterns during pattern growth for further pattern merge process     
     * @param prefix     
     */
    private void savePatternCandidate(SequentialPattern prefix){
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
     * Merge all candidate patterns of whole genome result from savePatternCandidate()
     * @param linkerOut
     * @param svRegionOut
     * @throws IOException 
     */
    private void wgsPatternMerge(String mergedPatternOut, BufferedWriter regionWriter, BufferedWriter susRegionWriter) throws IOException{
        System.out.println("\nStart pattern post-processing, total candidate patterns: " + patternCount);        

//        BufferedWriter regionWriter = new BufferedWriter(new FileWriter(svRegionOut));       
        BufferedWriter mergedWriter;
        
        if (mergedPatternOut == null){
            mergedWriter = null;
        }else{
            mergedWriter = new BufferedWriter(new FileWriter(mergedPatternOut));
        }                
        
        strMatcher = new stringMatcher();          
        
        int numChrs = patternCandidates.size();
        for (int i = 0; i < numChrs; i ++){
            List<pseudoSequentialPattern> allPatterns = patternCandidates.get(i);
            // For single chrom, others wont be processed.
            if (allPatterns.isEmpty()){
                continue;
            }
            Map<Integer, List<Integer>> indexMap = getPatternStartIndexMap(allPatterns);
            List<pseudoSequentialPattern> mergedPatterns = oneChromMerge(i, allPatterns, indexMap);
            oneChrPatternLinkageAnalysis(mergedWriter, regionWriter, mergedPatterns, susRegionWriter);
        }
        
        regionWriter.close();
        if (susRegionWriter != null){
            susRegionWriter.close();
        }        
        
        if (mergedWriter != null){
            mergedWriter.close();
        }        
    }
    
    private Map<Integer, List<Integer>> getPatternStartIndexMap(List<pseudoSequentialPattern> arpPatterns){
        Map<Integer, List<Integer>> indexMap = new HashMap<>();
        int numOfPatterns = arpPatterns.size();
        for (int i = 0; i < numOfPatterns ; i++){
            int patternLeftMostPos = arpPatterns.get(i).patternLeftMostPos;
            List<Integer> indexList = indexMap.get(patternLeftMostPos);
            if (indexList == null) {
                indexList = new ArrayList<>();
                indexMap.put(patternLeftMostPos, indexList);
            }
            indexList.add(i);
        }
        return indexMap;
    }
    
    /**
     * Aims at merge patterns from single chrom
     * @param chrom
     * @param patterns
     * @param patternStartIndexMap
     * @return 
     */
    private List<pseudoSequentialPattern> oneChromMerge(int chrom, List<pseudoSequentialPattern> patterns, Map<Integer, List<Integer>> patternStartIndexMap){
        String chr = chrIdxNameMap[chrom];
        System.out.println("\nProcess Chr: " + chr +" pattern before merge: " + patternCandidates.get(chrom).size());
        List<pseudoSequentialPattern> mergedPatternCandidates = new ArrayList<>();
//        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartAndIndexMap.entrySet());
        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartIndexMap.entrySet());
        Collections.sort(patternIndexEntrys, new Comparator<Entry<Integer, List<Integer>>>(){
            @Override
            public int compare(Entry<Integer, List<Integer>> o1, Entry<Integer, List<Integer>> o2){
                return o1.getKey().compareTo(o2.getKey());
            }
        
        });
        
        int entrysSize = patternIndexEntrys.size();
        Set<Integer> tracker = new HashSet<>();
        for (int i = 0; i < entrysSize - 1; i ++){
            int candidateSize = mergedPatternCandidates.size();
            Entry<Integer, List<Integer>> entry = patternIndexEntrys.get(i);
            int pos = entry.getKey();
            
            List<Integer> patternIndex = entry.getValue();

            Entry<Integer, List<Integer>> nextEntry = patternIndexEntrys.get(i + 1);                
            List<Integer> nextPatternIndex = nextEntry.getValue();
            int nextPos = nextEntry.getKey();
            
            pseudoSequentialPattern mergedPattern = mergePatternList(patterns, patternIndex);
            pseudoSequentialPattern nextMergedPattern = mergePatternList(patterns, nextPatternIndex);
    
//            System.out.println(entry.getKey() + ": " + patternIndex.toString() + "\t" + mergedPattern.toString(database));
//            System.out.println(nextEntry.getKey() + ": " + nextPatternIndex.toString() + "\t" + nextMergedPattern.toString(database));
            List<pseudoSuperItem> mergedSuperItems = mergedPattern.mergeTwoPattern(nextMergedPattern, database);
            if (!mergedSuperItems.isEmpty()){
                tracker.add(pos);
                tracker.add(nextPos);
                
                pseudoSequentialPattern newMergedPattern = new pseudoSequentialPattern(mergedSuperItems, database);
    
//                System.out.println("merged: " + newMergedPattern.toString(database));

                // the new pattern might be merged with the last pattern in the candidate list.
                if (!mergedPatternCandidates.isEmpty()){
                    
                    List<pseudoSuperItem> superitems = secondaryMerge(mergedPatternCandidates, newMergedPattern);
                    if (!superitems.isEmpty()){
                        pseudoSequentialPattern secondaryMergedPattern = new pseudoSequentialPattern(superitems, database);
//                        System.out.println("Added pattern: " + secondaryMergedPattern.toString(database));
                        mergedPatternCandidates.remove(candidateSize - 1);
                        mergedPatternCandidates.add(secondaryMergedPattern);
                    }else{                        
//                        System.out.println("Added pattern: " + newMergedPattern.toString(database));
                        mergedPatternCandidates.add(newMergedPattern);
                    }
                }else{
//                    System.out.println("Added pattern: " + newMergedPattern.toString(database));

                    mergedPatternCandidates.add(newMergedPattern);
                }
                               
            }else{
                if (!tracker.contains(pos)){
                    tracker.add(pos);
                    if (! mergedPatternCandidates.isEmpty()){
                        List<pseudoSuperItem> superitems = secondaryMerge(mergedPatternCandidates, mergedPattern);
                        if (!superitems.isEmpty()){
                            pseudoSequentialPattern secondaryMergedPattern = new pseudoSequentialPattern(superitems, database);
                            mergedPatternCandidates.remove(candidateSize - 1);
                            mergedPatternCandidates.add(secondaryMergedPattern);
                        }else{
                            mergedPatternCandidates.add(mergedPattern);
                        }
                    }else{
                        mergedPatternCandidates.add(mergedPattern);
                    }                                        
                }                                                               
            }
                                
        }

        System.out.println("pattern after merge: " + mergedPatternCandidates.size());
        return mergedPatternCandidates;
    }
    /**
     * This is used to merge a new pattern with the last pattern in the merged pattern list
     * @param mergedPatternList
     * @param aPattern
     * @return 
     */
    private List<pseudoSuperItem> secondaryMerge(List<pseudoSequentialPattern> mergedPatternList, pseudoSequentialPattern aPattern){
        int candidateSize = mergedPatternList.size();
        pseudoSequentialPattern lastSPInCandidates = mergedPatternList.get(candidateSize - 1);
        List<pseudoSuperItem> superitems = lastSPInCandidates.mergeTwoPattern(aPattern, database);
        return superitems;
    }
    /**
     * Merged patterns of a chrom will be processed to generate SV calls. 
     * Read-pair and consensus string matching are implemented to make a confident call
     * @param linkWriter
     * @param regionWriter
     * @param mergedPatterns
     * @throws IOException 
     */
    private void oneChrPatternLinkageAnalysis(BufferedWriter mergedPatternWriter, BufferedWriter regionWriter, 
            List<pseudoSequentialPattern> mergedPatterns, BufferedWriter susRegionWriter) throws IOException{
                
        
        Linkage linkageAnalyzer = new Linkage();        
        Collections.sort(mergedPatterns);            
        int patternNums = mergedPatterns.size();
        
        Map<Integer, List<int[]>> linkedPatternInfo = new HashMap<>();
        // pattern link by split align
        Map<Integer, Integer> splitLinkPatternBuffer = new HashMap<>();
        // pattern without any link
        Set<Integer> unLinkedPattern = new HashSet<>();

        for (int i = 0; i < patternNums; i ++){
            
            pseudoSequentialPattern pattern = mergedPatterns.get(i);
            

            
            if (hasMaskRegion && patternInMaskedRegion(pattern)){
                continue;
            }

            if (mergedPatternWriter != null){
                mergedPatternWriter.write(pattern.toString(database));
                mergedPatternWriter.newLine();
            }           

            int[] splitAlignCoords = pattern.splitAlignForBP(database); 
            pattern.splitAlignCheck(database, mergedPatterns, splitLinkPatternBuffer, splitAlignCoords);
            pattern.crossLinkBpEstimate(database);
            int splitMate = pattern.getSplitStatus(database);
            boolean isCrossSup = pattern.isCrossSup();
                        
            // Patterns do not have ARP SuperItems
            if (!pattern.hasArpSuperItem()){                  
                if (splitMate != -3 && pattern.getSplitSupCount() > 1){                    
                    if (!splitLinkPatternBuffer.containsKey(splitMate)){
                        splitLinkPatternBuffer.put(i, splitMate);
                    }                    
                }else{
                    unLinkedPattern.add(i);
                }                
            }
            else{
                boolean isSelfLinked = pattern.isSelfLinked(database);   
                                                
                if (isSelfLinked || isCrossSup){   
                    int[] linkinfo = new int[2];
                    linkinfo[0] = i;
                    linkinfo[1] = pattern.numOfLinkedEvidence;
                    List<int[]> oldVal = linkedPatternInfo.get(i);
                    if (oldVal == null){
                        oldVal = new ArrayList<>();
                        oldVal.add(linkinfo);
                        linkedPatternInfo.put(i, oldVal);
                    }else{
                        oldVal.add(linkinfo);
                        linkedPatternInfo.put(i, oldVal);
                    }
                }    
                
                boolean arpSpanned = false;
                // Search mate pattern if it exists
                Map<Integer, Integer> indexMap = new HashMap<>();
                List<pseudoSequentialPattern> removedPatternCandidates = linkageAnalyzer.minusItemAndCopyWithIndexMap(mergedPatterns, i, indexMap);
                int[] linkedMateInfo = searchMatePattern(removedPatternCandidates, pattern, linkageAnalyzer);
                
                // Mate pattern is found through ARPs
                if (linkedMateInfo[0] != -1 && linkedMateInfo[1] != -1){    
                    arpSpanned = true;
                    List<int[]> oldVal = linkedPatternInfo.get(i);
                    int orignialIndex = indexMap.get(linkedMateInfo[0]); 
                    linkedMateInfo[0] = orignialIndex;
                    if (oldVal == null){
                        oldVal = new ArrayList<>();
                        oldVal.add(linkedMateInfo);
                        linkedPatternInfo.put(i, oldVal);
                    }else{
                        oldVal.add(linkedMateInfo);
                        linkedPatternInfo.put(i, oldVal);
                    }                                                                                                                                                  
                }
                
                if (!isSelfLinked && !isCrossSup && !arpSpanned){
                    if (splitAlignCoords[0] > 0 && splitAlignCoords[1] > 0){
                        splitLinkPatternBuffer.put(i, i);
                    }else{
                        unLinkedPattern.add(i);
                    }                    
                }                                             
            }
        }        
        
        callSVFromLinked(linkedPatternInfo, splitLinkPatternBuffer, mergedPatterns, regionWriter);
        callSVFromUnlinked(unLinkedPattern, mergedPatterns, regionWriter, susRegionWriter);
        
    }
    /**
     * Call SVs from unlinkable patterns
     * @param unLinkedPattern
     * @param mergedPatterns
     * @param regionWriter
     * @throws IOException 
     */
    private void callSVFromUnlinked(Set<Integer> unLinkedPattern, List<pseudoSequentialPattern> mergedPatterns, 
            BufferedWriter regionWriter, BufferedWriter susRegionWriter) throws IOException{
        StringBuilder sb;
//        int localAlignedSV = 0;
        
        for (Integer id : unLinkedPattern){
            sb = new StringBuilder();
            pseudoSequentialPattern pattern = mergedPatterns.get(id); 
            int leftBound = pattern.patternLeftMostPos - 200;
            int rightBound = pattern.patternLeftMostPos + 200;
            
            int[] supEvi = new int[]{0,-1, -1};
            String chrName = pattern.getChrName(pattern.ChromId, chrIdxNameMap);

                                                
            if (!pattern.hasArpSuperItem() || pattern.patternRightMostPos - pattern.patternLeftMostPos <= 200){
                
                int[] arpBasedEstimateInfo = pattern.unlinkedArpPatternPosEstimate(database, 20);
                int[] crossLinkInfo = pattern.getCrossSupInfo();                
                int[] oemEstimatePos = pattern.oemPatternPosEstimate(database);

                ReferenceSequence seq = refSeqFile.getSubsequenceAt(chrName, leftBound, rightBound);
                String refStr = seq.getBaseString();
                                                                                
                // Do local alignment              
                List<svOutInfo> svFromLocalAlign = pattern.doLocalAlign(database, refStr, pattern.patternLeftMostPos - 200);                  
                                
                if (!svFromLocalAlign.isEmpty()){
                    for (svOutInfo sv : svFromLocalAlign){
                        sb = new StringBuilder();
                        sv.setSvInfo(pattern.toTypeString(database), pattern.getWeights(), pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                        sv.writeVariantsOutput(regionWriter, chrName, sb);
                    }
                }
                // Predict the sv region based on unlinked pattern with ARPs
//                else if (arpBasedEstimateInfo[0] > 0 && arpBasedEstimateInfo[1] > 0 && arpBasedEstimateInfo[2] >= 10){
//
//                    supEvi[1] = arpBasedEstimateInfo[2];
//                    supEvi[2] = 0;
//                    svOutInfo svInfo = new svOutInfo(arpBasedEstimateInfo[0], arpBasedEstimateInfo[1], pattern.toTypeString(database), 9, supEvi, pattern.getWeights(), 
//                            pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
//                    svInfo.writeVariantsOutput(regionWriter, chrName, sb);
//                }
                // Predict SV regions with "smalmultiClippedPatternPosEstimatel_insert" ARPs
                else if (arpBasedEstimateInfo[0] > 0 && arpBasedEstimateInfo[1] > 0 && arpBasedEstimateInfo[2] == -2){
                    supEvi[1] = arpBasedEstimateInfo[2];
                    supEvi[2] = 0;
                    svOutInfo svInfo = new svOutInfo(arpBasedEstimateInfo[0], arpBasedEstimateInfo[1], pattern.toTypeString(database), -2, supEvi, pattern.getWeights(), 
                            pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                    
                    svInfo.writeVariantsOutput(regionWriter, chrName, sb);
                }
                
                else if (crossLinkInfo[0] > 0 && crossLinkInfo[1] > 0 && crossLinkInfo[2] >= 15){
                    supEvi[1] = crossLinkInfo[2];
                    supEvi[2] = crossLinkInfo[3];
                    svOutInfo svInfo = new svOutInfo(crossLinkInfo[0], crossLinkInfo[1], pattern.toTypeString(database), 
                            8, supEvi, pattern.getWeights(), pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                    svInfo.writeVariantsOutput(regionWriter, chrName, sb);
                }
                else if (oemEstimatePos[0] > 0 || oemEstimatePos[1] > 0) {
                                    
                    supEvi[1] = oemEstimatePos[2];
                    
                    svOutInfo svInfo = new svOutInfo(oemEstimatePos[0], oemEstimatePos[1], pattern.toTypeString(database), -4, supEvi, pattern.getWeights(), 
                        pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());

//                    System.out.println(svInfo.toString());
                    svInfo.writeVariantsOutput(regionWriter, chrName, sb);
                    
                }
                else if (susRegionWriter != null && pattern.isComplexPattern(database)){                    
                    svOutInfo svInfo = new svOutInfo(pattern.patternLeftMostPos, pattern.patternRightMostPos, pattern.toTypeString(database),
                        pattern.getWeights(), pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());                    
                    svInfo.writeSusRegion(susRegionWriter, chrName, sb);
                }
            }                                                  
        }        
    }

    /**
     * Call SVs from patterns that are linked by ARP, cross or split align.
     * @param linkedPatternInfo
     * @param splitLinkPatternBuffer
     * @param mergedPatterns
     * @param regionWriter
     * @throws IOException 
     */
    private void callSVFromLinked(Map<Integer, List<int[]>> linkedPatternInfo, Map<Integer, Integer> splitLinkPatternBuffer, List<pseudoSequentialPattern> mergedPatterns, 
            BufferedWriter regionWriter) throws IOException{
        
//        BufferedWriter notCalledPatternWriter = new BufferedWriter(new FileWriter("/Users/jiadonglin/SV_data/NA19238/morphlingv1/wgs.notCalledLinkedPatterns.txt"));

        // Call SVs from split align linked patterns
        for (Entry<Integer, Integer> entry : splitLinkPatternBuffer.entrySet()){
            int linkType = 3;
            int idx = entry.getKey(); 
            
            int splitAlignStatus = entry.getValue();  
            
            StringBuilder sb = new StringBuilder();
            pseudoSequentialPattern pattern = mergedPatterns.get(idx);         
            
            int[] coords = pattern.getSplitAlignCoords();
            if (coords[0] > 0 && coords[1] > 0){
                int[] crossSupInfo = pattern.getCrossSupInfo();
                int[] supEvi = new int[]{0,-1, -1, pattern.getSplitSupCount(), pattern.getSplitReadMapQ()};
                if (crossSupInfo[0] > 0 && crossSupInfo[1] > 0 && splitAlignStatus != -2){ 
                    coords[0] = crossSupInfo[0];
                    coords[1] = crossSupInfo[1];
                    supEvi[1] = crossSupInfo[2];
                    supEvi[2] = crossSupInfo[3];                    
                    linkType = 7;
                }
                svOutInfo svInoInfo = new svOutInfo(coords[0], coords[1], pattern.toTypeString(database), linkType, supEvi, pattern.getWeights(), 
                        pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                svInoInfo.writeVariantsOutput(regionWriter, chrIdxNameMap[pattern.ChromId], sb);
            }
                       
        }
        // Call SVs from ARP linked patterns
        Set<Integer> linkedPatternCalled = new HashSet<>();
        for (Entry<Integer, List<int[]>> entry : linkedPatternInfo.entrySet()){
            int linkType = 0;
            StringBuilder sb;
            int patternIdx = entry.getKey();
            
            List<int[]> matePatternInfo = entry.getValue();
            int[] maxMateInfo = getMateWithMostSup(matePatternInfo);
            int mateIdx = maxMateInfo[0];
            
            int sup = maxMateInfo[1];
            pseudoSequentialPattern pattern = mergedPatterns.get(patternIdx);
            if (linkedPatternCalled.contains(patternIdx)){
                continue;
            }      

            if (pattern.patternLeftMostPos == 8185330){
                System.out.println(pattern.toString(database));
            }
            
            
            int[] splitAlignCoords = pattern.getSplitAlignCoords();
            int[] selfLinkedBP = pattern.selfLinkedPatternBP(database);
            
            int splitMate = pattern.getSplitStatus(database);
                       
            int[] crossSupInfo = pattern.getCrossSupInfo();
            int[] supEvi = new int[]{sup, -1, -1, pattern.getSplitSupCount(), pattern.getSplitReadMapQ()};
            
             // Self linked pattern (either ARP or split align)
//            if (mateIdx == patternIdx && pattern.selfLinkedSuperItemMapQCheck(20)){  
            if (mateIdx == patternIdx) {  
                // split align
                if (splitMate == -2){
                    linkType = 5;
                    // split align && cross linked
                    if (crossSupInfo[0] > 0 && crossSupInfo[1] > 0){
                        supEvi[1] = crossSupInfo[2];
                        supEvi[2] = crossSupInfo[3];
                        linkType = 6;
                    }
                    sb = new StringBuilder();
                    svOutInfo svInfo = new svOutInfo(splitAlignCoords[0], splitAlignCoords[1], pattern.toTypeString(database), linkType, supEvi, pattern.getWeights(), 
                            pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                    svInfo.writeVariantsOutput(regionWriter, chrIdxNameMap[pattern.ChromId], sb);
                    linkedPatternCalled.add(patternIdx);
                }
                // self only
                else if (selfLinkedBP[0] > 0 && selfLinkedBP[1] > 0){
                    linkType = -1; 
                    
                    if (pattern.selfLinkedSuperItemAfCheck()){
                        sb = new StringBuilder();
                        svOutInfo svInfo = new svOutInfo(selfLinkedBP[0], selfLinkedBP[1], pattern.toTypeString(database), linkType, supEvi, pattern.getWeights(), 
                                pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                        svInfo.setSelfLinkedInfo(pattern.getSelfLinkedItemMapQ(), pattern.getSelfLinkedItemAF(), pattern.getSelfLinkedItemType());
                        svInfo.writeVariantsOutput(regionWriter, chrIdxNameMap[pattern.ChromId], sb);
                        linkedPatternCalled.add(patternIdx);
                    }
                    
                }
                // cross linked only
                else if (crossSupInfo[0] > 0 && crossSupInfo[1] > 0 && crossSupInfo[2] > 20){   
                    
                    linkType = 4;   
                    supEvi[1] = crossSupInfo[2];
                    supEvi[2] = crossSupInfo[3];
                                        
                    sb = new StringBuilder();
                    svOutInfo svInfo = new svOutInfo(crossSupInfo[0], crossSupInfo[1], pattern.toTypeString(database), linkType, supEvi, pattern.getWeights(), 
                            pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                    svInfo.writeVariantsOutput(regionWriter, chrIdxNameMap[pattern.ChromId], sb);
                    linkedPatternCalled.add(patternIdx);
                    
                }                
            }
            else{
                pseudoSequentialPattern matePattern = mergedPatterns.get(mateIdx);
                
//                System.out.println(matePattern.toString(database));
                int[] estBps = pattern.arpLinkedPatternBp(matePattern, database);
                
                if (pattern.arpSpanUseSplit){
                    // arp linked two pattern with split align support
                    linkType = 2;
                }else{
                    linkType = 1;
                }
                
                if (estBps[0] > 0 && estBps[1] > 0){                   
                    linkedPatternCalled.add(patternIdx);
                    linkedPatternCalled.add(mateIdx);
                    
                    List<Integer> weights = pattern.getWeights();                  
                    weights.addAll(matePattern.getWeights());
                                        
                    
                    List<Double> weightRatio = pattern.getWeightRatio();
                    weightRatio.addAll(matePattern.getWeightRatio());
                    
                    List<String> oris = pattern.getOris();
                    oris.addAll(matePattern.getOris());
                    
                    String patternStr = pattern.toTypeString(database) + "<>" + matePattern.toTypeString(database);
                    svOutInfo svInfo = new svOutInfo(estBps[0], estBps[1], patternStr, linkType, supEvi, weights, pattern.getPatternSpanRegion(), weightRatio,oris);
                    svInfo.setArpSpanInfo(pattern.getArpSpanMapQ(), pattern.getArpSpanAF(), pattern.getArpSpanItemType());
                    sb = new StringBuilder();
                    svInfo.writeVariantsOutput(regionWriter, chrIdxNameMap[pattern.ChromId], sb);
                }                
                
            }            
        }                
    }
    
    /**
     * A pattern might find more than one mates, mate with most read pair support is used.
     * @param mateInfo
     * @return 
     */
    private int[] getMateWithMostSup(List<int[]> mateInfo){
        if (mateInfo.size() == 1){
            return mateInfo.get(0);
        }
        int maxSup = 0;
        int maxSupIdx = -1;
        for (int[] ele : mateInfo){            
            if (ele[1] > maxSup){
                maxSup = ele[1];
                maxSupIdx = ele[0];
            }
        }
        int[] maxInfo = new int[2];
        maxInfo[0] = maxSupIdx;
        maxInfo[1] = maxSup;
        return maxInfo;
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
    
    
    /**
     * If a target pattern has ARP SuperItem, we can use it to search its mate pattern.Otherwise, we need to do consensus matching of itself.
     * @param sortedPatterns
     * @param targetPattern
     * @param linker
     * @return 
     */
      
    private int[] searchMatePattern(List<pseudoSequentialPattern> sortedPatterns, pseudoSequentialPattern targetPattern, Linkage linker){        
        List<QueryInterval> targetPatternMateInterval = targetPattern.superitemMateInterval;
        int length = sortedPatterns.size();
        int startIdx = 0;
        int endIdx = length - 1;
        int mateIndex = -1;
        int noQueryInterval = targetPatternMateInterval.size();
        
//        int linkSup = -1;
        int targetPatternMatchSuperItemIdx = -1;
        int matchedMatePatternSuperItemIdx = -1;
        
        // Two values to return, one is mate index and another one is number of supported read-pairs
        int[] returnVals = new int[]{-1, -1};                
        
        for (int i = 0; i < noQueryInterval ; i++){
            QueryInterval interval = targetPatternMateInterval.get(i);
            while (startIdx <= endIdx){
                int midIdx = startIdx + (endIdx - startIdx) / 2;

                pseudoSequentialPattern pattern = sortedPatterns.get(midIdx);
                if (!pattern.hasArpSuperItem()){
                    continue;
                }
                List<QueryInterval> sortedIntervals = pattern.superitemInterval;

                int overlapAt = hasOverlap(interval, sortedIntervals);
                if (overlapAt != -1){
                    mateIndex = midIdx;                    
                    targetPatternMatchSuperItemIdx = i;
                    matchedMatePatternSuperItemIdx = overlapAt;
                    break;
                }
                if (isAfterInterval(interval, sortedIntervals)){
                    startIdx = midIdx + 1;                    
                }
                if (isAheadInterval(interval, sortedIntervals)){
                    endIdx = midIdx - 1;
                }
                else if(overlapAt == -1 && !isAfterInterval(interval, sortedIntervals) && !isAheadInterval(interval, sortedIntervals)){
                    startIdx = midIdx + 1;
                }
            }
            if (mateIndex != -1){
                SuperItem superitemOne = targetPattern.getSuperItemOfPatternAtPos(database, targetPatternMatchSuperItemIdx);            
                pseudoSequentialPattern matchedSequentialPattern = sortedPatterns.get(mateIndex);
                SuperItem superitemTwo = matchedSequentialPattern.getSuperItemOfPatternAtPos(database, matchedMatePatternSuperItemIdx);
                boolean isEnoughARPs = linker.twoSuperItemLinkCheck(superitemOne, superitemTwo);                
                if (isEnoughARPs){
                    returnVals[0] = mateIndex;
                    returnVals[1] = linker.getSupLink();;
                    break;
                }
            }
            // reset start and end for next interval match
            startIdx = 0;
            endIdx = length - 1;
            
        }        
       
        return returnVals;
    }
    
    
    private pseudoSequentialPattern mergePatternList(List<pseudoSequentialPattern> arpPatterns, List<Integer> patternIndex){
        List<pseudoSequentialPattern> patterns = new ArrayList<>();
        for (Integer idx : patternIndex){
            patterns.add(arpPatterns.get(idx));
        }
        int maxLength = 0;
        int maxLengthPatternIndex = 0;
        int patternsSize = patterns.size();
        for (int i = 0; i < patternsSize; i ++){                
            if (patterns.get(i).patternLength > maxLength){
                maxLength = patterns.get(i).patternLength;
                maxLengthPatternIndex = i;
            }
        }
        return patterns.get(maxLengthPatternIndex);
    }
    private int hasOverlap(QueryInterval targetInterval, List<QueryInterval> intervals){
        int overlapAtIdx = -1;
        for (int i = 0; i < intervals.size(); i++){
            QueryInterval interval = intervals.get(i);
            if (interval != null){
                if (reciprocalOverlap(targetInterval, interval)){
                    overlapAtIdx = i;
                    break;
                }
//                if (overlap(targetInterval, interval)){
//                    overlapAtIdx = i;
//                    break;
//                }
            }
        }
        return overlapAtIdx;
    }
    
    private boolean reciprocalOverlap(QueryInterval a, QueryInterval b){
        int aSize = a.end - a.start;
        int bSize = b.end - b.start;
        boolean isOverlapped = false;
        
        if (b.start < a.end && b.start >= a.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }else if (a.start < b.end && a.start >= b.start){
            int overlapSize = a.end - b.start;
            double aOverlapRatio = (double) overlapSize / aSize;
            double bOverlapRatio = (double) overlapSize / bSize;
            if (aOverlapRatio >= 0.1 && bOverlapRatio >= 0.1){
                isOverlapped = true;
            }
        }
        
        return isOverlapped;
    }
    private boolean isAheadInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAhead = false;
        QueryInterval leftMostInterval = intervals.get(0);
        if (targetInterval.end < leftMostInterval.start){
            isAhead = true;
        }
        return isAhead;        
    }
    
    private boolean isAfterInterval(QueryInterval targetInterval, List<QueryInterval> intervals){
        boolean isAfter = false;
        QueryInterval lastInterval = intervals.get(intervals.size() - 1);
        if (targetInterval.start > lastInterval.end){
            isAfter = true;
        }
        return isAfter;
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