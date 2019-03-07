/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import dataStructure.SequenceDatabase;
import dataStructure.SuperItem;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import utils.superItemLink;
import strMatch.stringMatcher;

/**
 *
 * @author jiadonglin 
 */
public class merge {
    
    final private double confAF;
    int patternCount;
    String[] chrIdxNameMap;
    public merge(double af, int num, String[] idxNameMap){
        confAF = af;
        patternCount = num;
        chrIdxNameMap = idxNameMap;
    }
    /**
     * Merge patterns from whole genome
     * @param mergedPatternOut
     * @param strMatcher
     * @param database
     * @param patternCandidates
     * @param regionWriter
     * @param susRegionWriter
     * @throws IOException 
     */
    public void wgsPatternMerge(String mergedPatternOut, stringMatcher strMatcher, SequenceDatabase database, 
            List<List<pseudoSequentialPattern>> patternCandidates, ReferenceSequenceFile refSeqFile, BufferedWriter regionWriter, BufferedWriter susRegionWriter) throws IOException{
        System.out.println("\nStart pattern post-processing, total candidate patterns: " + patternCount);        

//        BufferedWriter regionWriter = new BufferedWriter(new FileWriter(svRegionOut));       
        BufferedWriter mergedWriter;
        
        if (mergedPatternOut == null){
            mergedWriter = null;
        }else{
            mergedWriter = new BufferedWriter(new FileWriter(mergedPatternOut));
        }                        
        strMatcher = new stringMatcher();          
        estimateBPs estimator = new estimateBPs(confAF);
        
        int numChrs = patternCandidates.size();
        for (int i = 0; i < numChrs; i ++){
            List<pseudoSequentialPattern> allPatterns = patternCandidates.get(i);
            // For single chrom, others wont be processed.
            if (allPatterns.isEmpty()){
                continue;
            }
            Map<Integer, List<Integer>> indexMap = getPatternStartIndexMap(allPatterns, database);
            List<pseudoSequentialPattern> mergedPatterns = oneChromMerge(i, database, allPatterns, patternCandidates, indexMap);
           
            oneChrPatternLinkageAnalysis(database, estimator, mergedWriter, regionWriter, mergedPatterns, refSeqFile, susRegionWriter);
        }
        
        regionWriter.close();
        if (susRegionWriter != null){
            susRegionWriter.close();
        }        
        
        if (mergedWriter != null){
            mergedWriter.close();
        }        
    }
    
     /**
     * Merge patterns at single chromosome.
     * @param chrom
     * @param patterns
     * @param patternStartIndexMap
     * @return 
     */
    private List<pseudoSequentialPattern> oneChromMerge(int chrom, SequenceDatabase database, List<pseudoSequentialPattern> patterns, 
            List<List<pseudoSequentialPattern>> patternCandidates, Map<Integer, List<Integer>> patternStartIndexMap){
        String chr = chrIdxNameMap[chrom];
        System.out.println("\nProcess Chr: " + chr +" pattern before merge: " + patternCandidates.get(chrom).size());
        List<pseudoSequentialPattern> mergedPatternCandidates = new ArrayList<>();
//        List<Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartAndIndexMap.entrySet());
        List<Map.Entry<Integer, List<Integer>>> patternIndexEntrys = new ArrayList<>(patternStartIndexMap.entrySet());
        Collections.sort(patternIndexEntrys, new Comparator<Map.Entry<Integer, List<Integer>>>(){
            @Override
            public int compare(Map.Entry<Integer, List<Integer>> o1, Map.Entry<Integer, List<Integer>> o2){
                return o1.getKey().compareTo(o2.getKey());
            }
        
        });
        
        int entrysSize = patternIndexEntrys.size();
        Set<Integer> tracker = new HashSet<>();
        for (int i = 0; i < entrysSize - 1; i ++){
            int candidateSize = mergedPatternCandidates.size();
            Map.Entry<Integer, List<Integer>> entry = patternIndexEntrys.get(i);
            int pos = entry.getKey();
            
            List<Integer> patternIndex = entry.getValue();

            Map.Entry<Integer, List<Integer>> nextEntry = patternIndexEntrys.get(i + 1);                
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
                if (newMergedPattern.patternLeftMostPos == 181074103){
                    System.out.println(newMergedPattern.toString(database));
                }
//                System.out.println("merged: " + newMergedPattern.toString(database));

                // the new pattern might be merged with the last pattern in the candidate list.
                if (!mergedPatternCandidates.isEmpty()){
                    
                    List<pseudoSuperItem> superitems = secondaryMerge(database, mergedPatternCandidates, newMergedPattern);
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
                        List<pseudoSuperItem> superitems = secondaryMerge(database, mergedPatternCandidates, mergedPattern);
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
    private List<pseudoSuperItem> secondaryMerge(SequenceDatabase database, List<pseudoSequentialPattern> mergedPatternList, pseudoSequentialPattern aPattern){
        int candidateSize = mergedPatternList.size();
        pseudoSequentialPattern lastSPInCandidates = mergedPatternList.get(candidateSize - 1);
        List<pseudoSuperItem> superitems = lastSPInCandidates.mergeTwoPattern(aPattern, database);
        return superitems;
    }
    
     /**
     * Merged patterns of a chromosome will be processed to generate SV calls.      
     * @param linkWriter
     * @param regionWriter
     * @param mergedPatterns
     * @throws IOException 
     */
    private void oneChrPatternLinkageAnalysis(SequenceDatabase database, estimateBPs estimator, BufferedWriter mergedPatternWriter, BufferedWriter regionWriter, 
            List<pseudoSequentialPattern> mergedPatterns, ReferenceSequenceFile refSeqFile, BufferedWriter susRegionWriter) throws IOException{
                
        
        superItemLink linkageAnalyzer = new superItemLink();        
        Collections.sort(mergedPatterns);            
        int patternNums = mergedPatterns.size();
        
        Map<Integer, List<int[]>> linkedPatternInfo = new HashMap<>();
        // pattern link by split align
        Map<Integer, Integer> splitLinkPatternBuffer = new HashMap<>();
        // pattern without any link
        Set<Integer> unLinkedPattern = new HashSet<>();

        for (int i = 0; i < patternNums; i ++){
            
            pseudoSequentialPattern pattern = mergedPatterns.get(i);
            if (pattern.patternLeftMostPos == 43593626){
                System.out.println(pattern.toString(database));
            }

            if (mergedPatternWriter != null){
                mergedPatternWriter.write(pattern.toString(database));
                mergedPatternWriter.newLine();
            }           
            
            
            // Check if the pattern contains split aligned info, if it exist we check:
            // If Split aligned position is the position of another SuperItem in current pattern.
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
                int[] linkedMateInfo = searchMatePattern(database, removedPatternCandidates, pattern, linkageAnalyzer);
                
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
        
        
        
        estimator.callSVFromLinked(linkedPatternInfo, splitLinkPatternBuffer, mergedPatterns, database, chrIdxNameMap, regionWriter);
        estimator.callSVFromUnlinked(unLinkedPattern, mergedPatterns, database, chrIdxNameMap, refSeqFile, regionWriter, susRegionWriter);
        
    }   
     /**
     * If a target pattern has ARP SuperItem, we can use it to search its mate pattern.Otherwise, we need to do consensus matching of itself.
     * @param sortedPatterns
     * @param targetPattern
     * @param linker
     * @return 
     */
      
    private int[] searchMatePattern(SequenceDatabase database, List<pseudoSequentialPattern> sortedPatterns, pseudoSequentialPattern targetPattern, superItemLink linker){        
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
    /**
     * Merge patterns in a list to a new pattern.
     * @param arpPatterns
     * @param patternIndex
     * @return 
     */
    
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
    private Map<Integer, List<Integer>> getPatternStartIndexMap(List<pseudoSequentialPattern> arpPatterns, SequenceDatabase database){
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
}
