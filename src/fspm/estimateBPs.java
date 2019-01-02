/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import dataStructure.SequenceDatabase;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import utils.svOutInfo;

/**
 *
 * @author jiadonglin
 */
public class estimateBPs {
        
    final private double minAF;
    public estimateBPs(double af){        
        minAF = af;
    }
    /**
     * Call SVs from patterns that are linked by ARP, cross or split align.
     * @param linkedPatternInfo
     * @param splitLinkPatternBuffer
     * @param mergedPatterns
     * @param database
     * @param idxNameMap
     * @param regionWriter
     * @throws IOException 
     */
    public void callSVFromLinked(Map<Integer, List<int[]>> linkedPatternInfo, Map<Integer, Integer> splitLinkPatternBuffer, List<pseudoSequentialPattern> mergedPatterns, 
            SequenceDatabase database, String[] idxNameMap, BufferedWriter regionWriter) throws IOException{
        

        // Call SVs from split align linked patterns
        for (Map.Entry<Integer, Integer> entry : splitLinkPatternBuffer.entrySet()){
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
                svInoInfo.writeVariantsOutput(regionWriter, idxNameMap[pattern.ChromId], sb);
            }
                       
        }
        // Call SVs from ARP linked patterns
        Set<Integer> linkedPatternCalled = new HashSet<>();
        for (Map.Entry<Integer, List<int[]>> entry : linkedPatternInfo.entrySet()){
            int linkType = 0;
            StringBuilder sb;
            int patternIdx = entry.getKey();
            
            List<int[]> matePatternInfo = entry.getValue();
            int[] maxMateInfo = getMateWithMostSup(matePatternInfo);
            int mateIdx = maxMateInfo[0];
            
            int sup = maxMateInfo[1];
            pseudoSequentialPattern pattern = mergedPatterns.get(patternIdx);

            if (pattern.getPatternLeftMostPos() == 105054356){
                System.out.println(pattern.toString(database));
            }

            if (linkedPatternCalled.contains(patternIdx)){
                continue;
            }      
            
            int[] splitAlignCoords = pattern.getSplitAlignCoords();                        
            int splitMate = pattern.getSplitStatus(database);
                       
            int[] crossSupInfo = pattern.getCrossSupInfo();
            int[] supEvi = new int[]{sup, -1, -1, pattern.getSplitSupCount(), pattern.getSplitReadMapQ()};
            
             // Self linked pattern (either ARP or split align)
//            if (mateIdx == patternIdx && pattern.selfLinkedSuperItemMapQCheck(20)){  
            if (mateIdx == patternIdx) {  
                int[] selfLinkedBP = pattern.selfLinkedPatternBP(database);
                // split align only
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
                    svInfo.writeVariantsOutput(regionWriter, idxNameMap[pattern.ChromId], sb);

                    linkedPatternCalled.add(patternIdx);                    
                }
                // self only
                else if (selfLinkedBP[0] > 0 && selfLinkedBP[1] > 0){
                    linkType = -1; 
                    
                    if (pattern.selfLinkedSuperItemAfCheck(minAF)){
                        sb = new StringBuilder();
                        svOutInfo svInfo = new svOutInfo(selfLinkedBP[0], selfLinkedBP[1], pattern.toTypeString(database), linkType, supEvi, pattern.getWeights(), 
                                pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
                        svInfo.setSelfLinkedInfo(pattern.getSelfLinkedItemMapQ(), pattern.getSelfLinkedItemAF(), pattern.getSelfLinkedItemType());
                        svInfo.writeVariantsOutput(regionWriter, idxNameMap[pattern.ChromId], sb);
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
                    svInfo.writeVariantsOutput(regionWriter, idxNameMap[pattern.ChromId], sb);
                                                    
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
//                    if (patternIdx == 29676 || mateIdx == 29676){
//                        System.out.println(pattern.toString(database));
//                    }
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
                    svInfo.writeVariantsOutput(regionWriter, idxNameMap[pattern.ChromId], sb);
                }                
                
            }            
        }                
    }
   /**
    * Call SVs from unlinkable patterns
    * @param unLinkedPattern
    * @param mergedPatterns
    * @param database
    * @param idxNameMap
    * @param refSeqFile
    * @param regionWriter
    * @param susRegionWriter
    * @throws IOException 
    */
    public void callSVFromUnlinked(Set<Integer> unLinkedPattern, List<pseudoSequentialPattern> mergedPatterns, SequenceDatabase database, String[] idxNameMap, 
           ReferenceSequenceFile refSeqFile, BufferedWriter regionWriter, BufferedWriter susRegionWriter) throws IOException{
        StringBuilder sb;
//        int localAlignedSV = 0;
        
        for (Integer id : unLinkedPattern){
            sb = new StringBuilder();
            pseudoSequentialPattern pattern = mergedPatterns.get(id); 
            int leftBound = pattern.patternLeftMostPos - 200;
            int rightBound = pattern.patternLeftMostPos + 200;
            if (pattern.getPatternLeftMostPos() == 19670779){
                System.out.println(pattern.toString(database));
            }  
            int[] supEvi = new int[]{0,-1, -1};
            String chrName = pattern.getChrName(pattern.ChromId, idxNameMap);
            
            int[] arpBasedEstimateInfo = pattern.unlinkedArpPatternPosEstimate(database, 20);
            int[] crossLinkInfo = pattern.getCrossSupInfo();                
            int[] oemEstimatePos = pattern.oemPatternPosEstimate(database);
            
            if (arpBasedEstimateInfo[0] > 0 && arpBasedEstimateInfo[1] > 0){
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
            
//            else if (!pattern.hasArpSuperItem() || pattern.patternRightMostPos - pattern.patternLeftMostPos <= 200){                
//
//                ReferenceSequence seq = refSeqFile.getSubsequenceAt(chrName, leftBound, rightBound);
//                String refStr = seq.getBaseString();
//                                                                                
//                // Do local alignment              
//                List<svOutInfo> svFromLocalAlign = pattern.doLocalAlign(database, refStr, pattern.patternLeftMostPos - 200);                  
//                                
//                if (!svFromLocalAlign.isEmpty()){
//                    for (svOutInfo sv : svFromLocalAlign){
//                        sb = new StringBuilder();
//                        sv.setSvInfo(pattern.toTypeString(database), pattern.getWeights(), pattern.getPatternSpanRegion(), pattern.getWeightRatio(), pattern.getOris());
//                        sv.writeVariantsOutput(regionWriter, chrName, sb);
//                    }
//                }               
//            }
            
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
}
