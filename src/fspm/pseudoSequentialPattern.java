/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import strMatch.stringIdentity;
import strMatch.stringMatcher;
import dataStructure.SequenceDatabase;
import dataStructure.SuperItem;
import utils.*;
import htsjdk.samtools.QueryInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.HashSet;
import utils.svOutInfo;



/**
 *
 * @author jiadonglin
 */
public class pseudoSequentialPattern implements Comparable<pseudoSequentialPattern>{
    
    List<pseudoSuperItem> superitems = new ArrayList<>();
    List<pseudoSuperItem> ARPSuperItems = new ArrayList<>();
    List<pseudoSuperItem> MSSuperItems = new ArrayList<>();
    List<QueryInterval> superitemInterval;
    List<QueryInterval> superitemMateInterval;
    List<Integer> weights;
    List<Integer> postions;    
    List<Double> ratios;
    List<String> oris;   
    Map<QueryInterval, List<Integer>> indexMap = new HashMap<>();
    
    int ChromId = -1;    
    
    int patternLength;
    int numOfBPs;
    int patternLeftMostPos; // position of the first SuperItem in the pattern
    int patternRightMostPos; // position of the last SuperItem in the pattern
    int patternLeftMostIntervalStart; // read position of the left most interval of the pattern
    int[] linkedSuperItems = new int[2]; // for a self linked pattern, save its linked coord for boundary estimation
    int numOfLinkedEvidence = 0;
    int splitReadSup = 0;
    int splitReadMapQ = 0;
    
    int splitAlignLeftMatch = -1;
    int splitAlignRightMatch = -1;
    int splitAlignMate = -1;
    int splitAlignPairs; // number of SuperItem contains split-aligned reads
    int[] splitAlignCoords;
    int[] crossSupInfo;
    
    // Configurations of self-linked patterns
    int[] selfLinkedSuperItemMapQ; // linked ARP SuperItem mapping quality
    double[] selfLinkedSuperItemAF; // linked ARP SuperItem allele fraction
    String[] selfLinkedPatternBpItemType; // Breakpoint SuperItem type
    
    
    // Configurations of pair-wise pattern linkage
    int[] arpSpanBpMapQ;
    double[] arpSpanBpAF;
    String[] arpSpanBPItem;
    
   
    
    boolean arpSpanUseSplit = false;
    
    boolean hasSmallInsert = false;

    @Override
    public int compareTo(pseudoSequentialPattern other){                                          
        return patternLeftMostIntervalStart - other.patternLeftMostIntervalStart;
    }

    public List<QueryInterval> getSortedIntervals(){
        List<QueryInterval> sortedInterval = new ArrayList<>();
        for (QueryInterval val : superitemInterval){
            sortedInterval.add(val);
        }
        return sortedInterval;
    }

    public pseudoSequentialPattern(List<pseudoSuperItem> itemset, SequenceDatabase database, String clipped){
        superitems = itemset;
        patternLength = itemset.size();
        for (int i = 0; i < patternLength ; i ++){
            SuperItem superitem = superitems.get(i).getSuperItem(database);    
            if (ChromId == -1){
                ChromId = superitem.getChromIdx();
                break;
            }
        }  
    }

    public pseudoSequentialPattern (List<pseudoSuperItem> itemset, SequenceDatabase database){

        superitems = itemset;
        patternLength = itemset.size();
        superitemInterval = new ArrayList<>();
        superitemMateInterval = new ArrayList<>();
        postions = new ArrayList<>();
        weights = new ArrayList<>();
        ratios = new ArrayList<>();
        oris = new ArrayList<>();
        
        patternLeftMostPos = itemset.get(0).getSuperItem(database).getPos();
        patternRightMostPos = itemset.get(itemset.size() - 1).getSuperItem(database).getPos();

        for (int i = 0; i < patternLength ; i ++){
            SuperItem superitem = superitems.get(i).getSuperItem(database);   

//                System.out.println(superitem.toConciseString());
            if (superitem.getType().contains("SMALL")){
                hasSmallInsert = true;
            }
            weights.add(superitem.getWeight());            
            postions.add(superitem.getPos());
            ratios.add(superitem.getWeightRatio());
            oris.add(superitem.getReadOri());
            
            if (ChromId == -1){
                ChromId = superitem.getChromIdx();
            }
            if (superitem.isARPsuperitem() && ! superitem.getType().contains("OEM")){                      
                ARPSuperItems.add(superitems.get(i));                                         
                superitemInterval.add(superitem.getSuperitemRegion());
                superitemMateInterval.add(superitem.getSuperitemMateRegion());
            }else{
                MSSuperItems.add(superitems.get(i));                    
                numOfBPs += 1;
            }
        }   
        Collections.sort(ARPSuperItems);
        Collections.sort(superitemInterval);
        if (!superitemInterval.isEmpty()){
            patternLeftMostIntervalStart = superitemInterval.get(0).start;
        }

    } 
    public boolean hasSplitAlign(SequenceDatabase database){
        boolean hasSplit = false;
        if (splitAlignCoords[0] != 0 && splitAlignCoords[1] != 0){
            hasSplit = true;
        }
        return hasSplit;
    }
    
    public int[] getSplitAlignCoords(){
        return splitAlignCoords;
    }
    public int getPatternLeftMostPos(){
        return patternLeftMostPos;
    }
    public int getPatternRightMostPos(){
        return patternRightMostPos;
    }
    public boolean isCrossSup(){
        return crossSupInfo[0] > 0 && crossSupInfo[1] > 0 && crossSupInfo[2] > 30;
    }
    public int[] getCrossSupInfo(){
        return crossSupInfo;
    }
    public List<Integer> getIndex(QueryInterval aInterval){
        return indexMap.get(aInterval);
    }
    public List<Integer> getWeights(){
        return weights;
    }
    
    public List<Double> getWeightRatio(){
        return ratios;
    }
    public List<Integer> getPos(){
        return postions;
    }
    public List<String> getOris(){
        return oris;
    }
    public int getLinkSupport(){
        return numOfLinkedEvidence;
    }

    public int getNumOfBPs(){
        return numOfBPs;
    }
    public List<pseudoSuperItem> getSuperItems(){
        return superitems;
    }
    public boolean hasArpSuperItem(){
        return !ARPSuperItems.isEmpty();
    }
    public boolean hasClippedSuperItem(){
        return !MSSuperItems.isEmpty();
    }
    public pseudoSuperItem getSuperItemAt(int idx){
        return superitems.get(idx);
    }
    
    public int getSplitSupCount(){
        return splitReadSup;
    }
    
    public String getChrName(int dbIdx, String[] dbIdxToChr){
        return dbIdxToChr[dbIdx];
    }
    
    public int[] getSelfLinkedItemMapQ(){
        return selfLinkedSuperItemMapQ;
    }
    
    public double[] getSelfLinkedItemAF(){
        return selfLinkedSuperItemAF;
    }
    public String[] getSelfLinkedItemType(){
        return selfLinkedPatternBpItemType;
    }
    public int getSplitReadMapQ(){
        return splitReadMapQ;
    }
        
    public int[] getArpSpanMapQ(){
        return arpSpanBpMapQ;
    }
    public double[] getArpSpanAF(){
        return arpSpanBpAF;
    }
    public String[] getArpSpanItemType(){
        return arpSpanBPItem;
    }    
    
    /**
     * Return genome region that a pattern span
     * @return 
     */   
    public int[] getPatternSpanRegion(){
        return new int[]{patternLeftMostPos, patternRightMostPos};
    }
    
    
    public String toString(SequenceDatabase database){
        StringBuilder sb = new StringBuilder();
        
        for (pseudoSuperItem item : superitems){
            SuperItem superitem = item.getSuperItem(database);
            sb.append('(');
            sb.append(superitem.toConciseString());
            sb.append(')');
        }
        return sb.toString();
    }
    
    public String toTypeString(SequenceDatabase database){
        StringBuilder sb = new StringBuilder();
        for (pseudoSuperItem item : superitems){
            SuperItem superitem = item.getSuperItem(database);                                    
            sb.append(superitem.getType());                   
            sb.append(",");
        }

        return sb.substring(0, sb.length() - 1);
    }
    
    public List<SuperItem> getSuperItemsOfPattern(SequenceDatabase database){
        List<SuperItem> superItemsList = new ArrayList<>();
        for (pseudoSuperItem item : superitems){
            SuperItem superitem = item.getSuperItem(database);
            superItemsList.add(superitem);
        }
        return superItemsList;
    }
    

    
    public List<svOutInfo> doLocalAlign(SequenceDatabase database, String refSeqString, int refRegionLeft){
                
        stringMatcher strMatcher = new stringMatcher();
        List<String> mForStrings = new ArrayList<>(patternLength);
        StringBuilder sb;
        for (int i = 0; i < patternLength; i++){
            
            SuperItem superItem = superitems.get(i).getSuperItem(database);            
            sb = new StringBuilder();            
            if (superItem.getType().equals("MS")){
                sb.append(superItem.getMachedConsensus());
                sb.append(superItem.getClippedConsensus());
            }
            else if (superItem.getType().equals("SM")){
                sb.append(superItem.getClippedConsensus());
                sb.append(superItem.getMachedConsensus());
            }                        
            
            mForStrings.add(sb.toString());            
        }        
//        int[] alignCoords = new int[2];
        List<svOutInfo> svFromLocalAlign = new ArrayList<>();
        strMatcher.doStringAlign(mForStrings, superitems, refRegionLeft, refSeqString, svFromLocalAlign, database);               
        
        return svFromLocalAlign;
    }
    
    /**
     * pattern linked from split aligned reads
     * -3: inproper split align, split aligned pos either not in the pattern or in another pattern
     * -2: proper split align, where pos and split aligned pos are within the pattern spanned region
     * >0: split aligned pos is in another pattern
     * @param database
     * @return 
     */
    public int getSplitStatus(SequenceDatabase database){
        int splitStatus;    
        
        if (splitAlignLeftMatch >= 0 && splitAlignRightMatch >= 0){
            splitStatus = -2;
        }
        else if(splitAlignCoords[0] >= patternLeftMostPos && splitAlignCoords[1] <= patternRightMostPos){
            splitStatus = -2;
            if (splitAlignLeftMatch >= 0){
                pseudoSuperItem leftPSItem = superitems.get(splitAlignLeftMatch);
                String type = leftPSItem.getSuperItem(database).getType();                
                int closestToSplitAlign = Integer.MAX_VALUE;
                int rightSplitAlign = splitAlignCoords[1];               
                int rightClosestSuperItemPos = -1;
                for (int i = splitAlignLeftMatch; i < patternLength; i++){
                    pseudoSuperItem psItem = superitems.get(i);
                    SuperItem superItem = psItem.getSuperItem(database);
                    if (superItem.getType().equals(type)){
                        continue;
                    }
                    int dist = Math.abs(superItem.getPos() - rightSplitAlign);
                    if (dist < closestToSplitAlign){
                        closestToSplitAlign = dist;                        
                        rightClosestSuperItemPos = superItem.getPos();
                        
                    }
                }
                splitAlignCoords[1] = rightClosestSuperItemPos == -1 ? splitAlignCoords[1]:rightClosestSuperItemPos;
            }else if (splitAlignRightMatch >= 0){
                pseudoSuperItem rightPSItem = superitems.get(splitAlignRightMatch);
                String type = rightPSItem.getSuperItem(database).getType();
                int closestToSplitAlign = Integer.MAX_VALUE;
                int leftSplitAlign = splitAlignCoords[0];               
                int leftClosestSuperItemPos = -1;
                for (int i = splitAlignRightMatch; i > 0; i--){
                    pseudoSuperItem psItem = superitems.get(i);
                    SuperItem superItem = psItem.getSuperItem(database);
                    if (superItem.getType().equals(type)){
                        continue;
                    }
                    int dist = Math.abs(superItem.getPos() - leftSplitAlign);
                    if (dist < closestToSplitAlign){
                        closestToSplitAlign = dist;                        
                        leftClosestSuperItemPos = superItem.getPos();
                    }
                }
                splitAlignCoords[0] = leftClosestSuperItemPos == -1? splitAlignCoords[0]:leftClosestSuperItemPos;
            }
        }
        else if(splitAlignMate > 0){
            splitStatus = splitAlignMate;
        }else{
            splitStatus = -3;
        }               
        return splitStatus;
    }
    private boolean arpLinkedPatternSplitStatus(int[] splitCoords, pseudoSequentialPattern matePattern){
        int linkedPatternLeftPos = patternLeftMostPos < matePattern.patternLeftMostPos ? patternLeftMostPos : matePattern.patternLeftMostPos;
        int linkedPatternRightPos = patternRightMostPos > matePattern.patternRightMostPos ? patternRightMostPos : matePattern.patternRightMostPos;
        boolean splitStatus = false;
        if (splitCoords[0] > linkedPatternLeftPos && splitCoords[1] < linkedPatternRightPos){
            splitStatus = true;
        }
        return splitStatus;
    }
    /**
     * Check if coordinates of split aligned reads equals to Super-Item position in the pattern. 
     * If it does, record the index of such Super-Item
     * 
     * @param database
     * @param patternCandidates
     * @param patternIdx
     * @param splitAlignPos 
     */
    public void splitAlignCheck(SequenceDatabase database, List<pseudoSequentialPattern> patternCandidates, Map<Integer, Integer> patternIdx, int[] splitAlignPos){                         
        
        splitAlignCoords = splitAlignPos;
        
        for (int i = 0; i < superitems.size(); i++){            
            pseudoSuperItem psSuperItem = superitems.get(i);           
            SuperItem superItem = psSuperItem.getSuperItem(database);
            if (superItem.getPos() == splitAlignCoords[0]){
                splitAlignLeftMatch = i;
            }            
        }
        for (int i = superitems.size() - 1; i >= 0; i--){
            pseudoSuperItem psSuperItem = superitems.get(i);
            SuperItem superItem = psSuperItem.getSuperItem(database);

            if (superItem.getPos() == splitAlignCoords[1]){
                splitAlignRightMatch = i;
            }
        }
    }
    
    public boolean splitAlignLink(SequenceDatabase database, pseudoSequentialPattern matePattern, int[] expectedCoords){
//        boolean isSplitAlignLinked = false;               
        List<pseudoSuperItem> mateSuperItems = matePattern.getSuperItems();
        boolean onePartLinked = false;
        boolean otherPartLinked = false;
        if ((expectedCoords[1] > matePattern.patternLeftMostPos && expectedCoords[1] < matePattern.patternRightMostPos)
                || (expectedCoords[0] > matePattern.patternLeftMostPos && expectedCoords[0] < matePattern.patternRightMostPos)){
            onePartLinked = true;
        }        
        
        for (int i = 0; i < mateSuperItems.size(); i++){
            SuperItem item = mateSuperItems.get(i).getSuperItem(database);
            if (item.getSplitAlignPos() > expectedCoords[0] && item.getSplitAlignPos() < expectedCoords[1]){
                otherPartLinked = true;
            }
        }               
        return onePartLinked&&otherPartLinked;
    }
    /**
     * A pattern is self-linked by discordant read-pairs, indicates two ARP SuperItem linked.
     * Use ARP SuperItem for boundary estimation if there dose not exit split align or clipped reads
     * @param database
     * @return 
     */
    public boolean isSelfLinked(SequenceDatabase database){
        // Pattern linked by split alignment
        
        boolean linked = false;
        boolean hasEnoughARPs = false;
        superItemLink linker = new superItemLink();

        int Arps = ARPSuperItems.size();

        int curSuperItemIdx = -1;
        List<pseudoSuperItem> searchSpace;
        int machtedSuperItemIdx = -1;

        Set<Integer> matchedItem = new HashSet<>();
        int singleSuperItemARPcount = 0;
        for (int i = 0; i < Arps; i++){
            if (matchedItem.contains(i)){
                continue;
            }
            pseudoSuperItem target = ARPSuperItems.get(i);
            int siArpCount = singleSuperItemARPs(target, database);
            if (siArpCount > singleSuperItemARPcount){
                singleSuperItemARPcount = siArpCount;
            }
            
            Map<Integer, Integer> idxMap = new HashMap<>();
            searchSpace = minusSelf(ARPSuperItems, i, idxMap);

            int mateIndex = linker.mateSuperItemSearch(searchSpace, target, database);
            if (mateIndex != -1){    

                curSuperItemIdx = i;
                machtedSuperItemIdx = mateIndex;
                int originalIdx = idxMap.get(machtedSuperItemIdx);
                matchedItem.add(originalIdx);
                SuperItem superItemOne = ARPSuperItems.get(curSuperItemIdx).getSuperItem(database);
//                    SuperItem superItemTwo = searchSpace.get(machtedSuperItemIdx).getSuperItem(database);
                SuperItem superItemTwo = ARPSuperItems.get(originalIdx).getSuperItem(database);

                boolean isEnoughARPs = linker.twoSuperItemLinkCheck(superItemOne, superItemTwo);
                numOfLinkedEvidence = linker.getSupLink() > numOfLinkedEvidence ? linker.getSupLink():numOfLinkedEvidence;
                if (isEnoughARPs){
                    
                    hasEnoughARPs = true;
                    int superItemOnePos = superItemOne.getPos();
                    int superItemTwoPos = superItemTwo.getPos();
                    if (superItemTwoPos > superItemOnePos){
                        linkedSuperItems[0] = curSuperItemIdx;
                        linkedSuperItems[1] = idxMap.get(machtedSuperItemIdx);
                    }else{
                        linkedSuperItems[0] = idxMap.get(machtedSuperItemIdx);
                        linkedSuperItems[1] = curSuperItemIdx;
                    }
                }
            }
        }
        if (hasEnoughARPs){ 
            linked = true;              
        }
        else if (singleSuperItemARPcount > 2){
            linked = true;
            numOfLinkedEvidence = singleSuperItemARPcount;
        }
        
        return linked;
    }
    
    /**
     * Check if a SuperItem contain entire abnormal read-pair
     * @param ps
     * @param database
     * @return 
     */
    
    private int singleSuperItemARPs(pseudoSuperItem ps, SequenceDatabase database){
        SuperItem si = ps.getSuperItem(database);
        String[] qnames = si.getQNames();
        int readPairNums = readPairCounter(qnames);
        return readPairNums;
    }
    
    private List<pseudoSuperItem> minusSelf(List<pseudoSuperItem> psSuperItems, int idx, Map<Integer, Integer> idxMap){
        List<pseudoSuperItem> newSuperItems= new ArrayList<>();
        int length = psSuperItems.size();
        for (int i = 0; i < length; i ++){
            if (i != idx){
                newSuperItems.add(psSuperItems.get(i));
                idxMap.put(newSuperItems.size() - 1, i);
            }

        }
        return newSuperItems;
    }
    

    public SuperItem getSuperItemFromOriginal(SequenceDatabase database, int idx){
        pseudoSuperItem item = superitems.get(idx);
        return item.getSuperItem(database);
    }
    public SuperItem getSuperItemOfPatternAtPos(SequenceDatabase database, int idx){
        pseudoSuperItem item = ARPSuperItems.get(idx);
        return item.getSuperItem(database);
    }
    public List<pseudoSuperItem> mergeTwoPattern(pseudoSequentialPattern aPattern, SequenceDatabase database){          

        List<pseudoSuperItem> mergedSuperitems = new ArrayList<>();
        List<SuperItem> patternOneSuperItems = getSuperItemsOfPattern(database);
        List<SuperItem> patternTwoSuperItems = aPattern.getSuperItemsOfPattern(database);

        int lengthOne = patternOneSuperItems.size();
        int lengthTwo = patternTwoSuperItems.size();

        int matchedIndexAtPatternOne = -1;           
        int lastMatchedIndexAtPatternOne = lengthOne; 

        SuperItem patternTwoStartSuperItem = patternTwoSuperItems.get(0);
        SuperItem patternTwoLastSuperItem = patternTwoSuperItems.get(lengthTwo - 1);

        for (int i = 0; i < lengthOne ;i++){
            SuperItem patternOneSuperItem = patternOneSuperItems.get(i);
            if (patternTwoStartSuperItem.isEqual(patternOneSuperItem)){
                matchedIndexAtPatternOne = i;                    
            }
            if (patternTwoLastSuperItem.isEqual(patternOneSuperItem)){
                lastMatchedIndexAtPatternOne = i;
            }
        }

        if (matchedIndexAtPatternOne != - 1){
            if (lastMatchedIndexAtPatternOne == lengthOne){
                List<pseudoSuperItem> subListOfPatternOne = superitems.subList(0, matchedIndexAtPatternOne);
                mergedSuperitems.addAll(subListOfPatternOne);
                mergedSuperitems.addAll(aPattern.superitems);
            }else{
                mergedSuperitems = superitems;
            }                
        }

        return mergedSuperitems;            
    }
    /**
     * Estimate breakpoints of two arp linked patterns
     * @param matePattern
     * @param database
     * @return 
     */
    public int[] arpLinkedPatternBp(pseudoSequentialPattern matePattern, SequenceDatabase database){
        // First check the split alignment info on both side
        if (arpLinkedPatternSplitStatus(splitAlignCoords, matePattern)){
            arpSpanUseSplit = true;
            return splitAlignCoords;
        }
        int[] coords = new int[2];
        arpSpanBpMapQ = new int[2];
        arpSpanBpAF = new double[2];
        arpSpanBPItem = new String[2];
        
        // Process the patten ahead.
        for (int k = patternLength - 1; k >=0 ;k--){
            SuperItem superitem = getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                coords[0] = superitem.getPos();
                arpSpanBpMapQ[0] = superitem.getSumMapQ();
                arpSpanBpAF[0] = superitem.getRatio();
                arpSpanBPItem[0] = superitem.getType();

                arpLinkedBpHelper(false, 100, k, superitem.getPos(), superitem.getOri(), 
                        database, coords, arpSpanBpMapQ, arpSpanBpAF, arpSpanBPItem);                
            }
        }        

        // Process the mate pattern        
        int matePatternLength = matePattern.patternLength;
        for (int k = 0; k < matePatternLength; k ++){
            SuperItem superitem = matePattern.getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                coords[1] = superitem.getPos();
                arpSpanBpMapQ[1] = superitem.getSumMapQ();
                arpSpanBpAF[1] = superitem.getRatio();
                arpSpanBPItem[1] = superitem.getType();
                
                matePattern.arpLinkedBpHelper(true, 100, k, superitem.getPos(), superitem.getOri(), 
                        database, coords, arpSpanBpMapQ, arpSpanBpAF, arpSpanBPItem);                                
            }
        }
        if (coords[1] < coords[0]){
            int tmp = coords[0];
            coords[0] = coords[1];
            coords[1] = tmp;
        }
        return coords;
    }
    /**
     * Estimate ARP linked pattern bps.
     * @param isMate
     * @param expDist
     * @param arpItemIdx
     * @param arpItemPos
     * @param arpItemOri
     * @param database
     * @param coords
     * @param arpBpMapQ
     * @param aprBpAF
     * @param arpBpItem 
     */
    private void arpLinkedBpHelper(boolean isMate, int expDist, int arpItemIdx, int arpItemPos, String arpItemOri, 
            SequenceDatabase database, int[] coords, int[] arpBpMapQ, double[] aprBpAF, String[] arpBpItem){
        
        if (arpItemOri.equals("+")){
            if (arpItemIdx < patternLength){
                for (int i = arpItemIdx; i < patternLength; i ++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (superitem.getType().equals("MS")){
                        int dist = Math.abs(superitem.getPos() - arpItemPos);                        
                        if (dist < expDist){
                            expDist = dist;
                             
                            if (!isMate){
                                coords[0] = superitem.getPos(); 
                                arpBpMapQ[0] = superitem.getSumMapQ();
                                aprBpAF[0] = superitem.getRatio();
                                arpBpItem[0] = superitem.getType();
                            }else{
                                coords[1] = superitem.getPos(); 
                                arpBpMapQ[1] = superitem.getSumMapQ();
                                aprBpAF[1] = superitem.getRatio();
                                arpBpItem[1] = superitem.getType();
                            }                 
 
                                                      
                        }
                        break;
                    }
                }
            }
            if(arpItemIdx > 0){
                for (int i = arpItemIdx; i>=0 ; i--){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (superitem.getType().equals("MS")){
                        int dist = Math.abs(superitem.getPos() - arpItemPos);
                        if (dist < expDist){                            
                            if (!isMate){
                                coords[0] = superitem.getPos();
                                arpBpMapQ[0] = superitem.getSumMapQ();
                                aprBpAF[0] = superitem.getRatio();
                                arpBpItem[0] = superitem.getType();
                            }else{
                                coords[1] = superitem.getPos();
                                arpBpMapQ[1] = superitem.getSumMapQ();
                                aprBpAF[1] = superitem.getRatio();
                                arpBpItem[1] = superitem.getType();
                            }
                            break;
                        }
                    }
                }
            }
        }
        else{    
            if(arpItemIdx > 0){
                for (int i = arpItemIdx; i>=0 ; i--){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (superitem.getType().equals("SM")){
                        int dist = Math.abs(superitem.getPos() - arpItemPos);
                        if (dist < expDist){                                                                                    
                            expDist = dist;                            
                            if (!isMate){
                                coords[0] = superitem.getPos();
                                arpBpMapQ[0] = superitem.getSumMapQ();
                                aprBpAF[0] = superitem.getRatio();
                                arpBpItem[0] = superitem.getType();
                            }else{
                                coords[1] = superitem.getPos();
                                arpBpMapQ[1] = superitem.getSumMapQ();
                                aprBpAF[1] = superitem.getRatio();
                                arpBpItem[1] = superitem.getType();
                            }
                            break;
                        }
                    }
                }
            }            
            if (arpItemIdx < patternLength){
                for (int i = arpItemIdx; i < patternLength; i ++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (superitem.getType().equals("SM")){
                        int dist = Math.abs(superitem.getPos() - arpItemPos);                        
                        if (dist < expDist){
                            
                            if (!isMate){
                                coords[0] = superitem.getPos();
                                arpBpMapQ[0] = superitem.getSumMapQ();
                                aprBpAF[0] = superitem.getRatio();
                                arpBpItem[0] = superitem.getType();
                            }else{
                                coords[1] = superitem.getPos();
                                arpBpMapQ[1] = superitem.getSumMapQ();
                                aprBpAF[1] = superitem.getRatio();
                                arpBpItem[1] = superitem.getType();
                            }                     
                        }
                        break;
                    }
                }
            }
            
        }
        
    }
     
    
    public void crossLinkBpEstimate(SequenceDatabase database){    
        // info[0], info[1] start and end pos. 
        // info[2] shared string length
        // info[3] number of reads share string
        int[] info = new int[4];
        
        StringBuilder sb;
        stringMatcher strMatcher = new stringMatcher();
        List<String> mStrings = new ArrayList<>();
        List<String> sForwardStrings = new ArrayList<>();
        List<String> sReverseStrings = new ArrayList<>();

        int allStrLength = 0;
        for (pseudoSuperItem psItem : superitems){
            SuperItem superItem = psItem.getSuperItem(database);
            String matchedStr = superItem.getMachedConsensus();
            String clippedStr = superItem.getClippedConsensus();
            sb = new StringBuilder(clippedStr);
            String clippedRevereStr = sb.reverse().toString();

            allStrLength += matchedStr.length();
            allStrLength += clippedStr.length();
            allStrLength += clippedRevereStr.length();

            mStrings.add(matchedStr);
            sForwardStrings.add(clippedStr);
            sReverseStrings.add(clippedRevereStr);                
        }      
        
        strMatcher.strCrossMatch(mStrings, sForwardStrings, sReverseStrings);
                
        int numOfStrs = superitems.size();
        double avgLen = (double)allStrLength/superitems.size();
        double val = Math.log(numOfStrs * avgLen) / Math.log(4);
        int expect = (int) Math.ceil(val);
        
        int[] linkInfo = strMatcher.isCrossLinked();
        int infoSum = linkInfo[0] + linkInfo[1] + linkInfo[2];
        
        if (infoSum > 0){     
            int observedLen = strMatcher.isForwardExist() ? strMatcher.getForwardSharedStrLength():strMatcher.getReverseSharedStrLength();
            
            if (observedLen > expect){
//                System.out.println("Expect: " + expect + "\tObserve: " + observedLen); 
//                strMatcher.printSharedString(superitems, database);
//                System.out.println("\n");
                List<stringIdentity> identitys = strMatcher.isForwardExist()?strMatcher.getForwardSharedStrIdentitys():strMatcher.getReverseSharedStrIdentitys();
                if (validClipLink(identitys, database)){
                    int[] coords = strMatcher.getEstimateBp(superitems, database);
                    if (coords[0] != coords[1]){
                        info[0] = coords[0];
                        info[1] = coords[1];
                        info[2] = observedLen;
                        info[3] = infoSum;
                    }
                    
                }

            }            
        }
        crossSupInfo = info;
        
    }
    
    private boolean validClipLink(List<stringIdentity> strIdentitys, SequenceDatabase database){
        List<String> types = new ArrayList<>();
        boolean isValid = false;
        Collections.sort(strIdentitys, new Comparator<stringIdentity>(){
        @Override
            public int compare(stringIdentity o1, stringIdentity o2){
                return o1.getSeqId() - o2.getSeqId();
            }        
        });
        
        for (stringIdentity id : strIdentitys){
            pseudoSuperItem psItem = superitems.get(id.getSeqId());
            SuperItem superItem = psItem.getSuperItem(database);
            types.add(superItem.getType());
        }
        if (types.size() > 1){            
            String firstType = types.get(0);
            String lastType = types.get(types.size() - 1);
            if ((firstType.equals("MS") && lastType.equals("SM")) || (firstType.equals("SM") && lastType.equals("MS"))){
                isValid = true;
            }
            isValid = true;
        }        
        return isValid;
    }
    
    public int[] splitAlignForBP(SequenceDatabase database){
        int[] pos = new int[2];       
        for (int i = 0;i < MSSuperItems.size();i++){
            SuperItem superitem = MSSuperItems.get(i).getSuperItem(database);
//            String superItemType = superitem.getType();
            int superItemPos = superitem.getPos();
            int splitAlignPos = superitem.getSplitAlignPos();            
            if (splitAlignPos != -1 && superItemPos != splitAlignPos){
                                               
                if (Math.abs(superItemPos - splitAlignPos) > 1000000){
                    continue;
                }
                if (superitem.getSplitReadCount() > splitReadSup){
                    splitReadMapQ = superitem.getSplitReadMapQ();
                    splitReadSup = superitem.getSplitReadCount();
                    
                    if (superItemPos < splitAlignPos){
                        pos[0] = superItemPos;
                        pos[1] = splitAlignPos;
                    }else{
                        pos[0] = splitAlignPos;
                        pos[1] = superItemPos;
                    }                    
                    
                }
                
                // If we find the split align within the pattern defined range, return.
//                if (pos[0] >= patternLeftMostPos && pos[1] <= patternRightMostPos){                                        
//                    splitAlignPairs += 1;
//                }                   
                             
            }
            
        }
        
        return pos;
    }
    /**
     * Check self-linked pattern SuperItem allele fraction.
     * @return 
     */
    public boolean selfLinkedSuperItemAfCheck(double minAF){     
                
        if (selfLinkedSuperItemAF[0] >= minAF || selfLinkedSuperItemAF[1] >= minAF){
            return true;
        }        
        return false;
    }
    
    /**
     * Estimate breakpoints of a self-linked pattern
     * @param database
     * @return 
     */
    public int[] selfLinkedPatternBP(SequenceDatabase database){
        int[] pos = new int[]{-1,-1};
//        int[] splitAlignPos = splitAlignForBP(database);
        
        SuperItem leftARPItem = ARPSuperItems.get(linkedSuperItems[0]).getSuperItem(database);
        int leftARPpatternSuperItemPos = leftARPItem.getPos();
        SuperItem rigthARPItem = ARPSuperItems.get(linkedSuperItems[1]).getSuperItem(database);
        int rightARPpatternSuperItemPos = rigthARPItem.getPos();
        
        selfLinkedSuperItemMapQ = new int[2];
        selfLinkedSuperItemMapQ[0] = leftARPItem.getSumMapQ();
        selfLinkedSuperItemMapQ[1] = rigthARPItem.getSumMapQ();
        
        selfLinkedSuperItemAF = new double[2];
        selfLinkedSuperItemAF[0] = leftARPItem.getRatio();
        selfLinkedSuperItemAF[1] = rigthARPItem.getRatio();
        
        selfLinkedPatternBpItemType = new String[2];
        selfLinkedPatternBpItemType[0] = leftARPItem.getType();
        selfLinkedPatternBpItemType[1] = rigthARPItem.getType();
        
        int leftClosestToARP = -1;
        int rightClosestToARP = -1;
        int maxDistToLeftARP = 500;
        int maxDistToRightARP = 500;
        
        // Get the closest clipped SuperItem if it exist.
        for (int i = 0; i < patternLength; i++){
            SuperItem superitem = getSuperItemFromOriginal(database, i);
            
            boolean isARPSuperItem = superitem.isARPsuperitem();
            String SuperItemType = superitem.getType();
            int SuperItemPos = superitem.getPos();
            
            if (!isARPSuperItem && SuperItemType.equals("MS")){
                int dist = Math.abs(SuperItemPos - leftARPpatternSuperItemPos);
                if (dist < maxDistToLeftARP){
                    maxDistToLeftARP = dist;
                    leftClosestToARP = i;
                }
            }
            if (!isARPSuperItem && SuperItemType.equals("SM")){
                int dist = Math.abs(SuperItemPos - rightARPpatternSuperItemPos);
                if (dist < maxDistToRightARP){
                    maxDistToRightARP = dist;
                    rightClosestToARP = i;
                }
            }
        }
        
        if (leftClosestToARP != -1 && rightClosestToARP != -1){
            pos[0] = superitems.get(leftClosestToARP).getSuperItem(database).getPos();
            pos[1] = superitems.get(rightClosestToARP).getSuperItem(database).getPos();
            selfLinkedPatternBpItemType[0] = superitems.get(leftClosestToARP).getSuperItem(database).getType();
            selfLinkedPatternBpItemType[1] = superitems.get(leftClosestToARP).getSuperItem(database).getType();
            
        }        
        else if (leftClosestToARP != -1){            
            SuperItem leftSuperItem = getSuperItemFromOriginal(database, leftClosestToARP);
            pos[0] = leftSuperItem.getPos();
            selfLinkedPatternBpItemType[0] = leftSuperItem.getType();
            pos[1] = leftSuperItem.getSplitAlignPos();
            // Left clipped SuperItem dose not have supplementary alignments.
            if (pos[1] == -1 && rightClosestToARP > leftClosestToARP){
                SuperItem rightSuperItem = getSuperItemFromOriginal(database, rightClosestToARP);                
                pos[1] = rightSuperItem.getPos();
                int rightSuperItemSplitPos = rightSuperItem.getSplitAlignPos();
                
                // Refine the position if right clipped SuperItem has supplementary alignments and left
                // position is not decided yet.
                if (rightSuperItemSplitPos != -1 && pos[0] == 0){
                    pos[0] = rightSuperItemSplitPos;
                }
            }else if (rigthARPItem.getWeight() > 5 && pos[1] == -1){
                pos[1] = rightARPpatternSuperItemPos;
                
            }
        }
        else if (rightClosestToARP != -1){
            SuperItem rightSuperItem = getSuperItemFromOriginal(database, rightClosestToARP);
            pos[1] = rightSuperItem.getPos();
            selfLinkedPatternBpItemType[1] = rightSuperItem.getType();
            pos[0] = rightSuperItem.getSplitAlignPos();
            if (pos[0] == -1 && (leftClosestToARP != -1 && leftClosestToARP < rightClosestToARP)){
                SuperItem leftSuperItem = getSuperItemFromOriginal(database, leftClosestToARP);
                pos[0] = leftSuperItem.getPos();
                int leftSuperItemSplitPos = leftSuperItem.getSplitAlignPos();
                
                if (leftSuperItemSplitPos != -1 && pos[1] == 0){
                    pos[1] = leftSuperItemSplitPos;
                }
                
            }else if (leftARPItem.getWeight() > 5 && pos[0] == -1){
                pos[0] = leftARPpatternSuperItemPos;
            }
        }
        
        if (pos[0] == -1 || pos[1] == -1){
            
            pos[0] = leftARPpatternSuperItemPos;
            pos[1] = rightARPpatternSuperItemPos;                        
        }      
        
        return pos;
        
    }  
    /**
     * 
     * @param database
     * @param minQ
     * @return 
     */
    public int[] estimateBpFromArp(SequenceDatabase database, int minQ){
        // Use split-alignment to get the BP pos if there exist split align of a SuperItem.         
        int[] info = new int[]{-1, -1, -1};
        int arpsNum = ARPSuperItems.size();
        // Used to count number of whole read-pairs within a SuperItem.
        int maxReadPairIdx = -1;
        int maxReadPairNum = 2;
        if (arpsNum > 0){  
            for (int i = 0; i < arpsNum; i++){
                SuperItem si = ARPSuperItems.get(i).getSuperItem(database);
                if (si.getSumMapQ() < si.getWeight() * minQ && si.getRatio() > 0.2){
                    continue;
                }
                String[] qnames = si.getQNames();
                int readPairNums = readPairCounter(qnames);
                if (readPairNums > maxReadPairNum){
                    maxReadPairNum = readPairNums;
                    maxReadPairIdx = i;
                }
            }            
        }
        if (maxReadPairIdx != -1){
            QueryInterval interval = ARPSuperItems.get(maxReadPairIdx).getSuperItem(database).getSuperitemRegion();
            QueryInterval mateInterval = ARPSuperItems.get(maxReadPairIdx).getSuperItem(database).getSuperitemMateRegion();
            info[0] = interval.start;
            info[1] = mateInterval.end;
        }
       
        return info;
    }  
    
    private int readPairCounter(String[] qnames){
        int count = 0;
        Set<String> qSet = new HashSet<>();
        for (String q : qnames){
            if (qSet.contains(q)){
                count += 1;
            }else{
                qSet.add(q);
            }
        }
        return count;
    }
    
    /**
     * For patterns contain ARP superitem but do not have link info.
     * @param database
     * @param minQ
     * @return 
     */
    public int[] unlinkedArpPatternPosEstimate(SequenceDatabase database, int minQ){
        int[] info = new int[3];
        // If the pattern contains 'ARP_SMALL_INSERT', return the position as left and right pos of the pattern
        info[0] = patternLeftMostPos;
        info[1] = patternRightMostPos;
        if (hasSmallInsert){
            info[2] = -2;
        }
                
        if (!hasSmallInsert){
            info = estimateBpFromArp(database, minQ);
        }
        
        return info;
    }
    /**
     * OEM reads are important sign for potential SVs, but only OEM RPs are not confident enough
     * @param database
     * @return 
     */
    public int[] oemPatternPosEstimate(SequenceDatabase database){
        int[] info = new int[]{-1, -1, 0};
        for (int i = 0; i < superitems.size(); i++){
            SuperItem si = superitems.get(i).getSuperItem(database);
            if (si.getType().equals("ARP_OEM") && si.getWeightRatio() >= 0.2){
                info[2] = si.getWeight(); 
                int[] bp = oemBpSearch(i, si.getOri(), database);
                if (bp[0]!=-1 && bp[1]!=-1){
                    info[0] = bp[0];
                    info[1] = bp[1];
                }
                
            }
        }
        return info;
    }
    
    private int[] oemBpSearch(int oemIdx, String oemOri, SequenceDatabase database){
        int[] distToIdx;        
        int[] bp; 
        
        if(oemOri.equals("+")){            
            distToIdx = new int[superitems.size() - oemIdx];
            for (int i = oemIdx; i < superitems.size(); i++){
                SuperItem si = superitems.get(i).getSuperItem(database);                
                if (si.isClipped()){                     
                    distToIdx[i - oemIdx] = si.getPos() ;
//                    lastPos = si.getPos();                                       
                }else{
                    distToIdx[i - oemIdx] = -1;
                }
            }     
            
            bp = oemBpHelper(distToIdx, oemIdx, oemOri, database);
                        
        }else{
            distToIdx = new int[oemIdx + 1];
            for (int i = oemIdx; i >= 0; i--){
                SuperItem si = superitems.get(i).getSuperItem(database);                
                if (si.isClipped()){                     
                    distToIdx[i] = si.getPos();
//                    lastPos = si.getPos();                                       
                }else{
                    distToIdx[i] = -1;
                }
            }
            bp = oemBpHelper(distToIdx, oemIdx, oemOri, database);
            
        }
                                       
        return bp;
        
    }
    
    
    private int[] oemBpHelper(int[] distArray, int oemIdx, String oemOri, SequenceDatabase database){
        int bp[] = new int[]{-1, -1};
        if (oemOri.equals("+")){
            // ARP_OEM is the last superitem in the pattern
            if (oemIdx == superitems.size() - 1){
                return bp;
            }
            for (int i = 1; i < distArray.length;i++){
                if (distArray[i] == -1){
                    continue;
                }

                int closeNum = 1;
                int sumDist = 0;
                int leftBoundIdx = -1;

                for (int j = 1; i + j < distArray.length; j++){
                    if (distArray[i + j] == -1){
                        continue;
                    }
                    sumDist += distArray[i + j] - distArray[i];
                    if (sumDist < 100){
                        closeNum += 1;
                        leftBoundIdx = i + j;
                    }
                }
                if (closeNum > 1 && leftBoundIdx != -1){
                    bp[0] = superitems.get(oemIdx).getSuperItem(database).getPos();
//                    bp[0] = superitems.get(i + oemIdx).getSuperItem(database).getPos();
                    bp[1] = superitems.get(leftBoundIdx + oemIdx).getSuperItem(database).getPos();
                    
                }
            }
        }else{
            // ARP_OEM is the first superitem in the pattern
            if (oemIdx == 0){
                return bp;
            }
    
            for (int i = oemIdx - 1; i >= 0; i--){
                if (distArray[i] == -1){
                    continue;
                }

                int closeNum = 1;
                int sumDist = 0;
                int leftBoundIdx = -1;

                for (int j = 1; i - j >= 0; j++){
                    if (distArray[i - j] == -1){
                        continue;
                    }
                    sumDist += distArray[i] - distArray[i - j];
                    if (sumDist < 100){
                        closeNum += 1;
                        leftBoundIdx = i - j;
                    }
                }
                if (closeNum > 1 && leftBoundIdx != -1){
                    bp[0] = superitems.get(leftBoundIdx).getSuperItem(database).getPos();
                    bp[1] = superitems.get(oemIdx).getSuperItem(database).getPos();                    
                }
            }
        }
        
        return bp;
    }
    
    /**
     * For patterns cannot exact BPs
     * Check if it has multiple BP evidence, such as split-align, clipped reads and other signals, which might be a complex pattern
     * @param database
     * @return 
     */
    
    public boolean isComplexPattern(SequenceDatabase database){
        
        Map<String, Integer> siTypeCount = new HashMap<>();
        boolean hasArp = false;
        int clipNum = 0;
        for (pseudoSuperItem item : superitems){
            SuperItem si = item.getSuperItem(database); 
            String siType = si.getType();
            if (siType.contains("OEM")){
                continue;
            }
            
            if (!si.isARPsuperitem()){
                if (siTypeCount.containsKey(si.getType())){
                    clipNum += 1;
                    int count = siTypeCount.get(si.getType()) + 1;
                    siTypeCount.put(si.getType(), count);
                }else{
                    siTypeCount.put(si.getType(), 1);
                }                
            }else{               
                hasArp = siType.contains("ARP_RR") || siType.contains("ARP_RF") || siType.contains("ARP_FF");
                if (siTypeCount.containsKey(si.getType())){
                    int count = siTypeCount.get(si.getType()) + 1;
                    siTypeCount.put(si.getType(), count);
                }else{
                    siTypeCount.put(si.getType(), 1);
                }
            }
        }
//        double entropy = calPatternEntropy(siTypeCount);
        return hasArp && clipNum > 2;
        
    }
    
    private double calPatternEntropy(Map<String, Integer> map){      
        double ent = 0;
        int sum = 0;
        List<Integer> countList = new ArrayList<>();
        for (Entry<String, Integer> entry : map.entrySet()){
            countList.add(entry.getValue());
            sum += entry.getValue();
        }     
        for(int i = 0; i < countList.size(); i++){
            double prob = (double) countList.get(i) / (double) sum;
            double val = prob * (Math.log(prob)/Math.log(2));
            ent-= val;
        }        
        return ent;
    }
    
}
