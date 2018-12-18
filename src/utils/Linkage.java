/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

import fspm.pseudoSequentialPattern;
import fspm.pseudoSuperItem;
import dataStructure.SequenceDatabase;
import dataStructure.SuperItem;
import htsjdk.samtools.QueryInterval;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;

/**
 *
 * @author jiadonglin
 */
public class Linkage {
    
    int supLink;
    
    public Linkage(){
        
    }
    
    public int getSupLink(){
        return supLink;
    }
    
    public boolean twoSuperItemLinkCheck(SuperItem superItemOne, SuperItem superItemTwo){
        boolean supported = false;
        String[] qnameOneByteList = superItemOne.getQNames();
        String[] qnameTwoByteList = superItemTwo.getQNames();
        Set<String> uniqueQName = new HashSet<>();
        
        int supportedARPs = 0;
        for (String qname : qnameOneByteList){            
            uniqueQName.add(qname);
        }       
        for (String qname : qnameTwoByteList){            
            if (uniqueQName.contains(qname)){
                supportedARPs += 1;
                if (supportedARPs > 1){
                    supported = true;                   
                }
            }
            uniqueQName.add(qname);
        }
        supLink = supportedARPs;
        
        return supported;
    }
    
    public int mateSuperItemSearch(List<pseudoSuperItem> superItems, pseudoSuperItem targetSuperItem, SequenceDatabase database){
        QueryInterval targetInterval = targetSuperItem.getSuperItem(database).getSuperitemMateRegion();
        int length = superItems.size();
        int start = 0;
        int end = length - 1;
        int mateIndex = -1;
        if (length == 1){
            QueryInterval curInterval = superItems.get(0).getSuperItem(database).getSuperitemRegion();
            if (hasOverlap(curInterval, targetInterval)){
                return 0;
            }else{
                return -1;
            }
        }
        
        while (start <= end){
            int mid = start + (end - start) / 2;
            SuperItem curSuperItem = superItems.get(mid).getSuperItem(database);
            
            QueryInterval curInterval = curSuperItem.getSuperitemRegion();
            
            if (hasOverlap(curInterval, targetInterval)){
                               
                mateIndex = mid;
                break;
            }
            if (isAheadInterval(curInterval, targetInterval)){
                end = mid - 1;
            }
            else if (isAfterInterval(curInterval, targetInterval)){
                start = mid + 1;
            }
        }
        return mateIndex;
    }
        
    private boolean hasOverlap(QueryInterval aInterval, QueryInterval targetInterval){       
        return aInterval.overlaps(targetInterval);
    }
    private boolean isAheadInterval(QueryInterval aInterval, QueryInterval targetInterval){        
        return aInterval.start > targetInterval.end;
    }

    private boolean isAfterInterval(QueryInterval aInterval, QueryInterval targetInterval){
        return aInterval.end < targetInterval.start;
    }
    
    public List<pseudoSequentialPattern> minusItemAndCopy(List<pseudoSequentialPattern> patterns, int index){
        List<pseudoSequentialPattern> removedPatterns = new ArrayList<>();
        int length = patterns.size();
        for (int i = 0; i < length ; i ++){
            if (i != index){
                removedPatterns.add(patterns.get(i));                
            }
        }
        return removedPatterns;
    }
    public List<pseudoSequentialPattern> minusItemAndCopyWithIndexMap(List<pseudoSequentialPattern> patterns, int index, Map<Integer, Integer> indexMap){
        List<pseudoSequentialPattern> removedPatterns = new ArrayList<>();
        int length = patterns.size();
        for (int i = 0; i < length ; i ++){
            pseudoSequentialPattern pattern = patterns.get(i);
            if (i != index && pattern.hasArpSuperItem()){
                removedPatterns.add(pattern);  
                indexMap.put(removedPatterns.size() - 1, i);
            }
        }
        return removedPatterns;
    }
}
