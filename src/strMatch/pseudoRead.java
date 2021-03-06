/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package strMatch;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
/**
 *
 * @author jiadonglin
 */
public class pseudoRead {
    String dbId;
    int seqId;
    int firstBaseId;
    int seqSize;
    
    public pseudoRead(String dbId, int seqId, int firstBaseId, int size) {
        this.dbId = dbId;
        this.seqId = seqId;
        this.firstBaseId = firstBaseId;
        this.seqSize = size;
    }
    
    public pseudoRead(pseudoRead pseudo, int baseIdx){
        this.dbId = pseudo.dbId;
        this.seqId = pseudo.seqId;
        this.firstBaseId = pseudo.firstBaseId + baseIdx;
        this.seqSize = pseudo.seqSize - baseIdx;
        
    }
    
    public char getBase(pseudoRead pseudo, Map<String, List<char[]>> database){
        String dbId = pseudo.dbId;
        int baseIdx = pseudo.firstBaseId;
        int seqId = pseudo.seqId;
        char ch = database.get(dbId).get(seqId)[baseIdx];
        return ch;
    }
    
    public char getBaseLocalAlign(pseudoRead pseudo, Map<String, char[]> database){
        String dbId = pseudo.dbId;
        int baseIdx = pseudo.firstBaseId;        
        char ch = database.get(dbId)[baseIdx];
        return ch;
    }
    
    public int indexOf(String item, char[] sequence, int baseIdx){
        int len = sequence.length; 
        if (baseIdx > len){
            return -1;
        }
        String curChar = Character.toString(sequence[baseIdx]);
        if (curChar.equals(item)){
//                itemIdx.add(i);
            return baseIdx;
        }       
        return -1;
    }

}

