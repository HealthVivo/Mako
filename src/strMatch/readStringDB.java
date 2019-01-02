/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package strMatch;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

/**
 *
 * @author jiadonglin
 */
public class readStringDB {
    List<char[]> sequences = new ArrayList<>();
    String dbIdentity;
    int size;
    Map<Character, List<baseIdentity>> singleBaseTracker = new HashMap<>();
    List<pseudoRead> pseduoDB = new ArrayList<>();
    
    public readStringDB(List<String> strings, String dbName){
        dbIdentity = dbName;
        size = strings.size();
        for (int i = 0; i < size; i++){
            String str = strings.get(i);
            char[] charArray = str.toCharArray();
            int len = charArray.length;
            pseduoDB.add(new pseudoRead(dbIdentity, i, 0, len));
            for (int j = 0; j < len; j++){
                char ch = charArray[j];
                baseIdentity id = new baseIdentity(dbIdentity, i, j);
                List<baseIdentity> baseAppears = singleBaseTracker.get(ch);
                if (baseAppears == null){
                    baseAppears = new ArrayList<>();
                    baseAppears.add(id);                    
                }else{
                    baseAppears.add(id);
                }
            }
            sequences.add(charArray);            
        }
    }
    public char[] getSequenceAt(int idx){
        return sequences.get(idx);
    }       
}
