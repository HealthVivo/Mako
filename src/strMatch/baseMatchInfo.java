/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package strMatch;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author jiadonglin
 */
public class baseMatchInfo {
        String  strIdentity;
        String curMatchStr;
//        List<Integer> matchedStartAt = new ArrayList<>();
        int matchedStartAt = -1;
        int baseIdx;
        char base;
        boolean isForward;
        boolean isBoth;
        List<stringIdentity> baseStringIdentitys = new ArrayList<>();
        
        public baseMatchInfo(String base){
            this.curMatchStr = base;
        }
        
        public baseMatchInfo(String matchStr, int matchedAt){
            this.curMatchStr = matchStr;
            this.matchedStartAt = matchedAt;
        }
//        public baseMatchInfo(String strID, String matchStr, List<Integer> matchedAt){
//            this.curMatchStr = matchStr;
//            this.matchedStartAt.addAll(matchedAt);
//            this.strIdentity = strID;
//        }
        public baseMatchInfo(char base, boolean isForward, boolean isBoth){
            this.base = base;
            this.isForward = isForward;
            this.isBoth = isBoth;
        }
        
        public void setMatchedAt(int idx){
            matchedStartAt = idx;
        }
        public void setBase(char base){
            this.base = base;
        }

        public void setIsBoth(boolean isBoth){
            this.isBoth = isBoth;
        }
        public void setIsForward(boolean isForward){
            this.isForward = isForward;
        }
        public void setBaseStringIDs(List<stringIdentity> identitys){
            baseStringIdentitys = identitys;
        }
        public baseMatchInfo append(baseMatchInfo info){
            
            matchedStartAt = info.matchedStartAt;
            String newStr = curMatchStr + info.curMatchStr;
            return new baseMatchInfo(newStr, matchedStartAt);
        }
        public int curMatchedStrLen(){
            return curMatchStr.length();
        }
}
