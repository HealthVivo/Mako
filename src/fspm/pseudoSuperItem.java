/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fspm;

import dataStructure.SequenceDatabase;
import dataStructure.SuperItem;

/**
 *
 * @author jiadonglin
 */
public class pseudoSuperItem implements Comparable<pseudoSuperItem>{
    int sequenceId;
    protected int itemsetIdx;
    protected int itemIdx;   
    protected int leftMostPos; // start of the superitem interval. record position not mutational position.

    public pseudoSuperItem(int seqId, int itemsetIndex, int itemIndex){
        sequenceId = seqId;
        itemsetIdx = itemsetIndex;
        itemIdx = itemIndex;
    }
    public SuperItem getSuperItem(SequenceDatabase database){
        SuperItem si = database.getSequenceByID(sequenceId).superItemAtPos(itemsetIdx, itemIdx);
        return si;
    }
    public void setPsSuperitemLeftPos(SequenceDatabase database){
        leftMostPos = getSuperItem(database).getSuperitemRegion().start;
    }
    public int getLeftMostPos(){
        return leftMostPos;
    }
    @Override
    public int compareTo(pseudoSuperItem otherItem){
        return leftMostPos - otherItem.leftMostPos;
    }  
}
