    /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package channel;

import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import readers.CigarOps;
import readers.myCigarOp;

/**
 *
 * @author jiadonglin
 */


public class MutSignal implements Comparable<MutSignal>{
    
//    protected String queryName;
    protected String queryName;
    protected String mutSignalType;
    protected int mutSignalPos;
    protected String mutSignalOri;
    protected boolean isizeNormal = true;

    protected int signalRefIndex;
    protected String signalRefName;
    protected int insertSize;
    protected int recordPos;
    protected int recordMatePos;
    protected String cigarString;
    protected int mapQ;
    protected String signalMateRefName;
        
    protected boolean isARP = false;
    protected boolean isSplitAlign = false;
    protected boolean isInterChrom = false;
    protected int splitAlignPos;
    protected String splitAlignChr;
    protected String readSequence;
    protected int[] clippedStatus;
    protected int longestD;
    
    public MutSignal(){
        
    }

    @Override
    public int compareTo(MutSignal otherMutSignal){
        return mutSignalPos - otherMutSignal.getMutPos();
    }
    
    public MutSignal(SAMRecord record, String signalType, int pos, String ori){
        
        
        queryName = record.getReadName();        
        signalRefIndex = record.getReferenceIndex();
        signalRefName = record.getReferenceName();
        signalMateRefName = record.getMateReferenceName();                       
        insertSize = Math.abs(record.getInferredInsertSize());
        mutSignalType = signalType;
        mutSignalPos = pos;
        mutSignalOri = ori;
        recordPos = record.getAlignmentStart();
        recordMatePos = record.getMateAlignmentStart();
        cigarString = record.getCigarString();
        readSequence = record.getReadString();
        mapQ = record.getMappingQuality();
    
        CigarOps cigarOps = new CigarOps(record);
        clippedStatus = cigarOps.getClippedStatus();
        longestD = cigarOps.getLongestD(record.getCigar().getCigarElements());
        
        if (mutSignalType.contains("ARP")){
            isARP = true;
        }   
        if (record.isSecondaryAlignment()){           
            splitAlignPos = decodeSplitAlign(record.getAttribute("SA").toString(), mapQ);            
        }
        if (!signalRefName.equals(signalMateRefName)){
            isInterChrom = true;
        } 
    }
    @Override
    public boolean equals(Object obj){
        if (obj instanceof MutSignal){
            MutSignal mutSignal = (MutSignal) obj;
            return (mutSignal.getMutSignalType().equals(this.mutSignalType));
        }else{
            return false;
        }
    }
    @Override
    public int hashCode(){
        return mutSignalType.hashCode();
    }
    
    
    public boolean withinDistance(MutSignal otherMutSignal, int thresh){
        int disDiff = Math.abs(otherMutSignal.getMutPos() - mutSignalPos);
        int typePenalty = 0;
        if (!mutSignalType.equals(otherMutSignal.getMutSignalType())){
            typePenalty = 1;
        }
        return disDiff <= thresh / (1 + typePenalty);
    }
    
    @Override
    public String toString(){     
        StringBuilder sb = new StringBuilder();
        sb.append(queryName);
        sb.append("\t");
        sb.append(mutSignalType);
        sb.append("\t");
        sb.append(mutSignalPos);
        sb.append("\t");
        sb.append(recordPos);
        sb.append("\t");
        sb.append(signalRefName);
        sb.append("\t");
        sb.append(insertSize);
        sb.append("\t");
        sb.append(cigarString);
        return sb.toString();
    }
    public int getLongestD(){
        return longestD;
    }    
    public int[] getClippedStatus(){
        return clippedStatus;
    }
    public String getqName(){
        return new String(queryName);
    }
    public String getMutSignalType(){
        return mutSignalType;
    }
    public String getMutSignalOri(){
        return mutSignalOri;
    }
    public int getSignalChromIdx(){
        return signalRefIndex;
    }
    public String getSignalRef(){
        return signalRefName;
    }
    public int getRecordPos(){
        return recordPos;
    }
    public int getMateRecordPos(){
        return recordMatePos;
    }
    public boolean isARPSignal(){
        return isARP;
    }
    public boolean isIsizeNormal(){
        return isizeNormal;
    }
    public boolean isSplitAlign(){
        return isSplitAlign;
    }
    public boolean isInterChrom(){
        return isInterChrom;
    }
    public int getMutPos(){
        return mutSignalPos;
    }
    public int getMapQ(){
        return mapQ;
    }
//    public String getSplitAlignInfo(){
//        return splitAlignInfo;
//    }
    public int getSplitAlignPos(){
        return splitAlignPos;
    }
    public String getSplitAlignChr(){
        return splitAlignChr;
    }
    public String getReadSequence(){
        return readSequence;
    }
    public void setIsizeNormal(int isizeUpper, int isizeLower){
        if (this.insertSize > isizeUpper || this.insertSize < isizeLower){
            isizeNormal = false;
        }
    }
    
    private int decodeSplitAlign(String splitAlignString, int mapQ){
        
        String[] tokens = splitAlignString.split(",");        
        splitAlignChr = tokens[0];
        int pos = Integer.parseInt(tokens[1]);
        
        if (splitAlignChr.equals(signalRefName) && mapQ >= 20){
            isSplitAlign = true;
            List<myCigarOp> cigarOps = getSplitAlignPosFromCigarString(tokens[3]);                    
            myCigarOp firstOp = cigarOps.get(0);
            if (firstOp.getOp().equals("M")){
                pos += firstOp.getOpLength();
            }      
        }else{
            pos = -1;
        }          
        return pos;
    }
    public List<myCigarOp> getSplitAlignPosFromCigarString(String cigar){
       List<myCigarOp> cigarOps = new ArrayList<>();
       Pattern cigarPattern = Pattern.compile("[0-9]+[MIDNSHP]");
       Matcher cigarOpStrings = cigarPattern.matcher(cigar);

       while (cigarOpStrings.find()){
           String opString = cigarOpStrings.group();
           int strLen = opString.length();
           String op = opString.substring(strLen - 1, strLen);
           int opLength = Integer.parseInt(opString.substring(0, strLen - 1));

           myCigarOp cigarop = new myCigarOp(op, opLength);
           cigarOps.add(cigarop);

       }

       return cigarOps;
   }       
}

