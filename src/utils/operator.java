/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

/**
 *
 * @author jiadonglin
 */
public class operator {
    
    public String bamFile;
    public String fastaIndexFile;
    public String superitemOutFile;
    public String signalOutFile;
    public int readLen = -1;
    public int fragMean = -1;
    public int fragStd = -1;
    
    
    public void readConfigFile(){
        
    } 
    public String getHelpInfo(){
        StringBuilder sb = new StringBuilder();
        String bamConfigInfo = "config: BAM configuration file\n";
        String cutStdInfo = "cutStd: Fragment standared deviation cutoff (default 3)\n";
        String maxDistInfo = "maxD: Discordant read-pair clustering maximum distance (default mean-2*readlen)\n";
        String minQInfo = "minQ: Minimum mapping quality filter low quality alignment (default 20)\n";
        String chrInfo = "chrom: User specified chromsome to inspect (default None)\n";
        String regionSInfo = "start: User specified genome region start pos (default None)\n";
        String regionEInfo = "end: User sepcified genome region end pos (default None)\n";
        String bamFileInfo = "bamFile: The input BAM file to process\n";
        String fastaIndexInfo = "faiFile: The reference index file to use\n";
        String superitemOutInfo = "siOut: SuperItem output path\n";
        String signalOutInfo = "sigOut: Abnomral signals output path (optional)\n";
        
        return sb.toString();
    }
}
