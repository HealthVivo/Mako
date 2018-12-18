/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mako;
import dataStructure.SequenceDatabase;
import fspm.CFSPM;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import readers.fileReader;
import readers.SignalReader;


/**
 *
 * @author jiadonglin
 */
public class Mako {

      /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException{
    // TODO code application logic here

//        terminalParamLoader paramLoader = new terminalParamLoader();
//        paramLoader.loadParam(args);
//        if (paramLoader.hasParamInput){
//            if (!paramLoader.siFileMode){
//                System.out.println("fragMean: " + paramLoader.fragMean + " fragStd: " + paramLoader.fragStd 
//                + " readLen: " + paramLoader.readLen + " clustD: " + paramLoader.clusteringDist);
//                System.out.println("Working directory: " + paramLoader.workDir);
//                
//                SignalReader myReader = new SignalReader(paramLoader.fragMean, paramLoader.fragStd, paramLoader.cutStd, 
//                        paramLoader.readLen, paramLoader.clusteringDist, paramLoader.minMapQ);
//                
//                myReader.doWork(paramLoader.bamFile, paramLoader.fastaIndexFile, paramLoader.chr, paramLoader.givenRegionS, paramLoader.givenRegionE, paramLoader.superitemOut, paramLoader.abnormalSignalOut); 
//
//                
//                SequenceDatabase sequenceDatabase = new SequenceDatabase(); 
//                int minSup = 50;
//                
//                System.out.println("Super-Item generation completed!!\n\nSV path: " + paramLoader.svOut);
//                
//                BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(paramLoader.svOut));
//
//
//                sequenceDatabase.loadSequencesFromFile(paramLoader.superitemOut, svRegionWriter);
//
//                CFSPM algoContiguousFSPM = new CFSPM(minSup, paramLoader.patternMaxSpan);
//                algoContiguousFSPM.setParams(myReader.chromNameMap, paramLoader.regionMaskFile);
//                algoContiguousFSPM.runAlgorithm(sequenceDatabase, paramLoader.frequentPatternOut, paramLoader.mergedPatternOut, paramLoader.susRegionWriter, svRegionWriter, paramLoader.fastaFile);
//                algoContiguousFSPM.printAlgoStatistics();
//            }
//            else{
//                SequenceDatabase sequenceDatabase = new SequenceDatabase(); 
//                int minSup = 50;
//
//                BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(paramLoader.svOut));
//
//
//                sequenceDatabase.loadSequencesFromFile(paramLoader.superitemOut, svRegionWriter);
//
//                CFSPM algoContiguousFSPM = new CFSPM(minSup, paramLoader.fragMean);
//                fileReader reader = new fileReader();                
//                algoContiguousFSPM.setParams(reader.getIdxToChrName(paramLoader.fastaIndexFile), paramLoader.regionMaskFile);
//                algoContiguousFSPM.runAlgorithm(sequenceDatabase, paramLoader.frequentPatternOut, paramLoader.mergedPatternOut, paramLoader.susRegionWriter, svRegionWriter, paramLoader.fastaFile);
//                algoContiguousFSPM.printAlgoStatistics();
//            } 
//        }

            String workDir = "/Users/jiadonglin/SV_data/HG00514/";
            String bamFile = workDir + "5X/HG00514.5X.bam";

//        String workDir = "/Users/jiadonglin/SV_data/synthetics/svelter_simple/";
//        String bamFile = workDir + "simHomo.30X.sorted.rg.bam";
//            
            String superItemOut = workDir + "mako/5X/wgs.all.superitems.txt";
            String svOutPath = workDir + "mako/5X/Mako.vcf";
//            String susComplexRegion = workDir + "mako/CSVsRegion/susCSVs.region.vcf";

//        String superItemOut = workDir + "wgs.all.superitems.txt";
//        String svOutPath = workDir + "Mako.vcf";


            String fastaIndexFile = "/Users/jiadonglin/SV_data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai";
            String fastaFile = "/Users/jiadonglin/SV_data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";


//        String fastaIndexFile = "/Users/jiadonglin/SV_data/ref_genome/human_g1k_v37.fasta.fai";
//        String fastaFile = "/Users/jiadonglin/SV_data/ref_genome/human_g1k_v37.fasta";
//            
//            String excludableRegion = "/Users/jiadonglin/SV_data/ref_genome/grch38.exclude.bed";

        int fragMean = 566;
        int fragStd = 160;
        int readLen = 126;

//        int fragMean = 500;
//        int fragStd = 50;
//        int readLen = 150;
        
        int minMapQ = 10;
        int cutStd = 3;
        int ArpClusterDist = fragMean - readLen;     
        int patternGrowthLimit = fragMean;


        SignalReader myReader = new SignalReader(fragMean, fragStd, cutStd, readLen, ArpClusterDist, minMapQ);
        myReader.doWork(bamFile, fastaIndexFile, null, 0, 0, superItemOut, null);                        
        
        SequenceDatabase sequenceDatabase = new SequenceDatabase(); 

        int minSup = 50;

        BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(svOutPath));
//            BufferedWriter susRegionWriter = new BufferedWriter(new FileWriter(susComplexRegion));

        sequenceDatabase.loadSequencesFromFile(superItemOut, svRegionWriter);

        CFSPM algoContiguousFSPM = new CFSPM(minSup, patternGrowthLimit);
        fileReader reader = new fileReader();                
        algoContiguousFSPM.setParams(reader.getIdxToChrName(fastaIndexFile), null);
        algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, null, null, svRegionWriter, fastaFile);
        algoContiguousFSPM.printAlgoStatistics(); 
    }
   
}
