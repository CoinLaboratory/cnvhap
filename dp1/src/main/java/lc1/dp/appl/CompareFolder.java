package lc1.dp.appl;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.core.Sampler;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.util.Constants;

import org.apache.commons.cli.CommandLine;

public class CompareFolder {

    public static void main(String[] args){
        try{
        	 if(args.length==0){
                 args = new String[] {"--paramFile", "param_compare.txt", "--column", "1"};
             }
            String cols = Constants.getCols(args)[0];
            
             CommandLine para = Constants.parse(args, Integer.parseInt(cols),1, null);
             try{
             Thread.currentThread().wait(1000);
             }catch(Exception exc){
            	 
             }
          String[] inp =  para.getOptionValues("inputDir");
          System.err.println(Arrays.asList(inp));
          String chr = Constants.chrom0();
          final int[][] mid = Constants.mid();
          final int[] kb = Constants.restrictKb(0);
          for(int i=0; i<mid.length; i++){
        	  
          
          mid[i][0] = Math.max(0, mid[i][0] - kb[0]);
          mid[i][1] = Math.min(Integer.MAX_VALUE, mid[i][1]+kb[1]);
          }
            for(int i=1; i<inp.length; i++){
                File f1 = new File(inp[0], chr+".zip");
                File f2 = new File(inp[i], chr+".zip");
		        DataCollection original =  null;//new SimpleDataCollection(f1,(short)0, 2, mid, null);
		        DataCollection inferred =  null;//new SimpleDataCollection(f2,(short)0, 2, mid, null);
		       if(i==2){ Sampler sampler1 = new Sampler(inferred, new File(Constants.outputDir(), "samples"), false);
		       Integer[] is =new Integer[Constants.numRep];
		       for(int i1=0; i1<is.length; i1++){
		    	   is[i1] = i1;
		       }
		        sampler1.calcBestPathSampling(Arrays.asList(is));
		       }
		        if(!f1.exists())throw new RuntimeException("!!");
		        if(!f2.exists()) throw new RuntimeException("!!");
        //dc1.mix();
        ComparePhase cp = new ComparePhase( inferred, original, "inferred", "original");
        PrintWriter summ = new PrintWriter(new FileWriter(new File(inp[i].replace('/', '-')+"_phase_comp.txt")));
        PrintWriter pw = new PrintWriter(new FileWriter(new File(inp[i].replace('/', '-')+"_summary.txt")));
        PrintWriter[] logSumm = new PrintWriter[3] ;
        for(int k=0; k<logSumm.length; k++){
            logSumm[k] = new PrintWriter(new FileWriter(new File(inp[i]+"_"+k+"_summary.txt")));
        }
      //  cp.compare(0.0, pw, logSumm);
        for(int k=0; k<logSumm.length; k++){
            logSumm[k].close();
            
        }
        pw.close();
        File logDir  = new File(inp[i].replace('/', '-')+"_comparison");
        logDir.mkdir();
        cp.comparePhase(summ, logDir, original.loc);
        summ.close();
      }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
}
