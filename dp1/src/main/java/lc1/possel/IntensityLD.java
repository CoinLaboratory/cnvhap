package lc1.possel;

import java.io.File;

import lc1.dp.appl.CNVHap;
import lc1.dp.data.collection.DataCollection;
import lc1.util.Constants;
import lc1.util.Executor;

import org.apache.commons.cli.CommandLine;

public class IntensityLD {
	
	
	/** target is first */
	 public static void main(String[] args){
	        try{
	        	
	        	for(int ii =1; ii<=5; ii++){
	        	  //args = new String[] {"--paramFile", params[ii], "--column", "1:2"};
	        	int type1 = 4; //for reference; 4 is intensity
	        	boolean missing1 = false;
	        	boolean phased = false;
	        	int type2=1; //target
	        	boolean missing2 = false;
	        	  String[] cols  = new String[] {"1","2"};
	        	  CommandLine para = Constants.parse(args, Integer.parseInt(cols[0]), ii, null);
	        	   String dirF = Constants.baseDir();
	               if(dirF.equals(".")){
	               	 dirF = System.getProperty("user.dir");
	               }
	              File dir1 = new File(dirF);
	                 
	               DataCollection target =    CNVHap.read(dir1);
	          	 target.indiv();
	              for(int ij=1; ij<cols.length; ij++){
	               
	            	  Constants.parse(args, Integer.parseInt(cols[ij]), ii, null);
	            	 DataCollection reference = CNVHap.read(dir1); 
	            	 File outF = new File(new File(Constants.outputDir()), "results");
	            	 outF.mkdir();
	            	 DataCollection target1 = (DataCollection)target.clone();
	            	 target1.indiv();
	            	 reference.indiv();
	            
	            	 reference.restricToAlias(target.indiv());
	            	 target1.restricToAlias(reference.indiv());
	            	 File outF1 = new File(outF, Constants.experiment[0]+"_"+Constants.column()+"_"+ii);
	            	  HLALD hlald = new HLALD(outF1,null, reference, target1,type1, type2, missing1, missing2, phased,true,true,0,0);
	            	  hlald.run(false);
	              }
	        	}
	            // BaumWelchTrainer.es.shutdown();
	           //  ProbMultivariate.es.shutdown();
	       
	                                                                                               
	        }catch(Exception exc){
	        	Executor.shutdown();
	        	exc.printStackTrace();
	        }
}
	
}
