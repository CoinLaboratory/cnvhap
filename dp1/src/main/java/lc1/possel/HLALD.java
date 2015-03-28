package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.CGH.AberationFinder;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.HLADataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.states.EmissionState;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix1DProcedure;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Statistic;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.random.StudentT;
import edu.mit.wi.haploview.HaploView;

public class HLALD {
	static double rThresh =Constants.ld_r2_thresh();
	static double pThresh = 1.0;
	static double distThresh = Constants.ld_dist_thresh();
	static double printVals = 0.8;
	public static void main1(String[] args){
		try{
			Algebra a = new Algebra();
			 DoubleMatrix2D matrix = a.transpose(new DenseDoubleMatrix2D(new double[][] {{0,1,2,3,4,5},{3,2,6,4,7,4}}));
				DoubleMatrix2D matrix1 = Statistic.covariance(matrix);
		      	Statistic.correlation(matrix1);
		      	System.err.println(matrix1.get(0, 1));
		      	System.err.println(matrix.rows());
		      	StudentT t = new StudentT(matrix.rows()-2, null);
		      	System.err.println(t.cdf(statistic(matrix1.get(0,1), matrix.rows())));
		      	
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static double  statistic(double r, double N){
		return r*Math.sqrt((N-2)/(1-Math.pow(r, 2)));
	}
	 public static void main(String[] args1){
	        try{
	        	  File user = new File(System.getProperty("user.dir"));
	        	  File[] args = user.listFiles(new FileFilter(){

					public boolean accept(File pathname) {
						return pathname.isDirectory() && pathname.getName().startsWith("gg_1");
					}
	        	  });
	        	for(int ik=0; ik<args.length; ik++){
	        //	Constants.topBottom = false;
	        	 
	        	SimpleDataCollection sdt = 
	                SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader( args[ik],"phased2.txt"),Emiss.class, Emiss.getEmissionStateSpace(1));
	        	sdt.calcIndiv();
	        	HLADataCollection hla = new HLADataCollection(new File(user, "6.zip"),(short)0, 2, 
	        			new int[][] {new int[] {0,Integer.MAX_VALUE-1}}, 
	            		new File(user, "build35.txt"), null);
	        	hla.restricToAlias(sdt.indiv());
	        	
	        //	sdt.removeKeyIfStartsWith("NA");
	          //  sdt.calculateMaf(false);
	         
	            int len = sdt.length();
	           
	            File snpFile  = new File(user, "snp.txt");
	            if(snpFile.exists()){
	                sdt.snpid= new ArrayList<String>();
	                sdt.readPosInfo(snpFile, new int[] {0}, true, new List[] {sdt.snpid}, new Class[] {String.class});
	            }
	            File out = new File(user, "LD_results");
	            out.mkdir();
	            File outF = new File(out,args[ik].getName());
	            
	            HLALD hlald = new HLALD(outF,null, sdt, hla,3,3,false, false, false,true,true,0,0);
	          hlald.run(false);
	        	}
	           // writeEHH(sdt, new File(user, "sweep"),Constants.chrom0()+"", null, "");
	                                                                                               
	        }catch(Exception exc){
	        	exc.printStackTrace();
	        }
}
	
	 DataCollection dc_1; //source
	 DoubleMatrix2D matrix;
	 
	 //static int base = 65;
	 final PrintWriter pw;
	 final DataCollection dc_2;  //target
	
	 
	  class Column{
		 DataCollection dc;
		 int type;
		 int currentPos=-1;
		 List<Comparable> geno;
		// List<Comparable> poss;
		 boolean missing1;
	//	 Map<Double,Double> replace = null;
		 
		 public Column(DataCollection dc_1, int type1, boolean missing12) {
			this.dc = dc_1;
			this.type = type1;
			
			this.missing1 = missing12;
		}

		public void set(int i, boolean avg){
			if(currentPos==i) return;
			 this.currentPos = i;
			
			 if(!phased  || type==-4 || type==-5){
				 if(type==-4 || type==-5) geno =  dc.getIntensity(i, type==-4);
				 else if(type==0 || type==1 || type==2){ //type=1 is no A, type==2 is noB
					 geno = dc.getValues(i, type,avg);
				 }
				 else if (type==3){
					 geno = dc.getValues(i,0,avg);
					 List<Comparable> geno1 = dc.getValues(i,2,avg);
					 for(int ik=0; ik<geno.size(); ik++){
						 double cn = ((Double)geno.get(ik));
						 double v =  Double.isNaN(cn) || Math.abs(cn - 2.0)>0.01 ? Double.NaN : 
							 ((Double)geno1.get(ik));
						 geno.set(ik, 
								v);
					 }
					 
				 }
			 }
			 else{
				 geno = dc.getHaplotypes1(i);
			 }
			/* if(type==3){
				 Set<Comparable> s = new HashSet<Comparable>();
				 for(int j=0; j<geno.size(); j++){
					 ComparableArray comp = (ComparableArray)geno.get(j);
					 for(int k=0; k<comp.size(); k++){
						 s.add(comp.get(k));
					 }
				 }
				 poss = new ArrayList<Comparable>(s);
			 }
			 else{
				 poss = baseList;
			 }*/
		 }

		public Double get(int ik) {
			Comparable comp = geno.get(ik);
			Double res =  comp instanceof Double ? (Double)comp : noCopies((ComparableArray) geno.get(ik), type, missing1);
		/*	if(replace!=null){
				Double r1 = replace.get(res);
				if(r1!=null) return r1;
			}*/
			return res;
		}
		 
	 }
	 StudentT t;
	 int N;
	  boolean phased;
	 //type 0==no_copies; type 1 noA; type 2 noB ;type==3 all states type ==4 intensity type 5 = BAF
	
	 
	
	Column column1;
	Column column2;
	 
	 
	final File outDir; 
	public HLALD(File outDir, PrintWriter pw, DataCollection sdt,DataCollection hla, int type1, int type2, boolean missing1, boolean missing2, boolean phased,
			boolean toP,boolean avg, int index1, int index2) throws Exception{
		// this.geno2 =    new List[hla.loc.size()]; 
		this.phased = phased;
		 this.dc_1 = sdt;
		 this.avg = avg;
		 this.column1 = new Column(dc_1, type1, missing1);
		 this.outDir = outDir;
		 this.pw = pw;
		
		// this.phased = phased;
		 if(type1==-4 || type1==-5 || type2 ==-5 || type2==-4) this.phased = false;
		
//		 f.mkdir();
		 this.dc_2 = hla;
		 this.column2 = new Column(dc_2, type2, missing2);
		// this.outF = f;
		
	
		//   	pw.println(dc_1.loc.get(i1)+" "+dc_2.loc.get(i2)+
	      		//	(dc_1.snpid.size()>0 ? dc_1.snpid.get(i1):"")+" "+
	      			//dc_2.snpid.get(i2)+" "+r2+" "+p);
		// poss_2 = new List[geno2.length];
		 ty1 = getType(type1);
		 ty2 = getType(type2);
		
		// column1.set1(0);
		column1.set(index1, avg);
		 N = column1.geno.size();
		 matrix = new DenseDoubleMatrix2D(N, 2);
		 t = new StudentT(N-2, null);
	 }
	 
	 //type 0==no_copies; type 1 noA; type 2 noB ;type==3 all states type ==4 intensity type 5 = BAF
	public static String getType(int type22) {
		if(type22==0) return "noCopies";
		else if(type22==3) return "SNP allele";
		else if(type22==-4) return "LRR";
		else if(type22==-5) return "BAF";
		else if(type22 < Emiss.ems.length) return ""+Emiss.ems[type22];
		else return type22+"";
		
	}

	public List<Double> run(boolean keep) throws Exception{
		List<Double> l = new ArrayList<Double>();
		 for(int i=0; i<dc_1.loc.size(); i++){
			// this.set2(i);
			 int pos1 = dc_1.loc.get(i);
			 for(int j=0; j<dc_2.loc.size(); j++){
				 if(j==i && this.ty1.equals(this.ty2)) continue;
				 int pos2 = dc_2.loc.get(j);
				 if(Math.abs(pos2 - pos1) < distThresh ){
					if(keep) l.add( calcLD(i, j));
					else{
						calcLD(i,j);
					}
					 pw.flush();
				 }
			 }
		 }
		// pw.close();
		 return l;
	 }
	
	 int current2Pos =0;
	//final File outF;
	final String ty1, ty2;
	List<Comparable> baseList = Arrays.asList(new Comparable[] {-1});
	final boolean avg;
	/*2 is the 'target' */
	 public Double calcLD(int i1, int i2) throws Exception{
		 	column1.set(i1, avg);
		 	column2.set(i2, avg);
		 	Double R = null;
	      	List<Comparable> em1,genos1; 
	     // 	if(dc_1.snpid.get(i1).startsWith("cnv")){
	      //		Logger.global.info("h");
	     // 	}
	      
	      	 int cntNull=0;
	      	 for(int ik=0; ik<N; ik++){
			      
				      		Double c2 = column2.get(ik);
				      			
				      		Double c1 = column1.get(ik);
				      	//	if(c1!=2) num++;
				      		matrix.setQuick( ik,1, c2==null ? Double.NaN : c2);
				      		matrix.setQuick(ik,0,c1==null ? Double.NaN: c1);
				      	  	if(c1==null || c2==null || c1==null){
					      		cntNull++;
					      	}
				      	}
				      //	if(allEqual(matrix.)){
				      	DoubleMatrix2D matrix2 =  submatrix(matrix);
					      	DoubleMatrix2D matrix1 = Statistic.covariance(matrix2);
					      
					      	 if(zeroMat(matrix1)){
					    		 R = Double.NaN;
					    	 }
					    	 else {
					      	Statistic.correlation(matrix1);
					      	if(cntNull>0){
					      		System.err.println(matrix1.get(1, 0));
					      	}
					    	if(true || matrix1.getQuick(0, 1)>=0 ){
					    		double r = matrix1.getQuick(0, 1);
					    		
						      
						      double r2 = Math.pow(r, 2);
						      if(r2>1.0-1e-3){
						    	//  System.err.println("r2 out of bounds "+r2);
						    	  r2 = 1.0;
						      }
						    /*  if(r2<1.0){
						    	  System.err.println("h");
						      }*/
						      R = r2;
						    /* if(Double.isNaN(r) || Double.isInfinite(r) ){
						 	DoubleMatrix2D matx = Statistic.covariance(matrix2);
						  	 System.err.println("is nan");
						      }*/
						     /* else */
						     if(r2>rThresh){
						    	 double p = 2*(1-t.cdf(Math.abs(statistic(r, matrix2.rows()))));
						     if ( p<pThresh && pw!=null){
						    		int loc1 = dc_1.loc.get(i1);
						    		int loc2 = dc_2.loc.get(i2);
						      	pw.println(ty1+"\t"+ty2+"\t"+loc1+"\t"+loc2+"\t"
						      			+Math.abs(loc2 - loc1)	+"\t"+
						      			(dc_1.snpid.size()>0 ? dc_1.snpid.get(i1):"")+"\t"+
						      			dc_2.snpid.get(i2)+"\t"+String.format("%5.3g",r2)+"\t"+String.format("%5.3g", p));
						      	if(//dc_2.snpid.get(i2).equals(
						      		//	"rs1"
						      			//"A_14_P137822"
						      			r2>printVals
//						      			dc_1.snpid.get(i1).equals("rs10828906") && dc_2.snpid.get(i1).equals("rs11015065") || 
	//					      			dc_2.snpid.get(i1).equals("rs10828906") && dc_1.snpid.get(i1).equals("rs11015065")
						      			) {
						      		
						      		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(this.outDir, 
						      				Math.round(r2*100)+"_"+dc_1.snpid.get(i1)+"_"+dc_2.snpid.get(i2)+"_"+
						      				String.format("%5.3g", p).replaceAll("\\s+", "")))));
						      	//	pw.println(em.toString().replaceAll(Matcher.quoteReplacement(ch1.toString()), "***"));
						      	//	pw.println(str.toString().replaceAll(Matcher.quoteReplacement(ch.toString()), hla.conversion[j].get(ch).toString()));
						      		pw.println("id\t"+ty1+"_"+dc_1.snpid.get(i1)+"_"+dc_1.loc.get(i1)+"\t"+ty2+"_"+dc_2.snpid.get(i2)+"_"+dc_2.loc.get(i2));
						      		for(int kk=0; kk<N; kk++){
						    			String id =dc_2.indiv().get(phased ? (int) Math.floor((double)kk/2.0) : kk); 
						    			pw.print(id);
						    			for(int j1=0; j1<2; j1++){
						    				pw.print("\t"+matrix.getQuick(kk, j1));
						    			}
						    			pw.println();
						    		}
						      	pw.close();
						      	}
						      	
					    	}
						     }
						     }
			      //	}
	      	 }
					    	if(Double.isNaN(R)) {
							//	R = 0.0;
							}
					    	return R;
		 }
	
	 private boolean zeroMat(DoubleMatrix2D matrix1) {
		for(int i=0; i<matrix1.rows(); i++){
			for(int j=0; j<matrix.columns(); j++){
				if(matrix.get(i, j)>1e-8) return false;
			}
		}
		return true;
	}

	private DoubleMatrix2D submatrix(DoubleMatrix2D matrix2) {
		return matrix2.viewSelection(new DoubleMatrix1DProcedure(){

			public boolean apply(DoubleMatrix1D arg0) {
			   for(int i=0; i<arg0.size(); i++){
				   if(Double.isNaN(arg0.get(i))) return false;
			   }
			   return true;
			}
			
		});
	}

	private Double noCopies(ComparableArray array,  int type, boolean missing) {
		if(array==null) return Double.NaN;
		if(missing){
			for(int i=0; i<array.size(); i++){
				Comparable comp = array.get(i);
				if(comp==Emiss.N){
					return Double.NaN;
				}
			}
		}
		 if(type==0) return (double)  array.noCopies(true);
		 else if (type==6){
			 return array.noCopies(true)==array.size() ? array.noB(2) : Double.NaN; 
		 }
		 else if(type>0) return (double) array.noB(type);
		
		throw new RuntimeException("!!");/* int cnt =0;
			for(int i=0; i<array.size(); i++){
				Comparable string =  array.get(i);
				
				if(string.equals(c)) cnt++;
			}
			return (double) cnt;*/
	}
	

	public static int writeEHH(
		        
		        DataCollection sdt2, File f, String chr, Set<String> indiv, String name) throws Exception{
		        if(!f.exists()) {
		            f.mkdir();
		        }
		      
		        EmissionState maf1 = null;//sdt2.calculateMaf1();
		    //    sdt2.maf = maf1;
		        Collection<Integer> toD = maf1.getConstantPos();
		     // toD.addAll(maf1.getMultPos());
		        if(toD.size()>0) {
		            Logger.global.info
		           // throw new RuntimeException
		            ("multi deletions "+chr+" "+name);
		        }
		        File emFile  = new File(f, name+".emphase");
		        File snpFile = new File(f, name+".snp");
		        int max = sdt2.loc.get(sdt2.loc.size()-1);
		        int min = sdt2.loc.get(0);
		        //int lenR = (int) Math.round( (double)(max - sdt2.loc.get(mid[0]))/1000.0);
		        //int lenL = (int) Math.round( (double)(sdt2.loc.get(mid[1])-min)/1000.0);
		        File pngFile = new File(f, name+".emphase.LD.PNG");//sdt2.loc.get(0)+""+sdt2.loc.get(sdt2.loc.size()-1));
		        //File newPng = new File(f, name+
		          //      "-"+lenL+"-"+lenR+"kb"+
		            //    ".LD.PNG");//
		        
		        File track = new File(f, name+".track");
		        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(track)));
		       // pw.print(mid[0]+1);
		      //  for(int i=mid[0]+2; i<=mid[1]+1; i++){
		       //     pw.print(" "+i);
		       // }
		        pw.println();
		     //   pw.close();
		       // sdt2.snpid.set(mid[0], "deletion");
		        //sdt2.snpid.set(mid[1], "deletion");
		       //for(int i=mid[0]; i<mid[1]; i++){
		       //    toD.add(i);
		       //}
		     //   if(mid[0] < mid[1]){
		      //         toDrop.add(mid[0]);
		     //        mid[1]--;
		     //     }
		       String blocktype = "GAB";
		       if(true){ 
		    	   sdt2.writeDickFormat(emFile, false, indiv, toD);
		           sdt2.writeSNPFile(snpFile, chr, false, toD);
		      // Logger.global.info("hwe "+f+" "+chr+"_"+mid[0]+"_"+mid[1]);
		       HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
		                "-info", snpFile.getAbsolutePath(),
		    //            "-dprime",
		               // "-blocks", track.getAbsolutePath(),
		             //   "-blockoutput", blocktype,
		           //    "-blockCutHighCI", "1.0",
		            //   "-blockCutLowCI","1.0",
		              //  "-ldcolorscheme", "RSQ",
		            //  "-hwcutoff", "0",
		     // "-png" ,"-n"
		               });
		        //File newPng = new File(f, name+
		        //      "-"+lenL+"-"+lenR+"kb"+
		          //    ".LD.PNG");//
		     /* HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
		                "-info", snpFile.getAbsolutePath(),
		                "-dprime",
		               // "-blocks", track.getAbsolutePath(),
		                "-blockoutput", blocktype,
		                //"-blockCutHighCI", "1.0",
		                //"-blockCutLowCI","1.0"
		                "-ldcolorscheme", "RSQ",
		             //   "-hwcutoff", "0.1"
		             "-png" ,"-n"
		               });*/
		       }
		        
		        
		        sdt2.writeDickFormat(emFile, false, indiv, null);
		        sdt2.writeSNPFile(snpFile, chr, true, null);
		      //  Thread.sleep(1000);
		        String blockExt = blocktype.equals("GAB") ? ".GABRIELblocks" :".4GAMblocks";
		        File blockFile = new File(f, name+".emphase"+blockExt);
		        int cnt =0;
		        while(!blockFile.exists() && cnt < 3){
		            Thread.sleep(1000);
		            blockFile = new File(f, name+".emphase"+blockExt);
		            cnt++;
		        }
		        if(blockFile.exists()){
		        BufferedReader br = new BufferedReader(new FileReader(blockFile));
		        String st = "";
		       /* while((st = br.readLine())!=null){
		            if(st.startsWith("BLOCK")){
		                String[] str = st.split(":")[1].trim().split("\\s+");
		                int[] res = new int[str.length];
		                boolean contains = false;
		                for(int jk=0; jk<res.length; jk++){
		                    res[jk] = Integer.parseInt(str[jk]);
		                    if(res[jk]==mid[0] || res[jk]==mid[1]) contains = true;
		                }
		                if(contains) return sdt2.loc.get(res[res.length-1]) - sdt2.loc.get(res[0]);
		            }
		        }*/
		        }
		        return 0;
		       // if(true) Thread.currentThread().wait(1000000);
		       // String[] cmd = new String[] {"mv", f.getName()+"/"+pngFile.getName(), f.getName()+"/"+newPng.getName()};
		       // StringWriter err = new StringWriter();
		       // ProcessTools.exec(cmd, null, err, err);
		      //  System.err.println(err.toString());
		         // TODO Auto-generated method stub
		         
		        // TODO Auto-generated method stub
		    }

	public static void calcLD(DataCollection obj1, File outF2, boolean phased, boolean avg) {
		
		        try{
		        	
		        	DataCollection obj = obj1;
		        	File outF1 = new File(outF2,"ld");
		        	outF1.mkdir();
		           // File outF = new File(out,"ld_result.txt");
		            String[] ld1 = Constants.ld();
		            
		           File outDir = new File(outF1, "outDir");
		   		 outDir.mkdir();
		   		 File f = new File(outF1, "ld.txt");
		   		 PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f, false)));
		            
		   		  pw.println("A\tB\tloc_A\tloc_B\trsid_A\trisd_B\tr2\tp");
		   		   
		         //   int[] first = get(ld1[0]);
		        //   int[] second = get(ld1[1]);
		            for(int i=0; i<ld1.length; i++){
		           // 	for(int j=0; j<second.length; j++){
		            	String[] ld_ = ld1[i].split(";");
		            	int ld_1 = Integer.parseInt(ld_[0]);
		            	int ld_2 = Integer.parseInt(ld_[1]);
		            	File outF = new File(outF1,"ld_"+getType(ld_1)+
		            			"_"+getType(ld_2)+".txt");
		            
		            HLALD hlald = new HLALD(outDir,pw, obj, obj,ld_1,ld_2,false, false, phased,true, avg,0,0);
		          hlald.run(false);
		            	//}
		            }
		           // writeEHH(sdt, new File(user, "sweep"),Constants.chrom0()+"", null, "");
		             pw.close();                                                                                  
		        }catch(Exception exc){
		        	exc.printStackTrace();
		        }
		
	}
	
	static int[] get(String str){
		if(str.equals("all")){
			int alles = 3;//Constants.alleles;
			int[] res = new int[alles-1];
			for(int i=0; i<res.length; i++){
				res[i] = i+1;
			}
			return res;
		}
		else{
			String[] str1 = str.split(";");
			int[] res = new int[str1.length];
			for(int i=0; i<res.length; i++){
				res[i] = Integer.parseInt(str1[i]);
			}
			return res;
		}
	}
	public static Double calcLD(DataCollection dc1, DataCollection dc2, File outF2,boolean avg, int index1, int index2){
		try{
		File outF1 = new File(outF2,"ld");
    	outF1.mkdir();
       // File outF = new File(out,"ld_result.txt");
       // int[][] ld1 = new int[][]{ };
      
        	int[] ld =new int[] {0,0};
        	File outF = new File(outF1,"ld_"+getType(ld[0])+"_"+getType(ld[1])+".txt");
		  HLALD hlald = new HLALD(outF,null, dc1, dc2,ld[0],ld[1],false, false, ld[0]==6 || ld[1]==6,false, avg, index1,index2);
		//  hlald.column2.replace = replace;
         Double res = hlald.calcLD(index1, index2);
   /* if(Double.isNaN(res)){
    	System.err.println("is nana");
    }*/
         return res;
      
	}catch(Exception exc){
    	exc.printStackTrace();
    	return null;
    }
	}
}
