package ensj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipFile;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import conversion.Compressor;
import conversion.OptionBuild;
import conversion.SplitByLocation;
import conversion.Utils;

public class Regression1 {
	static boolean print = false;
	static  Algebra lg = new Algebra();
	static DoubleMatrix2D x0, Q;
	PrintWriter gc;
	static final int maxpow = 4;
	static boolean useAvgForPos =true; //whether to just regress avg LRR at each position, or include each individual
	static {
		int len = maxpow+1;
		double scale = 1e-2;
		x0 =  new DenseDoubleMatrix2D(len,1);
		Q = new DenseDoubleMatrix2D(len,len);
		for(int k=0; k<len; k++){
			Q.set(k, k, scale*Math.pow(1.2,k) );
		}
		
	}
	ZipFile current;
	File zf = null;
	String currentChr;
	int last_col;
	List<String>  header,header_snp,header_sample;
	double mint =-0.3;
	double maxt = 0.3;
	List<String> indiv;
	int num_indiv;
	static int start =0;//Integer.MAX_VALUE;
	static int end =  Integer.MAX_VALUE;
	public static void main(File dir1){
		try{
			
			// File dir1 = new File(args.length>0 ? args[0] : System.getProperties().getProperty("user.dir"));
	        	Regression1 reg = new Regression1(dir1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	class Inner{
		Set<String> toinclude = new HashSet<String>();
		DoubleMatrix2D X;
		DoubleMatrix2D Y;
		Set<Integer> na = new HashSet<Integer>();
		int cnt=0;
		int index=0;
		DoubleMatrix1D result;
		
		
	
		//double constant, gc_coeff, gc2_coeff;
		final boolean autosome;
		public Inner(File dir, boolean b) {
			this.autosome = b;
		
			 String[] name = dir.list(new ZipFileFilter(b));
				 for(int i=0; i<name.length; i++){
					 toinclude.add("chr"+name[i].split("\\.")[0]);
					
				 }
				
		}
		public boolean contains(String st){
			return this.toinclude.contains(st);
		}
		public void init(){
			this.X = new DenseDoubleMatrix2D(cnt,maxpow+1);
			this.Y = new DenseDoubleMatrix2D(cnt,1);
		}
		String currentChr ="";
		int chrcnt=0;
		double chrsum=0;
		int chr_start_index=0;
		
		private void correctAvgs(){
			if(chrcnt>0){
				chrsum = chrsum/(double) chrcnt;
				for(int k=chr_start_index; k<index; k++){
					Y.set(k, 0, Y.get(k,0) - chrsum); //subtract average
				}
			}
			chrsum=0;
			chr_start_index=index;
			chrcnt=0;
		}
		
		public void set(String chr, double gc, double lrr){
			if(!chr.equals(currentChr)){
				correctAvgs();
				currentChr = chr;
				
			}
			if(Double.isNaN(lrr)) this.na.add(index);
			else{
				chrsum+=lrr;
				chrcnt++;
			}
			/*else if(Regression1.this.gc!=null){
				
				
				Regression1.this.gc.print(lrr);
				for(int k=0; k<4; k++){
					Regression1.this.gc.print("\t"+Math.pow(gc,k+1));
				}
				Regression1.this.gc.println();
			}*/
			X.set(index,0, 1);
			for(int k=1; k<=maxpow; k++){
				X.set(index,k, Math.pow(gc,k));
			}
			
			Y.set(index, 0, lrr);
			index++;
			
		}
		public double calculate(double gcc){
			if(Double.isNaN(gcc)) return 0.0;
			double sum=result.get(0);
			for(int k=1; k<=maxpow; k++){
				sum+= result.get(k) *Math.pow(gcc, k);
			}
			return  sum;
			//constant +gc_coeff*gcc +gc2_coeff*Math.pow(gcc, 2);
		}
		public void solve(){
			this.correctAvgs();
			if(index!=cnt) {
				throw new RuntimeException("not set right number of variables");
			}
			else if(cnt==0){
				//constant =0;
				//gc_coeff =0;
				//gc2_coeff =0;
			}
			else{
				//X.viewSelection(arg0)
				DoubleMatrix2D X1 = X;
				DoubleMatrix2D Y1 = Y;
			    if(this.na.size()>0){
			    	int[] col = new int[X.columns()];
			    	for(int k=0; k<col.length; k++){
			    		col[k] = k;
			    	}
			    	int[] row = new int[X.rows()-na.size()];
			    	int k1 =0;
			    	for(int k=0; k<X.rows(); k++){
			    		if(!na.contains(k)){
			    			row[k1] =k;
			    			k1++;
			    		}
			    	}
			    	
			    	X1 = X.viewSelection(row, col);
			    	Y1 = Y.viewSelection(row, new int[]{0});
			    }
			   
			//	DoubleMatrix2D XT = lg.transpose(X1);
				try{
				 result = Regression1.solve(X1,Y1).viewColumn(0);
					//lg.solve(lg.mult(XT, X1),lg.mult(XT,Y1));
				 System.err.print("v = c("+result.get(0));
				for(int k=1; k<=maxpow; k++){
					System.err.print(","+result.get(k));
				}
				System.err.println(")");
				System.err.println(result);
				//constant = result.get(0, 0);
			//	gc_coeff = result.get(1, 0);
			//	if(!gconly) gc2_coeff =  result.get(2,0);
				}catch(Exception exc){
				   System.err.println("PROBLEM WITH REGRESSION");
//					System.err.println(X);
	//			System.err.println(Y);
				//	throw new RuntimeException("!!");
				}
				
				 // DoubleMatrix2D res1 = lg.mult(X, result);
				//  System.err.println(res1);
			}
			
		}
		public void print(PrintWriter pw_coeff) {
			pw_coeff.println(autosome ? "autosome" : "sex chroms");
			pw_coeff.println(result);
		}
	}
	Inner autosome,sexchrom;
	
	
	static class ZipFileFilter implements FilenameFilter{
		boolean autosome = false;
	
		ZipFileFilter(boolean autosome){
			this.autosome = autosome;
		}
		public boolean accept(File dir, String name) {
			boolean isautosome = (name.indexOf("X")>=0 || name.indexOf("Y")>=0);
			if(autosome) return name.endsWith("zip") && !isautosome;
			else return name.endsWith("zip") && isautosome;
		}
	}
	static String NA = "NaN";
	
	public void getCount(File buildF) throws Exception{
		 int cnt1 = 0;
		 BufferedReader br =Utils.getBufferedReader(buildF);
		 String st;
		 while((st=br.readLine())!=null && cnt1<end){
			 if(cnt1>=start){
				 int fi = st.indexOf('\t');
				 int li = st.lastIndexOf('\t');
				
				 if(!st.substring(li).equals(NA)){
				 if(!sexchrom.contains(st.substring(0,fi)) ){
					 autosome.cnt++;
				 }
				 else{
					 sexchrom.cnt++;
					 
				 }
				 }
			 }
			 cnt1++;
			 
		 }
		 br.close();
	}
	
	
	
	Map<String, Integer> centromere;
	
	/** this deals with _p, _q inconsistency in build file */
	private String process(String line){
		
		if(line==null) return null;
		int ind = line.indexOf('\t');
		String chr = line.substring(3,ind);
		if(centromere.containsKey(chr)){
			String p  = Integer.parseInt(line.substring(ind+1,line.indexOf('\t',ind+1))) < centromere.get(chr) ? "p" : "q";
			String res = line.replace("chr"+chr, "chr"+chr+""+p);
			return res;
		}
		else return line;
		
	}
	
	//used for calculating the average for each chrom
	
	
	public Regression1(File dir) throws Exception{
		 File buildF = new File(dir.getParentFile(), OptionBuild.build+".gc.txt");
		 this.dir = dir;
		
		 if(print) this.gc = new PrintWriter(new BufferedWriter(new FileWriter("gc.txt")));
		String[] f_ = dir.list(new FilenameFilter(){

			@Override
			public boolean accept(File dir, String name) {
				return name.indexOf("p.")>=0 || name.indexOf("q.")>=0;
			}
			
		});
	     centromere =
	    	 f_.length==0 ? new HashMap():
	    	 SplitByLocation.readCentromere(new File(dir.getParentFile(), "karyotypes_"+OptionBuild.build+".gc.txt"));
	     
	     for(Iterator<String> it = centromere.keySet().iterator(); it.hasNext();){
	    	 File f = new File(dir, it.next()+".zip");
	    	 if((f).exists()){
	    		 it.remove();
	    	 }
	     }
	     System.err.println("centromeres at "+centromere);
	     if(centromere.size()==0) System.err.println("centro file is "+new File(dir.getParentFile(), "karyotypes_"+OptionBuild.build+".gc.txt"));
		autosome = new Inner(dir, true);
		 sexchrom = new Inner(dir, false);
		// if(!buildF.exists()) throw new RuntimeException("build file does not exist "+buildF.getAbsolutePath());
		 BufferedReader br =Utils.getBufferedReader(buildF);
		 String st = process(br.readLine());
		 String[] str = st.split("\\s+");
		 last_col = str.length -1;
		 this.currentChr = str[0].substring(3);
		 
		 File currentF = new File(dir, currentChr+".zip");
		 //
		 System.err.println("file is "+currentF.getAbsolutePath());
		 while(!currentF.exists()){
			 st = process(br.readLine());
			  str = st.split("\\s+");
			 if(!str[0].substring(3).equals(currentChr)){
				 this.currentChr = str[0].substring(3);
				 currentF = new File( dir,currentChr+".zip");
				 System.err.println("file is "+currentF.getAbsolutePath());
			 }
			// last_col = str.length -1;
			 
			// File currentF = new File(dir, currentChr+".zip");
		 }
//		 if(currentF.exists()){
		 System.err.println("file is "+currentF.getAbsolutePath());
		 this.current = new ZipFile(currentF);
	//	 }
		 br.close();
		
		 List<String> headers=  Compressor.getIndiv(current, "Name", (int[])null);
         header =Arrays.asList( headers.get(0).split("\t"));
         header_snp = Arrays.asList(headers.get(1).split("\\s+"));
         header_sample = Arrays.asList(headers.get(2).split("\t"));
        this.snp_index = header_snp.indexOf("id");
        indiv = Compressor.getIndiv(current, "Samples", header_sample.indexOf("id"));
        if(indiv==null){
        	 indiv = Compressor.getIndiv(current, currentChr+"/Samples", header_sample.indexOf("id"));
        }
	       num_indiv = indiv.size();
	       System.err.println("num indiv is "+indiv.size());
	       int log_ind =-1;
	       for(int i=0; i<header.size(); i++){
	    	   if(header.get(i).indexOf("Log")>=0 && header.get(i).indexOf("R")>=0){
	    		   log_ind =i;
	    		   break;
	    	   }
	       }
        readData(log_ind, buildF, false);
       
		 sexchrom.init();
		 autosome.init();
		 
		 br = Utils.getBufferedReader(buildF);
	
         for(int i=0; i<header_snp.size(); i++){
        	 header_snp.set(i, header_snp.get(i).trim());
         }
		
		
		readData(log_ind, buildF, true);
		 if(gc!=null){
	        	gc.close();
	        }
		 PrintWriter pw_coeff = new PrintWriter(new BufferedWriter(new OutputStreamWriter((new FileOutputStream(new File(dir, "coefficients.txt"))) )));
		   
		 autosome.solve();
		  autosome.print(pw_coeff);
		 try{
		 sexchrom.solve();
		  sexchrom.print(pw_coeff);
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
		
	   
	     pw_coeff.close();
		// double[][] x = new double[][] {getOnes(10), getRandom(10)};//new double[][] {{1.0, 2.0, 3.0, 4.0}};
		// double[][]y = new double[][]{ getRandom(10)};
		 
				 // System.err.println(res1);
		this.printres(buildF);
	
			//  System.err.println(res1);
		    // System.err.println( result.getQuick(1, 0));
	 }
	
	public void printres(File buildF) throws Exception{
		 PrintWriter pw = null;
		 this.currentChr="";
		 BufferedReader br = Utils.getBufferedReader(buildF);
		 String st = "";
		 for(int row=0; (st = process(br.readLine()))!=null && row < this.end;row++){
			 if(row>=start){
				 String[]  str = st.split("\\s+");
				 String chr = str[0].substring(3);
				 if(!chr.equals(this.currentChr)){
					 currentChr = chr;
					 if(current!=null) current.close();
					 zf = new File(dir, currentChr+".zip");
					 if(zf.exists()){
						 try{
						 this.current = new ZipFile(zf);
						 }catch(Exception exc){
							 exc.printStackTrace();
							 current = null;
						 }
					 }
					 else if (currentChr.equals("X")){
						 File x_f = new File(dir, currentChr+"_F"+".zip");
						 File x_m = new File(dir, currentChr+"_M"+".zip");
						 zf = x_f.exists() ? x_f : x_m;
						 if(zf.exists()){
							 try{
							 current = new ZipFile(zf);
							 }catch(Exception exc){
								 exc.printStackTrace();
								 current = null;
							 }
						 }
						 else current = null;
					 }
					 else current = null;
					 if(pw!=null) pw.close();
					 File pwf = new File(dir, "SNPS_"+chr);
					// pwf.deleteOnExit()
					 pw =  new PrintWriter(new BufferedWriter(new OutputStreamWriter(
							(new FileOutputStream(pwf)) )));
				 }
					if(current!=null && (current.getEntry(str[this.snp_index])!=null ||current.getEntry(this.currentChr+"/"+str[this.snp_index]+".txt")!=null ) ){
					 double gc = Double.parseDouble(str[this.last_col]);
					 double base = sexchrom.contains(str[0])  ? sexchrom.calculate(gc) : autosome.calculate(gc);
					 pw.print(st);pw.print("\t");pw.println(base);
					}
			 }
			
		 }
		 pw.close();
		 br.close();
	}
	
	final File dir; 
	final int snp_index;

	private void readData(int log_ind, File buildF, boolean extract) throws Exception{
		 String st;
		 System.err.println("log ind is "+log_ind);
		 BufferedReader br = 	Utils.getBufferedReader(buildF);
		 String[] res = new String[num_indiv];
		// PrintWriter out_pw = null;//extract  ? new PrintWriter(new BufferedWriter(new FileWriter("test.txt"))) : null;
		 for(int cnt2=0;  (st = process(br.readLine()))!=null &&  cnt2 < end;cnt2++){
			 if(cnt2>=start){
				 String[] str = st.split("\\s+");
			 String chr = str[0].substring(3);
			 if(!chr.equals(this.currentChr)){
				 currentChr = chr;
				 System.err.println("chr is "+chr);
				if(current!=null) current.close();
				 zf = new File(dir, currentChr+".zip");
				 if(zf.exists()){
					 try{
					 this.current = new ZipFile(zf);
					 }catch(Exception exc){
						 exc.printStackTrace();
						 current = null;
					 }
					
				 }
				 else current = null;
			 }
			 if(str[this.last_col].equals(NA)) continue;
			 double gc = Double.parseDouble(str[this.last_col]);
			 String snpid = str[snp_index];
			 if(extract){
			 if(current!=null &&  Compressor.getIndiv1(current, snpid, log_ind, res, currentChr)){
			
			 double avg =0;
			 double cnt=0;
			 boolean sexch  = sexchrom.contains(str[0]);
			 { int i=0;
			 try{
				
			 for(i=0; i<res.length; i++){
				 double n = Double.parseDouble(res[i]);
				 
				 if(!Double.isNaN(n) && n>mint && n < maxt){
					 if(useAvgForPos){
						 cnt++;
						 avg+=n;
					 }
					 else{
						 if(sexch)  sexchrom.set(str[0], gc, n);
						 else autosome.set(str[0], gc, n);
					 }
					 
					
				 }
			 }
			 }catch(Exception exc){
				// System.err.println(res[i]+" ");//+res[i-1]);
				 exc.printStackTrace();
			 }
			 }
		//	 if(cnt==0) System.err.println("no probes matched");
			 avg = avg/(double)cnt;
			// if(out_pw!=null && !sexch){
				//	System.err.println("printing");
			//		out_pw.println(currentChr+"\t"+gc+"\t"+Math.pow(gc, 2)+"\t"+avg);
			//	}
			 if(useAvgForPos){
				 if(sexch)  sexchrom.set(str[0], gc, avg);
				 else autosome.set(str[0], gc, avg);
			 }
			// Y.setQuick(row, 1, );
		     }
			 }
			 else{
				 if(current!=null && (current.getEntry(snpid)!=null || current.getEntry(currentChr+"/"+snpid+".txt")!=null)){
					 if(useAvgForPos){
					 if( sexchrom.contains(str[0]))  sexchrom.cnt++;
					 else autosome.cnt++;
					 }
					 else if(Compressor.getIndiv1(current, snpid, log_ind, res, currentChr)) {
						 if( sexchrom.contains(str[0]))  sexchrom.cnt+= nonNa(res);
						 else autosome.cnt+= nonNa(res);
					 }
				 }
				 else{
					// System.err.println("not "+snpid);
				 }
			 }
			/* else{
				 System.err.println("error with "+snpid+" "+chr+" ");//+zf.getAbsolutePath());
				 double avg = Double.NaN;
				 if( sexchrom.contains(str[0]))  sexchrom.set(gc, avg);
				 else autosome.set(gc, avg);
			 }*/
		 }
		 }
	//	 if(out_pw!=null){
		//	 out_pw.close();
	//	 }
		 /*if(true && extract){
			 System.exit(0);
		 }*/
		 br.close();
		 if(current!=null && extract) current.close();
		
	}



	private int nonNa(String[] res) {
		int cnt=0;
		for(int k=0; k<res.length; k++){
			double n = Double.parseDouble(res[k]);
			if(!Double.isNaN(n) && n>mint && n < maxt) cnt++;
		}
		return cnt;
	}



	public static double[] getRandom(int len){
		 double[] vec = new double[len];
		 for(int i=0; i<vec.length; i++){
			 vec[i]  = Math.random();
		 }
		 return vec;
	 }
	 public static double[] getOnes(int len){
		 double[] vec = new double[len];
		 Arrays.fill(vec,1.0);
		 return vec;
	 }
	 
	 public static DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b
			  ) {
			
			   DoubleMatrix2D AT = lg.transpose(A);
			   DoubleMatrix2D prod = lg.mult(AT, A);
			   
			   for(int i=0; i<prod.rows(); i++){
				   for(int j=0; j<prod.columns(); j++){
					   prod.setQuick(i, j, prod.getQuick(i, j)+Q.getQuick(i, j));
				   }
			   }
			   DoubleMatrix2D b1 = b.copy();
			   DoubleMatrix2D Ax0 = lg.mult(A, x0);
			   for(int i=0; i<b1.rows(); i++){
				   b1.setQuick(i, 0, b1.getQuick(i, 0)-Ax0.getQuick(i, 0));
			   }
			   DoubleMatrix2D res =  lg.solve(prod,lg.mult(AT,b1));
			   for(int i=0; i<res.rows(); i++){
				   double v1 = res.getQuick(i, 0);
				  
				   res.setQuick(i, 0, v1+x0.getQuick(i, 0));
			   }
			   
			   return res;
		 }
}
