package lc1.dp.appl;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import lc1.util.CompressDir;
import lc1.util.Constants;

public class ConvertInversionCountsFileToZip {
 public static void main(String[] args){
	
	 try{
		 File f = new File(args[0]);
		// File f1 = new File(f,args[0]);
		// if(!f1.exists()) f1.mkdir();
	//	if((new File(f,"avg")).exists()){
		 main(f, false);
		//}else{
		//	mainImpute(f1,args[1], Integer.parseInt(args[2]));
			
	//	}
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
 }
	
 static int thresh = 1;
 

 final BufferedReader[] br;
 final File[] file;
 final String[][] st;
 //final String[][] str0,str1,str2; //first is files, second samples in file
 final int[] num;
 final int[] max;
 final int[] indiv;
 final double[] depth_per_indiv;
 
 public ConvertInversionCountsFileToZip(File f, final String pref) throws Exception{
	 file = f.listFiles(new FileFilter(){
	
		public boolean accept(File arg0) {
			return arg0.getName().startsWith(pref);
		}
		 
	 });
	 this.start = Integer.parseInt(pref.split("-")[0]);
	 br = new BufferedReader[file.length];
	 st = new String[file.length][];
	
	 
	 num = new int[file.length];
	 max = new int[file.length];
	 indiv = new int[file.length];
	 depth_per_indiv = new double[file.length];
	 this.indices = new List[file.length];
	// str0 = new String[file.length][];
	// str1 = new String[file.length][];
	 //str2 = new String[file.length][];
	refresh();
 }
 
 public void refresh() throws Exception{
	 for(int k=0; k<st.length; k++){
		 br[k] = new BufferedReader(new FileReader(file[k]));
	 }
 }
 public void close() throws Exception{
	 for(int k=0; k<st.length; k++){
		 br[k].close();
	 }
 }
 int minmax=0;
 public void getExtraInfo(){
	 for(int k=0; k<br.length; k++){
	 for(int j=0; j<num[k]; j++){
		 int v = Integer.parseInt(st[k][j]);
		  max[k]+= v;
		
		  if(v>0){
			  indiv[k]++;
			  indices[k].add(j);
		  }
	/*	 if( val> max[k]){
			 max[k] = val;
		 }*/
	 }
	 this.depth_per_indiv[k] = indiv[k]==0 ? Double.NaN :(double) max[k] /(double) indiv[k];
	 if(max[k] < minmax) minmax = max[k];
 }
 }
 
 public boolean readLine() throws Exception{
	
	 for(int k=0; k<br.length; k++){
		 String str = br[k].readLine();
		 st[k] = str==null ? null : str.split("\\s+");
		 minmax =Integer.MAX_VALUE;
		 indices[k] = new ArrayList<Integer>();
	//	 indices.
		 if(st[k]!=null){
			 num[k] = st[k].length;
			 max[k] =0;
			 indiv[k] =0;
			 for(int j=0; j<num[k]; j++){
				 int v = Integer.parseInt(st[k][j]);
				  max[k]+= v;
				//  indices[k].add(j);
				 // if(v>0) indiv[k]++;
			/*	 if( val> max[k]){
					 max[k] = val;
				 }*/
			 }
			// this.depth_per_indiv[k] = indiv[k]==0 ? Double.NaN :(double) max[k] /(double) indiv[k];
			 if(max[k] < minmax) minmax = max[k];
		 }
	//	 String st0 = st[k];
	//	 String st1 = br[k].readLine();
	//	 String st2 = br[k].readLine();
		 
		/* if(st[k]!=null){
			 str0[k] = st[k].split("\\s+");
			 if(!st[k].equals(st[0])) throw new RuntimeException("!!");
			 str1[k] = st1.split("\\s+");
			 str2[k] = st2.split("\\s+");
			 if(num[k]==0) num[k] = str1[k].length;
			 else if(str1[k].length!=num[k]) throw new RuntimeException("inconsistent lengths");
			// System.err.println(k+"_"+num[k]);
	 }*/
	 }
	 //System.err.println(max[0]+" "+max[1]);
	 if(st[0]==null) return false;
	 else return true;
 }
 
 final int start;
 final List<Integer>[] indices;
 static List<String> snps = new ArrayList<String>();
	static Map<Integer,String> snps1 = new HashMap<Integer,String>();
	
	public  List<Integer> findPeaks(int len) throws Exception{
		List<Integer> maxl = new ArrayList<Integer>();
		for(int kk=0;readLine();kk++){
			maxl.add(max[1]);
		}
		if(len>=50){
			return Arrays.asList(new Integer[] {Constants.getMax(maxl.toArray(new Integer[0]))});
		}
		int sz=  maxl.size();
		int lastAdded = 0;
		List<Integer> index = new ArrayList<Integer>();
		for(int k=0; k<sz; k++){	
			boolean gt = true;
			boolean gte = true;
			int va = maxl.get(k);
			for(int j=Math.max(0,k-len); j<Math.min(sz,k+len); j++){
				if(j==k) continue;
				int va1 = maxl.get(j);
				if(va < va1) {
					gt = false;
					gte = false;
				}else if(va==va1){
					gt = false;
				}
			}
			if(gt || gte && lastAdded < va){
				lastAdded = va;
				index.add(k);
			}
		}
		return index;
	}
	
 public  void mainTranspose(boolean firstTime, List<Integer> indices) throws Exception{
	  OutputStreamWriter pw  = null;
		if(firstTime){
			pw= 	compress.getWriter("Name",true);
			pw.write("countA\tcountB\n");
			pw.write("chr\tstart\tend\tsnpid\n");
			pw.write("id\tavgdepth\n");
			compress.closeWriter(pw);
		}
		//List<String>samples = new ArrayList<String>();
		
		
	//	compress.closeWriter(pw);
		
		boolean first = true;
	
		for(int kk=0;readLine();kk++){
			if(indices ==null || indices.contains(kk)){
			int kk1 = kk;
			String snp = "" +(start+"_"+kk1);//str0[0][0].substring(3)+"_"+str0[0][1];
			
			int st =(start+kk1);// Integer.parseInt(str0[0][1]);
			int end = st +20;
			
			
			try{
				this.getExtraInfo();
				List<Integer> indicesL = indiv[1] < indiv[0] ? this.indices[1] : this.indices[0];
				String string = chr+"\t"+st+"\t"+end+"\t"+snp+"\t"+max[0]+"\t"+max[1]+"\t"+indiv[0]+"\t"
				+indiv[1]+"\t"+String.format("%5.2g",depth_per_indiv[0])+"\t"+
				String.format("%5.2g",depth_per_indiv[1])+"\t"+indicesL;
				
				if(!snps1.containsKey(st)){
					snps.add(string);
					//System.err.println("duplicate\n "+snps1.get(st)+"\n"+string);
				//	if(true) System.exit(0);
				
			pw = compress.getWriter(snp,true);
			for(int j=0; j<num[0]; j++){
			 for(int k=0; k<this.st.length; k++){
			  	pw.write(this.st[k][j]);
				 pw.write(k<this.st.length-1 ? "\t":"\n");
			 }
			}
			compress.closeWriter(pw);
				}
				snps1.put(st,string);
			}catch(ZipException exc){
				System.err.println(exc.getMessage());
			}
			
			}
			/*if(first){
				first = false;
				for(int j=0; j<br.length; j++){
					 for(int k=0; k<str1[j].length; k++){
						 samples.add(file[j].getName().split("\\.")[0]+"_"+k);
					 }
				}
			}*/
		}
	
		/*pw = compress.getWriter("Samples",true);
		//new PrintWriter(new BufferedWriter(new FileWriter(new File(dir1,"Samples"))));
	    for(int i=0; i<samples.size(); i++){
		   pw.write(samples.get(i)+"\n");
	    }
	    compress.closeWriter(pw);
		*/
		
}
 
 public static void finish() throws Exception{
	 OutputStreamWriter pw = 	compress.getWriter("SNPS",true);
		for(int i=0; i<snps.size(); i++){
			pw.write(snps.get(i)+"\n");
		}
		compress.closeWriter(pw);
 }
static String chr;	
public static CompressDir compress=null; 


private static BufferedReader getBR(File f) throws Exception{
	return f.getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f)))) : 
		 new BufferedReader(new FileReader(f));
}

public static void main(File dir1, boolean remove) {
	 try{
		List<String> dir = Arrays.asList(dir1.list(new FilenameFilter(){
			public boolean accept(File arg0, String arg1) {
				return arg1.endsWith(".a");
			}
		 }));
		 
		Collections.sort(dir,new Comparator<String>(){

			@Override
			public int compare(String arg0, String arg1) {
			  Integer i1 =  Integer.parseInt(arg0.split("-")[0]);
			  Integer i2 =  Integer.parseInt(arg1.split("-")[0]);
			  return i1.compareTo(i2);
			}
			
		});
		 compress = (new CompressDir(dir1));
		 compress.delete = false;
		 for(int k=0; k<dir.size(); k++){
			  String str = dir.get(k).split("\\.")[0];
		 ConvertInversionCountsFileToZip crc = new ConvertInversionCountsFileToZip(dir1,str);
		 List<Integer> l = crc.findPeaks(50);
		 crc.close();
		 crc.refresh();
		 crc.mainTranspose(k==0,l);
		 crc.close();
		 }
		 ConvertInversionCountsFileToZip.finish();
		compress.close();
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
	
 }
}
