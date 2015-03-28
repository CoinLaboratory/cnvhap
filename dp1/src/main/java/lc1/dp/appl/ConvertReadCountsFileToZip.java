package lc1.dp.appl;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import lc1.util.CompressDir;

public class ConvertReadCountsFileToZip {
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
	
 
 public static void compressImpute(BufferedReader br, File dir,  String nme,  File karyfile, int max) throws Exception{
	 List<Integer> kary = new ArrayList<Integer>();
	 List<String> append = new ArrayList<String>();
		 getKary(karyfile,nme.substring(3), kary, append);
	
	 String st = br.readLine();
	outer: for(int ind =0; st!=null; ind++){
		int start =0;
	 outer1: for(int app1=0; st!=null && start < kary.get(ind); app1++){
	 // app1 = (int)Math.floor((double)kk/ (double)max);
	  File dir1 = new File(dir,nme+append.get(ind)+(app1==0 ? "" : app1));
	  int start_next = kary.get(ind);
	  compress = (new CompressDir(dir1));
	  OutputStreamWriter snps = compress.getWriter("SNPS", false);
		 String chr = dir1.getName();
		inner: for(int kk1=0; (st )!=null ; kk1++){
			String[] str = st.split(" ");
			//System.err.println(str[2]);
			 start = Integer.parseInt(str[2]);
			if(start > start_next || kk1>=max) break inner;
			snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"\t"+str[0]+"\n");
			OutputStreamWriter pw;
			try{
			pw= compress.getWriter(str[1], true);
			snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"\t"+str[0]+"\n");
			}catch(Exception exc){
				exc.printStackTrace();
				pw= compress.getWriter(str[1]+"_"+start, true);
				snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"_"+start+"\t"+str[0]+"\n");
			}
			for(int k=5;k<str.length; k+=3){
				StringBuffer sb = new StringBuffer();
				for(int j=0; j<3; j++){
					int pr = (int) Math.round(Double.parseDouble(str[k+j])*1000.0);
					sb.append(pr);
					if(j<2) sb.append(",");
				}
				sb.append("\n");
				pw.write(sb.toString());
			}
			compress.closeWriter(pw);
			st=br.readLine();
		}
		snps.close();
		OutputStreamWriter pw = 	compress.getWriter("Name",true);
		pw.write("geno\n");
		pw.write("chr\tstart\tend\tsnpid\timputed\n");
		pw.write("sample\n");
		pw.write("1000,0,0\n");
		compress.closeWriter(pw);
		compress.run();
	}
	}
		br.close();
		
		//compress.run();
 }
 
 private static void  getKary(File karyfile, String substring, List<Integer> res, List<String> app) throws Exception{
	String st = "";
	BufferedReader br = getBR(karyfile);
	while((st = br.readLine())!=null){
		String[] str = st.split("\\s+");
		if(str[0].equals(substring)){
			//Integer[] res = new Integer[str.length-1];
			for(int k=1; k<str.length; k++){
				res.add( Integer.parseInt(str[k]));
				if(k==1) app.add("p");
				else app.add(k+"");
			}
			res.add(Integer.MAX_VALUE-1);
			app.add("q");
			return;
		}
	}
	
	
}

 final BufferedReader[] br;
 final File[] file;
 final String[] st;
 final String[][] str0,str1,str2; //first is files, second samples in file
 final int[] num;
 
 public ConvertReadCountsFileToZip(File f) throws Exception{
	 file = f.listFiles(new FileFilter(){
		@Override
		public boolean accept(File arg0) {
			return arg0.getName().endsWith(".txt");
		}
		 
	 });
	 br = new BufferedReader[file.length];
	 st = new String[file.length];
	 num = new int[file.length];
	 str0 = new String[file.length][];
	 str1 = new String[file.length][];
	 str2 = new String[file.length][];
	 for(int k=0; k<st.length; k++){
		 br[k] = new BufferedReader(new FileReader(file[k]));
	 }
 }
 public boolean readLine() throws Exception{
	 for(int k=0; k<br.length; k++){
		 st[k] = br[k].readLine();
		 String st0 = st[k];
		 String st1 = br[k].readLine();
		 String st2 = br[k].readLine();
		 
		 if(st[k]!=null){
			 str0[k] = st[k].split("\\s+");
			 if(!st[k].equals(st[0])) throw new RuntimeException("!!");
			 str1[k] = st1.split("\\s+");
			 str2[k] = st2.split("\\s+");
			 if(num[k]==0) num[k] = str1[k].length;
			 else if(str1[k].length!=num[k]) throw new RuntimeException("inconsistent lengths");
			// System.err.println(k+"_"+num[k]);
		 }
	 }
	 if(st[0]==null) return false;
	 else return true;
 }
 
public  void mainTranspose() throws Exception{
		OutputStreamWriter pw = 	compress.getWriter("Name",true);
		pw.write("countA\tcountB\n");
		pw.write("chr\tstart\tend\tsnpid\n");
		pw.write("id\tavgdepth\n");
		compress.closeWriter(pw);
		List<String>samples = new ArrayList<String>();
		
		
		compress.closeWriter(pw);
		List<String> snps = new ArrayList<String>();
		boolean first = true;
		
		while(readLine()){
			String snp = str0[0][0].substring(3)+"_"+str0[0][1];
			int st = Integer.parseInt(str0[0][1]);
			int end = st +20;
			snps.add(str0[0][0]+"\t"+st+"\t"+end+"\t"+snp+"\t"+str0[0][2]);
			pw = compress.getWriter(snp,true);
			for(int j=0; j<br.length; j++){
			 for(int k=0; k<str1[j].length; k++){
			  	pw.write(str1[j][k]+"\t"+str2[j][k]);
				 pw.write("\n");
			 }
			}
			compress.closeWriter(pw);
			if(first){
				first = false;
				for(int j=0; j<br.length; j++){
					 for(int k=0; k<str1[j].length; k++){
						 samples.add(file[j].getName().split("\\.")[0]+"_"+k);
					 }
				}
			}
		}
	
		pw = compress.getWriter("Samples",true);
		//new PrintWriter(new BufferedWriter(new FileWriter(new File(dir1,"Samples"))));
	    for(int i=0; i<samples.size(); i++){
		   pw.write(samples.get(i)+"\n");
	    }
	    compress.closeWriter(pw);
		pw = 	compress.getWriter("SNPS",true);
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

public static void main(File dir, boolean remove) {
	 try{
		 ConvertReadCountsFileToZip crc = new ConvertReadCountsFileToZip(dir);
		 File dir1= new File(dir.getParentFile(),dir.getName().split("_raw")[0]);
		  compress = (new CompressDir(dir1));
		 crc.mainTranspose();
		compress.run();
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
	
 }
}
