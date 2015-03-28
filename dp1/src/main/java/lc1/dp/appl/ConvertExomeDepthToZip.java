package lc1.dp.appl;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import lc1.util.CompressDir;

public class ConvertExomeDepthToZip {
 public static void main(String[] args){
	//if(true) System.exit(0);
	 try{
		 File f = new File(args[0]);
		// File f1 = new File(f,args[0]);
		// if(!f1.exists()) f1.mkdir();
	//	if((new File(f,"avg")).exists()){
		 main(f, new File(args[1]), false);
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

 final BufferedReader br;
 //final File[] file;
  String st;
 String[] str0; //first is files, second samples in file
 //final int[] num;
 final Map<String, Double> avg_depth  = new HashMap<String, Double>();
 public ConvertExomeDepthToZip(File f, File avgdepth) throws Exception{
	String fname = f.getName();
	fname = fname.endsWith(".gz") ? fname.substring(0,fname.indexOf(".gz")) : fname;
	 String[] str =  fname.split("_");
	 int di = Arrays.asList(str).indexOf("depth");
	 this.start = Integer.parseInt(str[di+1]);
	 this.end = Integer.parseInt(str[di+2]);
	 this.chr = f.getParentFile().getName().substring(3);
	 br =this.getBR(f);// new BufferedReader(new FileReader(f));
	 
	 BufferedReader br1 = new BufferedReader(new FileReader(avgdepth));
	String st = br1.readLine();
	while((st = br1.readLine())!=null){
		String[] str1 = st.split("\\s+");
		try{
		avg_depth.put(str1[0], Double.parseDouble(str1[1]));
		}catch(Exception exc){
			avg_depth.put(str1[0], Double.NaN);
		}
	}
	br1.close();
 }
 public boolean readLine() throws Exception{
		 st = br.readLine();
		 if(st==null) return false;
			 str0 = st.split("\\s+");
	  return true;
 }
 final int start;
 final int end;
 int incr = 10;
 
 private  String getName(int i){
	 int st = this.start + i*incr;
		String nme = this.chr+"_"+st;
		return nme;
 }
 
 static List<String> samples_all = null;
 
 
public  void mainTranspose(File firstFile, File currentFile) throws Exception{
	OutputStreamWriter pw = null;
		List<String>samples = new ArrayList<String>();
		List<String> snps = new ArrayList<String>();
	//	List<String> snps = new 
	//	List<String> samples = new ArrayList<String>();
		boolean first = true;
		List<String>[] l=null;
		int[] alias = null;
	//	List<String> torem = new ArrayList<String>();
		for(int i=0; readLine(); i++){
			if(first){
				first = false;
				l = new List[str0.length-1];
				double incr1 =  ((double)this.end+1 - (double)this.start)/(double)l.length;
				incr =(int) Math.ceil(incr1);
				for(int k=0; k<l.length; k++){
					l[k] = new ArrayList<String>();
				}
				for(int k=0; k<l.length; k++){
					int st = this.start + k*incr;
					int end = st +incr;
					String nme = this.getName(k);
					snps.add("chr"+this.chr+"\t"+st+"\t"+end+"\t"+nme+"\n");
				}
			}
			String sample= str0[0];
			Double avgdepth = this.avg_depth.get(sample); 
			if(avgdepth==null){
				avgdepth = Double.NaN;
				System.err.println("could not find avg depth for "+sample+" will remove this sample");
			}else{
				samples.add(sample+"\t"+avgdepth);
			    for(int k=1; k<str0.length; k++){
				 // double reld = Double.parseDouble(str0[k])/avgdepth;
				  l[k-1].add(str0[k]);
			    }
			}
		}
	
		if(samples_all!=null){
			if(!samples_all.equals(samples)){
				alias = new int[samples.size()];
				int mismatch_index = -1;
				for(int i=0; i<samples.size(); i++){
					alias[i] = samples.indexOf(samples_all.get(i));
					if(mismatch_index<0 && alias[i]!=i) mismatch_index = i;
				}
				if(mismatch_index>=0) System.err.println("warning depth data not in same order, ie index of  "+
						samples_all.get(mismatch_index).split("\\s+")[0]
						 +" is at index "+mismatch_index+" and "+alias[mismatch_index]+
						 " between files "+firstFile+" vs "+currentFile);
			//	System.err.println("warning ")
				//throw new RuntimeException("mismatch of samples between depth files");
			}
			else alias = null;
		}else{
			alias = null;
			samples_all = samples;
			pw = compress.getWriter("Samples",true);
			//new PrintWriter(new BufferedWriter(new FileWriter(new File(dir1,"Samples"))));
		    for(int i=0; i<samples.size(); i++){
			   pw.write(samples.get(i)+"\n");
		    }
		    compress.closeWriter(pw);
		}
		
	    for(int k=0; k<l.length; k++){
	    	pw = 	compress.getWriter(this.getName(k),true);
	    	if(alias==null){
		    for(int i=0; i<l[k].size(); i++){
			 pw.write(l[k].get(i)+"\n");
		  }
	    	}
	    	else{
	    		for(int i=0; i<l[k].size(); i++){
	   			 pw.write(l[alias[k]].get(i)+"\n");
	   		  }
	    	}
		  compress.closeWriter(pw);
	    }
		OutputStreamWriter pw_snp = 	compress.getWriter("SNPS",false);
		for(Iterator<String> it = snps.iterator(); it.hasNext();){
			pw_snp.write(it.next());
		}
		compress.closeWriter(pw_snp);
}
final String chr;	
public static CompressDir compress=null; 


private static BufferedReader getBR(File f) throws Exception{
	return f.getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f)))) : 
		 new BufferedReader(new FileReader(f));
}

public static void main(File dir, File avgd, boolean remove) {
	 try{
		 File dir1= new File(dir.getParentFile(),dir.getName().substring(3));
		 if(dir1.exists() && dir1.length()>0){
			 System.err.println("exists");
			 //System.exit(0);
		 }
		  compress = (new CompressDir(dir1));
			OutputStreamWriter pw = 	compress.getWriter("Name",true);
			pw.write("depth\n");
			pw.write("chr\tstart\tend\tsnpid\n");
			pw.write("id\tavgdepth\n");
			compress.closeWriter(pw);
			
		  final String start ="depth"; 
		 File[] f =  dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File arg0) {
				// TODO Auto-generated method stub
				return arg0.getName().indexOf(start)>=0;
			}
			  
		  });
		 PrintWriter regions = new PrintWriter(new BufferedWriter(new FileWriter(new File(
				 dir.getParentFile(), "regions_"+dir1.getName()+".txt"))));
		 for(int i=0; i<f.length; i++){try{
			 ConvertExomeDepthToZip crc = new ConvertExomeDepthToZip(f[i], avgd);
			// crc.s
			 regions.println("chr"+crc.chr+"\t"+crc.start+"\t"+crc.end);
			 crc.mainTranspose(f[0], f[i]);
		 }catch(Exception exc){
			 exc.printStackTrace();
			 System.err.println("problem with "+f[i].getName());
			 
		 }
		 
		
		 }
		 regions.close();
		 compress.run();
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
	
 }
}
