package cnvtrans;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import assoc.Constants;
import data.util.CompressDir;


public class Assoc {
  
	public static void main(String[]args){
	//  System.err.println(args[0]);  
		
		try{
	    Constants.parse(args, Constants.class);
	 	File[] user1 = new File[Constants.inclDir().length];
	 	File[][] l = new File[user1.length][];
	 	File[] outFile=new File[user1.length];
	 	for(int k=0; k<user1.length; k++){
	 		user1[k] = new File(Constants.inclDir()[k]);
	 	}
	
	 	String combName = Constants.getCombName("cnv-",user1,"-");
	 	
	 	for(int k=0; k<user1.length; k++){
	 		File f1 = Constants.getFile(user1[k], Constants.experiment[0]);
		  String dir = f1.getAbsolutePath();//user1[k].getAbsolutePath()+"/"+Constants.dir(user1[k])[0];
		  //for (int j = 0; j < dir.length; j++) {
		
				int ind = dir.lastIndexOf('/');
				if(ind<0) {
					dir = "./"+dir;
					ind = dir.lastIndexOf('/');

				}
				/*
				File user = new File(user1 + "/"
						+ (ind < 0 ? "" : dir.substring(0, ind)));
				*/
				final String suff = dir.substring(ind);

				File avgDir = new File(user1[k], suff + "/avg");
				avgDir.mkdir();
				final String expt = suff.split("_")[0].split("/")[1];
				File cnvDir = new File(user1[k],suff+"/cnv");
				 l[k] = cnvDir.listFiles(new FileFilter() {

					@Override
					public boolean accept(File pathname) {
						return pathname.getName().endsWith(".txt");
								//&& pathname.getName()
								//		.startsWith(expt);
					}

				});
			
				 outFile[k]= new File(avgDir,combName);
				outFile[k].mkdir();
	 	}
		  
		
	 	int maxcnv=0;
		//File f = new File(cnvDir,expt+".txt");
	 	SortedSet<Integer> start_end= new TreeSet<Integer>();
		SortedMap<String, IndivCNV>[] m = new SortedMap[l.length];
		SortedSet<Integer> cnv_types = new TreeSet<Integer>();
		//PrintWriter[] pw = new PrintWriter[l.length];
		for(int k=0; k<l.length; k++){
		
			m[k] = new TreeMap();
			Iterator<String> iter = getBR(l[k]);
			//String[] header =iter.next().split("\\s+");
			process1(iter,m[k], k, cnv_types, start_end);
		}
		double[][] distr = new double[3][maxcnv+1];
		int multf = 100; //how much to multiply probs by
		int[][] distr100 = new int[3][maxcnv+1];
		int[] default_max_ind = new int[] {2,0,0};
		StringBuffer def_string = new StringBuffer();
		for(int k=0; k<default_max_ind.length; k++){
			Arrays.fill(distr100[k], 0);
			distr100[k][default_max_ind[k]]=multf;
			def_string.append(getString(distr100[k],-1));
			if(k<default_max_ind.length-1) def_string.append("\t");
		}
		for(int k=0; k<l.length; k++){
			PrintWriter pw = 	new PrintWriter(new BufferedWriter(new FileWriter(new File(outFile[k],"Name"))));
			pw.println("countAll\tstate.0\tstate.2");
			pw.println("chr\tstart\tend\tsnpid\tnosnps");
			pw.println("sample");
			pw.println(def_string.toString());
			pw.close();
		}
		for(int ind=0; ind<m.length; ind++){
			List<String> samples = readSamples(new File(outFile[ind].getParentFile().getParentFile(),"sample"));
			SortedSet<Integer> snpsset = readSNPS(outFile[ind].getParentFile().getParentFile());
			PrintWriter pw = 	new PrintWriter(new BufferedWriter(new FileWriter(new File(outFile[ind],"Samples"))));
			for(Iterator<String> it = samples.iterator(); it.hasNext();){
				pw.println(it.next());
			}
			pw.close();	
			String chr = 	"chr"+outFile[ind].getParentFile().getParentFile().getName().split("_")[2];
			Iterator<Integer> it = start_end.iterator();
			Integer[] startend = new Integer[] {  (it.next()),0};
		   String prev = "";
			PrintWriter snps = new  PrintWriter(new BufferedWriter(new FileWriter(new File(outFile[ind],"SNPS"))));
			for(int k=0; it.hasNext(); k++){
				int position = it.next();
				
				startend[1] =( (position));
				int nosnps = snpsset.tailSet(startend[0]).headSet(position+1).size();
				Iterator<String> it1= samples.iterator();
				double mid = ((double)startend[0]+(double)startend[1])/2.0;
				String rsid =String.format("R %8.6gmb", mid/(1000*1000)).replaceAll(" ", "");
			     if(rsid.equals("prev")) rsid =String.format("R %10.8gmb", mid/(1000*1000)).replaceAll(" ", "");
			     prev = rsid;
				PrintWriter pw1 = new  PrintWriter(new BufferedWriter(new FileWriter(new File(outFile[ind],rsid))));
				snps.println(chr+"\t"+startend[0]+"\t"+startend[1]+"\t"+rsid+"\t"+nosnps);
				
			
			//   System.err.println(mid);
				startend[0] = startend[1];
			
				int counter=0;
				for(counter=0;it1.hasNext();counter++){
					String indiv = it1.next();
					IndivCNV cnvs = m[ind].get(indiv);
					if(cnvs==null){
						pw1.println();
					}else{
					  for(int kj=0; kj<distr.length; kj++){
						Arrays.fill(distr[kj], 0.0);
					  }
					  if(cnvs.nocopies(position-1,distr,2)){
						  round(distr,distr100,multf);
						  for(int kk=0; kk<distr100.length; kk++){
							  String str  = getString(distr100[kk],default_max_ind[kk]);
							  pw1.print(str);
							  if(kk<distr100.length-1) pw1.print("\t");
						  }
					  }
					  pw1.println();
					}
					
				}
				if(counter<samples.size()){
					throw new RuntimeException("!!");
				}
				pw1.close();
			}
			snps.close();
		}
			
		for(int k=0; k<outFile.length; k++){
			(new CompressDir(outFile[k])).run();
		}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
    }

	private static void process1(Iterator<String> iter, Map<String, IndivCNV> m, int k, SortedSet<Integer> cnv_types,
			SortedSet<Integer> start_end) {
		String st ="";
		while(iter.hasNext()){
			st=iter.next();
			String [] str = st.split("\\s+");
			String indiv = str[0];
			IndivCNV cnvs = m.get(indiv);
			if(cnvs == null){
				cnvs = new IndivCNV(indiv,k);
				m.put(indiv, cnvs);
			}
			CNV cnv = new CNV(str,  new int[] {7,6,5,4,2,1},true);	
			cnv_types.add(cnv.type);
			cnvs.add(cnv);
			start_end.add(cnv.start);
			start_end.add(cnv.end);
		}	
		
	}

	private static String getString(int[] distr100, int default_max_ind) {
		if(default_max_ind>=0 && distr100[default_max_ind]==100) return "";
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<distr100.length; i++){
			if(distr100[i]>0) sb.append(distr100[i]);
			if(i<distr100.length-1) sb.append(",");
		}
		return sb.toString();
	}
	
	private static void round(double[][] distr, int[][] distr100, int mult) {
	  for(int j=0; j<distr.length; j++){
		int max_ind=0;
		double max = distr[j][0];
		int sum=0;
		for(int i=0; i<distr[j].length; i++){
			distr100[j][i] = (int) Math.round(mult*distr[j][i]);
			sum+=distr100[j][i];
			if(distr[j][i]>max){
				max_ind = i; max = distr[j][i];
			}
		}
		
		distr100[j][max_ind]= distr100[j][max_ind] + (mult-sum);
	  }
	}

	private static List<String> readSamples(File parentFile) throws Exception {
		List<String> l = new ArrayList<String>();
		File[] f = parentFile.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith("txt");
			}
			
		});
		for(int k=0; k<f.length; k++){
			BufferedReader br = new BufferedReader(new FileReader(f[k]));
			String st = "";
			while((st = br.readLine())!=null){
				l.add(st);
			}
			br.close();
		}
		return l;
	}
	
	private static SortedSet<Integer> readSNPS(File parentFile) throws Exception {
		SortedSet<Integer> l = new TreeSet<Integer>();
		File[] f = parentFile.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.getName().startsWith("snp") && pathname.getName().endsWith("txt");
			}
			
		});
		//for(int k=0; k<f.length; k++){
		int k=0;
		if(f.length==0){
			try{
				throw new RuntimeException("no snp files in "+parentFile.getAbsolutePath()
						+"\n"+Arrays.asList(parentFile.list()));
			}catch(Exception exc){
				exc.printStackTrace();
				System.exit(0);
			}
		}
			BufferedReader br = new BufferedReader(new FileReader(f[k]));
			String st = "";
			while((st = br.readLine())!=null){
				l.add(Integer.parseInt(st.split("\\s+")[1]));
			}
			br.close();
		//}
		return l;
	}

	public static Iterator<String> getBR(final File[] l) throws Exception{
		return new Iterator<String>(){
			int i =0;
			BufferedReader br = new BufferedReader(new FileReader(l[0]));
			String st = br.readLine();
			
			@Override
			public boolean hasNext() {
				return st!=null;
			}

			@Override
			public String next() {
			   String res = st;
			   try{
			   st = br.readLine();
			  i++;
				  for(;i<l.length && st==null; i++){
					  br.close();
					   br = new BufferedReader(new FileReader(l[i]));
					  st =  br.readLine();
					   st = br.readLine();
					//   if(st!=null){
						//   System.err.println(st);
					  // }
				   }
			  i--;
				// System.err.println(i); 
			   }catch(Exception exc){
				   exc.printStackTrace();
			   }
			   return res;
			}

			@Override
			public void remove() {
				// TODO Auto-generated method stub
				
			}
		};
	}
}
