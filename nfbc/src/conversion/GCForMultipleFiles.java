package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Sequence;
import org.ensembl.driver.SequenceAdaptor;

import ensj.GenesForRegion;

public class GCForMultipleFiles {
/* preconditions - each file must be sorted by chrom and then by pos  e.g. sort -k 2,2 -k 3,3g*/  
	static boolean getGC = true;
	final  int numRefs;
	
	
	
	public static void main(final String[] args){
		
		
		try{
			File  dir1 = new File(System.getProperty("user.dir"));
			String[]  ref = args.length==0 ? null : args[0].split(":");
		
			final List<String> ref1 = ref == null ? new ArrayList<String>(): Arrays.asList(ref);
			File[] in = dir1.listFiles(new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return 
					!pathname.isDirectory() &&
					!(pathname.getName().endsWith(".out") || pathname.getName().endsWith(".out.gz")) && (
					 ! ref1.contains(pathname.getName()));
				}
				
			});
			//if(ref ==null) ref = new String[] {in[0].
			GCForMultipleFiles gc = new GCForMultipleFiles(ref, in);
			gc.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
	
	 GenesForRegion example = new GenesForRegion();
	
	  CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
	
	GCForMultipleFiles(String[] ref, File[]in) throws Exception{
		 example.initialise();
			numRefs  = ref==null ? 0 : ref.length;
		
		
		//this.f1 = f1;
		if(ref!=null){
			f1= new File[in.length+ref.length];
		
			System.arraycopy(in, 0, f1, ref.length, in.length);
			for(int k=0; k<ref.length; k++){
				f1[k] = new File(in[0].getParentFile(), ref[k]);
				if(!f1[k].exists()) throw new RuntimeException("ref wrong "+f1[k].getAbsolutePath());
			}
		}
		else{
			f1 = in;
		}
		
		this.in = new BufferedReader[f1.length];
		this.out = new PrintWriter[f1.length];
		len = f1.length;
		currString = new String[len][];
		compPrev = new int[len];
		this.done = new LinkedList[len];
		List<String> chrs = Arrays.asList("1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:XY:Y:M".split(":"));
		File dir = in[0].getParentFile();
		ref_new = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "result.out"))));
		warning = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "warning.out"))));
		
		this.chrid = new int[len];
		this.snpid = new int[len];
		this.snpid1 = new int[len];
		this.snpid2 = new int[len];
		this.remid = new List[len];
		this.hasChr = new boolean[len];
		Arrays.fill(snpid, -1);
		Arrays.fill(snpid1, -1);
		Arrays.fill(snpid2, -1);
		Arrays.fill(chrid, -1);
		for(int k=0; k<len; k++){
			done[k] = new LinkedList();
			this.in[k] =Utils.getBufferedReader(f1[k]);
			this.out[k] = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, f1[k].getName().split("\\.")[0]+".out"))));
			List<String> l = Arrays.asList(this.in[k].readLine().split(sep));
			this.in[k].close();
			this.in[k] =Utils.getBufferedReader(f1[k]);
			remid[k] = new ArrayList<Integer>();
			System.err.println(k+" "+f1[k].getName());
			for(int j=0; j<l.size(); j++){
				boolean isint = true;
				try{
					int v = Integer.parseInt(l.get(j));
				}catch(Exception exc){
					isint = false;
				}
				
				if(chrid[k] < 0 && (l.get(j).startsWith("chr") || chrs.contains(l.get(j)))){
					chrid[k] = j;
					if(l.get(j).startsWith("chr")){
						hasChr[k] = true;
					}
					System.err.println("using as chrid: "+k+" "+l.get(j)+" "+hasChr[j]);
				}
				else if(snpid[k]<0 && isint){
					
						snpid[k] =j;
						System.err.println("using as snpid: "+k+" "+l.get(j));
					
				}
				else if(snpid1[k]<0 && isint){
						snpid1[k] =j;
						System.err.println("using as snpid1: "+k+" "+l.get(j));
				}
				else if(snpid2[k]<0){
					snpid2[k] = j;
					System.err.println("using as snpid2: "+k+" "+l.get(j));
				}
				else{
					remid[k].add(j);
				}
			}
			
		}
		
	}
	public void close() throws Exception{
		for(int k=0; k<len; k++){
			this.in[k].close();
			this.out[k].close();
		}
	}
	
	final BufferedReader[] in;
	final PrintWriter[] out;
	final PrintWriter ref_new;
	final String[][] currString;
	final int[] compPrev;
	final int len;
	
	String[] current;
	
	final String sep = "\t";
	final int[] chrid;
	final boolean[] hasChr;
	final int[] snpid;
	final int[] snpid1;
	final int gc_id=4;
	public void run() throws Exception{
		
		while(readLine()){
			double gc = Double.NaN;
			inner: for(int k=0; k<numRefs; k++){
				if(this.compPrev[k]==0){
				String[] curr0 = currString[k];
					if(curr0!=null){// && gc_id <currString[0].length){
						double d = Double.parseDouble(curr0[curr0.length-1]);
						if(!Double.isNaN(d) && d<1 && d>0){
							gc = d;
							break inner;
						}
					}
				}
			}
			
			 if(Double.isNaN(gc) && getGC && false){
				String chr = current[chrid[currentId]];
				if(hasChr[currentId]) chr = chr.substring(3);
				int pos = 
					snpid1[currentId] <0 ? Integer.parseInt(current[snpid[this.currentId]]) : 
					(Integer.parseInt(current[snpid[this.currentId]]) +  Integer.parseInt(current[snpid1[currentId]]))/2;
				 Location chromosomeLoc = new Location("chromosome:"+chr+":"+(pos-150)+"-"+(pos+150));
					
				  SequenceAdaptor se = example.coreDriver.getSequenceAdaptor();
				  Sequence seq =  se.fetch(chromosomeLoc);
				  String str1 =seq.getString();
				  int gcc =0;
				  int at=0;
				  for(int i=0; i<str1.length(); i++){
					  char ch = str1.charAt(i);
					  if(ch=='G' || ch=='C') gcc++;
					  else if(ch=='A' || ch=='T') at++;
				  }
				   gc = ((double)gcc)/ ((double)gcc + (double)at);
				   System.err.println("getting GC for "+Arrays.asList(current)+" "+f1[currentId].getName());
			}
			for(int k=0; k<len; k++){
				if(compPrev[k]==0){
					print(out[k], currString[k], k==0  || !getGC? null : gc,k);
				}
			}
			print(this.ref_new, this.current,currentId==0 || !getGC ? null : gc,this.currentId);
		}
		this.close();
		warning.close();
	}
	
	private void print(PrintWriter pw, String[] strings, Double double1, int id) {
		if(!hasChr[id]) pw.print("chr");
		pw.print(strings[this.chrid[id]]+"\t");
		pw.print(strings[this.snpid[id]]+"\t");
	
		pw.print((snpid1[id]>=0? strings[this.snpid1[id]] : (40+Integer.parseInt(strings[snpid[id]])))+"\t");
		pw.print(strings[this.snpid2[id]]);
		if(remid[id].size()>0){
			for(int k=0; k<remid[id].size(); k++){
				pw.print("\t"+strings[remid[id].get(k)]);
			}
			if(!remid[id].contains(strings.length-1)){
				pw.print("\t"+strings[strings.length-1]);
			}
		}
		pw.println(double1!=null ? "\t"+String.format("%5.3g", double1) :"");
		pw.flush();
	}
	final int[] snpid2;
	final List<Integer>[]remid;

	/*Comparator<String[]> compara = new Comparator<String[]>(){

		@Override
		public int compare(String[] o1, String[] o2) {
			int k = o1[chrid].compareTo(o2[chrid]);
			if(k==0){
				int i1 = Integer.parseInt(o1[snpid]);
				int i2 = Integer.parseInt(o2[snpid]);
				if(i1==i2) return 0;
				else return i1< i2 ? -1 : 1;
			}
			return k < 0 ? -1 : 1;
		}
		
	};*/
	final File[] f1;
	public int compare(String[][] str, int i1, int i2){
		String[] o1 = str[i1];
		String[] o2 = str[i2];
		String chr1 = o1[chrid[i1]];
		String chr2 = o2[chrid[i2]];
		if(hasChr[i1]) chr1 = chr1.substring(3).replace("Mt", "M");
		if(hasChr[i2]) chr2 = chr2.substring(3).replace("Mt", "M");
	
		int k = chr1.compareTo(chr2);
		if(k==0){
			int ii1 = Integer.parseInt(o1[snpid[i1]]);
			int ii2 = Integer.parseInt(o2[snpid[i2]]);
			if(ii1==ii2){
				return o1[snpid2[i1]].compareTo(o2[snpid2[i2]]);
			}
			else return ii1< ii2 ? -1 : 1;
		}
		else{
			System.err.println("warning comparing diff chroms "+chr1+" "+chr2+" "+f1[i1].getName()+" "+f1[i2].getName());
		}
		return k < 0 ? -1 : 1;
	}
	
	
	final PrintWriter warning;
	public boolean readLine() throws Exception{
	   int newid =-1;
		for(int k=0; k<len; k++){
		   if(compPrev[k]<=0){
			   String str = in[k].readLine();
			   if(str!=null){
				   currString[k] = str.split(sep);
			   }
			   else currString[k]=  null;
		   }
		   if(currString[k]!=null && ( newid<0 || compare(currString, k, newid)<0)){
				   newid =k;
		   }
	   }
		if(newid<0) return false;
		current = currString[newid];
		currentId = newid;
		String id = current[this.snpid2[newid]];
		for(int k=0; k<len; k++){
			if(k==newid) compPrev[k] =0;
			else{
				compPrev[k] = 
			
				currString[k] == null ? +1:
				compare(currString,k,newid);
				if(done[k].contains(id)){
					warning.println(" not using same build file - prob with "+id+" "+f1[k]+" vs "+f1[newid]);
				}
			}
			//if(!id.equals)
		}
		
		done[newid].add(id);
		if(done[newid].size()>500) done[newid].removeFirst();
		
		return true;
	}
	
	final LinkedList<String>[] done ;
	int currentId;
	
}
