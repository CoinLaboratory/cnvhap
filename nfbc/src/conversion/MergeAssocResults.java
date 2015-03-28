package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

public class MergeAssocResults {

	public static void main(String[] args){
		try{
			File dir = new File(".");
			String prefix = args[0];
			String suffix = args[1];
			File dirr = new File(dir, "result");
			dirr.mkdir();
			File res = new File(dirr, suffix);
		MergeAssocResults mar = new MergeAssocResults(dir, prefix, suffix, res);
			mar.read();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	String suffix;
	
	String prefix;
	
	File[] in;
	BufferedReader[] br;
	
	Map<Integer, Integer> best = new HashMap<Integer, Integer>();
	
	final PrintWriter pw;
	int chr_ind, pos_ind;
	
	public MergeAssocResults(File dir, final String prefix, final String suffix, File res) throws Exception{
		File[] f = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && pathname.getName().startsWith(prefix);
			}
			
		});
		FileFilter ff = new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.getName().startsWith("result") && pathname.getName().endsWith(suffix);
			}
			
		};
		this.in = new File[f.length];
		String st = "";
		this.br = new BufferedReader[f.length];
		for(int k=0; k<f.length; k++){
		File avg = new File(f[k], "avg");
//			System.err.println(avg.getAbsolutePath()+" "+Arrays.asList(avg.list()));
			File[] f1 = avg.listFiles(ff);
			if(f1.length==0){
				System.err.println("no file in "+avg.getAbsolutePath());
			}
			else{
			in[k] = f1[0];
			this.br[k] = new BufferedReader(new FileReader(this.in[k]));
			st = br[k].readLine();
			}
		}
		List<String> str = Arrays.asList(st.split("\t"));
		this.chr_ind = str.indexOf("chrom");
		this.pos_ind = str.indexOf("loc");
		this.pw = new PrintWriter(new BufferedWriter(new FileWriter(res)));
		pw.println(st);
	}
	
	class Snp implements Comparable{
		String st;
		int mindist;
		//int start,  end;
		//String name;
		public Snp(File f, String st, String[] str) {
			this.st = st;
			this.mindist = Integer.parseInt(str[pos_ind]);
		//	this.name = f.getName();
		}
		
		public void setEnd(int st, int end) {
			if(mindist < st || mindist > end){
				throw new RuntimeException("!!");
			}
			//this.start = st;
			//this.end = end;
			this.mindist = Math.min(mindist-st , end - mindist);
			
		}

		@Override
		public int compareTo(Object o) {
			int d1 = ((Snp)o).mindist;
			if(mindist==d1) return 0;
			return (mindist < d1) ? 1 : -1;
		}
		public String toString(){
			return //start+"\t"+end+"\t"+
			st;
		}
	}
	
	Map<String, List<Snp>> m = new TreeMap<String, List<Snp>>();
	
	public void read() throws Exception{
		List<String> snpid[] = new List[br.length];
		for(int k=0; k<br.length; k++){
			if(in[k]!=null){
				snpid[k] = new ArrayList<String>();
				read(k);
			}
			
		}
		for(Iterator<Entry<String, List<Snp>>> ents = m.entrySet().iterator(); ents.hasNext();){
			Entry<String, List<Snp>> ent = ents.next();
			List<Snp> vals = ent.getValue();
			Collections.sort(vals);
			Snp snp = vals.get(0);
			pw.println(snp.st);
		}
		pw.close();
		
	}
	
	public void read(int k) throws Exception{
		String st = "";
		Integer start = null;
		List<Snp> snps_ = new ArrayList<Snp>();
		String[] str=null;
		while((st  = br[k].readLine())!=null){
			
			if(st.startsWith("NaN")) continue;
			str = st.split("\t");
			if(start==null){
				start = Integer.parseInt(str[pos_ind]);
			}
			String nme = str[chr_ind]+"_"+str[pos_ind];
			List<Snp> snps = m.get(nme);
			if(snps==null){
				m.put(nme, snps = new ArrayList<Snp>());
			}
			Snp snp = new Snp(this.in[k].getParentFile().getParentFile(), st, str);
			snps.add(snp);
			snps_.add(snp);
		}
		br[k].close();
		if(str!=null){
			Integer end = Integer.parseInt(str[pos_ind]);
		for(int k1=0; k1<snps_.size(); k1++){
			snps_.get(k1).setEnd(start, end);
		}
		}
	}
	
}
