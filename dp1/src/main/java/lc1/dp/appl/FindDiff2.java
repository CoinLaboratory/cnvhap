package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class FindDiff2 {

	
	public static void main(String[] args){
		try{
			String comp1 = "1M-200910161213:1M-200910211125";
			String types = "0.1.3.4.5.6.";
				//"1M_hap-20091010.txt:1M-20091010_1.txt";
			//1M-20090909_
			//"244k_agilent1-20090804.txt";//"1M_PennCNV.txt";
//			"1M-20090804:1M_PennCNV-20090804:1M-20091009:1M_hap-20091010:1M_20091009_:1M-20091010:1M-20091010_1";
			String ref =  "244k_agilent1-20090804";//"244k-2009101318";
			String pn = "negative";//"positive";
			String[] posneg = pn.split(":");
				//"244k_agilent1-20090804";
			String[] args1 = comp1.split(":");
			//for(int k=0; k<args1.length; k++){
			int k = 0;
				for(int k1=k+1; k1<args1.length; k1++){
			boolean[] tf = new boolean[]{true, false};
			String[] type = types.split(":");
			for(int i=0; i<tf.length; i++){
				for(int j=0; j<type.length; j++){
					for(int k2=0; k2<posneg.length; k2++){
				FindDiff2 fd = new FindDiff2(tf[i], type[j],ref,args1[k], args1[k1], posneg[k2]);
				fd.run();
					}
				}
			//}
			}
			}
		
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	BufferedReader[] br= new BufferedReader[2];
	Map<Locs, Locs>[]map = new Map[2];
	PrintWriter[] out = new PrintWriter[2];
	
	String comp1 = null;//"1M-20090804.txt";//"244k_agilent1-20090804.txt";//"1M_PennCNV.txt";
	String comp2 = null;//"1M-20090928.txt";//"1M_cnvHap.txt";
	FindDiff2(boolean all, String type,
			String ref,
			String comp1, String comp2, String posneg) throws Exception{
		String base = "roc_base_";
		this.comp1 = comp1;
		this.posneg = posneg;
		this.comp2  = comp2;
		base = base +(all ? "all_" : "max_");
	base = base +type;//el ? "0.1." :"3.4.5.6.");
		base = base +"_index_"+ref+"_line_";
		String[] base1 = new String[] { base+comp1,base+comp2};
	 

		
		for(int i=0; i<br.length; i++){
			File f = new File(base1[i]);
			System.err.println(f);
			diff_l[i] = new TreeSet<Locs>();
			this.br[i] = new BufferedReader(new FileReader(base1[i]+".txt"));
			map[i] = new TreeMap<Locs, Locs>();
			File outf = new File(base1[i]+"."+posneg);
			System.err.println(outf.getAbsolutePath());
			this.out[i] = new PrintWriter(new BufferedWriter(new FileWriter(outf)));
		}
		
		
		
		
	}
	
	static class Locs implements Comparable{
		String chr="";
		int start;
		int end;
		int no;
		String id;
		public Locs( String str, boolean pos){
			String[] st = str.split("_");
			
			if(pos){
				chr = st[0].split("\\.")[0];
				id = st[1];
			start = Integer.parseInt(st[2]);
			end = Integer.parseInt(st[3]);
			no = Integer.parseInt(st[4]);
			}
			else{
				id = st[1];
				start =  Integer.parseInt(st[0]);
				end = start;
			}
		}
		@Override
		public int compareTo(Object o) {
			if(this==o) return 0;
			Locs loc1 = (Locs)o;
			int c1 = chr.compareTo(loc1.chr);
			int c2 = id.compareTo(loc1.id);
			if(c1!=0) return c1;
			if(c2!=0) return c2;
			if(start !=loc1.start) return start < loc1.start ?  -1 :1;
			if(end !=loc1.end) return end < loc1.end ? -1 :1;
			
			return 0;
		}
		public String toString(){
			return "--mid "+chr+":"+start+":"+end+"  #"+id;
		}
	}
	
	static Comparator<Locs> compa1 = new Comparator<Locs>(){

		@Override
		public int compare(Locs o1, Locs loc1) {
			int c1 = o1.chr.compareTo(loc1.chr);
			if(c1!=0) return c1;
			if(o1.start !=loc1.start) return o1.start < loc1.start ?  -1 :1;
			if(o1.start !=loc1.start) return o1.end < loc1.end ? -1 :1;
			return 0;
		}
		
	};
	public void run() throws Exception{
		read(0); read(1);
		findDiff(0); findDiff(1);
	}
	
	Set<Locs>[] diff_l = new Set[2];
	private void findDiff(int comp) {
		for(Iterator<Locs> it = this.map[comp].keySet().iterator(); it.hasNext();){
			Locs nxt = it.next();
			if(!map[1-comp].containsKey(nxt)){
				//String[] str = map[comp].get(nxt).split("\\s+")[3].split("_");
				//String
				this.out[comp].println(map[comp].get(nxt)+"\t--restrictKb 50kb:50kb");
			}
			else{
				System.err.println("both "+nxt);
			}
		}
		out[comp].close();
	}
	final String posneg;
	public void read(int i) throws Exception{
		String st = br[i].readLine();
		System.err.println(st);
		while((st = br[i].readLine())!=null ){
		//	System.err.println(st);
			st = st.trim();
			if(st.startsWith(posneg)){
				String[] str = st.split("\\s+");
				double pp = Double.parseDouble(str[2]);
				if(pp>0.5){
					Locs loc = new Locs(str[3], this.posneg.equals("positive"));
					map[i].put(loc,loc);
				}
			}
		}
		br[i].close();
	}
	
	
}
