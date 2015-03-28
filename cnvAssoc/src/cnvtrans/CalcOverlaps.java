package cnvtrans;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

/* Program used for analysis in BGI paper on SVs from denovo assembly */
public class CalcOverlaps {

	static int overlThresh =0;
	  static int minSnpLength = 2;
	  //static int cnv_seq_rest = 2;
	  
	  static int minCNVLength = 0;
	  static int maxCNVLength = 900*1000*1000;
	  double includeFlanking = 0.5;
	public static int inclSnpThresh = 0;  //how far to look left or right for snp
	
	Map<String, IndivCNV> map;
	Map<String, IndivCNV> flanking;
	Map<String, SortedSet<Integer>> snps;
	String nme;
	TreeSet<Integer> s;
	boolean expand = true;
	
	public CalcOverlaps(String string2, int[] is, boolean b, String string)  throws Exception{
	 if(string!=null) snps = readSNPS(new File(string));
		TreeSet<Integer> s = new TreeSet<Integer>();
		this.map = CalcOverlaps.read(new File(string2),s,is,b);
		this.nme = string2; 
		if(includeFlanking>0 && string!=null) calcFlanking(snps,includeFlanking);
		//this.flanking = null;
	}
	public CalcOverlaps(String f2, int[] is, boolean b,Map<String, SortedSet<Integer>> snps)  throws Exception{
			TreeSet<Integer> s = new TreeSet<Integer>();
			this.snps = snps;
			this.map = CalcOverlaps.read(new File(f2),s,is,b);
			this.nme = f2;
			if(includeFlanking> 0 && snps!=null) calcFlanking(snps,includeFlanking);
			//this.flanking = null;
		}
	
	private Map<String, IndivCNV> calcFlanking(Map<String, SortedSet<Integer>> snps, double expand) {
		Map<String, IndivCNV> res = new HashMap<String,IndivCNV>();
		for(Iterator<String> it = map.keySet().iterator(); it.hasNext();){
			String key = it.next();
			res.put(key, map.get(key).getFlanking(snps.get(key),expand));
		}
		return res;
		
	}

	public CalcOverlaps(Map<String, IndivCNV> m, String string) {
	this.map = m;
	this.nme = string;
	}

	public static void main(String[] args){
		try{
			String[] str = args[0].split(":");
		
		
			File[] f = new File[str.length];
			f[0] = new File(str[0]);
			CalcOverlaps cnvs_cgh = new CalcOverlaps(args[2],  new int[] {0,2,3,4,5},false,"build36_cgh.txt");
			
			CalcOverlaps cnvs = new CalcOverlaps(str[0],  new int[] {7,6,5,4,2,1},true,"build36.gc.txt");
			cnvs.nme  = "Illumina1M";
			CalcOverlaps cnvs_merged = cnvs;
			CalcOverlaps cnvs_sequence = new CalcOverlaps(args[1],  new int[] {4,5,6,-1,1},false,(String)null);
			Map<String, Map<String,String>>  topr = new TreeMap<String,Map<String,String>>();
			
		    cnvs_sequence.nme="Assembly";
		    cnvs_cgh.nme = "CGH";
		    
		    findOverlapsBoth(cnvs_cgh,cnvs_sequence,topr);
			/*for(int k=1; k<str.length; k++){
				CalcOverlaps cnvs_ = new CalcOverlaps(str[k],  new int[] {7,6,5,4,2,1},true, cnvs.snps);
				cnvs_.snps = cnvs.snps;
				 cnvs_merged = findOverlapsBoth(cnvs_merged,cnvs_,topr);
			}
			cnvs_merged.nme = "1M_predictions";
			/*
			findOverlapsBoth(cnvs,cnvs_cgh,topr);
			 */
			
			CalcOverlaps mergedSeqCGH = findOverlapsBoth(cnvs_merged,cnvs_cgh,topr);
			mergedSeqCGH.nme = "Combined 1M/CGH";
			
//			System.err.println("cgh with seq");
			
	//		System.err.println("1M with seq");
			findOverlapsBoth(cnvs,cnvs_sequence,topr);
			findOverlapsBoth(cnvs_merged,cnvs_sequence,topr);
		//	System.err.println("all with seq");
			findOverlapsBoth(mergedSeqCGH, cnvs_sequence,topr);
		  
			for(Iterator<Map<String,String>> it = topr.values().iterator(); it.hasNext();){
				for(Iterator<String> it1 = it.next().values().iterator(); it1.hasNext();){
					System.err.println(it1.next());
				}
					
			}
		  
		   //findOverlapsBoth(cnvs,cnvs3);
		  
		//	
		  //}
		//	System.err.println("h");
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	private Map<String, IndivCNV> restrict(Map<String, SortedSet<Integer>> snps, int num) {
		Map<String, IndivCNV> res = new HashMap<String, IndivCNV>();
		for(Iterator<String> it = map.keySet().iterator(); it.hasNext();){
			String key = it.next();
			res.put(key, map.get(key).annontateNum(snps.get(key),num));
		}
		return res;
		
	}

	
	

	public static Map<String,SortedSet<Integer>> readSNPS(File f) throws Exception{
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String st = "";
		boolean shorten = f.getName().indexOf("cgh")<0;
		int[] ids = shorten ? new int[] {0,1}: new int[] {1,2};
		Map<String, SortedSet<Integer>> m = new HashMap<String, SortedSet<Integer>>();
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			String chr = shorten ? str[ids[0]].substring(3) : str[ids[0]];
			int pos = Integer.parseInt(str[ids[1]]);
			SortedSet<Integer> l = m.get(chr);
			if(l==null) m.put(chr, l = new TreeSet<Integer>());
			l.add(pos);
		}
		return m;
	}
	
	
	public static SortedMap<String, IndivCNV> read(File f, Set<Integer> set, int[] inds, boolean reverse) throws Exception{
		int maxcnv =0;
		SortedMap<String, IndivCNV> m = new TreeMap();
		Iterator<String> iter = Assoc.getBR(new File[] {f});
		String[] header =iter.next().split("\\s+");
		String st ="";
		int maxSize =0;
		int cnt=0;
		while(iter.hasNext()){
			st=iter.next();
			String [] str = st.split("\\s+");
			CNV cnv = new CNV(str,inds,reverse);	
			//String indiv = str[0];
			IndivCNV cnvs = m.get(cnv.chrom);
			if(cnvs == null){
				cnvs = new IndivCNV(cnv.chrom,0);
				m.put(cnv.chrom, cnvs);
			}
			//CNV cnv = new CNV(str,inds,reverse);
			int cnvmax = cnv.type.max();
			if(cnvmax>maxcnv) maxcnv=cnvmax;
			if(cnv.length()>maxSize){
				maxSize = cnv.length();
			//	System.err.println("max size "+maxSize);
			}
			if(cnv.nosnp>=minSnpLength
					
					//&& cnv.length()>=minCNVLength
					// && cnv.length() <=maxCNVLength
					//&& cnv.type<2
					){
				/*if(cnvs.cnvs.contains(cnv)){
					//CNV cnv1 = cnvs.cnvs.headSet(cnv).last();
					CNV cnv2 = cnvs.cnvs.tailSet(cnv).first();
					
					System.err.println(cnv.start+"-"+cnv.end);
					//System.err.println(cnv1.start+"-"+cnv1.end);
					System.err.println(cnv2.start+"-"+cnv2.end);
					throw new RuntimeException("!!");
				}*/
				cnvs.add(cnv);
				cnt++;
			}
			set.add(cnv.start);
			set.add(cnv.end);
		}
		/*for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String key = it.next();
			System.err.println(key+" "+m.get(key).cnvs.size());
		}*/
		return m;
	}
	
	public static CalcOverlaps findOverlapsBoth(CalcOverlaps m1, CalcOverlaps m2, Map<String, Map<String,String>> topr){
		
		CalcOverlaps merge =findOverlaps(m1,m2,topr);
		 findOverlaps(m2,m1,topr);
		/*if(m1.flanking!=null) {
			System.err.println(m1.nme+"_flanking which are found in "+m2.nme);
			findOverlaps(m1.flanking,m2.map,m1.nme+"_f_"+m2.nme);
			System.err.println(m2.nme+"_flanking which are found in "+m1.nme);
			findOverlaps(m2.map,m1.flanking,m2.nme+"_f_"+m1.nme);
		}
		if(m2.flanking!=null){
			System.err.println(m2.nme+"_flanking which are found in "+m1.nme);
		
			findOverlaps(m2.flanking,m1.map,m2.nme+"_f_"+m2.nme);
			System.err.println(m1.nme+"_flanking which are found in "+m2.nme);
			findOverlaps(m1.map,m2.flanking,m1.nme+"_f_"+m2.nme);
		}*/
		return merge;
	}
	
	
	//how many m1 could be found in m2
	public static CalcOverlaps findOverlaps(CalcOverlaps o1, CalcOverlaps o2,Map<String, Map<String,String>> topr){
		int[] res1 = new int[] {0,0};
		Map<String, IndivCNV> m1,m2;
		m2 = o2.map; 
		String nme = o1.nme+"_"+o2.nme;
		if(o2.snps!=null) {
			m1 = o1.restrict(o2.snps, CalcOverlaps.minSnpLength);
		}
		else m1 = o1.map;
		int cnt =0;
		Map<String, IndivCNV> m = new HashMap<String, IndivCNV>();
		for(Iterator<String> it = m1.keySet().iterator(); it.hasNext();){
			String key = it.next();
			IndivCNV cnv1 = m1.get(key);
			cnt+=cnv1.cnvs.size();
			IndivCNV cnv2 = m2.get(key);
			if(cnv2!=null){
			 m.put(key, findOverlaps(cnv1,cnv2,res1));
			}
		}
		
		Map<String, String> topr1 = topr.get(o1.nme);
		if(topr1==null) topr.put(o1.nme, topr1  = new TreeMap<String, String>());
		topr1.put(o2.nme,o1.nme+","+o2.nme+","+res1[1]+","+res1[0]+","+String.format("%5.2g",(double)res1[0]/(double)res1[1]));
		CalcOverlaps res =  new CalcOverlaps(m,"inters_"+nme);
		res.snps = merge(o2.snps,o1.snps);
		return res;
	}
	private static Map<String, SortedSet<Integer>> merge(
			Map<String, SortedSet<Integer>> snps1,
			Map<String, SortedSet<Integer>> snps2) {
		Map<String, SortedSet<Integer>> res = new HashMap<String, SortedSet<Integer>>();
		Set<String> keys  = new HashSet<String>();
		if(snps1!=null)	keys.addAll(snps1.keySet());
		if(snps2!=null) keys.addAll(snps2.keySet());
		for(Iterator<String> it = keys.iterator(); it.hasNext();){
			String key = it.next();
			SortedSet<Integer> s1 = snps1==null || !snps1.containsKey(key) ? new TreeSet<Integer>() : snps1.get(key);
			SortedSet<Integer> s2 = snps2==null || !snps2.containsKey(key) ? new TreeSet<Integer>() : snps2.get(key);
			SortedSet<Integer> s3 = new TreeSet<Integer>(s1);
			s3.addAll(s2);
			res.put(key, s3);
		}
		return res;
	}
	//how many cnv1 could be found in cnv2
	public static IndivCNV  findOverlaps(IndivCNV cnv1, IndivCNV cnv2, int[] overl){
		int overlap =0;
		int cnt=0;
		IndivCNV indiv = new IndivCNV(cnv1.id, cnv1.index);
		for(Iterator<CNV> it = cnv1.cnvs.iterator(); it.hasNext();){
			CNV cnv = it.next();
			CNV cnv_2 = cnv2.hasOverlap(cnv);
			if(cnv_2!=null){
				overlap++;
				indiv.add(cnv.merge(cnv_2));
			}
			cnt++;
		}
		overl[0]+=overlap;
		overl[1]+=cnt;
		return indiv;
	}
	//}
	//double[][]
	
}
