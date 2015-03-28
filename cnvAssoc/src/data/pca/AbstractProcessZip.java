package data.pca;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.zip.GZIPOutputStream;

import lc1.util.ApacheCompressor;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipFile;

import conversion.Utils;

public abstract class AbstractProcessZip {
	 PrintWriter missingpw;
	static String limit = null;
	static boolean compress = false;
	static int limitDone = Integer.MAX_VALUE-1;
	protected static String regionsF = "regions_all.txt";
	protected static String buildF = "build37.txt";
	protected static boolean singleEntry = true;
	protected static boolean overwrite = true;
	protected String chr_out;
	boolean calcWeights = false;
	 int[][] alias,aliasRev;
	 boolean[][] include;
	 File[] outf;//, outf1;
		PrintWriter[] out;
	protected static class Location implements Comparable{
			String chrom;
			String chrom_out;
			int start;
			int end;
			boolean all = false;
			final String locString;
			public String toString(){
				return locString;
			}
			public Location(String[] str, String chr_out, boolean all) {
				this.chrom = str[0].substring(3);
				this.chrom_out = chr_out;
				this.all = all;
				this.start = Integer.parseInt(str[1]);
				this.end = Integer.parseInt(str[2]);
				locString = chrom+"_"+start+"_"+end;
				}
			@Override
			public int compareTo(Object o) {
				Location l1 = (Location)o;
				if(l1.chrom.equals(this.chrom)){
					if(l1.start == this.start) return 0;
					else return start < l1.start ? -1 : 1;
				}
				else return l1.chrom.compareTo(chrom);
			}
			public boolean contains(String chrom, int pos) {
				boolean res =  pos>=this.start && pos <=this.end && (all || this.chrom.equals(chrom));
				if(res){
					return true;
				}else{
					return false;
				}
			}
			
			public String line() {
				return "chr"+chrom_out+"\t"+start+"\t"+end+"\t"+locString;
			}
			int cnt=0;
			public void addCount() {
			cnt++;
				
			}
		}
	int[] offset = null;
	List<Double> depth = new ArrayList<Double>();
	int depth_index=-1;
	protected final double thresh;
	protected static boolean keepZipOpen = true;
	static boolean randomOffset = false;
	static boolean checkOrder = false;
	static double threshNA = 0.1;
	protected static double var_thresh = 0.25;
	protected static int maxIndiv = 1000*1000;
	protected static boolean tocache = true;
	Set<String> toexcl = new HashSet<String>();
	protected final int thin_min;
	protected final int thin_mult;
	protected final int thin_start;
	protected final boolean[] readFile;
	protected ZipFile[] zf = null;
	protected ZipArchiveEntry[] ze = null;
	protected File[] f1 = null;
	protected ProcessSNP[] l1 = null;
	protected boolean[] type;
	protected File[] in;
	protected final File[] dir;
	String chr1 = "";
	protected List<String>[] sample_id;
	protected final String[] lrrst;
	int[][] lrr_id = null;
	Boolean[][] gl = null;
	protected Location currentRegion = new Location(new String[] {"all",0+"",9*1000*1000*1000+""},"all",true);
	protected SortedSet<Location> regions = null;
	Map<String, ZipFile[]> m = new HashMap<String, ZipFile[]>();
	protected final SortedMap<String, String> chroms;
	public Map<String, Integer> karyo;
	protected int[] sze;

	public void readKaryotypes(File dirp) throws IOException {
		
			
			File[] karyofiles = dirp.listFiles(new FileFilter(){
	
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().indexOf("karyo")>=0;
				}
				
			});
			if(karyofiles.length>0){
				karyo = new HashMap<String, Integer>();
			BufferedReader br = new BufferedReader(new FileReader(karyofiles[0]));
			String st = "";
			while((st = br.readLine())!=null){
				String[] str = st.split("\t");
				karyo.put(str[0], Integer.parseInt(str[1]));
			}
			br.close();
		}	
		
		
	}

	protected String skip(String subdir, String string, String suff, File[] f1) {
		String nme = string;
		for(int k=0; k<dir.length; k++){
			File f = new File(dir[k],subdir+"/"+string+".zip");
			//note next line will probably break other applications of this code!!! IT is a quick hack
			f = new File(dir[k],string+".zip");
			if(!f.exists() && this.karyo!=null){
				nme = string+suff;
				f = new File(dir[k],nme+".zip");
			}
			if(!f.exists() || f.length()==0) {
				return null;
			}
			else {
				f1[k] = f;
			}
		}
		return nme;
	}

	protected void finish() {
		try{
			this.missingpw.close();
		for(Iterator<ZipFile[]> it = this.m.values().iterator(); it.hasNext();){
			ZipFile[] nxt = it.next();
			for(int k=0; k<nxt.length; k++){
				nxt[k].close();
			}
		}
		for(int i=0; i<out.length; i++){
			out[i].close();
		}
		/*if(snp_list.size()>0){
			PrintWriter pw_snps = new PrintWriter(new FileWriter(new File("snp_weights.txt")));
			for(int i=0; i<this.snp_list.size(); i++){
				pw_snps.println(snp_list.get(i));
				
			}
			pw_snps.close();
		}*/
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}


   abstract void makeRes();
  
	
	public void reInitialise(String nme1) throws Exception {
		
		zf = this.getZip(nme1, f1);
		
		ze = new ZipArchiveEntry[zf.length];
		
		
		if(lrr_id==null){
			lrr_id = new int[dir.length][lrrst.length];
			gl = new Boolean[dir.length][lrrst.length];
			//sample_id = new List[dir.length];
		
			for(int jk=0; jk < dir.length; jk++){
			
			  List<String > l = Arrays.asList(ApacheCompressor.readZipFrom(zf[jk], "Name").get(0).split("\\t+"));
			  for(int i=0; i<l.size(); i++){
				  l.set(i, l.get(i).trim());
			  }
			  for(int k=0; k<lrr_id[jk].length; k++){
				lrr_id[jk][k] = l.indexOf(lrrst[k]);
				gl[jk][k] = lrrst[k].equals("GL");
				if(lrr_id[jk][k] <0) {
					throw new RuntimeException(l+"\n"+lrrst[k]);
				}
			}
			}
			if(l1!=null){
			for(int k=0; k<lrr_id[0].length; k++){
				if(l1[k]!=null) if(l1[k]!=null) l1[k].setSize(sze);
			}
			}
			
			//l1_sex.setSize(sze);
			makeRes();
		}
		
		
		
	}

	protected ZipFile[] getZip(String nme1, File[] file)  {
		if(!keepZipOpen && m.size()>0) m.clear();
		ZipFile[] zf = m.get(nme1);
		if(zf==null){
			try{
			System.err.println("opening "+Arrays.asList(file));
			m.put(nme1, zf = new ZipFile[file.length]);
			for(int k=0; k<file.length; k++){
				zf[k] = new ZipFile(file[k]);
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		return zf;
	}

	public abstract void makeHists(String nme, String[] lrrst);

	

	protected String[] getChroms(File pdir, String[] chr) {
		File f = new File(pdir,chr[0]+".zip");
		if(f.exists()) return chr;
		else{
			String[] split = chr[0].split("\\*");
			List<String> res = new ArrayList<String>();
			File[] f1 = pdir.listFiles(Utils.dirFilter("chr",""));
			
			for(int k=0; k<f1.length; k++) {
				res.addAll(Arrays.asList(f1[k].list(Utils.dirFilter(split[0], split[1]+".zip"))));
			}
			
			String[] resu =  res.toArray(new String[0]);
			for(int i=0; i<resu.length; i++) resu[i] = resu[i].replace(".zip", "");
			return resu;
		}
	}

	private int match(String[] chroms2, String name) {
		for(int k=0; k< chroms2.length; k++){
			if(chroms2[k].startsWith("pcs_"+name.substring(3))) {
				return k;
			}
		}
		return -1;
	}

	/*private Iterator<String[]> getBuildReaders() throws Exception {
	
		final Iterator<String>[] res = new Iterator[in.length];
		for(int i=0; i<res.length; i++){
			res[i] = getBuildReader(i);
		}
		return Util.getI
	}*/
	private Iterator<String> getBuildReader(int i) throws Exception {
		BufferedReader br = Utils.getBufferedReader(this.in[i]);
		Entry<String, String> ent = chroms.entrySet().iterator().next();
		if(br==null){
			
			File nxtF = new File(dir[0],ent.getKey()+".zip");
			File pdir = in[i].getParentFile();
			if(!nxtF.exists()){
				return Utils.getStringIterator(Utils.listFiles(pdir, chroms));
			}
		 else if(chroms.size()==1){
			 ZipFile[] zf1  = getZip(ent.getKey(), new File[] {nxtF});
				 br = ApacheCompressor.getBufferedReader(zf1[0],  "SNPS");
			}
		}
			return Utils.getStringIterator(br);
		
		
		
	}

	protected void filter() {
		try{
			Iterator<String> br = getBuildReader(0);
			
			for(int i=0; br.hasNext(); i++){
				String st = br.next();
				String[] str = st.split("\t");
				try{
			
				Location loc = new Location(str, this.chr_out,false);
				Location loc2 = new Location(str, this.chr_out,false);
				loc.start = loc.start+1;
				loc2.start = loc2.start - 50000;
				SortedSet<Location> ss = this.regions.headSet(loc);
				if(ss.size()>0){
					ss = ss.tailSet(loc2);
					
					inner: for(Iterator<Location > it = ss.iterator(); it.hasNext();){
						Location loc1 = it.next();
						if(loc1.contains(loc.chrom, loc.start)){
							loc1.addCount();
							break inner;
						}
					}
				}
				}catch(Exception exc){
					System.err.println("problem at line "+i);
					System.err.println("problem with "+Arrays.asList(str));
					exc.printStackTrace();
				}
			}
			for(Iterator<Location> it = this.regions.iterator(); it.hasNext();){
				Location loc = it.next();
				if(loc.cnt==0){
					System.err.println("removing "+loc);
					it.remove();
				}
			}
			int sze = regions.size();
			//br.close();
			//if(zf1!=null) zf1.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	protected File getFile(File file, String string) {
		File f1 = new File(file, string);
		if(!f1.exists()) f1 = new File(file, string+".gz");
		return f1;
	}

	int subdir_ind =-1;
	public void runInner(int index, int ii) throws Exception {
		int thin = (int)Math.floor(Math.max((double)this.thin_min, (double)this.thin_start/Math.pow((double)this.thin_mult,ii)));
		int thin_offset = randomOffset? (int) Math.floor(Math.min(thin-1, Math.random()*Math.min(thin, 10))) : 0;
		System.err.println("thin "+ii+" "+thin+" "+thin_offset);
	chr1 = "";
	boolean skip = false;
	int cnt=0;
	int cntDone=0;
	
	Integer centromerePos = Integer.MAX_VALUE;
	Iterator<String> br = this.getBuildReader(0);
	for(int kk=0; br.hasNext() && cntDone < limitDone   ; kk++){
	  try{
		  String st = br.next();
	  String[] str = st.split("\\s+");
	  String subdir = ".";
	  if(str.length==3){
		  String[] str1 = new String[4];
		  str1[0] = "chrpcs_out0";
		  str1[2] = str[2]; str1[1] = str[1];
		  str1[3] = str[0].substring(3)+"_"+str[1]+"_"+str[2]+"PC";
		  str = str1;
	  }
	  String chr = str[0].startsWith("chr") ? str[0].substring(3) : str[0];
	  if(subdir_ind>0 && str.length>subdir_ind){
		  subdir = str[subdir_ind];
		  chr = subdir.startsWith("chr") ? subdir.substring(3) : subdir;
	  }
	 /// System.err.println(kk+" "+this.l1[0].cnt);
	
	  
	  int pos = Integer.parseInt(str[1]);
	  if(!chr.equals(chr1) || (centromerePos!=null && pos > centromerePos)){
		if(this.chroms.keySet().contains(chr)){
			String suff = "p";
			centromerePos = karyo==null ? (Integer)(Integer.MAX_VALUE) : this.karyo.get(chr);
			if(centromerePos!=null && pos >  centromerePos){
				centromerePos= Integer.MAX_VALUE;
				suff = "q";
			}
		chr1 = chr;
		skip = false;
		
	//	File[] f1 = getFile(dir,chr+".zip");
		String nme =skip(subdir, chr,suff, f1); 
		if(nme==null) skip =true;
		else{
			
				reInitialise(nme);
				this.makeP();
		}
		}else skip =true;
	  }
	  if(skip) {
		  continue;
	  }
	//  int pos = Integer.parseInt(str[1]);
	  String id = str[3];
	  int rem = cnt %thin;
	  if(Math.abs(rem-thin_offset)<0.01 && (limit==null || id.indexOf(limit)>=0)){
	//  System.err.
		//  System.err.println(id);
		//  System.err.println(cnt+" "+thin_offset);
		  boolean complete = true;
		  boolean allNull = true;
		  boolean excl = false;
		
		  if(!currentRegion.contains(chr,pos) || toexcl.contains(id)){
			  complete=false;
			  excl = true;
		  }else{
			  complete = getEntries(zf,ze,id, readFile);
			  //if(complete){
				  for(int k=0; k<zf.length; k++){
					  if(readFile[k]) readData( id, complete,k);
					  //complete = complete && readFile[k];
		              allNull = allNull && !readFile[k];
				  }
			  //}
		  }
		if(allNull || excl) {
		//	System.err.println("did not find "+id);
			missingpw.println(id);
			//missingpw.flush();
			continue;
		}
		{
			
			try{
			 process(id, pos);
			}catch(Exception exc){
				exc.printStackTrace();
			}
			 cntDone++;
			
		  }
	}
	
	cnt++;
	 
	  }catch(Exception exc){
	  exc.printStackTrace();
		 // out[0].println(st+"\t"+Double.NaN);
	  }
	  out[0].close();
	}
	}

	
	private boolean getEntries(ZipFile[] zf2, ZipArchiveEntry[] ze2, String id, boolean[] read) {
		 // Arrays.fill(read, false);
		  boolean complete = true;
		for(int k=0; k<zf2.length; k++){
			ze2[k] = zf2[k].getEntry(id);
			read[k] = ze2[k]!=null ;
			complete = complete && read[k];
		}
		return complete;
	}

	public abstract void process(String id, int pos) ;


	abstract boolean readData(String id, boolean complete, int k);
	

	

	protected void makeP() {}

	public abstract void run() throws Exception;

	private void close(ZipFile[] zf) throws IOException {
		for(int k=0; k<zf.length; k++)zf[k].close();
		
	}

	public AbstractProcessZip(File pdir, String[] dir, String[] lrrst, String[] strtype, String[] chroms1_, 
			String[] thin, String thresh, boolean mergeRegions, Integer no_pcs) {
		this.lrrst = lrrst;
		try{
			File resultsDir = new File("results");
			resultsDir.mkdir();
			FileOutputStream  missingf = new FileOutputStream(new File(resultsDir,chroms1_[0]+"_missing.txt"+
					(compress ? ".gz":"")));
			this.missingpw =  new PrintWriter(new BufferedWriter(new OutputStreamWriter(
					compress  ? new GZIPOutputStream(missingf) : missingf
					)));
			}catch(Exception exc){
				exc.printStackTrace();
			}
		//this.start = Integer.parseInt(startend[0]);
		//this.end = Integer.parseInt(startend[1]);
		this.thresh = Double.parseDouble(thresh);
		File regionsFile = new File(regionsF);
		if(!regionsFile.exists() || lrrst.length>1) singleEntry = false;
		String[] chroms_ =  getChroms(new File(pdir,dir[0]), chroms1_);
		this.chroms = new TreeMap<String, String>();
	if(chroms.size()>10){
		System.err.println("not keeping zips in buffer, too many files "+chroms);
		keepZipOpen = false;
	}
		for(int k=0; k<chroms_.length; k++){
			String nxt = chroms_[k];
			String nxt1 = "chr"+nxt.replaceAll("pcs_", "").replaceAll("_out", "").split("_")[0];
			this.chroms.put(nxt,nxt1);
		}
		if(regionsFile.exists() && regionsFile.length()>0){
			calcWeights = false;
			this.regions = new TreeSet<Location>();
			try{
			BufferedReader br = new BufferedReader(new FileReader(regionsFile));
			String st = "";
			while((st = br.readLine())!=null){
				if(st.startsWith("#") || st.length()==0) continue;
				String[] str = st.split("\\s+");
				if(this.chroms.keySet().contains(str[0].substring(3))){
					Location loc = new Location(str,this.chr_out, false);
					if(regions.contains(loc)){
						System.err.println("duplicate "+loc+" "+regions.tailSet(loc).first());
					}
				this.regions.add(loc);
				}
			}
			}catch(Exception exc){exc.printStackTrace();}
		}
		this.type = new boolean[lrrst.length];
		this.thin_start = Integer.parseInt(thin[0]);
		this.thin_min =Math.max(1, Integer.parseInt(thin[1]));
		this.thin_mult = Math.max(1, Integer.parseInt(thin[2]));
		for(int k=0; k<type.length; k++){
			type[k] = (strtype[k].indexOf("str")>=0);
		}
		
		
		
		
		this.dir = new File[dir.length];
		
		this.f1 = new File[dir.length];
		sample_id = new List[dir.length];
		
		this.include = new boolean[dir.length][];
		sze = new int[dir.length];
		this.readFile = new boolean[dir.length];
		this.alias = new int[dir.length][];
		this.aliasRev = new int[dir.length][];
		offset = new int[dir.length];
	
		
		
		int totsze =0;
		
		try{
			
			if(chroms.size()==0) throw new RuntimeException("chroms size is 0");
			Entry<String, String> ent = this.chroms.entrySet().iterator().next();
			for(int k=0; k<dir.length; k++){
				String subdir = ".";
				
				this.dir[k] = new File(pdir, dir[k]);
				File pdir1 = (new File(dir[0]));
				File f2 = new File(pdir1,ent.getKey()+".zip");
				if(!f2.exists()){
					File pdir2 = new File(pdir1, ent.getValue());
					subdir = ent.getValue();
				}
				f1[k] = new File(dir[k], subdir+"/"+ent.getKey()+".zip");
			}
			zf = this.getZip(ent.getKey(), f1);
			
			for(int k=0; k<dir.length; k++){
			    BufferedReader br2;
				ZipFile zf = this.zf[k];
				processSampleInfo(zf,k);
				 // if(zf!=null)zf.close();
				  //now read in previously defined pcs
				  
				
				  
			} 
		
			 this.readKaryotypes(pdir);
		// in =getFile(this.dir[0], OptionBuild.build+".gc.txt");
		 in = getFile(this.dir, buildF);
		if(regions!=null){
			filter();
			if(mergeRegions) regions = null;
		}
		
		 String nme = in[0].getName();
		 nme = nme.substring(0,nme.length()-4);
		makeHists(nme, lrrst);
		outf = new File[dir.length];
		out = new PrintWriter[dir.length];
		for(int k=0; k<dir.length; k++){
	     outf[k] = new File(dir[k], nme+".var.txt");
	     out[k] = new PrintWriter(new BufferedWriter(new FileWriter(outf[k],false)));
		}
		 	  
	   
	    //  BufferedReader br1 =  conversion.Utils.getBufferedReader(outf);
	      
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
	}

	private File[] getFile(File[] dir2, String buildF2) {
		File[] res = new File[dir2.length];
		for(int k=0; k<res.length; k++)res[k] = getFile(dir2[k], buildF2);
		return res;
	}

	void processSampleInfo(ZipFile zf, int k)  throws Exception{
		File var = new File(dir[k], "variance.txt");
		int var_index = -1;
	//	int depth_index = -1;
		BufferedReader br2;
		if(var.exists()){
	    	br2 = new BufferedReader(new FileReader(var));
		br2.readLine(); 
		var_index = 1;
		
		}
		else{
			
			 br2 = ApacheCompressor.getBufferedReader(zf,  "Samples");
			 BufferedReader br1 = ApacheCompressor.getBufferedReader(zf, "Name");
			 br1.readLine();br1.readLine();
			 List<String> l = Arrays.asList(br1.readLine().split("\\s+"));
			 var_index = l.indexOf("variance");
			 depth_index = l.indexOf("avgdepth");
			 br1.close();
			 
		}
		String st = "";
		this.sample_id[k] = new ArrayList<String>();
		List<Integer> toincl = new ArrayList<Integer>();
		int kj =0;
		List<Double> vars = new ArrayList<Double>();
		{
			List<String> sample_info = new ArrayList<String>(); 
			for(kj=0;(st = br2.readLine())!=null; kj++){
				String[] str = st.split("\t");
				double v=  var_index<0 ? 0 : Double.parseDouble(str[var_index]);
				depth.add( depth_index<0 ? 1 : Double.parseDouble(str[depth_index]));
				vars.add(v);
				sample_info.add(str[0]+"\t"+v+"\n");
				if(v < var_thresh && sample_id[k].size() < maxIndiv){
					this.sample_id[k].add(str[0]);
					toincl.add(kj);
					
				}
			}
			saveSampleInfo(sample_info,k);
		}
		
		alias[k] = new int[toincl.size()];
		include[k] = new boolean[kj];
		aliasRev[k] = new int[kj];
		Arrays.fill(include[k], false);
		Arrays.fill(aliasRev[k], -1);
		for(int kk=0; kk<toincl.size(); kk++){
			alias[k][kk] = toincl.get(kk);
		    include[k][toincl.get(kk)] = true;
		    aliasRev[k][toincl.get(kk)] = kk;
		    
		}
	
		System.err.println("before "+aliasRev.length+" after "+alias.length);
		 sze[k] = sample_id[k].size();
		  if(k>0) offset[k-1] = sze[k];
		
		//totsze+=sze[k];
			
	
		  
	}

	 abstract void saveSampleInfo(List<String> sample_info, int k)  ;
	
}