package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class ConvertFinalReport {
	
	static  boolean writeToFile = true;
	static final int countmax = 9999999;
	
	static boolean matches(Pattern p, String st){
		Matcher m = toexcludefinal.matcher(st);
		m.useAnchoringBounds(true);
		return m.find();
	}
	
	public class ReadRun  {
		final File input;
		BufferedReader bread;
		String st = "";
		int offset=0;
		int target = 0;
		int start =0;
		int cnt;
		String sample;
		public ReadRun(File input, BufferedReader br1, int cnt){
			this.bread = br1;
			this.input = input;
			this.cnt = cnt;
		}
		void setTarget(int pos){
			this.target = pos;
		}
		boolean active(){
			return st!=null;
		}
		public void runInner() throws Exception{
			//if(!active()) return;
		
			for(int k=start; k<=target+offset; k++){
			//	bread.
				try{
				st = bread.readLine();
				}catch(java.lang.OutOfMemoryError exc){
					throw new RuntimeException(this.input.getName()+":  out of memory at line "+ cnt);
				}
				cnt++;
				if(st==null) return;
				if(sample==null){
					
					sample = st.split(splitstr)[reports[0].sampleid];
					System.err.print("new sample :");
					System.err.println(sample+ " "+cnt+ " "+start+ " "+offset);
				}
			//	System.err.println("st "+st);
				if(toexcludefinal!=null && matches(toexcludefinal, st)){
					offset++;
				}
				
			}
			start = target+1+offset;
			
			// TODO Auto-generated method stub

		}
		public String closeInner() throws Exception {
			//if(!active()) return null;
			this.target = reports[0].numsnps-1;
			this.runInner();
			if(st==null) bread.close();
			//else{
				//System.err.println("next line");
				//System.err.println(st);
			//}
			target=0;
			offset=0;
			start =0;
			String sample1 = sample;
			
			sample = null;
			return sample1;
		}

	}



	//static Pattern toexcludesnp = null;
	static Pattern toexcludefinal = null;
	//Name    GenomeBuild     Chr     MapInfo
	public static void main( String[] args1){
		if(args1.length==2) {
			System.err.println("merging");
			mainMerge(args1);
			return;
		}
		final String[] args = args1;
		System.err.println("here!!!");
		
		//snpfile loc numperrun:maxthreads finalreportexclude:snpexclude
//		HumanOmni2.5-4v1_H_SNPlist.txt  1:0:1000000  1:2  snp kgp18461489,|kgp22734341,:null
		try{
		
			File f = new File( System.getProperty("user.dir"));
		
			
			final File snpf = new File(f,args[0]);
			final File target = new File(f,  "target");
			System.err.println("target file : "+target.getAbsolutePath());
			//final  String[] locs = args[1].split(":");
			String[] params = args[2].split(":");
			final int numperrun = Integer.parseInt(params[0]);
			int maxthreads = Integer.parseInt(params[1]);
			String prefix = null;
			if(!args[3].equals("null")){
				prefix = args[3];
			}
			final  File out = new File(f, prefix+"output");
			String[] toexcl = args[4].split(":");
			System.err.println("printing args "+args.length+" "+args1.length);
			System.err.println(Arrays.asList(args).toString());
			System.err.println(Arrays.asList(toexcl).toString());
			System.err.println("done");
			final Pattern excludefinalreport = toexcl[0].equals("null")? null : Pattern.compile(toexcl[0]) ;
			final Pattern toexcludesnps = toexcl[1].equals("null")? null : Pattern.compile(toexcl[1]) ;
			boolean removeInternal = Boolean.parseBoolean(args[5]);
			writeToFile = Boolean.parseBoolean(args[6]);
			out.mkdir();
			final List<SNP> snps = new ArrayList<SNP>();
			String[] locs = read(new File(f, args[1]));
			final String[] chroms = new String[locs.length];
			final int[][] pos = new int[locs.length][];
			
			getSNPS(snpf, target, toexcludesnps, locs, snps, chroms,pos);
			System.err.println("SNPS: ");
			System.err.println(snps.size());
			File[][] reports = getReports(f, prefix, numperrun);
			if(reports.length==0){
				Logger.global.info(f.getAbsolutePath());
				Logger.global.info(prefix);
				throw new RuntimeException("found no reports");
			}
			Thread[] t = new Thread[reports.length];
			ThreadGroup tg = new ThreadGroup("threads");
			System.err.println("running "+t.length+" threads");
			for(int k=0; k<t.length; k++){
				System.err.println("thread  "+k+" "+reports[k].length);
			}
			final List<File>[] outf = new List[chroms.length];//
			String[] outnme = new String[chroms.length];
			for(int k=0; k<outf.length; k++){
				outf[k] = new ArrayList<File>();
				outnme[k] = chroms[k]+"_"+pos[k][0]+"_"+pos[k][1];
			}
			for(int k=0; k<reports.length; k++){
				final int k1 = k;
				//final File out1 = new File(out, outnme+"."+numperrun+"."+k);
				
				//outf.add(out1);//new File(out, out1.getName()+".zip");
				//if(!out1.exists() || out1.length()<20){
					final File[] reportsk = reports[k];
					Runnable run = new Runnable(){
						@Override
						public void run() {
							try{
								System.err.println("Running batch "+k1);
								System.err.println(reportsk[0].getName());
								System.err.println(reportsk.length);//+":"+Arrays.asList(reportsk));
								
								List<FinalReport> reports =  open(reportsk, vals);
								//
								//Arrays.fill(active, true);
								//boolean active = true;
								int count=0;
								ConvertFinalReport cfr = null;
								 
								try{
								while(reports.size()>0 && count < countmax){
									if(count==0 || !writeToFile){
										if(cfr!=null) cfr.close();
										String suffix = "."+numperrun+"."+k1+"."+count;
										cfr = new ConvertFinalReport( out, suffix, snps, chroms,pos, reports.toArray(new FinalReport[0]));
										for(int i=0; i< outf.length; i++){
											outf[i].add(cfr.outf[i]);
										}
									}
									cfr.extract(excludefinalreport);
									boolean[] active = new boolean[reports.size()];
									
									cfr.reset(active);
									for(int k=active.length-1; k>=0; k--){
										if(!active[k]){
											reports.remove(k);
										}
									}
									count++;
								}
								}catch(Exception exc){
									exc.printStackTrace();
								}
								cfr.close();
							}catch(Exception exc){
								exc.printStackTrace();
							}
						}
					};
					t[k] = new Thread(tg,run);
				//}
				
			}
			int kk=0;
			while(kk < t.length){
				int tc = tg.activeCount();
				for(int kj=0; kj <maxthreads-tc && kk < t.length; kj++)
				{
					if(t[kk]!=null) t[kk].start();
					kk++;
				}
				//kk = kk+maxthreads - tc;
					Thread.sleep(1000);
			}
			while(tg.activeCount()>0) Thread.sleep(1000);
			System.err.println("merging");
			for(int k=0; k<outf.length; k++){
			ConvertFinalReport.merge(new File(out, outnme[k]), outf[k], removeInternal,10);
			}
			 if(args.length>2) toexcludefinal = Pattern.compile(args[2]);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static void mainMerge(String[] args){
		try{
			File f = new File( System.getProperty("user.dir"));
			final String prefix = args[0];
			String[] locs = read(new File(f, args[1]));
			List<File>[] outf = new List[locs.length];
			String[] outnme = new String[locs.length];
			final String[] outnme1 = new String[locs.length];
			for(int j=0; j<locs.length; j++){
				String[] locsj = locs[j].split(":");
				outnme[j] = locsj[0]+"_"+locsj[1]+"_"+locsj[2];
				outnme1[j] = outnme[j]+".zip";
			}
			File[] dirs = f.listFiles(new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().startsWith(prefix) && Arrays.asList(pathname.list()).contains(outnme1[0]);
				}
				
			});
			for(int j=0; j<locs.length; j++){
				outf[j] = new ArrayList<File>();
				for(int k=0; k<dirs.length; k++){
					outf[j].add(new File(dirs[k],outnme1[j]));
				}
			}
			boolean removeInternal = true;
			final  File out = new File(f, "output");
			out.mkdir();
			for(int k=0; k<outf.length; k++){
				ConvertFinalReport.merge(new File(out, outnme[k]), outf[k], removeInternal,10);
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	private static String[] read(File file) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String st = "";
		List<String> locs = new ArrayList<String>();
		while((st = br.readLine())!=null){
			locs.add(st);
		}
		return  locs.toArray(new String[0]);
	}

	private void extract(Pattern p)throws Exception{
		toexcludefinal = p;
		for(Iterator<SNP> it = this.snps.iterator(); it.hasNext();){
			try{
			SNP snp = it.next();
			CompressDir compress = allcompress[snp.which];
		
				//System.err.println(snp);
			for(int k=0; k< this.run.length; k++){
				if(run[k].active()){
					this.run[k].setTarget(snp.setind);
					System.err.println(snp.setind);
					this.run[k].runInner();
					
				}
			}
			OutputStreamWriter os= compress.getWriter(snp.id, !writeToFile);
			for(int k=0; k< this.run.length; k++){
				if(run[k].st!=null){
					String[] str = this.run[k].st.split(splitstr);
				//Logger.global.info(Arrays.asList(str).toString());
					if(!str[reports[k].id_index].equals(snp.id)){
							throw new RuntimeException("mismatch "+snp.toString()+" vs \n  "+Arrays.asList(str));
					}
					snp.write(os, str, this.reports[k].inds);
				}
			}
			compress.closeWriter(os);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		if(reports[0].numsnps==0) throw new RuntimeException("!!");
		/*for(int i=0; i< this.run.length; i++){
			run[i].setTarget(numsnps-1);
			run[i].run();
			
		}*/
		
		
	}

	 
	static File[][] getReports(File dir, final String prefix, int numpergroup){
		List<File> res;
		FileFilter ff = new FileFilter(){
			@Override
			public boolean accept(File arg0) {
				return arg0.getName().indexOf("FinalReport")>=0;
			}
		};
		
		if(prefix!=null && !prefix.equals("null")){
			FileFilter ff1 = new FileFilter(){
				@Override
				public boolean accept(File arg0) {
					return arg0.isDirectory() && arg0.getName().startsWith(prefix);
				}
			};
			res = new ArrayList<File>();
			File[] dirs = dir.listFiles(ff1);
			for(int k=0; k<dirs.length; k++){
				res.addAll(Arrays.asList(dirs[k].listFiles(ff)));
			}
		}
		else res = Arrays.asList(dir.listFiles(ff));
		
		int num = (int) Math.ceil((double)res.size()/(double) numpergroup);
		File[][] ress = new File[num][];
		for(int k=0; k<ress.length; k++){
			
			ress[k] = res.subList(k*numpergroup, Math.min(res.size(), (k+1)*numpergroup)).toArray(new File[0]);
		}
		return ress;
	}
	
	static class SNP {
		String id;
		final Integer pos, index, setind, which;
		
		SNP(int pos, String name,  int k, int setind, int which){
			index = k;
			id = name;
			this.pos = pos;
			this.setind = setind;
			this.which = which;
		}
		public SNP(SNP snp) {
			this.index = snp.index;
			this.id = snp.id;
			this.pos = snp.pos;
			this.setind = snp.setind;
			this.which  = snp.which;
			// TODO Auto-generated constructor stub
		}
		public void write(OutputStreamWriter os, String[] str, int[] is) throws Exception{
			for(int k=0; k<is.length; k++){
				os.write(str[is[k]]);
				os.write(k<is.length-1?"\t":"\n");
			}
			
		}
		
		public String toString(String chrom){
			return chrom+"\t"+toString();
		}
		public String toString(){
			return pos+"\t"+(pos+20)+"\t"+id+"\t"+index + (setind==null ? "": "\t"+setind);
		}
	}
	
	static class SNPComp implements Comparator<SNP>{
		final boolean usePos;
		public SNPComp(boolean usePos){
			this.usePos = usePos;
		}
		//boolean index = true;
		@Override
		public int compare(SNP arg0, SNP arg1) {
		//	if(index) return arg0.index.compareTo(arg1.index);
			if(usePos) return arg0.pos.compareTo(arg1.pos);
			else return arg0.setind.compareTo(arg1.setind);
		}
		
	};
	
	
	
	List<SNP> snps = new ArrayList<SNP>();
	
	String[] chroms;

	FinalReport[] reports;
	ReadRun[] run;
	String[] samples;
	
	
	final static String splitstr = "\t";
	static void  getSNPS(File snpsfile, File target, Pattern toexcludesnp,String[] locs, List<SNP> snps, String[] chroms, int[][]from) throws Exception{
		List<String> set = null;
		
		if(target.exists()){
			set = new ArrayList<String>();
			BufferedReader br = new BufferedReader(new FileReader(target));
			String st = "";
			while((st = br.readLine())!=null){
				set.add(st.split(splitstr)[0]);
			}
			System.err.println("total "+set.size()+" first: "+set.get(0));
			br.close();
		}else{
			System.err.println("no target file");
		}
		String[] chrom = new String[locs.length];
	//	int[][] from = new int[locs.length][];
		
		for(int k=0; k<locs.length; k++){
			String[] lock = locs[k].split(":");
			chrom[k] = lock[0];
			chroms[k] = lock[0];
			from[k] = (new int[] {Integer.parseInt(lock[1]),
		Integer.parseInt(lock[2])});
		}
	
	
		BufferedReader br = new BufferedReader(new FileReader(snpsfile));
		List<String> head = Arrays.asList(br.readLine().split("\\s+"));
		String st = "";
		int chromind =head.indexOf("Chr");
		int posind = head.indexOf("MapInfo");
		int nameind = head.indexOf("Name");
		for(int k=0; (st = br.readLine())!=null ; k++){
			if(toexcludesnp==null || !matches(toexcludesnp, st)){
				String[] str  = st.split("\\s+");
				inner: for(int kk=0; kk<chrom.length; kk++){
				if(str[chromind].equals( chrom[kk])){
					int pos = Integer.parseInt(str[posind]);
					if(pos>=from[kk][0] && pos <from[kk][1]){
						int setind = set==null ? k : set.indexOf(str[nameind]);
						if(setind>=0){
							SNP snp = new SNP(pos, str[nameind], k, setind, kk);
							snps.add(snp);
							continue inner;
						}
					//	System.err.println(snp);
					}
				}
				}
			}else{
				System.err.println("matched "+st);
			}
		}
		br.close();
		if(set != null){
			Collections.sort(snps, new ConvertFinalReport.SNPComp(false));
		}
		System.err.println("snps ");
		for(int k=0; k<snps.size(); k++){
		System.err.println(snps.get(k));
		}
		
	}
	
	//int sampleid = -1;
	//int numsnps =0;
	static String val = "B Allele Freq,Log R Ratio";
	static String[] vals = val.split(",");
	static class FinalReport{
		final BufferedReader br;
		int cnt=0;
		final File in;
		List<String> header;
		int numsnps;
		int sampleid;
		final int[]inds;
		int id_index;
		//Name    GenomeBuild     Chr     MapInfo

		
		
		String snpname = "SNP Name";
		FinalReport(File reports, String[] vals) throws Exception{
			inds = new int[vals.length];
			this.in = reports;
				if(reports.getName().endsWith(".gz")){
				   br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(reports))));
				}else{
					br =  new BufferedReader(new InputStreamReader((new FileInputStream(reports))));
				}
				
				String str = "";
				while(!(str = br.readLine()).startsWith("[Data]")){	
					cnt++;
					if(str.toLowerCase().startsWith("num snps")){
						int numsnps1 = Integer.parseInt(str.split(splitstr)[1]);
						if(numsnps>0 && numsnps1 !=numsnps) throw new RuntimeException("diff num snps");
						numsnps = numsnps1;
					}
					//System.err.println(str);
				}
				List<String> header = Arrays.asList(br.readLine().split(splitstr));
				cnt++;
				for(int kk=0; kk<header.size(); kk++){
					if(header.get(kk).toLowerCase().indexOf("sample")>=0) sampleid = kk;
				}
				Logger.global.info("header:");
				Logger.global.info(header.toString());
				for(int j=0; j<vals.length; j++){
					inds[j] = header.indexOf(vals[j]);
					if(inds[j]<0){
						Logger.global.info("problem finding "+vals[j]);
					}
				}
				this.id_index = header.indexOf(this.snpname);
				if(id_index<0){
					throw new RuntimeException("did not find index "+this.snpname+" in "+header);
				}
			}
	}
	
	static List<FinalReport> open(File[] reports, String[] vals) throws Exception{
		List<FinalReport> res = new ArrayList<FinalReport>();
		for(int k=0; k<reports.length; k++){
			
		res.add( new FinalReport(reports[k], vals));	
		}
		return res;
	}
	
	ConvertFinalReport( File outdir, String suffix, List<SNP>snps1, String[] chroms, int[][] pos, FinalReport[]reports) throws Exception{
		Logger.global.info("input ");
	//	for(int k=0; k<reports.length; k++){
	//		Logger.global.info(reports[k].getAbsolutePath());
	//	}
		this.reports = reports;
		Logger.global.info("output "+outdir.getAbsolutePath());
		this.allcompress =new CompressDir[chroms.length];
		this.outf = new File[chroms.length];
		this.samplespw = new OutputStreamWriter[chroms.length];
		this.chroms = chroms;
		this.snps = new ArrayList(snps1.size());
		for(int k=0; k<snps1.size(); k++){
			this.snps.add(new SNP(snps1.get(k)));
		}
		for(int k=0; k<chroms.length; k++){
			// final String outnme = args[1].replace(':', '_');
			 outf[k] = new File(outdir, chroms[k]+"_"+pos[k][0]+"_"+pos[k][1]+suffix);
			
			allcompress[k] = new CompressDir(outf[k]);
			CompressDir compress = allcompress[k];
		OutputStreamWriter pw = compress.getWriter("Name", true);
		
		for(int kk=0; kk<this.vals.length; kk++){
			pw.write(vals[kk]);
			pw.write(kk<vals.length-1 ? "\t":"\n");
		}
		pw.write("chr\tstart\tend\tsnpid\tindex\n");
		pw.write("Sample\n");
		compress.closeWriter(pw);
		 samplespw[k] = compress.getWriter("Samples", false);
		}{
	     //this.id_index = new int[reports.length];
	    // this.br = new BufferedReader[reports.length];
	     this.run = new ReadRun[reports.length];
	     //this.inds = new int[reports.length][vals.length];
		
		for(int k=0; k<reports.length; k++){
			
			BufferedReader br1 = reports[k].br;
		
			
			run[k] = new ReadRun(reports[k].in, br1, reports[k].cnt);
		}
		
		}
	}
	final OutputStreamWriter[] samplespw;
	
	
	
	void reset(boolean[] active) throws Exception{
		//boolean active = false;
		
		for(int i=0; i< this.run.length; i++){
			
			for(int kk=0; kk<allcompress.length; kk++){
			//	CompressDir compress = allcompress[kk];
				if(run[i].sample!=null){
					samplespw[kk].write(this.run[i].sample);
					samplespw[kk].write("\n");
					Logger.global.info("Writing sample "+run[i].sample+" "+run[i].input.getName());
				}else{
					Logger.global.warning("Sample was null for "+run[i].input.getName());
					}
			
			}
			
			
		}
		for(int i=0; i<run.length; i++){
			 run[i].closeInner();
			 active[i] =  run[i].active();
		}
		
		
			
		// return active;
	
	}
	public void close() throws Exception{
		Collections.sort(snps, new SNPComp(true));
		OutputStreamWriter[] pws = new OutputStreamWriter[allcompress.length];
		for(int kk=0; kk<allcompress.length; kk++){
			allcompress[kk].closeWriter(samplespw[kk]);
			 pws[kk] = allcompress[kk].getWriter("SNPS", true);
		}
		for(Iterator<SNP> it = snps.iterator(); it.hasNext();){
			SNP snp = it.next();
				pws[snp.which].write(snp.toString(this.chroms[snp.which]));
				pws[snp.which].write("\n");
		}
		for(int kk=0; kk<allcompress.length; kk++){
			CompressDir compress = allcompress[kk];
			compress.closeWriter(pws[kk]);
			compress.run();
			compress.close();
	}
	}
	
final 	CompressDir[] allcompress;
final File[] outf;



public static void merge(File nme, List<File> zips, boolean removeInternal, int numperbin) throws Exception{
	int len = zips.size();
	System.err.println(zips);
	for(int k=0; zips.size()>numperbin; k++){
		System.err.println(k);
		zips = mergecombine(zips, removeInternal, numperbin);
		System.err.println(zips);
		if(zips.size()>=len) throw new RuntimeException("this should be decreasing "+zips.size()+" "+len);
	}
	System.err.println(zips);
	 mergeinternal(nme, zips, removeInternal);
}
//takes in zip files, and merges them to fewer based on numperbin
private static List<File> mergecombine(List<File> in, boolean removeInternal, int numperbin) throws Exception{
	int len = (int) Math.ceil((double) in.size()/(double)numperbin);
	File[] out = new File[len];
	for(int k=0; k<out.length; k++){
		List<File> zips_ = in.subList(k*numperbin, Math.min(in.size(),(k+1)*numperbin));
		out[k] =  new File(zips_.get(0).getAbsolutePath()+".i");
		mergeinternal(out[k], zips_, removeInternal);
	}
	return Arrays.asList(out);
}

static void mergeinternal(File nme, List<File> zips, boolean removeInternal) throws Exception{
	CompressDir dir = new CompressDir(nme);
	ZipFile[] zf = new ZipFile[zips.size()];
	for(int k=0; k<zf.length; k++){
		File zipf = zips.get(k).getName().endsWith(".zip") ? new File(zips.get(k).getAbsolutePath()): new File(zips.get(k).getAbsolutePath()+".zip");
		Logger.global.info("opening "+zipf.getAbsolutePath());
		zf[k] = new ZipFile(zipf);
		if(removeInternal)zipf.deleteOnExit();
	}
	List<String> nomerge = Arrays.asList("Name:SNPS".split(":"));
	Enumeration en= zf[0].entries();
	while(en.hasMoreElements()){
		ZipEntry ze = (ZipEntry)en.nextElement();
		OutputStreamWriter osw = dir.getWriter(ze.getName(), true);
		for(int k=0; k<zf.length; k++){
			if(k==0 || ! nomerge.contains(ze.getName())){
				BufferedReader br = new BufferedReader(new InputStreamReader(zf[k].getInputStream(zf[k].getEntry(ze.getName()))));
				String st = "";
				while((st=br.readLine())!=null){
					osw.write(st);
					osw.write("\n");
				}
			}
		}
		dir.closeWriter(osw);
	}
	dir.run();
	dir.close();
}


}
