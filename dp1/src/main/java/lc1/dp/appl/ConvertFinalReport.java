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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class ConvertFinalReport {
	
	static boolean matches(Pattern p, String st){
		Matcher m = toexcludefinal.matcher(st);
		m.useAnchoringBounds(true);
		return m.find();
	}
	
	public class ReadRun implements Runnable {
		BufferedReader bread;
		String st = "";
		int offset=0;
		int target = 0;
		int start =0;
		public ReadRun(BufferedReader br1){
			this.bread = br1;
		}
		void setTarget(int pos){
			this.target = pos;
		}
		@Override
		public void run() {
			try{
			for(int k=start; k<=target+offset; k++){
				st = bread.readLine();
				if(toexcludefinal!=null && matches(toexcludefinal, st)){
					offset++;
				}
			}
			start = target+1+offset;
			}catch(Exception exc){
				exc.printStackTrace();
			}
			// TODO Auto-generated method stub

		}
		public void close() throws Exception {
			this.bread.close();
			
		}

	}



	//static Pattern toexcludesnp = null;
	static Pattern toexcludefinal = null;
	//Name    GenomeBuild     Chr     MapInfo
	public static void main( String[] args1){
		final String[] args = args1;
		//snpfile loc numperrun:maxthreads finalreportexclude:snpexclude
//		HumanOmni2.5-4v1_H_SNPlist.txt  1:0:1000000  1:2  snp kgp18461489,|kgp22734341,:null
		try{
		
			File f = new File( System.getProperty("user.dir"));
			final  File out = new File(f, "output");
			final File snpf = new File(f,args[0]);
			final  String[] locs = args[1].split(":");
			String[] params = args[2].split(":");
			int numperrun = Integer.parseInt(params[0]);
			int maxthreads = Integer.parseInt(params[1]);
			String prefix = null;
			if(!args[3].equals("null")){
				prefix = args[3];
			}
			String[] toexcl = args[4].split(":");
		
			final Pattern excludefinalreport = toexcl[0].equals("null")? null : Pattern.compile(toexcl[0]) ;
			final Pattern toexcludesnps = toexcl[1].equals("null")? null : Pattern.compile(toexcl[1]) ;
			boolean removeInternal = Boolean.parseBoolean(args[5]);
			out.mkdir();
			final   List<SNP> snps =  getSNPS(snpf, toexcludesnps, locs);
			System.err.println("SNPS: ");
			System.err.println(snps.size());
			File[][] reports = getReports(f, prefix, numperrun);
			Thread[] t = new Thread[reports.length];
			ThreadGroup tg = new ThreadGroup("threads");
			String outnme = args[1].replace(':', '_');
			File[] outf = new File[reports.length];
			for(int k=0; k<reports.length; k++){
				final int k1 = k;
				final File out1 = new File(out, outnme+"."+numperrun+"."+k);
				outf[k] = out1;//new File(out, out1.getName()+".zip");
				if(!outf[k].exists() || outf[k].length()<20){
					final File[] reportsk = reports[k];
					Runnable run = new Runnable(){
						@Override
						public void run() {
							try{
								System.err.println("Running batch "+k1);
								System.err.println(reportsk[0].getName());
								System.err.println(reportsk.length);//+":"+Arrays.asList(reportsk));
							ConvertFinalReport cfr = new ConvertFinalReport(out1, snps, locs[0], reportsk);
							cfr.extract(excludefinalreport);
							cfr.close();
							}catch(Exception exc){
								exc.printStackTrace();
							}
						}
					};
					t[k] = new Thread(tg,run);
				}
				
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
			ConvertFinalReport.merge(new File(out, outnme), Arrays.asList(outf), removeInternal,10);
			 if(args.length>2) toexcludefinal = Pattern.compile(args[2]);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	private void extract(Pattern p)throws Exception{
		toexcludefinal = p;
		for(Iterator<SNP> it = this.snps.iterator(); it.hasNext();){
			SNP snp = it.next();
			OutputStreamWriter os= compress.getWriter(snp.id, true);
				//System.err.println(snp);
			for(int k=0; k< this.run.length; k++){
				this.run[k].setTarget(snp.index);
				this.run[k].run();
				String[] str = this.run[k].st.split(",");
				if(!str[id_index[k]].equals(snp.id)){
					throw new RuntimeException("mismatch "+snp.id+" "+str[id_index[k]]);
				}
				snp.write(os, str, this.inds[k]);
			}
			compress.closeWriter(os);
		}
		for(int i=0; i< this.run.length; i++){
			run[i].close();
		}
		
	}

	 
	static File[][] getReports(File dir, final String prefix, int numpergroup){
		List<File> res;
		FileFilter ff = new FileFilter(){
			@Override
			public boolean accept(File arg0) {
				return arg0.getName().indexOf("_FinalReport_")>=0;
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
		Integer pos, index;
		
		SNP(int pos, String name,  int k){
			index = k;
			id = name;
			this.pos = pos;
		}
		public SNP(SNP snp) {
			this.index = snp.index;
			this.id = snp.id;
			this.pos = snp.pos;
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
			return pos+"\t"+(pos+20)+"\t"+id+"\t"+index;
		}
	}
	
	class SNPComp implements Comparator<SNP>{
		//boolean index = true;
		@Override
		public int compare(SNP arg0, SNP arg1) {
		//	if(index) return arg0.index.compareTo(arg1.index);
			 return arg0.pos.compareTo(arg1.pos);
		}
		
	};
	
	
	
	List<SNP> snps = new ArrayList<SNP>();
	
	String chrom;

	BufferedReader[] br;
	ReadRun[] run;
	int[][]inds;
	int[] id_index;
	//Name    GenomeBuild     Chr     MapInfo

	String val = "X Raw,Y Raw,B Allele Freq,Log R Ratio";
	String[] vals = val.split(",");
	String snpname = "SNP Name";
	
	static List<SNP> getSNPS(File snpsfile, Pattern toexcludesnp,String[] locs) throws Exception{
		
		List<SNP> snps= new ArrayList<SNP>();
		String chrom = locs[0];
		int from = Integer.parseInt(locs[1]);
		int to = Integer.parseInt(locs[2]);
		BufferedReader br = new BufferedReader(new FileReader(snpsfile));
		List<String> head = Arrays.asList(br.readLine().split("\\s+"));
		String st = "";
		int chromind =head.indexOf("Chr");
		int posind = head.indexOf("MapInfo");
		int nameind = head.indexOf("Name");
		for(int k=0; (st = br.readLine())!=null ; k++){
			if(toexcludesnp==null || !matches(toexcludesnp, st)){
				String[] str  = st.split("\\s+");
				if(str[chromind].equals( chrom)){
					int pos = Integer.parseInt(str[posind]);
					if(pos>=from && pos <to){
					
						SNP snp = new SNP(pos, str[nameind], k);
						snps.add(snp);
					//	System.err.println(snp);
					}
				}
			}else{
				System.err.println("matched "+st);
			}
		}
		br.close();
		return snps;
	}
	
	ConvertFinalReport(File outdir, List<SNP>snps1, String chrom, File[]reports) throws Exception{
		this.compress = new CompressDir(outdir);
		this.chrom = chrom;
		this.snps = new ArrayList(snps1.size());
		for(int k=0; k<snps1.size(); k++){
			this.snps.add(new SNP(snps1.get(k)));
		}
		{
		OutputStreamWriter pw = this.compress.getWriter("Name", true);
		
		for(int k=0; k<this.vals.length; k++){
			pw.write(vals[k]);
			pw.write(k<vals.length-1 ? "\t":"\n");
		}
		pw.write("chr\tstart\tend\tsnpid\tindex\n");
		pw.write("Sample\n");
		compress.closeWriter(pw);
		}{
	     this.id_index = new int[reports.length];
	     this.br = new BufferedReader[reports.length];
	     this.run = new ReadRun[reports.length];
	     this.inds = new int[reports.length][vals.length];
		OutputStreamWriter pw = this.compress.getWriter("Samples", true);
		for(int k=0; k<reports.length; k++){
			pw.write(reports[k].getName());
			pw.write("\n");
			BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(reports[k]))));
			
			String str = "";
			while(!(str = br1.readLine()).startsWith("[Data]")){			}
			List<String> header = Arrays.asList(br1.readLine().split(","));
			for(int j=0; j<vals.length; j++){
				inds[k][j] = header.indexOf(vals[j]);
			}
			this.id_index[k] = header.indexOf(this.snpname);
			run[k] = new ReadRun(br1);
		}
		compress.closeWriter(pw);
		}
	}
	
	void close() throws Exception{
		OutputStreamWriter pw = this.compress.getWriter("SNPS", true);
		Collections.sort(snps, new SNPComp());
		for(Iterator<SNP> it = snps.iterator(); it.hasNext();){
			pw.write(it.next().toString(this.chrom));
			pw.write("\n");
		}
		compress.closeWriter(pw);
		compress.run();
		 compress.close();
	}
	
final 	CompressDir compress;




public static void merge(File nme, List<File> zips, boolean removeInternal, int numperbin) throws Exception{
	int len = zips.size();
	for(int k=0; zips.size()>numperbin; k++){
		zips = mergecombine(zips, removeInternal, numperbin);
		if(zips.size()>=len) throw new RuntimeException("this should be decreasing "+zips.size()+" "+len);
	}
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
		File zipf = new File(zips.get(k).getAbsolutePath()+".zip");
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
