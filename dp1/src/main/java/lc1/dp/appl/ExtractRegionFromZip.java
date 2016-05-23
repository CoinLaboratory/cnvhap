package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class ExtractRegionFromZip {

	public static void main( String[] args){
	//	if(true) System.exit(0);
	//	infile outfile locs
		try{
		final File dir = new File(args[0]);
		final File outdir = new File(args[1]);
		outdir.mkdir();
		final File regions = new File(args[2]);
		BufferedReader br1 = new BufferedReader(new FileReader(regions));
		String str = "";
		Map<String, List<String>> ls = new HashMap<String, List<String>>();
		boolean groupByChrom = false;
		while((str= br1.readLine())!=null){
			String chr = str.split(":")[0];
			List<String>l = ls.get(chr);
			if(l==null) ls.put(chr, l = new ArrayList<String>());
			l.add(str);
		}
		
			for(Iterator<String> it = ls.keySet().iterator(); it.hasNext();){
				String chr = it.next();
				File in = new File(dir, chr+".zip");
				List<String> locs = ls.get(chr);
				int[] from = new int[locs.size()];
				int[] to = new int[locs.size()];
				for(int i=0; i<locs.size(); i++){
					String[] loc = locs.get(i).split("\\s+")[0].split(":");
					from[i] = Integer.parseInt(loc[1]);
					to[i] = Integer.parseInt(loc[2]);
					if(!groupByChrom) extractregion(in, outdir, new int[] {from[i]}, new int[] {to[i]}, chr);
				}
				if(groupByChrom) extractregion(in, outdir, from, to, chr);
			}
		
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
public static void extractregion(File in, File outdir,int[] from, int[]to, String chrom ){
	try{
		
		String outname = in.getName().split("\\.")[0];
		if(!outdir.exists()){
			outdir.mkdirs();
		}
		
		outname = outname +(in.getName().equals("all") ? "_"+chrom: "")+"_"+from[0]+"_"+to[to.length-1];
		File outf = new File(outdir,outname);
		CompressDir compress = new CompressDir(outf);
		ZipFile zf = new ZipFile(in);
	//	String chrom = loc[0];
		String chrom1 = "chr"+chrom;
		List<String> header_snp = Arrays.asList(lc1.util.Compressor.readZipFrom(zf, "Name").get(1).toLowerCase().split("\\s+"));
		int snpid_ = header_snp.indexOf("id");
	    if(snpid_<0) snpid_ = header_snp.indexOf("snpid");
	    int chrid = header_snp.indexOf("chr");
	    int startid = header_snp.indexOf("start");
		BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
		String st = "";
		OutputStreamWriter os = compress.getWriter("SNPS", false);
		compress.copy(zf, "Name");
		compress.copy(zf, "Samples");
		int cnt=0;
		while((st = snps.readLine())!=null){
			String[] str = st.split("\\s+");
			if(str[chrid].equals(chrom) ||str[chrid].equals(chrom1) ){
				int pos = Integer.parseInt(str[startid]);
				inner: for(int j=0; j<from.length; j++){
				if(pos<=to[j] && pos>=from[j]){
					cnt++;
					compress.copy(zf, str[snpid_]);
					os.write(st);
					os.write("\n");
					break inner;
				}
				}
			}
		}
		compress.closeWriter(os);
		compress.run();
		compress.close();
		if(cnt==0) compress.outFile.deleteOnExit();		
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
