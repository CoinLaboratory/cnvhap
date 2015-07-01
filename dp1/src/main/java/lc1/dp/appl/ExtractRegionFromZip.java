package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class ExtractRegionFromZip {

	public static void main( String[] args){
	//	if(true) System.exit(0);
	//	infile outfile locs
		try{
		final File in = new File(args[0]);
		final File outdir = new File(args[1]);
		String locs = args[2];
			extractregion(in, outdir, locs);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
public static void extractregion(File in, File outdir, String locs ){
	try{
		String[] loc = locs.split(":");
		String outname = in.getName().split("\\.")[0];
		if(!outdir.exists()){
			outdir.mkdirs();
		}
		outname = outname +(in.getName().equals("all") ? "_"+loc[0]: "")+"_"+loc[1]+"_"+loc[2];
		File outf = new File(outdir,outname);
		CompressDir compress = new CompressDir(outf);
		ZipFile zf = new ZipFile(in);
		String chrom = loc[0];
		String chrom1 = "chr"+chrom;
		int from = Integer.parseInt(loc[1]);
		int to = Integer.parseInt(loc[2]);
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
		while((st = snps.readLine())!=null){
			String[] str = st.split("\\s+");
			if(str[chrid].equals(chrom) ||str[chrid].equals(chrom1) ){
				int pos = Integer.parseInt(str[startid]);
			
				if(pos<=to && pos>=from){
					compress.copy(zf, str[snpid_]);
					os.write(st);
					os.write("\n");
				}
			}
		}
		compress.closeWriter(os);
		compress.run();
		compress.close();
				
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
