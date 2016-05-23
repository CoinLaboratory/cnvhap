package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class ExtractSamplesFromZip {
	public static void main( String[] args){
		//	if(true) System.exit(0);
		//	infile outfile locs
			try{
			final File dir = new File(args[0]);
			final File outdir = new File(args[1]);
			outdir.mkdir();
			final File samplesF = new File(args[2]);
			BufferedReader br1 = new BufferedReader(new FileReader(samplesF));
			String str = "";
			List<String> samples = new ArrayList<String>();
			while((str= br1.readLine())!=null){
				
				samples.add(str);
			}
				File[] ins = dir.listFiles(new FilenameFilter(){

					@Override
					public boolean accept(File dir, String name) {
						return name.endsWith(".zip");
					}
					
				});
				for(int i=0; i<ins.length; i++){
				 extractsamples(ins[i], outdir, samples);
				}
			
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		
		
	public static void extractsamples(File in, File outdir,List<String> samplesToSelect ){
		try{
			String outname = in.getName().split("\\.")[0];
			if(!outdir.exists()){
				outdir.mkdirs();
			}
			outname = outname +"_"+samplesToSelect.size();
			File outf = new File(outdir,outname);
			CompressDir compress = new CompressDir(outf);
			ZipFile zf = new ZipFile(in);
		//	String chrom = loc[0];
			List<String> samples = lc1.util.Compressor.readZipFrom(zf, "Samples");
			boolean[] include = new boolean[samples.size()];
			Arrays.fill(include, false);
			for(int i=0; i<include.length; i++){
				if(samplesToSelect.contains(samples.get(i).split("\\s+")[0])) include[i] = true;
			}
			List<String> header_snp = Arrays.asList(lc1.util.Compressor.readZipFrom(zf, "Name").get(1).toLowerCase().split("\\s+"));
			int snpid_ = header_snp.indexOf("id");
		    if(snpid_<0) snpid_ = header_snp.indexOf("snpid");
		    int chrid = header_snp.indexOf("chr");
		    int startid = header_snp.indexOf("start");
			BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
			String st = "";
			OutputStreamWriter os = compress.getWriter("SNPS", false);
			compress.copy(zf, "Name");
			compress.copy(zf, "Samples", include);
			while((st = snps.readLine())!=null){
				String[] str = st.split("\\s+");
				compress.copy(zf, str[snpid_], include);
				os.write(st);
				os.write("\n");
			}
			compress.closeWriter(os);
			compress.run();
			compress.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	}


