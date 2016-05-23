package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.zip.ZipFile;

public class ConvertZipToVCF {
	// work in progress
	public static void main( String[] args){
		//	if(true) System.exit(0);
		//	infile outfile locs
			try{
			final File dir = new File(args[0]);
			final File outdir = new File(args[1]);
			outdir.mkdir();
			final List<String> columns = Arrays.asList(args[2].split(":"));
			for(int i=0; i<columns.size();i++) columns.set(i, slug(columns.get(i)));
				File[] ins = dir.listFiles(new FilenameFilter(){

					@Override
					public boolean accept(File dir, String name) {
						return name.endsWith(".zip");
					}
					
				});
				File outf = new File(outdir,"all.vcf");
				Arrays.sort(ins, new Comparator<File>(){

					@Override
					public int compare(File o1, File o2) {
						String[] n1 = o1.getName().split("_");
						String[] n2 = o2.getName().split("_");
						for(int i=0; i<n1.length; i++){
							if(i==0){
								Integer i1 = Integer.parseInt(n1[0]);
								Integer i2 = Integer.parseInt(n2[0]);
								if(i1.compareTo(i2)!=0) return i1.compareTo(i2);
							}else{
							if(n1[i].compareTo(n2[i])!=0) return n1[i].compareTo(n2[i]);
							}
						}
						return 0;
					}

					
					
				});
				for(int i=0; i<ins.length; i++){
				 convertToVCF(ins[i], outf, columns);
				}
				out.close();
			
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
	public static String slug(String input){
		return input.toLowerCase().replaceAll("[^a-z0-9-]", "");
	}
	public static void write(List<String> vals, PrintWriter out){
		for(int i=0; i<vals.size(); i++) {
			out.write(vals.get(i));
			if(i<vals.size()-1) out.write("\t");
			
		}
	}
	static List<String> samples = null;
	static PrintWriter out = null;
	public static void convertToVCF(File in, File outf,List<String> columns ){
		try{
			String outname1 = in.getName().split("\\.")[0];
			ZipFile zf = new ZipFile(in);
		//	String chrom = loc[0];
			List<String> header_geno = Arrays.asList(lc1.util.Compressor.readZipFrom(zf, "Name").get(0).toLowerCase().split("\t"));
		
			boolean[] include = new boolean[header_geno.size()];
			for(int i=0; i<include.length; i++){
				if(columns.contains(slug(header_geno.get(i)))){
					include[i] = true;
				}
			}
			List<String> header_snp = Arrays.asList(lc1.util.Compressor.readZipFrom(zf, "Name").get(1).toLowerCase().split("\\s+"));
			int snpid_ = header_snp.indexOf("id");
		    if(snpid_<0) snpid_ = header_snp.indexOf("snpid");
		    int chrid = header_snp.indexOf("chr");
		    int startid = header_snp.indexOf("start");
			BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
			String st = "";
			
			List<String> samples1 = lc1.util.Compressor.readZipFrom(zf, "Samples", new boolean[] {true}, ":");
			if(samples==null){
				 out = new PrintWriter(new FileWriter(outf));
				samples = samples1;
				out.write("#chr\tregion\tpos\tsnpid\t");
				write(samples, out);
				out.write("\n");
			}else{
				if(!samples1.equals(samples)) throw new Exception ("sample files diff");
			}
			while((st = snps.readLine())!=null){
				String[] str = st.split("\\s+");
				List<String> vals = lc1.util.Compressor.readZipFrom(zf, str[snpid_], include, ":");
				out.write(str[chrid].replace("chr", "")+"\t"+outname1+"\t"+str[startid]+"\t"+str[snpid_]+"\t");
				write(vals, out);
				out.write("\n");
			}
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	}


