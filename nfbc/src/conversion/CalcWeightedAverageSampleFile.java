package conversion;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import lc1.util.ApacheCompressor;

import org.apache.commons.compress.archivers.zip.ZipFile;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

public class CalcWeightedAverageSampleFile {
	
	public static void main(String[] args){
		try{
			CalcWeightedAverageSampleFile cwas = new CalcWeightedAverageSampleFile(new File("."));
			cwas.calc();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	File[] f;
	int[] len;
	List<String[]> samples_info =new ArrayList<String[]>();
	List<RealMatrix> samples = new ArrayList<RealMatrix>();
	List<String[]>headers = new ArrayList<String[]>();
	List<List<String>> original = new ArrayList<List<String>>();
	
	List<Double> snpsize = new ArrayList<Double>();
	int nocol,norow;
	File out;
	CalcWeightedAverageSampleFile(File dir) throws Exception{
		out = new File(dir, "Samples_avg");
		final Pattern p = Pattern.compile("\\d");
		
		f = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File arg0) {
				Matcher mat = p.matcher(arg0.getName().substring(0, 1));
				
				mat.useAnchoringBounds(false);
				boolean b = mat.find();
				int st = mat.regionStart();
				return st==0 && b && arg0.getName().endsWith("zip");
			//return  && arg0.getName().startsWith("\\d");
			}
			
		});
		double totsnps =0;
		System.err.println(Arrays.asList(f));
		for(int k=0; k<f.length; k++){
			ZipFile zf = new ZipFile(f[k]);
			headers.add(ApacheCompressor.getEntries(zf, "Name").get(2).split("\\s+"));
			
			List<String> origk = lc1.util.ApacheCompressor.getEntries(zf, "Samples");
			this.original.add(origk);
			samples.add(process(origk, headers.get(k)));
			String[] res = new String[samples.get(k).getRowDimension()];
			lc1.util.ApacheCompressor.readZip(zf, "Samples", res, 0);
			samples_info.add(res);
			snpsize.add((double) lc1.util.ApacheCompressor.getEntries(zf, "SNPS").size());
			totsnps+=snpsize.get(k);
			if(k>0){
				if(headers.get(k).length!=headers.get(0).length){
					System.err.println("problem with "+this.f[k]+" "+this.f[0]);
					System.err.println(Arrays.asList(headers.get(k))+"\n"+
							Arrays.asList(headers.get(0)));
					throw new RuntimeException("!!");
				}
				//if(samples.get(k).getRowDimension()!=samples.get(0).length) throw new RuntimeException("!!");
				if(!Arrays.asList(res).equals(Arrays.asList(samples_info.get(0))))  throw new RuntimeException("!!");
			}else{
				nocol = headers.get(0).length;
				norow = samples_info.get(0).length;
			}
			zf.close();
		}
		double sum =0;
		for(int k=0; k<f.length; k++){
			double v = snpsize.get(k)/totsnps;
			sum+=v;
			snpsize.set(k, v);
		}
		System.err.println(snpsize);
		System.err.println(sum);
	}
	
	
	
	public void calc() throws Exception{
		RealMatrix rm = new Array2DRowRealMatrix(samples.get(0).getRowDimension(), samples.get(0).getColumnDimension());
		double sum =0; double sum1 = 0;
		for(int k=0; k<samples.size(); k++){
			double v = samples.get(k).getEntry(1, 1);
			double v1 = snpsize.get(k);
			sum1+=v1*v;
			sum+=v1;
			System.err.println(f[k]+" "+v+" "+v1);
			RealMatrix rm1 = samples.get(k).scalarMultiply(v1);
			
			rm = rm.add(rm1);
		}
		System.err.println("final "+rm.getEntry(1,1));
		System.err.println(sum+" "+sum1);
		List<String> orig = this.original.get(0);
		PrintWriter pw = new PrintWriter(new FileWriter(this.out));
		String[] samps = this.samples_info.get(0);
		boolean[] num = getNumCodes(this.headers.get(0));
		for(int k=0; k<samps.length; k++){
			String[] origk = orig.get(k).split("\\t+");
			pw.print(samps[k]);
			for(int j=1; j<this.nocol; j++){
				pw.print("\t"+(num[j] ?  String.format("%5.3g",rm.getEntry(k, j)).trim() : origk[j]));
			}
			pw.println();//"\t"+orig.get(k).replace('\t', '_'));
		}
		pw.close();
	}

	boolean[] getNumCodes(String[] headers){
		boolean[] num = new boolean[headers.length];
		Arrays.fill(num, true);
		for(int k=0; k<headers.length; k++){
			if(headers[k].toLowerCase().equals("plate")) num[k] = false;
		}
		num[0] = false;
		return num;
	}
	private RealMatrix process(List<String> entries,String[]  headers) {
		RealMatrix rm = new Array2DRowRealMatrix(entries.size(),headers.length);
		boolean[] num = getNumCodes(headers);
		for(int i=0; i<entries.size(); i++){
			String[] str = entries.get(i).split("\\s+");
			for(int k=0; k<str.length; k++){
				rm.setEntry(i, k, num[k] ? Double.parseDouble(str[k]) : Double.NaN);
			}
		}
		return rm;
	}
	
}
