package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipFile;

public class ExtractOverallVariance {
	 public static String dir = System.getProperty("user.dir");
	 public static String prefix = "";
	  public static String[] chr;// = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22".split(":");
public static void main(String[] args){
	try{
		chr = (new File(dir)).list(new FilenameFilter(){

			@Override
			public boolean accept(File dir, String name) {
			return name.endsWith(".zip") && ! name.startsWith("M") && !name.startsWith("X") && !name.startsWith("Y")&& !name.startsWith("p");
			}
			
		});
		for(int k=0; k<chr.length; k++){
			chr[k] = chr[k].substring(0, chr[k].indexOf('.'));
		}
		run(new File(dir));
	/*	File [] f = (new File(dir)).listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.getName().startsWith(prefix) && pathname.isDirectory();
			}
			
		});
		for(int i=0; i<f.length; i++){
			run(f[i]);
		}*/
	
		
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
public static void run(File base) throws Exception{
	PrintWriter res = new PrintWriter(new BufferedWriter(new FileWriter(new File(base,"variance.txt"))));
	double[] num   = new double[chr.length];
	int var_id = 1;
	{
		File f1 = new File(base, chr[0]+".zip");
		System.err.println("opening "+f1.getAbsolutePath());
		ZipFile zf = new ZipFile(f1);
		var_id = Arrays.asList(Compressor.read(zf,"Name").get(2).split("\\s+")).indexOf("variance");
		zf.close();
	}
	double sum=0;
	for(int i=0; i<chr.length; i++){
		File f1 = new File(base, chr[i]+".zip");
		if(!f1.exists() || f1.length()==0) continue;
		ZipFile zf = new ZipFile(f1);
		List<String>snps = new ArrayList<String>();
		Compressor.read(zf, "SNPS", snps, 0);
		num[i] = snps.size();
		sum+=num[i];
		zf.close();
	}
	
	List<String>indiv =get(chr[0], "Samples",0, base);
	double[] sums = new double[indiv.size()];
	Arrays.fill(sums, 0);
	for(int i=0; i<chr.length; i++){
		num[i] = num[i]/sum;
		File f = new File(base, chr[i]+".zip");
		if(!f.exists()) continue;
		ZipFile zf = new ZipFile(f);
		List<String>snps = new ArrayList<String>();
		Compressor.read(zf, "Samples", snps,var_id);
		for(int k=0; k<sums.length; k++){
			sums[k] +=num[i]*Double.parseDouble(snps.get(k));
		}
		zf.close();
	}
	res.println("PATIENT\tVariance");
	for(int i=0; i<indiv.size(); i++){
		res.println(indiv.get(i).split("#")[0]+"\t"+sums[i]);
	}
	res.close();
}
private static List<String> get(String chr, String header, int col, File base) throws Exception{
	List<String> l = new ArrayList<String>();
		ZipFile zf = new ZipFile(new File(base, chr+".zip"));
		Compressor.read(zf, header, l, col);
	return l;
}
}
