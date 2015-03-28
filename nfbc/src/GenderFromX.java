import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipFile;

import conversion.Compressor;


public class GenderFromX {
	public static void main(String[]args){
		try{
			ZipFile f = new ZipFile(new File("X.zip"));
			List<String> snps =Compressor.read(f, "SNPS");
			List<String> samples =Compressor.read(f, "Samples");
			String[] res = new String[samples.size()];
			int[] cnt = new int[res.length];
			Arrays.fill(cnt,0);
			for(int i=0; i<snps.size(); i++){
				Compressor.readZip(f, snps.get(i), res, 0, "\t");
				for(int i1=0; i1<res.length; i1++){
					char[] c = res[i1].toCharArray();
					if(c[0]!=c[1]) cnt[i1]++;
				}
			}
			
			
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("out.txt")));
			Set<String > done = new HashSet<String>();
			for(int i=0; i<samples.size(); i++){
				double d = ((double)cnt[i])/(double)snps.size();
				if(!done.contains(samples.get(i))){
				pw.println(samples.get(i)+" "+ (d< 0.05 ? "Male" : "Female"));
					done.add(samples.get(i));
				}
			}
			pw.close();
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
}
