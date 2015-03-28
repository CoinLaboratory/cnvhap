package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;
import java.util.zip.ZipFile;

public class Uncompress {
	public static void main(String[] args){
		try{
		//	if(true) return;
			String[] str = //"11:12:13:14:15:16:17:18:19:20:21:22"
			"10:20:21:22:X:XY:M"
				.split(":");
			File dir = 
			new File(System.getProperties().getProperty("user.dir"));
			 
			for(int k=0; k<str.length; k++){
				Uncompress unc=  new Uncompress(dir,str[k]);
				unc.writeAll();
				(new CompressDir1(unc.dir)).run();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public Uncompress(File dir, String name)throws Exception{
		this.zf = new ZipFile(new File(dir, name+".zip"));
		this.dir = new File(dir, name);
		this.dir.mkdir();
	}
	final File dir;
	
	final ZipFile zf;
	public List<String> write(String name) throws Exception
	{
		List<String> l = Compressor.read(zf, name);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name))));
		for(int i=0; i<l.size(); i++){
			pw.println(l.get(i));
		}
		pw.close();
		return l;
	}
	
	public List<String> write1(String name) throws Exception
	{
		List<String> l = Compressor.read(zf, name);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name))));
		for(int i=0; i<l.size(); i++){
			pw.println(l.get(i).replaceFirst("\t", ""));
		}
		pw.close();
		return l;
	}
	public List<String> write2(String name) throws Exception
	{
		List<String> l = Compressor.read(zf, name);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name))));
		String st1 = l.get(0);
		int ind = st1.indexOf('\t',st1.indexOf('\t')+1);
		
		pw.println("Genotype"+st1.substring(ind));
		for(int i=1; i<l.size(); i++){
			pw.println(l.get(i));
		}
		pw.close();
		return l;
	}
	
	public void writeAll() throws Exception{
		write2("Name");
		write("Samples");
		List<String> snps = write("SNPS");
		for(int k=0; k<snps.size(); k++){
			write1(snps.get(k).split("\t")[3]);
		}
		
		
	}
}
