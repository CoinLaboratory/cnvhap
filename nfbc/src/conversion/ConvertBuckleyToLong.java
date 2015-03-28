package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

public class ConvertBuckleyToLong {

	public static void main(String[] args){
		try{
			File dir = new File( System.getProperty("user.dir"));
			String[] names = dir.list(new FilenameFilter(){

				@Override
				public boolean accept(File dir, String name) {
					return name.startsWith("ids");
				}
				
			});
			for(int i=0; i<names.length; i++){
				ConvertBuckleyToLong cbtl = new ConvertBuckleyToLong(dir, names[i].substring(names[i].indexOf('.')+1));
				cbtl.run();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	/** replaces file name with converted file for next stage */
	public static void main(File[] file){
		try{
			File[] res = new File[file.length];
			for(int i=0; i<file.length; i++){
			File dir = file[i].getParentFile();
			String name = file[i].getName();
			ConvertBuckleyToLong cbtl = new ConvertBuckleyToLong(dir, name.substring(name.indexOf('.')+1));
			cbtl.run();
			file[i] = cbtl.outf;
			}
			
		}catch(Exception exc){
			exc.printStackTrace();
		
		}
	}
	final File outf;
	List<String> snps = new ArrayList<String>();
	BufferedReader likelihoods;
	PrintWriter out;
	String header;
	public ConvertBuckleyToLong(File dir, String id) throws Exception{
		System.err.println(id);
		BufferedReader br = Utils.getBufferedReader(new File(dir, "ids."+id));
		String st = "";
		this.outf = new File(dir, "long."+id+".gz");
		out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf)))));
		
		while((st=br.readLine())!=null){
			this.snps.add(st);
		}
		br.close();
		this.likelihoods = new BufferedReader(new FileReader(new File(dir, "likelihoods."+id)));
		this.id = likelihoods.readLine();
		 header = likelihoods.readLine();
		out.println("Sample ID\tChr\tPosition\tSNP Name\t"+header);
		
	}
	String id;
	public void run() throws Exception{
		while(id!=null){
			getIndiv();
			id = likelihoods.readLine();
			likelihoods.readLine();
		}
		likelihoods.close();
		out.close();
	}
	public void getIndiv() throws Exception{
		System.err.println("indiv "+id);
		for(int i=0; i<snps.size(); i++){
			out.print(id);
			out.print("\t");
			out.print(snps.get(i));
			out.print("\t");
			out.println(likelihoods.readLine());
		}
		
		
	}
	
}
