package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;

import lc1.util.CompressDir;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipFile;

public class MergeZip {

	public static void main(String[] args){
		try{
			//if(true) System.exit(0);
			File f = new File( System.getProperty("user.dir"));
			String outf ="";
		
			File[] fall = new File[args.length];
			for(int k=0; k<args.length; k++){
				outf = outf+args[k]+"_";
				fall[k] = new File(f, args[k]+"/all");
			}
			final  File out = new File(f, outf);
			out.mkdir();
			mergeinternal(new File(out,"all"), Arrays.asList(fall), false);
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}

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
	System.err.println(nme.getAbsolutePath());
	ZipFile[] zf = new ZipFile[zips.size()];
	System.err.println(zf.length);
	for(int k=0; k<zf.length; k++){
		System.err.println(zips.get(k).getAbsolutePath()+".zip");
		File zipf = new File(zips.get(k).getAbsolutePath()+".zip");
	
		zf[k] = new ZipFile(zipf);
		System.err.println("done opening "+k);
		if(removeInternal)zipf.deleteOnExit();
	}
	System.err.println("here");
	List<String> nomerge = Arrays.asList("Name:SNPS".split(":"));
	Enumeration en= zf[0].getEntries();
	while(en.hasMoreElements()){
		ZipArchiveEntry ze = (ZipArchiveEntry)en.nextElement();
		boolean append = ze.getName().equals("Samples");
	//	System.err.println(ze.getName());
		OutputStreamWriter osw = dir.getWriter(ze.getName(), true);
		for(int k=0; k<zf.length; k++){
			//System.err.println(k);
			if(k==0 || ! nomerge.contains(ze.getName())){
				String nme1 = ze.getName();
				ZipArchiveEntry entry = zf[k].getEntry(nme1);
				if(entry==null){
					System.err.println(zips.get(k).getAbsolutePath()+".zip not contained "+nme1);
				}
				BufferedReader br = new BufferedReader(new InputStreamReader(zf[k].getInputStream(entry)));
				String st = "";
				while((st=br.readLine())!=null){
					if(append) st = zips.get(k).getParentFile().getName()+"."+st;
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
