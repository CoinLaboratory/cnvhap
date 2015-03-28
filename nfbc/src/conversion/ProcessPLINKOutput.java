package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public class ProcessPLINKOutput {
	public static void main(String[] args){
		try{
			final File user = new File(System.getProperty("user.dir"));
			File plink = new File(user, "plink");
			File[] todo =plink.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.isDirectory();
				}
				
			});
			List l = new ArrayList();
			for(int i=0; i<todo.length; i++){
				final String tod = todo[i].getName();
				Callable call = new Callable(){
					
				public Object call(){
					try{
				ProcessPLINKOutput ppo = new ProcessPLINKOutput(user, 
						tod
						);
				ppo.run();
					}catch(Exception exc){
						exc.printStackTrace();
					}
				return null;
				}
				};
				l.add(call);
			}
			es.invokeAll(l);
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
		es.shutdown();
	}
	
	
	
	private static String[] read(File file) throws Exception {	
		BufferedReader br = new BufferedReader(new FileReader(file));
		List<String> l = new ArrayList<String>();
		String st = "";
		while((st = br.readLine())!=null){
			l.add(st);
		}
		return l.toArray(new String[0]);
	}



	ZipFile[] files;
	
	File outputdir;
	
	BufferedReader br;//map
	String line; //current line
	String[] el;
	String chrom;
	PrintWriter pw;
	List<String> header = new ArrayList<String>();
	int chr_index = 0;
	int id_index =1;
	int pos_index = 3;
	int hwd_index, maf_index, call_index;
	int[] p_index;
	String[][] zipFileHeader;
	String headerString;
	String formatString;
	//String formatString1;
	 File[] files2;
	 static boolean minimal = true;
	ProcessPLINKOutput(File f, final String str) throws Exception{
		 File outputdir1 = new File(f, "summary");
		 outputdir1.mkdir();
		 outputdir = new File(outputdir1, str);
		 outputdir.mkdir();
		 log = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputdir, "log.txt"))));
		// log1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputdir, "log1.txt"))));
		 File plink = new File(f, "plink");
		 File[] files1 = plink.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				if(pathname.isFile()) return false;
				if(str==null ) return true;
				
					if(str.equals(pathname.getName())) return true;
			
				
				return false;
			}
			 
		 });
		
		files2 = makeZIPFiles(files1);
		 br =new BufferedReader(new FileReader(new File(f,"join0.txt")));
		  files = new ZipFile[files2.length];
		 for(int i=0; i<files.length; i++){
			 files[i] = new ZipFile(files2[i]);
		 }
		 line = br.readLine();
		 el = line.split("\t");
		List l1 = Arrays.asList(el);
		chr_index = l1.indexOf("CHR");
		pos_index = l1.indexOf("POSITION");
		id_index = l1.indexOf("MARKER");
		hwd_index = l1.indexOf("P_HWD");
		maf_index = l1.indexOf("MAF");
		call_index = l1.indexOf("CALLRATE");
		 StringBuffer formatStr = new StringBuffer();
		// StringBuffer formatStr1 = new StringBuffer();
		 header.add("snpid");
		// header.add("loc");
		// header.add("summ_summ_MAF");
		// header.add("summ_summ_CALLRATE");
		// header.add("summ_summ_HWE");
		 start = 1;
		 formatStr.append("%7s");
		// formatStr1.append("%7s\t%7s");
		 p_index = new int[files.length];
		 zipFileHeader = new String[files.length][];
		 for(int i=0; i<zipFileHeader.length; i++){
			 zipFileHeader[i] =  Compressor.read(files[i],"MARKER").get(0).split("\\s+");
			 p_index[i] = Arrays.asList(zipFileHeader[i]).indexOf("P");
			 header.add(files1[i].getName().replace('_', '-')+"_P");
			 formatStr.append("\t%7s");
		//	 formatStr1.append("\t%7.3f");
		 }
		 formatString = formatStr.toString();
		 headerString = String.format(formatString, header.toArray());
		 toPrint = new Object[header.size()];
	}
	public void run() throws Exception{
		while(line!=null){
			//log1.println(line.substring(0,10));
			runInner();
		}
		for(int i=0; i<files.length; i++){
			files[i].close();
		}
		log.close();
	//	log1.close();
	}
	final int start;
	final Object[] toPrint;
	final PrintWriter log;
	public void runInner() throws Exception{
		//log1.println(line.substring(0,10));
		 chrom = el[chr_index];
		 File outp1 = new File(outputdir, chrom);
		 outp1.mkdir();
		 pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(outp1, "assoc.txt"))));
		 pw.println(headerString);
		 inner: while(el[chr_index].equals(chrom)){
			toPrint[0] = el[id_index];
			//toPrint[1] = el[pos_index];
			//toPrint[2] = el[this.maf_index];
		//	toPrint[3] = el[this.call_index];
		//	toPrint[4] = el[this.hwd_index];
			for(int i=0;i<files.length; i++ ){
				List<String> l1 = Compressor.read(files[i],el[id_index]);
				if(l1==null){
					log.println("WARNING DID NOT FIND SNP "+el[id_index] +" in "+files2[i]);
					toPrint[i+start] = "NA";
					/* line = br.readLine();
					 if(line==null) break inner;
					 el = line.split("\t");
					continue inner;*/
				}
				else toPrint[i+start] =  l1.get(0).split("\\s+")[this.p_index[i]];
			}
			pw.println(String.format(formatString, toPrint));
			 line = br.readLine();
			 if(line==null) break inner;
			 el = line.split("\t");
		}
		pw.close();
	}
	final static ExecutorService es  = Executors.newFixedThreadPool(4);;;
	private File[] makeZIPFiles(File[] files2) throws Exception{
		File[] res = new File[files2.length];
		List l  = new ArrayList();
		System.err.println(Arrays.asList(files2));
		for(int i=0; i<files2.length; i++){
			
			final File outp = new File(files2[i], "assoc.zip");
			res[i] = outp;
			File[] ass = files2[i].listFiles(new FileFilter(){
				
				public boolean accept(File pathname) {
					return pathname.getName().startsWith("plink.assoc");
				}
				 
			 });
			if(ass.length==0){
				//deleteNonZip(files2[i]);
				continue;
			}
			
			final  File assoc = ass[0];
			
			
			final File fi = files2[i];
			Callable call = new Callable(){
				public Object call(){
					try{
					
					
					
						makeZip(outp, assoc);
					deleteNonZip(fi);
					}catch(Exception exc){
						exc.printStackTrace();
					}
					return null;
				}
			};
			l.add(call);
			
		}
		es.invokeAll(l);
		return res;
		
	}

	private void deleteNonZip(File file) {
		File[] f = file.listFiles();
		for(int i=0; i<f.length; i++){
			if(!f[i].getName().endsWith("zip")) f[i].delete();
		}
		
	}

	static class WriteZip{
		 FileOutputStream dest;
		 CheckedOutputStream checksum;
		    ZipOutputStream outS;
		    OutputStreamWriter osw;
		String curr;
		String st;
		BufferedReader br;
		List<String> toWrite = new ArrayList<String>();
		String[] str;
		 List<String> l= new ArrayList<String>();
		WriteZip(File outp, File assoc) throws Exception{
			dest = Compressor.getOS(outp);
			  checksum = new   CheckedOutputStream(dest, new Adler32());
	           outS = new 
	          ZipOutputStream(new 
	            BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	         br = new BufferedReader(new FileReader(assoc));
	         st = br.readLine();
	         str = st.split("\\s+");
	         curr =  str[0];
		}
		public void writeList() throws Exception{
			 ZipEntry headings = new ZipEntry(curr);
	           outS.putNextEntry(headings);
	           for(int i=0; i<1; i++){
	        	   osw.write(l.get(i));
	        	   osw.write("\n");
	           }
	           osw.flush();
	           outS.closeEntry();
			   l.clear();
		}
		public void run() throws Exception{
			while(st!=null){
				runInner();
				writeList();
			}
			 outS.close();
			 
		}
		public void runInner()throws Exception{
				curr = str[0];
				while(curr.equals(str[0])){
					l.add(st);
			        st =  br.readLine();
			        if(st==null){
			        	return;
			        }
			        str = st.split("\\s+");
				}
				
		}
	}
	
	private void makeZip(File outp, File assoc) throws Exception{
		System.err.println("making"+outp+" "+assoc);
		WriteZip wz = new WriteZip(outp, assoc);
		wz.run();
		
	}

	public void process(){
		
	}
	
}
