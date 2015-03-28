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
import java.util.zip.ZipOutputStream;


public class ProcessPLINKOutput1 {
	public static void main(String[] args){
		try{
			File user = new File(System.getProperty("user.dir"));
			File[] todo =(new File(user, "summary")).listFiles(new java.io.FileFilter(){

				public boolean accept(File pathname) {
					return pathname.isDirectory();
				}
				
			});
			
			ProcessPLINKOutput1 ppo = new ProcessPLINKOutput1(user, 
					todo
					);
			ppo.run();
			
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



	BufferedReader[] files;
	
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
	 final File f;
	ProcessPLINKOutput1(File f, final File[] str) throws Exception{
		 outputdir = new File(f, "output");
		 outputdir.mkdir();
		 this.f = f;
		 
		
		 this.files2 = str;
		 log = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputdir, "log.txt"))));
		 log1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputdir, "log1.txt"))));
		
		
		//files2 = makeZIPFiles(files1);
		 br =new BufferedReader(new FileReader(new File(f,"join0.txt")));
		 
		 line = br.readLine().trim();
		 el = line.split("\t");
		List l1 = Arrays.asList(el);
		chr_index = l1.indexOf("CHR");
		pos_index = l1.indexOf("POSITION");
		id_index = l1.indexOf("MARKER");
		hwd_index = l1.indexOf("P_HWD");
		maf_index = l1.indexOf("MAF");
		call_index = l1.indexOf("CALLRATE");
		 files = new BufferedReader[str.length];
		 for(int i=0; i<files.length; i++){
			 files[i] = new BufferedReader(new FileReader(new File(new File(str[i], el[chr_index]), "assoc.txt")));
		 }
		 StringBuffer formatStr = new StringBuffer();
		// StringBuffer formatStr1 = new StringBuffer();
		 header.add("snpid");
		 header.add("loc");
		 header.add("summ_summ_MAF");
		 header.add("summ_summ_CALLRATE");
		 header.add("summ_summ_HWE");
		 start = 5;
		 formatStr.append("%7s\t%7s\t%7s\t%7s\t%7s");
		// formatStr1.append("%7s\t%7s");
		 p_index = new int[files.length];
		 zipFileHeader = new String[files.length][];
		 for(int i=0; i<zipFileHeader.length; i++){
			
			 zipFileHeader[i] =  files[i].readLine().trim().split("\\s+");
			 p_index[i] = Arrays.asList(zipFileHeader[i]).indexOf(str[i].getName().replace('_', '-')+"_P");
			 if(p_index[i] <0){
				String[] st = zipFileHeader[i];
				 throw new RuntimeException("!!");
			 }
			 header.add(str[i].getName().replace('_', '-')+"_P");
			 formatStr.append("\t%7s");
		//	 formatStr1.append("\t%7.3f");
		 }
		 formatString = formatStr.toString();
		// formatString1 = formatStr1.toString();
		 headerString = String.format(formatString, header.toArray());
		 toPrint = new Object[header.size()];
		
		 
	}
	public void run() throws Exception{
		boolean first = true;
		while(line!=null){
			if(!first){
				 for(int i=0; i<files.length; i++){
					 files[i].close();
					 files[i] = new BufferedReader(new FileReader(new File(new File(files2[i],el[chr_index]), "assoc.txt")));
					 files[i].readLine().trim();
				 }
			}
			log1.println(line.substring(0,10));
			runInner();
			first = false;
		}
		br.close();
		log.close();
		log1.close();
	}
	final int start;
	final Object[] toPrint;
	final PrintWriter log,log1;
	public void runInner() throws Exception{
		log1.println(line.substring(0,10));
		 chrom = el[chr_index];
		 File outp1 = new File(outputdir, chrom);
		 outp1.mkdir();
		 pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(outp1, "assoc.txt"))));
		 pw.println(headerString);
			
		inner: while(el[chr_index].equals(chrom)){
			toPrint[0] = el[id_index];
			toPrint[1] = el[pos_index];
			toPrint[2] = el[this.maf_index];
			toPrint[3] = el[this.call_index];
			toPrint[4] = el[this.hwd_index];
			for(int i=0;i<files.length; i++ ){
				String st1 = files[i].readLine();
				//if(i==0) System.err.println(st1);
				String  l1 = st1.trim();//el[id_index]);
			    toPrint[i+start] =  l1.split("\\s+")[this.p_index[i]];
			}
			pw.println(String.format(formatString, toPrint));
			 line = br.readLine().trim();
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
	         st = br.readLine().trim();
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
			        st =  br.readLine().trim();
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
