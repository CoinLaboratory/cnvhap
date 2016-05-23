package lc1.util;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipException;
import java.util.zip.ZipOutputStream;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;



public class CompressDirApache {
	/** Class used to write compressed file 
	* allows writing direct to zip, or to file, which is appeneded at end
	 * */
	public static void main(String[] args){
		try{
			//if(true) System.exit(0);
			 File dir1 = new File(System.getProperties().getProperty("user.dir"));
	        	CompressDirApache.compress(dir1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	 FileOutputStream dest;
	   CheckedOutputStream checksum;
	    ZipArchiveOutputStream outS;
	    OutputStreamWriter osw;
	    
	    public File inDir;
	    int len;
	    
	    public CompressDirApache(File f1) throws Exception{
	    	//if(f1.exists()) throw new RuntimeException("!!");
	    	String nme = f1.getName();
	    	File f = nme.endsWith("zip") ? new File(f1.getParentFile(),nme.substring(0,nme.length()-4))  : f1;
	    	this.inDir = f;
	    	inDir.mkdir();
	    	len = inDir.getAbsolutePath().length()+1;
	    	 dest = Compressor.getOS(new File(inDir.getParentFile(), inDir.getName()+".zip"));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipArchiveOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	    }
	    public  boolean delete = true;
	    
	    public void close() throws Exception{
	    	this.outS.close();
	    }
	    public void run() throws Exception{
	    	try{
	    	File[] f = inDir.listFiles();
	    	for(int i=0; i<f.length; i++){
	    	this.writeHeader(f[i]);
	    	this.delete(f[i]);
	    	}
	    	this.delete(inDir);
	    	this.outS.close();
	    	}catch(Exception exc){
	    		System.err.println("problem with "+inDir);
	    		exc.printStackTrace();
	    	}
	    }
	    public void delete(File f) throws Exception{
	    	if(!delete) return;
	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			delete(f1[i]);
	    		}
	    	}
	    	f.delete();
	    }
	    public void writeHeader(File f) throws Exception{
	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			writeHeader(f1[i]);
	    		}
	    	}
	    	else{
	    		OutputStreamWriter osw1 = this.getWriter(f.getAbsolutePath().substring(len), true);
	    		BufferedReader br = new BufferedReader(new FileReader(f));
	            String str = "";
	            while((str = br.readLine())!=null){
	        	  osw.write(str);osw.write("\n");
	           }
	          this.closeWriter(osw1);
	           br.close();
	    	}
	    }
	    
	    
	    
	    public void closeWriter(OutputStreamWriter osw) throws Exception{
	    	if(osw==this.osw){
	    	   osw.flush();
	           outS.closeArchiveEntry();
	           outS.close();
	    	}else{
	    		osw.close();
	    	}
	    }
	    
	    public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip, String comment)throws ZipException, IOException{
	    	if(writeDirectToZip){
	    		
	    	ZipArchiveEntry headings = new ZipArchiveEntry(entry);
	    	headings.setComment(comment);
		    outS.putArchiveEntry(headings);
		    return osw;
	    	}else{
	    		return new OutputStreamWriter((new FileOutputStream(new File(inDir, entry),true)));
	    	}
		        
	    }
	    public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip)throws Exception{
	    	if(writeDirectToZip){
	    	ZipArchiveEntry headings = new ZipArchiveEntry(entry);
		    outS.putArchiveEntry(headings);
		    return osw;
	    	}else{
	    		return new OutputStreamWriter((new FileOutputStream(new File(inDir, entry),true)));
	    	}
		        
	    }

	    public static void compress(File dir) {
			try{
						(new CompressDirApache(dir)).run();
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		public void copy(ZipFile zf, String name) throws Exception {
			ZipArchiveEntry e = zf.getEntry(name);
			OutputStreamWriter os = this.getWriter(name, false, e.getComment());
			BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(e)));
			String st = "";
			while((st = snps.readLine())!=null){
				os.write(st);
				os.write("\n");
			}
			snps.close();
			this.closeWriter(os);
		}
}
