package data.util;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;



public class CompressDir {
	/** Class used to write compressed file 
	* allows writing direct to zip, or to file, which is appeneded at end
	 * */
	public static void main(String[] args){
		try{
			//if(true) System.exit(0);
			 File dir1 = new File(System.getProperties().getProperty("user.dir"));
	        	CompressDir.compress(dir1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	 FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter osw;
	    
	    File inDir;
	    int len;
	    
	    public CompressDir(File f) throws Exception{
	    	this.inDir = f;
	    	inDir.mkdir();
	    	len = inDir.getAbsolutePath().length()+1;
	    	 dest = new FileOutputStream(new File(inDir.getParentFile(), inDir.getName()+".zip"));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
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
	           outS.closeEntry();
	    	}else{
	    		osw.close();
	    	}
	    }
	    
	    public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip)throws Exception{
	    	if(writeDirectToZip){
	    	ZipEntry headings = new ZipEntry(entry);
		    outS.putNextEntry(headings);
		    return osw;
	    	}else{
	    		return new OutputStreamWriter((new FileOutputStream(new File(inDir, entry))));
	    	}
		        
	    }

	    public static void compress(File dir) {
			try{
						(new CompressDir(dir)).run();
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
}
