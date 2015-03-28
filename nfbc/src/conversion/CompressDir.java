package conversion;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import lc1.util.Compressor;

public class CompressDir {
	/** Class used to write compressed file 
	* allows writing direct to zip, or to file, which is appeneded at end
	 * */
	
	
	    FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter1 osw_zip;
	    
	    File inDir;
	    int len;
	    
	    public CompressDir(File f) throws Exception{
	    	this.inDir = f;
	    	if(!inDir.getParentFile().exists()){
	    		inDir.getParentFile().mkdir();
	    	}
	    	inDir.mkdir();
	    	len = inDir.getAbsolutePath().length()+1;
	    	 dest = Compressor.getOS(new File(inDir.getParentFile(), inDir.getName()+".zip"));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw_zip =  new OutputStreamWriter1(outS);
	        	
	         outS.setMethod(ZipOutputStream.DEFLATED);
	    }
	    public  boolean delete = true;
	    public boolean writeToZipAsYouGo=false;
	    
	    public void close() throws Exception{
	    	File[] f = inDir.listFiles();
	    	for(int k=0; k<f.length; k++){
	    		writeToZip(f[k], null);
	    	}
	    	this.delete(inDir);
	    	this.outS.close();
	    }
	    
	    private void writeToZip(File f, String comment) throws IOException {
			String entry = f.getName();
			ZipEntry headings = new ZipEntry(entry);
	    	if(comment!=null) headings.setComment(comment);
		    outS.putNextEntry(headings);
			BufferedReader br = new BufferedReader(new FileReader(f));
            String str = "";
            while((str = br.readLine())!=null){
        	  osw_zip.write(str);osw_zip.write("\n");
           }
		    osw_zip.flush();
		    outS.closeEntry();
		    delete(f);
			
		}
	   

	    public void delete(File f) throws IOException{
	    	if(!delete) return;

	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			delete(f1[i]);
	    		}
	    	}
	    	f.delete();
	    }
	   
	    
	   
	    public class OutputStreamWriter1 {
	    	OutputStreamWriter osw;
	    	File f;
	    	String comment;
	    	private int curr_index=0;
	    	private final String name;
	    	public String toString(){
	    		return name;
	    	}
	    	public int curr_index(){
	    		return curr_index;
	    	}

	    	public OutputStreamWriter1(ZipOutputStream outS) {
	    		 osw = new OutputStreamWriter(outS);
	    		 f = null;
	    		 name = "zipfile";
			}
			public OutputStreamWriter1(File f, String comment) throws IOException {
				osw = new OutputStreamWriter((new FileOutputStream(f)));
				this.f = f;
				this.comment=comment;
				this.name = f.getName();
			}
			public  void printLine(String[] str, int[] cols) throws IOException{
	    			osw.write(str[cols[0]]);
	    		for(int k1=1; k1<cols.length && cols[k1] < str.length; k1++){
	    			osw.write("\t"+str[cols[k1]]);
	    		}
	    		osw.flush();
	    		osw.write("\n");
	    		curr_index=curr_index+1;
	    	}
	    	public void printLineSNP(String[] str, int end, int[] cols) throws IOException{
	    		osw.write(str[cols[0]]);
	    	for(int k1=1; k1<cols.length && cols[k1] < str.length; k1++){
	    		if(cols[k1]<0)osw.write("\t"+end);
	    		else osw.write("\t"+str[cols[k1]]);
	    	}
	    	osw.flush();
	    	osw.write("\n");
	    	curr_index=curr_index+1;
    		
	    }
	    	
	    	public void write(String st) throws IOException{
	    		osw.write(st);
	    	
	    	}
	    	public void close() throws IOException{
	    		if(f==null){
	    		   osw.flush();
	 	           outS.closeEntry();
	    		}else{
	    			osw.close();
	    			if(writeToZipAsYouGo){
		    			writeToZip(f,comment);
	    			}
	    		}
	    	}
			
			public void flush() throws IOException {
				osw.flush();
				
			}
			public void println(String string) throws IOException {
				osw.write(string);osw.write("\n");
				curr_index=curr_index+1;
				
			}
	    }
	    
	    public OutputStreamWriter1 getWriter(String entry, boolean writeDirectToZip, String comment)throws Exception{
	    	if(writeDirectToZip){
	    	 ZipEntry headings = new ZipEntry(entry);
	    	 if(comment!=null) headings.setComment(comment);
		     outS.putNextEntry(headings);
		     return osw_zip;
	    	}else{
	    		File f = new File(inDir, entry);
	    		return new OutputStreamWriter1(f,comment);
	    		
	    	}
		        

	    }
	    public OutputStreamWriter1 getWriter(String entry, boolean writeDirectToZip)throws Exception{
	      return this.getWriter(entry, writeDirectToZip,null);
		        

	    }

	   
}