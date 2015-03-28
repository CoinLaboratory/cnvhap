

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;



public class CompressDir {
	
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
	    
	    CompressDir(File f) throws Exception{
	    	this.inDir = f;
	    	len = inDir.getAbsolutePath().length()+1;
	    	 dest = getOS(new File(inDir.getParentFile(), inDir.getName()+".zip"));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	    }
	    public static FileOutputStream getOS(File f) throws Exception{
	    	 // if(f.exists() && f.length()>0) throw new RuntimeException("!!" +f.getAbsolutePath());
	    	    return new FileOutputStream(f);
	    	}
	    public void run() throws Exception{
	    	File[] f = inDir.listFiles();
	    	for(int i=0; i<f.length; i++){
	    	this.writeHeader(f[i]);
	    	this.delete(f[i]);
	    	}
	    	this.delete(inDir);
	    	this.outS.close();
	    }
	    public static void delete(File f) throws Exception{
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
	    	else if(f.exists() && !f.getName().startsWith(".")){
	        ZipEntry headings = new ZipEntry(f.getAbsolutePath().substring(len));
	        BufferedReader br = new BufferedReader(new FileReader(f));
	        outS.putNextEntry(headings);
	          String str = "";
	          while((str = br.readLine())!=null){
	        	  osw.write(str);osw.write("\n");
	          }
	            osw.flush();
	           outS.closeEntry();
	           br.close();
	    	}
	    	
	           
	    }

		
		
		public static void compress(File dir) {
			try{
			
			File[] f = dir.listFiles(new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory();
				}
				
			});
			List l = new ArrayList();
			for(int i=0; i<f.length; i++){
				final File fi = f[i];
				Callable call = new Callable(){

					@Override
					public Object call() throws Exception {
						(new CompressDir(fi)).run();
						return null;
					}
					
				};
				l.add(call);
				
			}
			for(int i=0; i<l.size(); i++){
				((Callable)l.get(i)).call();
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
}
