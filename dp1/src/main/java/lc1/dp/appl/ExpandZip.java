package lc1.dp.appl;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class ExpandZip {
	
	public static void main(String[] args){
		File user = new File( System.getProperty("user.dir"));
		final String date  = args[0];
		if(false)main(user);
		else{File[] f = user.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.isDirectory() && pathname.getName().endsWith(date);
			}
			
		});
		for(int i=0; i<f.length; i++){
			main(f[i]);
		}
		
	   	     //    VirtualDataCollection.main(new String[] {"1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22"});
	   	     //    Regression.es.shutdown();
		}
	}
	
	public static void main(File dirF){
		try{
//			File dirF =new File( System.getProperty("user.dir"));
			File[] f = dirF.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.getName().endsWith("zip");
				}
				
			});
			for(int i=0; i<f.length; i++){
				ZipFile zf = new ZipFile(f[i]); 
				 String[] resEntry = findEntry(zf, "res");
				 
	               File dirOut = new File(dirF, f[i].getName().substring(0, f[i].getName().length()-4));
	               if(dirOut.exists() && numFiles(dirOut)>0 && false) continue;
	               dirOut.mkdir();
	               for(int j=0; j<resEntry.length; j++){
	                  
	                  
	                   expand(dirOut, zf,  resEntry[j]);
	               }
				//expand(dirF, zf, "res");
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	private static int numFiles(File dirOut) {
		if(dirOut.isDirectory()){
			File[] f = dirOut.listFiles();
			int num=0;
			for(int i=0; i<f.length; i++){
				num+=numFiles(f[i]);
			}
			return num;
		}
		else return 1;
	}

	public static String[] findEntry(ZipFile zf, String string) {
    	List<String> l = new ArrayList<String>();
    	
        for(Enumeration en = zf.entries(); en.hasMoreElements();){
           ZipEntry ent =  (ZipEntry) en.nextElement();
           String nme = ent.getName();
           if(nme.indexOf(string)>=0 ){//&& nme.length() > nme1.length() ){
        	   String[] spl = nme.split("/");
        	   
        	   if( spl.length>=3 && nme.indexOf("geno")>=0 && nme.endsWith("zip"))
               //if(!spl[0].startsWith(spl[2]) || spl[2].indexOf('_')<0) 
            	   l.add(nme);// nme1 = nme;
           }
        }
        return l.toArray(new String[0]);
      }
	
	
	public static File expand(File dir, ZipFile zf,  String resEntry) throws Exception {
        // TODO Auto-generated method stub
        String[] f = resEntry.split("/");
        File par = dir;
        for(int i=1; i<f.length-3; i++){
            par = new File(par, f[i]);
            par.mkdir();
        }
        String name = f[f.length-3];
        if(name.indexOf("_")>0){
        	//String[] str1 = name.split("_");
//        	String[] str2 = name.split("\\.");
        	name = name.split("_")[0];//+"."+name.split("\\.")[1];
        }
        File outF = new File(par, name+".zip");
        if(outF.exists() && outF.length()>0) return outF;
        BufferedInputStream in = new BufferedInputStream(zf.getInputStream(zf.getEntry(resEntry)));
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outF));
        byte[] buf = new byte[1024];
        int len = 0;
        int reallen = 0;

        //System.out.println(file+":"+getLocalPath()+outfile);
        while(true)
        {
            len = in.read(buf);

            //System.out.print(".");
            if(len == StreamTokenizer.TT_EOF)
            {
                break;
            }

            out.write(buf, 0, len);
            reallen += len;

        }
        in.close();
        out.close();
       return  outF;
    }
}
