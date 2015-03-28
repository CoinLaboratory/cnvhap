package data.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipFile;
import org.apache.commons.math.linear.RealMatrix;

public class ApacheUtil {
 public static  List<Integer> readPosInfo(File f, int index, boolean header) throws Exception{
        
        List<Integer> res = new ArrayList<Integer>();
       readPosInfo(f, new int[] {index}, header, new List[] {res}, new Class[] {Integer.class});
        return res;
    }
 static   class ZipBufferedReader extends BufferedReader {
 	BufferedReader br;
 	String st = "";
 	Enumeration<ZipArchiveEntry> en;
 	String filter;
 	ZipFile zf;
 	int head;
	int cnt=0;	
 	public ZipBufferedReader(ZipFile zf,  Enumeration en, String filter, int head) {
			super(getNextReader(zf,en, filter));
			try{
				
				for(int k=0; k<head; k++){
					st = this.readLineS();
				}
			}catch(Exception exc){
				exc.printStackTrace();
			}
			this.filter = filter;
			this.head = head;
			this.en = en;
			this.zf = zf;
			br = this;
			try{
			st = readLineS();
			}catch(Exception exc){
				exc.printStackTrace();
			}
			// TODO Auto-generated constructor stub
		}
		private String readLineS() throws IOException{
			// TODO Auto-generated method stub
			return super.readLine();
		}
		@Override
		public String readLine() throws IOException{
			
			String res = st;
			if(br==null) st =  null;
			else{
			st = br==this ?super.readLine(): br.readLine();
			if(st==null){
				if(br==this)  this.closeS();
				else br.close();
				cnt=0;
				Reader r = getNextReader(zf,en, filter);
				if(r!=null){
					br = new BufferedReader(r);
					for(int k=0; k<=head; k++){
						st = br.readLine();
					}
					//System.err.println(st);
				}
				else{
					br = null;
				}
			}
			}
			cnt++;
			return res;
		}
		public void closeS() throws IOException{
			super.close();
		}
		public void close() throws IOException{
			if(br==this)super.close();
			else br.close();
			this.zf.close();
		}
 	
 }
 
 static String filter = "";
 static int header = -1;
    public static BufferedReader getBufferedReader(File f) throws Exception{
        if(f.exists() && f.length()>0){
        	if(f.getName().endsWith(".gz")){
        		 return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
        	}
        	else if(f.getName().endsWith("zip")){
        		ZipFile zf = new ZipFile(f);
        		return new ZipBufferedReader(zf, zf.getEntries(),filter,header);
        		
        	}
            return new BufferedReader(new FileReader(f));
        }
        else{
            File f1 = new File(f.getAbsolutePath()+".gz");
            if(f1.exists())
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f1))));
            else return null;
        }
    }
    static Reader getNextReader(ZipFile zf, Enumeration<ZipArchiveEntry> en, String st) {
    	try{
    	String currentEntry = "";
    	while(en.hasMoreElements()){
    		ZipArchiveEntry ent = en.nextElement();
    		System.err.println(ent.getName());
    		currentEntry = ent.getName().split("\\.")[0];
    		if(ent.getName().indexOf(st)>=0){
    			return new InputStreamReader(zf.getInputStream(ent));
    		}
    	}
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    	return null;
    }
 
    
    public static  void readPosInfo(File f, int[] index, boolean header, List[] res, Class[] cl) throws Exception{
        if(f.exists() && f.length()>0)
        readPosInfo(getBufferedReader(f), index, header, res, cl);
    }
    
    public static  void readPosInfo(BufferedReader br, int[] index, boolean header, List[] res, Class[] cl) throws Exception{
        readPosInfo(br, index, header, res, cl, "\\s+");
    }
    public static  void readPosInfo(BufferedReader br, int[] index, boolean header, List[] res, Class[] cl, String spl) throws Exception{
        if(header) br.readLine();
        String st = "";
      
        while((st = br.readLine())!=null){
            String st1 = st.trim();
            String[] str = st1.split(spl);
         //   System.err.println(st1);
            for(int i=0; i<index.length; i++){
                if(str.length>index[i]){
                    try{
                        if(cl[i].equals(String.class)){
                            res[i].add(str[index[i]]);
                        }
                        else{
                            res[i].add(
                                    cl[i].getConstructor(new Class[] {String.class}).newInstance(new Object[] {str[index[i]]}));
                        }
                    }catch(Exception exc){
                        System.err.println(Arrays.asList(str));;
                        exc.printStackTrace();
                        System.exit(0);
                    }
                }
            }
        }
        br.close();
    }
    public static void copy(File in, File out)  throws Exception{
      BufferedReader br = new BufferedReader(new FileReader(in));
      PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
      String st = "";
      while((st = br.readLine())!=null){
          pw.println(st);
      }
      br.close();
      pw.close();
    }
    public static void delete(File file) {
      //  if(file.getName().startsWith("WG") || file.getName().endsWith("zip")) throw new RuntimeException("!!");
       if(file.isDirectory()){
           File[] f = file.listFiles();
           for(int i=0; i<f.length; i++){
               delete(f[i]);
           }
       }
       file.delete();
       
    }
	public static BufferedReader[] getBufferedReader(File[] f) throws Exception {
		BufferedReader[] res = new BufferedReader[f.length];
		for(int i=0; i<res.length;  i++){
			res[i] = getBufferedReader(f[i]);
		}
		return res;
	}
	public static BufferedReader getBufferedReader(ZipFile f, String string) throws Exception{
	  ZipArchiveEntry entry =  f.getEntry(string);
	  if(entry==null) return null;
       return   new BufferedReader(new InputStreamReader(
               f.getInputStream(entry)));
         
   }
	public static List<String> read(ZipFile pci, String string, int i) {
		List<String>res = new ArrayList<String>();
		try{
		
		BufferedReader br = getBufferedReader(pci, string);
		String st = "";
		while((st = br.readLine())!=null){
			res.add(st.split("\\s+")[i]);
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		return res;
	}
	public static void read(ZipFile pci, String string, RealMatrix mat) {
	
		try{
		
		BufferedReader br = getBufferedReader(pci, string);
		String st = "";
		for(int i=0; (st = br.readLine())!=null; i++){
			String[] str = st.split("\\s+");
			for(int j=0; j<str.length; j++){
				mat.setEntry(i, j, Double.parseDouble(str[j]));
			}
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	public static List<String> read(ZipFile pci, String string) {
		List<String>res = new ArrayList<String>();
		try{
		
		BufferedReader br = getBufferedReader(pci, string);
		String st = "";
		while((st = br.readLine())!=null){
			res.add(st);
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		return res;
	}
}
