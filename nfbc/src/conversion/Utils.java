package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.util.Compressor;

public class Utils {
 public static FileFilter zipFilter = new FileFilter(){
	@Override
	public boolean accept(File arg0) {
		return arg0.getName().endsWith("zip");
	}
	 
 };
 
 public static FileFilter dirFilter(final Set<String> s){
	 return new FileFilter(){
 
		@Override
		public boolean accept(File arg0) {
			return arg0.isDirectory() && s.contains(arg0.getName());
		}
		 
	 };
 }
 
 public static FilenameFilter dirFilter(final String pref, final String suff){
	 return new FilenameFilter(){
 
		@Override
		public boolean accept(File arg0, String nme) {
			return nme.startsWith(pref) && nme.endsWith(suff);
		}
		 
	 };
 }
 

public static  List<Integer> readPosInfo(File f, int index, boolean header) throws Exception{
        
        List<Integer> res = new ArrayList<Integer>();
       readPosInfo(f, new int[] {index}, header, new List[] {res}, new Class[] {Integer.class});
        return res;
    }
public static  List<String> readStringInfo(File f, int index, boolean header) throws Exception{
    
    List<String> res = new ArrayList<String>();
   readPosInfo(f, new int[] {index}, header, new List[] {res}, new Class[] {String.class});
    return res;
}
 static   class ZipBufferedReader extends BufferedReader {
 	BufferedReader br;
 	String st = "";
 	Enumeration<ZipEntry> en;
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
    public static BufferedReader getBufferedReader(File f) throws Exception{
        if(f.exists() && f.length()>0){
        	if(f.getName().endsWith(".gz")){
        		 return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
        	}
        	else if(f.getName().endsWith("zip")){
        		ZipFile zf = new ZipFile(f);
        		return new ZipBufferedReader(zf, zf.entries(), ConvertLong.filter,ConvertLong.header);
        		
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
    static Reader getNextReader(ZipFile zf, Enumeration<ZipEntry> en, String st) {
    	try{
    	
    	while(en.hasMoreElements()){
    		ZipEntry ent = en.nextElement();
    		System.err.println(ent.getName());
    		//ConvertLong.currentEntry = ent.getName().split("\\.")[0];
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
	  ZipEntry entry =  f.getEntry(string);
	  if(entry==null) return null;
       return   new BufferedReader(new InputStreamReader(
               f.getInputStream(entry)));
         
   }
	
	public static Iterator<String> getStringIterator(final File f) throws Exception{
		return getStringIterator(new BufferedReader(new FileReader(f)));
	}
	
	public static Iterator<String> getStringIterator(final File[] f) throws Exception{
		return new Iterator<String>(){
			int index =-1;
			ZipFile zf;
			BufferedReader br;
			String st = null;
			String suffix = null;
			{
				nextFile();
			}
			public boolean nextFile(){
				try{
				if(br!=null){
					br.close();
					zf.close();
				}
				index++;
				if(index < f.length){
					System.err.println("opening "+f[index]+".zip");
					zf = new ZipFile(f[index]+".zip");
					br = Compressor.getBufferedReader(zf,  "SNPS");
					suffix = "\t"+f[index].getParentFile().getName();
					st = br.readLine();
					return true;
				}
				else{
					zf = null;
					return false;
				}
				}catch(Exception exc){
					exc.printStackTrace();
					return false;
				}
				
			}
			@Override
			public boolean hasNext() {
				if(st==null) {
					return nextFile();
				}else return true;
			}

			@Override
			public String next() {
				String st1 = st+suffix;
				try{
				 st = br.readLine();
				}catch(Exception exc){
				 exc.printStackTrace();	
				}
				// TODO Auto-generated method stub
				return st1;
			}

			@Override
			public void remove() {
				// TODO Auto-generated method stub
				
			}
			
		};
	}
	
	public static Iterator<String> getStringIterator(final BufferedReader br) throws Exception{
		return new Iterator<String>(){
			String st = br.readLine();
			public boolean hasNext(){
				if(st==null) {
					try{br.close();}catch(Exception exc){exc.printStackTrace();}
				}
				return st!=null;
			}
			@Override
			public String next() {
				String st1 = st;
				try{
				st = br.readLine();
				}catch(Exception exc){
					exc.printStackTrace();
					st = null;
				}
				return st1;
			}
			@Override
			public void remove() {
				// TODO Auto-generated method stub
				
			}
		};
	}


	public static File[] listFiles(File dir, Map<String, String> chroms) {
		List<File>f = new ArrayList<File>();
		for(Iterator<Map.Entry<String, String>> it = chroms.entrySet().iterator(); it.hasNext();){
			Map.Entry<String, String> nxt = it.next();
			f.add(new File(dir, nxt.getValue()+"/"+nxt.getKey()));
			}
		return f.toArray(new File[0]);
	}

/*	public static Iterator<String[][]> getIterator(final Iterator<String>[] res) {
		final int len = res.length;
		return new Iterator<String[][]>(){
			String[][] resu = new String[len][];
			@Override
			public boolean hasNext() {
				return res[0].hasNext();
			}
			@Override
			public String[] next() {
				for(int k=0; k<len; k++){
					if(!res[k].hasNext()) return null;
					resu[k] = res[k].next();
				}
				return resu;
			}
			@Override
			public void remove() {}
		};
	}*/
	
}
