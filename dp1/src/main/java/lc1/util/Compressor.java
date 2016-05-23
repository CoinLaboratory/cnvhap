package lc1.util;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.Adler32;
import java.util.zip.CheckedInputStream;
import java.util.zip.CheckedOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;


/**
 * @author ebellos
 *
 */
public class Compressor {
    
    
    static int idColumn = 0;
    static int[] snpColumns = new int[] {3,2,1};
    static int[] resultColumns = new int[] {4,5};
    static String[] resultColumnHeader = new String[] {"","","","", "Log R","B allele"}; 
    static String split = "\t";
    final File in;
    final File out; 
    static boolean header =true;
    
    
    public   static void main(String[] args){
       try{
        File user = new File(System.getProperty("user.dir"));
        File[] f = user.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
               return pathname.getName().endsWith("data1M.txt");
            }
            
        });
        for(int i=0; i<f.length; i++){
            Compressor compr = new Compressor(f[i]);
            compr.compress();
          //  if(true) System.exit(0);
        }
       }catch(Exception exc){
           exc.printStackTrace();
       }
    }
   
   private void writeLine(OutputStreamWriter osw, String[] str) throws Exception{
       osw.write(str[resultColumns[0]]);
       for(int i=1; i<resultColumns.length; i++){
           osw.write("\t");
           int j = resultColumns[i];
           osw.write(str[j]);
       }
       osw.write("\n");
   }
   
   private void writeLine(OutputStreamWriter osw, List<String> ids) throws Exception{
       for(int i=0; i<ids.size(); i++){
           osw.write(ids.get(i)+"\n");
       }
   }
 
    static final int BUFFER = 2048;
    
    private String getEntryName(String[] str){
        StringBuffer sb = new StringBuffer(str[snpColumns[0]]);
        for(int i=1; i<snpColumns.length; i++){
            sb.append("_");
            sb.append(modify(str[snpColumns[i]]));
        }
        return sb.toString();
    }
    
    
    
    public  void compress () {
        if(out.exists())return ;
       try {
          List<String> ids = new ArrayList<String>();
          BufferedReader origin = new BufferedReader(new FileReader(in));
          if(header) origin.readLine();
          FileOutputStream dest = new      FileOutputStream(out);
          CheckedOutputStream checksum = new   CheckedOutputStream(dest, new Adler32());
          ZipOutputStream out = new 
            ZipOutputStream(new 
              BufferedOutputStream(checksum));
          OutputStreamWriter osw = new OutputStreamWriter(out);
          out.setMethod(ZipOutputStream.DEFLATED);
          String st = "";
          String comp = null;
          
          int index =   snpColumns[0];
          ZipEntry headings = new ZipEntry("Name");
          out.putNextEntry(headings);
          writeLine(osw, resultColumnHeader);
          osw.flush();
         for(int ik=0; (st = origin.readLine())!=null; ik++){
              String[] str = st.split(split);
              if(comp==null || !comp.equals(str[index])){
                  osw.flush();
                  ZipEntry entry = new ZipEntry(getEntryName(str));
                  out.putNextEntry(entry);
                  comp = str[index];
                  ik=0;
              }
              if(ik>=ids.size()) ids.add(str[0]);
              else if(!ids.get(ik).equals(str[0])) throw new RuntimeException("!!");
              writeLine(osw, str);
          }
          osw.flush();
          ZipEntry samples = new ZipEntry("Samples");
          out.putNextEntry(samples);
          writeLine(osw,ids);
          osw.flush();
          out.close();
          System.out.println("checksum: "
            +checksum.getChecksum().getValue());
       } catch(Exception e) {
          e.printStackTrace();
       }
    }
    
    public static void readZip(final ZipFile f, String st, String[] res, int col){
        try{
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(st))));
          String str = "";
          for(int i=0; (str = br.readLine())!=null && i<res.length;i++){
              String[] stri = str.split("\t");
              res[i]=stri[col];
          }
          br.close();
        }catch(Exception exc){
            exc.printStackTrace();
            //return null;
        }
    }
    
    public static void readZip(final ZipFile f, String st, String[][] res, int[] col){
        try{
        	int len = res[0].length;
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(st))));
          String str = "";
          for(int i=0; (str = br.readLine())!=null && i<len;i++){
              String[] stri = str.split("\t");
              for(int k=0; k<col.length; k++){
            	  res[k][i]=stri[col[k]];
              }
          }
          br.close();
        }catch(Exception exc){
            exc.printStackTrace();
            //return null;
        }
    }
    public static boolean readZip(ZipFile f, String st, double[][] res,
			int[] col, int offset, boolean[] strtype, double threshNA, boolean[] include) {
    	try{
    		  for(int k=0; k<col.length; k++){
    			  Arrays.fill(res[k],Double.NaN);
    		  }
        	int len = res[0].length;
        	ZipEntry ent = f.getEntry(st);
        	int cntNa =0;
        	
        	if(ent==null){
        		return false;
        	}
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(ent)));
          String str = "";
         
          int i2=0;
          int i=0;
          for( i=0; (str = br.readLine())!=null && i2<len;i++){
        	  if(include[i]){
	              String[] stri = str.split("\t");
	              int len1 = stri.length;
	              for(int k=0; k<col.length; k++){
	            	  int i1 = col[k];
	            	  if(i1 >= len1) {
	            		  return false;
	            	  }
	            	  String strv = stri[i1];
	            	  if(strv.indexOf('N')>0){
	            		  res[k][i2+offset] = Double.NaN;
	            		  cntNa++;
	            	  }
	            	  else res[k][i2+offset]= strtype[k] ? (strv.charAt(0)=='B' ? 2 :strv.charAt(1)=='B' ? 1 : 0 ) : Double.parseDouble(strv);
	              }
	              i2++;
        	  }
          }
          br.close();
          if(cntNa>threshNA * i2) {
        	  return false;
          }
        }catch(Exception exc){
            exc.printStackTrace();
            //return null;
        }
		return true;
	}
    public static String reverse(String st){
    	int len1 = st.length()-1;
    	char[] ch = new char[len1+1];
    	for(int i=0; i<ch.length; i++){
    		ch[i] = st.charAt(len1-i);
    	}
    	return new String(ch);
    }
    
    /* res is individual index by genotype index */
    public static String readZip(ZipFile f, ZipEntry ent, String st, double[][] res,
			boolean[] include, int[] col, String comment) throws Exception {
    	String comm="";
    	try{
    		  for(int k=0; k<col.length; k++){
    			  Arrays.fill(res[k],Double.NaN);
    		  }
        	int len = res.length;
        	//ZipEntry ent = f.getEntry(st);
        	int cntNa =0;
        	if(ent==null) return null;
        	 comm = ent.getComment();
        	boolean swtch = false;
        	if(comm!=null && comment.length()>0){
        		char[] ch = comm.toCharArray();
        		if(!comm.equals(comment)){
        			swtch = true;
        			if(!reverse(comment).equals(comm)) throw new RuntimeException("!! alleles don't match "+comm+" "+new String(comment));
        		}
        	}
        	
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(ent)));
          String str = "";
         
        
          for(int i=0; (str = br.readLine())!=null;i++){
        	  
        	String stri1 =   str.startsWith("./") ? str : str.split("\t")[col[0]];
        	  if(!stri1.startsWith(".")){
	              String[] stri = stri1.split(",");
	              int len1 = stri.length;
	              for(int k=0; k<stri.length; k++){
//	            	  int i1 = col[k];
	//            	  if(i1 >= len1) {
	  //          		  return false;
	    //        	  }
	            	  String strv = stri[k];
	            	res[i][k]=  Math.exp(Double.parseDouble(strv));
	              }
        	  }else{
        		  Arrays.fill(res[i], Double.NaN);
        		  include[i] = false;
        	  }
          }
          br.close();
          if(swtch){
      		System.err.println("switching!!!!!! "+st);
      		for(int i=0; i<res.length; i++){
      			double tmp = res[i][0];
      			res[i][0] = res[i][2];
      			res[i][2] = tmp;
      		}
      	}
        /*  if(cntNa>threshNA * i2) {
        	  return false;
          }*/
    	
        }catch(RuntimeException exc){
        	System.err.println(st+" "+exc.getMessage());
        	return null;
        }
    	catch(Exception exc){
            exc.printStackTrace();
            return null;
        }
		return comm;
	}
    
    
   public static BufferedReader readZip(final ZipFile f, final List<String> l ){
       Enumeration<InputStream> inputStreams = new Enumeration<InputStream>(){
          int i=0;
          public boolean hasMoreElements() {
              return i < l.size();
          
          }

          public InputStream nextElement() {
              try{
                  System.err.println(l.get(i));
                InputStream nxt = f.getInputStream(f.getEntry(l.get(i)));
              i++;
                return nxt;
              }catch(Exception exc){
                  exc.printStackTrace();
                  return null;
              }
          }
          
      };
     
     return new BufferedReader(new InputStreamReader(new SequenceInputStream(inputStreams)));
   }
   
   public static List<String> getIndiv(ZipFile f, String entryName, Integer column) throws Exception{
	   ZipEntry ent = f.getEntry(entryName);
	   if(ent==null){
		   ent = f.getEntry(entryName.replace('$', ':'));
		   if(ent==null) {
			   Logger.global.info("no entry "+entryName);
			   return null;
		   }
	   }
       BufferedReader  nxt ;
       if(entryName.endsWith(".gz")){
    	   nxt = new BufferedReader(new InputStreamReader(
    	          new GZIPInputStream( f.getInputStream(ent))));
       }
/*       else if (entryName.endsWith(".bin"))
       {
    	   return getIndiv(f, ent, column);    	   
       }*/
       else{
    	   nxt = new BufferedReader(new InputStreamReader(
           f.getInputStream(ent)));
       }
       return getIndiv(nxt, column);
   }
/*   
	static byte toByte(int number)
	{
		return (byte) (number-128);
	}
	
	static int toInt(byte number)
	{
		return ((int) number)+128;
	}   
   
   public static List<String> getIndiv(ZipFile f, ZipEntry ent, Integer column) {
       List<String> indiv = new ArrayList<String>();
       
       
	   	BufferedInputStream bis = new BufferedInputStream(f.getInputStream(ent));
		int size;
		byte[] buffer = new byte[2048];
		BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream((new File(path)).getParent()+File.separator+(new File(path)).getName()+".unzipped"), buffer.length);
		int lineCount=0;
		String line="";
		while ((size = bis.read(buffer, 0, buffer.length)) != -1) {
			for (int currByte=0;currByte<size;currByte+=2)
			{
				if ((toInt(buffer[currByte])==0) && (toInt(buffer[currByte+1])==0))
				{
					//System.out.println(lineCount);
					lineCount=0;
					indiv.add(arg0)
					line="";
					bos.write("\n".getBytes());
				}
				else
				{
					line+=Integer.toString(toInt(buffer[currByte]));
					lineCount+=toInt(buffer[currByte+1]);
					bos.write((Integer.toString(toInt(buffer[currByte]))+";"+Integer.toString(toInt(buffer[currByte+1]))+",").getBytes());
				}
			}
		}
       
       
       
       String st = "";
       while((st = nxt.readLine())!=null){
           String str = st.trim();
           if(column==null){
               indiv.add(str);
           }
           else{
        	   String[] str1 = str.split(Constants.splString());
               indiv.add(str1[column]);
           }
       }
      nxt.close();
       return indiv;
}
*/
public static List<String> getIndiv(BufferedReader nxt,  Integer column) throws Exception{
	 
      
       List<String> indiv = new ArrayList<String>();
       String st = "";
       for(int i=0; (st = nxt.readLine())!=null; i++){
           String str = st.trim();
           if(column==null){
               indiv.add(str);
           }
           else{
        	   String[] str1 = str.split(Constants.splString());
               indiv.add(str1[column]);
             //  if(i==44){
            //	   System.err.println("h");
              // }
           }
       }
      nxt.close();
       return indiv;
   }
   
   public static List<String> getIndiv(BufferedReader nxt,  Integer column,String spl) throws Exception{
		 
	      
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           String str = st.trim();
           if(column==null){
               indiv.add(str);
           }
           else{
        	   String[] str1 = str.split(spl);
               indiv.add(str1[column]);
           }
       }
      nxt.close();
       return indiv;
   }
   
   public static List<String> getIndiv(File f,  Integer column) throws Exception{
	  
       BufferedReader  nxt = 
           new BufferedReader(new FileReader(f));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           String str = st.trim();
           if(column==null){
               indiv.add(str);
           }
           else{
        	   String[] str1 = str.split(Constants.splString());
               indiv.add(str1[column]);
           }
       }
      nxt.close();
       return indiv;
   }
   public static List<String> getIndiv(ZipFile f, String entryName, Integer column, String spl) throws Exception{
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           String str = st.trim();
           if(column==null){
               indiv.add(str);
           }
           else{
        	   String[] st_ = str.split(spl);
        	   if(st_.length>column){
        		   String res = st_[column];
        		   if(res.length()>0)
        			   indiv.add(res);
        	   }
           }
       }
      nxt.close();
       return indiv;
   }
   public static boolean getIndiv(ZipFile f, String entryName, Integer column, String[] indiv) throws Exception{
	   ZipEntry ent = f.getEntry(entryName);
	   if(ent==null) {
		   return false;
	   }
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(ent)));
     //  List<String> indiv = new ArrayList<String>();
       String st = "";
       int i=0;
      for(; (st = nxt.readLine())!=null; i++){
           String str = st.trim();
           if(column==null){
               indiv[i] = str;
           }
           else{
               indiv[i] = (str.split("\t")[column]);
           }
       }
      
      nxt.close();
      return i==indiv.length;
    //   return indiv;
   }
   public static BufferedReader getBufferedReader(ZipFile f, String string) throws Exception{
	  ZipEntry entry =  f.getEntry(string);
	  if(entry==null) return null;
       return   new BufferedReader(new InputStreamReader(
               f.getInputStream(entry)));
         
   }
   private static void read(ZipFile f, String string, List<String> indiv,
           int i) throws Exception {
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(string))));
       String st = "";
       while((st = nxt.readLine())!=null){
           indiv.add(st.split("\t")[i]);
       }
       nxt.close();
       
   }
   private static String conc(String[] str, int from){
       StringBuffer sb = new StringBuffer(str[from]);
       for(int i=from+1; i<str.length; i++){
           sb.append(str[i]);
       }
       return sb.toString();
   }
   
  public static void readBuildFile(BufferedReader br, String chrom, final int from , final int to
          ,List<Integer> loc, List<String> chr,  List<String> snpid, 
          List<Character> majorAllele, List<Character> minorAllele,
          int loc_index, int chr_index, int snp_index, int[] maf_index) throws Exception{
      List<String> l = new ArrayList<String>();
   
      String st = "";
      while((st = br.readLine())!=null){
          String[] str = st.split("\\s+");
          if(str[chr_index].equals(chrom)){
              int no = Integer.parseInt(str[loc_index]);
              if(no>=from){
                  if(no <= to){
                     l.add(str[loc_index]);
                      loc.add(no);
                      chr.add(str[chr_index].substring(3));
                      snpid.add(str[snp_index]);
                      if(maf_index!=null &&maf_index[0]>=0 &&  str.length>maf_index[0]){
                          majorAllele.add(str[maf_index[0]].charAt(0));
                          minorAllele.add(str[maf_index[1]].charAt(0));
                      }
                   //   System.err.println(no);
                  }
                  else break;
              }
          }
      }
     // return readZip(zf, l);
      
  }
    public static BufferedReader readZip(final ZipFile f, final int from, final int to
            ,List<Integer> loc, List<String> chr,  List<String> snpid) throws Exception{
       List<String> l = new ArrayList<String>();
      Enumeration entries = f.entries();
           while(entries.hasMoreElements()){
               ZipEntry entry = (ZipEntry) entries.nextElement();
               String name = entry.getName();
               if(name.equals("Name") || name.equals("Sample")) continue;
               try{
                   String[] names = name.split("_");
               int no = Integer.parseInt(names[0]);
            //   System.err.println(no);
              if(no>=from){
                  if(no <= to){
                      l.add(name);
                      loc.add(no);
                      chr.add(names[1]);
                      snpid.add(conc(names,2));
                   //   System.err.println(no);
                  }
                  else break;
              }
               }catch(NumberFormatException exc){
                   continue;
               }
           }
          // read(f, l.get(0), indiv, 0);
           return readZip(f, l);
       }
    
    
    public static List<String> readZipFrom(ZipFile f, String name) throws Exception{
        BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(name))));
        String st;;// = br.readLine();
        List<String> l = new ArrayList<String>();
        while((st = br.readLine())!=null){
            l.add(st);
        }
        return l;
    }
    public static String join(String[] vals, String join){
    	StringBuffer r = new StringBuffer();
    	for(int i=0; i<vals.length; i++){
    		r.append(vals[i]+(i<vals.length-1? join : "")); 
    	}
    	return r.toString();
    }
    public static List<String> readZipFrom(ZipFile zf, String name,	boolean[] include,String join) throws Exception {
    	 BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(name))));
         String st;;// = br.readLine();
         List<String> l = new ArrayList<String>();
         int toinc = 0;
         for(int i=0; i<include.length; i++) if(include[i]) toinc++;
         String[] v1 = new String[toinc];
         while((st = br.readLine())!=null){
        	 String[] vals = st.split("\t");
        	 int k=0;
        	 for(int i=0; i<include.length; i++){
        		 if(include[i])v1[k] = vals[i];
        	 }
             l.add(join(v1, join));
         }
         br.close();
         return l;
	}
    
    public static List<String> readZipFrom(ZipFile f, String name, boolean[] include) throws Exception{
        BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(name))));
        String st;;// = br.readLine();
        List<String> l = new ArrayList<String>();
        for(int j=0;(st = br.readLine())!=null; j++){
        	if(include[j]){
        		l.add(st);
        	}
        }
        return l;
    }
    public static BufferedReader readZipFrom(ZipFile f, final int from, final int to
            ,List<Integer> loc, List<String> chr,  List<String> snpid
    ) throws Exception{
     //  final  ZipFile f = new ZipFile(in);
        final Enumeration entries =  f.entries();
        List<String> l = new ArrayList<String>();
        for(int i=0; i<from; i++){
            entries.nextElement();
        }
        for(int i=0; i<to; i++){
            String name =((ZipEntry)entries.nextElement()).getName(); 
            l.add(name);
            String[] names = name.split("_");
            loc.get(Integer.parseInt(names[0]));
            chr.add(names[1]);
            snpid.add(conc(names,2));
        }
     //   read(f, l.get(0), indiv, 0);
        return readZip(f, l);
        
        
     }
   

    public static void unzip(File in) {
        try {
           final int BUFFER = 2048;
           BufferedOutputStream dest = null;
           FileInputStream fis = new 
         FileInputStream(in);
           CheckedInputStream checksum = new 
             CheckedInputStream(fis, new Adler32());
           ZipInputStream zis = new 
             ZipInputStream(new 
               BufferedInputStream(checksum));
           ZipEntry entry;
           while((entry = zis.getNextEntry()) != null) {
              System.out.println("Extracting: " +entry);
              int count;
              byte data[] = new byte[BUFFER];
              // write the files to the disk
              FileOutputStream fos = new 
                FileOutputStream(entry.getName());
              dest = new BufferedOutputStream(fos, 
                BUFFER);
              while ((count = zis.read(data, 0, 
                BUFFER)) != -1) {
                 dest.write(data, 0, count);
              }
              dest.flush();
              dest.close();
           }
           zis.close();
           System.out.println("Checksum:  "+checksum.getChecksum().getValue());
            
        } catch(Exception e) {
           e.printStackTrace();
        }
     }
    
    public String modify(String st){
        return st.replace(':',',');
    }
      /*  if(i==0){
          
        }
        else return st;
        }*/
    
    
   
    Compressor(File f) throws Exception{
        this.in = f;
         out = new File(f.getAbsolutePath()+".zip");
      
        //unzip( out);
       // if(true) System.exit(0);
    }

    public static List<Integer> readPositions(File file) throws Exception {
      ZipFile f = new ZipFile(file);
      List<Integer> l = new ArrayList<Integer>();
      for(Enumeration en = f.entries(); en.hasMoreElements();){
          ZipEntry ent = (ZipEntry) en.nextElement();
          if(ent.getName().equals("Samples") || ent.getName().equals("Name")) continue;
          String nm = ent.getName().split("_")[0];
          l.add(Integer.parseInt(nm));
      }
      return l;
    }

	public static BufferedReader getBufferedReader(ZipFile zf, File file, String string) throws Exception {
		File f1 = new File(file, string);
		if(zf==null || f1.exists()){
			return new BufferedReader(new FileReader(f1));
		}
		else{
			return Compressor.getBufferedReader(zf, file.getName()+"/"+string);
		}
	}

		public static FileOutputStream getOS(File f) throws Exception{
		  //  if(f.exists() && f.length()>0) throw new RuntimeException("!!");
		    return new FileOutputStream(f);
		}

		public static List<String> getEntries(ZipFile originalZip, String string)throws Exception {
			List<String> l = new ArrayList<String>();
			BufferedReader br = Compressor.getBufferedReader(originalZip, string);
			String st = "";
			while((st = br.readLine())!=null){
				l.add(st);
			}
			br.close();
			return l;
		}

		
	}
