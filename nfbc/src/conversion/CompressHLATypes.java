package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class CompressHLATypes {

    public static void main(String[] args){
        try{
          
                CompressHLATypes chmp = new CompressHLATypes(new File(args[0]));
                chmp.run();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
 
    File dir;
    String prefix;
    FileOutputStream dest;
    CheckedOutputStream checksum;
    ZipOutputStream outS;
    OutputStreamWriter osw;
    ///String chr;
    
  String[] header =   new String[] {
            "Genotype",
            "chr\tstart\tend\tid",
            "id"};
  File f;
  
  
  
    CompressHLATypes(File f) throws Exception{
        this.dir = f.getParentFile();
        this.f =f;
        prefix = f.getName().substring(0, f.getName().indexOf(".txt"));
       // chr = prefix.split("_")[1].substring(3);
        dest = Compressor.getOS(new File(dir, prefix+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
    }
    
    public void writeHeader() throws Exception{
        ZipEntry headings = new ZipEntry("Name");
        outS.putNextEntry(headings);
            for(int i=0; i<header.length; i++){
                osw.write(header[i]);
                osw.write("\n");
            }
            osw.flush();
           outS.closeEntry();
    }
    
    public void writeSnps() throws Exception{
        ZipEntry headings = new ZipEntry("SNPS");
        outS.putNextEntry(headings);
            for(int i=0; i<snps.size(); i++){
                osw.write(snps.get(i));
                osw.write("\n");
            }
            osw.flush();
           outS.closeEntry();
    }
    List<String> snps = new ArrayList<String>();
    public void writeSamples(List<String> ids) throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
        outS.putNextEntry(headings);
            for(int i=2; i<ids.size(); i++){
                    osw.write(ids.get(i)+"\n");
            }
            osw.flush();
           outS.closeEntry();
    }
 
  public String writeSNP(String[] str) throws Exception{
      String snp_id = str[1];
      ZipEntry headings = new ZipEntry(snp_id);
      outS.putNextEntry(headings);
      for(int i=2; i<str.length; i++){
          osw.write(str[i]+"\n");
  }
      osw.flush();
      outS.closeEntry();
      return "chr"+str[0]+"\t"+str[1]+"\t"+str[1]+"\t"+str[1];
  }
    
  BufferedReader br;
  
  public List<String> readColumn(int col) throws Exception{
	  List<String> l = new ArrayList<String>();
	  BufferedReader br =Utils.getBufferedReader(f);
	  String st = "";
	  while((st = br.readLine())!=null){
         l.add(st.trim().split("\\s+")[col]);
	  }
	  br.close();
      return l;
  }
  
    public  void run() throws Exception{
       writeHeader();
       
       BufferedReader br =Utils.getBufferedReader(f);
      int no_cols =  br.readLine().trim().split("\\s+").length;
      br.close();
       List<String> st = readColumn(0);
       int len = st.size();
       writeSamples(st);
       
      for(int i1 =1; i1<no_cols; i1+=2){
    	  List<String> st1 = readColumn(i1);
    	  List<String> st2 = readColumn(i1+1);
    	  String snpid = st1.get(0);
    	  Map<String, Character> m = new HashMap<String, Character>();
    	  int latest = 65;
    	  ZipEntry headings = new ZipEntry(snpid);
	      outS.putNextEntry(headings);
    	  for(int i=2; i<st1.size(); i++){
    		  String str1 = st1.get(i);
    		  String str2 = st2.get(i);
    		  Character c1 = m.get(str1);
    		  if(c1==null){
    			  System.err.println(str1+" "+latest+" "+m.keySet().size());
    			  m.put(str1,c1 = (char)latest);
    			  latest++;
    		  }
    		  Character c2 = m.get(str2);
    		  if(c2==null){
    			  System.err.println(str2+" "+latest+" "+m.keySet().size());
    			  m.put(str2,c2 = (char)latest);
    			  latest++;
    		  }
    	          osw.write(c1.toString()+c2.toString()+"\n");
    	    
    	  
    		  
    	  }
    	    osw.flush();
  	      outS.closeEntry();
    	  
    	 
    	  snps.add(prefix+"\t"+st1.get(1)+"\t"+st2.get(1)+"\t"+snpid+getString(m));
    	
    	  //writeSNP
      }
      
      
        writeSnps();
       
        outS.close();
    }

	private String getString(Map<String, Character> m) {
		StringBuffer sb = new StringBuffer();
		for(Iterator<Entry<String, Character>> it = m.entrySet().iterator(); it.hasNext();){
			Entry<String, Character> ent = it.next();
			sb.append("\t"+ent.getKey()+"="+ent.getValue());
		}
		return sb.toString();
	}
}
