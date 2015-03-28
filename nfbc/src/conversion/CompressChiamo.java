package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class CompressChiamo {
    public static void main(String[] args){
        try{
            File dir1 = new File(System.getProperties().getProperty("user.dir"));
            File[] files = dir1.listFiles(new FileFilter(){

                public boolean accept(File arg0) {
                  return arg0.getName().indexOf(".txt")>=0;
                }
                
            });
            for(int i=0; i<files.length; i++){
                CompressChiamo chmp = new CompressChiamo(files[i]);
                System.err.println("running "+i);
                chmp.run();
            }
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
    
  String[][] header =   new String[][] {
            new String[] {"Genotype"},//,tScore",
            new String[] {"id"},
            new String[]{"id"}};
  File f;
  
  final List<String> file_header;
  final int id_index, sample_id;
  final int[] geno_id;
    CompressChiamo(File f) throws Exception{
        this.dir = f.getParentFile();
        this.f =f;
       br  =Utils.getBufferedReader(f);
    
      file_header = Arrays.asList( br.readLine().split("\t"));
      id_index = file_header.indexOf("MARKER");
      sample_id = file_header.indexOf("PATIENT");
      geno_id = new int[] {file_header.indexOf("ALLELES")};
    //   String[] header = br.readLine().split();
        prefix = f.getName().substring(0, f.getName().indexOf(".txt"));
       // chr = prefix.split("_")[1].substring(3);
        dest = Compressor.getOS(new File(dir, prefix+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        System.err.println("new output stream");
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
    }
    
    
    
    public void writeHeader() throws Exception{
        ZipEntry headings = new ZipEntry("Name");
        outS.putNextEntry(headings);
            for(int i=0; i<header.length; i++){
            	for(int j=0; j<header[i].length; j++){
            		osw.write(header[i][j]);
            		if(j<header[i].length-1) osw.write("\t");
            	}
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
    
    public void writeSamples() throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
        outS.putNextEntry(headings);
            for(int i=0; i<ids.size(); i++){
                    osw.write(ids.get(i)+"\n");
            }
            osw.flush();
           outS.closeEntry();
    }
   // List<Integer> loc1 = new ArrayList<Integer>();
    List<String> snps = new ArrayList<String>();
  //  List<String> major = new ArrayList<String>();
 //   List<String> minor =  new ArrayList<String>();
  // List<String> list1 = new ArrayList<String>();
   //List<String> list2 = new ArrayList<String>();
   List<String> ids = new ArrayList<String>();
   
   final BufferedReader br;
    public  void run() throws Exception{
       writeHeader();
     
       String st = "";
       String snp_id ="";
       boolean first = true;
       int sample_cnt =0;
       String[] genos=null;
       while((st = br.readLine())!=null){
           String[] str = st.split("\t");
           if(!snp_id.equals(str[this.id_index])){
               snp_id = str[id_index];
             
               this.snps.add(snp_id);
              
               if(!first){
            	   if(genos==null){
                	   genos =  new String[sample_cnt];
                   }
                   else{
                	   for(int k=0; k<genos.length; k++){
                		   osw.write(genos[k]);
                	   }
                	   Arrays.fill(genos, null);
                   }
            	   
                   osw.flush();
                   outS.closeEntry();
               }
               else{
                   first = false;
               }
               sample_cnt =0;
               ZipEntry headings = new ZipEntry(snp_id);
               outS.putNextEntry(headings);
           }
           if(sample_cnt>=ids.size()){
               ids.add(str[sample_id]);
               for(int k=0; k<geno_id.length; k++){
            	   if(str.length <= geno_id[k]) osw.write("NC");
            	   else osw.write(str[geno_id[k]]);
            	   if(k<geno_id.length-1) osw.write("\t");
            	   else osw.write("\n");
               }
           }
           else{
        	   int pos_index = ids.indexOf(str[sample_id]);
        	   StringBuffer sb = new StringBuffer();
        	   for(int k=0; k<geno_id.length; k++){
            	   if(str.length <= geno_id[k]) sb.append("NC");
            	   else sb.append(str[geno_id[k]]);
            	   if(k<geno_id.length-1) sb.append("\t");
            	   else sb.append("\n");
               }
        	   genos[pos_index] = sb.toString();
           }
           sample_cnt++;
          
       }
       if(genos!=null){
	       for(int k=0; k<genos.length; k++){
			   osw.write(genos[k]);
		   }
       }
       osw.flush();
       outS.closeEntry();
        writeSnps();
       writeSamples();
       System.err.println("closing");
        outS.close();
    }
}
