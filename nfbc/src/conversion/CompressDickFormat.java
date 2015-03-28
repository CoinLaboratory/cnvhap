package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class CompressDickFormat {

    public static void main(String[] args){
        try{
          
                CompressDickFormat chmp = new CompressDickFormat(new File(args[0]));
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
    CompressDickFormat(File f) throws Exception{
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
    public void writeSamples(String[] ids) throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
        outS.putNextEntry(headings);
            for(int i=2; i<ids.length; i++){
                    osw.write(ids[i]+"\n");
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
    
    public  void run() throws Exception{
       writeHeader();
       BufferedReader br =Utils.getBufferedReader(f);
       writeSamples( br.readLine().trim().split("\\s+"));
       
       String st = "";
      
       while((st = br.readLine())!=null){
          snps.add( writeSNP( st.trim().split("\\s+")));
          
       }
        writeSnps();
       
        outS.close();
    }
}
