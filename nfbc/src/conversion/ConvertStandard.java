package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class ConvertStandard {
    public static void main(String[] args){
        try{
            File dir1 = new File(System.getProperties().getProperty("user.dir"));
            File[] files = dir1.listFiles(new FileFilter(){

                public boolean accept(File arg0) {
                  return arg0.getName().indexOf(".txt")>=0;
                }
                
            });
            for(int i=0; i<files.length; i++){
                ConvertStandard chmp = new ConvertStandard(files[i]);
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
    String chr;
    
  String[] header =   new String[] {
            "Genotype",
            "chr\tbuild35_start\tbuild35_end\tid\tA\tB",
            "i"};
  File file;
    ConvertStandard(File f) throws Exception{
        this.dir = f.getParentFile();
        this.file = f;
        prefix = f.getName().substring(0, f.getName().indexOf(".txt"));
        chr = prefix.split("_")[1].substring(3);
        dest = Compressor.getOS(new File(dir, chr+".zip"));
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
        for(int i=0; i<loc1.size(); i++){
            int loc = loc1.get(i);
            osw.write("chr"+chr+"\t"+loc+"\t"+(loc+40)+"\t"+snps.get(i)+"\t"+major.get(i)+"\t"+minor.get(i));
            osw.write("\n");
        }
        osw.flush();
       outS.closeEntry();
}
public void writeSamples() throws Exception{
    ZipEntry headings = new ZipEntry("Samples");
    outS.putNextEntry(headings);
        for(int i=0; i<ids.size(); i++){
            osw.write(ids.get(i));
            osw.write("\n");
          //  osw.write(ids.get(i)+"_2");
           // osw.write("\n");
        }
        osw.flush();
       outS.closeEntry();
}
List<Integer> loc1 = new ArrayList<Integer>();
List<String> snps = new ArrayList<String>();
List<String> major = new ArrayList<String>();
List<String> minor =  new ArrayList<String>();
List<String> list = new ArrayList<String>();
List<String> ids = new ArrayList<String>();

public  void run() throws Exception{
    writeHeader();
    BufferedReader br = new BufferedReader(new FileReader(file));
    String[] header = br.readLine().trim().split("\\s+");
    for(int i=4; i<header.length; i++){
        ids.add(header[i]);
    }
    String st;
    while((st = br.readLine())!=null){
        String[] str = st.trim().split("\\s+");
        this.snps.add(str[0]);
        this.loc1.add(Integer.parseInt(str[1]));
        this.major.add(str[2]);
        this.minor.add(str[3]);
        ZipEntry headings = new ZipEntry(str[0]);
        outS.putNextEntry(headings);
        for(int i=4; i<header.length; i++){
            osw.write(str[i]+"\n");
        }
        osw.flush();
        outS.closeEntry();
       
    }
     
     
     writeSnps();
    writeSamples();
     outS.close();
 }


}
