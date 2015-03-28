package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class ConvertHapMapPhased {

    public static void main(String[] args){
        try{
            File dir1 = new File(System.getProperties().getProperty("user.dir"));
            File[] files = dir1.listFiles(new FileFilter(){

                public boolean accept(File arg0) {
                  return arg0.getName().indexOf("phased")>=0;
                }
                
            });
            for(int i=0; i<files.length; i++){
                ConvertHapMapPhased chmp = new ConvertHapMapPhased(files[i]);
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
            "id"};
    ConvertHapMapPhased(File f) throws Exception{
        this.dir = f.getParentFile();
        prefix = f.getName().substring(0, f.getName().indexOf("phased"));
        int indofChr = prefix.indexOf("chr");
        chr = prefix.substring(indofChr+3).substring(prefix.indexOf('_', indofChr));
        dest = Compressor.getOS(new File(dir, chr+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
    }
    boolean trio = true;
    
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
    
    static boolean single =true;
    public void writeSamples() throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
        outS.putNextEntry(headings);
            for(int i=0; i<ids.size(); i++){
                if(!single){
                    osw.write(ids.get(i)+"\n");
                }
                else{
                    osw.write(ids.get(i)+"_1");
                    osw.write("\n");
                    osw.write(ids.get(i)+"_2");
                    osw.write("\n");
                }
            }
            osw.flush();
           outS.closeEntry();
    }
    List<Integer> loc1 = new ArrayList<Integer>();
    List<String> snps = new ArrayList<String>();
    List<String> major = new ArrayList<String>();
    List<String> minor =  new ArrayList<String>();
   List<String> list1 = new ArrayList<String>();
   List<String> list2 = new ArrayList<String>();
   List<String> ids = new ArrayList<String>();
   
    
    public  void run() throws Exception{
       writeHeader();
        File[] f = new File[] {
           new File(dir, prefix+"legend.txt"),
           new File(dir, prefix+"phased"),  
           new File(dir, prefix+"sample.txt"),  
        };
       
        
        Utils.readPosInfo(f[0], new int[] {0,1, 2, 3}, true, new List[] {snps, loc1, major, minor}, new Class[] {String.class, Integer.class, String.class, String.class});
        Utils.readPosInfo(new File(dir, prefix+"sample.txt"), new int[] {0}, false, new List[] {ids}, new Class[] {String.class});
        if(trio){
            int noFam = ids.size();
            int noParents = (int) Math.round((double)noFam*(2.0/3.0));
            if(noParents*(3.0/2.0)!=noFam) throw new RuntimeException("!!");
            ids = ids.subList(0, noParents);
        }
       
        
        
        writeSnps();
       writeSamples();
        BufferedReader br =Utils.getBufferedReader(f[1]);
        String st = ""; String st1 = "";
        
        for(int i=0; (st = br.readLine())!=null; i++){
            list1.add( st.replaceAll("\\s+", "").replace('0','A').replace('1','B'));
            list2.add( br.readLine().replaceAll("\\s+", "").replace('0','A').replace('1','B'));
        }
        for(int i=0; i<snps.size(); i++){
            ZipEntry headings = new ZipEntry(snps.get(i));
            outS.putNextEntry(headings);
            for(int j=0; j<this.ids.size(); j++){
                osw.write(list1.get(j).charAt(i));
                if(single) osw.write("\n");
                osw.write(list2.get(j).charAt(i)+"\n");
            }
            osw.flush();
            outS.closeEntry();
        }
      //  osw.close();
        outS.close();
    }
}
