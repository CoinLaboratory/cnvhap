package conversion;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public class ModifyHeader {
    
    File in;
    ZipOutputStream outS = null;
    OutputStreamWriter osw = null;
    CheckedOutputStream checksum;
    
   ModifyHeader(File in){
       this.in = in;
   }
    
    public static void main(final String[] args){
        try{
           //if(true) throw new RuntimeException("!!");
            File user = new File(System.getProperty("user.dir"));
            File[] fs = user.listFiles(new FilenameFilter(){

                public boolean accept(File arg0, String arg1) {
                   return arg1.endsWith("zip");
                }
                
            });
            for(int i=0; i<fs.length; i++){
                String nme = fs[i].getName();
                nme = nme.substring(0,nme.indexOf(".zip"));
                ModifyHeader comp = new ModifyHeader(fs[i]);
                comp.copy(map1 , false);
//                        new String[] {
  //                      "Log R\tAgilentPred",
    //                    "chr\tbuild35_start\tbuild35_end\tid\tgc\tloess",
      //                  "id\tmedianCorrection"});
                 
                        
            }
            
            
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    static Map<String, String[]> map1 = new HashMap<String, String[]>();
    static{
      //  map1.put("Name", new String[] {"population"});
      //  map1.put("Samples", new String[] {"CEU"});
        map1.put( "Name", 
                new String[] {"Genotype",
                "chr\tbuild35_start\tbuild35_end\tid\tA\tB", 
                "id\tpopulation"});
    }
           
   
    
    
public  void copy(Map<String, String[]> rep, boolean addition) throws Exception{
    //List<String> excs = Arrays.asList(exception);
    File f1 = this.in;
    ZipFile zf = new ZipFile(f1);
   
    String chr = f1.getName().substring(0,f1.getName().indexOf(".zip"));
    File tmp = new File(f1.getParentFile(), chr+System.currentTimeMillis()+".zip");
    
    FileOutputStream dest = Compressor.getOS(tmp);
    checksum = new   CheckedOutputStream(dest, new Adler32());
    outS = new 
     ZipOutputStream(new 
       BufferedOutputStream(checksum));
   osw = new OutputStreamWriter(outS);
   outS.setMethod(ZipOutputStream.DEFLATED);
    for(Enumeration en = zf.entries(); en.hasMoreElements();){
        ZipEntry ent = (ZipEntry) en.nextElement();
       
        List<String> res = Compressor.getIndiv(zf, ent.getName());
           ZipEntry headings = new ZipEntry(ent.getName());
           outS.putNextEntry(headings);
           String[] replacement = rep.get(ent.getName());
           if(replacement!=null){
               if(addition){
                   for(int i=0; i<res.size(); i++){
                       osw.write(res.get(i));
                       osw.write("\t");
                       osw.write(replacement[0]);
                       osw.write("\n");
                   }
               }
               for(int i=0; i<replacement.length; i++){
                   osw.write(replacement[i]);
                   osw.write("\n");
               }
           }
           else{
               for(int i=0; i<res.size(); i++){
                   osw.write(res.get(i));
                   osw.write("\n");
               }
           }
           osw.flush();
           outS.closeEntry();
    }
    outS.close();
    zf.close();
   boolean deleted = f1.delete();
   if(deleted) tmp.renameTo(f1);
   else{
       throw new RuntimeException("could not move");
   }
}



    
}
