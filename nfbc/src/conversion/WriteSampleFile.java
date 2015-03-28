package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

public class WriteSampleFile {
    
    
    public WriteSampleFile(File file) {
        this.dir = file;
    }

    public static void main(String[] args){
        try{
            File user = new File(System.getProperty("user.dir"));
            File[] f = user.listFiles(new FileFilter(){

                public boolean accept(File pathname) {
                   return pathname.getName().startsWith("WG"); 
                 //  &&   !pathname.getName().endsWith("zip")
                }
                
            });
            sampleInfo = readSampleInfo(new File(user, "list300k_31Aug2007.txt"));
            for(int i=0; i<f.length; i++){
                WriteSampleFile wsw = new WriteSampleFile(f[i]);
                wsw.run();
            }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
static Map<String, String[]> sampleInfo;

static Map<String, String[]> readSampleInfo(File f) throws Exception{
          
        BufferedReader br = new BufferedReader(new FileReader(f));
        HashMap<String, String[]> cc = new HashMap<String, String[]>();
        if(br!=null){
            String st = br.readLine();
            while((st = br.readLine())!=null){
                String[] str = st.split("\\s+");
                String id = str[1];
               
                 cc.put(id, str);
            }
            br.close();
            //if(cc.size()!=this.dataL.size()) throw new RuntimeException("!!!");
            
        }
        return cc;

}
     
    File dir;
    
   
    
    public void run() throws Exception {
        ZipFile zf = new ZipFile(new File(dir, "1.zip"));
        List<String> res = Compressor.getIndiv(zf, "Samples");
        PrintWriter sam = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Samples.txt"))));
        sam.println("id\tT2Diabetes\tPopulation\tSub");
        boolean dummyP = Math.random() < 0.5;
        for(int i=0; i<res.size(); i++){
            String val = res.get(i).split("\\s+")[0];
            String[] vals = sampleInfo.get(val.split("#")[0]);
            String pop = vals[2];
           
            String cc = pop.equals("Desir") ? "Control" :"Case";
            System.err.println(dir.getName()+" "+cc);
            sam.println(val+"\t"+cc+"\t"+pop+"\t"+pop+(dummyP?"a":"b"));
        }
        sam.close();
    }
    
}
