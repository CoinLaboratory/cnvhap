package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ModifyGCFile {
    
    public static void main(String[] args){
        boolean reverse = true;;
        try{
            File dir = new File(System.getProperty("user.dir"));
            File outF  = new File("gc_build36.txt");
            File[] f = (dir).listFiles(new FilenameFilter(){

                public boolean accept(File arg0, String arg1) {
                    return arg1.indexOf("1k")>=0;
                }
                
            });
            if(reverse){
                File outdir = new File(dir, "build36");
                outF = new File(dir, "gc_build36.txt");
                for(int i=0; i<f.length; i++){
                    f[i] = new File(outdir, f[i].getName());
                }
            }
            ModifyGCFile mb = new ModifyGCFile(f,  outF);
           if(reverse) mb.reverse();
           else mb.run();
            
        }catch(Exception exc){
            exc.printStackTrace();
        }
        
        
        
    }
    
    
    private void reverse() throws Exception{
        BufferedReader in = new BufferedReader(new FileReader(outF));
        String st = "";
        String chr = "";
        PrintWriter pw = null;
        while((st = in.readLine())!=null){
            String[] str = st.split("\t");
            if(!str[0].equals(chr)){
                 if(pw!=null) pw.close();
                chr = str[0];
                File out = m.get(chr);
                if(out.exists()) throw new RuntimeException("!!");
                pw = new PrintWriter(new FileWriter(out));
                
            }
            pw.print(str[1]);
            for(int i=2; i<str.length; i++){
                pw.print("\t"+str[i]);
            }
            pw.println();
        }
        in.close();
    }
    private void run() throws Exception {
        out = new PrintWriter(new FileWriter(outF));
      for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
          String chr = it.next();
         
          readGC(chr);
       
      }
        
    }
    PrintWriter out;
    File outF;
    Map<String, File> m;
    ModifyGCFile(File []f, File outF) throws Exception{
       
        m = new HashMap<String, File>();
       for(int i=0; i<f.length; i++){
           m.put("chr"+f[i].getName().split("_")[0], f[i]);
                  
       }
       this.outF = outF;
     
     
    }
    
    private void readGC(String chr2) throws Exception {
        File f = m.get(chr2);       
      
        BufferedReader br1 = new BufferedReader(new FileReader(f));
        String stri = "";
        while((stri = br1.readLine())!=null){
//            String[] str = stri.split("\\s+");
            out.println(chr2+"\t"+stri);
        }
        br1.close();
    }
}
