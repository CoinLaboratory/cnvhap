package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

public class ModifyBuild {
   public static void main(String[] args){
        try{
            ModifyBuild mb = new ModifyBuild();
           
           mb.run();
            
        }catch(Exception exc){
            exc.printStackTrace();
        }
        
    }
    Map<String, File> m;
    String chr = "";
    PrintWriter out;
    BufferedReader br;
    
    ModifyBuild() throws Exception{
        br = new BufferedReader(new FileReader("build35.txt"));
        File[] f = (new File("../gc_data/")).listFiles();
        m = new HashMap<String, File>();
       for(int i=0; i<f.length; i++){
           m.put("chr"+f[i].getName().split("_")[0], f[i]);
                  
       }
      out = new PrintWriter(new FileWriter("build35_mod.txt"));
    }
    String st;
   SortedMap<Integer, Double> gc = new TreeMap<Integer, Double>();
    public void run() throws Exception{
        int i=0; 
        while((st = br.readLine())!=null){
            String[] str = st.split("\t");
            if(!str[0].equals(chr)){
                chr = str[0];
                readGC(chr);
            }
            int pos = Integer.parseInt(str[1]);
            double gc = getGC(pos);
            out.println(st+"\t"+Math.round(gc));
            //System.err.println(i);
            i++;
        }
        out.flush();
        out.close();
    }
    private double getGC(int pos) {
       SortedMap<Integer, Double> gt = gc.tailMap(pos);
       SortedMap<Integer, Double> lt = gc.headMap(pos+1);
       Integer fgt = gt.size()==0 ? null : gt.firstKey();
       Integer llt = lt.size()==0 ? null : lt.lastKey();
       if(fgt==null){
           return lt.get(llt);
       }
       else if(llt==null){
           return gt.get(fgt);
       }
       else if(Math.abs(pos - fgt) < Math.abs(pos -llt)){
           
           return gt.get(fgt);
           
       }
       else{
           return lt.get(llt);
       }
    }
    private void readGC(String chr2) throws Exception {
        File f = m.get(chr2);       
        this.gc.clear();
        BufferedReader br1 = new BufferedReader(new FileReader(f));
        String stri = "";
        while((stri = br1.readLine())!=null){
            String[] str = stri.split("\\s+");
            gc.put((int) Math.round((Integer.parseInt(str[0])+Integer.parseInt(str[1]))/2.0), Double.parseDouble(str[2]));
        }
        br1.close();
    }
}
