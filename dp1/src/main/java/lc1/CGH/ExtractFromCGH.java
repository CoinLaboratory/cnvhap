package lc1.CGH;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ExtractFromCGH {
   
    List<String> chrs = new ArrayList<String>();
    {
        for(int i=1; i<=22; i++){
            chrs.add(i+"");
        
        }
        chrs.add("X");
        chrs.add("Y");
        
    }
    public static void main(String[] args){
       try{
       ExtractFromCGH eh = new ExtractFromCGH();
       eh.run();
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
   
    File[] files;
    Map<String, PrintWriter> output  = new HashMap<String, PrintWriter>();
    int[] indices = new int[] {11,12,16,17};
    Map<String, String> arrays;
    File outputParent;
    File[] outputDirs;
    BufferedReader[] in;
    PrintWriter log;
    ExtractFromCGH() throws Exception{
        File user = new File(System.getProperty("user.dir"));
        outputParent = user.getParentFile().getParentFile();
        outputDirs =  new File[chrs.size()];
        log= new PrintWriter(new BufferedWriter(new FileWriter("log.txt")));
        for(int i=0; i<outputDirs.length; i++){
            outputDirs[i] = new File(user, chrs.get(i));
            if(!outputDirs[i].exists()) outputDirs[i].mkdir();
        }
        arrays = readArrays();
        files = user.listFiles(new FileFilter(){
            public boolean accept(File f){
                return get(f.getName())!=null && f.getName().endsWith(".txt");
//                return f.getName().startsWith("US") && f.getName().endsWith("txt");
            }
        });
        in = new BufferedReader[files.length];
      if(in.length< arrays.size()){
            Set<String> set = new HashSet<String>();
            for(int i=0; i<in.length; i++){
                int index = files[i].getName().indexOf("_S0");
                System.err.println(files[i].getName());
                set.add(files[i].getName().substring(0, index));
            }
                if(set.size()!=in.length) {
                    throw new RuntimeException("!!");
                }
                else{
                    Set<String> set1 = arrays.keySet();
                    set1.removeAll(set);
                    for(Iterator<String> it = set1.iterator(); it.hasNext();){
                        log.println(it.next());
                    }
                    log.close();
                    System.err.println(set);
                }
                System.exit(0);
            throw new RuntimeException(in.length+" "+arrays.size());
        }
        StringBuffer sb = new StringBuffer();
        StringBuffer sb_h = new StringBuffer();
        String[] headRow = new String[files.length*2];
       
        for(int i=0; i<in.length; i++){
            in[i] = new BufferedReader(new FileReader(files[i]));
            sb.append("\t%10.3f\t%10.3f");
            sb_h.append("\t%10s\t%10s");
            String name = files[i].getName();
            String alt = get(name);
            if(alt==null){
                System.err.println("not contained "+name);
            }
            else{
                System.err.println("cont "+name);
            }
            headRow[2*i] = alt+"_R";
            headRow[2*i+1]= alt+"_RU";
        }
        formatString = sb.toString();
        formatString_h = sb_h.toString();
        for(int i=0; i<outputDirs.length; i++){
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputDirs[i], "cghdata.txt"))));
            output.put(outputDirs[i].getName(),pw);
            pw.println("Name\tStart\tEnd"+String.format(formatString_h, headRow));
        }
       
       
    }
    public void finish() throws Exception{
        for(int i=0; i<in.length; i++){
            in[i].close();
        }
        for(Iterator<PrintWriter> it = output.values().iterator(); it.hasNext();){
            it.next().close();
        }
    }
    public String get(String st){
       // return arrays.get(st);
       for(Iterator<String> it = arrays.keySet().iterator(); it.hasNext(); ){
           String key = it.next();
           if(st.startsWith(key)){
               return arrays.get(key);
           }
        }
        return null;
    }
    public Map<String, String> readArrays() throws Exception{
        Map<String, String> arrays = new HashMap<String, String>();
        File f = new File("aCGH_barcodes.csv");
        BufferedReader br = new BufferedReader(new FileReader(f));
        String st = br.readLine();
        while((st = br.readLine())!=null){
            String[] str = st.split("\\t");
         arrays.put(str[0]
                        //+"_S01_CGH-v4_91.txt"
                        , str[1]);   
        }
        return arrays;
    }
    String formatString;
    String formatString_h;
    private void checkHeaders(String st){
        if(st.indexOf("LogRatio")>=0 && st.indexOf("ProbeName")>=0){
            String[] str = st.split("\\t");
            if(!str[15].equals("LogRatio") || !str[16].equals("LogRatioError")) throw new RuntimeException("!!");
        }
    }
    
    String[] forb = new String[] {"Corner", "NC2", "SRN", "SM", "PC"};
    public void run() throws Exception{
//        String[] str = new String[in.length];
        Double[] res = new Double[in.length*2];
        String[] str;
        String st="";
        String name="";
        String chr="";
        int start=0;
        int end=0;
        for(int ik=0; ik<13; ik++){
            for(int i=0; i<in.length;i++){
                st = in[i].readLine();
                if(i==0)System.err.println("skipping "+st);
               checkHeaders(st);
            }
        }
        System.err.println(st);
        outer: while(true){
            
            for(int i=0; i<in.length;i++){
                 st = in[i].readLine();
                 if(st==null) break outer;
                 str = st.split("\\t");
                 try{
                 
                 if(i==0){
                     for(int j=0; j<forb.length; j++){
                     if(st.indexOf("A_")<0){
                         checkHeaders(st);
                         for(int i1=1; i1<in.length; i1++){
                             st = in[i1].readLine();
                             if(st.indexOf("A_")>=0) throw new RuntimeException("");
                         }
                         continue outer;
                     }
                     
                     }
                     name = str[10];
                     String[] st1 = str[12].split(":");
                     chr = st1[0].substring(st1[0].indexOf("chr")+3);
                //     System.err.println(Arrays.asList(st1));
                     String[] st2 = st1[1].split("-");
                     start = Integer.parseInt(st2[0]);
                     end = Integer.parseInt(st2[1]);
                     
                 }
                 else if(!str[10].equals(name)) {
                     throw new RuntimeException("!!");
                 }
                 }catch(Exception exc){
                     exc.printStackTrace();
                     System.err.println("prob with "+str[12]+"\n"+st);
                     System.exit(0);
                 }
                 String sc1 = str[15];
                 String sc2 = str[16];
                 res[2*i] = Double.parseDouble(sc1);
                 res[2*i+1] = Double.parseDouble(sc2);
            }
            PrintWriter outp = this.output.get(chr);
            outp.println(name+"\t"+start+"\t"+end+String.format(formatString, res));
            outp.flush();
        }
        this.finish();
    }
    
    
}
