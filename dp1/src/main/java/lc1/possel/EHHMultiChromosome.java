package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import lc1.CGH.AberationFinder;
import lc1.stats.TrainableNormal;
import lc1.util.Constants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class EHHMultiChromosome {
    final File user;
    Map<Double, TrainableNormal> normal;
   
  //  File baseDir;
   public static boolean overl = true;
   public static boolean graphAll = false;
    int noCop;
   PrintWriter score;
    public void setNoCop(int  noCop) throws Exception{
        File nF = new File(user, "norm.txt");
        scoring = nF.exists() && nF.length()>0;
          normal = read(nF);
        File out = new File(user, "score_"+noCop+"_"+"0.txt");
        for (int i=0; out.exists() && out.length()>0; i++){
            out = new File(user, "score_"+noCop+"_"+i+".txt");
        }
        System.err.println("out file is "+out);
        score=  new PrintWriter(new BufferedWriter(new FileWriter(out)));
        String fs1 = "%-5s\t%15s\t%15s\t%10s\t%5s\t%5s\t%5s";
        Object[] toPrintH = new Object[] {"Name", "snp_id",  "snp_id", "loc", "derF", "ratio", "pval"};
       
        score.println(String.format(fs1, toPrintH));
        score.flush();
        this.noCop = noCop;
    }
    public EHHMultiChromosome(File user, File baseDir, File input, int num) throws Exception{
        this.user = user;
     
       
     
      //  this.baseDir = baseDir;
        File[] f;
        if(input==null){
            File[] f1 = user.listFiles(new FileFilter(){
    
                public boolean accept(File pathname) {
                  return pathname.getName().endsWith(".tar.gz") && !pathname.getName().startsWith("X") && !pathname.getName().startsWith("Y");
                }
                   
               });
            f  =  new File[f1.length];
            for(int i=0; i<f.length; i++){
                f[i] = f1[i];
            }
        }
        else{
            throw new RuntimeException("!!");
/*            List<File> l = new ArrayList<File>();
            BufferedReader br = new BufferedReader(new FileReader(input));
            String st = br.readLine();
            while((st = br.readLine())!=null){
                String[] str = st.split("\\s+");
                File f1 = new File(user, str[1]+".tar.gz");
                l.add(new File(f1,Integer.parseInt(str[3]),num ));
                if(!str[1].equals(str[2])) {
                    File f2 = new File(user, str[2]+".tar.gz");
                    int i1 = Integer.parseInt(str[3]);
                    File s = new File(f2, i1, num);
                    l.add(s);
                   
                }
            }
            f =l.toArray(new File[0]);*/
        }
       /* for(int i=0; i<f.length; i++){
            String chr = 
                split(f[i]).get(0);
            if(!overl) chr = chr+"_"+Constants.nextInt(1000);
                List<File> l = fileMap.get(chr);
                if(l==null) fileMap.put(chr, l =  new ArrayList<File>());
                l.add(f[i]);
            
        }*/
    }
   
    
     static Map<Double, TrainableNormal> read(File file) throws Exception {
        Map<Double, TrainableNormal> normal = new TreeMap<Double, TrainableNormal>();
        if(file.exists()){
           BufferedReader br = new BufferedReader(new FileReader(file));
           String st = "";
           while((st = br.readLine())!=null){
               String[] str = st.split("\\t");
               normal.put(Double.parseDouble(str[0]), new TrainableNormal(str[2], 100));
               //c(0.0,0.7949197906775308,0);
             
           }
           return normal;
        }
        return new HashMap<Double, TrainableNormal>();
     }
    static Comparator<File> cc = new Comparator<File>(){

        public int compare(File o1, File o2) {
          int i1 = Integer.parseInt(split(o1).get(1));
          int i2 = Integer.parseInt(split(o2).get(1));
          if(i1!=i2) return i1<i2 ? -1: 1;
          else return 0;
        }
        
    };
    static List<String> split(File f){
        String[] st = f.getName().split("\\.")[0].split("_");
        return Arrays.asList(st);//.subList(st.length-3,st.length);
    }
    boolean scoring;
    Class clazz = LDInner.class; //EHH.class
 /** object is grouped by chromosome;  [file, start_pos_of_next_file; end_pos_of_prev_file*/
    public static Iterator<Iterator<Object[]>> getIterator(File user){
        final Map<String, List<File>> fileMap = new TreeMap<String, List<File>>(new Comparator<String>(){

            public int compare(String o1, String o2) {
                String[] oo1 = o1.split("_");
                String[] oo2 = o2.split("_");
                for(int i=0; i<oo1.length; i++){
                    int i1 = Integer.parseInt(oo1[i]);
                    int i2 = Integer.parseInt(oo2[i]);
                    if(i1!=i2) return i1<i2 ? -1 :1;
                }
                return 0;
            }
            
        });
        File[] f = user.listFiles(new FileFilter(){
            
            public boolean accept(File pathname) {
              return pathname.getName().endsWith(".tar.gz") && !pathname.getName().startsWith("X") && !pathname.getName().startsWith("Y");
            }
               
           });
        
        for(int i=0; i<f.length; i++){
            String chr = 
                split(f[i]).get(0);
           // if(!overl) chr = chr+"_"+Constants.nextInt(1000);
                List<File> l = fileMap.get(chr);
                if(l==null) fileMap.put(chr, l =  new ArrayList<File>());
                l.add(f[i]);
            
        }
        return new Iterator<Iterator<Object[]>>(){
            Iterator<List<File>> it1 = fileMap.values().iterator();
            public boolean hasNext(){
                return it1.hasNext();
            }
            public void remove(){}
            public Iterator<Object[]> next(){
                final List<File> nxt = it1.next();
                Collections.sort(nxt, cc);
                return new Iterator<Object[]>(){
                    int i=0; 
                    int end_prev = 0;
                    Object[] res = new Object[3];
                    public boolean hasNext() {
                       return i<nxt.size();
                    }

                    public Object[] next() {
                        res[0] = nxt.get(i);
                        res[2] = end_prev;
                      
                        if(i<nxt.size()-1){
                            try{
                            File f_i =  nxt.get(i+1);
                            BufferedReader br = AberationFinder.getBufferedReader(f_i, "snp.txt");
                            br.readLine();
                            res[1] = Integer.parseInt(br.readLine().split("\\s+")[2]);
                            br.close();
                            }catch(Exception exc){
                                exc.printStackTrace();
                            }
                        }
                        else res[1] = Integer.MAX_VALUE;
                        try{
                            end_prev =  getStartEnd(nxt.get(i))[1];
                            }catch(Exception exc){
                                exc.printStackTrace();
                            }
                            i++;
                        return res;
                    }
                    public void remove() {}
                    
                };
            }
        };
      
    }
    
    public static int[] getStartEnd(File f)  throws Exception{
        int[] res = new int[2];
        BufferedReader br = AberationFinder.getBufferedReader(f, "snp.txt");
        br.readLine();
        res[0] = Integer.parseInt(br.readLine().split("\\s+")[2]);
        String st;
        String st1 = "";
        while( (st =  br.readLine())!=null){
            st1 = st;
            
        }
        res[1] = Integer.parseInt(st1.split("\\s+")[2]);
        br.close();
        return res;
    }
    
    public void run(int noCop) throws Exception{
        this.setNoCop(noCop);
        for(Iterator<Iterator<Object[]>> it = getIterator(user); it.hasNext();){
            Iterator<Object[]> chr  = it.next();
            try{
            //    System.err.println(chr);
            RunEHH rehh = new RunEHH(user,normal, score, this.scoring, clazz);
        
           // List<File> l1 = fileMap.get(chr);
            boolean first = true;
            while(chr.hasNext()){
                Object[] obj = chr.next();
                File f = (File) obj[0];
                if(first){
                    first = false;
                    rehh.setChrom(f.getName().split("_")[0]);
                }
             //   System.err.println("running with "+l.get(i));
                int nxtStart = (Integer) obj[1];
                rehh.run(f, nxtStart,noCop);
                
            }
            System.err.println("calculating scores!");
            rehh.calcScores();
            rehh.score.flush();
        //   this.printNormal();
         // System.err.println("done");
            }catch(Exception exc){
             
                exc.printStackTrace();
            }
        }
        score.close();
    }
    public void printNormal() throws Exception{
        File out = new File(user, "norm_"+noCop+"_0.txt");
        for (int i=0; out.exists() && out.length()>0; i++){
            out = new File(user, "norm_"+noCop+"_"+i+".txt");
        }
        System.err.println("normal output is "+out);
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
        for(Iterator<Double> it = normal.keySet().iterator(); it.hasNext();){
            Double  key = it.next();
            TrainableNormal norm = normal.get(key);
            if(norm.obsx.size()>0) norm.maximise(0, 0,0,0);
            pw.println(key+"\t"+norm.size()+"\t"+norm.toString());//+"\t"+norm.observations);
        }
        pw.close();
       
    }
    public static void main(String[] args){
        try{
            Options opt = Constants.OPTIONS;
            opt.addOption( OptionBuilder.withLongOpt( "hap" ).withDescription( "hap").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "baseDir" ).withDescription( "baseDir").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "overl" ).withDescription( "overl").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "cn" ).withDescription( "cn").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "input" ).withDescription( "input").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "graphAll" ).withDescription( "graphAll").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "extra" ).withDescription( "extra").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "chr" ).withDescription( "chr").withValueSeparator( ':' ).hasArgs().create());
            opt.addOption( OptionBuilder.withLongOpt( "start" ).withDescription( "start").withValueSeparator( ':' ).hasArgs().create());
           // opt.addOption( OptionBuilder.withLongOpt( "end" ).withDescription( "end").withValueSeparator( ':' ).hasArgs().create());
            Parser parser= new PosixParser();
            final CommandLine params = parser.parse(opt, args, false);
         
            EHHFinder.thresh = 0;
            EHHFinder.lim=  Integer.MAX_VALUE;
            EHHFinder.thresh1 = 0.0;
            File user = new File(System.getProperty("user.dir"));
            File baseDir = new File(user, params.getOptionValue("baseDir"));
          File input = null;
            if(params.hasOption("input")){
              input = new File(user, params.getOptionValue("input"));
          }
            if(params.hasOption("hap")){
            	String[] hap = params.getOptionValues("hap");
           	 RunEHH rehh = 
           		 new RunEHH(user, read(   new File(user, "norm.txt")), null, true, EHH.class);
           	 rehh.setChrom(params.getOptionValue("chr"));
           	 rehh.run(input, 
           			 Integer.parseInt(params.getOptionValue("start")), 
           		//	Integer.parseInt(params.getOptionValue("end")), 
           			hap[0], hap[1]);
           	return;
           }
            EHHMultiChromosome.overl = Boolean.parseBoolean(params.getOptionValue("overl", "false"));
            EHHMultiChromosome.graphAll = Boolean.parseBoolean(params.getOptionValue("graphAll", "false"));
         //  if(scoreDel){
            EHHMultiChromosome rehh = new EHHMultiChromosome(user, baseDir, input,Integer.parseInt(params.getOptionValue("extra", "0")));
            String[] cn = params.getOptionValues("cn");
           
            for(int i=0; i<cn.length; i++){
                System.err.println("running with "+ i +" cnv");
                    rehh.run(Integer.parseInt(cn[i]));
                    rehh.printNormal();
                }
     //   System.err.println("running with 2 cnv");
      //  rehh.run(1);
          
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   
}
