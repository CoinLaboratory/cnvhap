package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import lc1.CGH.AberationFinder;
import lc1.CGH.Location;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.sequenced.Convert;
import lc1.sequenced.Deletion;

public class EHHFinder1 implements Runnable {
    final static int column=1;
    static boolean restrict =true;
    static{
        try{
        ConsoleHandler handler = new ConsoleHandler();
        FileHandler handlerF = new FileHandler("stderr_EHHfinder_"+column+"_"+restrict, false);
        Logger.global.addHandler(handlerF);
        Formatter formatter = 
        new Formatter(){
            public String format(LogRecord record){
                return record.getSourceClassName()+":\n"+record.getMessage()+"\n";
            }
        };
        handler.setFormatter(formatter);
        handlerF.setFormatter(formatter);
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
        }
    }
    SimpleDataCollection sdt;
    File cgh;
   // Deletion[] dels;
    File[] files;
    public  EHHFinder1(File dir, File[] user, String[] args, PrintWriter out, 
            Map<String, Number> sum, String chromosome, Location reg) throws Exception{
        EHHFinder.thresh = 2;
        EHHFinder.lim=  Integer.MAX_VALUE;
        EHHFinder.thresh1 = 0.0;
        this.cgh = dir;
      //  dels = Convert.getDeletions(new File(cgh, "deletion_samples.txt"),
     //          null, null);
        this.files = user;
    }

    public void run() {
        try{
            File delF = new File(this.cgh, "deletion_samples.txt");
           
           
    //    Map<String, String> map1= new HashMap<String, String>();
        List<String> probs=  new ArrayList<String>();
    
        File sweepDir = new File( "sweep");
        if(sweepDir.exists()){
       // File[] ss = sweepDir.listFiles();
       // for(int ik=0; ik<ss.length; ik++){
         //   ss[ik].delete();
        //}
        }
        else sweepDir.mkdir();
        PrintWriter pw_out = new PrintWriter(new BufferedWriter(new FileWriter(new File(sweepDir, "sweep.many"))));
        try{      Logger.global.info("files "+Arrays.asList(files));
          for(int i=0;
          i<
        //1;
             files.length; 
          i++){ 
              try{
                  String[] del = EHHFinder.getMid(AberationFinder.getBufferedReader(files[i], "log_out.txt"+column), "mid");
                  String[] chr = EHHFinder.getMid(AberationFinder.getBufferedReader(files[i], "log_out.txt"+column), "chrom");
               //   Location delL = new Location(chr[0], Long.parseLong(del[0]), Long.parseLong(del[1]));
                  Location[] dels = Convert.getDeletions(delF,Integer.parseInt(chr[0]), new int[]{Integer.parseInt(del[0]), Integer.parseInt(del[1])});
            if(!chr[0].equals("1")) return;
              System.err.println("doing "+i);
              Logger.global.info("file  "+i+" of "+files.length);
              SimpleDataCollection sdt = 
                  SimpleDataCollection.readFastPhaseOutput((AberationFinder.getBufferedReader(files[i], "phased2.txt_"+column)),Emiss.class, Emiss.getEmissionStateSpace(1));
              sdt.removeKeyIfStartsWith("NA");
                  BufferedReader  snpFile  = AberationFinder.getBufferedReader(files[i], "snp.txt");
          //  if(snpFile.exists()){
                sdt.snpid= new ArrayList<String>();
                List<Integer> snp_orig = new ArrayList<Integer>();
                sdt.readPosInfo(snpFile, new int[] {0,4}, true, new List[] {sdt.snpid, snp_orig}, new Class[] {String.class, Integer.class});
                snpFile.close();
                Map<Integer, List<Integer>> summary = summary(snp_orig);
             //  if(!restrict) sdt.drop(summary.get(3));
            //}
            File outDir = files[i].getParentFile();
            if(outDir.getName().endsWith(".gz")){
                String nm  = outDir.getName();
                outDir = new File(outDir.getParentFile(), nm.substring(0, nm.indexOf(".tar.gz")));
                if(!outDir.exists())outDir.mkdir();
            }
         
              //files[i].getName().split("\\.")[0].split("_");
            //Deletion del = findDeletion(dels, sdt);
        //    map1.put(del.toStringShort(), files[i].getName());
        //    System.err.println("MAP "+map1);
          int st1 =sdt.loc.indexOf((int) dels[0].min);//sdt.loc.indexOf(Integer.parseInt(del[0]));
          int end1 = sdt.loc.indexOf((int) dels[0].max);//sdt.loc.indexOf(Integer.parseInt(del[1]));
          int sub =restrict ?  sdt.restrict(st1, end1, 50, 50, true) : 0;
         st1 -=sub;
        end1-=sub;
          String name = files[i].getName();
          if(name.indexOf('.')>=0){
              name = name.substring(0, name.indexOf("."));
          }
          if(restrict){
              Logger.global.info("hwe "+chr[0]+"_"+del[0]+"_"+del[1]);
           int blocksize = EHHFinder.writeEHH(sdt, sweepDir,chr[0], null, name);  
           
           Logger.global.info("block size is "+blocksize+" for "+files[i].getName());
          }
         //  sdt.calculateMaf(false);
           pw_out.println(name+".emphase\t"+name+".snp");
         //  if(true) continue;
            sdt.split();
        
            File hapF = new File("hap_"+column);
            if(!hapF.exists()) hapF.mkdir();
            
            PrintWriter pw = !restrict ? new PrintWriter(new BufferedWriter(new FileWriter(new File( hapF, "hap_summary"+files[i].getName()+".txt")))) : null;
            EHHFinder ehh = new EHHFinder(sdt);
            
          //  sdt.calculateMaf();
        if(restrict){  
            double[] maxR =  EHHFinder.calcLD(sdt, st1, end1, files[i].getName(), false,10, true);
       //  double maxR = max[Constants.getMax(max)];
       //  double maxd = max[1][Constants.getMax(max[1])];
          Logger.global.info("maxDprime is "+maxR[0]+" "+maxR[1]+" "+maxR[2]+" for "+files[i].getName());
          maxR =  EHHFinder.calcLD(sdt, st1, end1, files[i].getName(), false,1000, false);
          //  double maxR = max[Constants.getMax(max)];
          //  double maxd = max[1][Constants.getMax(max[1])];
             Logger.global.info("maxr2 is "+maxR[0]+" "+maxR[1]+" "+maxR[2]+" for "+files[i].getName());
        } 
        else{
          Set<String> coreHaps =  sdt.getAllHaplotypes(st1, end1).keySet();
           SortedMap<Integer, SortedMap<String, Object[]>> first = null;
          for(Iterator<String> it = coreHaps.iterator(); it.hasNext();){
          
              String st = it.next();
              Logger.global.info("HAPLOTYPE "+st);
              pw.println("HAPLOTYPES FOR "+st);
              ehh.base = st;
              System.err.println(st);
              boolean fir = st.startsWith("_") && st.endsWith("_");
              if(!fir) continue;
          //   if(fir ) ehh.set =getDel(dels, Integer.parseInt(del[0])).noObs;
            // else {
              //   ehh.master = first;
                //ehh.set = null;
            // }
            ehh.extend(st, st1, end1, pw);
            if(fir){
                first = ehh.map;
            }
            pw.println("##############################################");
            pw.flush();
            System.err.println("done "+i);
          }
          }
              }catch(Exception exc){
                  exc.printStackTrace();
                  probs.add(files[i].getName());
                  System.err.println("error with file "+files[i]);
              }
//        List<String> indiv = sdt.getNames();
              }
              
         // System.err.println("MAP FINAL"+map1);
          System.err.println("problems "+probs);
    }catch(Exception exc){
        exc.printStackTrace();
    }
    pw_out.close();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }

    private Map<Integer, List<Integer>> summary(List<Integer> snp_orig) {
       Map<Integer, List<Integer>> m = new HashMap<Integer, List<Integer>>();
       for(int i=0; i<snp_orig.size(); i++){
          Integer val = snp_orig.get(i);
          List<Integer> l = m.get(val);
          if(l==null){
              m.put(val, l = new ArrayList<Integer>());
          }
          l.add(i);
       }
       return m;
    }

    private Location getDel(Location[] dels, int mid) {
     for(int i=0; i<dels.length; i++){
         if((int) dels[i].min == mid) return dels[i];
     }
     return null;
    }

    private Deletion findDeletion(Deletion[] dels2,DataCollection dt) {
        List<Deletion>d = new ArrayList<Deletion>();
        int start = dt.loc.get(0);
        int end = dt.loc.get(dt.length()-1);
      for(int i=0; i<dels2.length; i++){
          if(dels2[i].start>= start && dels2[i].start <= end
                   && dels2[i].end>= start && dels2[i].end <= end
          ){
              if(dt.loc.indexOf(dels2[i].start)>=0){
              d.add(dels2[i]);
              }
          }
      }
      if(d.size()!=1)
      throw new RuntimeException( "!!WRONG NUMBER"+d);
      return d.get(0);
    }
}
