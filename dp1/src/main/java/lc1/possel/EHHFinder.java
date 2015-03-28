/**
 * 
 */
package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import lc1.CGH.Aberation;
import lc1.CGH.AberationFinder;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

import edu.mit.wi.haploview.HaploView;


public class EHHFinder{
    
    static{
        try{
        ConsoleHandler handler = new ConsoleHandler();
        FileHandler handlerF = new FileHandler("stderr_EHHfinder", false);
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
    Collection<String> set;
    //String base;
    SortedMap<Integer, SortedMap<String, Object[]>> map;
    SortedMap<Integer, SortedMap<String, Object[]>> master;
    
    public boolean masterContains(String st){
        if(master==null) return true;
        for(Iterator<SortedMap<String, Object[]>> it = master.values().iterator(); it.hasNext();){
            SortedMap<String, Object[]> nxt =  it.next();
            for(Iterator<String> it1 = nxt.keySet().iterator(); it1.hasNext();){
                if(it1.next().indexOf(st)>=0) return true;
//                if(it.next().containsKey(st)) return true;
            }
        }
        return false;
    }
   
    public static void main(String[] args){
        try{
            Options opt = Constants.OPTIONS;
            opt.addOption( OptionBuilder.withLongOpt( "hap" ).withDescription( "hap").withValueSeparator( ':' ).hasArgs().create());
            Parser parser= new PosixParser();
            final CommandLine params = parser.parse(opt, args, false);
            int col = Integer.parseInt(params.getOptionValue("column"));
            EHHFinder.thresh = 0;
            EHHFinder.lim=  Integer.MAX_VALUE;
            EHHFinder.thresh1 = 0.0;
            File user = new File(System.getProperty("user.dir"));
           String[] md =  params.getOptionValues("mid");
             
           final int[] mid = new int[md.length];
           String chr = md[0];
           for(int i=0; i<md.length; i++){
               mid[i] = Integer.parseInt(md[i+1]);
           }
           File[] f = user.listFiles(new FilenameFilter(){
            public boolean accept(File arg0, String arg1) {
                String[] str = arg1.split("_");
                int start = Integer.parseInt(str[0]);
                int end = Integer.parseInt(str[1]);
                if(mid[0]>start && mid[1]<end)
               return arg1.endsWith(".zip");
                else return false;
            }
           });
        
        SimpleDataCollection sdt = 
            SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader(user, "phased2.txt_"+col),Emiss.class, Emiss.getEmissionStateSpace(1));
        sdt.removeKeyIfStartsWith("NA");
      //  sdt.calculateMaf(false);
        int len = sdt.length();
       
        File snpFile  = new File(user, "snp.txt");
        if(snpFile.exists()){
            sdt.snpid= new ArrayList<String>();
            sdt.readPosInfo(snpFile, new int[] {0}, true, new List[] {sdt.snpid}, new Class[] {String.class});
        }
        if((sdt.loc.get(sdt.loc.size()-1) -sdt.loc.get(0))>100000){
            System.err.println(sdt.loc.size());
    int sub =    sdt.restrict(sdt.loc.indexOf(mid[0]), sdt.loc.indexOf(mid[1]), 50, 50,true);
        }
 //       st1 -=sub;
   //     end1-=sub;
        int[] nmid = new int[] {
                find(sdt.loc,mid[0]),
                find(sdt.loc,mid[1])};
        if(nmid[0]==-1){
            nmid[0] = nmid[1] = sdt.loc.indexOf(Constants.mid()[0]);
        }
        if(nmid[0]==-1 || nmid[1]==-1) throw new RuntimeException("!!");
     
       
       
        if(true){  
            SimpleDataCollection sdt1 = sdt.clone();
            sdt1.restrict(sdt.loc.indexOf(mid[0]), sdt.loc.indexOf(mid[1]), 50, 50, true);
            writeEHH(sdt1, new File(user, "sweep"),Constants.chrom0()+"", null, mid[0]+"-"+mid[1]
                                                                                              );   
            int st1 =nmid[0];
            int end1 = nmid[1];
       
         if(false){   double[] maxR =  EHHFinder.calcLD(sdt1, st1, end1, "", false,1000, false);
       //  double maxR = max[Constants.getMax(max)];
       //  double maxd = max[1][Constants.getMax(max[1])];
          Logger.global.info("maxR/D is "+maxR[0]+" "+maxR[1]+" "+maxR[2]+" for "+"");}
        } 
        sdt.split();
        Set<String> coreHaps =  sdt.getAllHaplotypes(nmid[0], nmid[1]).keySet();
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(user, "hap_summary.txt"))));
        EHHFinder ehh = new EHHFinder(sdt);
        
        ehh.set = new HashSet<String>();
        List<String > ll= new ArrayList<String>();
        for(Iterator<String> it = ll.iterator(); it.hasNext();){
            ehh.set.add(it.next().split("_")[0]);
        }
        System.err.println("set is "+ehh.set);
  
   String hap = params.getOptionValue("hap");
  // double[] ld = calcLD(sdt, nmid[0], nmid[1], mid[0]+"-"+mid[1], true, 100);
           // coreHaps =  sdt.getAllHaplotypes(nmid[0], nmid[1]).keySet();
           // System.err.println("haps "+coreHaps);
          // SortedMap<Integer, SortedMap<String, Object[]>> first = null;
           SortedMap<Integer, SortedMap<String, Object[]>> first = null;
          for(Iterator<String> it = coreHaps.iterator(); it.hasNext();){
              String st = it.next();
              pw.println("HAPLOTYPES FOR "+st);
              ehh.base = st;
              System.err.println(st);
              int l = hap.length();
              boolean fir = st.startsWith(hap.substring(0,1)) && st.endsWith(hap.substring(l-1, l));
            //  if(fir){
              //    ehh.set = getLoc(user, mid);
           //   }
           //   else {
            //      ehh.master = first;
            //     ehh.set = null;
           //   }
             ehh.extend(st, nmid[0], nmid[1], pw);
             if(fir){
                 first = ehh.map;
             }
           
          //  ehh.extend(st, st1, end1, pw);
     
            pw.println("##############################################");
            pw.flush();
          }
//            List<String> indiv = sdt.getNames();
      
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    private static int find(List<Integer> loc, int pos) {
        double min = Integer.MAX_VALUE;
        int min_i = 0;
        for(int i=0; i<loc.size(); i++){
            double dist = Math.abs(loc.get(i) - pos);
            if(dist < min){
                min_i = i;
                min = dist;
            }
        }
        return min_i;
    }
    public static String[] getMid(File file, String tag) throws Exception{
        BufferedReader br=  new BufferedReader(new FileReader(file));
        return getMid(br, tag);
    }
    static String[] getMid(BufferedReader br, String tag) throws Exception{
      String st = "";
      String tg = "--"+tag;
      while((st = br.readLine())!=null){
          String[] str = st.trim().split("\\s+");
          if(str[0].equals(tg)){
              str = str[1].split(":");
              return str;
//              int[] res = new int[2];
 //             for(int i=0; i<res.length; i++){
  //                res[i] = Integer.parseInt(str[i]);
   //           }
    //          br.close();
     //         return res;
          }
      }
      br.close();
      return null;
    }
public static int writeEHH(
        
        DataCollection sdt2, File f, String chr, Set<String> indiv, String name) throws Exception{
        if(!f.exists()) {
            f.mkdir();
        }
      
        EmissionState maf1 =null;// sdt2.calculateMaf1();
        Collection<Integer> toD = maf1.getConstantPos();
        toD.addAll(maf1.getMultPos());
        if(toD.size()>0) {
            Logger.global.info
           // throw new RuntimeException
            ("multi deletions "+chr+" "+name);
        }
        File emFile  = new File(f, name+".emphase");
        File snpFile = new File(f, name+".snp");
        int max = sdt2.loc.get(sdt2.loc.size()-1);
        int min = sdt2.loc.get(0);
        //int lenR = (int) Math.round( (double)(max - sdt2.loc.get(mid[0]))/1000.0);
        //int lenL = (int) Math.round( (double)(sdt2.loc.get(mid[1])-min)/1000.0);
        File pngFile = new File(f, name+".emphase.LD.PNG");//sdt2.loc.get(0)+""+sdt2.loc.get(sdt2.loc.size()-1));
        //File newPng = new File(f, name+
          //      "-"+lenL+"-"+lenR+"kb"+
            //    ".LD.PNG");//
        
        File track = new File(f, name+".track");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(track)));
       // pw.print(mid[0]+1);
      //  for(int i=mid[0]+2; i<=mid[1]+1; i++){
       //     pw.print(" "+i);
       // }
        pw.println();
        pw.close();
       // sdt2.snpid.set(mid[0], "deletion");
        //sdt2.snpid.set(mid[1], "deletion");
       //for(int i=mid[0]; i<mid[1]; i++){
       //    toD.add(i);
       //}
     //   if(mid[0] < mid[1]){
      //         toDrop.add(mid[0]);
     //        mid[1]--;
     //     }
       String blocktype = "GAB";
       if(true){ 
    	   sdt2.writeDickFormat(emFile, false, indiv, toD);
           sdt2.writeSNPFile(snpFile, chr, false, toD);
      // Logger.global.info("hwe "+f+" "+chr+"_"+mid[0]+"_"+mid[1]);
       HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
                "-info", snpFile.getAbsolutePath(),
    //            "-dprime",
               // "-blocks", track.getAbsolutePath(),
                "-blockoutput", blocktype,
               "-blockCutHighCI", "1.0",
               "-blockCutLowCI","1.0",
                "-ldcolorscheme", "RSQ",
            //  "-hwcutoff", "0",
     // "-png" ,"-n"
               });
        //File newPng = new File(f, name+
        //      "-"+lenL+"-"+lenR+"kb"+
          //    ".LD.PNG");//
     /* HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
                "-info", snpFile.getAbsolutePath(),
                "-dprime",
               // "-blocks", track.getAbsolutePath(),
                "-blockoutput", blocktype,
                //"-blockCutHighCI", "1.0",
                //"-blockCutLowCI","1.0"
                "-ldcolorscheme", "RSQ",
             //   "-hwcutoff", "0.1"
             "-png" ,"-n"
               });*/
       }
        
        
        sdt2.writeDickFormat(emFile, false, indiv, null);
        sdt2.writeSNPFile(snpFile, chr, true, null);
      //  Thread.sleep(1000);
        String blockExt = blocktype.equals("GAB") ? ".GABRIELblocks" :".4GAMblocks";
        File blockFile = new File(f, name+".emphase"+blockExt);
        int cnt =0;
        while(!blockFile.exists() && cnt < 3){
            Thread.sleep(1000);
            blockFile = new File(f, name+".emphase"+blockExt);
            cnt++;
        }
        if(blockFile.exists()){
        BufferedReader br = new BufferedReader(new FileReader(blockFile));
        String st = "";
       /* while((st = br.readLine())!=null){
            if(st.startsWith("BLOCK")){
                String[] str = st.split(":")[1].trim().split("\\s+");
                int[] res = new int[str.length];
                boolean contains = false;
                for(int jk=0; jk<res.length; jk++){
                    res[jk] = Integer.parseInt(str[jk]);
                    if(res[jk]==mid[0] || res[jk]==mid[1]) contains = true;
                }
                if(contains) return sdt2.loc.get(res[res.length-1]) - sdt2.loc.get(res[0]);
            }
        }*/
        }
        return 0;
       // if(true) Thread.currentThread().wait(1000000);
       // String[] cmd = new String[] {"mv", f.getName()+"/"+pngFile.getName(), f.getName()+"/"+newPng.getName()};
       // StringWriter err = new StringWriter();
       // ProcessTools.exec(cmd, null, err, err);
      //  System.err.println(err.toString());
         // TODO Auto-generated method stub
         
        // TODO Auto-generated method stub
    }

Comparator comp= new Comparator<String>(){
       public int compare(String o1, String o2) {
          int l1 = o1.length();
          int l2 = o2.length();
          if(l1!=l2) return l1> l2 ? -1 : 1;
          else return o1.compareTo(o2);
       }
         
     };
     static int thresh =3 ;
     static double thresh1 = 0.2;
    final SimpleDataCollection sdt;
    public EHHFinder(SimpleDataCollection sdt2) {
        this.sdt = sdt2;
    }
    private Set<String> getIndividuals(List<PIGData> value) {
        Set<String> res = new TreeSet<String>();
        for(int i=0; i<value.size(); i++){
            res.add(value.get(i).getName());
        }
        return res;
      }
    int orig_size;
    public void extend(Aberation ab, PrintWriter pw){
         String[] nme = ab.name.split("_");
         PIGData haplotype = sdt.get(nme[0]).split()[Integer.parseInt(nme[1])];
          String string = haplotype.getStringRep(ab.start, ab.end);
          extend(string, ab.start, ab.end, pw);
    }
    String base;
    
    public void extend(String string, int start_, int end_,  PrintWriter pw){
     
        List<PIGData> dat = sdt.getAllHaplotypes(start_, end_).get(string);
       /// check(dat);
        if(dat==null) return;
        if(set!=null){
            SortedSet<String> hapl_names = new TreeSet<String>();
            for(Iterator<PIGData> it = dat.iterator(); it.hasNext();){
                PIGData dat_ = it.next();
                String name = dat_.getName().split("_")[0].split("\\|")[0];
               // if(name.startsWith("NA")) it.remove();
                
                 if(!set.contains(name)){
                 //   CompoundScorableObject dat1_ = (CompoundScorableObject) dat_.clone();
                  //  dat1_.restrictSites(start_-5, end_+5);
                    pw.println("deletion_samples.txt does not contain "+name+" "+dat_.getStringRep(Math.max(0, start_-5), Math.min(dat_.length()-1, end_+5)));
                   // it.remove();
                 }
                else hapl_names.add(name);
            }
            for(Iterator<String>  it = set.iterator(); it.hasNext();){
                String key = it.next();
              //  if(!hapl_names.contains(key)){
              //      pw.println("phased deletion does not contain "+key);
              //  }
            }
        }
      orig_size = dat.size();
      pw.println("original size is "+orig_size);
     map = new TreeMap<Integer, SortedMap<String, Object[]>>();
   //   System.err.println("HEEEERE");
   //  check(dat);
      extend(start_, end_, string, dat, map, 0, string.length());
      pw.println("aberation "+string+" "+sdt.loc.get(start_)+" "+sdt.loc.get(end_));
   
      int cnt =0;
      for(Iterator<Integer> it = map.keySet().iterator(); it.hasNext() ;cnt++){
          Integer key1 = it.next();
          SortedMap<String, Object[]> m1 = map.get(key1);
          List<PIGData> vals = null;
          for(Iterator<String> it1 = m1.keySet().iterator(); it1.hasNext();){
          String key = it1.next();
          Object[] val = map.get(key1).get(key);
          List<PIGData> val0 = (List<PIGData>) val[0];
          if(//val0.size() < 0.5*(double)orig_size && 
                  vals!=null && val0.equals(vals)) continue;
          else vals = val0;
          int start = (Integer) val[1];
          int end = (Integer) val[2];
        //  String st = nxt.getKey();
          PIGData da =SimpleScorableObject.make("", Arrays.asList(new String[] {key}),     Emiss.getEmissionStateSpace(1),(short)-1);
          StringWriter sw = new StringWriter();
          PrintWriter pw1 = new PrintWriter(sw);
        //  da.print(pw1, false, true, true, null, null, null);
          String st1 = sw.toString();
          Set<String> indv = getIndividuals(val0);
          if(this.master==null && val0.size() < 0.5 *(double) orig_size) continue;
       //  if(cnt<5 || val0.size()>=4) {
             pw.println(base+" : "+da.length()+"\tr"+val0.size()+" of "+orig_size+"\tl "+(sdt.loc.get(end) - sdt.loc.get(start))+"\t"+indv);
             pw.println(st1.substring(0, st1.indexOf("\n")-1)+"\t"+base+"\t"+val0.size()+"\t"+((double)val0.size() / (double)orig_size)+"\t"+((double)(sdt.loc.get(end) - sdt.loc.get(start))/1000.0));
       //  }
            /* if(key.length()>1){
                 String keyL = key.substring(1);
                 String keyR = key.substring(0, key.length()-1);
                inner: for(int key11 = key1; ; key11++){
                     SortedMap<String, Object[] > m11 = map.get(key11);
                     if(m11==null) break inner;
                     if(m11.containsKey(keyL)){
                         Set<String> valL = getIndividuals((List<PIGData>) m11.get(keyL)[0]);
                         valL.removeAll(indv); 
                         pw.println("\t\tleft lo "+valL);
                         break inner;
                     }
                }
                 inner: for(int key11 = key1; ; key11++){
                     SortedMap<String, Object[] > m11 = map.get(key11);
                     if(m11==null) break inner;
                     if(m11.containsKey(keyR)){
                         Set<String> valR = getIndividuals((List<PIGData>) m11.get(keyR)[0]);
                         valR.removeAll(indv); 
                         pw.println("\t\tright lo "+valR);
                         break inner;
                     }
                } 
                  
             }*/
          }
          //+map+"\n");
      }
    //    System.err.println(string+" "+dat.size()+" "+sdt.getKeys().size()+" "+left.getValue().size()+" "+right.getValue().size());

    }
   private void print(SortedMap<Integer, SortedMap<String, Object[]>> map, Integer key1, String key) {
        // TODO Auto-generated method stub
        
    }
   
   public static double[] calcLD(SimpleDataCollection sdt2, int st1, int end1, String name, boolean b, int numPerm, boolean dprime) throws Exception {
       LDInner inner = new LDInner(sdt2, b, numPerm, dprime);
       inner.setCore(st1, end1);
       return inner.run();
    }
/* public static void check(Collection<PIGData> dat){
        for(Iterator<PhasedIntegerGePXJ4*4*v
        notypeData> it = dat.iterator(); it.hasNext();){
            if(it.next().length()<35) throw new RuntimeException("!!");
        }
    }*/
   static Integer lim = 1;
    public void extend(int start, int end, String key, List<PIGData> dat,
            SortedMap<Integer, SortedMap<String, Object[]>>map, int base_index0, int base_index1){
         if(dat.size() < thresh || dat.size() < thresh1 * orig_size) return;
      
     //   else if(key.length()>10){
            SortedMap<String, Object[]> map1  = map.get(dat.size());
             if(map1==null) map.put(dat.size(), map1 = new TreeMap<String, Object[]>(comp));
             char[] ch = key.toCharArray();
             Arrays.fill(ch, base_index0, base_index1, '*');
             String key1 =new String(ch);
            if( map1.containsKey(key1)) return;
            if(!this.masterContains(key1)) return;
           // if(base_index1 - base_index0>1){
             //   throw new RuntimeException("!!");
          //  }
            map1.put(key1, new Object[] {dat, start, end});
       // }
     //   System.err.println(key+" "+dat.size());
        if(end < dat.get(0).length()-1){
            int k=0;
            for(Iterator<Entry<String, List<PIGData>>> it =
                sdt.getHaplotypes(sdt.getHaplotypes(start,end+1, dat)).iterator();
            it.hasNext() && k<lim;k++){
               Entry<String, List<PIGData>> right= it.next();
               extend(start, end+1, right.getKey(), right.getValue(), map, base_index0, base_index1);
            }
        }
        if(start-1>0){
            int k=0;
            for(Iterator<Entry<String, List<PIGData>>> it =sdt.getHaplotypes(sdt.getHaplotypes(start-1,end, dat)).iterator();
            it.hasNext()&& k<lim;k++){
            Entry<String, List<PIGData>> left = it.next();
            extend(start-1, end, left.getKey(), left.getValue(), map, base_index0+1, base_index1+1);
            }
           
        }
    
    }
    
}