package lc1.CGH;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.logging.Logger;

import lc1.possel.EHHFinder1;
import lc1.possel.ExtendedHaplotypeFinder;
import lc1.possel.ExtractHalpotypeStructure;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class AberationFinder {
    
    public static double lowIntens=0.11;
    
    public static int noOfSnps =2;
    public static int[] maf = new int[] {1,1,1};
    public static int overlapThresh =2;
    public static int requiredNumberOfProbesForFP =2;
    public static double overlapFracForTrueMatch =0;
    public static double threshold = -1;
    public static double overlapThresh1 =0.0;
   public static int mode =0;
    //should we merge across individuals
    public static double merge = 0.5;
    public static boolean includeCN  = true;
    public static boolean includeName =true;
    public static boolean add185KProbes = false;
    public static boolean add244KProbes = true;
    public static boolean addADM1Probes = false;
    public static int[] reg = null;
    public static int mult = 0;
    public static int endEffects = 0;
    public static boolean includeFrench = false;
    public static String[] type = new String[] {"cnv.txt"};//, "FrenchSamples.txt"};
   public static boolean includeSingle = true;
    public static boolean randomLoc =false;
    public static boolean limit = true;
    
   static final Options OPTIONS  = new Options(){
       {
           Field[] f = AberationFinder.class.getFields();
           for(int i=0; i<f.length; i++){
               if(Modifier.isStatic(f[i].getModifiers())){
                   this.addOption( OptionBuilder.withLongOpt( f[i].getName() ).withDescription( f[i].getName()).withValueSeparator( ':' ).hasArgs().create());
               }
           }
       }
   };
   
   
   public static double merge(){
       return merge;
   }
   public static boolean includeCN(){
       return includeCN;
   }
   public static boolean includeName() {
       return includeName;
   }
   
  public static  FileFilter ff_top = new FileFilter(){
       public boolean accept(File pathname) {
          boolean acc =  pathname.isDirectory() && pathname.getName().indexOf("subm")<0 && pathname.getName().indexOf("CGH")<0
           && pathname.getName().indexOf("X")<0;
          if(acc){
              try{
              Integer.parseInt(pathname.getName());
          }
          catch(Exception exc){
              return false;
          }
       }
          return acc;
       }
   };
   
  static FileFilter ff_mid  = new FileFilter(){

       public boolean accept(File pathname) {
         if(!pathname.isDirectory()){
           if(pathname.getName().endsWith(".tar.gz") || pathname.getName().endsWith(".tar")) return true;
         }
         return false;
       }
       
   };
  
 
   
   public static Map<String, List<File>> getDirs1(File user){
       File[] files  = user.listFiles(ff_mid);
       Map<String, List<File>> f = new TreeMap<String, List<File>>();
       for(int i=0; i<files.length; i++){
           
           String nm = files[i].getName().split("_")[0];
           List<File> l = f.get(nm);
           if(l==null){
               f.put(nm, l = new ArrayList<File>());
           }
           l.add(files[i]);
       }
       return f;
   }
    public static File[] getDirs(File user){
        return user.listFiles(ff_top);
    }
   public static File[] getFiles(File user){
      
       
          File[] f =  user.listFiles(ff_mid);
          return f;
       
    }
  
    public static Object[] getChromDir(File user){
        int chrom;
        File dir= user.getParentFile();;
        try{
            chrom= Integer.parseInt(user.getName());
            File cgh = new File(dir, "CGH");
            if(!cgh.exists()) throw new RuntimeException("!!");
         }
        catch(Exception exc){
            chrom = Integer.parseInt(user.getParentFile().getName());
            dir = dir.getParentFile();
        }
        return new Object[] {chrom, dir};
    }
    public static Runnable makeAbFinder(File par, File[] user, String[] type, PrintWriter out, Map<String, Number> sum,
            String chrom, Location region) throws Exception{
      if(mode==0) 
            return  null;//new AbFinder(par,user , type, out, sum,chrom, region);
       else if(mode==1)
           return new ExtractHalpotypeStructure(par,user , type, out, sum,chrom, region);
       else if(mode==2)
           return new ExtendedHaplotypeFinder(par,user , type, out, sum,chrom, region);
       else if(mode==3) return new EHHFinder1(par,user , type, out, sum,chrom, region);
       else return null;//new ErrorRate(par, user, type, out, sum, chrom, region);
    }
    public static void main(String[] args){
       
        try{ 
            parse(args);
            Map<String, Number> sum = new TreeMap<String, Number>();
            if(!extract && overlapThresh1 > 0.0) throw new RuntimeException("!!");
            else if (extract && overlapFracForTrueMatch > 0.0)  throw new RuntimeException("!!");
            File user = new File(System.getProperty("user.dir"));
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File("abResults.txt"))));
            printOptions(out);
            Location region = reg!=null ? new Location(reg[0]+"", reg[1], reg[2]) : null;
            if(mult==3){
                File cgh = new File(user, "CGH");
                while(!cgh.exists()){
                    cgh = new File(cgh.getParentFile().getParent(), "CGH");
                    
                }
                Map<String, List<File>> files = getDirs1(user);
                List keys = new ArrayList<String> (files.keySet());
                Collections.reverse(keys);
                for(Iterator<String> it =keys.iterator(); it.hasNext();){
                    String key = it.next();
                    if(key.equals("X")) continue;
                    File[] fi = files.get(key).toArray(new File[0]);
                    int chrom = Integer.parseInt(key);
                    try{
                       // out.println("directory: "+files[j].getName());
                     //  Logger.global.info("files from "+files[j].getName());
                      //  if(chrom!=1) continue;
                        for(int i=0; i<fi.length; i++){
                        Runnable abF = makeAbFinder(cgh, new File[] {fi[i]}, type, out, sum, chrom+"", null);
                    abF.run();
                        }
                    }catch(Exception exc){
                        exc.printStackTrace();
                        //System.exit(0);
                    }
                }
            }
            else if(mult==2){
               File cgh = new File(user, "CGH");
               if(!cgh.exists()){
                   cgh = new File(user.getParent(), "CGH");
                   
               }
            File[] files = getDirs(user);
//                getFiles(user);
          
          // Arrays.fill(sum, 0.0);
                for(int j=0; j<files.length; j++){
                    int chrom = Integer.parseInt(files[j].getName());
                    try{
                        out.println("directory: "+files[j].getName());
                        Logger.global.info("files from "+files[j].getName());
                        Runnable abF = makeAbFinder(cgh, getFiles(files[j]), type, out, sum, chrom+"", null);
                    abF.run();
                    }catch(Exception exc){
                        exc.printStackTrace();
                        //System.exit(0);
                    }
                }
                
            }
            else if(mult==1){
                File cgh ;
                int chrom=0;
                if(mode!=3){
                    Object[] obj = getChromDir(user);
                     chrom = (Integer)obj[0];
                    File dir= (File) obj[1];
                    cgh = new File(dir, "CGH");
                }
                else{
                    cgh = new File(user, "CGH");
                    if(!cgh.exists()){
                        cgh = new File(user.getParent(), "CGH");
                        
                    }
                    if(!cgh.exists()){
                        cgh = new File(user.getParentFile().getParentFile(), "CGH");
                    }
                }
                try{
                    Runnable abF = makeAbFinder(cgh, getFiles(user), type, out, sum, chrom+"", null);
                    abF.run();
                }catch(Exception exc){
                    exc.printStackTrace();
                    //System.exit(0);
                }
            }
            else{
                
                String chrom =  user.getName();
                if(!chrom.equals("X")){
                    try{
                        Integer.parseInt(chrom);
                    }
                    catch(Exception exc){
                        chrom = user.getParentFile().getName();
                    }
                }
                File par = user.getParentFile();
                while(!(new File(par, "CGH")).exists()){
                    par = par.getParentFile();
                }
                Runnable abF = makeAbFinder(new File(par, "CGH"), new File[] {user}, type, out, sum,chrom, region);
                abF.run();
            }
            System.err.println("printing summary");
            out.println("summary");
            for(Iterator<Entry<String, Number>> it = sum.entrySet().iterator(); it.hasNext();){
                Entry<String,Number> nxt = it.next();
                out.println(String.format("%-40s : %20i", new Object[] {nxt.getKey(),nxt.getValue()}));
            }
            
       out.close();
           
        }catch(Exception exc){
            exc.printStackTrace();
            //System.exit(0);
        }
    }
    
    public static List<Integer> readPosInfo(File f, int index, boolean header) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        if(header) br.readLine();
        String st = "";
        List<Integer> res = new ArrayList<Integer>();
        while((st = br.readLine())!=null){
            String st1 = st.trim();
            String[] str = st1.split("\\s+");
         //   System.err.println(st1);
            if(str.length>index){
                try{
                res.add(Integer.parseInt(str[index]));
                }catch(Exception exc){
                    System.err.println(Arrays.asList(str));;
                    exc.printStackTrace();
                }
            }
        }
        br.close();
        return res;
    }
  
    private static Integer sum(SortedSet<Entry<Integer, Integer>> name) {
       int sum =0;
       for(Iterator<Entry<Integer, Integer>> it = name.iterator(); it.hasNext();){
           sum+= it.next().getValue();
       }
       return sum;
    }
    static long gap = 0;//(long)1e9;
   // Set<String> chromosomes = new HashSet<String>();

    public static boolean extract = true;
    
   
    public static BufferedReader getBufferedReader(File dir, String name1) throws Exception{
        BufferedReader br;
        /*if(dir.getName().endsWith(".tar.gz")){
            TarInputStream dir1 = new TarInputStream(new GZIPInputStream(new FileInputStream(dir)));
            // br = dir1;
             TarEntry nxt =  dir1.getNextEntry();
             while(nxt!=null){
                 String name = nxt.getName();
                 if(name.indexOf(name1)>=0){
                     break;
                 }
                 nxt = dir1.getNextEntry();
             }
             br = new BufferedReader(new InputStreamReader(dir1));
            // loc = nxt.getFile();
        }
        else if(dir.getName().endsWith(".tar")){
            TarInputStream dir1 = new TarInputStream(new FileInputStream(dir));
           // br = dir1;
            TarEntry nxt =  dir1.getNextEntry();
            while(nxt!=null){
                String name = nxt.getName();
                if(name.indexOf(name1)>=0){
                    break;
                }
                nxt = dir1.getNextEntry();
            }
            br = new BufferedReader(new InputStreamReader(dir1));
           // loc = nxt.getFile();
        }
        else{*/
           File loc = new File(dir,name1 ); 
           br   = new BufferedReader(new FileReader(loc));
        //}
        return br;
    }
    public static void printOptions(PrintWriter pw){
        try{
            Calendar cal = new GregorianCalendar();
            pw.println(cal.getTime());
            Field[] f = AberationFinder.class.getFields();
            for(int i=0; i<f.length; i++){
                if(Modifier.isStatic(f[i].getModifiers())){
                    if(f[i].getType().equals(Random.class)) continue;
                    else if(f[i].getType().equals(int[].class)){
                        int[] val =(int[]) f[i].get(null);
                        if(val!=null){
                            Integer[] d = new Integer[val.length];
                            StringBuffer sb = new StringBuffer();
                            for(int ik=0; ik<d.length; ik++){
                                d[ik] = val[ik];
                                sb.append("%5i ");
                            }
                            pw.println(f[i].getName()+" "+String.format(sb.toString(), d));
                        }
                    }
                    else{
                        pw.println(f[i].getName()+" "+f[i].get(null));
                    }
                }
            }
           
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }


    public static void parse(String[] args) throws Exception{
        Parser parser= new PosixParser();
        final CommandLine params = parser.parse(OPTIONS, args, false);
        Option[] options = params.getOptions();
        for(int i=0; i<options.length; i++){
            String argName = options[i].getLongOpt();
          Field field =  AberationFinder.class.getField(argName);
          Class type = field.getType();
          if(type.equals(double.class)){
              field.set(null, Double.parseDouble(options[i].getValue()));
          }
          else if(type.equals(double[].class)){
              field.set(null, parse(options[i]));
          }
          else if(type.equals(int[].class)){
              field.set(null, parse1(options[i]));
          }
          else if(type.equals(int.class)){
              field.set(null, Integer.parseInt(options[i].getValue()));
          }
          else if(type.equals(String.class)){
              String val = options[i].getValue();
            
              field.set(null, val);
          }
          else if(type.equals(boolean.class)){
              field.set(null, Boolean.parseBoolean(options[i].getValue()));
          }
          else if(type.equals(Boolean.class)){
              field.set(null, new Boolean(Boolean.parseBoolean(options[i].getValue())));
          }
          else if(type.equals(URL.class)){
              String val = options[i].getValue();
              if(argName.equals("dir")){
               String[] str = options[i].getValues();
               StringBuffer sb = new StringBuffer();
               for(int ik=0; ik<str.length; ik++){
                   sb.append(str[ik]);
                   if(ik<str.length-1) sb.append(":");
               }
               val = sb.toString();
               if(!val.startsWith("file") && !val.startsWith("http")){
                   File f = new File(val);
                   field.set(null, f.toURL());
               }
               else{
                   field.set(null, new URL(val));
               }
              }
          }
          else{
              throw new RuntimeException("unknown type "+argName);
          }
          System.err.println("Set "+field.getName()+" : "+field.get(null));
        }
}

    
    private static int[] parse1(Option params){
        String[] cop = params.getValues();
       int[]  noCopies = new int[cop.length];
        for(int i=0; i<cop.length; i++){
            noCopies[i] =Integer.parseInt(cop[i]);
        }
        return noCopies;
}
    private static int[] parse1(String[] cop, String def){
        //String[] cop = params.getValues();
        if(cop==null) cop = def.split(":");
       int[]  noCopies = new int[cop.length];
        for(int i=0; i<cop.length; i++){
            noCopies[i] =Integer.parseInt(cop[i]);
        }
        return noCopies;
}
    
    private static double[] parse(Option params){
        String[] cop = params.getValues();
       double[]  noCopies = new double[cop.length];
        for(int i=0; i<cop.length; i++){
            noCopies[i] = Double.parseDouble(cop[i]);
        }
        return noCopies;
}
    public static boolean limit() {
        if(AberationFinder.includeName()) return true;
        else return limit;
    }
}
