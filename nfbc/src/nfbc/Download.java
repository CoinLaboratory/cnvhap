package nfbc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JOptionPane;
import javax.swing.JPasswordField;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;


public class Download {
    
    public static final Options OPTIONS  = new Options(){
        {
            Field[] f = Download.class.getFields();
            for(int i=0; i<f.length; i++){
                if(Modifier.isStatic(f[i].getModifiers())){
                    this.addOption( OptionBuilder.withLongOpt( f[i].getName() ).withDescription( f[i].getName()).withValueSeparator( ':' ).hasArgs().create());
                }
                else{
                    System.err.println("excluded "+f[i]);
                }
            }
        }
    };
    
    public static String[] user;
    public static String[] password;
    public static String[] fields;
    public static String[] experiment;
    
    public static void setOptions(String[] args){
            try{
            Parser parser= new PosixParser();
            final CommandLine params = parser.parse(OPTIONS, args, false);
            Field[] f = Download.class.getFields();
            for(int i=0; i<f.length; i++){
                if(Modifier.isStatic(f[i].getModifiers()) && !f[i].getName().equals("OPTIONS")){
                        try{    
                            String[] val = params.getOptionValues(f[i].getName());
                            f[i].set(null, params.getOptionValues(f[i].getName()));
                        }catch(Exception exc){
                            exc.printStackTrace();
                        }
                }
               
            }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    };
    
    File[] f1;
    File dir;
    
    Map<String, PrintWriter> output = new HashMap<String, PrintWriter>();
    static FileFilter assoc = new FileFilter(){

        public boolean accept(File arg0) {
          return  arg0.getName().endsWith(".dat");
        }
        
    };
    static FileFilter cover = new FileFilter(){

        public boolean accept(File arg0) {
          return  arg0.getName().startsWith("cover");
        }
        
    };
    
    
    public static void ftp() {
        boolean error = false;
        FTPDownload ftp=null;
        try{
            String user =  (String) JOptionPane.showInputDialog(null,
                    "user name", "user name", JOptionPane.QUESTION_MESSAGE, null, null,
                   "");
            
            JPasswordField pwd = new JPasswordField(10);
            int action = JOptionPane.showConfirmDialog(null, pwd,"Enter Password",JOptionPane.OK_CANCEL_OPTION);
            
            ftp = new FTPDownload("fh-bcsnp3.sm.med.ic.ac.uk", user,new String(pwd.getPassword()));
            experiment =  new String[] { (String) JOptionPane.showInputDialog(null,
                    "experiment", "experiment", JOptionPane.QUESTION_MESSAGE, null, null,
            "")};
           for(int i=0; i<experiment.length; i++){
               ftp.getFiles(experiment[i]);
           }
            ftp.disconnect();
        }catch(Exception exc){
            exc.printStackTrace();
            if(ftp!=null)  ftp.disconnect();
        }
        
    }
    public static void main(String[] args){
        try{
            if(true){
                setOptions(args);
                ftp();
            }
            //else{
                File f = new File(  System.getProperty("user.dir"));
                File[] fi = f.listFiles();
                for(int i=0; i<fi.length; i++){
                    Download dc = new Download(fi[i]);
                        dc.join();
                }
           //}

        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   File outputDir;
   
   public void join() throws Exception{
    
       File[] f = dir.listFiles(new FileFilter(){
           public boolean accept(File arg0) {
               return arg0.isDirectory();
           }
          });
       for(int i=0; i<f.length; i++){
           merge(f[i], true);
       }
       List<String[]> header = new ArrayList<String[]> ();
       List<BufferedReader> l = new ArrayList<BufferedReader>();
    
       String[] proto_header = null;
       for(int i=0; i<f.length; i++){
           BufferedReader br = new BufferedReader(new FileReader(new File(f[i], "summ.txt")));
           String[] st = br.readLine().split("\t");
           if(proto_header==null || proto_header.length<st.length){
               proto_header = st;
           }
           header.add(st);
           l.add(br);
         
       }
       PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "summ.txt"))));
       for(int i=0; i<proto_header.length; i++){
           String toP = proto_header[i];
           if(toP.equals("SNP")){
              toP = "MARKER";
          }
           if(toP.startsWith("phen")){
               toP = "TRAIT";
           }
           
           pw.print(toP+(i==proto_header.length-1 ? "\n" :"\t"));
       }
       for(int i=0; i<f.length; i++){
           int[] alias = getAlias(proto_header,header.get(i));
           String st = "";
           BufferedReader bri = l.get(i);
           while((st= bri.readLine())!=null){
               String[] str_i = st.split("\t");
               for(int j=0; j<alias.length; j++){
                   if(alias[j]<0) pw.print("");
                   else pw.print(str_i[alias[j]]);
                   if(j==alias.length-1) pw.println();
                   else pw.print("\t");
               }
           }
       }
       
       
       pw.close();
       for(int i=0; i<l.size(); i++){
          l.get(i).close();
       }
   }
   private int whichIndex(String st, String[] str){
       for(int i=0; i<str.length; i++){
           if(str[i].equals(st)) return i;
       }
       return -1;
   }
   
   private int[] getAlias(String[] proto_header, String[] strings) {
       int[] res = new int[proto_header.length];
       for(int i=0; i<proto_header.length; i++){
           res[i] = whichIndex(proto_header[i], strings );
       }
       return res;
}
/** merges files*/
   public void merge(File dir, boolean avoidRepeat) throws Exception{
       File output = new File(dir, "summ.txt");
       if(avoidRepeat && output.exists() && output.length()>0) {
           System.err.println("already done "+dir);
           return;
           
       }
       System.err.println("doing "+dir);
       PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(output)));
       File[] f = dir.listFiles(new FileFilter(){
        public boolean accept(File arg0) {
            return arg0.getName().endsWith("dat");
        }
       });
       List<String[]> header = new ArrayList<String[]> ();
       List<BufferedReader> l = new ArrayList<BufferedReader>();
       List<List<Integer>> cols = new ArrayList<List<Integer>>();
       List<String> headers = new ArrayList<String>();
       for(int i=0; i<f.length; i++){
           BufferedReader br = new BufferedReader(new FileReader(f[i]));
           String[] st = br.readLine().split("\t");
           if(st[0].equals("SNP")){
               l.add(br);
               header.add(st);
           }
           else{
               br.close();
           }
       }
       for(int i=0; i<header.size(); i++){
           String[] h = header.get(i);
           List<Integer> col_i = new ArrayList<Integer>();
           cols.add(col_i);
           for(int j=0; j<h.length; j++){
              if(!headers.contains(h[j])){
                  col_i.add(j);
                  headers.add(h[j]);
              }
           }
       }
       String[] nxtString = new String[l.size()];
       pw.print("TRAIT \t");
       for(int i=0; i<headers.size(); i++){
           pw.print(headers.get(i)+(i==headers.size()-1 ? "\n" :"\t"));
       }
       String snp_id= null;
       while(getNextString(l, nxtString)){
           pw.print(dir.getName().trim()+"\t");//+dir.getParentFile().getName().trim()+"\t");
       for(int i=0; i<nxtString.length; i++){
           
           String[] nString = nxtString[i].split("\t");
           if(i==0) snp_id = nString[0];
           else if(!snp_id.equals(nString[0])){
               throw new RuntimeException("!!");
           }
           for(Iterator<Integer> col = cols.get(i).iterator(); col.hasNext();){
               Integer column = col.next();
               pw.print(nString[column]+(col.hasNext() || i<nxtString.length-1 ? "\t": "\n"));
           }
       }
     //  pw.println();
       }
       pw.close();
       for(int i=0; i<l.size(); i++){
          l.get(i).close();
       }
   }
    private boolean getNextString(List<BufferedReader> l, String[] nxtString) throws Exception {
        boolean isNull = false;
        boolean allNull = true;
        for(int i=0; i<l.size(); i++){
            nxtString[i] = l.get(i).readLine();
            isNull = isNull || nxtString[i]==null;
            allNull = allNull && nxtString[i]==null;
        }
        if(isNull && !allNull) throw new RuntimeException("diff lenghts");
        if(isNull) return false;
        else return true;
    }
    public void run() throws Exception{
        Set<String> pheno = new HashSet<String>();
        PrintWriter err = new PrintWriter(new BufferedWriter(new FileWriter("stderr.txt")));
        PrintWriter out1 = new PrintWriter("summary.txt");
       Set<String> done = new HashSet<String>();
        for(Iterator<String> it = this.assocFiles.keySet().iterator(); it.hasNext();){
            String str = it.next();
          
            File cover = coverFiles.get(str);
            File assoc = assocFiles.get(str);
            File perm = permFiles.get(str);
            String[] info = getInfo(cover);
            String nme = info[0].trim()+"_"+info[1].trim()+".txt";
            if(done.contains(nme)) throw new RuntimeException("!!");
            done.add(nme);
            File outF = new File(outputDir, nme.replace(' ', '_'));
            if(outF.exists()) throw new RuntimeException("!!");
            PrintWriter out2 = new PrintWriter(new BufferedWriter(new FileWriter(outF)));
            err.println("found "+Arrays.asList(info));
            PrintWriter out = this.output.get(info[0]);
            BufferedReader br = new BufferedReader(new FileReader(assoc));
            BufferedReader br1 = new BufferedReader(new FileReader(assoc));
            BufferedReader br2 = new BufferedReader(new FileReader(assoc));
            String st = br.readLine();
            String st1 = br1.readLine();
            String st2 = br2.readLine();
            out2.print("phen\t"); out2.println(st2);
            if(out==null){
                String outname = info[0].replace(' ', '_');
                output.put(info[0], out = new PrintWriter(new BufferedWriter(new FileWriter(outname))));
                out.print("phen\t"); out.println(st);
                out1.print("phen\t"); out1.println(st1);
              
               // out2.close();
               // System.exit(0);
            }
            append(br, out, info[1]);
            append(br1, out1, info[1]);
            append(br2, out2, info[1]);
            br.close();
            br1.close();
            out2.close();
        }
        for(Iterator<PrintWriter> it = this.output.values().iterator(); it.hasNext();){
            it.next().close();
        }
        out1.close();
        err.close();
        
    }
    
    Map<String, File> assocFiles = new HashMap<String, File>();
    Map<String, File> permFiles = new HashMap<String, File>();
    Map<String, File> coverFiles = new HashMap<String, File>();
 
    Download(File f) throws Exception{
        dir = f;
      //  outputDir = new File(f, "output");
    /*    if(!outputDir.exists()) outputDir.mkdir();
       //  out = new PrintWriter(new BufferedWriter(new FileWriter()));
        List<File> datFiles = new ArrayList<File>(Arrays.asList(f.listFiles(assoc)));
        File[] coverFiles = f.listFiles(cover);
        for(int i=0; i<coverFiles.length; i++){
            String name = coverFiles[i].getName();
            name = name.substring(name.indexOf("cover.")+6);
            List<File> df = new ArrayList<File>();
            for(int j=0; j<datFiles.size(); j++){
                if(datFiles.get(j).getName().startsWith(name)){
                    df.add(datFiles.remove(j));
                }
            }
            int[] ind = getIndices(df);
            this.coverFiles.put(name, coverFiles[i]);
            if(ind[1]>=0) this.permFiles.put(name, df.get(ind[1]));
            else{
                System.err.println("no permutation for "+name);
            }
            if(ind[0]>=0) 
                
            this.assocFiles.put(name, df.get(ind[0]));
            else{
                System.err.println("no assoc for "+name);
            }
        }
        //names = new String[f1.length];
       */
        
    }
    private int[] getIndices(List<File> df) throws Exception {
        int[] res = new int[] {-1,-1};
        for(int i=0; i<df.size(); i++){
            BufferedReader br = new BufferedReader(new FileReader(df.get(i)));
            String st =  br.readLine();
            int len =  st.split("\t").length;
            if(len==5) res[1] = i;
            else if(len ==8) res[0] = i;
            br.close();
        }
        return res;
    }

    private void append(BufferedReader br, PrintWriter out, String phen) throws Exception{
       String st = "";
        while((st = br.readLine())!=null){
            out.print(phen+"\t");
            out.println(st);
        }
        out.flush();
    }
    private String[] getInfo(File cover2) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(cover2));
        String st = "";
        while((st = br.readLine())!=null){
           
            if(st.indexOf("Job description")>=0){
                String[] str = st.split(":")[1].split("-");
                return str;
            }
        }
        return null;
    }
    
}
