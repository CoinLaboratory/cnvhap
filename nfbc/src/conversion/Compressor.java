package conversion;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.Adler32;
import java.util.zip.CheckedInputStream;
import java.util.zip.CheckedOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;


public class Compressor {
    
    
    static int idColumn = 0;
    static int[] snpColumns =// new int[] {0,1};
                         new int[] {0}; //rs, chrom, position
    static int[] sortOrder = new int[] {2,3,1};
    static int[] resultColumns =new int[] {2,3};// new int[] {4,5};
    static int id_index = 1;
    static String[] resultColumnHeader = 
        new String[] {"Log R", "AgilentPred"}; 
    static String split = "\t";
    final File in;
    File infinium;
    File out; 
    static boolean header =true;
    
   
    
    public static void main1(final String[] args){
        try{
            Compressor comp = new Compressor();
            File f = new File("Desir");
            File[] fs = f.listFiles(new FilenameFilter(){

                public boolean accept(File arg0, String arg1) {
                   return arg1.endsWith("zip");
                }
                
            });
            for(int i=0; i<fs.length; i++){
                String nme = fs[i].getName();
                nme = nme.substring(0,nme.indexOf(".zip"));
                comp.join(new File[] {f, 
                            new File("DNID"), new File("Corbeil")},
                            new File( "diab_all"), nme);
            }
            
            
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    /** args 0 is the string the file should include */
    public   static void main(final String[] args){
       try{
        File user = new File(System.getProperty("user.dir"));
        Compressor.unzip(new File(user, args[0]));
       }catch(Exception exc){
           exc.printStackTrace();
       }
    }
   
    private static File[] getOutFiles( File user, List<String> types ){
        File [] out = new File[types.size()];
        for(int i=0; i<out.length; i++){
            out[i] = new File(user, types.get(i).replaceAll("\\s+", "_"));
        }
        return out;
    }
    
   private void writeLine(OutputStreamWriter osw, String[] str) throws Exception{
       osw.write(str[resultColumns[0]]);
       for(int i=1; i<resultColumns.length; i++){
           osw.write("\t");
           int j = resultColumns[i];
           osw.write(""+str[j]);
       }
       osw.write("\n");
   }
 
   
 public static void writeLine(OutputStreamWriter osw2, String[] str,
           int[] second, String string) throws Exception{
       osw2.write(str[second[0]]);
       for(int i=1;  i< second.length; i++){
         osw2.write("\t");
         osw2.write(str[second[i]]);
       }
       osw2.write(string);
       
   }
  public static void writeLine1(OutputStreamWriter osw, String[] str) throws Exception{
       osw.write(str[0]);
       for(int i=1; 
     i< (str.length);
       i++){
           osw.write("\t");
         osw.write(str[i]);
        
       }
       osw.write("\n");
   }
  public static void writeLine2(OutputStreamWriter osw, String[] str) throws Exception{
      osw.write(str[0]);
      for(int i=1; 
    i< (str.length);
      i++){
          osw.write("\t");
        osw.write(str[i]);
       
      }
      osw.write("\t");
  }
   private void writeLine(OutputStreamWriter osw, List<String> ids) throws Exception{
       for(int i=0; i<ids.size(); i++){
           osw.write(ids.get(i)+"\n");
       }
   }
 
    static final int BUFFER = 2048;
    
    private String getEntryName(String[] str){
        StringBuffer sb = new StringBuffer(str[snpColumns[0]]);
        for(int i=1; i<snpColumns.length; i++){
            sb.append("_");
            sb.append(modify(str[snpColumns[i]]));
        }
        return sb.toString();
    }
    
    /** gets all rs numbers which lie between start and end */
    public List<String[]> getRS(ZipFile f, int start, int end, String chr) throws Exception{
       BufferedReader br = new BufferedReader(new InputStreamReader( f.getInputStream(f.getEntry("SNPS"))));
       String st = "";
       List<String[]> l = new ArrayList<String[]>();
       while((st = br.readLine())!=null){
           String[] str = st.split("\t");
           if(!chr.equals(str[1])) throw new RuntimeException("!!");
           int sta = Integer.parseInt(str[2]);
           if(str.length>=4){
               sta =(int) Math.round( ((double)sta+ (double) Integer.parseInt(str[3]))/2.0);
           }
           if(sta<=end && sta>=start){
              l.add(str);
           }
       }
       br.close();
       return l;
    }
    public void extract(ZipFile f, List<String[]> rs, Map<String, List<String[]>> results) throws Exception{
       
        for(int i=0; i<rs.size(); i++){
            String rs_i = rs.get(i)[0];
         
            List<String[]> results_i = results.get(rs_i);
            if(results_i==null){
                results.put(rs_i, results_i = new ArrayList<String[]>());
                BufferedReader br = new BufferedReader(new InputStreamReader( f.getInputStream(f.getEntry(rs_i))));
                String st = "";
                while((st = br.readLine())!=null){
                    String[] st1 = st.split("\t");
                    String[]st2 = new String[st1.length+1];
                    System.arraycopy(st1, 0, st2, 0, st1.length);
                    st2[st1.length] = "0";
                    results_i.add(st2);
                    
                }
                br.close();
            }
        }
    }
    
    public void write(Map<String, List<String[]>> map, Map<String, String[]> detMap, List<List<String>> details) throws Exception{
       Map<Integer, String> mm = new TreeMap<Integer, String>();
       for(Iterator<String> it = map.keySet().iterator(); it.hasNext();){
           String key = it.next();
           int pos =Integer.parseInt(detMap.get(key)[2]);
           mm.put(pos, key);
           
       }
        for(Iterator<String> it = mm.values().iterator(); it.hasNext();){
            String key = it.next();
            details.add(Arrays.asList(detMap.get(key)));
          //  System.err.println("writing "+key);
            ZipEntry entry = new ZipEntry(key);
            outS.putNextEntry(entry);
            for(Iterator<String[]> it1 = map.get(key).iterator(); it1.hasNext();){
                writeLine1(osw, it1.next());
            }
            osw.flush();
        }
    }
    public double convertToE(double baseFrom, double y){
        return Math.log(Math.pow(baseFrom, y));
    }
    public void compressSpreadsheet(){
        try{
        BufferedReader origin = new BufferedReader(new InputStreamReader(
                        new FileInputStream(new File(out.getParentFile(), "AllProbeIntervals.txt"))));
       // BufferedReader single = new BufferedReader(new InputStreamReader(
      //                          new FileInputStream(new File(out.getParentFile(), "SingleProbeIntervals.txt"))));
        Map<String, List<String[]>> map = new HashMap<String, List<String[]>>();
        out = new File(out.getParentFile(), out.getName().substring(0, out.getName().indexOf(".txt.zip")));
        File indir = new File(out.getParentFile(), "all_cgh");
        resultColumnHeader = new String[] {"logR",  "CN"};
        origin.readLine();
        Map<String, Integer> idMap = new HashMap<String, Integer>();
        String st = "";
        String comp = null;
        int index =   snpColumns[0];
        String chr = "";
        int chr_index = 9;
        ZipFile inF=null;
        Map<String, String[]> detMap = new HashMap<String, String[]>();
        for(int ik=0; (st = origin.readLine())!=null; ik++ ){
             String[] str = st.split(split);
             String id = str[6];
             
             if(!str[chr_index].equals(chr)  || ik==20){
               System.err.println(str[chr_index]);
               File fIn = new File(indir, str[chr_index]+".zip");
                 inF = new ZipFile(fIn);
                // if(fIn.exists()) throw new RuntimeException("!!");
                 if(chr.length()>0){
                     write(map, detMap, details);
                     this.finishChrom();
                     map.clear();
                     detMap.clear();
                     
                     
                 }
                 else{
                     this.ids = this.getIndiv(inF, "Samples");
                     for(int ij=0; ij<ids.size(); ij++){
                         idMap.put(ids.get(ij),ij);
                     }
                 }
                 chr= str[chr_index];
                 FileOutputStream dest = getOS(new File(out, chr+".zip"));
                 checksum = new   CheckedOutputStream(dest, new Adler32());
                 outS = new 
                 ZipOutputStream(new 
                   BufferedOutputStream(checksum));
                   osw = new OutputStreamWriter(outS);
                   outS.setMethod(ZipOutputStream.DEFLATED);
                   ZipEntry headings = new ZipEntry("Name");
                   outS.putNextEntry(headings);
                   writeLine1(osw, resultColumnHeader);
                   osw.flush();
               
                
                 //}
             }
         //    System.err.println(Arrays.asList(str));
             List<String[]> affected =  this.getRS(inF, Integer.parseInt(str[10]), Integer.parseInt(str[11]), chr);
             for(int i=0; i<affected.size(); i++){
                 detMap.put(affected.get(i)[0], affected.get(i));
             }
             this.extract(inF, affected, map);
             for(Iterator<String[]> aff_it = affected.iterator(); aff_it.hasNext();){
                 List<String[]> str1 = map.get(aff_it.next()[0]);
                 int id_ind = idMap.get(id);
                 str1.get(id_ind)[1] =String.format("%5.3g", ( Double.parseDouble(str[14])/1.0));
             }
          
         }
        write(map, detMap, details);
        this.finishChrom();
        }catch(Exception exc){
            exc.printStackTrace();
        }
        
    }
    public Compressor(){
        in = null;
    }
    
public List<Integer> getMap(List<String> l1, List<String> l2){
    List<Integer> m = new ArrayList<Integer>();
   /* for(int i=0; i<l1.size(); i++){
        l1.set(i, l1.get(i).split("#")[0]);
    }
    for(int i=0; i<l2.size(); i++){
        l2.set(i, l2.get(i).split("#")[0]);
    }*/
    for(int i=0; i<l1.size(); i++){
        int ind = l2.indexOf(l1.get(i));
        if(ind<0) {
            throw new RuntimeException("!!");
        }
        m.add(ind);
    }
    return m;
}

public static FileOutputStream getOS(File f) throws Exception{
 // if(f.exists() && f.length()>0) throw new RuntimeException("!!" +f.getAbsolutePath());
    return new FileOutputStream(f);
}

public static FileOutputStream getOS1(File f) throws Exception{
 // if(f.exists() && f.length()>0) throw new RuntimeException("!!" +f.getAbsolutePath());
    return new FileOutputStream(f);
}

public static void split(File in, File[] out, final Map<String, Integer> indices, List<String> types, boolean dosample) throws Exception{
    FileOutputStream[] dest = new FileOutputStream[out.length];
    CheckedOutputStream[] checksum = new CheckedOutputStream[out.length];
    ZipOutputStream[] outS = new ZipOutputStream[out.length];
    OutputStreamWriter[] osw = new OutputStreamWriter[out.length];
    List[] resu = new List[out.length];
    ZipFile zf = new ZipFile(in);
    if(in.length()==0) return;
    List<String> snps = Compressor.getIndiv(zf, "SNPS");
    List<String> indiv = Compressor.getIndiv(zf, "Samples");
    List<String> names = Compressor.getIndiv(zf, "Name");
    String head = in.getName();
    Set<Integer> poss_vals = new HashSet<Integer>();
        for(int j=0; j<indiv.size(); j++){
            poss_vals.add(indices.get(indiv.get(j).split("#")[0]));
        }
    head = head.substring(0, head.indexOf(".zip"));
    for(int i=0; i<out.length; i++){
        if(poss_vals.contains(i)){
            if(!out[i].exists()) out[i].mkdir();
            File outF = new File(out[i], head+".zip");
            dest[i] = getOS(outF);
            checksum[i] = new CheckedOutputStream(dest[i], new Adler32());
            outS[i] = new 
             ZipOutputStream(new 
               BufferedOutputStream(checksum[i]));
           osw[i] = new OutputStreamWriter(outS[i]);
           outS[i].setMethod(ZipOutputStream.DEFLATED);
           writeEntry("Name", names, outS[i], osw[i]);
           writeEntry("SNPS", snps, outS[i], osw[i]);
           resu[i] = new ArrayList();
        }
    }
    List<Integer> alias = new ArrayList<Integer>(indiv.size());
    for(int i=0; i<indiv.size(); i++){
        String ind = indiv.get(i);
        int index = indices.get(ind.split("#")[0]);
        System.err.println(ind+" "+index+" "+types.get(index));
        alias.add(index);
    }
    
    
    for(Enumeration en = zf.entries(); en.hasMoreElements();){
        ZipEntry ent = (ZipEntry) en.nextElement();
        if(ent.getName().startsWith("Name") || ent.getName().startsWith("SNPS")) continue;
        List<String> res = getIndiv(zf, ent.getName());
        Compressor.subList(res, alias, resu);
        for(int i=0; i<out.length; i++){
            if(outS[i]==null || resu[i].size()==0) continue;
         //   System.err.println(out[i].getName()+" "+types.get(i)+" "+resu[i].size());
            writeEntry(ent.getName(), resu[i],outS[i], osw[i]);
            resu[i].clear();
        }
    }
    String sample_header = null;
    if(dosample){
        List<String> samples = new ArrayList<String>();
    
    BufferedReader br =Utils.getBufferedReader(new File(in.getParentFile(), "Samples.txt"));
    
   sample_header =br.readLine(); 
    Utils.readPosInfo(br,new int[] { 0}, false, new List[] {samples}, new Class[] {String.class}, "________________");
    Compressor.subList(samples, alias, resu);
    br.close();
    }
    for(int i=0; i<out.length; i++){
        if(outS[i]!=null) outS[i].close();
        zf.close();
        if(sample_header!=null && outS[i]!=null){
            PrintWriter sample_pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(out[i], "Samples.txt"))));
            sample_pw.println(sample_header);
            for(int j=0; j<resu[i].size(); j++){
                sample_pw.println(resu[i].get(j));
            }
            sample_pw.close();
        }
    }
    System.err.println("done "+in);
    
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#readCaseControl(java.io.File)
 * 
 *  
    
  desir are the controls!
 */
public static Map<String, Integer> readCaseControl(File cca, 
        List<Integer> toExcl,
        String val, List<String>types) throws Exception{
    BufferedReader br = new BufferedReader(new FileReader(cca));
    HashMap<String, Integer> cc = new HashMap<String, Integer>();
    if(br!=null){
        String st = br.readLine();
        while((st = br.readLine())!=null){
            String[] str = st.split("\\s+");
            String id = str[1];
            if(toExcl.contains(Integer.parseInt(id))) continue;
            int col=0;
            inner: for(; col<str.length; col++){
                if(str[col].startsWith(val)) break inner;
            }
            String type = str[col];//.split("-DNA")[0];
            int index = types.indexOf(type);
            if(index<0){
                index = types.size();
                types.add(type);
            }
             cc.put(id, index);
        }
        br.close();
        //if(cc.size()!=this.dataL.size()) throw new RuntimeException("!!!");
        
    }
    return cc;

}

public  void copy(String[] exception, int[][] alias) throws Exception{
    List<String> excs = Arrays.asList(exception);
    File f1 = this.in;
    ZipFile zf = new ZipFile(f1);
   
    String chr = f1.getName().substring(0,f1.getName().indexOf(".zip"));
    File tmp = new File(f1.getParentFile(), chr+System.currentTimeMillis()+".zip");
    
    FileOutputStream dest = getOS(tmp);
    checksum = new   CheckedOutputStream(dest, new Adler32());
    outS = new 
     ZipOutputStream(new 
       BufferedOutputStream(checksum));
   osw = new OutputStreamWriter(outS);
   outS.setMethod(ZipOutputStream.DEFLATED);
    for(Enumeration en = zf.entries(); en.hasMoreElements();){
        ZipEntry ent = (ZipEntry) en.nextElement();
        List<String> res = this.getIndiv(zf, ent.getName());
           ZipEntry headings = new ZipEntry(ent.getName());
           outS.putNextEntry(headings);
           int index = excs.indexOf(ent.getName());
           if(index>=0){
               for(int i=0; i<res.size(); i++){
                   osw.write(applyAlias(res.get(i), alias[index]));
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
    tmp.renameTo(f1);
}

    private String applyAlias(String string, int[] is) {
        String[]str = string.split("\t");
        StringBuffer sb = new StringBuffer();
        sb.append(str[is[0]]);
        for(int i=1; i<str.length; i++){
           sb.append("\t");
            sb.append(str[is[1]]);
           
        }
        return "B allele\tLog R";
//"b allele\tlogR\tGenotype\tBases";
     //   return sb.toString();
    }
    
    public static void subList(List<String> res, List<Integer> alias, List[] resu){
        for(int i=0; i<res.size(); i++){
            resu[alias.get(i)].add(res.get(i));
        }
    }
    public static void writeEntry(String name, List<String> res, ZipOutputStream outS, OutputStreamWriter osw) throws Exception{
        ZipEntry headings = new ZipEntry(name);
        outS.putNextEntry(headings);
        for(int i=0; i<res.size(); i++){
            osw.write(res.get(i));
            osw.write("\n");
        }
        osw.flush();
        outS.closeEntry();
    }
    
    public static void writeEntry(String name, List<String> res, List<String>[] res1,  ZipOutputStream outS, OutputStreamWriter osw) throws Exception{
        ZipEntry headings = new ZipEntry(name);
        outS.putNextEntry(headings);
        for(int i=0; i<res.size(); i++){
        	String st = res.get(i);
        	int ind = st.lastIndexOf("\t",st.lastIndexOf("\t")-1);
            if(ind>=0) {
            	osw.write(st.substring(0,ind));
            
	            for(int j=0; j<res1.length; j++){
	            	osw.write("\t"+res1[j].get(i).substring(ind+1));
	            }
            }
            else{
            	osw.write(st);
            }
            osw.write("\n");
        }
        osw.flush();
        outS.closeEntry();
    }
    
    public static void writeEntry(String name, List<String> res, File[] names, int loess_ind, ZipOutputStream outS, OutputStreamWriter osw) throws Exception{
        ZipEntry headings = new ZipEntry(name);
        outS.putNextEntry(headings);
        for(int i=0; i<res.size(); i++){
        	String st = res.get(i);
        
            osw.write(st.substring(0,st.lastIndexOf("\t")));
            for(int j=0; j<names.length; j++){
            	osw.write("\t"+names[i].getName()+"_loess");
            }
            osw.write("\n");
        }
        osw.flush();
        outS.closeEntry();
    }
    public void merge(File dir1, File dir2, File out, String chr) throws Exception{
         int[] first = new int[] {0};
         int[] second = new int[] {1};
        ZipFile zf1 = new ZipFile(new File(dir1, chr+".zip"));
        ZipFile zf2 = new ZipFile(new File(dir2, chr+".zip"));
        String[] names1 = this.getIndiv(zf1, "Name").get(0).split("\t");
        String[] names2 = this.getIndiv(zf2, "Name").get(0).split("\t");
        List<String> samples1 = this.getIndiv(zf1, "Samples");
        List<String> samples2 = this.getIndiv(zf2, "Samples");
        List<Integer> conv = getMap(samples1, samples2);
        for(int i=0; i<samples1.size(); i++){
            if(!samples2.get(conv.get(i)).equals(samples1.get(i))) throw new RuntimeException("!!");

        }out.mkdir();
        FileOutputStream dest = getOS(new File(out, chr+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
         ZipOutputStream(new 
           BufferedOutputStream(checksum));
       osw = new OutputStreamWriter(outS);
       outS.setMethod(ZipOutputStream.DEFLATED);
       //if(header){
           ZipEntry headings = new ZipEntry("Name");
           outS.putNextEntry(headings);
           writeLine1(osw, resultColumnHeader);
           osw.flush();
       // for(int i=0; i<samples1.size(); i++){
        //    if(!samples1.get(i).equals(samples2.get(i))) throw new RuntimeException("!!");
       // }
        for(Enumeration en = zf1.entries(); en.hasMoreElements();){
            
            ZipEntry ent1 =(ZipEntry) en.nextElement();
           String ent_name = ent1.getName();
           if(ent_name.startsWith("Name") || ent_name.startsWith("SNP") || ent_name.startsWith("Sample")) continue;
           List<String> geno1 = this.getIndiv(zf1, ent_name);
          
           List<String> geno2 = null;
           if( zf2.getEntry(ent_name)!=null) geno2 = this.getIndiv(zf2, ent_name);
          
           ZipEntry ze = new ZipEntry(ent_name);
           outS.putNextEntry(ze);
           for(int i=0; i<geno1.size(); i++){
              
               this.writeLine(osw,  geno1.get(i).split("\t"), first, "\t");
              if(geno2!=null)
               this.writeLine(osw, geno2.get(conv.get(i)).split("\t"), second, "\n");
              else{
                  osw.write("0");
                  for(int ik=1; ik<second.length; ik++){
                      osw.write("\t0");
                  }
                  osw.write("\n");
              }
              
           }
           osw.flush();
           outS.closeEntry();
         
        }
         ZipEntry ze = new ZipEntry("Samples");
         outS.putNextEntry(ze);
         for(int i=0; i<samples1.size(); i++){
             osw.write(samples1.get(i)+"\n");
         }
         osw.flush();
        outS.closeEntry();
       
        ZipEntry ze1 = new ZipEntry("SNPS");
        outS.putNextEntry(ze1);
        List<String> snps = this.getIndiv(zf1, "SNPS");
        for(int i=0; i<snps.size(); i++){
            osw.write(snps.get(i)+"\n");
        }
        osw.flush();
       outS.closeEntry();
      
       osw.close();
    }
    
    /** rather than merge, simply joins together multiple samples on same chip */
    public void join(File[] dir1,  File out, String chr) throws Exception{
        int[] first = new int[] {0};
        int[] second = new int[] {1};
       ZipFile[] zf1 = new ZipFile[dir1.length];//(new File(dir1, chr+".zip"));
     //  String[][] names1 = new String[dir1.length][];
       List<String>[] samples = new List[zf1.length];
       for(int i=0; i<zf1.length; i++){
           zf1[i] = new ZipFile(new File(dir1[i], chr+".zip"));
       //   names1[i] = this.getIndiv(zf1[i], "Name").get(0).split("\t");
          samples[i] = this.getIndiv(zf1[i], "Samples");
       }
       List<String> names = this.getIndiv(zf1[0], "Name");
       out.mkdir();
       FileOutputStream dest = getOS(new File(out, chr+".zip"));
       checksum = new   CheckedOutputStream(dest, new Adler32());
       outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
      osw = new OutputStreamWriter(outS);
      outS.setMethod(ZipOutputStream.DEFLATED);
      //if(header){
       //   ZipEntry headings = new ZipEntry("Name");
       //   outS.putNextEntry(headings);
          this.writeEntry("Name", names, this.outS, this.osw);
       //   writeLine1(osw, resultColumnHeader);
       //   osw.flush();
      // for(int i=0; i<samples1.size(); i++){
       //    if(!samples1.get(i).equals(samples2.get(i))) throw new RuntimeException("!!");
      // }
       for(Enumeration en = zf1[0].entries(); en.hasMoreElements();){
           
           ZipEntry ent1 =(ZipEntry) en.nextElement();
          String ent_name = ent1.getName();
          if(ent_name.startsWith("Name") || ent_name.startsWith("SNP") ) continue;
          List<String>[] geno1 = new List[dir1.length];//
          for(int i=0; i<zf1.length; i++){
             geno1[i] =  this.getIndiv(zf1[i], ent_name);
          }
          ZipEntry ze = new ZipEntry(ent_name);
          outS.putNextEntry(ze);
          for(int k=0; k<geno1.length; k++){
              for(int i=0; i<geno1[k].size(); i++){
                  osw.write(geno1[k].get(i)+"\n");
//                  this.writeLine(osw,  geno1[k].get(i).split("\t"), first, "\n");
              }
          }
          osw.flush();
          outS.closeEntry();
        
       }
      
      
       ZipEntry ze1 = new ZipEntry("SNPS");
       outS.putNextEntry(ze1);
       List<String> snps = this.getIndiv(zf1[0], "SNPS");
       for(int i=0; i<snps.size(); i++){
           osw.write(snps.get(i)+"\n");
       }
       osw.flush();
      outS.closeEntry();
     
      osw.close();
   }
  
    private String[] changeBase(String[] geno1St) {
      for(int i=0; i<geno1St.length; i++){
          geno1St[i] = ""+this.convertToE(2,Double.parseDouble(geno1St[i]));
      }
        return geno1St;
    }
    public  void compress () {
        out = new File(out.getParentFile(), out.getName().substring(0, out.getName().indexOf(".txt.zip")));
        out.mkdir();
        int chr_index = -1;
        header = false;
    //    if(out.exists())return ;
       try {
       
          ids = new ArrayList<String>();
          BufferedReader origin = new BufferedReader(new InputStreamReader(
                 // new GZIPInputStream(
                          new FileInputStream(in)
                   //       )
                  ));
          if(header) origin.readLine();
        
          String st = "";
          String comp = null;
          
          int index =   snpColumns[0];
          
        String chr = "";
        int count = 0;
        String prev_id = "";
        
         for(int ik=0; (st = origin.readLine())!=null; ik++){
              String[] str = st.split(split);
       
              if(prev_id.equals(str[id_index])) {
              //    System.err.println(st);
                 System.err.println("warning should not equal "+prev_id+" "+str[id_index]);
                  ik--;
                 continue;
              }
            
              prev_id = str[id_index];
              if((chr.length()==0 && chr_index==-1) || (chr_index>=0 && !str[chr_index].equals(chr))){
                 
                  if(chr.length()>0) {
                    //  System.err.println("finishing chromosome "+chr);
                    //  System.err.println("no entries "+count);
                      this.finishChrom();
                    //  if(true) System.exit(0);
                  }
                  chr= chr_index < 0 ? "res" : str[chr_index];
               System.err.println("doing chr "+chr);
                  FileOutputStream dest = getOS(new File(out, chr+".zip"));
                   checksum = new   CheckedOutputStream(dest, new Adler32());
                   outS = new 
                    ZipOutputStream(new 
                      BufferedOutputStream(checksum));
                  osw = new OutputStreamWriter(outS);
                  outS.setMethod(ZipOutputStream.DEFLATED);
                  //if(header){
                      ZipEntry headings = new ZipEntry("Name");
                      outS.putNextEntry(headings);
                      writeLine1(osw, resultColumnHeader);
                      osw.flush();
                  //}
                  
              }
              
              if(comp==null || !comp.equals(str[index])){
                 osw.flush();
                 count++;
                  ZipEntry entry = new ZipEntry(proc(str[this.snpColumns[0]]));
                  //System.err.println(entry.getName());
                  outS.putNextEntry(entry);
                  comp = str[index];
                  ik=0;
                  List<String> det = getSubString(str);
                  details.add(det);
              }
            //  System.err.println(st+" "+ik);
              if(ik>=ids.size()) {
                  ids.add(str[id_index]);
              }
              else if(!ids.get(ik).equals(str[id_index])){
                  throw new RuntimeException("is "+str[id_index]+" should be "+ids.get(ik)+ " at "+ik+
                          "\n"+ids.subList(0, ik+1));
              }
              writeLine(osw, str);
          }
         System.err.println("finishing chromosome "+chr);
         this.finishChrom();
       } catch(Exception e) {
          e.printStackTrace();
       }
    }
    ZipOutputStream outS = null;
    OutputStreamWriter osw = null;
    CheckedOutputStream checksum;
    List<String> ids   = new ArrayList<String>();
    List<List<String>> details = new ArrayList<List<String>>();
    
    public void finishChrom() throws Exception{
        osw.flush();
        System.err.println("writing samples "+this.ids.size());
        ZipEntry samples = new ZipEntry("Samples");
        outS.putNextEntry(samples);
        writeLine(osw,ids);
        osw.flush();
        outS.closeEntry();
        ZipEntry snps = new ZipEntry("SNPS");
        outS.putNextEntry(snps);
       // Collections.sort(details,compa);
        for(Iterator<List<String>> it = details.iterator(); it.hasNext();){
           writeLine1(osw, it.next().toArray(new String[0]));
        }
        osw.flush();
        outS.closeEntry();
        outS.finish();
        osw.close();
        outS.close();
        System.out.println("checksum: "
          +checksum.getChecksum().getValue());
        details.clear();
    }
    
 
    public  void compressInfinium () {
      // if(out.exists())return ;
     //  else{
           this.infinium =null;// new File(out.getParentFile(), "inf_sorted.txt.gz");
           out = new File(out.getParentFile(), out.getName().substring(0, out.getName().indexOf(".txt.zip")));
           out.mkdir();
     //  }
       try {
           int chr_index = 1;
           header = true;
              resultColumnHeader = new String[] {
                                         "Genotype", 
                                          "B allele", 
                                          "Log R",
                                          "Top allelles"};
//                                          "Log RU"};
              int start = 2;
              
           int rep = resultColumnHeader.length;
        
          BufferedReader origin = new BufferedReader(new InputStreamReader(
                  //new GZIPInputStream(
                          new FileInputStream(in)
                    //      )
                          ));
          BufferedReader infiniumB = 
              infinium==null ? null:
              new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(this.infinium))));
          if(infinium!=null){
              infiniumB.readLine();
          }
          if(header){
              String[] st = origin.readLine().split(split);
                  //process( origin.readLine().split(split));
              for(int i= start; i<st.length; i+=rep){
                  ids.add(st[i].split("_")[0]);
              }
          }
         
          String st = "";
         
     
       String[] toPrint = new String[rep];
       String chr = "";
      int count =0;
       String[] inf =   infiniumB==null ? null : infiniumB.readLine().split(split);
         for(int ik=0; (st = origin.readLine())!=null; ik++){
            
              String[] str = process(st.split(split));
            //  System.err.println(str[0]);
              if( !str[chr_index].equals(chr)){
                  if(chr.length()>0){
                      
                      finishChrom();
                     /* ZipFile zf = new ZipFile(new File(out,1+".zip"));
                      ZipEntry ent = zf.getEntry("Samples");
                     List<String> samples =  this.getIndiv(zf, "Samples");
                     System.err.println(samples);
                      System.exit(0);*/
                      
                  }
                  chr = str[chr_index];
                  System.err.println("doing "+chr);
                  if(inf!=null){while(!chr.equals(inf[2])){
                      inf =   infiniumB.readLine().split(split);
                  }
                  }
                  FileOutputStream dest = getOS(new File(out, chr+".zip"));
                  checksum = new   CheckedOutputStream(dest, new Adler32());
                  outS = new 
                    ZipOutputStream(new 
                      BufferedOutputStream(checksum));
                  osw = new OutputStreamWriter(outS);
                  outS.setMethod(ZipOutputStream.DEFLATED);
                  if(header){
                      ZipEntry headings = new ZipEntry("Name");
                      outS.putNextEntry(headings);
                      writeLine1(osw, resultColumnHeader);
                      osw.flush();
                      outS.closeEntry();
                  }
              }
              count++;
              List<String> det = getSubString(str);
              if(infiniumB!=null){
                 
                //  System.err.println(det.get(0)+" "+inf[1]);
                  int pos = Integer.parseInt(det.get(2));
                  while(!det.get(0).equals(inf[1]) || inf[0].indexOf("v1")<0){
                      if(Integer.parseInt(inf[3])>pos) throw new RuntimeException("!!");
                      inf = infiniumB.readLine().split(split);
                     
                  }
//                  if(inf[0].indexOf("v1")<0) System.err.println("!! "+Arrays.asList(inf)+" vs "+det);
                  det.add(inf[4]); det.add(inf[5]);
                  String in = infiniumB.readLine();
                  if(in!=null) inf = in.split(split);
                  
              }
              
               details.add(det);
               ZipEntry entry = new ZipEntry(proc(str[snpColumns[0]]));
               
               outS.putNextEntry(entry);
              for(int j = start; j<str.length; j+=rep){
                  for(int k=0; k<rep; k++){
                      toPrint[k] = str[j+k];
                  }
                  writeLine1(osw, toPrint);
              }
              
              osw.flush();
              outS.closeEntry();
             
          }
          finishChrom();
         // if(true) System.exit(0);
       } catch(Exception e) {
          e.printStackTrace();
       }
    }
   
    
    public String[] process(String[] str){
        if(true) return str;
        String[] st1 = new String[str.length+1];
        System.arraycopy(str, 0, st1, 1    , str.length);
        if(str[0].equals("Name")){
            st1[0] ="chr";
        }
        else{
            st1[0] = str[0].split("_")[0];
            st1[1] = str[0].split(":")[1];
        }
        return st1;
    }
    public List <String> getSubString(String[] str){
       // String[] str = process(str1);
        String[] res = new String[this.snpColumns.length];
        for(int i=0; i<snpColumns.length; i++){
            res[i] = proc(str[snpColumns[i]]);
        }
        return new ArrayList<String>(Arrays.asList(res));
    }
   private String proc(String string) {
       return string.replace(":", "$");
    }
public static BufferedReader readZip(final ZipFile f, final List<String> l ){
       Enumeration<InputStream> inputStreams = new Enumeration<InputStream>(){
          int i=0;
          public boolean hasMoreElements() {
              return i < l.size();
          
          }

          public InputStream nextElement() {
              try{
           //       System.err.println(l.get(i));
                InputStream nxt = f.getInputStream(f.getEntry(l.get(i)));
              i++;
                return nxt;
              }catch(Exception exc){
                  exc.printStackTrace();
                  return null;
              }
          }
          
      };
     
     return new BufferedReader(new InputStreamReader(new SequenceInputStream(inputStreams)));
   }
   
   public static List<String> getIndiv(ZipFile f, String entryName) throws Exception{
       if(f.getEntry(entryName)==null) return null;
//    	   throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           indiv.add(st.trim());
       }
      nxt.close();
       return indiv;
   }
   
   
   
   public static List<String>[] getIndiv(ZipFile[] f, String entryName) throws Exception{
      List[] res = new List[f.length];
      for(int i=0; i<res.length; i++){
    	  res[i] = getIndiv(f[i], entryName);
      }
      return res;
   }
   
   public static List<String> getIndiv(ZipFile f, String entryName, int id) throws Exception{
       if(f.getEntry(entryName)==null){
    	 
    	   return null;
       }
    	   //throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
    	   String[] str = st.trim().split("\t");
    	   if(id <str.length){
           indiv.add(str[id]);
    	   }
    	   else{
    		   indiv.add(null);
    	   }
       }
      nxt.close();
       return indiv;
   }
   
   public static List<String> getIndiv(ZipFile f, String entryName, int[] alias) throws Exception{
       if(f.getEntry(entryName)==null){
    	 
    	   return null;
       }
    	   //throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
     
       while((st = nxt.readLine())!=null){
    	   String[] str = st.trim().split("\t");
    	   StringBuffer sb = new StringBuffer(); 
    	   for(int k=0; k<alias.length; k++){
    		 if(k>0) sb.append("\t");
    		 sb.append(alias[k]<0 ? "NA" : str[alias[k]]);
    	   }
    	   indiv.add(sb.toString());
    	   
       }
      nxt.close();
       return indiv;
   }
   
   public static List<String> getIndiv(File f, int id) throws Exception{
   
    	   //throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
       BufferedReader  nxt = Utils.getBufferedReader(f);
           
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           indiv.add(st.trim().split("\t")[id]);
       }
      nxt.close();
       return indiv;
   }
   
   public static List<String> getIndiv(File f) throws Exception{
	   
	   //throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
   BufferedReader  nxt = Utils.getBufferedReader(f);
       
   List<String> indiv = new ArrayList<String>();
   if(nxt==null) return indiv;
   String st = "";
   while((st = nxt.readLine())!=null){
       indiv.add(st.trim());
   }
  nxt.close();
   return indiv;
}
   public static boolean getIndiv(ZipFile f, String entryName, int id,
			String[] res) throws Exception{
	   if(f.getEntry(entryName)==null){
		   
		   return false;
		  // throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
	   }
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       for(int i=0; (st = nxt.readLine())!=null; i++){
           res[i]=(st.trim().split("\t")[id]);
       }
      nxt.close();
      return true;
	}
   
   public static boolean getIndiv1(ZipFile f, String entryName, int id,
			String[] res, String chr) throws Exception{
	   if(f.getEntry(entryName)==null){
		   entryName = chr+"/"+entryName+".txt";
		   if(f.getEntry(entryName)==null)
		   return false;
		  // throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
	   }
      BufferedReader  nxt = 
          new BufferedReader(new InputStreamReader(
          f.getInputStream(f.getEntry(entryName))));
      //System.err.println(entryName);
   //   List<String> indiv = new ArrayList<String>();
      String st = "";
      for(int i=0; (st = nxt.readLine())!=null; i++){
    	  String[] str_ = st.trim().split("\t");
    	  if(id>=str_.length){
    		  res[i] = "NaN";
    	  }
    	  else res[i]=(str_[id]);
      }
     nxt.close();
     return true;
	}
   public static List<String> getIndiv(ZipFile f, String entryName, String toadd) throws Exception{
       if(f.getEntry(entryName)==null) throw new RuntimeException("does not contain entry in \n"+f.getName()+" \n "+entryName);
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(entryName))));
       List<String> indiv = new ArrayList<String>();
       String st = "";
       while((st = nxt.readLine())!=null){
           indiv.add(st+"\t"+toadd);
       }
      nxt.close();
       return indiv;
   }
   
   public static void read(ZipFile f, String string, List<String> indiv,
           int i) throws Exception {
       BufferedReader  nxt = 
           new BufferedReader(new InputStreamReader(
           f.getInputStream(f.getEntry(string))));
       String st = "";
       while((st = nxt.readLine())!=null){
           indiv.add(st.split("\t")[i]);
       }
       nxt.close();
       
   }
   private static String conc(String[] str, int from){
       StringBuffer sb = new StringBuffer(str[from]);
       for(int i=from+1; i<str.length; i++){
           sb.append(str[i]);
       }
       return sb.toString();
   }
   public static BufferedReader readZip(final ZipFile f, String ent) throws Exception{
	 return new BufferedReader(new InputStreamReader( f.getInputStream( f.getEntry(ent))));
   }
   
    public static BufferedReader readZip(final ZipFile f, final int from, final int to
            ,List<Integer> loc, List<String> chr,  List<String> snpid) throws Exception{
       List<String> l = new ArrayList<String>();
      Enumeration entries = f.entries();
           while(entries.hasMoreElements()){
               ZipEntry entry = (ZipEntry) entries.nextElement();
               String name = entry.getName();
               if(name.equals("Name") || name.equals("Sample")) continue;
               try{
                   String[] names = name.split("_");
               int no = Integer.parseInt(names[0]);
            //   System.err.println(no);
              if(no>=from){
                  if(no <= to){
                      l.add(name);
                      loc.add(no);
                      chr.add(names[1]);
                      snpid.add(conc(names,2));
                   //   System.err.println(no);
                  }
                  else break;
              }
               }catch(NumberFormatException exc){
                   continue;
               }
           }
          // read(f, l.get(0), indiv, 0);
           return readZip(f, l);
       }
    
    
    public static List<String> readZipFrom(File in, String name) throws Exception{
        final  ZipFile f = new ZipFile(in);
        BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(name))));
        String st;;// = br.readLine();
        List<String> l = new ArrayList<String>();
        while((st = br.readLine())!=null){
            l.add(st);
        }
        return l;
    }
    public static BufferedReader readZipFrom(ZipFile f, final int from, final int to
            ,List<Integer> loc, List<String> chr,  List<String> snpid
    ) throws Exception{
     //  final  ZipFile f = new ZipFile(in);
        final Enumeration entries =  f.entries();
        List<String> l = new ArrayList<String>();
        for(int i=0; i<from; i++){
            entries.nextElement();
        }
        for(int i=0; i<to; i++){
            String name =((ZipEntry)entries.nextElement()).getName(); 
            l.add(name);
            String[] names = name.split("_");
            loc.get(Integer.parseInt(names[0]));
            chr.add(names[1]);
            snpid.add(conc(names,2));
        }
     //   read(f, l.get(0), indiv, 0);
        return readZip(f, l);
        
        
     }
   

    public static void unzip(File in) {
        try {
           final int BUFFER = 2048;
           BufferedOutputStream dest = null;
           FileInputStream fis = new 
         FileInputStream(in);
           CheckedInputStream checksum = new 
             CheckedInputStream(fis, new Adler32());
           ZipInputStream zis = new 
             ZipInputStream(new 
               BufferedInputStream(checksum));
           ZipEntry entry;
           while((entry = zis.getNextEntry()) != null) {
              System.out.println("Extracting: " +entry);
              int count;
              byte data[] = new byte[BUFFER];
              // write the files to the disk
              FileOutputStream fos = 
                  getOS(new File(entry.getName()));
              dest = new BufferedOutputStream(fos, 
                BUFFER);
              while ((count = zis.read(data, 0, 
                BUFFER)) != -1) {
                 dest.write(data, 0, count);
              }
              dest.flush();
              dest.close();
           }
           zis.close();
           System.out.println("Checksum:  "+checksum.getChecksum().getValue());
            
        } catch(Exception e) {
           e.printStackTrace();
        }
     }
    
    public String modify(String st){
        return st.replace(':',',');
    }
      /*  if(i==0){
          
        }
        else return st;
        }*/
    
    
   
    Compressor(File f) throws Exception{
        this.in = f;
        String o = f.getName();
        int ind = o.indexOf(".gz");
        if(ind>=0){
            o = o.substring(0,ind);
        }
         out = new File(f.getParentFile(), o+".zip");
         int ind1 = o.indexOf("_data.txt");
         if(ind1>=0){
             String o1 = o.substring(0,ind1);
             this.infinium = new File(f.getParentFile(), o1+"_infinium.txt.gz");
         }
         else infinium = null;
    
    }
    
    public static void readZip(final ZipFile f, String st, String[] res, int col, String spl){
        try{
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(st))));
          String str = "";
          int i=0;
          for(; (str = br.readLine())!=null;i++){
              String[] stri = str.split(spl);
              res[i]=stri[col];
          }
          if(i!=res.length) throw new RuntimeException("prob with lengths");
          br.close();
        }catch(Exception exc){
            exc.printStackTrace();
            //return null;
        }
    }
    
public static List<String> readSamples(ZipFile f) throws Exception {
    List<String> l = new ArrayList<String>();
    ZipEntry ent = f.getEntry("Samples");
    if(ent==null) return null;
    BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(ent)));
    String st = "";
    while((st = br.readLine())!=null){
        l.add(st.split("\\t")[0].trim());
    }
    return l;
}

public static List<String> read(ZipFile f, String head) throws Exception {
    List<String> l = new ArrayList<String>();
    ZipEntry ent = f.getEntry(head);
    if(ent==null) return null;
    BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(ent)));
    String st = "";
    while((st = br.readLine())!=null){
        l.add(st);
    }
    br.close();
    return l;
}
public static List<String> read(ZipFile f, ZipEntry ent) throws Exception {
    List<String> l = new ArrayList<String>();
   // ZipEntry ent = f.getEntry(head);
    if(ent==null) return null;
    BufferedReader br = new BufferedReader(new InputStreamReader(f.getInputStream(ent)));
    String st = "";
    while((st = br.readLine())!=null){
        l.add(st);
    }
    br.close();
    return l;
}
    public static List<Integer> readPositions(File file) throws Exception {
      ZipFile f = new ZipFile(file);
      List<Integer> l = new ArrayList<Integer>();
      for(Enumeration en = f.entries(); en.hasMoreElements();){
          ZipEntry ent = (ZipEntry) en.nextElement();
          if(ent.getName().equals("Samples") || ent.getName().equals("Name")) continue;
          String nm = ent.getName().split("_")[0];
          l.add(Integer.parseInt(nm));
      }
      return l;
    }
	
}
