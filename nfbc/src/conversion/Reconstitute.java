package conversion;


import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class Reconstitute {
    
    public static final Options OPTIONS  = new Options(){
        {
                    this.addOption( OptionBuilder.withLongOpt( "r_home" ).withDescription( "r_home").withValueSeparator( ':' ).hasArgs().create());
                    
        }
    };    

    File f;
  
    static Logger logger = Logger.getAnonymousLogger();
    static boolean hasgc = false; 
    static int chr_index = 0;
    static int pos_index = 1;
    static int id_index = 3;
 
      
    
  
      
      
        
    
 
public static void main( File dir, File outF, String todo, String index){
      String end = ".txt";
      try{
       // 
          File outDir = new File(outF, dir.getName());
          outDir.mkdir();
         /* File sampleFile = new File(dir, "Samples.txt");
          if(sampleFile.exists() && sampleFile.length()>0){
              Utils.copy(sampleFile, new File(outDir,sampleFile.getName() ));
          }*/
          File f = new File(dir, todo+".zip");
           
             //outer: for(int i=0; i<f.length; i++){
                
                 try{
                     File outf = new File(outDir, todo);
                     /*  if(outf.exists() && outf.length()>0) {
                    	   System.err.println("file exists "+outf.getAbsolutePath());
                    	   System.err.println("skipping");
                    	   return ;
                       }*/
                       try{
                    	 
                         
                           File summ = new File(f.getParentFile(),f.getName().substring(0,f.getName().indexOf(".zip"))+ "/summary"+index+".out");
                           if(!summ.exists() || summ.length()==0){
                        	   throw new RuntimeException("cannot find summary file "+summ.getAbsolutePath()+"  this is probably due to " +
                        	   		"not having limma accesible by R.  It needs to be installed in the R distribution, or alternatively you can specify" +
                        	   		"--rlib param which is an extra r library which will be searched");
                        	  // return;
                           }
                           Reconstitute cw = new Reconstitute(f, new String[] {"Log R"}, new File(dir, ConvertToWide.build+".txt"), end, outDir, index);
                     cw.run();
                       }catch(Exception exc){
                    	   exc.printStackTrace();
                       }
                 
                 }catch(Exception exc){
                     exc.printStackTrace();
                 }
            // }
            
             //System.err.println(err.getBuffer().toString());
             
        // }
         
             }catch(Exception exc){
          exc.printStackTrace();
      }
      
  }
  
 
    String[] samples;
    Map<String,List<String>> sampleInfo = new HashMap<String, List<String>>();
    List sampleInfoHeader = new ArrayList<String>();
    
    final String chr;
     List<String> rs_header;
     List<String> snp_header;
     List<String> sample_header;
    int[] indices;
    List<String> snps = new ArrayList<String>();
    List<String> loess = new ArrayList<String>();
    List<String> id = new ArrayList<String>();
    File outDir;
    FileOutputStream dest;
    CheckedOutputStream checksum;
    ZipOutputStream outS;
    OutputStreamWriter osw;
    File join;
    ZipFile zf;
    
    private void readSummaryFile(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        String[] header = br.readLine().split("\\t");
        List[] vals = new List[header.length];
        for(int i=0; i<header.length; i++){
        ////	int ind = header[i].lastIndexOf('_');
          //  header[i] = header[i].substring(0,ind);
            sampleInfo.put(header[i],vals[i] =  new ArrayList());
        }
        String st = "";
        while((st = br.readLine())!=null){
            if(st.startsWith("overall")) continue;
            else{
                sampleInfoHeader.add(st);
                String[] res = br.readLine().trim().split("\\s+");
                
                for(int i=0; i<header.length; i++){
                	if(res[i].length()==0){
                		System.err.println("h");
                	}
                    vals[i].add(res[i]);
                }
            }
        }
    }
    
    Reconstitute(File f_in, String[] feat, File buildFile, String ending, File outDir, String index) throws Exception{
        this.f = f_in;
      
        this.chr = f.getName().substring(0,f.getName().indexOf(".zip"));
        File SNPSF = new File(f.getParentFile(), "SNPS_"+chr);
        join  =  new File(f.getParentFile(),chr);
        File[] summ = join.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
			return pathname.getName().startsWith("summ");
			}
        	
        });
        if(summ.length==0){
        	throw new RuntimeException("!!");
        }
                this.zf = new ZipFile(f);
              
          this.outDir =outDir;
          outDir.mkdir();
          System.err.println("join - outdir is "+outDir+" "+chr+index+".zip");
          dest = Compressor.getOS(new File(outDir, chr+index+".zip"));
          checksum = new   CheckedOutputStream(dest, new Adler32());
          outS = new 
          ZipOutputStream(new 
            BufferedOutputStream(checksum));
         osw = new OutputStreamWriter(outS);
         outS.setMethod(ZipOutputStream.DEFLATED);
         File loessF = new File(join, "loess"+index);
      for(int i=0; i<summ.length; i++){
    	  readSummaryFile(summ[i]);
      }
        
             read(loessF, loess);
         readSNPS(SNPSF.exists()  ? Utils.getBufferedReader(SNPSF) : Utils.getBufferedReader(zf, "SNPS"));
   hasgc = SNPSF.exists();
      samples =  Compressor.readSamples(zf).toArray(new String[0]);
      List<String> nameInf = Compressor.read(zf, "Name");
      this.rs_header = Arrays.asList(nameInf.get(0).split("\\t")); 
      this.snp_header = new ArrayList<String>(Arrays.asList(nameInf.get(1).split("\\t"))); 
      this.sample_header = new ArrayList<String>();
      sample_header.add(Arrays.asList(nameInf.get(2).split("\\t")).get(0));
      if(sample_header.size()>1) sample_header.remove(sample_header.size()-1);
      int snp_head = this.snps.get(0).split("\t").length-2; //number of cols to keep in snp_header
      if(hasgc && snp_head < snp_header.size()) snp_header = snp_header.subList(0, snp_head);
      if(hasgc){
    	  snp_header.add("gc_content");
          snp_header.add("gc_correction");
      }
      snp_header.add("loess");
      List<Integer> snp_headerToRemove = new ArrayList<Integer>();
          sample_header.addAll(sampleInfoHeader);
    }
   
   
   
//final boolean replaceLoessColumn;
  
    public void run() throws Exception{
           ZipEntry headings = new ZipEntry("Name");
           outS.putNextEntry(headings);
           Compressor.writeLine1(osw, rs_header.toArray(new String[0]));
           Compressor.writeLine1(osw, snp_header.toArray(new String[0]));
           Compressor.writeLine1(osw, sample_header.toArray(new String[0]));
           osw.flush();
           outS.closeEntry();
           
           //now write all snp entries
           for(int ij=0; ij<this.id.size(); ij++){
              String ent_name = id.get(ij);
              List<String> geno1 =    Compressor.getIndiv(zf, ent_name);
              ZipEntry ze = new ZipEntry(ent_name);
              outS.putNextEntry(ze);
              for(int i=0; i<geno1.size(); i++){
                 if(sampleInfo.containsKey(samples[i])){
                	 osw.write(geno1.get(i)+"\n");
                 }
                 else{
                	 System.err.println("WARNING sample info did not contain "+samples[i]+" "+this.chr+" "+this.f.getName());
                 }
              }
              osw.flush();
              outS.closeEntry();
           }
           
           
           ZipEntry ze = new ZipEntry("Samples");
           outS.putNextEntry(ze);
           for(int i=0; i<samples.length; i++){
        	   List<String> vals = this.sampleInfo.get(samples[i]);
	            if(vals!=null){
	               osw.write(samples[i]);
	               for(int k2 =0; k2<vals.size(); k2++){
	                   osw.write("\t"+String.format("%5.3g", Double.parseDouble(vals.get(k2))));
	               }
	               osw.write("\n");
               }
           }
           osw.flush();
          outS.closeEntry();
          ZipEntry ze1 = new ZipEntry("SNPS");
          outS.putNextEntry(ze1);
          List<String> geno1 =   this.snps;// Compressor.getIndiv(zf, "SNPS");
          for(int i=0; i<snps.size(); i++){
           //   String str1 = snps.get(i);
              String str = geno1.get(i);
              if(!id.get(i).equals(str.split("\\s+")[3])) throw new RuntimeException("!!");
             /* if(replaceLoessColumn){
                  String[] val = str.split("\t");
                  val[0] = "chr"+chr;
                  val[val.length-1] = loess.get(i);
                Compressor.writeLine1(osw, val);
                  //osw.write(str+"\n");
              }
              else{*/
            	  String[] val = str.split("\t");
            	  val[0] = "chr"+chr;
            	  Compressor.writeLine2(osw, val);
                  //osw.write(str+"\t");
                  osw.write(loess.get(i)+"\n");
              //}
             
             
          }
          osw.flush();
         outS.closeEntry();
         outS.close();
        if(zf!=null) zf.close();
    }
    
   
    
    private void read(File file, List<String> res) throws Exception{
        BufferedReader br =  new BufferedReader(new FileReader(file));
        String st = "";
        res.clear();
        while((st = br.readLine())!=null){
            res.add(st);
        }
    }

    private void readSNPS(BufferedReader br) throws Exception{
     
        String st = "";
        String start = "chr"+chr;
        while((st = br.readLine())!=null){
            if(!st.startsWith(start)) continue;
            String[] str = st.split("\\s+");
            if(!str[chr_index].equals(start)) continue;
            snps.add(st);
            id.add(str[id_index]);
        }
        br.close();
    }
   
 
   
        
    
}
