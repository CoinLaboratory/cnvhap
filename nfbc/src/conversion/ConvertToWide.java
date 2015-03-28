package conversion;


import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

import conversion.Compressor;
import conversion.OptionBuild;
import conversion.Reconstitute;
import conversion.RunLoess;
import conversion.Utils;
import exec.Exec;

public class ConvertToWide {
    
    public static final Options OPTIONS  = new Options(){
        {
                    this.addOption( OptionBuilder.withLongOpt( "todo" ).withDescription( "todo").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "build" ).withDescription( "build").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "rlib" ).withDescription( "rlib").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "prefix" ).withDescription( "prefix").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "mode" ).withDescription( "mode").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "runR" ).withDescription( "runR").withValueSeparator( ':' ).hasArgs().create());
                    
        }
    };    

    File f;
    
    
    static void setOptions(boolean runR1, String build1){
    	runR = runR1;
    	build = build1;
    	snp_header = 
    	        hasgc ? 
    	        new String[] {"chr", build+"_start", build+"_end", "id", "gc", "loess"} :
    	            new String[] {"chr", build+"_start",build+"_end", "id", "loess"} ;
    }
    static double median_thresh = -0.1;
    static boolean runR = false;
    static String build = "build36";
    static Logger logger = Logger.getAnonymousLogger();
    static boolean hasgc = false; 
    static int chr_index = 0;
    static int pos_index = 1;
    static int id_index = 3;
    static String[] snp_header = 
        hasgc ? 
        new String[] {"chr", build+"_start", build+"_end", "id", "gc", "loess"} :
            new String[] {"chr", build+"_start",build+"_end", "id", "loess"} ; //"A", "B",
   static String[] sample_header = new String[] {"id", "medianCorrection"};
      static  boolean split = true; //if true, then we are splitting the zip into the outputs, otherwise we are joining
      
      static FileFilter ff = new FileFilter(){

          public boolean accept(File pathname) {
            return pathname.getName().startsWith("1M") && pathname.isDirectory();
          }
            
        };
        static FileFilter ff1 = new FileFilter(){

            public boolean accept(File pathname) {
              return  pathname.isDirectory();
            }
              
          };
      
      public static void main1(String[] args){
          try{
          File dir1 =  new File(System.getProperty("user.dir"));
          File[] dirs = dir1.listFiles(ff);
          for(int i=0; i<dirs.length; i++){
              File[] dirs1 = dirs[i].listFiles(ff1);
              for(int k=0; k<dirs1.length; k++){
                  System.err.println("deleting "+dirs1[k]);
                  
                  Utils.delete(dirs1[k]);
              }
          }
          }catch(Exception exc){
              exc.printStackTrace();
          }
      }
      public static void main(String[] args){
          try{
         Parser parser= new PosixParser();
         final CommandLine params = parser.parse(OPTIONS, args, false);
          File dir1 =  new File(System.getProperty("user.dir"));
          final String prefix = params.getOptionValue("prefix", "");
          build = params.getOptionValue("build", "build35");
         // rlib = params.getOptionValue("rlib", null);
          todo = params.getOptionValues("todo");
          runR = Boolean.parseBoolean(params.getOptionValue("runR", "false"));
          //r_home = r_home.substring(1,r_home.length()-1);
        //  Writer err = new BufferedWriter(new OutputStreamWriter(System.out));
          File[] dirs = dir1.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
              return pathname.getName().startsWith(prefix) && pathname.isDirectory();
            }
              
          });
          File outdir = new File(dir1, "modified");
          File pdf = new File(dir1, "pdf");
          pdf.mkdir();
          outdir.mkdir();
          List<String> mode = Arrays.asList(params.getOptionValues("mode"));
          for(int i=0; i<dirs.length; i++){
                  File dir = dirs[i];
              for(int k=0; k<todo.length; k++){
                  if(mode.contains("0")) main(true, dir, todo[k]);
                  if(mode.contains("1")){
                	  if(todo[k].equals("X")){
                		  runR(dir, todo[k], pdf, true, false, OptionBuild.varThresh);
                		
                	  }
                	  else{
                		  runR(dir, todo[k], pdf, false, false, OptionBuild.varThresh);
                	  }
                  }
                  if(mode.contains("2")){
                	  if(todo[k].equals("X") && false){
                		  Reconstitute.main(dir, outdir, todo[k], "_M");
                		  Reconstitute.main(dir, outdir, todo[k], "_F");
                	  }
                	  else{
                		  Reconstitute.main(dir, outdir, todo[k], "");
                	  }
                  }
                // Utils.delete(new File(dir, todo[k]));
              }
              
          }
          }catch(Exception exc){
              exc.printStackTrace();
          }
      }
   
     
      
      //static String rlib = null;
      
      public static void runR( File dir, String todo, File pdf_out, boolean stratify, boolean plotPdf, double var_thresh) throws Exception {
    	  RunLoess loess = new RunLoess(new File(dir, todo), var_thresh);
  			loess.run();
      }
      
  public static void runROld( File dir, String todo, File pdf_out, boolean stratify, boolean plotPdf, double var_thresh) throws Exception {
      File[] f = dir.listFiles(new FileFilter(){

        public boolean accept(File pathname) {
           return pathname.isDirectory() && !pathname.getName().startsWith("out");
        }
          
      });
      int i=0;
      Writer err = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "err"+todo+".txt"))));
   
     // for(; i<f.length; i++){
      if(!stratify){
          File f1 = new File(todo, "summary.out");
          if(f1.exists() && f1.length()>0){
        	  System.err.println("no need to run R for "+dir.getName());
              return;
          }
      }
      else{
    	  File f1 = new File(dir, todo+ "/summary_F.out");
    	  File f2 = new File(dir, todo+"/summary_M.out");
          if(f1.exists() && f1.length()>0 && f2.exists() && f2.length()>0){
        	  System.err.println("no need to run R for "+dir.getName());
              return;
          }
          
      }
          
     // }
     
      {
         
          ConvertToWide cw = new ConvertToWide();
         
         //  ClassLoader.getResourceAsStream("gcplot.R");
          File tmp = new File("gcplot_tmp+"+todo+"_"+System.currentTimeMillis()+".R");
       if(false) tmp.deleteOnExit();
        	 
          PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(tmp)));
          String st = "";
          File pdf = new File(pdf_out, dir.getName());
          pdf.mkdir();
          pw.println("pdfpath = \""+pdf.getAbsolutePath().replace('\\','/')+ "\"");
          pw.println("path = \""+dir.getAbsolutePath().replace('\\','/')+ "\"");
          pw.println("files1 = c(\""+todo+"\")");
          pw.println("doplot = "+(plotPdf ? "T": "F"));
          pw.println("var_thresh = "+var_thresh);
          if(!stratify){
        	   pw.println("med_thresh = c("+(-1000)+");");
               pw.println("lessThan = c(FALSE);");
               pw.println("index = c(\"\")");
          }
         
          else {
        	  pw.println("med_thresh = c("+(median_thresh)+","+(median_thresh)+")");
              pw.println("lessThan = c(FALSE, TRUE)");
              pw.println("index = c(\"_F\",\"_M\" )");
          }
          if(OptionBuild.rlib==null){
              pw.println("library(\"limma\")");
          }
          else{
              pw.println("library(\"limma\", lib.loc = \""+OptionBuild.rlib+"\")");
          }
          {
              BufferedReader br =  new BufferedReader(new InputStreamReader(cw.getClass().getResourceAsStream("gcplot")));
	          while((st = br.readLine())!=null){
	              pw.println(st);
	          }
	          pw.close();
	          br.close();
          }
          String[] exec = new String[] {"R", "--vanilla", "<", tmp.getAbsolutePath()};
          String ex_str = Arrays.asList(exec).toString().replaceAll(",", "");
          ex_str = ex_str.substring(1, ex_str.length()-1);
          System.err.println(ex_str);
       String[] exec1 = new String[] {"R", "--vanilla", "--slave"};
        BufferedReader br1 = new BufferedReader(new FileReader(tmp));
       if(runR) {
    	   try
           {            
    		   try
    	        {     
    			
    				  Exec.exec(exec1, br1,err,err);
    	        } catch (Throwable t)
    	          {
    	            t.printStackTrace();
    	          }

    		
           } catch (Throwable t)
             {
               t.printStackTrace();
             }

    	   
    	   
    	 //  
       }
        //  err.flush();
        
          
      }
      
      
        
    }
  
  
  
  static String[] todo = new String[0];
public static void main(boolean spl, File dir, String chr){
      String end = ".txt";
      ConvertToWide.split = spl;
      try{
       // 
          File f = new File(dir, chr+".zip");
        
          if(!f.exists()) return;
          System.err.println("opening "+f);
          ZipFile zf = new ZipFile(f);
     //   int  numSamples = (int) Math.ceil(Compressor.readSamples(zf).size()/(double)  ConvertToWide.lim);
             
                 try{
                    
                        
                         File outputdir = new File(f.getParentFile(), chr);
                         System.err.println(chr);
                       
                    Logger.global.info("running for "+f);
                    //for(Iterator<Map.Entry<String, List<String>>> entries = plate.entrySet().iterator(); entries.hasNext();){
                   
                     ConvertToWide cw = new ConvertToWide(f, new String[] {OptionBuild.feat}, new File(dir, build+".txt"), end, 0);
                     cw.run();
                    //}
                 
                 }catch(Exception exc){
                     exc.printStackTrace();
                 }
          
            
             //System.err.println(err.getBuffer().toString());
             
        // }
         
             }catch(Exception exc){
          exc.printStackTrace();
      }
      
  }
  
 
    File[][] files;
    PrintWriter[][] pw;
    int[][]cnt;
    BufferedReader[][] buff;
    String[] samples;
    double[][] medianCorrection;
     String chr;
    List<String> name;
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
    File SNPSF;
 public  static  int lim = 500;
    ConvertToWide(File f_in, String[] feat, File buildFile, String ending, int pos) throws Exception{
        if( this.split){
            this.f = f_in;
            this.zf = new ZipFile(f);
          
            this.chr = f_in.getName().substring(0, f.getName().indexOf(".zip"));
            SNPSF = new File(f.getParentFile(), "SNPS_"+(chr));
            join = new File(f.getParentFile(), chr);
            
           // this.checksum = f_in.getName();
        }
        else{
            File z_f = new File(f_in.getAbsolutePath()+".zip");
            if(z_f.exists() && z_f.length()>0){
                this.f = z_f;
                this.zf = new ZipFile(z_f);
                this.chr = f.getName().substring(0,f.getName().indexOf(".zip"));
                join  =  new File(f.getParentFile(),chr);
            }
            else{
                this.f = null;
                this.join = f_in;
                this.chr = f_in.getName();
            }
        }
      if(!join.exists()) join.mkdir();
      if(!split){
          outDir = new File ( f_in.getParentFile(), "outF");
          outDir.mkdir();
          dest = Compressor.getOS(new File(outDir, chr+".zip"));
          checksum = new   CheckedOutputStream(dest, new Adler32());
          outS = new 
          ZipOutputStream(new 
            BufferedOutputStream(checksum));
         osw = new OutputStreamWriter(outS);
         outS.setMethod(ZipOutputStream.DEFLATED);
         File loessF = new File(join, "loess");
         if(loessF.exists()){
             read(loessF, loess);
         }
         else loess = null;
      }
    BufferedReader br1;
    if(SNPSF.exists() && SNPSF.length()>0 ){
    	hasgc = true;
    	br1 = Utils.getBufferedReader(SNPSF);
    }
    	else{
    		br1  = Utils.getBufferedReader(zf, "SNPS");
    }
    if(br1==null){
    	br1 = Utils.getBufferedReader(zf, this.chr+"/SNPS");
    }
    if(br1==null){
    	br1 = Utils.getBufferedReader(buildFile);
    }
    if(br1==null){
    	throw new RuntimeException("no SNPS file for "+f_in);
    }
    if(zf!=null){
        samples =  Compressor.readSamples(zf).toArray(new String[0]);
        name = Arrays.asList(Compressor.read(zf, "Name").get(0).split("\\t+"));
        List<String>snp_ids = trim(Arrays.asList(Compressor.read(zf, "Name").get(1).split("\\s+")));
       id_index=   snp_ids.indexOf("id");
      }
      else{
          samples = getSampleList(f_in);
          name = Arrays.asList(new String[] {OptionBuild.feat});
      }
      readSNPS(br1);
      Logger.global.info("read "+this.snps.size());
    
    indices = new int[feat.length];
    this.medianCorrection = new double[indices.length][samples.length];
  
    for(int i=0; i<feat.length; i++){
        indices[i] = getIndex(name, feat[i]);
        if(indices[i]<0) indices[i] = getIndex(name, feat[i].replace(' ', '_'));
        if(indices[i]<0) {
        	
        	throw new RuntimeException(" cannot find "+feat[i]+" in "+name);
        }
    }

    pw = new PrintWriter[feat.length][samples.length];
    cnt = new int[feat.length][samples.length];
    buff = new BufferedReader[feat.length][samples.length];
    files = new File[feat.length][samples.length];
    
    if(!split) readMedianCorrection(new File(join, "summary.out"));
    for(int k=0; k<feat.length; k++){
        for(int i=pos*lim; i<Math.min(samples.length,(pos+1)*lim); i++){
        	System.err.println(i);
            files[k][i] = new File(join,(samples[i]));//+"_"+feat[k]).replaceAll("\\s+", "_")+ending);
            if(split){
              //  if(files[k][i].exists() && files[k][i].length()>0) throw new RuntimeException("!! "+files[k][i]+"\n"+Arrays.asList(files[k]));
                pw[k][i] = new PrintWriter(new BufferedWriter(new FileWriter(files[k][i])));
            }
            else{
                buff[k][i] = new BufferedReader(new FileReader(files[k][i]));
            }
        }
    }
    }
   
   
    private List<String> trim(List<String> asList) {
		for(int j=0; j<asList.size(); j++){
			asList.set(j, asList.get(j).trim());
		}
		return asList;
	}
	private int getIndex(List<String> name2, String string) {
    	int ind = name2.indexOf(string);
    	if(ind>=0) return ind;
		for(int i=0; i<name2.size(); i++){
			if(name2.get(i).indexOf(string)>=0) return i;
		}
		return -1;
	}
	public ConvertToWide() {
		// TODO Auto-generated constructor stub
	}
	private void readMedianCorrection(File file) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(file));
       List<String >ids = Arrays.asList(br.readLine().split("\t"));
       if(ids.size()!=this.samples.length){
           throw new RuntimeException( "!!");
       }
      /* for(int i=0; i<ids.size(); i++){
           ids.set(i, ids.get(i).split("_")[0]);
       }*/
        if(!br.readLine().equals("variance")) throw new RuntimeException("!!");
        String[] vari = br.readLine().split("\t");
        if(!br.readLine().equals("median")) throw new RuntimeException("!!");
        String[] median = br.readLine().split("\t");
        String[] overall_m = br.readLine().split("\\s+");
        if(!overall_m[0].equals("overall")) throw new RuntimeException("!!");
        double overallMedian = Double.parseDouble(overall_m[2]);
        double[] corr = new double[median.length];
        for(int i=0; i<corr.length; i++){
            corr[i] = overallMedian - Double.parseDouble(median[i]) ;
        }
        for(int i=0; i< this.samples.length; i++){
            this.medianCorrection[0][i] = corr[ids.indexOf(samples[i])];
        }
    }

    private String[] getSampleList(File f_in) {
       String[] res = f_in.list(new FilenameFilter(){

        public boolean accept(File arg0, String arg1) {
         return arg1.endsWith("R.txt");
        }
           
       });
       for(int i=0; i<res.length; i++){
           res[i] = res[i].split("_")[0];
       }
       return res;
    }

    public void workOutMedianCorrection() throws Exception{
            for(int i=0; i<this.id.size(); i++){
                ZipEntry ent = (ZipEntry) zf.getEntry(id.get(i));
                if(ent==null) throw new RuntimeException("!! did not contain "+id.get(i));
                BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
               for(int k1=0; k1<samples.length; k1++){
                    String[] in = br.readLine().split("\t");
                    for(int k2=0; k2<indices.length; k2++){
                        String[] in1 = this.buff[k2][k1].readLine().split("\t");
                        if(!snps.get(i).split("\\s+")[pos_index].equals(in1[0])){
                           logger.warning("loc mismatch "+snps.get(i)+"\n"+Arrays.asList(in1));
                           // throw new RuntimeException("!!");
                        }
                        double val_mod = Double.parseDouble(in1[1]);
                        double val_orig = Double.parseDouble(in[indices[k2]]);
                      //  System.err.println("orig "+ent.getName()+" "+id.get(i)+" "+i+" "+val_orig);
                      //  System.err.println(i+" "+k1+" "+samples[k1]+" "+id.get(i)+" "+val_mod+" "+val_orig);
                        if(i==0){
                            this.medianCorrection[k2][k1] = val_mod-val_orig;
                       //     System.err.println(samples[k1]+" "+medianCorrection[k2][k1]);
                        }
                        else{
                            double diff1 = val_mod-val_orig;
                            if(Math.abs(diff1 - this.medianCorrection[k2][k1])>0.00001) throw new RuntimeException(samples[k1]+"  "+diff1+" "+this.medianCorrection[k2][k1]);
                        }
                    }
       
                }
            
              }
    }

  
    public void write() throws Exception{
           ZipEntry headings = new ZipEntry("Name");
           outS.putNextEntry(headings);
           Compressor.writeLine1(osw, name.toArray(new String[0]));
           Compressor.writeLine1(osw, snp_header);
           Compressor.writeLine1(osw, sample_header);
           File association = new File(join, "association.txt");
           if(association.exists()){
               String[] ass = readAssociationFile(association);
               osw.write(ass[0]+"\n");
               osw.write(ass[1]+"\n");
           }
           osw.flush();
           outS.closeEntry();
           
           //now write all snp entries
           for(int ij=0; ij<this.id.size(); ij++){
              String ent_name = id.get(ij);
              List<String> geno1 = zf==null ?
                  readGeno(Integer.parseInt(this.snps.get(ij).split("\\s+")[1])):
                  Compressor.getIndiv(zf, ent_name);
              ZipEntry ze = new ZipEntry(ent_name);
              outS.putNextEntry(ze);
              for(int i=0; i<geno1.size(); i++){
                 osw.write(geno1.get(i)+"\n");
              }
              osw.flush();
              outS.closeEntry();
            
           }
           
           
           ZipEntry ze = new ZipEntry("Samples");
           outS.putNextEntry(ze);
           for(int i=0; i<samples.length; i++){
               for(int k2 =0; k2<medianCorrection.length; k2++){
               osw.write(samples[i]+"\t"+String.format("%5.3g", medianCorrection[k2][i])+"\n");
               }
           }
           osw.flush();
          outS.closeEntry();
          ZipEntry ze1 = new ZipEntry("SNPS");
          outS.putNextEntry(ze1);
          for(int i=0; i<snps.size(); i++){
              osw.write(snps.get(i)+"\t");
              if(loess!=null) osw.write(loess.get(i)+"\n");
              else osw.write("\n");
          }
          osw.flush();
         outS.closeEntry();
         outS.close();
        if(zf!=null) zf.close();
    }
    
    private List<String> readGeno(Integer pos) throws Exception {
       List<String> res = new ArrayList<String>();
       for(int k=0; k<this.buff[0].length; k++){
           String[] st = buff[0][k].readLine().split("\t");
           String str = st[1];
           Double pos1 = Double.parseDouble(st[0]);
           if(Math.abs(pos1-pos)>1.0){
              Logger.global.warning(chr+" warning mismatch "+pos1+" "+pos+"\n"+str);
           }
           for(int k2 = 1; k2<this.buff.length; k2++){
               str = str+"\t"+buff[k2][k].readLine().split("\t")[1];
           }
           res.add(str);
       }
       return res;
    }

    private String[] readAssociationFile(File file) throws Exception{
      BufferedReader br =  new BufferedReader(new FileReader(file));
      String[] res = new String[] {br.readLine(), br.readLine()};
        br.close();
        return res;
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
        String start = chr;//chr.split("_")[0];
       boolean first = true; 
    //   for(int i=0; (st = br.readLine())!=null ; i++){}
        for(int i=0; (st = br.readLine())!=null; i++){
        	if(first){
        		if(st.startsWith("chr")){
        			start = "chr"+start;
        		}
        		first = false;
        	}
        		
        	
        	
            if(!st.startsWith(start)) continue;
            String[] str = st.split("\\s+");
            if(!str[chr_index].equals(start)) continue;
            snps.add(st);
            id.add(str[id_index]);
        }
        if(id.size()==0){
        	throw new RuntimeException("no snps for ");
        }
        br.close();
    }
   
 
    public void run() throws Exception{
        if(split) {
            String st = "";
            int k=0;
            for(int i=0; i<this.id.size(); i++){
                ZipEntry ent = (ZipEntry) zf.getEntry(id.get(i));
                if(ent==null){
                	ent = (ZipEntry) zf.getEntry(chr+"/"+id.get(i)+".txt");
                	if(ent==null){
                		throw new RuntimeException(this.f+" did not contain "+id.get(i));
                	}
                }
                BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
                for(int k1=0; k1<samples.length; k1++){
                    String[] in = br.readLine().split("\t");
                    for(int k2=0; k2<indices.length; k2++){
                    	if(pw[k2][k1]!=null){
                    		String[] str = snps.get(i).split("\\s+");
                    		double res1=Double.NaN;
                    		double res2=Double.NaN;
                    		try{
                    		
                    		 res1 =  Double.parseDouble(in[indices[k2]]);
                    		 res2 = hasgc ? Double.parseDouble(str[str.length-1]) : 0;
                    		// 
                    		}catch(Exception exc){
                    		//	Exception exc1  = new Exception("problem with "+Arrays.asList(str));
                    			//throw exc1;
                    			System.err.println(exc.getMessage());
                    		
                    			
                    			//exc.printStackTrace();
                    		}
                    		this.pw[k2][k1].println(str[1]+"\t"+(res1-res2));
                    		
                    		
                    	}
                    }
           
                }
                br.close();
            }
            for(int k2=0; k2<pw.length; k2++){
                for(int k1 =0; k1<pw[k].length; k1++){
                	if(pw[k2][k1]!=null)
                    pw[k2][k1].close();
                }
            }
    }
        else{
           // this.workOutMedianCorrection();
            this.write();
            
        }
    }
        
    
}
