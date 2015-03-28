package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

/** merges 1.zip, 2.zip ... into a single file */
public class MergeChroms {
	
	
	public static int threads = 1;
	public static ExecutorService es;
       
        public static void involeTasks(List l, boolean seq) throws Exception{
        	
           if(!seq && threads>1) es.invokeAll(l);
            
           else{
           for(Iterator it = l.iterator(); it.hasNext();){
               try{
                ((Callable)it.next()).call();
               }catch(Exception exc){
                   exc.printStackTrace();
                   System.exit(0);
               }
            }
           }
        }
       
	public static void main(final String[] args){
		es = Executors.newFixedThreadPool(threads);;
		List l = new ArrayList();
		 final File user = new File(System.getProperty("user.dir"));
	        //main(user,args);
	        File[] f = user.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.isDirectory()&& !pathname.getName().startsWith("comparison") && !pathname.getName().startsWith("res")
					&& (args.length==0 || pathname.getName().indexOf(args[0])>=0);
				}
				
			});
	        //if(!all){
			for(int i=0; i<f.length; i++){
				File res = new File(f[i], "res");
				if(!res.exists()) continue;
				File[] f1 = res.listFiles(new FileFilter(){

					public boolean accept(File pathname) {
						return pathname.isDirectory();
					}
					
				});
				for(int j=0; j<f1.length; j++){
					l.add(main(f1[j]));
				}
			}
			try{
			involeTasks(l, true);
			}catch(Exception exc){
				exc.printStackTrace();
			}
			es.shutdown();
	}
	
	public static void main1(String args[]){
		 File dir1 = new File(System.getProperties().getProperty("user.dir"));
		 Callable call = main(dir1);
		 try{
		 call.call();
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
	}
	
    public static Callable main(final File dir1){
    	return new Callable(){
    		public Object call(){
		        try{
		        //	 File dir_ = new File(System.getProperties().getProperty("user.dir"));
		        //	 File dir1 = new File(dir_, args[0]);
		        	final File[] dirs =  dir1.listFiles(new FileFilter(){
		
							public boolean accept(File pathname) {
							   return  pathname.getName().endsWith("zip") && ! pathname.getName().startsWith("all") && !pathname.getName().startsWith("X");
							}
		             	
		             });
		        Arrays.sort(dirs, chromComparator);
		        	//File output = new File(dir1.getParent(), "res1");
		        //	output.mkdir();
		        	//File outputdir = new File(output, args[0]);
		//        	outputdir.mkdir(); 
		        boolean checkZip = false;
		        if(dirs.length>0){
		        	 File res = new File(dir1, "all.zip");
		        	 boolean conti=res.exists();
		        	 if(conti && checkZip){
		        		 
		        	
		        		 try{
		        		 ZipFile zf = new ZipFile(res);
		        		 zf.entries();
		        		 zf.close();
		        	 }catch(Exception exc){
		        		 conti = false;
		        	 }
		        	 }
		        	 if(!conti){
		        	try{
			        MergeChroms chmp = new MergeChroms(dir1, dirs);
			       
			        	chmp.run();
		        	}catch(Exception exc){
		        		exc.printStackTrace();
		        	}
			        }
			        
		        	 
		        	 else{
				        	System.err.println("skipping "+res.getAbsolutePath());
				        }
		        }
		        
		        }catch(Exception exc){
		            exc.printStackTrace();
		        }
		        return null;
    		}
    	};
    }
    static Comparator<File> chromComparator = new Comparator<File>(){

		public int compare(File o1, File o2) {
			return getNumber(o1.getName()).compareTo(getNumber(o2.getName()));
		}

		
		 
	 };
    
   
    public static Integer getNumber(String o) {
    	String o1 = o.substring(0, o.indexOf(".zip"));
		if(o1.equals("X")) return 23;
		else if(o1.equals("XY")) return 24;
		else if(o1.equals("Y")) return 25;
		else if(o1.equals("M")) return 26;
		else return Integer.parseInt(o1);
	}
    public static String getString(int o1){
    	if(o1==23) return "X";
		else if(o1==24) return "Y";
	//	else if(o1==25) return "";
	//	else if(o1==26) return "M";
		else return o1+"";
    }
   final  File dir;
    FileOutputStream dest;
   final  CheckedOutputStream checksum;
    final ZipOutputStream outS;
    final OutputStreamWriter osw;
    ///String chr;
    
 String[] header;
 
  ZipFile f;
  File[] f1;
 // int id_index;
//  int index = 0;
 
 // BufferedReader buildF;
  List<String> build_header;
  int chr_index;
  int pos_index;
  int id_index;
  String currentChr = "";
  PrintWriter snps = null;
  File snpsF ;
  File excl;
    MergeChroms(File dir, File[] f) throws Exception{
        this.dir = dir;
        this.f1 = f;
       snpsF = new File(dir, "SNPS");
       snps = new PrintWriter(new BufferedWriter(new FileWriter(snpsF)));
       
     //   id_index = Arrays.asList(header).indexOf("snpid");
       // chr = prefix.split("_")[1].substring(3);
       File res = new File(dir, "all.zip");
        dest = Compressor.getOS1(res);
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
       outS.setMethod(ZipOutputStream.DEFLATED);
       this.f = new ZipFile(f[0]);
       writeEntry(this.f.getEntry("Name"));
  	   writeEntry(this.f.getEntry("Samples"));
  	 
    }
    public void copySNPS(File f) throws Exception{
      	 ZipEntry headings = new ZipEntry(f.getName());
           outS.putNextEntry(headings);
          
      	BufferedReader br = new BufferedReader(new FileReader(f));
      	String st = "";
      	while((st = br.readLine())!=null){
      		 osw.write(st);
               osw.write("\n");
      	}
      br.close();
      	 osw.flush();
           outS.closeEntry();
      }
  
    
   
   
    public void writeHeader(File[] f) throws Exception{
        ZipEntry headings = new ZipEntry("Name");
        outS.putNextEntry(headings);
            for(int i=0; i<f.length; i++){
                osw.write(f[i].getName());
                osw.write(i<f.length-1 ? "\t" : "\n");
            }
            osw.flush();
           outS.closeEntry();
    }
 
  
   
  
   /**tmp storage */
 //  final ZipEntry ents_;
  BufferedReader br;
   String str;
   
   public void  writeSNPS(int index ) throws Exception{
	   ZipEntry ents_ = f.getEntry("SNPS");
	   this.toexcl.clear();
	  // PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(snpsF, index>0)));
	   InputStream is =   f.getInputStream(ents_);
		  br =  new BufferedReader(new InputStreamReader(is));
	   String st = "";
		 
	   while((st = br.readLine())!=null){
			
			/* if(st.startsWith(exclChrom) && false){
				 String[] str = st.split("\\s+");
				int pos = Integer.parseInt(str[1]);
				if(pos>= startExcl && pos <=endExcl){
					this.toexcl.add(str[3]);
					System.err.println("EXCLUDING "+str);
				}
			 }*/
			  this.snps.write(st);
			  this.snps.write("\n");
		 }
		snps.flush();
		br.close();
		 
   }
   
   public boolean writeEntry(ZipEntry ents_) throws Exception{
	   
		  String name = ents_.getName();
			//ZipEntry   ents_ = f.getEntry(name);
			  
			   InputStream is =   f.getInputStream(ents_);
				  br =  new BufferedReader(new InputStreamReader(is));
		  
		   ZipEntry headings = new ZipEntry(name);
		  // headings.setComment(det);
		   outS.putNextEntry(headings);
//		   if(det!=null) osw.write(det+"\n");
		   String st = "";
			 while((st = br.readLine())!=null){
				 
				  osw.write(st);
				  osw.write("\n");
			 }
			 br.close();
			 is.close();
			 osw.flush();
			  outS.closeEntry();
			  return true;
	   }
   
    public void run() throws Exception{
    	for(int i=0; i<f1.length; i++){
    		run(i);
    	}
    	this.f.close();
    	this.snps.close();
    	copySNPS(snpsF);
    	this.osw.close();
    	outS.close();
    	snpsF.delete();
    }
    
  Set<String> toexcl = new HashSet<String>();
  int startExcl = 42000000; String exclChrom = "chr7"; int endExcl = 82000000;
    public void run(int i) throws Exception{
    	//String st = "";
    	System.err.println(f1[i].getAbsolutePath());
    	if(this.f!=null) f.close();
    	this.f = new ZipFile(this.f1[i]);
    	
    	this.writeSNPS(i);
    	String chrom = f1[i].getName().split("\\.")[0];
    	inner: for(Enumeration ent = f.entries(); ent.hasMoreElements();){
    		ZipEntry ent1 = (ZipEntry)ent.nextElement();
    		String nme = (ent1).getName();
    		if(toexcl.contains(nme)) continue inner;
    		if(!nme.startsWith("SNP") && ! nme.startsWith("Sam") && ! nme.startsWith("excl") && !nme.startsWith("Nam")){
    			writeEntry(ent1);
    		}
    		
    	}
    	
    }

	
}
