package lc1.assoc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import lc1.dp.data.collection.DataCInterface;
import lc1.dp.data.representation.PIGData;
import lc1.stats.Regression;
import lc1.util.Compressor;
import lc1.util.Constants;

public class VirtualDataCollection implements DataCInterface {
 
	 private static int substr = 1;
    
    int[] mid = new int[] {0,Integer.MAX_VALUE};
    int[] kb = new int[] {0,0};
    public final String chrom;
 
    public static boolean check(File out) {
		if(out.exists()){
			//ZipFile zf=null;
    		try{
    			if(out.isDirectory()) return false;
    			ZipFile zf = new ZipFile(out);
    			List<String> l = Compressor.getIndiv(zf,"SNPS", null);
    			zf.close();
    			if(l.size()>0) {
    				return true;
    			}
    		
    		}catch(Exception exc){
    			//if(zf!=null) zf.close();
    			//exc.printStackTrace();
    		}
    	}
		return false;
	}
    
    public static void main(String[] args){
    	boolean all = true;
        final File user = new File(System.getProperty("user.dir"));
       final  String date = args[0];
       final String chroms = args.length==1 ? "": args[1];
        //main(user,args);
        File[] f = user.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.isDirectory()&& pathname.getName().endsWith(date);
			}
			
		});
        //if(!all){
		for(int i=0; i<f.length; i++){
			if(
					//f[i].getName().indexOf("-3-2")>=0
					//||
					f[i].getName().indexOf("20")>=0 
					){

				main(f[i], chroms, 1);
			}
		}
        //}
       /* else{
        	main(user, args, 0);
        }*/
		//main(user, args);
    }
    
    public static void main(final File user, String chroms, int substring){
        try{
        	String[] probeset = getProbes(user);
           //Constants.modelCNP = 10;
            Constants.plot = 2;
            Constants.numThreads = 1;
         //   Constants.topBottom = false;
            Constants.median_correction = new boolean [] {false};
            Constants.loess = new boolean [] {false};
        substr = substring;
            final String[] chrom = chrom(user, false, chroms);
            final File outDir = new File(user, "res");
            outDir.mkdir();
           String[] str = new String[] {""};//, "0", "1", "2", "3", "4", "5", "6", "7","8","9"};
          //  for(int j=0; j<str.length; j++){
           List l = new ArrayList<Callable>();
           
            for(int i=0; i<chrom.length;
            //chrom.length; 
            i++){
            	for(int j=0; j<probeset.length; j++){
            		final String probset = probeset[j];
            	final int i1 = i;
            	//for(int k=0; k<type.length; k++){
            		//for(int l=0; l<sexP.length; l++){
            	Callable call = new Callable(){
            		public Object call(){
            		try{
            			//String out_name =  "assoc.txt";
            			//String probeset = "1M";
            			File outDir1 = new File(outDir, probset);
            			outDir1.mkdir();
                    	File out = new File(outDir1, chrom[i1]+".zip");
                    	if(check(out)) return null;
                    
		                VirtualDataCollection vdc =  new VirtualDataCollection(user,1,chrom[i1], probset);
		               // vdc.pheno = pheno;
		             //  File outF = new File(outDir, chrom[i1]);
		            //   outF.mkdir();
		                if(vdc.result.length>0) vdc.copy(out);
		            //	CompressDC cdc = new CompressDC(outDir,vdc);
		             //  cdc.run();
                    	
               
            		}catch(Exception exc){
            			exc.printStackTrace();
            		}
            		return null;
            		}

					
            };
            l.add(call);
            	//}
            	//}
         //   }
          /* AssociationFrame app = new AssociationFrame(vdc, chrom, "regrP", "assoc"+str[j]+".txt");//"armitage");
           app.pack();
           app.setVisible(true);*/
            	}
         
           }
            if(Constants.numThreads==1){
            	for(int i=0; i<l.size(); i++){
            		((Callable)l.get(i)).call();
            	}
            }
            else{
            	Regression.es.invokeAll(l);
            }
           // app.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        }catch(Exception exc){
            exc.printStackTrace();
        }
      // Executor.shutdown();
    }
    
    private static String[] getProbes(File user2) {
		Set<String> probes = new HashSet<String>();
		FilenameFilter filter = new FilenameFilter(){

			@Override
			public boolean accept(File f, String name) {
				return name.endsWith("zip");
			}
			
		};
		File[] f = user2.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && ((new File(pathname, "res"))).exists();
			}
			
		});
		
		for(int i=0; i<f.length; i++){
			probes.addAll(Arrays.asList(new File(f[i], "res").list(filter)));
		}
		String[]res = new String[probes.size()];
		Iterator<String> it = probes.iterator(); 
		for(int i=0; i<res.length; i++){
			String st = it.next();
			res[i] = st.substring(0, st.length()-4);
		}
		return res;
	}

	ZipOutputStream outS = null;
    OutputStreamWriter osw = null;
    CheckedOutputStream checksum;
    public  void copy(ZipFile zf, String id) throws Exception{
    
 	   copy( Compressor.getIndiv(zf, id, null),id);
 	  
    }
    public void copy(List<String> l, String id) throws Exception{
    	ZipEntry headings = new ZipEntry(id);;
  	   outS.putNextEntry(headings);
  	   for(int i1=0; i1<l.size(); i1++){
  		osw.write(l.get(i1)+"\n") ;  
  	   }
  	   osw.flush();
  	  outS.closeEntry();
    }
    public  void copy(File out) throws Exception{
        //List<String> excs = Arrays.asList(exception);
       for(int i=0; zf==null; i++){
    	   try{
    	   zf = new ZipFile(result[i]);
    	   }catch(Exception exc){
    		   System.err.println("Warning"+exc.getMessage());
    		   //exc.printStackTrace();
    	   }
       }
        FileOutputStream dest = Compressor.getOS(out);
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
         ZipOutputStream(new 
           BufferedOutputStream(checksum));
       osw = new OutputStreamWriter(outS);
       outS.setMethod(ZipOutputStream.DEFLATED);
      String[] toCopy = new String[] {"Samples", "Name"};
       for(int i=0; i<toCopy.length;i++){
    	   
    	   copy(zf,toCopy[i]);
       }
       for(int i=0; i<fs.length; i++){
    	   try{
    		   zf = new ZipFile(result[i]);
    		
              List<String > snps =Compressor.getIndiv(zf, "SNPS", null);
                
               
                inner: for(int j=0; j<snps.size() ; j++){
                	String[] str = snps.get(j).split("\t");
                	int pos = Integer.parseInt(str[1]);
                	if(pos>=this.start[i]){
                		if(pos >=this.end[i]) break inner;
                		this.snpid.add(snps.get(j));
                		copy(zf, str[3]);;
                	}
                }

    	   }catch(Exception exc){
    		   System.err.println("warngin "+exc.getMessage());
    	   }
       }
       copy(this.snpid,"SNPS");
        outS.close();
        
      
    }
    public static String[] chrom(File f, boolean lng, final String names){
    	final String[] name = names.split(":");
        File[] fs =  f.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
              return !pathname.getName().startsWith("res") && (pathname.isDirectory())&& !pathname.getName().startsWith("res")
              		&& !pathname.getName().startsWith("compa");
              //            		  pathname.getName().endsWith("zip") && startsWith(name, pathname.getName());
            }

			private boolean startsWith(String[] name, String name2) {
				for(int i=0; i<name.length; i++){
					if(name2.substring(substr).startsWith(name[i]) ) return true;
				}
				return false;
			}
            
        });
        Set<String> s = new HashSet<String>();
        for(int i=0; i<fs.length; i++){
            String[] str = fs[i].getName().substring(substr).split("_");
            if(str.length<=1) continue;
            s.add(str[0]+(lng ? "_"+str[1] : ""));
        }
        return s.toArray(new String[0]);
    }
    public String name = "-";
    public Integer length;
   public List<String> snpid = new ArrayList<String>();
  //  public List<Integer> loc = new ArrayList<Integer>();
   // public List<Character>majorAllele = new ArrayList<Character>();
    //public List<Character> minorAllele = new ArrayList<Character>();
  //  public List<Character> majorAllele = new ArrayList<Character>();
  //  public List<Character> minorAllele= new ArrayList<Character>();
  //  List<List<String>>[] info;
   // DataCollection[] dc;
    File[] fs;
    int col;
    
   
 private boolean contains(Collection<String> allowed, String string) {
		for(Iterator<String> it = allowed.iterator(); it.hasNext();){
			if(string.indexOf(it.next())>=0) return true;
		}
		return false;
	}
/*   Map<Integer, Integer> leftPos = new HashMap<Integer, Integer>();
    Map<Integer, Integer> rightPos = new HashMap<Integer, Integer>();
    Map<Integer, Integer> mid = new HashMap<Integer, Integer>();*/
   
    File[] result;
    //DataCollection[] dc;
    
    
    public File[] getFiles(File f, final String chrom, final String probeset){
    	Set<String> s= new HashSet();
    	File[] fs =  f.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
              if((pathname.isDirectory() || pathname.getName().endsWith(".zip")) && pathname.getName().substring(substr).startsWith(chrom+"_")){
            	  File f1 = new File(pathname,"res/"+probeset+".zip") ;
            	  return f1.exists();//check(pathname);
              }
            	  else return false;
             
            
            }
            
        });
    	for(int i=0; i<fs.length; i++){
    		s.add(fs[i].getName().split(".zip")[0]);
    		
    	}
    	File[]res = new File[s.size()];
    	Iterator<String> it = s.iterator();
    	for(int i=0; i<res.length; i++){
    		res[i] = new File(f, it.next());
    	}
    	return res;
    }
    
    static Comparator<File> COMPA = new Comparator<File>(){

		public int compare(File o1, File o2) {
			String[] st1 = o1.getName().substring(substr).split("_");
			String[] st2 = o2.getName().substring(substr).split("_");
			if(!st1[0].equals(st2[0])) throw new RuntimeException("!!");
			int i1 = Integer.parseInt(st1[1]);
			int i2 = Integer.parseInt(st2[1]);
			if(i1==i2)return 0;
			else return i1<i2  ? -1 : 1;
		}
    	
    };
    //String sexV;
  //  int[] startIndex;//which position we start
   // int[] endIndex;//which position we end in this file
    List<Integer> fileIndex = new ArrayList<Integer>();
    List<Integer>posIndex = new ArrayList<Integer>();
    ZipFile zf;
    public VirtualDataCollection(File f, int col, final String chrom, final String probeset) throws Exception{
        Constants.mid = 
        	new String[][]{
        	new String[] {chrom.split("_")[0], "0", ""+(Integer.MAX_VALUE -1)}};
       fs = getFiles(f, chrom, probeset);
       this.start = new int[fs.length];
       this.end = new int[fs.length];
      // this.startIndex = new int[fs.length];
       //this.endIndex = new int[fs.length];
       Arrays.sort(fs, COMPA);
      // this.sexV = sex;
       this.col = col;
       this.user = f;
       this.chrom = chrom.split("_")[0];
       //this.dc = new DataCollection[fs.length];
       this.result = new File[fs.length];
     
     //  this.info = new List[fs.length];
       for(int i=0; i<fs.length; i++){
    	   try{
           String nme1 = fs[i].getName();
          if(nme1.indexOf("core")>=0){
        	  	String[] str = nme1.split("core_")[1].split("_");
        	  	start[i] = Constants.convert(str[0]);
        	  	end[i] = Constants.convert(str[1]);
          }
          else{
           String[] spl = nme1.substring(substr).split("_");
          
           
           if(i==0)start[i] = Constants.convert(spl[1]);
           else{
        	   int mind = (int) Math.round((Double.parseDouble(spl[1])+(double)end[i-1])/2.0);
        	   end[i-1] = mind;
        	   start[i] = mind;
           }
           end[i] =Constants.convert(spl[2]);
          }
         File zipFile = new File(fs[i].getAbsolutePath()+".zip");
           ZipFile zf;
           if(zipFile.exists()){
        	   zf = new ZipFile(zipFile) ;
           }
           else{
        	   zf = null;
           }
         //  File nme = new File(nme1+"/"+assoc);
           
           result[i] = new File(fs[i],"res/"+probeset+".zip") ;
         if(!result[i].exists()){
          	 File zipFile1 = new File(fs[i].getAbsolutePath()+".zip");
          	 ZipFile zf1 = new ZipFile(zipFile1);
          	 String resEntry = CalcAssociation.findEntry(zf1, fs[i].getName().split("_")[0]+".zip");
               File dirOut = new File(user, fs[i].getName());
               if(resEntry!=null){
                  
                   dirOut.mkdir();
                   result[i] = CalcAssociation.expand(dirOut, zf,  resEntry);
               }
         }
          //  dc[i]= new LightWeightDataCollection(result[i],(short) 0, 2, mid,  null);
           //dc[i] =  new SimpleDataCollection() ;
         //  System.err.println(dc[i].data.size());
           
          
          
          
              
              
       }catch(Exception exc){
    	   exc.printStackTrace();
    	  // System.err.println("prob with "+fs[i]);
    	  // System.exit(0);
       }
       }
       this.snpid = new ArrayList<String>();
      
       System.err.println("snpid loc");
    }
    
    private String conc(String[] strn, int i, String join) {
       StringBuffer sb = new StringBuffer();
       for(int k=i; k<strn.length; k++){
           sb.append(strn[k]);
           if(k<strn.length-1) sb.append(join);
       }
       return sb.toString();
    }

   

    
   
static final Comparator<Map.Entry<String, List<PIGData>>> compa = new Comparator<Map.Entry<String, List<PIGData>>>(){

    public int compare(Entry<String, List<PIGData>> o1, Entry<String, List<PIGData>> o2) {
         int sz1 = o1.getValue().size();
         int sz2 = o2.getValue().size();
         if(sz1==sz2) return 0;
         else return sz1 < sz2 ? 1 :-1;
    }
    
};
public  static final Comparator< Entry<String, List<int[]>>> compa1 = new Comparator<Entry<String, List<int[]>>>(){

    public int compare(Entry<String, List<int[]>> o1,
            Entry<String, List<int[]>> o2) {
        // TODO Auto-generated method stub
      Integer i1 = o1.getValue().get(0)[0];
      Integer i2 = o2.getValue().get(0)[0];
       return i1.compareTo(i2);
    }
    
};
 File user;
int extra = 30000;

final int[] start;
final int[] end;
public double getAverageNumberCN(){
	return 0;
}
   
    public List<Integer> loc() {
       return null;
    }

    public String chrom(int i) {
    	return chrom;
//       return this.dc[0].chrom;
    }

	public String chrom() {
		return chrom;
	}

	public String getCompressedString(String key, int i, boolean b, boolean b1) {
		return null;
	//	int pos = this.fileIndex.get(i);
//		if(dc[pos]==null)return "---";//
//	return this.dc[pos].getCompressedString(key, this.posIndex.get(i), b);
	}

	public String headSNP() {
		return null;
		// TODO Auto-generated method stub
		//return this.dc[0].headSNP();
	}

	public String head_snp() {
		return null;
		// TODO Auto-generated method stub
		//return this.dc[0].head_snp();
	}

	public List<Character> alleleA() {
		return null;//return majorAllele;
	}

	public List<Character> alleleB() {
		return null;//return minorAllele;
	}

	public String name() {
		// TODO Auto-generated method stub
		return this.chrom;
	}

	public List<String> snpid() {
		// TODO Auto-generated method stub
		return snpid;
	}

	public List<String> indiv() {
		return null;//this.dc[0].indiv();
	}
    
}
