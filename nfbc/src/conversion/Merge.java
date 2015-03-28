package conversion;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

public class Merge {
	
public static void main(final String[] args){
	if(args.length>=2){
	//File user = new File(args[1]);
	//split(user, args[0], args.length>2 ? new File(args[2]) :  new File(System.getProperties().getProperty("user.dir"))
	//   );
	//}
	//else{
		chromPref = args[1];
	}
		File user = new File(args[0]);
	System.err.println(Arrays.asList(args));
		split(user, "", new File(System.getProperties().getProperty("user.dir")));
	//}
	//ConvertIllumina.es.shutdown();;
}

public static String chromPref=null;
public static void split(File user, final String prefix,final File out){
	try{ 
		out.mkdir();
   final  File[] f = user.listFiles(new FileFilter(){

        public boolean accept(File pathname) {
           return pathname.isDirectory() && pathname.getName().startsWith(prefix) ;
         //  &&   !pathname.getName().endsWith("zip")
        }
        
    });
   FilenameFilter fnf = new FilenameFilter(){

       public boolean accept(File pathname, String name) {
       	System.err.println("checking "+name);
          return (chromPref==null || name.equals(chromPref+".zip")) && name.endsWith("zip");// && name.startsWith("X_M");// && pathname.getName().startsWith(prefix);
       }
       
   };
    final String[] f1 = f[0].list(fnf);
    for(int k=1; k<f.length; k++){
    	String[] f2 = f[k].list(fnf);
    	if(f2.length!=f1.length){
    		throw new RuntimeException("discrepancy in number of zip files "+Arrays.asList(f1)+"\n"+Arrays.asList(f2));
    	}
    }
    System.err.println("prefix "+prefix);
    System.err.println(Arrays.asList(f));
    System.err.println(Arrays.asList(user.list()));
    System.err.println(user.getAbsolutePath());
    System.err.println(chromPref);
    System.err.println(Arrays.asList(f1));
    List<String> types = new ArrayList<String>();
  
	//  Map<String, Integer> m =  readCaseControl(new File(user, "Samples.txt"),toExcl,  "obesity",types);
    final File[] outF = new File[f.length];
    
    ExecutorService es = Executors.newFixedThreadPool(1);;
    List l = new ArrayList();
    for(int i1=0; i1<f1.length; i1++){
    	final int i = i1;
    	Callable call = new Callable(){
    		public Object call(){
    			try{
    	
    	for(int k=0; k<outF.length; k++){
    		outF[k] = new File(f[k], f1[i]);;
    	}
         split(new File(out, f1[i]), outF);
    		
    		}catch(Exception exc){
    			exc.printStackTrace();
    		}
    		return null;}
    	};
    	call.call();
//    	l.add(call);
    }
    
    //es.invokeAll(l);
    es.shutdown();
	}catch(Exception exc){
		exc.printStackTrace();
	}
}


public static void split(File in, File[] out ) throws Exception{
    FileOutputStream dest = new FileOutputStream(in);
    CheckedOutputStream checksum = new CheckedOutputStream(dest, new Adler32());
    ZipOutputStream outS = new ZipOutputStream(new BufferedOutputStream(checksum));
    OutputStreamWriter osw = new OutputStreamWriter(outS);
    List[] resu = new List[out.length];
    ZipFile[] zf = new ZipFile[out.length];
   List<String> indiv = new ArrayList<String>();
 //  PrintWriter plate = new PrintWriter(new BufferedWriter(new FileWriter(new File(in.getParentFile(), in.getName().split("\\.")[0]+".plate"))));
    for(int i=0; i<out.length; i++){
    	try{
    		System.err.println("opening "+out[i]);
    	zf[i] = new ZipFile(out[i]);
    	List<String> ind =Compressor.getIndiv(zf[i], "Samples",out[i].getParentFile().getName()) ;
    	indiv.addAll(ind );
    	//for(int ij=0; ij<ind.size(); ij++){
    	//	plate.println(ind.get(ij)+" "+out[i].getParentFile().getName());
    	//}
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    }
  //  plate.close();
   // if(in.length()==0) return;
    List<String>[] snps = Compressor.getIndiv(zf, "SNPS");
   
    List<String>[] names = Compressor.getIndiv(zf, "Name");
    int[][] alias = new int[names.length][];
    String nme_out = getCommonHeaders(names,0, alias);
   // int loes_ind = Arrays.asList(names[0].get(1).split("\t")).indexOf("loess");
    String head = out[0].getName();
    
   //Compressor.writeEntry("Name", names, out, outS, osw);
   {
   ZipEntry headings = new ZipEntry("Name");
   outS.putNextEntry(headings);
       osw.write(nme_out); osw.write("\n");
       String str1 = names[0].get(1);
       int li = str1.lastIndexOf("\t");
   
       String[] snp_names = names[0].get(1).split("\t");
      
       boolean hasgc = false;
       for(int k=0; k<snp_names.length; k++){
    	   if(snp_names[k].indexOf("gc")>=0){
    		   hasgc = true;
    	   }
       }
       if(hasgc) li = str1.lastIndexOf("\t",li-1);
       String towrite = names[0].get(1).substring(0, li);
       osw.write(towrite);
       for(int j=0; j<names.length; j++){
    	   if(hasgc){
          	osw.write("\t"+out[j].getParentFile().getName()+"_gc");
    	   }
          	osw.write("\t"+out[j].getParentFile().getName()+"_loess");
          }
       osw.write("\n");
       osw.write(names[0].get(2));
       osw.write("\t"); osw.write("Plate"); osw.write("\n");
       osw.flush();
       outS.closeEntry();
   }
  if(snps.length>0 && snps[0]!=null){
   Compressor.writeEntry("SNPS", snps[0],snps,  outS, osw);
  }
   Compressor.writeEntry("Samples", indiv, outS, osw);
    head = head.substring(0, head.indexOf(".zip"));
   
    
    for(Enumeration en = zf[0].entries(); en.hasMoreElements();){
        ZipEntry ent = (ZipEntry) en.nextElement();
        String name = ent.getName();
        if(ent.getName().startsWith("Name") || ent.getName().startsWith("SNPS") || name.startsWith("Sample")) continue;
        List<String> res = new ArrayList<String>();
        for(int i=0; i<out.length; i++){
        	if(zf[i]!=null){
        	try{
        	 res.addAll(Compressor.getIndiv(zf[i], ent.getName(), alias[i]));
        	}catch(Exception exc){
        		exc.printStackTrace();
        	}
        	}
        }
      Compressor.writeEntry(name, res, outS, osw);
        
    }
    String sample_header = null;
   
    osw.flush();outS.close();
    System.err.println("done "+in);
    
}

private static String getCommonHeaders(List<String>[] names, int i,
		int[][] alias) {
	List<String> res = new ArrayList<String>();
	StringBuffer sb = new StringBuffer();
	for(int k=0; k<names.length; k++){
		List<String> str = Arrays.asList(names[k].get(i).split("\t"));
		for(int j=0; j<str.size(); j++){
			if(!res.contains(str.get(j))){
				res.add(str.get(j));
				sb.append("\t"+str.get(j));
			}
		}
	}
	for(int k=0; k<names.length; k++){
		alias[k] = new int[res.size()];
		List<String> str = Arrays.asList(names[k].get(i).split("\t"));
		for(int j=0; j<res.size(); j++){
			alias[k][j] = str.indexOf(res.get(j));
		}
	}
	return sb.toString().substring(1);
}


private static File[] getOutFiles( File user, List<String> types ){
    File [] out = new File[types.size()];
    for(int i=0; i<out.length; i++){
        out[i] = new File(user, types.get(i).replaceAll("\\s+", "_"));
    }
    return out;
}
}
