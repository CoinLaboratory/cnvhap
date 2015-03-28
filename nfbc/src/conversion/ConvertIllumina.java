  package conversion; 

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import ensj.GenesForRegion;
import ensj.Regression1;


public class ConvertIllumina {
	
	/*public static final Options OPTIONS  = new Options(){
        {
                    this.addOption( OptionBuilder.withLongOpt( "todo" ).withDescription( "todo").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "build" ).withDescription( "build").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "output" ).withDescription( "output").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "mode" ).withDescription( "mode").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "runR" ).withDescription( "runR").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "phenFile" ).withDescription( "phenFile").withValueSeparator( ':' ).hasArgs().create());
                    
        }
    };  */  
	public final static boolean delete = false;
///	public static Map<String, List<String>> plate = null;
public static class ZipFileFilter implements FileFilter{
public ZipFileFilter(String prefix){
	pref = prefix;
}
String pref;
	@Override
	public boolean accept(File pathname) {
		boolean dir = pathname.isDirectory() && (
				pref.length()==0 || 
				pathname.getName().equals(pref));
		if(dir){
			String[] f = pathname.list();
			for(int i=0; i<f.length; i++){
				if(f[i].endsWith(".zip")) return true;
			}
		}
		return false;
	}
	
}
	public static void  split(final List<String> mode, final String todo, final File dir, final File outdir, final File pdf) throws Exception{
		 // for(int k1=0; k1<todo.length; k1++){
			 
              if(mode.contains("qc")){
            	  if(todo.equals("X")){
            		 ConvertToWide.runR(dir, todo, pdf, false, false, 0.35);
            	
            	  }
            	  else{
            		  System.err.println("runR "+dir+" "+todo);
            		  ConvertToWide.runR(dir, todo, pdf, false, false, 0.35);
            	  }
              }
              if(mode.contains("join")){
            	  System.err.println("join on "+todo+" "+dir.getAbsolutePath()+" output "+outdir.getAbsolutePath()+" "+pdf.getAbsolutePath());
            	  if(todo.equals("X") && false){
            		  Reconstitute.main(dir, outdir, todo, "_M");
            		  Reconstitute.main(dir, outdir, todo, "_F");
            	  }
            	  else{
            		  Reconstitute.main(dir, outdir, todo, "");
            	  }
              }
					
				  
			
      
              
  }
		  
		  public static void split( final String chr, final File dir, final File outdir, final File pdf ) throws Exception{
				
		            	  ConvertToWide.main(true, dir, chr);
		         
		  }
	public static int threads = 4;
	public static ExecutorService es;
       
        public static void involeTasks(List l, boolean seq) throws Exception{
        	
           if(!seq && threads>1 && es!=null) es.invokeAll(l);
            
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
        
        public static void involeTasks(Collection l, boolean seq) throws Exception{
        	
            if(!seq && threads>1 && es!=null) es.invokeAll(l);
             
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
	public static void main(String[] args){
		try{
			OptionBuild.getOptions(args);
		//	File plate = OptionBuild.plateFile;
			//if(plate!=null && plate.exists()){
				//ConvertIllumina.plate = readPlateFile(plate);
			//}
			es = threads== 0 ? null:   Executors.newFixedThreadPool(OptionBuild.numThreads);;
			   File dir1 =  OptionBuild.output[0].getParentFile();
		         final  File[] outputDir = OptionBuild.output;
		    //    	if(dir1==null) delete =false;
		         final String build = OptionBuild.build;
		         final String phenFile = OptionBuild.phenFile;
		       //  final String[] todo = OptionBuild.todo;
		         final List<String>mode = Arrays.asList( OptionBuild.mode);
		         Boolean runR = Boolean.parseBoolean(OptionBuild.runR);
		        // for(int ik=0; ik<outputDir.length; ik++){
		         File parent = outputDir[0].getParentFile();
		         if(parent!=null)parent.mkdir();
		        Collection<Callable> l = new ArrayList();
		         
		         if(mode.contains("likelihood")){
		        	ConvertBuckleyToLong.main(OptionBuild.file);
		        	return;
		         }
		         if(mode.contains("compress")){
		        	
		        	 if(outputDir.length>1) {
		        		 throw new RuntimeException("!!");
		        	 }
		            // for(int k=todo.length-1; k>=0;k--){
		            	final  File nme =  outputDir[0];
		            //	final String tod = todo[k];
		            	 Callable call = new Callable(){
		            		 public Object call(){
				            	 try{
				                 if(OptionBuild.format.equals("long")){
				                	 CompressIlluminaTableLong.run(OptionBuild.file[0],nme, build, phenFile);
				                 }
				                 else{
				                	 if(OptionBuild.format.equals("widesorted")){
				                		 CompressIlluminaTable.run(OptionBuild.file,nme, build, phenFile);

				                	 }
				                	 else
				                	 CompressIlluminaTable.run(OptionBuild.file,nme, build, phenFile);
				                 }
				            	 }catch(Exception exc){
				            		// System.err.println("done "+k);
				            	 }
				            	 return null;
		            		 }	
		            	 };
		            	 l.add(call);
		            // }
		        // }
		            	  involeTasks(l, true);
		 		         l.clear();
		 		        
		         }
		      
		   //   if(mode.size()>1 || !mode.get(0).equals("compress")) 
		     // 
		       if(mode.contains("build")){
		        	  File buildF = BuildConversion.makeBuildFile( outputDir[0],build);
		          }
		      
		       if(mode.contains("gc")){
		        	 GenesForRegion.main( outputDir[0]);
		          }
		         if(mode.contains("sep")){
		        	
		        	 if(outputDir.length>1) throw new RuntimeException(" output dir length should be zero if we are separating!"+Arrays.asList(outputDir));
		        	 String[] todo = getTodo( outputDir[0]);
		        	// System.err.println("running sep "+Arrays.asList(todo));
		        	 for(int k=0; k<todo.length; k++){
		        		 final  File nme =  outputDir[0];
			            	final String tod = todo[k];
		        		 Callable call = new Callable(){
		            		 public Object call(){
				            	 try{
				            		 return SplitByPheno.split(nme, OptionBuild.plateFile, tod, OptionBuild.plateId, OptionBuild.sampleId, OptionBuild.sampleId1);
				            	 }catch(Exception exc){
				            		 return null;
					            	 }
					            	
			            		 }	
			            	 };
			            	 l.add(call);
		        	 }
		        	 for( Iterator<Callable> it = l.iterator(); it.hasNext();){
		        		 it.next().call();
		        	 }
//		        	  List<Future> fut = es.invokeAll(l);
				       //  outputDir = (File[]) fut.get(0).get();
		         }
		       
		         l.clear();
		       //  if(todo.contains("qc")){
		        final  File outdir = new File(dir1, "modified");
		         final File pdf = new File(dir1, "pdf");
		         pdf.mkdir();
		         outdir.mkdir();
		        File[] dir = outputDir;
		        ConvertToWide.setOptions(true, build);
		        if(mode.contains("splitpq")){
		        	for(int kk=0; kk<dir.length; kk++){
			        	 File[] dirs = dir[kk].listFiles(new ZipFileFilter(OptionBuild.prefix));
			        	 for(int ik=0; ik<dirs.length; ik++){
			        		 String[] todo = getTodo(dirs[ik]);
			         inner: for(int k=0; k<todo.length; k++){
			        	 final String tod = todo[k];
			        	 if(tod.endsWith("q") || tod.endsWith("p") || tod.indexOf("XY")>=0 || tod.indexOf("M")>=0|| tod.indexOf("N")>=0) continue inner;
			        	 final File di = dirs[ik];
			        	 Callable call = new Callable(){
		            		 public Object call(){
				            	 try{
								       SplitByLocation.split( tod, di,   new File(outputDir[0], OptionBuild.build+".gc.txt"),
								    		   new File(outputDir[0], "karyotypes_"+OptionBuild.build+".gc.txt"),
								    		   true);
				            	 }catch(Exception exc){
				            		exc.printStackTrace();
					            	 }
				            	 return null;
			            		 }	
			            	 };
			            	 l.add(call);
				        
			         }
			        	 }
			        }
			        involeTasks(l, true);
			         l.clear();
		        }
		        if(mode.contains("regr")){
		        	  for(int kk=0; kk<dir.length; kk++){
		        		  System.err.println("running regr");
		        		  File[] dirs = dir[kk].listFiles(new ZipFileFilter(OptionBuild.prefix));
		        		  if(dirs.length==0) throw new RuntimeException("no files to run it on "+OptionBuild.prefix);
		        		  for(int i=0; i<dirs.length; i++){
		        			final   File plateDir = dirs[i];
		        			  Callable call = new Callable(){
		 	            		 public Object call(){
		 			            	 try{
		 							      Regression1.main(plateDir);
		 			            	 }catch(Exception exc){
		 			            		
		 				            	 }
		 			            	 return null;
		 		            		 }	
		 		            	 };
		 		            	 l.add(call);
		        		  }
		        	  }
		        	   involeTasks(l, true);
				         l.clear(); 
		        }
		       
		       if(mode.contains("split")){
		        for(int kk=0; kk<dir.length; kk++){
		        	 File[] dirs = dir[kk].listFiles(new ZipFileFilter(OptionBuild.prefix));
		        	 for(int ik=0; ik<dirs.length; ik++){
		        		 String[] todo = getTodo(dirs[ik]);
		         for(int k=0; k<todo.length; k++){
		        	 final String tod = todo[k];
		        	 final File di = dirs[ik];
		        	 Callable call = new Callable(){
	            		 public Object call(){
			            	 try{
							       split( tod, di, outdir, pdf);
			            	 }catch(Exception exc){
			            		
				            	 }
			            	 return null;
		            		 }	
		            	 };
		            	 l.add(call);
			        
		         }
		        	 }
		        }
		        involeTasks(l, true);
		         l.clear();
		       }
		         for(int kk=0; kk<dir.length; kk++){
		        	 File[] dirs = dir[kk].listFiles(new ZipFileFilter(OptionBuild.prefix));
		        	 for(int ik=0; ik<dirs.length; ik++){
		        		 String[] todo = getTodo(dirs[ik]);
			         for(int k=0; k<todo.length; k++){
			        	 final String tod = todo[k];
			        	 final File di = dirs[ik];
			        	 Callable call = new Callable(){
		            		 public Object call(){
				            	 try{
								        split(mode, tod, di, outdir, pdf);
				            	 }catch(Exception exc){
				            		exc.printStackTrace();
					            	 }
				            	 return null;
			            		 }	
			            	 };
			            	 l.add(call);
				        
			         }
		        	 }
			        }
			        involeTasks(l,false);
			         l.clear();
			         if(mode.contains("clean") ){
			        	if(pdf!=null ){
			        	// File pdf = new File(parent, "pdf");
			        	 Utils.delete(pdf);
			        	}
			        	  for(int kk=0; kk<dir.length; kk++){
			        		  if(mode.contains("sep")){
			        			  String[] todo = getTodo(dir[kk]);
			        			  for(int kkk=0; kkk<todo.length; kkk++){
			        				  File todel = new File(dir[kk], todo[kkk]+".zip");
			        				  todel.delete();
			        			  }
				        		}
					        	 File[] dirs = dir[kk].listFiles(new ZipFileFilter(OptionBuild.prefix));
					        	for(int ik=0;ik<dirs.length; ik++){
					        	 String[] todo = getTodo(dirs[ik]);
					        	 for(int kkk=0; kkk<todo.length; kkk++){
					        		if(mode.contains("join")){
					        			File todel = new File(dirs[ik],todo[kkk]);
					        			File todel1 = new File(dirs[ik],"err"+todo[kkk]+".txt");
					        			System.err.println("deleting "+todel.getAbsolutePath());
					        			System.err.println("deleting "+todel1.getAbsolutePath());
					        			Utils.delete(todel);Utils.delete(todel1);
					        		}
					        		if(mode.contains("regr")){
					        			File todel = new File(dirs[ik],"SNPS_"+todo[kkk]);
					        			System.err.println("deleting "+todel.getAbsolutePath());
					        			Utils.delete(todel);
					        		}
					        		
					        	 }
					        	
					        	}
					        	/* if(mode.contains("join") && outdir!=null){
					        		 File buildF = new File(dir[kk],OptionBuild.build+".gc.txt");
					        		 buildF.renameTo(new File(outdir, OptionBuild.build+".txt"));
					        	  }*/
			        	  }
			        	  
				       }
		                   // Utils.delete(new File(dir, todo[k]));
			         if( mode.contains("join") || mode.contains("qc")){
			        	 if(delete){
		          Utils.delete(outputDir[0]);
		          File mod = new File(outdir,outputDir[0].getName());
		          mod.renameTo(outputDir[0]);
			        	 }
			         }
		         
		       //   (new File(dir1,"modified/")).delete();
		         outdir.delete();
		         if(dir1!=null){
		          System.err.println(dir1.getAbsolutePath());
		          File[] gc = dir1.listFiles(new FilenameFilter(){

					public boolean accept(File dir, String name) {
						return name.startsWith("gcplot_tmp");
					}
		        	  
		          });
		          for(int i=0; i<gc.length; i++){
		        	  gc[i].delete();
		          }
		         }
		        // }
		        
		         es.shutdown();
		        // }       
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}

	private static Map<String, List<String>> readPlateFile(File plate2) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(plate2));
		String st = br.readLine();
		List<String> head = Arrays.asList(st.split("\t"));
		int id_ind = head.indexOf("PATIENT");
		int plate_ind = head.indexOf("PLATE");
		Map<String, List<String>> m= new HashMap<String,List<String> >();
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			List<String >v = m.get(str[plate_ind]);
			if(v==null){
				m.put(str[plate_ind], v = new ArrayList<String>());
			}
			v.add(str[id_ind]);
		}
		return m;
	}

	private static String[] getTodo(File file) {
		System.err.println("list for "+file.getAbsolutePath());
	final 	List<String> pref = OptionBuild.todo == null ? null : Arrays.asList(OptionBuild.todo);
		String[] f =  file.list(new FilenameFilter(){

				public boolean accept(File pathname, String name) {
					return name.endsWith("zip") && (pref == null || pref.contains(name.substring(0, name.length()-4)));
				}
	        	 
	         });
		for(int i=0; i<f.length; i++){
			f[i] = f[i].substring(0,f[i].length()-4);
		}
		System.err.println("dir: "+file.getName()+" to do "+Arrays.asList(f));
		return f;
	}
}

/* if(mode.contains("phenotype")){
			        	 BufferedReader br = new BufferedReader(new FileReader(new File(dir1, phenFile)));
			        	 ZipFile zf = new ZipFile(new File(dir[kk], todo[0]+".zip"));
			        	 List<String> samples = Compressor.readSamples(zf);
			        	 String[] order = new String[samples.size()];
			        	 String st = "";
			        	 String header_st = br.readLine();
			        	 String[] header = header_st.split("\t");
			        	 Set<String>[] s = new Set[header.length-1]; 
			        	 for(int k=1; k<header.length; k++){
			        		 s[k-1] = new TreeSet<String>(new Comparator<String>(){
	
								public int compare(String o1, String o2) {
									return -1*o1.compareTo(o2);
								}
			        			 
			        		 });
			        	 }
			        	 int id_index = 0;
	//		        	 int phen_index = 1;
			        	 while((st = br.readLine())!=null){
			        		st =  st.replaceAll("Unknown", "null");
			        		 String[] str = st.split("\t");
			        		 for(int k=1; k<header.length; k++){
				        		 s[k-1].add(str[k].trim());
				        	 }
			        		 int ind =samples.indexOf(str[id_index]);
			        		 if(ind>=0) order[ind] = st;
			        		 else{
			        			 System.err.println("warning - no sample in zip file for "+str[id_index]);
			        		 }
			        	 }
			        	 PrintWriter pw = new PrintWriter(new FileWriter(new File(outputDir[kk],"Samples.txt")));
			        	 PrintWriter pw1 = new PrintWriter(new FileWriter(new File(outputDir[kk],"include.txt")));
			        	  
			        	 for(int k=1; k<header.length; k++){
			        		 pw1.print(header[k]+"\t");
			        		for(Iterator<String> it = s[k-1].iterator();it.hasNext();){
			        			String next = it.next();
			        			if(next.equals("null")) continue;
			        			pw1.print(next+(it.hasNext() ? "\t":"\n"));
			        		}
			        	 }
			        	 pw1.close();
			        	 pw.println(header_st);
			        	 for(int i=0; i<order.length; i++){
			        		 if(order[i]==null) throw new RuntimeException("!!");
			        		 pw.println(order[i]);
			        	 }
			        	 pw.close();
			         }*/
