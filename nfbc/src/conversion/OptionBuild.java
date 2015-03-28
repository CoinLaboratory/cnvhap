package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class OptionBuild {
	static{
		try{
		   javax.swing.UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
		}catch(Exception exc){
			//exc.printStackTrace();
		}
	}
	 static JFileChooser filechooser = new JFileChooser();
	 static File user =  new File( System.getProperty("user.dir"));
	public static  Boolean orderedByChrom = true;
	 public static String[] todo = null;//"1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:XY".split(":");
	 public static File[] file;
	 public static String build = "build36";
	 public static String[] mode = new String[] {"compress", "build", "sep", "split",  "qc", "join"};
	 public static String runR = "true";
	 public static File[] output = new File[] {new File(user, "output_dir")};  //dir with zip files after compression, but before qc
	/* static{
		 for(int i=0; i<output.length; i++){
		if(!output[i].exists()){
			output[i].mkdir();
		}
		 }
	 }*/
	 public static Integer numThreads = 3;
	 public static String splStr = "\t";//"\t"
	 static String phenFile = "phenotype.txt";//params.getOptionValue("phenFile", "phenotype.txt");
	public static String feat = "Log R";
	public static String plateId = "PLATE";
	public static String sampleId = "PATIENT";
	public static String[] sampleId1 = "SentrixBarcode:SentrixPosition_A".split(":");
	public static File plateFile = new File(user, "plate.txt");
	    static final Options OPTIONS  = new Options(){
	         {
	             Field[] f = OptionBuild.class.getFields();
	             for(int i=0; i<f.length; i++){
	                 if(Modifier.isStatic(f[i].getModifiers()) && Modifier.isPublic(f[i].getModifiers())){
	                	 System.err.println("field "+f[i].getName());
	                     this.addOption( OptionBuilder.withLongOpt( f[i].getName() ).withDescription( f[i].getName()).withValueSeparator( ':' ).hasArgs().create());
	                 }
	                 else{
	                     System.err.println("excluded "+f[i].getName());
	                 }
	             }
	         }
	     };
		public static final String SNPS = "SNPS.txt";
		public static  double varThresh = 0.5;
	     
	     public static void getVal(String stri, String header1, CommandLine params) throws Exception{
	    	 String header = stri+":"+header1;
	    	 Field f_i = OptionBuild.class.getField(stri);
	    	 if(params.hasOption(stri)) return;
	    	 if(Modifier.isStatic(f_i.getModifiers())){
            	 Class clazz = f_i.getType();
            	 Object valo =  f_i.get(null);
            	 if(valo instanceof String[]){
            		 valo = Arrays.asList((String[])valo).toString().replaceAll("\\s+","").replace(',', ':');
            		 valo = ((String)valo).substring(1,((String)valo).length()-1);
            	 }
            	 String nme = f_i.toString().substring(f_i.toString().lastIndexOf('.')+1);
            	 if(nme.equals("dir") || nme.equals("plateFile")){
            		 filechooser.setDialogTitle(header);
            		 filechooser.setCurrentDirectory(user);
            		 filechooser.setMultiSelectionEnabled(true);
            		 filechooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
            		 filechooser.showOpenDialog(null);
            		 filechooser.setName(header);
            		 File[] f = filechooser.getSelectedFiles();
            		
            		 if(nme.equals("dir")) f_i.set(null, f);
            		 else f_i.set(null, f[0]);
            	 }
            	 else if(nme.equals("output")){
            		 filechooser.setDialogTitle(header);
            		 filechooser.setCurrentDirectory(user);
            		 filechooser.setMultiSelectionEnabled(true);
            		 filechooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            		 filechooser.showOpenDialog(null);
            		 filechooser.setName(header);
            		 File[] f = filechooser.getSelectedFiles();
            		 f_i.set(null, f);
            	 }
            	 else{
	            	 String val =(String) JOptionPane.showInputDialog(null,
	                         header, nme, JOptionPane.QUESTION_MESSAGE, null, null,
	                valo);
	            	if(val==null) return;
	            	 Object val1 ;
	            	 if(clazz.equals(String[].class)){
	            		 val1 = val.split(":");
	            	 }
	            	 else val1= clazz.getConstructor(new Class[] {String.class}).newInstance(new Object[] {val});
	            	 f_i.set(null, val1);
            	 }
             }
             else{
                 System.err.println("excluded "+f_i);
             }
	     }
	     
	    public static String format = "long"; //can be long or wide
	     public static void getVals(CommandLine params) throws Exception{
	    	 getVal("mode","What steps do you want to run", params );
	    	 List<String> modeL = Arrays.asList(mode);
	    	 if(modeL.contains("compress")){
	    		 getVal("file", "data files", params);
	    		 user = file[0].getParentFile();
	    		 System.err.println("user is "+file[0].getName());
	    		 getVal("output","Name of output directory for zip files", params );
	    		 getVal("splStr","What is the separator", params );
	    	 }
	    	 else{
	    		 
	    		 getVal("output","Name of directories with .zip files", params );
//	    		 getVal("dir", "Name of directory with compressed zip files");
	    		
	    		 user = output[0].getParentFile();
	    		
	    	 }
	    	 if(modeL.contains("sep")){
	    		 getVal("plateFile", "The plate file", params);
	    		 if(plateFile!=null){
	    			 getVal("sampleId", "What is the header in the file for sample", params);
	    			 getVal("plateId", "What is the header in the file for plate", params);
	    		 }
	    		
	    		
	    	 }
	    	// getVal("todo", "Chromosomes to do, separated by ':'", params);
	    	 getVal("numThreads", "how many cpus you want to run this on", params);
	    	 getVal("build", " Which build ", params);
	    	
	    	 getVal("feat","Name feature for loess correction" , params);
	    	 getVal("varThresh","variance threshold" , params);
	    	// getVal("runR","What steps do you want to run" );
	         
	     }
	     private static void delete(File file) {
			if(file.exists()){
				if(file.isDirectory()){
					File [] f = file.listFiles();
					for(int i=0; i<f.length; i++){
						delete(f[i]);
					}
				}
				file.delete();
			}
			
		}

		private static String[] getHeader(String string) {
			File f = new File(user, string);
			try{
				BufferedReader br = new BufferedReader(new FileReader(f));
				String[] res =  br.readLine().split("\t");
				br.close();
				return res;
			}catch(Exception exc){
				exc.printStackTrace(); 
				return null;
			}
		} 
		private static boolean useFileChooser = true;
		private static boolean hardCoded = false;
		public static String default_type = "Log R";
		public static String prefix = "";
		public static String rlib = null;
		
		 static void setVals(CommandLine params) throws Exception{
    		 Field[] f = OptionBuild.class.getFields();
             for(int i=0; i<f.length; i++){
            	 if(Modifier.isStatic(f[i].getModifiers()) && Modifier.isPublic(f[i].getModifiers())){
            		 if(params.hasOption(f[i].getName())){
            			 if(f[i].getType().isArray()){
            				 Class inn = f[i].getType().getComponentType();
            				 String[] val = params.getOptionValues(f[i].getName());
            				 Object newObj = Array.newInstance(inn, val.length);
            				 System.err.println("Setting  "+f[i].getName()+" "+Arrays.asList(val));
            				 for(int k=0; k<val.length; k++){
            					 if(val[k].equals("null")) {
            						 Array.set(newObj, k,null);
            					 }
            					 else 
            					 Array.set(newObj, k, inn.getConstructor(new Class[] {String.class}).newInstance(
	            						 new Object[] {val[k]}));
            				 }
            				 f[i].set(null, newObj);
            				 
            			 }
            			 
            			 else{
            				 String val = params.getOptionValue(f[i].getName());
            				 if(val.equals("null")) f[i].set(null, null);
            				 else{
	            				 Object newObj = f[i].getType().getConstructor(new Class[] {String.class}).newInstance(
	            						 new Object[] {val});
	            				 f[i].set(null, newObj);
            				 }
            				 System.err.println("Setting "+f[i].getName()+" "+val);
            			 }
            		 }
            	 }
             }
			 
	          
	        // user = new File(dir);
	         
	         
		
		
     }
		private static void printParams(PrintWriter pw2) throws Exception {
	    	 Field[] f = ConvertIllumina.class.getFields();
	         for(int i=0; i<f.length; i++){
	             if(Modifier.isStatic(f[i].getModifiers())){
	            	 Object val = f[i].get(null);
	            	 if(val instanceof String[]){
	            		 pw2.println("field "+f[i].getName()+" "+Arrays.asList((String[])f[i].get(null)));
	            	 }
	            	 else
	            	 pw2.println("field "+f[i].getName()+" "+f[i].get(null));
	             }
	             
	         }
	 		
	 	}

		public static void getOptions(String[] args) throws Exception {
			Parser parser= new PosixParser();
			 final CommandLine params = parser.parse(OPTIONS, args, false);
			 setVals(params);
			
	     
	          //  getVals(params);
	          
	       
			
		}
}
