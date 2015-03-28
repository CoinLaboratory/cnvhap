package lc1.assoc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Arrays;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import lc1.dp.data.collection.Phenotypes;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class OptionBuild {
	
	 static Boolean hardCoded = false;
	public static Boolean hardCoded(){
		return hardCoded;
	}
	 static JFileChooser filechooser;
	 static File user;
	 public static Boolean probs = true;//
	 public static Integer contextsnps = 5;
	 public static Boolean recalculate = false;  // calculate or show pvalues
	 public static Boolean plot =true;  //show plot
	 public static Boolean append = false;
	 public static Boolean deletePrevious = false;
	 public static Double plotUpToPValue = 1.0; //show all points
	 public static String[] chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:XY".split(":");
	 public static String dir = System.getProperty("user.dir");
	 public static String[] input = new String [0];
	 public static    String[] type = new String[] {"armitage", "chisq"};//regrP"};//,"regrP"};
	public static String[] experiment = new String[]{"all"};// {"assoc-all.txt","assoc-igt.txt"};
	 
	 public static String[] excl = new String[0];
     public static String[] exclList = new String[0];
     //public static String []prefix ="1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y:XY".split(":");
     public static String[] tag_pheno =   new String[] {"regress","regrP" };
    
    // public static String experiment = "all1";
     public static Boolean regr=false;
   //  public static String dir = System.getProperty("user.dir");
	public static Integer numThreads = 3;
	public static Boolean doqq = true;
	 private static boolean useFileChooser = false;
	public static Boolean plotGenes = true;;
	 public static Double plotUpToPValue(){
		 return plotUpToPValue;
	 }
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
	     
	     public static File user(){
	    	 return user;
	     }
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
				 
		          
		         user = new File(dir);
		         
		         
			
			
	     }
	     public static String[] experiment(){
	    	 return experiment;
	     }
	     public static void getVal(String stri, String head, CommandLine params) throws Exception{
	    	 String header = stri+": "+head;
	    	 if(params.hasOption(stri)) return;
	    	 Field f_i = OptionBuild.class.getField(stri);
	    	 if(Modifier.isStatic(f_i.getModifiers())){
            	 Class clazz = f_i.getType();
            	 Object valo =  f_i.get(null);
            	 if(valo instanceof String[]){
            		 valo = Arrays.asList((String[])valo).toString().replaceAll("\\s+","").replace(',', ':');
            		 valo = ((String)valo).substring(1,((String)valo).length()-1);
            	 }
            	 String nme = f_i.toString().substring(f_i.toString().lastIndexOf('.')+1);
            	 if(nme.equals("dir")){
            		 String val1;
            		 if(hardCoded){
            			File f = new File((String)valo);
            			File[] fs = f.listFiles(new java.io.FileFilter(){

							public boolean accept(File pathname) {
								return pathname.getName().endsWith("zip");
							}
            				
            			});
            			if(fs.length>0) val1 = fs[0].getAbsolutePath();
            			 val1 =(String) JOptionPane.showInputDialog(null,
    	                         header, nme, JOptionPane.QUESTION_MESSAGE, null, null,
    	                valo);
            		 }
            		 else{
            		if(useFileChooser){
	            		 if(filechooser==null)  filechooser= new JFileChooser();
	            		 filechooser.setCurrentDirectory(new File((String)valo));
	            		 filechooser.setMultiSelectionEnabled(false);
	            		 filechooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
	            		 filechooser.showOpenDialog(null);
	            		 filechooser.setName(header);
	            		  val1 = filechooser.getSelectedFile().getAbsolutePath();
	            		 }
	            		 
	            		 else{
	            			 val1 =(String) JOptionPane.showInputDialog(null,
	    	                         header, nme, JOptionPane.QUESTION_MESSAGE, null, null,
	    	                valo);
	            		 }
            		 }
            		 f_i.set(null, val1);
            		 
            	 }	
            	 else if(nme.equals("exclList")){
            		 if(useFileChooser){
	            		 if(filechooser==null)  filechooser= new JFileChooser();
	            		 filechooser.setDialogTitle("Choose files with list of samples to exclude");
	            		 filechooser.setCurrentDirectory(user);
	            		 filechooser.setMultiSelectionEnabled(true);
	            		 filechooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
	            		 filechooser.showOpenDialog(null);
	            		 filechooser.setName(header);
	            		 File[] f = filechooser.getSelectedFiles();
	            		 String[] str = new String[f.length];
	            		 for(int i=0; i<f.length; i++){
	            			 str[i] = f[i].getName();
	            		 }
	            		 f_i.set(null, str);
            		 }
            		 else{
            			 String str =(String) JOptionPane.showInputDialog(null,
    	                         header, nme, JOptionPane.QUESTION_MESSAGE, null, null,
    	                valo);
            			 f_i.set(null, new String[] {str});
            		 }
            	 }
            	 else{
	            	 String val =(String) JOptionPane.showInputDialog(null,
	                         header, nme, JOptionPane.QUESTION_MESSAGE, null, null,
	                valo);
	            	
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
	     
	     public static void getVals(CommandLine p) throws Exception{
	    	
	    	 getVal("dir", "Results directory",p);
	    	//  getVal("chr", "Chromosomes to do, separated by ':'",p);
	    	 if(!hardCoded) {
	    		 getVal("recalculate", "true if you want to re-calculate the pvalues, false just to display previously calculated values",p);
	    	 }
	    	 getVal("numThreads", "how many cpus you want to run this on",p);
	    	 if(recalculate){
	    		 input = (new File(dir)).list(new FilenameFilter(){

					public boolean accept(File dir, String name) {
						File f1 = new File(dir, name);
						return 
						
						(input.length==0 || input[0].equals("*") || name.startsWith(input[0])) && 
						!name.equals("assoc") && (f1).isDirectory()
						 && (f1.list(new FilenameFilter(){

							public boolean accept(File dir, String name) {
								return name.endsWith("zip");
							}
							 
						 })).length>0;
					}
	    			 
	    		 });
	    		 getVal("input", "Directory(s) of cnv predictions together with Sample files", p);
	    		
	    		 getVal("experiment", "Names of experiments, separated by ':'",p);
	    		// ("excl");
	    		 
	    		 boolean resExists = false;
	    		 File assoc = new File(user, "assoc");
	    		 if(assoc.exists()){
		    		 inner: for(int i=0; i<experiment.length; i++){
		    			 File f_i = new File(assoc, experiment[i]);
		    			 if(f_i.exists()){
			    			 String[] fs = f_i.list();
							  if(fs!=null &&  (fs.length>0)){
								  resExists =true;
								  break inner;
							  }
		    			 }
					  }
		    		 if(resExists){
			    		 getVal("deletePrevious","do you want to delete all previous association results?",p);
			    		 if(deletePrevious){
			    			 deletePrevious = false;
			    			 getVal("deletePrevious","Confirmation: delete all previous association results?",p);
			    			 if(deletePrevious){
			    				
			    				  for(int i=0; i<experiment.length; i++){
			    					  delete(new File(user, experiment[i]));
			    				  }
			    			 }
			    		 }
		    		 }
	    		 }
	    		 
	    		 if(!p.hasOption("excl")){
	    			 excl = new String[experiment.length];
	    			 exclList = new String[experiment.length];
	    		 for(int k=0; k<experiment.length; k++){
	    			 File include = new File(new File(dir), "include.txt");
	    			 if(!include.exists()) throw new RuntimeException("include file does not exist "+Arrays.asList(user.list()));
	    			
		    		 Phenotypes pheno = new Phenotypes(include);
		    		 StringBuffer sb = new StringBuffer();
		    		
		    		
		    		 for(int i=0; i<pheno.size(); i++){
		    			String phe =  pheno.phen.get(i);
			    			if(pheno.type()[i]!=0){
			    				String excval =(String)  JOptionPane.showInputDialog(null, 
			    						"Experiment "+k+", "+experiment[k]+": exclude based on phenotype: "+phe+" ? ",
			    						"Exclusions",
			    						JOptionPane.QUESTION_MESSAGE,
			    						null, 
			    						pheno.phenVals[i].keySet().toArray(new String[0]), null);
			    				if(excval!=null) sb.append(phe+","+excval+":");
			    				
			    				
			    				
			    			}
			    			else{
			    				String excval =(String)  JOptionPane.showInputDialog(null, 
			    						"Exclude based on phenotype: "+phe+". Choose with vals less than? ",
			    						"Exclusions",
			    						JOptionPane.QUESTION_MESSAGE,
			    						null, 
			    						null, null);
			    				if(excval!=null) sb.append(phe+","+excval+";");
			    			}
			    		 }
			    		 String str = sb.toString();
			    		 if(str.length()>0) str = str.substring(0, str.length()-1);
		    		 
		    		 excl[k] = str;
		    		 exclList[k] = (String)  JOptionPane.showInputDialog(null, 
	    						"Experiment "+k+", "+experiment[k]+": exclude list:  ? ",
	    						"Exclusions",
	    						JOptionPane.QUESTION_MESSAGE,
	    						null, 
	    						null, null);
		    		 if(k>0 && excl[k].equals(excl[0])) throw new RuntimeException("experiments are the same!");
		    		// experiment = str.length()==0 ? new String[] {"all"} : new String[] {str.replace(':', ';')};
		    		 inner: for(int j=0; j<chr.length; j++){
		    			 File out = new File(user, "assoc-"+experiment[k]+"/"+chr[j]);
		    			 if(out.exists()){
		    				 getVal("append", "append (true) or overwrite(false)",p);
		    				 break inner;
		    			 }
		    		 }
	    		 }
	    		 }
	    		 getVal("regr", "true to do regression, false to do armitage and chisq",p);
	    		 if(regr.equals(false)){
		        	   tag_pheno = "chisq_state:armitage_state:odds_state:cases_state:controls_state".split(":");
		           }
		           else{
		        	   tag_pheno =  new String[] {"regress","regrP" };
		           }
	    		 
	    		 
	    	 }
	    	 else{
	    		 getVal("plot", "plot the results, or print significant results to file",p);
	    		
	    		if(!hardCoded) getVal("probs", "plot probs or chisq values",p);
	    		 getVal("plotUpToPValue", 
	    				 OptionBuild.probs ? 
	    		 "p value threshold for plotting or reporting results.  Set < 1e-7 for reporting results":
	    			 "chisq threshold for plotting or reporting results.  Set >28  for reporting results",p	
	    		 
	    		 );
	    		 getVal("doqq", "plot QQ plot as well",p);
	    		if(!hardCoded) getVal("type", "type of results to report, separated by ':', can be armitage, chisq, regrP",p);
	    		else type = new String[] {"P"};
	    		String chrom = chr[0];
	    		System.err.println(user.getAbsolutePath());
	    		File ass = new File(user, "assoc");
	    		if(ass.exists()){
	    			experiment = ass.list();
	    		}
	    		else  experiment =  (new File(user, chrom)).list();
	    		if(!hardCoded){ getVal("experiment", "name of experiments, separated by ':'",p);
	    		 getVal("contextsnps", "context for reporting haplotype structures",p);
	    		}
	    	 }
	        
	    	
	         
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

		private static void printParams(PrintWriter pw2) throws Exception {
	    	 Field[] f = CalcAssociation1.class.getFields();
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
			
	            //Constants.modelCNP = 10;
	        //    System.err.println("prefix "+chr);
	            getVals(params);
	            if(recalculate){
		            final PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(user, "log_"+experiment[0]+".txt"))));
		            printParams(pw);
		            pw.close();
	            }
	             user = new File(dir);
	          //   List<String> t=  Arrays.asList(type);
	          //   scoreChi = t.contains("armitage") || t.contains("chiSq");
	          //   scoreRegression = t.contains("regrP");
	            /*if(experiment.length==0 || experiment[0].length()==0) experiment = new String[] {".txt"};
	            final  List<String> exp = Arrays.asList(experiment);
	             {
	              	experiment = (new File(user, chr[0])).list(new FilenameFilter(){
	  					public boolean accept(File dir, String name) {
	  						for(int i=0; i<exp.size(); i++){
	  							if(name.indexOf(exp.get(i))>=0) return true;
	  						}
	  						return false;
	  					}
	              		
	              	});
	              }*/
			
		}

		//private static boolean scoreRegression = false;
		//private static boolean scoreChi = false;
		public static boolean scoreRegression() {
			return regr;
		}
		public static boolean scoreChi() {
			return !regr;
		}
		public static String[] input() {
			return input;
		}
}
