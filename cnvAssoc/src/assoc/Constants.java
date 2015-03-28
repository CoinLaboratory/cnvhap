package assoc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;


public class Constants {
	public static String[] experiment = null;
	public static String[] type = new String[] {"_countAll"};//,"_countB"};//, "_countAll"};
	public static String[] mid = null;
	public static String restrict  = null;
	public static String probe = "rs800729";
	public static String src = "2";
	public static String target = "2";
	public static String printNaN = "true";
	public static String covariates = "covar.txt";
	public static String distance = "10kb";
	public static String results = "QTL";
	public static String baseDir = ".";
	public static String pheno = "pheno.txt";
	public static String residuals = "true";
	public static String[] dir = null;
	public static String minCount = "0";
	public static String thresh = "0.5";
	public static String[] inclDir = new String[] {"."};
	
	
	public static String[] experiment(){
		if(experiment==null){
			experiment = new String[dir.length];
			for(int k=0; k<dir.length; k++){
				int ind = dir[k].lastIndexOf('/');
				
			experiment[k]  = ind >=0 ?dir[k].substring(ind).split("_1_")[0]: dir[k].split("_1_")[0];
			}
		}
		return experiment;
	}
	public static String mid(){
		return mid();
	}
	
	public static Options makeOptions(final Class clazz){
		 return 	new Options(){
				{
					Field[] f = clazz.getFields();
					for(int i=0; i<f.length; i++){
						if(Modifier.isStatic(f[i].getModifiers())){
							this.addOption(OptionBuilder.withLongOpt(f[i].getName()).withValueSeparator(':').hasArgs().create());
						}
					}
				}
			};
		
		}
	
	public static final Options OPTIONS = makeOptions(Constants.class);

	
	public static void parse(File params, Class clazz, int col, int row) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(params));
		String st = "";
		Map<String, Integer> cnts = new HashMap<String, Integer>();
		while((st=br.readLine())!=null){
			String st1 = st.trim();
			if(!st1.startsWith("#") && st1.length()>0 && st1.startsWith("--")){
				String[] str = st.substring(2).split("\\s+");
				Integer cnt = cnts.get(str[0]);
				if(cnt==null) cnt = 1;
				else cnt++;	
				cnts.put(str[0], cnt);
				if(cnt<=row){
					set(clazz, str[0], str[Math.min(col, str.length-1)]);
				}
			}
		}
		br.close();
	}
	
	public static String[] inclDir(){
		return Constants.inclDir;
	}
	
	   private static void set(Class clazz, String string, String string2) throws Exception {
		   Field f = clazz.getField(string);
			Class type_ = f.getType();
			System.err.println("set "+type_.getName()+" "+string+" "+string2+" "+f.getName());
			if(type_.isArray()){
				String[] str = string2.split(":");
				Class compt = type_.getComponentType();
				Object val = Array.newInstance(compt, str.length);
				
				for(int i=0; i<str.length; i++){
				//	if(compt.)
					
					Object value = compt.getConstructor(String.class).newInstance(new Object[] {str[i]});
					Array.set(val, i, value);
				}
				f.set(null, val);
			}
			
			else {
				try{
				f.set(null, type_.getConstructor(String.class).newInstance(new Object[] {string2}));
				}catch(Exception except){
					System.err.println("prob with "+string);
				}
			}
			
	}
	public static CommandLine parse(String[] args, Class clazz) throws Exception{
		Parser parser = new PosixParser();
		CommandLine cli  = parser.parse(OPTIONS, args);
		Option[] o = cli.getOptions();
		for(int i=0; i<o.length; i++){
			Field f = clazz.getField(o[i].getLongOpt());
			Class type = f.getType();
			System.err.println(f.getName());
			if(!type.isArray()){
				f.set(null, type.getConstructor(String.class).newInstance(o[i].getValue()));
			}
			else {
				
				Object arr = Array.newInstance(type.getComponentType(), o[i].getValues().length);
				for(int k=0; k< o[i].getValues().length; k++){
					
					Array.set(arr, k, type.getComponentType().getConstructor(String.class).newInstance(o[i].getValues()[k]));
				}
				f.set(null, arr);
			}
		}
		
		return cli;
	}
	   
	   public static int process(String txt){
	  		
	  		int ind = txt.indexOf("mb");
	  		if(ind>=0){
	  			return (int) Math.round(1000000*Double.parseDouble(txt.substring(0,ind)));
	  		}
	  		 ind = txt.indexOf("kb");
	  		if(ind>=0){
	  			return (int) Math.round(1000*Double.parseDouble(txt.substring(0,ind)));
	  		}
	  		 ind = txt.indexOf("gb");
	  		if(ind>=0){
	  			return (int) Math.round(1000000000*Double.parseDouble(txt.substring(0,ind)));
	  		}
	  		
	  		else return (int) Math.round(Double.parseDouble(txt));
	  	}

	public static boolean residuals() {
	return residuals.startsWith("t");
	}
	public static String[] dir(File user) {
	
		if(dir==null){
		dir = new String[Constants.experiment.length];	
		for(int j=0; j<dir.length; j++){
		dir[j] = Constants.experiment[j]+"_1_"+Constants.mid[0];
		for(int k=1; k<Constants.mid.length; k++){
			dir[j]+="_"+
			//Constants.process(
					Constants.mid[k]
					              //)
					              ;
		}
		}
		}
		else{
			for(int j=0; j<dir.length; j++){
				final String str = dir[j];
				File[] fs = 
				 user.listFiles(new FileFilter(){

					@Override
					public boolean accept(File pathname) {
						// TODO Auto-generated method stub
						return pathname.getName().indexOf(str)>=0;
					}
					
				});
				dir[j] = fs[0].getName();
				File f = new File(user,dir[j]);
				
			}
		}
		return dir;
	}
	
	public static String chrom() {
		if(mid!=null) return mid[0];
		else return Constants.dir[0].split("_")[2];
	}
	public static Boolean printNa=null;
	public static boolean printNaN() {
		if(printNa==null) printNa = printNaN.toLowerCase().startsWith("t");
	return printNa;
	}
	public static int[][] cod = null;
	public static String[] coding = null;
	public static String mult = "false";
	public static String step = "true";
	public static int[][] coding(){
		if(cod==null || coding!=null){
			String[] str = coding;
			if(str!=null){
			 cod = new int[str.length][];
			for(int k=0; k<cod.length; k++){
				String[] str1 = str[k].split(";");
				cod[k] = new int[str1.length];
				for(int k1=0; k1<str1.length; k1++)
				{
				cod[k][k1] = Integer.parseInt(str1[k1]);
				}
			}
			}
		}
		return cod;
	}
	public static boolean mult() {
		return mult.toLowerCase().startsWith("t");
	}
	public static double maf_thresh() {
		// TODO Auto-generated method stub
		return maf_thresh;
	}
	public static Double maf_thresh=1e-5;
	public static String[] snpsToInclude;
	public static List<String> snpsToInclude(){
		if(snpsToInclude==null) return null;
		List l = new ArrayList<String>(Arrays.asList(snpsToInclude));
		//l.add("snpid");
		return l;
	}
	
	public static List<String> snpsToExclude(){
		if(snpsToExclude==null) return null;
		List l = new ArrayList<String>(Arrays.asList(snpsToExclude));
		return l;
	}
	public static String[] snpsToExclude;

	public static Integer buffer = -1;
    public static Integer limit = Integer.MAX_VALUE;


    
    public static Integer buffer(){ return buffer;}
	public static void write(File in, File out) throws Exception{
		// TODO Auto-generated method stub
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
		BufferedReader br = new BufferedReader(new FileReader(in));
		String st = "";
		while((st = br.readLine())!=null){
			pw.println(st);
		}
		pw.close();br.close();
	}
	public static Double round = 0.05;
	public static double fisherThresh = 5; //min counts when we do chisq instead of fisher
	public static double round() {
		// TODO Auto-generated method stub
		return round;
	}
	public static File getFile(File dir, final String toincl){
		try{
		File[] f = dir.listFiles(new FileFilter(){
	
			@Override
			public boolean accept(File f) {
				return f.getName().indexOf(toincl)>=0;
			}
			
		});
		if(f.length!=1) {
			
			throw new RuntimeException("!!");
		}
		return f[0];
		}catch(Exception exc){
			
			System.err.println(Arrays.asList(dir.list()));
			System.err.println(toincl);
			exc.printStackTrace();
			return null;
		}
	}
	public static String getCombName(String pref,File[] user1, String string) {
		StringBuffer sb = new StringBuffer(pref);
	 	for(int k=0; k<user1.length; k++){
	 			
	 		sb.append(user1[k].getName().split("-")[0]);
	 		if(k<user1.length-1) sb.append("-");
	 	}
		return sb.toString();
	}
	public static int sum(int[] sze) {
		int n = 0;
		for(int k=0; k<sze.length; k++) n+=sze[k];
		return n;
	}
	
}
