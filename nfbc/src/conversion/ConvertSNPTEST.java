package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class ConvertSNPTEST {
	static int index = 0;
	int len1 = 0; // set to 3 if it is chr3 rather than just 3
	static String[] prefix = new String[] {"WED11_female_DIAS","withCOUNTS_femaleSYS", "WED11_male_DIAS", "withCOUNTS_maleSYS", "SUNf", "SUNm", "BMI"};
	static String[] pref = new String[] {"DIAS", "SYS","DIAS", "SYS", "", "", "BMI" };
	static int[] index1 = new int[] {0,0,0,0,1,1, 0};
	static String[] f = new String[] {"femaleDIAS", "withCounts_female_SYS", "maleDIAS", "withCounts_male_SYS", "HTN_female", "HTN_male", "BMI"};
	final int prefixLength;
	public static void main(String[] args){
		try{
			for(int i=0; i<args.length; i++){
			index = Integer.parseInt(args[i]);
			  File dir1 = new File(System.getProperties().getProperty("user.dir"));
			//for(int i=0; i<todo.length; i++){
				ConvertSNPTEST cst = new ConvertSNPTEST(new File(dir1, f[index]));
				cst.run();
			}
			//}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	//public static String[] todo = 
	
	//new String[] {"males_dystolic", "females_dystolic", "males_systolic", "females_systolic"};
	//new String[] {"_females", "_males"};
	//String stri;
	
	
	File[] in;
	
	
	
	PrintWriter pw;
	int num_ = 1;
	ConvertSNPTEST(File dir) throws Exception{
		this.prefixLength =   prefix[index].length();
		setDet();
		final File outF = new File(dir, "out.txt");
		in = dir.listFiles(new FileFilter(){
			public boolean accept(File pathname) {
				return  !outF.equals(pathname);
			}
			
		});
		pw = new PrintWriter(new BufferedWriter(new FileWriter(outF)));
	//	int pl = prefixLength;
		Arrays.sort(in, new Comparator<File>(){

			public int compare(File o1, File o2) {
				
				String nme1 = o1.getName();
				String nme2 = o2.getName();
				
				try{
				Integer i1 = Integer.parseInt(split(nme1));
				Integer i2 = Integer.parseInt(split(nme2));
				return i1.compareTo(i2);
				}catch(Exception exc){
					System.err.println("prob "+nme1+" "+nme2);
					exc.printStackTrace();
					System.exit(0);
					return 0;
				}
			}
		});
		
		
	}
	
	List<String> snptest;
	
	String chrom = "";
	 String split(String str){
		String[] st1 = str.substring(prefixLength).split("_");
		String[] st11 = st1[st1.length-1].substring(len1).split("\\.");
		return st11[0];
	}
	public void run() throws Exception{
		this.print(det[index1[index]].outH);
		for(int i =0; i<in.length; i++){
			chrom = split(in[i].getName());//.substring(prefixLength).split("_")[0].substring(len1).split("\\.")[0];
			BufferedReader br = new BufferedReader(new FileReader(in[i]));
			snptest = Arrays.asList(br.readLine().split("\\s+"));
			alias = new int[det[index1[index]].outH.length];
			for(int k=0; k<det[index1[index]].outH.length; k++){
				alias[k] = det[index1[index]].equiv[k]==null ? -1 : snptest.indexOf(det[index1[index]].equiv[k]);
			}
			String str="";
			while((str = br.readLine())!=null){
				this.convert(str);
				this.print(toPrint);
			}
		}
		pw.close();
	}
	public void convert(String str){
		String[] stri = str.split("\\s+");
		for(int i=0; i<alias.length; i++){
			if(alias[i]>=0){
				toPrint[i] = stri[alias[i]];
			}
		}
		toPrint[1] = chrom;
		toPrint[5] = "+";
		toPrint[toPrint.length-1] = stri[0].equals("---") ? "N" : "Y";
	}
	
	static String sep = ",";
	
	static class Details{
		String[] outH;
		String[] equiv;
		
	}
	
	static Details[] det = new Details[] {new Details(), new Details()}; //first is qtl, second cc
	public void setDet(){
	  det[0].outH =  new String[] {
			"SNPID", "chr", "position", "coded_all",
			"noncoded_all", "strand_genome", "beta", "SE", 
			"pval", "obs_var_imp", "info", 
			"genotype count coded-coded",
			"genotype count coded-noncoded",
			"genotype count noncoded-noncoded",
			"exp genotype count coded-coded",
			"exp genotype count coded-noncoded ",
			"exp genotype count noncoded-noncoded",
			"directly genotyped"
	};
	  det[0].equiv = 
			new String[] {
				"rsid", null, "pos", "allele_B",
				"allele_A", null, pref[index]+"_frequentist_add_proper_beta_1", pref[index]+"_frequentist_add_proper_se_1", 
				pref[index]+"_frequentist_add_proper", "info", pref[index]+"_frequentist_add_proper_info", 
				"all_BB",
				"all_AB",
				"all_AA",
				"all_BB_exp",
				"all_AB_exp",
				"all_AA_exp",
				null
		};
		det[1].outH = new String[] {
			"SNPID", "chr", "position", "coded_all",
			"noncoded_all", "strand_genome", "beta", "SE", 
			"pval", "obs_var_imp", "info", 
			"control genotype count coded-coded",
			"control genotype count coded-noncoded",
			"control genotype count noncoded-noncoded",
			"case genotype count coded-coded",
			"case genotype count coded-noncoded",
			"case genotype count noncoded-noncoded",
			"directly genotyped"
		};
	
	det[1].equiv = 
		new String[] {
			"rsid", null, "pos", "allele_B",
			"allele_A", null, "frequentist_add_cov_all_proper_beta_1", "frequentist_add_cov_all_proper_se_1", 
			"frequentist_add_cov_all_proper", "info", "frequentist_add_cov_all_proper_info", 
			"controls_BB",
			"controls_AB",
			"controls_AA",
			"cases_BB",
			"cases_AB",
			"cases_AA",
			null
	};
		toPrint =  new String[det[index1[index]].outH.length];
	}
	String[] toPrint; 
	int[] alias;
	public void print(String[] str){
		pw.print(str[0]);
		for(int i=1; i<str.length; i++){
			pw.print(sep+str[i]);
		}
		pw.println();
	}
	
}
