package lc1.assoc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.List;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.IlluminaRDataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.PhasedDataState;
import lc1.stats.IlluminaDistribution;
import lc1.util.Constants;

import org.biojava.utils.ProcessTools;

public class IntensityAssoc {

	public static void main(String[] args){
		try{
			extractToTempFile();
			runassoc();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	private static void runassoc() throws Exception{
		File dir = new File(System.getProperty("user.dir"));
		 Writer err = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "err.txt"))));
		   File tmp = new File(dir, "assoc.R");
	     String[] exec = new String[] {"R", "--vanilla", "<", tmp.getAbsolutePath()};
         String ex_str = Arrays.asList(exec).toString().replaceAll(",", "");
         ex_str = ex_str.substring(1, ex_str.length()-1);
         System.err.println(ex_str);
         String[] exec1 = new String[] {"R", "--vanilla"};
         BufferedReader br1 = new BufferedReader(new FileReader(tmp));
      ProcessTools.exec(exec1, br1, err, err);
		
	}

	static String chr = "11";
	static int start = 27345080-1;
	static int end =  27723114+1;
	static String pref = "";
	static void extractToTempFile() throws Exception{
		   Constants.plot = 2;
	        Constants.var_thresh= new double[] {0.1};
	        Constants.numThreads = OptionBuild.numThreads;
	      //  Constants.topBottom = false;
	        Constants.loess = new boolean [] {true};
	        Constants.median_correction = new boolean [] {true};
			File user = new File(System.getProperty("user.dir"));
			String suff = ".txt";//+System.currentTimeMillis();
			File r = new File(user, "r"+suff);
		//	File b = new File(user, "b"+suff);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(r)));
			//PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(b)));
			File[] f =  user.listFiles(new FileFilter(){
	
				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().startsWith(pref);
				}
				
			});
			int[] mid = new int[] {start, end};
			DataCollection[] dc = new DataCollection[f.length];
			for(int i=0; i<f.length; i++){
				dc[i] = null;//	new IlluminaRDataCollection(new File(f[i],chr+".zip"),(short) 0, 2,
						//new int[][]{mid},  null, null);
			}
			boolean genotype =true;
			List<String> snps = dc[0].snpid();
				pw.print("id");//pw1.print("indiv");
				List<String> phen = dc[0].pheno.phen;
				for(int k1=0; k1<phen.size(); k1++){
					pw.print("\t"+phen.get(k1));
				}
				for(int i=0; i<snps.size(); i++){
					if(genotype) pw.print("\t"+dc[0].snpid().get(i));
					else{
						pw.print("\t"+dc[0].snpid().get(i)+"_r");
						pw.print("\t"+dc[0].snpid().get(i)+"_b");
					}
				}
				pw.println();
				pw.print("id");//pw1.print("indiv");
				//List<String> phen = dc[0].pheno.phen;
				for(int k1=0; k1<phen.size(); k1++){
					pw.print("\t"+phen.get(k1));
				}
				for(int i=0; i<snps.size(); i++){
					if(genotype){
						pw.print("\t"+dc[0].loc.get(i));
					}
					else{
						pw.print("\t"+dc[0].loc.get(i));
						pw.print("\t"+dc[0].loc().get(i));
					}
				}
				pw.println();
				for(int k=0; k<dc.length; k++){
					for(int j=0; j<dc[k].indiv().size(); j++){
						String key =dc[k].indiv().get(j); 
						pw.print(key);
						//pw1.print(key);
						IlluminaNoBg hes = (IlluminaNoBg) dc[k].dataL.get(key);
						PhasedDataState hes1 = (PhasedDataState) dc[k].data.get(key);
						Double[] phenv = hes.phenValue();
						for(int k1=0; k1<phenv.length; k1++){
							pw.print("\t"+(phenv[k1]==null ? "NA" : phenv[k1]));
						}
						for(int i=0; i<snps.size(); i++){
							Number  r_ = ((IlluminaDistribution) hes.emissions[i]).r();
							Double b_ = ((IlluminaDistribution) hes.emissions[i]).b();
							if(genotype){
								int kk = hes1.emissions[i].fixedInteger();
								ComparableArray comp = (ComparableArray)hes1.getEmissionStateSpace().get(kk);
								char[] ch = (comp.get(0).toString()+comp.get(1).toString()).toCharArray();
								Arrays.sort(ch);
								String st = new String(ch);
								if(st.indexOf('_')>=0){
									pw.print("\tNA");
								}
								else if(st.equals("AA")){
									pw.print("\t0");
								}
								else if(st.equals("AB")){
									pw.print("\t1");
								}
								else if(st.equals("BB")){
									pw.print("\t2");
								}
								else throw new RuntimeException(st);
							}
							else{
								pw.print("\t"+(r_==null ? "NA" : r_));
								pw.print("\t"+(b_==null ? "NA" : b_));
							}
							
						//	pw.print(dc[k].
						}
						pw.println();
						//pw1.println();
					}
				}
				pw.close(); //pw1.close();
		
		
	}
}
