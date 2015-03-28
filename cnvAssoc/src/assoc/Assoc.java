package assoc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RList;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

import cnvtrans.LoopCall;

import cern.colt.matrix.linalg.Algebra;

public class Assoc {

	public static void main() {
		try {
			String path = System.getProperty("java.library.path");
			System.err.println("java.library.path = " + path);
		
			String[] inclDir = Constants.inclDir();
			String[] dir = new String[inclDir.length];
			for(int k=0; k<dir.length; k++){
				File user1 = new File(System.getProperty("user.dir")+"/"+inclDir[k]);
			  dir[k] = user1.getAbsolutePath()+"/"+Constants.dir(user1)[0];	
			}
			Assoc[][] assoc = new Assoc[Constants.type.length][dir.length];
List<String> phenos = new ArrayList<String>();
			for (int k = 0; k < Constants.type.length; k++) {

				for (int j = 0; j < dir.length; j++) {
					int ind = dir[j].lastIndexOf('/');
					if(ind<0) {
						dir[j] = "./"+dir[j];
						ind = dir[j].lastIndexOf('/');

					}
					
					File user =new File(dir[j].substring(0, ind));
					//+ "/"+ (ind < 0 ? "" : 
					                      //  );
					File pheno = new File(user, Constants.pheno);
					File limit = new File(user, "limit.txt");
					File exclude = new File(user, "exclude.txt");
					
					String suff = dir[j].substring(ind);
					if(!user.exists()) {
						throw new RuntimeException("!!");
					}
					File suff1 = new File(user,suff);
					if(!suff1.exists()) {
						throw new RuntimeException("!!");
					}
					File avgDir = new File(suff1, "avg");
					File[] l = avgDir.listFiles(new FileFilter() {

						@Override
						public boolean accept(File pathname) {
							boolean acc =  pathname.getName().endsWith(".txt")
									&& !pathname.getName()
											.startsWith("result_");
							return acc;
						}

					});
					File avg = l.length == 1 ? l[0] : new File(avgDir,
							Constants.experiment()[j] + ".txt");

					assoc[k][j] = new Assoc(avg, pheno, limit, exclude,
							Constants.type[k]);
					if (!assoc[k][j].loc.equals(assoc[0][j].loc)) {
						throw new RuntimeException("assume same locs");
					}
					phenos.addAll(assoc[k][j].pheno);
				}
			}
			boolean merge = Constants.mult();
			if (merge) {
				Assoc[][] ass1 = new Assoc[1][Constants.dir.length];
				for (int k = 0; k < dir.length; k++) {
					ass1[0][k] = assoc[0][k];
					for (int j = 1; j < Constants.type.length; j++) {
						ass1[0][k].merge(assoc[j][k]);
					}
					ass1[0][k].removeNaRows(true, false);
					ass1[0][k].removeNaCols(true);
					ass1[0][k].removeNaRows(true, true);
					// ass1[0][k].removeNaCols(false);
				}
				assoc = ass1;
			}

			
				for (int k = 0; k < assoc.length; k++) {
					List<String> snpids = new ArrayList<String>(
							assoc[k][0].probeids);
					PrintWriter pw = new PrintWriter(new BufferedWriter(
							new FileWriter(new File(System
									.getProperty("user.dir")
									+ "/resultMeta_"
									+ Constants.type[k]
									+ ".txt"))));
					List[] beta = new List[phenos.size()];
					List[] se = new List[phenos.size()];
					List[] pv = new List[phenos.size()];
					for(int kj = 0; kj<phenos.size(); kj++){
						 beta[kj] = new ArrayList<double[]>();
							se[kj] = new ArrayList<double[]>();
							pv[kj] = new ArrayList<double[]>();
					}
					
					for (int j = 0; j < assoc[k].length; j++) {
						Assoc assoc_ = assoc[k][j];
						for(int kj = 0; kj<phenos.size(); kj++){
							beta[kj].add(new double[assoc_.probeids.size()]);
							se[kj].add(new double[assoc_.probeids.size()]);
							pv[kj].add(new double[assoc_.probeids.size()]);
						}
						
						
						snpids.retainAll(assoc_.probeids);
						while(assoc_.pos<assoc_.genoline.size()){
							int start = assoc_.pos-1;
							if(Constants.buffer()>=0) assoc_.refreshGenotypes();
							for (int kj=0; kj< phenos.size();kj++) {
								String phen = phenos.get(kj);
							int ind = assoc_.pheno.indexOf(phen);
							if (ind >= 0) {
								assoc_.prepare(ind);
								assoc_.run(ind);
								
								 System.arraycopy(assoc_.beta, 0, beta[kj].get(j), start, assoc_.beta.length);
								 System.arraycopy(assoc_.se, 0, se[kj].get(j), start, assoc_.se.length);
								 System.arraycopy(assoc_.pv, 0, pv[kj].get(j), start, assoc_.pv.length);
								
								//pw[k].close();
							}
							}
						}
						for(int kk=0; kk<assoc_.pw.length; kk++){
							assoc_.pw[kk].close();
						}
					}
					int[][] alias = new int[assoc[k].length][];
					for (int j = 0; j < assoc[k].length; j++) {
						alias[j] = assoc[k][j].getPhenoInd(
								assoc[k][j].probeids, snpids, null);
					}
					for(int kj = 0; kj<phenos.size(); kj++){
						metaAnalyse(beta[kj], se[kj], pv[kj], pw, snpids, alias);
					}
					pw.close();
				}
				// assoc.run();
			//}

			System.exit(0);
		} catch (Exception exc) {
			exc.printStackTrace();
			System.exit(0);
		}
	}

	private void setNaToAvg() {
		// TODO Auto-generated method stub

	}

	private static void metaAnalyse(List<double[]> beta2, List<double[]> se2,
			List<double[]> pv2, PrintWriter pw2, List<String> snpid,
			int[][] alias) {
		int len = snpid.size();
		Meta meta = new Meta(beta2.size());
		double[] beta = new double[beta2.size()];
		double[] se = new double[beta.length];
		double[] pv = new double[beta.length];
		for (int i = 0; i < len; i++) {
			boolean notNa = false;
			for (int j = 0; j < beta.length; j++) {
				int i1 = alias[j][i];
				beta[j] = beta2.get(j)[i1];
				se[j] = se2.get(j)[i1];
				pv[j] = pv2.get(j)[i1];
				if (!Double.isNaN(pv[j])) {
					notNa = true;
				}
			}
			if (notNa) {
				try{
				meta.calc(beta, se, pv);
				}catch(Exception exc){exc.printStackTrace();}
				if (!Double.isNaN(meta.pFix)) {
					pw2.println(snpid.get(i) + "\t" + meta.betafixed + "\t"
							+ meta.sebetafixed + "\t" + meta.pFix + "\t"
							+ meta.countstudy);
				}
			}
		}
		pw2.close();

	}

	public static void main(String[] args) {
		try {
			Constants.parse(args, Constants.class);
			main();

		} catch (Exception exc) {
			exc.printStackTrace();
		}
	}

	static Rengine re = new Rengine(new String[] { "--vanilla" }, true,
			new LoopCall());
			
	static {
		re.eval("library(MASS)");
		re.eval("vecm = c(2,1,0)");
		re.eval("allFreq<-function(y) (vecm[1:length(y)] %*% y) / (2*(sum(y)))");
		re.eval("expDist<-function(y,p)   c(p^2, 2*p*(1-p), (1-p)^2)[1:length(y)]*sum(y)");
		re.eval("hwestat<-function(obs,exp) sum((obs-exp)^2/exp) ");
        re.eval("hwep<-function(y) if(length(which(y>0))<=1) 1.0 else 1-pchisq(hwestat(y,expDist(y,allFreq(y))),1) ");
        re.eval("hwepall<-function(y) min(hwep(y[1:min(3,length(y))]), hwep(y[3:min(5,length(y))]),na.rm=TRUE) ");
  /*  double pt0 = re.eval("hwep(c(25,50,25))").asDouble();
       double pt = re.eval("hwep(c(0,18,1403))").asDouble();
       double af = re.eval("allFreq(c(0,18,1403))").asDouble();
       double[] ed = re.eval("expDist(c(0,18,1403),"+af+")").asDoubleArray();
       double pt1 = re.eval("hwep(c(2808,0))").asDouble();
       double pt3 = re.eval("hwepall(c(2806,3))").asDouble();
       System.err.println("h");*/
		//re.eval("library(stepPlr)");
	}
	final PrintWriter[] pw;
	final double[][] phenotypes; // by phenotype, and then indiv
	final double[][] covariates; // by covariate, then indiv
	double[][] genotypes; // by snp, then indiv
	String[] pheno_name;

	List<String> indiv;
	final int[] snp_col_ind;
	BufferedReader cnvs;
	int[] pheno_alias; // maps pheno row to correct row in pheno matrix

	private void getSubset(List<String> split, String suff, List<String> res,
			List<String> res_all, List<String> exclude) {
		List<Integer> inds = new ArrayList<Integer>();
		for (int k = 0; k < split.size(); k++) {
			int ind = split.get(k).lastIndexOf(suff);
			if (ind >= 0) {
				String indv = split.get(k).substring(0, ind);
				if (!exclude.contains(indv)) {
					res.add(indv);
					res_all.add(indv);
					// inds.add(k);
				} else {
					res_all.add(split.get(k));
				}
			} else {
				res_all.add(split.get(k));
			}
		}

		// return inds;
	}

	private int lastIndexOf(String string, String[] suff) {
		for (int k = 0; k < suff.length; k++) {
			int i = string.lastIndexOf(suff[k]);
			if (i >= 0)
				return i;
		}
		return -1;
	}

	final List<String> covar, pheno, familyList;

	boolean transpose = false;
	// This is constructor

	List<String> probeids;
	List<String> loc = new ArrayList<String>();
	String chrom;

	Algebra alg = new Algebra();
	// String header =
	// "rsid    chrom   loc     pheno   N_case  N_control       beta    se      pvalue";
	String header = "rsid\tchrom\tloc\tpheno\tN_control\tN_case\thwe_control\thwe_case\tbeta\tse\tpvalue\tpvalue_fisher\tcontrol_CN\tcase_CN\tP(case|CN)";

	// first index of geno, second of pheno

	String assocline = "y ~ x";
	String assocCov = "y ~ ";

	public double[] fitCovar(int k, String family) {
		double beta = Double.NaN;
		double se = Double.NaN;
		double pv = Double.NaN;

		double[] pheno = this.phenotypes[k];
		double[] offset = new double[pheno.length];
		//
		// String family =this.getType(phenotypes[k], null, casesControls);
		if (family != null) {
			// if(pheno.length==0 || geno.length==0) throw new
			// RuntimeException("!!");
			re.assign("y", pheno);// new double[] {0,1,2,3});//pheno);
			// re.assign("x",geno);// new double[] {0,0,1,1});//geno);
			List<String>[] levels = new List[covariates.length];
			for (int k1 = 0; k1 < covariates.length; k1++) {
				re.assign(covar.get(k1) + "__", covariates[k1]);
				if (this.covar_type.get(this.covar.get(k1)).toLowerCase()
						.endsWith("factor")) {
					re.eval(covar.get(k1) + "__=" + "as.factor("
							+ covar.get(k1) + "__)");
					levels[k1] = Arrays.asList(re.eval(
							"levels(" + covar.get(k1) + "__)").asStringArray());
				}
			}
			String line = "summ= summary(glm(" + assocCov + ", fami=" + family
					+ "))";
			re.eval(line);
			REXP summ = re.eval("summ");
			// System.err.println(line);
			RList l1 = re.eval("dimnames(summ$coefficients)").asList();
			String[] names = l1.at(0).asStringArray();
			String[] names1 = l1.at(1).asStringArray();
			double[][] coeff = summ.asList().at("coefficients").asMatrix();
			for (int k1 = 0; k1 < offset.length; k1++) {
				double off = coeff[0][0];
				for (int j = 1; j < names.length; j++) {
					String[] str = names[j].split("__");
					double v = covariates[this.covar.indexOf(str[0])][k1];
					if (str.length == 1) {
						off += v * coeff[j][0];
					} else {
						if (Math.abs(v - Double.parseDouble(str[1])) < 0.001) {
							off += coeff[j][0];
						}
					}

				}
				offset[k1] = off;
			}
		}
		/*
		 * if(true){ re.assign("offs",offset); String line
		 * ="summ= summary(glm("+assocCov+", fami="+family+", offset=offs))";
		 * re.eval(line); REXP summ = re.eval("summ"); //
		 * System.err.println(line); RList l1 =
		 * re.eval("dimnames(summ$coefficients)").asList(); String[] names=
		 * l1.at(0).asStringArray(); String[] names1= l1.at(1).asStringArray();
		 * double[][]coeff = summ.asList().at("coefficients").asMatrix();
		 * System.err.println(coeff); }
		 */
		 return offset;
	}

	public void applyFilter(double diff, int cntThresh, double base) {
		for (int k = 0; k < this.indiv.size(); k++) {
			for (int j = 0; j < genotypes.length; j++) {
				int cnt = 0;
				// List<Double> l = new ArrayList<Double>();
				if (Math.abs(genotypes[j][k] - base) > diff) {
					cnt++;
					// l.add(genotypes[j][k]);
					inner1: for (int j1 = j - 1; j1 >= 0; j1--) {
						if (Math.abs(genotypes[j1][k] - base) > diff) {
							cnt++;
							// l.add(genotypes[j1][k]);
						} else
							break inner1;
					}
					inner1: for (int j1 = j + 1; j1 < genotypes.length; j1++) {
						if (Math.abs(genotypes[j1][k] - base) > diff) {
							cnt++;
							// l.add(genotypes[j1][k]);
						} else
							break inner1;
					}
				}
				if (cnt < cntThresh) {
					genotypes[j][k] = base;
				}
				/*
				 * else{ System.err.println("kept "+genotypes[j][k]); }
				 */
			}
		}
	}

	// double beta, se, pv;

	// final double[] beta, se, pv;

	public void removeNaRows(boolean all, boolean setNaAsAvg) {
		boolean[] include = new boolean[genotypes.length];
		// double[] weights = new double[genotypes[0].length];
		int cntInc = 0;
		for (int i = 0; i < include.length; i++) {
			double sum = 0;
			double cnt = 0;
			boolean hasNA = false;
			boolean allNA = true;
			for (int j = 0; j < genotypes[0].length; j++) {
				if (!Double.isNaN(genotypes[i][j])) {
					sum += genotypes[i][j];
					cnt++;
					allNA = false;
				} else {
					hasNA = true;
				}

			}
			double mean = sum / cnt;
			cnt = 0;
			sum = 0;
			for (int j = 0; j < genotypes[0].length; j++) {
				if (!Double.isNaN(genotypes[i][j])) {
					sum += Math.pow(genotypes[i][j] - mean, 2);
					cnt++;
				} else if (setNaAsAvg) {
					genotypes[i][j] = mean;
				}

			}
			double var = sum / cnt;
			include[i] = var <= 0 || (all && allNA || !all && hasNA) ? false
					: true;
			if (include[i])
				cntInc++;
		}
		int[] alias = new int[cntInc];
		cntInc = 0;
		for (int k = 0; k < include.length; k++) {
			if (include[k]) {
				alias[cntInc] = k;
				cntInc++;
			}
		}
		this.probeids = Arrays.asList((String[]) sublist(probeids
				.toArray(new String[0]), alias));
		this.loc = Arrays.asList((String[]) sublist(loc.toArray(new String[0]),
				alias));
		this.genotypes = (double[][]) sublist(genotypes, alias);
		this.maf = (double[]) sublistd(maf, alias);
	}

	public void removeNaCols(boolean all) {
		boolean[] include = new boolean[genotypes[0].length];
		// double[] weights = new double[genotypes[0].length];
		int cntInc = 0;
		for (int i = 0; i < include.length; i++) {
			// double sum =0;
			// double cnt=0;
			boolean hasNA = false;
			boolean allNA = true;
			for (int j = 0; j < genotypes.length; j++) {
				if (!Double.isNaN(genotypes[j][i])) {
					// sum+=genotypes[j][i];
					// cnt++;
					allNA = false;
				} else {
					hasNA = true;
				}

			}
			/*
			 * double mean = sum/cnt; cnt=0; sum=0; for(int j=0;
			 * j<genotypes.length; j++){ if(!Double.isNaN(genotypes[j][i])){
			 * sum+=Math.pow(genotypes[j][i]-mean,2); cnt++; }
			 * 
			 * } double var = sum/cnt;
			 */
			include[i] = (all && allNA || !all && hasNA) ? false : true;
			if (include[i])
				cntInc++;
		}
		int[] alias = new int[cntInc];
		cntInc = 0;
		for (int k = 0; k < include.length; k++) {
			if (include[k]) {
				alias[cntInc] = k;
				cntInc++;
			}
		}
		this.indiv = Arrays.asList((String[]) sublist(indiv
				.toArray(new String[0]), alias));
		for (int j = 0; j < phenotypes.length; j++) {
			this.phenotypes[j] = (double[]) sublistd(this.phenotypes[j], alias);
		}
		for (int j = 0; j < this.covariates.length; j++) {
			this.covariates[j] = (double[]) sublistd(this.covariates[j], alias);
		}
		for (int k = 0; k < genotypes.length; k++) {
			this.genotypes[k] = (double[]) sublistd(genotypes[k], alias);
		}
	}

	private Object sublist(Object genotypes2, int[] alias) {
		Object obj = Array.newInstance(Array.get(genotypes2, 0).getClass(),
				alias.length);
		for (int k = 0; k < alias.length; k++) {
			Array.set(obj, k, Array.get(genotypes2, alias[k]));
		}
		return obj;
	}

	private double[] sublistd(double[] genotypes2, int[] alias) {
		double[] obj = new double[alias.length];
		for (int k = 0; k < alias.length; k++) {
			obj[k] = genotypes2[alias[k]];// Array.get(genotypes2, alias[k]));
		}
		return obj;
	}

	public void runMult(int k, String family, int[] casesControls, boolean step) {
		Arrays.fill(beta, Double.NaN);
		Arrays.fill(se, Double.NaN);
		Arrays.fill(pv, Double.NaN);
		Arrays.fill(pvtab, Double.NaN);
		// double[] geno = this.genotypes[i];
		double[] pheno = this.phenotypes[k];
		// Arrays.fill(casesControls,-1);
		if (family != null) {
			if (pheno.length == 0)
				throw new RuntimeException("!!");
			re.assign("y", pheno);// new double[] {0,1,2,3});//pheno);
			StringBuffer assocline1 = new StringBuffer(assocline.substring(0,
					assocline.length() - 1));
			String[] names = new String[genotypes.length];
			for (int i = 0; i < genotypes.length; i++) {
				// if(include[i]){
				names[i] = "x_" + this.probeids.get(i).replace('.', '_');
				// String var =
				assocline1.append("+" + names[i]);
				re.assign(names[i], genotypes[i]);// new double[]
													// {0,0,1,1});//geno);
				// }
			}
			String line = assocline1 + ", family=" + family + "()";
			if (offset != null) {
				re.assign("offs", offset);
				line += " ,offset=offs";
			} else {
				for (int k1 = 0; k1 < covariates.length; k1++) {
					re.assign(this.covar.get(k1) + "__", covariates[k1]);
				}
			}
			String line1 = "l <-glm(" + line + ") ";
			re.eval(line1);
			if (step) {
				String line3 = "l<-stepAIC(l, steps=100" 
//						+",scope = list(lower = ~1), direction=\"forward\"" 
								+")";
				re.eval(line3);
			}
			String line2 = "summ= summary(l)";
			re.eval(line2);
			REXP summ = re.eval("summ");
			// System.err.println(line);

			RList names1 = re.eval("dimnames(summ$coefficients)").asList();
			RList suml = summ.asList();

			StringBuffer form = new StringBuffer("y~1");
			String form1 = new String("y~1");

			double[][] coeff = suml.at("coefficients").asMatrix();

			for (int i = 0; i < genotypes.length; i++) {
				// if(include[i]){
				int indrow = Arrays.asList(names1.at(0).asStringArray())
						.indexOf("x_" + this.probeids.get(i).replace('.', '_'));
				// String[] names1 = names.at(1).asStringArray();
				if (indrow >= 0 && indrow < coeff.length) {
					beta[i] = coeff[indrow][0];

					se[i] = coeff[indrow][1];
					pv[i] = coeff[indrow][3];
					form.append("+" + names[i]);
				}
				if (Constants.printNaN() || !Double.isNaN(pv[i])) {
					this.pw[k].print(this.probeids.get(i) + "\t");
					this.pw[k].print(this.chrom + "\t");
					this.pw[k].print(this.loc.get(i) + "\t");
					this.pw[k].print(this.pheno.get(k) + "\t");
					this.pw[k].print(this.casesControls[0] + "\t");
					this.pw[k].print(this.casesControls[1] + "\t");
					this.pw[k].print("NA\t");
					this.pw[k].print( "NA\t");
					this.pw[k].print(String.format("%5.3g", beta[i]) + "\t");
					this.pw[k].print(String.format("%5.3g", se[i]) + "\t");
					this.pw[k].print(String.format("%5.3g", pv[i]) + "\n");
					this.pw[k].print(String.format("%5.3g", pvtab[i]) + "\n");
				}
			}
			re.eval("hyper.lm = glm(as.formula(" + form
					+ ") ,fami=\"binomial\")");
			re.eval("hyper.lm1 = glm(as.formula(" + form1
					+ ") ,fami=\"binomial\")");
			// print(summary(hyper.lm))
			REXP su = re.eval("anova(hyper.lm1,hyper.lm,test =\"F\" )");
			String[] keys = su.asList().keys();
			// double[]d = su.asDoubleArray();
			REXP pr = su.asList().at(keys[5]);
			double pv_ = pr == null ? Double.NaN : pr.asDoubleArray()[1];
			// print()

			// double[] fst= suml.at("fstatist").asDoubleArray();
			// double pv_ =
			// re.eval("pf(summ$fst[1],summ$fst[2],summ$fst[3],lower.tail=F)").asDouble();
			this.pw[k].print("overall" + "\t");
			this.pw[k].print(this.chrom + "\t");
			this.pw[k].print("\t");
			this.pw[k].print(this.pheno.get(k) + "\t");
			this.pw[k].print(this.casesControls[0] + "\t");
			this.pw[k].print(this.casesControls[1] + "\t");
			// this.pw[k].print(String.format("%5.3g", )+"\t");
			// this.pw[k].print(String.format("%5.3g", se[i])+"\t");
			this.pw[k].print(String.format("%5.3g", pv_) + "\n");
		}
		// }
	}

	public void run(int i, int k, double[] offset, String family,
			int[] casesControls) {
		//System.err.println("i,k " +i+","+k);
		beta[i] = Double.NaN;
		se[i] = Double.NaN;
		pv[i] = Double.NaN;
		pvtab[i] = Double.NaN;
		String cnstring_case = null, cnstring_control = null, cont_prob = null;
		double[] geno = this.genotypes[i];
		double[] pheno = this.phenotypes[k];
		double hwe_control = Double.NaN;
		double hwe_case = Double.NaN;
		// Arrays.fill(casesControls,-1);
		if (family != null && (this.maf[i] > Constants.maf_thresh()) && (1-this.maf[i] > Constants.maf_thresh())) {
			if (pheno.length == 0 || geno.length == 0)
				throw new RuntimeException("!!");
			re.assign("y", pheno);// new double[] {0,1,2,3});//pheno);
		
			re.assign("x", geno);// new double[] {0,0,1,1});//geno);
			String line = assocline + ", family=" + family + "()";
			if (offset != null) {
				re.assign("offs", offset);
				line += " ,offset=offs";
			} else {
				for (int k1 = 0; k1 < covariates.length; k1++) {
					re.assign(this.covar.get(k1) + "__", covariates[k1]);
				}
			}
			String line1 = "summ= summary(glm(" + line + "))";
			re.eval(line1);
			REXP summ = re.eval("summ");
			// System.err.println(line);

			RList names = re.eval("dimnames(summ$coefficients)").asList();
			double[][] coeff = summ.asList().at("coefficients").asMatrix();
			int indrow = Arrays.asList(names.at(0).asStringArray())
					.indexOf("x");
			
			// String[] names1 = names.at(1).asStringArray();
			if (indrow >= 0 && indrow < coeff.length) {
				beta[i] = coeff[indrow][0];

				se[i] = coeff[indrow][1];
				pv[i] = coeff[indrow][3];
				
				re.eval("ind1 = !is.na(x)");
				if(family.startsWith("binom")){
				re.eval("indCont = as.factor(y)==" + 0);
				re.eval("indCase = as.factor(y)==" + 1);
				}else{
					re.eval("indCont = rep(TRUE,length(y))");
					re.eval("indCase = rep(FALSE,length(y))");
				}
				// re.eval("n_control = length(which(ind1 & indCont))");
				// re.eval("n_case = length(which(ind1 & indCase))");
			//	double[] d = re.eval("spl").asDoubleArray();
				re
						.eval("cn_case = hist(as.numeric(x[indCase & ind1]),br = spl,plot=F)$count");
				re
						.eval("cn_control = hist(as.numeric(x[indCont & ind1]),br = spl,plot=F)$count");
				//double[] d = re.eval("cn_case").asDoubleArray();
				 hwe_control = re.eval("hwepall(cn_control)").asDouble();
				 hwe_case = re.eval("hwepall(cn_case)").asDouble();
				re.eval("cnstring_case = paste(cn_case,collapse=\";\")");
				re.eval("cnstring_control = paste(cn_control,collapse=\";\")");
				re
						.eval("cont_prob = paste(round(cn_case/(cn_control+cn_case)*100)/100, collapse=\";\")");
				cnstring_case = re.eval("cnstring_case").asString();
				cnstring_control = re.eval("cnstring_control").asString();
				cont_prob = re.eval("cont_prob").asString();
				re.eval("inds = cn_case>0 | cn_control>0");
				re.eval("tab = rbind(cn_case[inds],cn_control[inds])");
				re.eval("vec1 = apply(tab,1,sum)");
				re.eval("vec2 = apply(tab,2,sum)");
				re.eval("tab = tab[vec1>0, vec2>0]");
				
				double minv = 	re.eval("min(tab)").asDouble();
			    double dim = re.eval("max(dim(tab))").asDouble();
				if(minv<Constants.fisherThresh && dim <=2){
				 pvtab[i] = 	re.eval("fisher.test(tab)$p").asDouble();
				}else{
				  pvtab[i] = 	re.eval("chisq.test(tab)$p.value").asDouble();
				}
				//double[][] tab = re.eval("tab").asDoubleMatrix();
				//System.err.println("h");
			}
		}else{
			System.err.println("excluded");
		}
			if (Constants.printNaN() || !Double.isNaN(pv[i])) {
				this.pw[k].print(this.probeids1.get(i) + "\t");
				this.pw[k].print(this.chrom + "\t");
				this.pw[k].print(this.loc1.get(i) + "\t");
				this.pw[k].print(this.pheno.get(k) + "\t");
				this.pw[k].print(this.casesControls[0] + "\t");
				this.pw[k].print(this.casesControls[1] + "\t");
				this.pw[k].print(String.format("%5.3g", hwe_control) + "\t");
				this.pw[k].print(String.format("%5.3g", hwe_case) + "\t");
				this.pw[k].print(String.format("%5.3g", beta[i]) + "\t");
				this.pw[k].print(String.format("%5.3g", se[i]) + "\t");
				this.pw[k].print(String.format("%5.3g", pv[i]) + "\t");
				this.pw[k].print(String.format("%5.3g", pvtab[i]) + "\t");
				this.pw[k].print(cnstring_control + "\t");
				this.pw[k].print(cnstring_case + "\t");
				this.pw[k].print(cont_prob + "\n");
				
			//}
		}
	}

	public double max(double[] d) {
		double m = Double.NEGATIVE_INFINITY;
		for (int k = 0; k < d.length; k++) {
			if (d[k] > m) {
				m = d[k];
			}
		}
		return m;
	}

	public double min(double[] d) {
		double m = Double.POSITIVE_INFINITY;
		for (int k = 0; k < d.length; k++) {
			if (d[k] < m) {
				m = d[k];
			}
		}
		return m;
	}

	int[] casesControls;

	final double[][]offsets;
	public void prepare(int k) {

		Arrays.fill(casesControls, -1);
		family = this.getType(phenotypes[k], null, casesControls);
		if(familyList.size()>0 && !familyList.get(k).equals(family)) throw new RuntimeException("problem with pheno file");
		if (Constants.residuals() && this.covar.size() > 0) {
			if(offsets[k]==null){
		  	offsets[k] = this.fitCovar(k, family);
			}
			this.offset = offsets[k];
		}
	}

	double[] offset = null;

	String family;

	double[] beta, se, pv, pvtab;

	public void run(int k) {
		// String[] family = new String[phenotypes.length];
		

		// for(int k=0; k<phenotypes.length; k++){

		if (!Constants.mult.toLowerCase().startsWith("t")) {
			for (int i = 0; i < genotypes.length; i++) {
				try {
					
					run(i, k, offset, family, casesControls);

				} catch (Exception exc) {
					exc.printStackTrace();
				}
			}

		} else {
			this.runMult(k, family, casesControls, Constants.step.toLowerCase()
					.startsWith("t"));
		}
		// }
		// for(int k=0; k<phenotypes.length; k++){
	
		// }
	}

	private void applyRecoding(int[][] coding1) {
		int len = genotypes.length;

		String[] pi = new String[len * (1 + coding1.length)];
		String[] loc = new String[len * (1 + coding1.length)];
		double[] maf = new double[len * (1 + coding1.length)];

		double[][] gen1 = new double[len * (1 + coding1.length)][genotypes[0].length];
		System.arraycopy(this.probeids.toArray(new String[0]), 0, pi, 0, len);
		System.arraycopy(this.loc.toArray(new String[0]), 0, loc, 0, len);
		System.arraycopy(this.maf, 0, maf, 0, len);
		System.arraycopy(this.genotypes, 0, gen1, 0, len);
		for (int k1 = 0; k1 < coding1.length; k1++) {
			int[] coding = coding1[k1];
			StringBuffer cod = new StringBuffer();
			for (int kj = 0; kj < coding.length; kj++) {
				cod.append("." + coding[kj]);
			}
			for (int j = 0; j < this.genotypes.length; j++) {
				int j1 = (k1 + 1) * len + j;
				pi[j1] = probeids.get(j) + "." + cod;
				loc[j1] = this.loc.get(j);
				maf[j1] = this.maf[j];
				for (int k = 0; k < genotypes[j].length; k++) {
					double gen = genotypes[j][k];
					if (!Double.isNaN(gen)) {
						int lb = (int) Math.floor(gen);
						int ub = (int) Math.ceil(gen);
						if (ub == lb) {
							gen1[j1][k] = coding[ub];
						} else {
							double w1 = (gen - lb);
							double w2 = ub - gen;
							gen1[j1][k] = (coding[lb] * w1 + coding[ub] * w2);
						}
					}
					else{
						gen1[j1][k] = Double.NaN;
					}
				}
			}
		}

		this.genotypes = gen1;
		;
		this.probeids = new ArrayList<String>(Arrays.asList(pi));
		this.loc = new ArrayList<String>(Arrays.asList(loc));
		this.maf = maf;
	}

	private static String getType(double[] ds, double[] geno, int[] caseControl) {

		SortedMap<Double, Integer> m = new TreeMap<Double, Integer>();
		// SortedMap<Double, Integer> m_geno = new TreeMap<Double, Integer>();
		for (int k = 0; k < ds.length; k++) {
			if (!Double.isNaN(ds[k])
					&& (geno == null || !Double.isNaN(geno[k]))) {
				Integer cnt = m.get(ds[k]);
				m.put(ds[k], cnt == null ? 1 : cnt + 1);
				// Integer cnt1 = m_geno.get(ds[k]);
				// m_geno.put(geno[k], cnt1==null ? 1 : cnt1+1);

			}
		}
		if (m.size() <= 1) {
			return null;
		} else if (m.size() == 2) {
			caseControl[0] = m.get(m.firstKey());
			caseControl[1] = m.get(m.lastKey());
			double first = m.firstKey();
			double last = m.lastKey();
			if (first != 0 || last != 1) {
				for (int k = 0; k < ds.length; k++) {
					if (ds[k] == first)
						ds[k] = 0;
					else if (ds[k] == last)
						ds[k] = 1;
				}
			}
			return "binomial";
		} else {
			return "gaussian";
		}
	}

	final Map<String, Integer>[] pheno_decoding;
	final String[] pheno_family;
	Map<String, String> covar_type = new HashMap<String, String>();
	Map<String, Integer>[] covar_decoding;

	public void merge(Assoc assoc) {
		int len1 = probeids.size();
		int len2 = assoc.probeids.size();
		for (int i = 0; i < assoc.probeids.size(); i++) {
			this.probeids.add(assoc.probeids.get(i) + "_" + assoc.suff);
			this.loc.add(assoc.loc.get(i));
		}
		if (!this.indiv.equals(assoc.indiv) || !this.pheno.equals(assoc.pheno))
			throw new RuntimeException("!!");
		double[] maf = new double[len1 + len2];
		System.arraycopy(this.maf, 0, maf, 0, this.maf.length);
		System.arraycopy(assoc.maf, 0, maf, len1, assoc.maf.length);
		double[][] genotypes1 = new double[len1 + len2][genotypes[0].length];
		for (int i = 0; i < genotypes1[0].length; i++) {
			for (int j = 0; j < len1; j++) {
				genotypes1[j][i] = genotypes[j][i];
			}
			for (int j = 0; j < len2; j++) {
				genotypes1[j + len1][i] = assoc.genotypes[j][i];
			}
		}
		this.genotypes = genotypes1;
		this.maf = maf;
	}

	String suff;
int pos = 0;
//final int buffer = Constants.buffer;
	 public Assoc(File avgFile, File phenoFile, File limitFile,
			File excludeFile, String suff) throws Exception {
		covar = new ArrayList<String>();
		this.suff = suff;
		pheno = new ArrayList<String>();
		familyList = new ArrayList<String>();
		chrom = Constants.chrom();
		this.avgFile = avgFile;
		this.casesControls = new int[2];
		cnvs = new BufferedReader(new FileReader(avgFile));
		List<String> avg_file_indiv = new ArrayList<String>();
		List<String> avg_file_indiv_all = new ArrayList<String>();
		// new ArrayList<String>();
		List<String> cnvline = Arrays.asList(cnvs.readLine().split("\t"));
		transpose = !cnvline.get(1).equals("loc");

		String st = "";
		 genoline = new ArrayList<String>();
		genoline.add("snpid");
		for (int i = 0; (st = cnvs.readLine()) != null; i++) {
		//	if(st.replaceAll(" ", "").length()>0){
			String[] str_ = st.split("\t");
			if(str_[0].length()>0){
			genoline.add(str_[0]);
			if (transpose) {
				if (i == 0) {
					List<String> l1 = Arrays.asList(str_);
					loc.addAll(l1.subList(1, l1.size()));
				}
			} else {
				loc.add(str_[1]);
			}
			}
		//	}
		}
		cnvs.close();

		if (transpose) {
			List<String> tmp = genoline;
			genoline = cnvline;
			cnvline = tmp;
		}
		// cnvline = cnvline.subList(1, cnvline.size());
	
		getSubset(cnvline, suff, avg_file_indiv, avg_file_indiv_all,
				read(excludeFile));

		Map<String, String> subsPhenoToAvg = new HashMap<String, String>();
		BufferedReader br1 = new BufferedReader(new FileReader(phenoFile));
		List<String> phen_head = Arrays.asList(br1.readLine().split("\t")); // header

		read(limitFile, covar, "covar", phen_head, this.covar_type, null);
		read(limitFile, pheno, "pheno", phen_head, this.covar_type, this.familyList);
		if (pheno.size() == 0) {
			pheno.addAll(phen_head.subList(1, phen_head.size()));
		}
		this.offsets = new double[pheno.size()][];
		int[] covar_ind = getPhenoInd(covar, phen_head, null);
		int[] pheno_ind = getPhenoInd(pheno, phen_head, null);

		List<String> pheno_indiv = new ArrayList<String>();
		for (; (st = br1.readLine()) != null;) {
			String[] str = st.split("\t");
			pheno_indiv.add(str[0]);
		}

		br1.close();
		List<String> avg_file_indiv_mod = this.reconcile(pheno_indiv,
				avg_file_indiv, subsPhenoToAvg);

		indiv = new ArrayList<String>(pheno_indiv);// reconcile(pheno_indiv,
													// avg_file_indiv,
													// subsAvgToPheno);
		indiv.retainAll(avg_file_indiv_mod);
		int[] pheno_indiv_ind = getPhenoInd(indiv, pheno_indiv, null);
		this.snp_col_ind = getPhenoInd(indiv, avg_file_indiv_all,
				subsPhenoToAvg);
		this.covariates = new double[covar.size()][indiv.size()];
		this.phenotypes = new double[pheno.size()][indiv.size()];
		pheno_decoding = new Map[pheno.size()];
		covar_decoding = new Map[covar.size()];
		pheno_family = new String[pheno.size()];
		String[] covar_family = new String[covar.size()];
		readMatrix(covar_ind, pheno_indiv_ind, phenoFile, covar.size(),
				covariates, covar_decoding, pheno_family, true);
		readMatrix(pheno_ind, pheno_indiv_ind, phenoFile, pheno.size(),
				phenotypes, pheno_decoding, covar_family, true);
		if(Constants.limit < genoline.size()) {
			genoline = genoline.subList(0, Constants.limit+1);
			this.loc = loc.subList(0, Constants.limit);
		}
		probeids = genoline.subList(1, genoline.size());
		
		pos = 1;
	//	refreshGenotypes();
    //    pos = buffer +1;
		pw = new PrintWriter[this.pheno.size()];
		for (int k = 0; k < pw.length; k++) {
			pw[k] = new PrintWriter(new BufferedWriter(new FileWriter(new File(
					avgFile.getParentFile(), "result_" + pheno.get(k) + "_"
							+ suff + ".txt"))));
			pw[k].println(header);
			pw[k].flush();
		}
		for (int k = 0; k < covariates.length; k++) {
			if (Constants.residuals()) {
				this.assocCov += "+" + covar.get(k) + "__";
			} else {
				assocline += "+" + covar.get(k) + "__";
				;
			}
		}
		
		if(Constants.buffer()<0) this.refreshGenotypes(); 	
	}

	 List genoline;
	 final File avgFile;
	List<String> probeids1;
	List<String> loc1;
	 BufferedReader avgFileReader=null;
	private void refreshGenotypes() throws Exception{
		int start = pos;
		int end = Constants.buffer()<0 ? genoline.size() : Math.min(this.genoline.size(), Constants.buffer()+start);
		
		int[] geno_col_ind = getPhenoInd(probeids, genoline, null);
		probeids1 = genoline.subList(start, end);
		loc1 = this.loc.subList(start -1, end-1);
		
			beta = new double[this.loc1.size()];
			se = new double[this.loc1.size()];
			pv = new double[this.loc1.size()];
			pvtab = new double[this.loc1.size()];
		if (transpose) {
			if(Constants.buffer()>=0 && Constants.buffer() < probeids1.size()) throw new RuntimeException("!! bufferring does not work unless transposed");
			this.genotypes = new double[probeids1.size()][indiv.size()];
			;
			readMatrix(geno_col_ind, snp_col_ind, avgFile,
					this.probeids.size(), genotypes, null, null, false);
		} else {
			if(avgFileReader==null){
				avgFileReader = new BufferedReader(new FileReader(avgFile));
				avgFileReader.readLine();
			}
			this.genotypes = new double[probeids1.size()][indiv.size()];
			readMatrixT(snp_col_ind, geno_col_ind, avgFileReader, this.indiv.size(),
					genotypes, false, pos, end);
			//genotypes = transpose(genotypes1);//
		}
		pos = end;
		if(Constants.round()>0) round(genotypes,Constants.round());
		if(pos>=this.genoline.size()&& avgFileReader!=null) avgFileReader.close(); 
		int max = (int) Math.ceil(this.getMax(genotypes));
		spl = new double[max + 2];
		for (int i = 0; i < spl.length; i++) {
			spl[i] = i - 0.5;
		}
		this.re.assign("spl", spl);
		this.maf = new double[this.genotypes.length];
		for (int k = 0; k < genotypes.length; k++) {
			double[] hist = getHist(genotypes[k]);
			double sum = sum(hist);
			if (suff.endsWith("countAll")) {
				maf[k] = hist.length <= 2 ? 1 : 1 - hist[2] / sum;
			}else{
				maf[k] = 1-hist[0]/sum;
			}
			/*else {
				maf[k] = 1.0;
			}*/
			// if(maf[k]>0){
			// double[] g = genotypes[k];
			// System.err.println("h");
			// }

		}
		double thresh = Double.parseDouble(Constants.thresh);
		int cnt = Integer.parseInt(Constants.minCount);
		if (cnt > 0) {
			this.applyFilter(thresh, cnt, suff.endsWith("countAll") ?2 : 0);
		}
		if (Constants.coding() != null) {
			this.applyRecoding(Constants.coding());
		}
		
	}

	private void round(double[][] genotypes2, double round) {
		for(int i=0; i<genotypes.length; i++){
			for(int j=0; j<genotypes[i].length; j++){
				double rnd = Math.round(genotypes[i][j]);
				if(Math.abs(rnd - genotypes[i][j])<round){
					genotypes[i][j] = rnd;
				}
			}
		}
		
	}

	double[] spl;
	double[] maf;

	private double sum(double[] hist) {
		double sum = 0;
		for (int k = 0; k < hist.length; k++) {
			sum += hist[k];
		}
		return sum;
	}

	private double getMax(double[][] geno) {
		double m = Double.NEGATIVE_INFINITY;
		int ik  =0;
		int ij = 0;
		for (int k = 0; k < geno.length; k++) {
			for (int j = 0; j < geno[k].length; j++) {
				if (geno[k][j] > m) {
					m = geno[k][j];
					ij = j;
					ik = k;
				}
			}
		}
		return m;
	}

	private double[] getHist(double[] ds) {
		re.assign("x", ds);
		re.eval("ind1 = !is.na(x)");
		// re.eval("indCont = as.factor(y)=="+0);
		// re.eval("indCase = as.factor(y)=="+1);
		// re.eval("n_control = length(which(ind1 & indCont))");
		// re.eval("n_case = length(which(ind1 & indCase))");
		// re.eval("cn_case = hist(as.numeric(x[ind1]),br = spl,plot=F)$count");
		re.eval("cn_control = hist(as.numeric(x[ind1]),br = spl,plot=F)$count");

		// re.eval("cnstring_case = paste(cn_case,collapse=\";\")");
		// re.eval("cnstring_control = paste(cn_control,collapse=\";\")");
		// re.eval("cont_prob = paste(round(cn_case/(cn_control+cn_case)*100)/100, collapse=\";\")");
		// cnstring_case = re.eval("cnstring_case").asString();
		return re.eval("cn_control").asDoubleArray();

		// cont_prob = re.eval("cont_prob").asString();
	}

	private List<String> reconcile(List<String> pheno_indiv,
			List<String> avg_file_indiv, Map<String, String> phenoToAvg) {
		SortedSet<String> pheno = new TreeSet<String>(pheno_indiv);
		SortedSet<String> avg = new TreeSet<String>(avg_file_indiv);
		List<String> res = new ArrayList<String>();

		for (Iterator<String> it = pheno.iterator(); it.hasNext();) {
			String nxt = it.next();
			SortedSet<String> tail = avg.tailSet(nxt);
			boolean assg = false;
			if (tail.size() > 0) {
				String first = tail.first();
				if (first.equals(nxt)) {
					res.add(nxt);
					assg = true;
					avg.remove(first);
				} else if (accept(first, nxt, phenoToAvg)) {
					res.add(nxt);
					assg = true;
					avg.remove(first);
				}

			}
			if (!assg) {
				SortedSet<String> head = avg.headSet(nxt);
				if (head.size() > 0) {
					String last = head.last();
					if (accept(last, nxt, phenoToAvg)) {
						res.add(nxt);
						assg = true;
						avg.remove(last);
					}
				}
			}
			if (!assg) {
				System.err.println("not found avg for " + nxt);
			}

		}
		System.err.println("not found " + avg);
		return res;
		// pheno_indiv.retainAll(res);
		// return pheno_indiv;
	}

	private boolean accept(String avg, String pheno, Map<String, String> subs) {
		if (avg.startsWith(pheno) || pheno.startsWith(avg)) {
			if (subs.containsKey(avg))
				throw new RuntimeException("!!");
			System.err.println("accepted substitute " + avg + "->" + pheno);
			if (!pheno.equals(avg)) {
				subs.put(avg, pheno);
			}
			return true;
		}
		subs.put(pheno, "NOT FOUND");
		return false;
	}

	private double[][] transpose(double[][] m) {
		double[][] res = new double[m[0].length][m.length];
		for (int k = 0; k < m.length; k++) {
			for (int j = 0; j < res.length; j++) {
				res[j][k] = m[k][j];
			}
		}
		return res;
	}

	private void readMatrix(int[] col_ind, int[] row_ind, File phenoFile,
			int numcols, double[][] phenotypes,
			Map<String, Integer>[] decoding, String[] type, boolean header)
			throws Exception {

		for (int k = 0; k < phenotypes.length; k++) {
			Arrays.fill(phenotypes[k], Double.NaN);
		}
		// DoubleMatrix2D phenotypes = new DenseDoubleMatrix2D(indiv.size(),
		// numcols);
		if (numcols == 0)
			return;
		BufferedReader br1 = new BufferedReader(new FileReader(phenoFile));

		String st = header ? br1.readLine() : "";
		// boolean[] nonNumeric = new boolean[phenotypes.length];
		for (int i = 0; (st = br1.readLine()) != null; i++) {
			String[] str = st.split("\t");
			if (i == 0 && decoding != null) {
				for (int k = 0; k < str.length; k++) {
					if (col_ind[k] >= 0) {
						int col_ = col_ind[k];
						String v = str[k];
						try {
							Double.parseDouble(v);
						} catch (Exception exc) {
							decoding[col_] = new HashMap<String, Integer>();
						}
					}
				}
			}

			if (row_ind[i] >= 0) {
	
				for (int k = 0; k < str.length; k++) {
					if (col_ind[k] >= 0) {
						int col_ = col_ind[k];
						String v = str[k];
						if (v.equals("NA"))
							v = "NaN";
						Map<String, Integer> decoding_ = decoding == null ? null
								: decoding[col_];
						phenotypes[col_][row_ind[i]] = decoding_ != null ? decode(
								decoding_, v)
								: Double.parseDouble(v);

					}
				}

			}
		}
		br1.close();
	
	}
	
	private void readMatrixT(int[] col_ind, int[] row_ind, BufferedReader br1,
			int numcols, double[][] phenotypes,
			boolean header, int start, int end)
			throws Exception {

		for (int k = 0; k < phenotypes.length; k++) {
			Arrays.fill(phenotypes[k], Double.NaN);
		}
		// DoubleMatrix2D phenotypes = new DenseDoubleMatrix2D(indiv.size(),
		// numcols);
		if (numcols == 0)
			return;
	//	BufferedReader br1 = new BufferedReader(new FileReader(phenoFile));
		int offset = start -1;
		String st = header ? br1.readLine() : "";
		// boolean[] nonNumeric = new boolean[phenotypes.length];
		for (int i = start; i<end && (st = br1.readLine()) != null; i++) {
			String[] str = st.split("\t");
			int len = st.length();
			
			if (i>=start && row_ind[i] >= 0) {
	
				for (int k = 0; k < str.length; k++) {
					if (col_ind[k] >= 0) {
						int col_ = col_ind[k];
						String v = str[k];
						if (v.equals("NA")  || v.length()==0)
							v = "NaN";
						phenotypes[row_ind[i]-offset][col_] =Double.parseDouble(v);

					}
				}

			}
		}
	//	br1.close();
	
	}

	private double decode(Map<String, Integer> map, String v) {
		Integer d = map.get(v);
		if (d == null) {
			map.put(v, d = map.size());
		}
		return d.doubleValue();
	}

	private int[] getPhenoInd(List<String> subsample, List<String> list,
			Map<String, String> m) {
		int[] res = new int[list.size()];
		for (int k = 0; k < res.length; k++) {
			String st = list.get(k);
			String st1 = st;
			// if(st.equals("loc")) continue;
			if (m != null && m.containsKey(st)) {
				st1 = m.get(st);
			}

			res[k] = subsample.indexOf(st1);
			// if(k>0 && res[k]==res[k-1] && res[k]>=0 ){
			// System.err.println("h");
			// }
			if (res[k] < 0) {
				System.err.println(st1);
			} else {
				// System.err.println(st1);
			}
		}
		return res;
	}

	private List<String> read(File excludeFile) throws Exception {
		List<String> res = new ArrayList<String>();
		if (excludeFile.exists()) {
			BufferedReader br = new BufferedReader(new FileReader(excludeFile));
			String st = br.readLine();
			int pat_ind =  Arrays.asList(st.split("\t")).indexOf("PATIENT");
			while ((st = br.readLine()) != null) {
				res.add(st.split("\t")[pat_ind]);
			}
			br.close();
		}

		return res;
	}

	private void read(File excludeFile, List<String> res, String match,
			List<String> phens, Map<String, String> type, List<String> family) throws Exception {

		if (excludeFile.exists()) {
			BufferedReader br = new BufferedReader(new FileReader(excludeFile));
			String st = "";
			while ((st = br.readLine()) != null) {
				if (!st.startsWith("#")) {
					String[] str = st.split("\t");
					if (str[1].startsWith(match) && (phens.contains(str[0]))) {
						res.add(str[0]);
						type.put(str[0], str[1]);
						if(family!=null && str.length>2){
							family.add(str[2]);
						}
					}
					else if (str[0].startsWith(match) && (phens.contains(str[1]))) {
						res.add(str[1]);
						type.put(str[1], str[0]);
						if(family!=null && str.length>2){
							family.add(str[2]);
						}
					}
				}
			}
			br.close();
		}

	}

}
