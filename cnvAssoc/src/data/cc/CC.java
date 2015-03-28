package data.cc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.regression.OLSMultipleLinearRegression;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.Rengine;

import cnvtrans.LoopCall;

import data.Dirichlet;

public class CC {
	
	static double thresh = 0.001; //threshold for iteration to find maximally correlated vector
	static double u = 10;// parameter for initialise random weights
	public static Boolean plotCC = false;
	static boolean normalise = false;
	
	public Double[] q = new Double[] { 1.0,1.0, 5.0, 4.0, 4.0,4.0};

	
	//public static Boolean reportAnova =true;
//	public static Integer stepAIC = 0;
	static boolean allowBinary = false;
	static boolean addOne = true;  //do we add column of ones
	static boolean centralise = true;  //do we centralise prior to CC
//	static double threshold = 0.2;        
	public static File dir2 = new File( "results");


//this has a list of all the maximally correlated vectors.  	
	List<RealVector>[] correlatedVectors;
	
	//this has a list of the betas chosen at each iteration
	List<RealVector>[] combinations;
	
	//this has a list of  the correlation matrices at each iteration
	List<RealMatrix> covar;;
	
	//these are the regularisation matrices
	RealMatrix[] Q;

	
	/** this makes the Q matrices */
	public void makeGamma(Double[] q,  int[] band) {
		Q = new RealMatrix[mats.length];
		
		for (int k = 0; k < mats.length; k++) {
			if(q[k]==0 || band[k]==0) continue;
			if (mats[k] == null)
				continue;
			double[] array = new double[mats[k].getColumnDimension()];
			Arrays.fill(array, q[k]);
			Q[k] = MatrixUtils.createRealDiagonalMatrix(array);
		   System.err.println("h");
/*			RealMatrix gamma = new OpenMapRealMatrix(this.mats[k]
					.getColumnDimension(), this.mats[k].getColumnDimension());
			for (int j = 0; j < gamma.getColumnDimension(); j++) {
				for (int i = 0; i < gamma.getColumnDimension(); i++) {
					int bnd = Math.abs(i - j);
					if (bnd <= band[k]) {
						gamma.setEntry(j, i, q[k]);
					} 
				}

			}
			Q[k] = gamma.transpose().multiply(gamma);
			*/
		}
	}

	
	public static Integer perms = 0;
	public static String base = ".";
	public static String[] args_;
	

	private void permutePhen() {
		List<Integer> s = new ArrayList<Integer>();
		int rows = data[data.length-1].data.getRowDimension();
		for(int k=0; k<rows; k++){
			s.add(k);
		}
		double sze = s.size();
		List<Integer> s1 = new ArrayList<Integer>();
		for(int k=0; sze>0; k++){
			int ind = s.remove((int) Math.floor(Math.random()*sze));
			s1.add(ind);
//			indiv.set(k, indiv_.get(ind));
			sze = s.size();
		}
		
		 ((Data)this.data[data.length-1]).reorder(s1);
		initialise();
		
			centralise(true, true, true);
		updateY();
		
		
		//this.data[data.length-1];
		
	}

	

	private List<String> baf_snps() {
		return this.data[baf_index].probes;
	}

	private List<String> lrr_snps() {
		return this.data[lrr_index].probes;
	}

	int baf_index, lrr_index;

	//public static Boolean joinVntr = false;

	private static void deletePdfs(File dir) {
		System.err.println("deleting pdfs from "+dir.getAbsolutePath());
		if(!dir.exists()) return;
		File[] pdfs = (dir).listFiles(new FileFilter() {

			@Override
			public boolean accept(File arg0) {
				String n = arg0.getName();
				return n.endsWith(".pdf")
						&& (n.startsWith("cor") || n.startsWith("out_") || n
								.startsWith("cca"));

			}

		});
		if(pdfs!=null){
		for (int k = 0; k < pdfs.length; k++) {
			pdfs[k].delete();
		}
		}

	}

	private static Data getData(File[] files, String[] object, String chr,
			int[] start, boolean modify, boolean c, String[] asList,
			boolean only, String type_name) throws Exception {
		boolean remZero = files.length == 1;
		Data d = null;
		boolean correct = type_name.equals("lrr");
		for (int k = 0; k < files.length; k++) {
			try {
				Data d1 = new Data(files[k], object, files[k].getName().substring(0,files[k].getName().length()-4), start, modify, c,
						asList, only, type_name, remZero, correct);
				if (d == null)
					d = d1;
				else
					d.append(d1);
			} catch (Exception exc) {
				exc.printStackTrace();
				// System.err.println(exc.getMessage());
			}
		}
		return d;
	}

	public static int[] thin(AbstractData data, List<String> probes) {
		List<Integer> l1 = new ArrayList<Integer>();
		List<Integer> l2 = new ArrayList<Integer>();
		getAlias(data.probes, probes, l1, l2);
		int[] alias = get(l1);

		data.probes = sub(data.probes, alias);
		if (data.loc.size() > 0) {
			data.loc = subI(data.loc, alias);
		}
		data.data = data.data.getSubMatrix(getRows(data.indiv.size()), alias);
		return alias;

	}

	public void thin(int k, List<String> probes) {
		int[] rows = getRows(indiv.size());
		int[] alias = thin(data[k], probes);
		this.matsOrig[k] = this.matsOrig[k].getSubMatrix(rows, alias);
		this.mats[k] = matsOrig[k].copy();
		this.cols[k] = new int[matsOrig[k].getColumnDimension()];
		for(int j=0; j<cols[k].length; j++){
			cols[k][j] = j;
		}
	}

	double[][] means;

	public void centralise(boolean setNaAsMean, boolean centralise,
			boolean standardise) {
		means = new double[mats.length][];
		for (int k = 0; k < this.mats.length; k++) {
			if (mats[k] == null)
				continue;
			// if(matsOrig[k]!=null){
			// na[k] =new boolean[ matsOrig[k].getRowDimension()][
			// matsOrig[k].getColumnDimension()];
			if (k != this.pheno_index && !this.binary[k]) {
				means[k] = new double[matsOrig[k].getColumnDimension()];
				this.centralise(this.data[k].type_name, this.data[k].probes,
						this.matsOrig[k], setNaAsMean, centralise
								&& k < 2, standardise && k < 2, means[k],
						this.alias[k]);
				this.mats[k] = matsOrig[k].copy();
			}
			// }
		}
		// this.updateY();
	}

	public void centralise(String nme, List<String> probes, RealMatrix m,
			boolean setNaAsMean, boolean centralise, boolean std,
			double[] means, int[] alias

	) {

		for (int k = 0; k < m.getColumnDimension(); k++) {
			if (probes.get(k).equals("ones"))
				continue;
			double sum = 0;
			double cnt = 0;
			double cntNA = 0;
			Set<Double> vals = new HashSet<Double>();
			for (int j = 0; j < m.getRowDimension(); j++) {
				double v = m.getEntry(j, k);
				if (!Double.isNaN(v)) {
					sum += v;
					cnt++;
					vals.add(v);
					// na[j][k] = true;
				} else if (alias[j] >= 0) {
					cntNA++;
				}
				// else{
				// na[j][k] = false;
				// }
			}
			if (setNaAsMean && cntNA / (cnt + cntNA) > 0.1) {

				throw new RuntimeException("!! "
						+ (cntNA + " " + (cnt + cntNA)) + " " + nme + " "
						+ probes.get(k));
			}
			if (vals.size() <= 2) {
				continue;
			}
			double mean = sum / cnt;
			means[k] = mean;
			double var = 1;

			if (std) {
				sum = 0;
				for (int j = 0; j < m.getRowDimension(); j++) {
					double v = m.getEntry(j, k);
					if (!Double.isNaN(v)) {
						sum += Math.pow(v - mean, 2);
					}
				}
				var = Math.sqrt(sum / cnt);
				if (sum == 0 || cnt == 0)
					var = 1.0;
			}
			if (Double.isNaN(mean)) {
				throw new RuntimeException("exc");
			}
			for (int j = 0; j < m.getRowDimension(); j++) {
				double v = m.getEntry(j, k);
				if (!Double.isNaN(v)) {
					// sum += Math.pow(d[j] - mean, 2);
					// cnt++;
					if (centralise) {
						m.setEntry(j, k, (v - mean) / var);
					}
				} else if (setNaAsMean && alias[j] >= 0) {

					m.setEntry(j, k, centralise ? 0 : mean);
				}

			}
		}
	}

	public void set(CC data) {
	
		List<RealVector>[] combinations = data.combinations;
		this.combinations = new List[combinations.length];

		RealMatrix[] matsN = new RealMatrix[combinations.length];
		int[] rows = getRows(this.indiv.size());
		int[] rows_d = new int[] { 0 };
		int[][] alias = new int[data.mats.length][];
		int[][] alias_d = new int[data.mats.length][];
		for (int k = 0; k < data.mats.length; k++) {
			if (this.data[k] == null || data.data[k] == null)
				continue;
			this.combinations[k] = new ArrayList<RealVector>();
			List<Integer> l1 = new ArrayList<Integer>();
			List<Integer> l2 = new ArrayList<Integer>();
			getAlias(this.data[k].probes, data.data[k].probes, l1, l2);
			alias[k] = get(l1);
			alias_d[k] = get(l2);
			if(this.data[k].type_name.equals("baf")){
			   for(int j=0; j<l1.size(); j++){
				   int j1 = l1.get(j);
				   int j2 = l2.get(j);
				   String probe1 = this.data[k].probes.get(j1);
				   String probe2 = data.data[k].probes.get(j2);
				   if(!probe1.equals(probe2)) throw new RuntimeException("!! prob with probes "+probe1+" "+probe2 );
				   checkAndFix(this.data[k].data, data.data[k].data,j1,j2, probe1);
				  
				   
			   }
			}
			this.data[k].probes = sub(this.data[k].probes, alias[k]);
		//	this.data[k].data = 	this.data[k].data.getSubMatrix(rows, alias[k]);
			if (this.data[k].loc.size() > 0)
				this.data[k].loc = subI(this.data[k].loc, alias[k]);
		}
		for (int k = 0; k < correlatedVectors.length; k++) {
			if (this.correlatedVectors[k] == null || data.data[k] == null)
				continue;
			this.correlatedVectors[k].clear();
		}
		this.covar.clear();
		for (int j = 0; j < combinations[0].size(); j++) {
			for (int k = 0; k < correlatedVectors.length; k++) {
				if (this.correlatedVectors[k] == null)
					continue;
				RealVector v = combinations[k].get(j);

				RealMatrix ma = this.matsOrig[k].getSubMatrix(rows, alias[k]);
				RealMatrix m = (new Array2DRowRealMatrix(v.getData()))
						.getSubMatrix(alias_d[k], rows_d);
				this.combinations[k].add(m.getColumnVector(0));
				matsN[k] = ma;
				// this.na[k] = getSub()
				RealVector corre = ma.multiply(m).getColumnVector(0);

				this.correlatedVectors[k].add(corre);
				this.y.setColumn(k, corre.getData());
				// this.y[k].
				// ols.

			}
			RealMatrix cor = computeCorrelationMatrix(this.y);
			// cor = corr.getEntry(0,1);
			this.covar.add(cor);
		}
		this.matsOrig = matsN;
	 
		// this.combinations = combinations;
	}

	public static void checkAndFix(RealMatrix data1, RealMatrix data2, int j1, int j2, String probe1) {
		 double freq = getFreq(data1.getColumn(j1));
		   double freq1 = getFreq(data2.getColumn(j2));
		   double small = Math.min(freq, freq1);
		   double big = Math.max(freq, freq1);
		   double diff = Math.abs(freq-freq1);
		   if(diff>0.2 && (Math.abs(freq-(1-freq1))<diff)){
			   data2.setColumnVector(j2, data2.getColumnVector(j2).mapMultiply(-1).mapAdd(2));
			   try{
			   throw new RuntimeException("big freq diff at "+probe1+" "+freq1+" "+freq);
			   }catch(RuntimeException exc){
				   exc.printStackTrace();
			   }
		   }
		
	}

	public static double getFreq(double[] column) {
		double sum=0;
		double cnt=0;
		for(int k=0; k<column.length; k++){
			if(!Double.isNaN(column[k])){
				sum+=column[k];
				cnt+=2.0;
			}
		}
		return sum/cnt;
	}

	public static List<String> sub(List<String> probes, int[] is) {
		List<String> s = new ArrayList<String>();
		for (int k = 0; k < is.length; k++) {
			s.add(probes.get(is[k]));
		}
		return s;
	}

	public static List<Integer> subI(List<Integer> probes, int[] is) {
		List<Integer> s = new ArrayList<Integer>();
		for (int k = 0; k < is.length; k++) {
			if (is[k] < probes.size())
				s.add(probes.get(is[k]));
		}
		return s;
	}

	public static int[] getRows(int size) {
		int[] res = new int[size];
		for (int k = 0; k < res.length; k++) {
			res[k] = k;
		}
		return res;
	}

	public static void getAlias(List<String> probes, List<String> probes2,
			List<Integer> l1, List<Integer> l2) {

		for (int k = 0; k < probes.size(); k++) {
			int index = probes2.indexOf(probes.get(k));
			if (index >= 0) {
				l1.add(k);
				l2.add(index);
			}
		}

	}

	public static int[] get(List<Integer> l1) {
		int[] res = new int[l1.size()];
		for (int k = 0; k < res.length; k++) {
			res[k] = l1.get(k);
		}
		return res;
	}

	final AbstractData[] data;

	public void setDefaultCorrVectors(int k) {
		this.correlatedVectors[k].clear();
		for (int j = 0; j < this.data[k].probes.size(); j++) {
			if (!this.data[k].probes.get(j).equals("ones")) {
				correlatedVectors[k].add(this.matsOrig[k].getColumnVector(j));
			}
		}
	}

    public List<String>  makeGreaterThan(int k, int max_ind){
    	//RealMatrix or = matsOrig[k];
    	this.mats[k] = matsOrig[k].copy();
    	RealMatrix mat = mats[k];
    	 
    	for(int i=1; i<mat.getRowDimension(); i++){
    		double[] d = mat.getRow(i);
    		double sum = 0;
    		for(int j= max_ind-1; j>=1; j--){
    			sum+=d[j];
    			if(sum>2){
    				System.err.println("greater than");
    			}
    			mat.setEntry(i, j, sum);
    		}
    		sum = 0;
    		for(int j =max_ind; j<d.length;j++){
    			sum+=d[j];
    			if(sum>2){
    				System.err.println("h");
    			}
    			mat.setEntry(i, j, sum);
    		}
//    		double[] d1 = mat.getRow(i);
  //  		System.err.println("h");
    	}
    	List<String> pr = new ArrayList<String>(data[k].probes);
    	String n1 = pr.get(1).split("__")[0];
    	String nme_ = max_ind == data[k].probes.size()-1 ? "Inf"  : "<="+(pr.get(max_ind).split("__")[1].split("_")[0].replaceAll("<","").replaceAll("=",""));
    	String nme1_ = max_ind == data[k].probes.size()-1 ? "Inf" : (pr.get(max_ind).split("__")[1].split("_")[0].replaceAll("<","").replaceAll("=",""))+"<";
    	for(int i=max_ind+1; i<this.data[k].probes.size(); i++){
    		String nme = i == pr.size()-1 ? "Inf" : pr.get(i).split("__")[1].split("_")[2];
    		pr.set(i,n1+"__"+nme1_+"_C_<="+nme);
    	}
    	for(int i=1; i<max_ind; i++){
    		String nme = pr.get(i).split("__")[1].split("_")[0];
    	 	pr.set(i,n1+"__"+nme+"_C_<="+nme_);
    	}
    	return pr;
    }
    
    
	
	static class Res implements Comparable{
		Double d;
		String st;
		int min,max;
		public Res(Double pv, String string, int min,int max) {
			this.d = pv;
			this.st = string;
			this.min = min;
			this.max = max;
		}
		@Override
		public int compareTo(Object arg0) {
		    return d.compareTo(((Res)arg0).d);
		}
		public void println(PrintWriter pw){
			pw.println(st+"\t"+d);
		}
	}

	
	
  

	private void plot(String num, RealVector realVector,
			RealVector realVector2, RealVector phen) {
		re.assign("x1", realVector.getData());
		String name = "\"out_" + num + ".pdf\"";

		re.assign("y", phen.getData());
		;
		re.eval("ind1 = y==0");
		re.eval("ind2 = y==1");
		REXP rex = re.eval("ind1");
		re.assign("x2", realVector2.getData());
		re.eval("pdf(file=" + name + ")");
		re.eval("plot(x1[ind1],x2[ind1], col=\"red\")");
		;
		re.eval("lines(x1[ind2],x2[ind2], type=\"p\",col=\"green\")");
		;
		re.eval("dev.off()");
	}

	final File dir;

	private void plotCCs(String name, int j) {
		if(!plotCC) return;
		// for(int j=0; j<this.corr.size(); j++){
		// String name = "cca_j+_"+String.format("%5.3g",this.corr.get(j));;

		double ht = 2.5 * (this.mats.length - 1);
		re.eval("pdf(file=\"" + dir1.getAbsolutePath() + "/" + name + ".pdf\", height="
				+ ht + ",width=7)");

		re.eval("par(mfrow = c(" + 2 + ",1),cex=0.3)");
		for (int k = 0; k < mats.length - 1; k++) {
			if (mats[k] == null)
				continue;
			double[] data = this.combinations[k].get(j).getData();
			double min = 0;
			for (int kk = 0; kk < data.length; kk++) {
				if (!Double.isNaN(data[kk]) && data[kk] < min) {
					min = data[kk];
				}
			}
			double[] mins = new double[data.length];
			Arrays.fill(mins, 0);
			re.assign("mins", mins);
			// re.eval("mins[is.na(mins)]=0");
			String[] probes = this.data[k].probes.toArray(new String[0]);
			for (int k1 = 0; k1 < probes.length; k1++) {
				probes[k1] = probes[k1].replaceAll("VNTR", "");
				if (k1 < this.data[k].loc.size()) {
					probes[k1] += "_" + this.data[k].loc.get(k1);
				}
				if (Math.abs(data[k1]) < 0.01)
					probes[k1] = "";
				/*
				 * if(probes[k1].indexOf('<')>=0){ probes[k1] =
				 * probes[k1].split("<")[0]; }
				 */
			}
			/*
			 * for(int kj=0; kj<data.length; kj++){ if(Double.isNaN(data[kj]))
			 * data[kj]=0; }
			 */
			re.assign("weights", data);
			re.eval("weights[is.na(weights)]=0");
			re.assign("probes", probes);
			re.eval("pos = 1:" + data.length + "-0.5");

			re.eval("plot(pos,weights,type=\"p\",xaxt=\"n\",main=\""
					+ this.data[k].type_name + "\")");
			re.eval("text(pos,weights,probes,srt=0)");
		}
		re.eval("dev.off()");
		// }
	}

	public void plotCC1(String name, int j, int partition) {
		// for(int j=0; j<this.corr.size(); j++){
		// String name = "cca_j+_"+String.format("%5.3g",this.corr.get(j));;
		if (re == null)
			this.startEng();
		double ht = (this.correlatedVectors[j].size() - 1) * 0.5;
		re.eval("pdf(file=\"" + name + "\")"
		// +", height="+ht+",width="+ht+")"
				);
		// int[][] inds = null;

		re.eval("par(mfrow = c(" + 2 + ",1),cex=0.3)");
		for (int k = 0; k < correlatedVectors[j].size() - 1; k++) {
			double[] data = this.correlatedVectors[j].get(k).getData();
			double[] data1 = this.correlatedVectors[j].get(k + 1).getData();

			re.assign("x", data);
			re.assign("y", data1);
			re.eval("x[is.na(x)]=0");
			re.eval("y[is.na(y)]=0");
			if (partition < 0) {
				String st = "plot(x,y,col=" + 1 + " ,type=\"p\",xlab=\"corr_"
						+ k + "\",ylab=\"corr_" + (k + 1) + "\")";
				re.eval(st);

			} else if (partition >= 0) {
				RealMatrix m_ = this.matsOrig[partition];

				for (int kk = 1; kk < m_.getColumnDimension(); kk++) {
					double[] d = m_.getColumnVector(kk).mapAdd(
							this.means[partition][kk]).getData();
					re.assign("tmp", d);
					re.eval("tmp[is.na(tmp)]=0");
					re.eval("inds= tmp!=0");
					String st = "plot(x[inds],y[inds],col=" + (kk)
							+ " ,type=\"p\",xlab=\"corr_" + k
							+ "\",ylab=\"corr_" + (k + 1) + "\")";
					re.eval(kk == 1 ? st : st.replace("plot", "lines"));
				}
			}
		}
		re.eval("dev.off()");
		// }
	}

	//final PrintWriter pw;

	
	final String header  = "namek\tnamej\tbeta\tse\tpv";

	static Rengine re;

	public void startEng() {
		re = new Rengine(new String[] { "--vanilla" }, true, 
			//	null
				new LoopCall()
				);
		
		File f = new File("funcs.R");
		if(false && f.exists()){
			try{
			BufferedReader br = new BufferedReader(new FileReader(f));
			String st = "";
			while((st = br.readLine())!=null){
				re.eval(st);
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		re.eval("library(MASS)");
	if(false)	re.eval("library(stepPlr)");
	
	}
	public double[] getRow(REXP r, int[] dim, double nonZero, double thresh){
		//RealMatrix r1 = new Array2DRowRealMatrix(dim[0], dim[1]);
		double[] r1=null;
		int cnt = 0;
		double[] lambda = re.eval("fit$lambda").asDoubleArray();
		double[] r_prev =  null;
		for(int k=0; k<dim[1]; k++){
			 cnt=0;
			int k1 = k+1;
			r_prev = r1;
			 r1 = re.eval("fit$beta[,"+k1+"]").asDoubleArray();
			for(int j=0; j<r1.length; j++){
				if(Math.abs(r1[j])>thresh){
					cnt++;
				}
			}
			if(cnt>=nonZero){
				if(true){
					return r1;
				}
				if(r_prev!=null) {
					return r_prev;
				}
				return r1;
			}
		}
//		int k1 = dim[1];
		// r1 = re.eval("fit$beta[,"+k1+"]").asDoubleArray();
		return r1;
	}

	//final double[] q;
	//final double[] q1;

	/* This method calculates the CC decomposition between the matrices mat[0] and mat[1], and also
	 * calculates the linear combs maximally correlated with the combination derived from mat[0]
	 */
	public void decompose(boolean assoc, double thresh1) throws Exception {

		for (int i = 0;; i++) {//iteration continues to find new correlated directions
			
			while (reg() > thresh) {}  //reg() finds maximal correlations between mat[0] and mat[1], 
			for (int k = 2; k < this.mats.length; k++) {
				if (data[k] == null) continue;
				//once we have found that, we should regress the other matrices (2 ... length) on the first matrix
				this.regress(k, 0);  
			}
			//System.err.println(cor);
			if (Double.isNaN(cor.getEntry(0, 1)))
				break;
			    cor = computeCorrelationMatrix(this.y);
			    updateVectors();
			//	
			if (assoc) {
			

				this.plotCCs("corr_" + q , covar.size() - 1);

			}
			RealMatrix corr1 = p1.covarianceToCorrelation(covar.get(covar.size() - 1));
			double corrd = corr1.getEntry(0, 1);
			if (corrd < thresh1 || Double.isNaN(corrd))
				break;
			System.err.println("correlation "+i+" "+corr1.getEntry(0, 1));
			//System.err.println(this.covar.get(covar.size() - 1).getEntry(0, 1));

			adjust();

		}

	}
/** this method used when we want to keep the independent dimensions of mat[k] fixed */
	private void decomposeFixFirst(boolean assoc) throws Exception {
		double[] d = new double[mats[0].getColumnDimension()];
		for (int i = 0; i < this.mats[0].getColumnDimension(); i++) {
			if (this.data[0].probes.get(i).equals("ones"))
				continue;
			Arrays.fill(d, 0);
			d[i] = 1;
			combs[0] = new Array2DRowRealMatrix(d);
			y.setColumn(0, mats[0].multiply(combs[0]).getColumn(0));
			try {
				for (int k = 1; k < this.mats.length; k++) {
					if (mats[k] != null)
						this.regress(k, 0);
				}
				cor = computeCorrelationMatrix(this.y);
				updateVectors();
				if (assoc) {
					
					this.plotCCs("corr_" + q , covar.size() - 1);
				}
			} catch (Exception exc) {
				exc.printStackTrace();
			}
			// if(corr.get(corr.size()-1).getEntry(1)<0.2) break;
			// System.err.println(this.corr.get(corr.size()-1).getEntry(0, 1));
			// adjust();

		}

	}

	public void associate() {
		for (int k = 0; k < covar.size(); k++) {
			try{
			//String name = associate(k);
			this.plotCCs("corr_" + q[0], k);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
	}

	private void adjust() {
		for (int k = 0; k < mats.length; k++) {
			adjust(k);
		
			combs[k] =this.getRandomCombs(mats[k].getColumnDimension(), this.data[k].probes.get(0).equals("ones"));
		}
		this.updateY();
		// TODO Auto-generated method stub

	}

	final boolean[] modified;

	private void updateVectors() {

		this.covar.add(this.cor);

		for (int k = 0; k < this.mats.length; k++) {

			if (mats[k] == null)
				continue;
			RealVector rv = this.y.getColumnVector(k).copy();
		//double norm1 = rv.getNorm();
			//check that rv is genuinely orthogonal to all other vectors
			for (int j = 0; j < correlatedVectors[k].size(); j++) {
				double dp = rv.dotProduct(correlatedVectors[k].get(j));
				if ( k < 2 && dp > 1e-3) {
					throw new RuntimeException("!!");
				}
				// System.err.println("dot product "+dp);
			}
			this.correlatedVectors[k].add(rv);
			boolean onesk = this.data[k].probes.get(0).equals("ones");
			
		  //now, we have found the linear combination, but we need to know the combination in terms
		// or original matrix, not updated matrix.  next section calculates that
			RealVector v1;
			if (!modified[k]) {
				RealVector v = this.combs[k].getColumnVector(0);
				double norm = normalise ? v.getNorm() : 1.0;
				v1 = v.mapDivide(norm);
			} else {
			//	if(true)throw new RuntimeException("!!!");
				RealMatrix A_ = matsOrig[k];
				RealMatrix b_ = y.getColumnMatrix(k);
				RealMatrix A, b;
				int[] col_A = this.cols[k];
				int[] col_b = single;
				int[] rows = null;
				if (col_A != null) {
					rows = CC.getNonNARow(A_, b_);
					A = A_.getSubMatrix(rows, col_A);
					b = b_.getSubMatrix(rows, col_b);
				} else {
					A = A_;
					b = b_;
				}

				RealMatrix beta1 = solveExact(binary[k],onesk, A, b, 1e-10);
				RealVector resid = A.multiply(beta1).subtract(b).getColumnVector(0);
				//double resnorm = resid.getNorm();
				double resnorm1 = resid.getNorm();
				if(resnorm1>1e-3){
					throw new RuntimeException("did not find direction in orig matrix "+resnorm1);
				}
				RealVector v = beta1.getColumnVector(0);
				double norm = normalise ? v.getNorm() : 1.0;
				v1 = v.mapDivide(norm);
			}
			this.combinations[k].add(v1);
		}

	}

	private double reg() throws Exception {
		double corst = p1.covarianceToCorrelation(cor).getEntry(0, 1);
		for (int k = 0; k < 2; k++) {
			this.regress(k, 1 - k);
			// System.err.println(this.cor);
		}
		cor = computeCorrelationMatrix(this.y);
		double corend = p1.covarianceToCorrelation(cor).getEntry(0, 1);
		double diff = corend - corst;
		if (Double.isNaN(diff))
			return 0;
		//System.err.println(corend);
		return diff;
	}

	private RealMatrix computeCorrelationMatrix(RealMatrix y2) {
		RealMatrix res = new Array2DRowRealMatrix(y2.getColumnDimension(), y2
				.getColumnDimension());
		int[] col = new int[] { 0, 0 };
		for (int k = 0; k < y2.getColumnDimension() - 1; k++) {
			RealMatrix res1 = y2.getColumnMatrix(k);
			col[0] = k;
			for (int j = k + 1; j < y2.getColumnDimension(); j++) {
				col[1] = j;
				RealMatrix res2 = y2.getColumnMatrix(j);
				int[] row = getNonNARow(res1, res2);
				if (row.length > 0) {  
					RealMatrix sub = y2.getSubMatrix(row, col);
					Covariance p2 = new Covariance(sub);
					RealMatrix corr = p2.getCovarianceMatrix();//p2.computeCovarianceMatrix(sub);//
					   // p1.computeCorrelationMatrix(sub);
					// this.checkNA(corr);
					for (int jk = 0; jk < 2; jk++) {
						for (int jk2 = 0; jk2 < 2; jk2++) {
							res.setEntry(col[jk], col[jk2], corr.getEntry(jk,
									jk2));
						}
					}
				} else {
					for (int jk = 0; jk < 2; jk++) {
						for (int jk2 = 0; jk2 < 2; jk2++) {
							res.setEntry(col[jk], col[jk2], Double.NaN);
						}
					}
				}
			}
		}
		return res;
	}

	public void adjust(int k) {
		RealMatrix mat = mats[k];
		RealVector basis = this.correlatedVectors[k].get(correlatedVectors[k]
				.size() - 1);
		
		double denom = dotProduct(basis, basis);
		
		for (int j = 0; j < mat.getColumnDimension(); j++) {
			if (!this.data[k].probes.get(j).equals("ones")) {
				RealVector vec = mat.getColumnVector(j);
				// for(int k1=0; k1<basis.size(); k1++){
				// vec.d
				double num = dotProduct(basis, vec);
				// }
				double ratio = num / denom;
				RealVector yind = mat.getColumnVector(j).subtract(
						basis.mapMultiply(ratio));
				this.modified[k] = true;
				// for(int k1=0; k1<basis.size(); k1++){
				// yind = yind
				// }
				// double norm = yind.getNorm();
				mat.setColumnVector(j, yind);// .mapDivide(yind.getNorm()));
			}
		}
		// this.centralise(this.data[k].probes, this.mats[k],true,true);
		// }
	}

	private static double dotProduct(RealVector basis, RealVector basis2) {
		double res = 0;
		for (int k = 0; k < basis.getDimension(); k++) {
			double v1 = basis.getEntry(k);
			double v2 = basis2.getEntry(k);
			if (!Double.isNaN(v1) && !Double.isNaN(v2)) {
				res += v1 * v2;
			}
		}
		return res;
	}

//	static double thresh = 0.0001;
	public static Double lowerB = null;

	List<String> indiv;
	RealMatrix[] mats;
	RealMatrix[] matsOrig;
	// boolean[][][] na;
	RealMatrix[] combs;
	RealMatrix y;

	int pheno_index;
 int[][] cols;
	final int[] single = new int[] { 0 };

	int[][] alias;

	

	//boolean[] restrictIndiv;

	
	public static boolean restrictToCommonIds = true;
	public void initialise(){
		this.correlatedVectors = new List[data.length];
		this.combinations = new List[data.length];
		this.covar = new ArrayList<RealMatrix>();
		Arrays.fill(modified, false);
		//int i = 0;
		//for (; data[i] == null; i++) {
		//}
		indiv = new ArrayList<String>();
		for (int i1 = 0; i1 < data.length; i1++) {

			if (data[i1] != null)
			// && !restrictIndiv[i1])
			{
				for (Iterator<String> it = data[i1].indiv.iterator(); it
						.hasNext();) {
					String nxt = it.next();
					if (!indiv.contains(nxt)) {
						indiv.add(nxt);
					}
				}
			}
		}
		if (restrictToCommonIds) {
			for (int i1 = 0; i1 < 2; i1++) { //restrict indiv list to indiv common to first two datasets
				if (data[i1] != null ) {
					System.err.println("before " + indiv.size());
					Set<String> l1 = new HashSet<String>(indiv);
					indiv.retainAll(data[i1].indiv);
					l1.removeAll(new ArrayList<String>(data[i1].indiv));
					if (l1.size() > 0) {
						System.err.println(data[i1].name + " not contained "
								+ l1);
					}
					System.err.println("after " + indiv.size());
				}
			}
		}
		alias = new int[data.length][indiv.size()];
		mats = new RealMatrix[data.length];
		this.cols = new int[mats.length][];
		this.matsOrig = new RealMatrix[data.length];
		combs = new RealMatrix[data.length];

		// int[][] alias = getCombinedList(data,this.indiv = new
		// ArrayList<String>());
		
		for (int k = 0; k < mats.length; k++) {
			if (data[k] == null)
				continue;

			correlatedVectors[k] = new ArrayList<RealVector>();
			combinations[k] = new ArrayList<RealVector>();

			for (int j = 0; j < indiv.size(); j++) {
				alias[k][j] = data[k].indiv.indexOf(indiv.get(j));
			}
			int[] col = new int[data[k].data.getColumnDimension()];
			for (int j = 0; j < col.length; j++) {
				col[j] = j;
			}
			cols[k] = col;

			matsOrig[k] =
			restrictToCommonIds ?
			 data[k].data.getSubMatrix(alias[k], col) :
			expand(data[k].data, alias[k], col);
			mats[k] = // matsOrig[k];
			//
			matsOrig[k].copy();
			
			combs[k] = getRandomCombs(mats[k].getColumnDimension(),this.data[k].probes.get(0).equals("ones"));
			
			// .viewSelection(alias, col);
			int ind = -1;
			if (mats[k].getColumnDimension() == 0)
				ind = 0;
			else if (mats[k].getColumnDimension() == 2
					&& data[k].probes.get(0).equals("ones"))
				ind = 1;
			if (ind >= 0) {
				Set<Double> s = new HashSet<Double>();
				for (int j = 0; j < mats[k].getRowDimension(); j++) {
					double v = mats[k].getEntry(j, ind);
					if (!Double.isNaN(v)) {
						s.add(v);
					}
				}
				binary[k] = s.size() <= 2;
			} else {
				binary[k] = false;
			}
		}

		y = new Array2DRowRealMatrix(indiv.size(), data.length);
		this.updateY();
		this.makeGamma(q,  band);
	}
	/*get random initial combinations, close to uniform, but not exactly uniform*/
	private RealMatrix getRandomCombs(int columnDimension,boolean ones) {
		double[] d = new double[columnDimension];
		Arrays.fill(d, 1);
		if (ones) {
			d[0] = 0;
		}
		try{
		 d = (new Dirichlet(d,u)).sample();
		}catch(MathException exc){
			exc.printStackTrace();
		}
		
		RealMatrix combsk = new Array2DRowRealMatrix(d);
		double norm = combsk.getNorm();
		combsk.scalarMultiply(1 / norm);
		return combsk;
	}


	File dir1;
	public CC(AbstractData[] data, File dir, 
			final Double[] q,  final int[] band)
			throws Exception {
		this.data = data;
		//this.q1 = q1;
		this.q = q;
		this.modified = new boolean[data.length];
		this.band = band;
		
		this.dir = dir;
		this.pheno_index = -1;
		this.baf_index = -1;
		this.binary = new boolean[data.length];
		for (int k = 0; k < data.length; k++) {
			if (data[k] == null)
				continue;
			if (data[k].type_name == "baf") {
				this.baf_index = k;
			}
			if (data[k].type_name == "pheno") {
				this.pheno_index = k;
			}
			if (data[k].type_name == "lrr") {
				this.lrr_index = k;
			}
		}
		// this.band = band;
		dir2.mkdir();
		dir1 = new File(dir2,dir.getName());
		dir1.mkdir();
		this.deletePdfs(dir1);
		
		this.initialise();
	}

	private static RealMatrix expand(RealMatrix data2, int[] alias, int[] col) {
		if (alias.length == data2.getRowDimension())
			return data2;
		RealMatrix data = new Array2DRowRealMatrix(alias.length, col.length);
		for (int k = 0; k < col.length; k++) {
			for (int j = 0; j < alias.length; j++) {
				if (alias[j] < 0) {

					data.setEntry(j, col[k], col[k] == 0 ? 1 : Double.NaN);
				} else {
					data.setEntry(j, col[k], data2.getEntry(alias[j], col[k]));
				}
			}
		}
		// TODO Auto-generated method stub
		return data;
	}

	private static RealMatrix expand1(RealMatrix data2, int[] alias, int[] col,
			int len) {
		// if(alias.length==data2.getRowDimension()) return data2;
		RealMatrix data = new Array2DRowRealMatrix(len, col.length);

		for (int k = 0; k < col.length; k++) {
			for (int j = 0; j < len; j++) {
				data.setEntry(j, col[k], Double.NaN);
			}
			for (int j = 0; j < alias.length; j++) {
				data.setEntry(alias[j], col[k], data2.getEntry(j, col[k]));
			}
		}
		// TODO Auto-generated method stub
		return data;
	}

	private static RealVector expand1(RealVector data2, int[] alias, int len) {
		// if(alias.length==data2.getRowDimension()) return data2;
		RealVector data = new ArrayRealVector(len);

		// for(int k=0; k<col.length; k++){
		for (int j = 0; j < len; j++) {
			data.setEntry(j, Double.NaN);
		}
		for (int j = 0; j < alias.length; j++) {
			data.setEntry(alias[j], data2.getEntry(j));
		}
		// }
		// TODO Auto-generated method stub
		return data;
	}

	/*
	 * private RealMatrix getOneMatrix(RealMatrix matsOrig2) { RealMatrix data =
	 * new Array2DRowRealMatrix(matsOrig2.getRowDimension(),
	 * matsOrig2.getColumnDimension()+1); for(int k=1;
	 * k<data.getColumnDimension(); k++){ data.setColumnMatrix(k,
	 * matsOrig2.getColumnMatrix(k-1)); } for(int j=0;
	 * j<matsOrig2.getRowDimension(); j++){ matsOrig2.setEntry(j, 0, 1.0); }
	 * return data; }
	 */

	public void updateY() {
		for (int k = 0; k < combs.length; k++) {
			if (mats[k] != null) {
				y.setColumn(k, mats[k].multiply(combs[k]).getColumn(0));
			}
		}
		cor = computeCorrelationMatrix(this.y);
		// cor = corr.getEntry(0,1);
		// System.err.println("h "+cor);
	}

	RealMatrix cor;
	PearsonsCorrelation p1 = new PearsonsCorrelation();
	
	NormalDistribution dist = new NormalDistributionImpl();

	public void regress1(int k) throws Exception {
		OLSMultipleLinearRegression ols = new OLSMultipleLinearRegression();
		double[] y1 = y.getColumn(1 - k);

		ols.newSampleData(y1, mats[k].getData());
		double[] params = ols.estimateRegressionParameters();
		double[] beta = ols.estimateRegressionParameters();
		double[] se = ols.estimateRegressionParametersStandardErrors();
		double[] pv = new double[se.length];
		RealMatrix beta1 = new Array2DRowRealMatrix(beta);
		RealMatrix y = mats[k].multiply(beta1);
		double[] resid = ols.estimateResiduals();
		for (int j = 0; j < beta.length; j++) {
			dist.setStandardDeviation(se[j]);
			pv[j] = 2 * (1 - dist.cumulativeProbability(Math.abs(beta[j])));

		}
		this.combs[k] = beta1;
		this.y.setColumn(k, y.getColumn(0));
		// this.y[k].
		// ols.
		cor = computeCorrelationMatrix(this.y);
		// cor = corr.getEntry(0,1);
		// System.err.println("h "+cor);
	}

	// final int[] band;
	final boolean[] binary;
	final int[] band;

	/*This method uses the k1 column of y (representing the current linear combination of mat k1)
	*then finds combinations of mat[k] which are maximally correlated
	*/
	//This uses the fixed columns of mat[k] to define a variable, and then regresses mat[k1] on this. 
	public void regress(int k, int k1) throws Exception {
		RealMatrix beta1;
		boolean onesk = this.data[k].probes.get(0).equals("ones"); //is the first column just a 'ones' column
		RealMatrix A_ = mats[k];
		RealMatrix b_ = y.getColumnMatrix(k1);
		RealMatrix A, b;
		int[] col_A = this.cols[k];
		int[] col_b = single;
		int[] rows = null;
		if (col_A != null) {
			rows = CC.getNonNARow(A_, b_);
			A = A_.getSubMatrix(rows, col_A);
			b = b_.getSubMatrix(rows, col_b);
		} else {
			A = A_;
			b = b_;
		}
		if (onesk && mats[k].getColumnDimension() <= 2 ) {
			//this if there is only one column in mats[k] in which case answer is trivial
			beta1 = new Array2DRowRealMatrix(new double[][] {
					new double[] { 0 }, new double[] { 1 } });
		} else if(!onesk && mats[k].getColumnDimension()<=1){
			beta1 = new Array2DRowRealMatrix(new double[][] {
					 new double[] { 1 } });
		}else {
			beta1 = 	solve(binary[k1], onesk, A, b, k);
			RealVector v = beta1.getColumnVector(0);
			if (onesk) {
				v.setEntry(0, 0);
			}
			if(true){
				double norm = v.getNorm();
				if (norm > 0) {
					v.mapDivideToSelf(norm);
				}
			}
			//int m = getMax(v.getData());
			//double max = v.getEntry(m);
			beta1.setColumnVector(0, v);
		}
		RealMatrix y_ = A.multiply(beta1);
		double norm = y_.getNorm();
		if (Double.isNaN(norm)) {
			checkNA(beta1);
			checkNA(mats[k]);
			// throw new RuntimeException("!!");
		}
		// y_ = y_.scalarMultiply(1.0/norm);
		this.combs[k] = beta1; //set the combinations of k
		if (rows != null) {
			// this would be if there were NA rows
			y_ = expand1(y_, rows, this.single, this.y.getRowDimension());
		}
		this.y.setColumn(k, y_.getColumn(0));  //update y
	}

	private int getMax(double[] data2) {
		int max_ind = 0;
		for(int k=0; k<data2.length; k++){
			if(data2[k]>data2[max_ind]){
				max_ind = k;
			}
		}
		return max_ind;
	}

	private void checkNA(RealMatrix beta1) {
		for (int j = 0; j < beta1.getRowDimension(); j++) {
			for (int k = 0; k < beta1.getColumnDimension(); k++) {
				double v = beta1.getEntry(j, k);
				if (Double.isNaN(v)) {
					throw new RuntimeException("!!");
				}
			}
		}

	}

	private boolean binary(RealMatrix realMatrix, int i) {
		Double v1 = null;
		Double v2 = null;
		for (int k = 0; k < realMatrix.getRowDimension(); k++) {
			double v = realMatrix.getEntry(k, i);
			if (!Double.isNaN(v)) {
				if (v1 == null)
					v1 = v;
				else if (Math.abs(v1.doubleValue() - v) > 0.01) {
					if (v2 == null)
						v2 = v;
					else if (Math.abs(v2.doubleValue() - v) > 0.01) {
						return false;
					}
				}
			}
		}
		return true;
	}

	

	
	
	
	
    /*just returns the rows of the matrix which are NA */
	public static int[] getNonNARow(RealMatrix A_, RealMatrix b_) {
		List<Integer> l = new ArrayList<Integer>();
		outer: for (int k = 0; k < A_.getRowDimension(); k++) {
			for (int j = 0; j < A_.getColumnDimension(); j++) {
				if (Double.isNaN(A_.getEntry(k, j)))
					continue outer;
			}
			for (int j = 0; j < b_.getColumnDimension(); j++) {
				if (Double.isNaN(b_.getEntry(k, j)))
					continue outer;
			}
			l.add(k);
		}
		int[] rows = new int[l.size()];
		for (int k = 0; k < rows.length; k++) {
			rows[k] = l.get(k);
		}
		return rows;
	}

	

	
	
	/*basic method for calculating the regression */
	public  RealMatrix solve(boolean binaryPhenotype, boolean onesk, RealMatrix A, RealMatrix b,int k
			) {
		RealMatrix Q = this.Q[k];
		RealMatrix AT = A.transpose();
		RealMatrix prod = AT.multiply(A);//.add(Q);
		if(Q!=null) prod = prod.add(Q);
		else{
			for(int i=0; i<prod.getColumnDimension(); i++){
				prod.addToEntry(i, i, q[k]);
			}
		}
		//RealMatrix b1 = b.copy();
		LUDecompositionImpl lu = new LUDecompositionImpl(prod);
		DecompositionSolver solver = lu.getSolver();
		RealMatrix res = solver.solve(AT.multiply(b));
		return res;
	}

	
	public  RealMatrix solveExact(boolean binaryPhenotype, boolean onesk, RealMatrix A, RealMatrix b,double q
	) {

RealMatrix AT = A.transpose();
RealMatrix prod = AT.multiply(A);//.add(Q);
/*RealMatrix Q = this.Q[k];
if(Q!=null) prod = prod.add(Q);
else{*/
	for(int i=0; i<prod.getColumnDimension(); i++){
		prod.addToEntry(i, i, q);
	}

//RealMatrix b1 = b.copy();
LUDecompositionImpl lu = new LUDecompositionImpl(prod);
DecompositionSolver solver = lu.getSolver();
RealMatrix res = solver.solve(AT.multiply(b));
return res;
}

	

}
