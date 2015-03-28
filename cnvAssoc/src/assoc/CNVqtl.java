package assoc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.regression.OLSMultipleLinearRegression;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class CNVqtl {
	//comment this out if you don't have JRI installed
	
	final String src_chrom, target_chrom;
	public static void main(String[] args){
		try{
			Constants.parse(args, Constants.class);
			File dir = new File(System.getProperty("user.dir"));
			File user = new File(Constants.baseDir);
			String probe = Constants.probe;
			File src = new File(user, Constants.src+".zip");
			File target = new File(user,Constants.target+".zip");
			File covariates = new File(dir,Constants.covariates);
			File output = new File(dir,probe+"_"+Constants.target+"_"+Constants.distance+"_"+Constants.results+".txt");
			CNVqtl cnvqtl = new CNVqtl(src,target, covariates, output);
			cnvqtl.setSrcId(probe);
			List<String> probes =  new ArrayList<String>();
			List<Integer> loc =  new ArrayList<Integer>();
			cnvqtl.srcLoc = 	cnvqtl.getProbesNearSrc(Constants.process(Constants.distance), probes, loc);
			for(int k=0; k<probes.size(); k++){
				cnvqtl.calcAssociation(probes.get(k), loc.get(k));
			}
			cnvqtl.close();
			System.exit(0);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	



	String probeSrc;
	int srcLoc;
	final ZipFile srcF;
	
	final ZipFile targetFile;
	
	final DoubleMatrix2D src; //contains LRR
	final DoubleMatrix2D target; 
	//first column is all 1s to capture the intercept
	//second column is genotypes
	//remaining columns are covariates

	
	final int srcColumnInd, targetColumnInd;
	
	final Set<Integer> nonNa_indSrc = new TreeSet<Integer>();
	final Set<Integer> nonNa_indTarget = new TreeSet<Integer>();
	
	final int[] cols_target, cols_src;  
	
	final List<String> covars;
	final List<String> indiv;
	String assocline = "y ~ x";
	boolean isGenotypeTarget;
	final PrintWriter pw;
	//Assume first column of covariates is patient id, but also assume it is in same order as 'Samples' entry in zip
	
	public void close() throws Exception{
		this.srcF.close();
		this.targetFile.close();
		this.pw.close();
	}
	
	public CNVqtl(File inSrc, File inTarget, File covariates, File output) throws Exception{
		this.srcF = new ZipFile(inSrc);
		this.src_chrom = inSrc.getName().split("\\.")[0];
		this.target_chrom = inTarget.getName().split("\\.")[0];
		this.pw = new PrintWriter(new BufferedWriter(new FileWriter(output)));
		this.targetFile = new ZipFile(inTarget);
		srcColumnInd = getIndex(srcF,"Name",1,"Log R Ratio:Log R".split(":") );
		targetColumnInd = getIndex(srcF,"Name",1,"B Allele Freq".split(":"));
				//"GType:Genotype".split(":") );
		isGenotypeTarget = false;
		
		indiv = read(srcF, "Samples",0);
		if(!indiv.equals(read(targetFile, "Samples",0))){
			throw new RuntimeException("problem -Sample files do not agree");
		}
		this.src = new DenseDoubleMatrix2D(1,indiv.size());
		
		if(covariates!=null && covariates.exists()){
		BufferedReader br = new BufferedReader(new FileReader(covariates));
		List<String> head = Arrays.asList(br.readLine().split("\t"));
		covars = head.subList(1, head.size());
			for(int k=0; k<covars.size(); k++){
				assocline+="+cov_"+k;
			}
			br.close();
		}
		else{
			covars = new ArrayList<String>();
		}
		this.target = new DenseDoubleMatrix2D(indiv.size(), 2+covars.size() );
		for(int k=0; k<target.rows(); k++){
			target.set(k, 0, 1);
		}
		cols_target = getColIndex(target.columns());
		cols_src = getColIndex(src.rows());
		pw.println(header);
		
	}
	
	private int[] getColIndex(int columns) {
		int[] res= new int[columns];
		for(int k=0; k<res.length; k++){
			res[k] = k;
		}
		return res;
	}

	
	
	private int getIndex(ZipFile zf, String probe, int i, String[] string2) throws Exception {
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		for(int k=0; k<i; k++){
			st = br.readLine();
		}
		br.close();
		int ind =-1;
		List<String> str = Arrays.asList( st.split("\t"));
		for(int ij=0; ij<str.size(); ij++){
			str.set(ij, str.get(ij).trim());
		}
		for(int k=0; k<string2.length; k++){
			ind = str.indexOf(string2[k]);
			if(ind>=0) return ind;
		}
		return ind;
		
	}
	
	private int getProbesNearSrc(int len, List<String> res, List<Integer> loc) throws Exception {
		String probe = "SNPS";
		int pos =-1;
		{
			ZipFile zf = this.srcF;
			BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
			String st = "";
		
			for(int k=0;(st = br.readLine())!=null; k++){
				String[] str = st.split("\\s+");
				if(str[3].equals(this.srcId)){
					 pos = Integer.parseInt(str[1]);
					 break;
				}
			}
			br.close();
		}
		if(pos>=0){
			ZipFile zf = this.targetFile;
			BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
			String st = "";
			
			for(int k=0;(st = br.readLine())!=null; k++){
				String[] str = st.split("\\s+");
				if(Math.abs(Integer.parseInt(str[1])-pos)<len){
					res.add(str[3]);
					loc.add(Integer.parseInt(str[1]));
				}
			}
			br.close();
		}
		
		return pos;
	}
	
	private List<String> read(ZipFile zf,  String probe, int i) throws Exception {
		List<String> res = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st.split("\t")[0]);
		}
		br.close();
		return res;
		
	}

	public void fillMatrix(DoubleMatrix2D mat, ZipFile zf, String probe, int columnIndInZip, int colIdMatrix,
			Collection<Integer> nonNa_ind, boolean isGenotype, boolean transpose
	) throws Exception{
		ZipEntry entry = zf.getEntry(probe);
		if(entry==null) throw new RuntimeException("no entry for "+probe);
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(entry)));
		//BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		for(int i=0; (st = br.readLine())!=null; i++){
			String str = st.split("\t")[columnIndInZip];
			
			double v = 
				isGenotype ?  ( str.indexOf('N')>=0 ? Double.NaN : count(str.toCharArray(),'B')) : 
				Double.parseDouble(str);
			
			if(transpose)
				mat.set( colIdMatrix,i, v);
			else
				mat.set(i, colIdMatrix, v);
			if(!Double.isNaN(v)){
				nonNa_ind.add(i);
			}
		}
		br.close();
	}
	private double count(char[] charArray, char c) {
		int cnt=0;
		for(int k=0; k<charArray.length; k++){
			if(c==charArray[k]) cnt++;
		}
		return cnt;
	}





	String srcId = "";
	//call this to change the src matrix of lrr values
	public void setSrcId(String rsid) throws Exception{
		nonNa_indSrc.clear();
		this.srcId = rsid;
		fillMatrix(src, srcF, rsid, srcColumnInd, 0,nonNa_indSrc, false, true);
	}
	
	public double[] solve(DoubleMatrix2D Y_, DoubleMatrix2D X_) throws Exception{
		DoubleMatrix2D AT = alg.transpose(X_);
		DoubleMatrix2D prod = alg.mult(AT, X_);
		DoubleMatrix2D B =  alg.solve(prod,alg.mult(AT,Y_));
	  
	  	DoubleMatrix2D Yha = alg.mult(X_, B);
    	double S = 0;
    	
    	//double mean = B.get(0, 0);
    	
    	//double var0 =0;
    	for(int i1=0; i1<Y_.rows(); i1++){
    		double y1 = Y_.getQuick(i1, 0);
    		double y2 = Yha.getQuick(i1, 0);
    		double p1 = Math.pow(y1-y2,2);
	    	  S+=p1;
	    }
    	 
    	
         double p = X_.columns();
         double n =  X_.rows();
         double var = Math.sqrt(S/(double)n);
         double   se = Math.sqrt((S/(n-p)) * alg.inverse(prod).getQuick(1,1));
         double  beta = B.getQuick(1, 0);
         dist.setStandardDeviation(se);
     	double pv = 2*(1 - dist.cumulativeProbability(Math.abs(beta)));
         return new double[] {beta,se,pv};
	}
	OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
	NormalDistribution dist = new NormalDistributionImpl();
	public double[] solve1(DoubleMatrix2D Y_, DoubleMatrix2D X_) throws Exception{
		double[][] y = Y_.toArray();
		double[][] x = X_.toArray();
	
		regression.newSampleData(y[0], x);
		double beta = regression.estimateRegressionParameters()[1];
		double se = regression.estimateRegressionParametersStandardErrors()[1];
		dist.setStandardDeviation(se);
		double pv = 2*(1 - dist.cumulativeProbability(Math.abs(beta)));
		return new double[] {beta,se,pv};
	}
	
	/*
	 * This runs the R
	Rengine re=new Rengine(new String[] {"--vanilla"}, false, null);
public double[] solve2(DoubleMatrix2D Y_, DoubleMatrix2D X_){
	double[][] y = Y_.toArray();
	double[][] x = alg.transpose(X_).toArray();
		re.assign("y", y[0]);
		re.assign("x", x[1]);
		for(int k1=2; k1<x.length; k1++){
			re.assign("cov_"+k1, x[k1]);
		}
		re.eval("summ= summary(lm("+assocline+"))");
        REXP summ = re.eval("summ");
        RList names= re.eval("dimnames(summ$coefficients)").asList();
        double[][]coeff = summ.asList().at("coefficients").asMatrix();
        int indrow = Arrays.asList(names.at(0).asStringArray()).indexOf("x");
 //       String[] names1 = names.at(1).asStringArray();
        if(indrow>=0 && indrow<coeff.length){
        	double beta   = coeff[indrow][0];
       
      double se = coeff[indrow][1];
      double pv =  coeff[indrow][3];
      return new double[] {beta, se,pv};
        }
        return new double[0] ;//{beta, se};
	}
	*/
	public void calcAssociation(String rsid, int rsLoc) throws Exception{
		nonNa_indTarget.clear();
		//na_indTarget.addAll(na_indSrc);
		fillMatrix(target, targetFile, rsid, targetColumnInd, 1, nonNa_indTarget, isGenotypeTarget,false);
		nonNa_indTarget.retainAll(nonNa_indSrc);
	  	int[] row_ind = getInd(nonNa_indTarget);
	  	DoubleMatrix2D Y_ = src.viewSelection(cols_src,row_ind);  //row matrix
	  	DoubleMatrix2D X_ = target.viewSelection(row_ind, cols_target);
		try{
		double[] res1 = solve1(Y_,X_);
		//two different ways to solve
		//double[] res = solve(alg.transpose(Y_), X_);
		//double[] res2 = solve2(Y_,X_);
		pw.println(this.srcId+"\t"+rsid+"\t"+this.src_chrom+"\t"+this.target_chrom+"\t"+srcLoc+"\t"+rsLoc+"\t"+String.format("%5.3g",res1[0]).trim()+"\t"+String.format("%5.3g",res1[1]).trim()+"\t"+String.format("%5.3g",res1[2]).trim());
		}catch(Exception exc){
			System.err.println("prob at "+rsid+" "+exc.getMessage());
		}
		//System.err.println("h");
	}
	String header = "CNV_Probe_id\tGeno_probe_id\tCNV_Probe_chrom\tGeno_chrom\tCNV_loc\tGeno_loc\tbeta\tse\tpvalue";
	private int[] getInd(Collection<Integer> inds){
		int[] res = new int[inds.size()];
		int k=0;
		for(Iterator<Integer> it = inds.iterator(); it.hasNext();k++){
			res[k] = it.next();
		}
		return res;
	}
	static Algebra alg = new Algebra();
	
	
	
	//remove NA rows
	
	
	
	
	
}
