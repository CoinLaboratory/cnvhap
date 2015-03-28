package lc1.sfs;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public abstract class SFS {

	public static final boolean CHECK = true;
	
	
	protected AbstractMultiDM AFS;
	protected int[] pop,  lim;
	protected RealMatrix[] hweprobs;

	protected ArrayIt it;
	List<String> nme;
	final int[] colsToInclude = new int[] {0,1,2};

	//public static boolean useMaxGeno = false;
	

	ProbMap[] pmap;

	
	public class Likelihood{ //this calculates P(D | p)  for each population
	//	RealVector prob;
	//	RealVector counts;
		
		RealVector lhood;  //P(D | p)  = lhood * exp(logScale)
	
		
		double logScale = 0.0;
		
		Likelihood(int len){
			lhood = new ArrayRealVector(len);
			lhood.set(1);
		}
		
		Likelihood(Likelihood lh){
			this.lhood = lh.lhood.copy();
			this.logScale = lh.logScale;
		}
		public Likelihood clone(){
			return new Likelihood(this);
		}
		
		public void setData(RealMatrix lh1, RealMatrix hweprobs12,int[] rowsToInclude){
			RealMatrix subm = lh1.getSubMatrix(rowsToInclude, colsToInclude);
			
			RealMatrix  rm1 = subm.multiply(hweprobs12); //each row is P(D_i | G)
			lhood.set(1);
			logScale = 0;
			for(int k=0; k<rm1.getRowDimension(); k++){
				RealVector rv = rm1.getRowVector(k);
				//RealVector rv = likeli.ebeMultiply(prob);
				double norm = rv.getL1Norm();
				lhood = lhood.ebeMultiply(rv);
				double sc = lhood.getLInfNorm();
				this.logScale+=Math.log(sc);
				lhood.mapDivideToSelf(sc); //keeps probs within suitable range (not too small)
				rv.mapDivideToSelf(norm);
				//this.addCounts(rv);
			}
//			System.err.println("LOG SCALE "+logScale);
		}
		
		

	
		
	}
	
	Likelihood[] priors;
//	RealVector[] prob; // stores P(D | G)
	
	File sfs_in, sfs_out;//  this stores the trained probs, or has input probs;
	
	final double[] v_het, v; //this keeps the current het frequency;
	public void print(String suffix) {
		try{
			PrintWriter pw = new PrintWriter(new FileWriter(new File(sfs_out,suffix)));
			for(int k=0; k<pop.length; k++){
				pw.print("p"+k+"\thobs"+k+"\t");
				
			}
			pw.println("prob");
			it.reset();
			while(true){
				int[] inds = it.inds;
				for(int k=0; k<inds.length; k++){
					int j = inds[k];
					double[] re = hweprobs[k].getColumn(j);
					double tot = re[0]+re[1]+re[2];
					double p = (re[1] + 2*re[0]) / (2*tot);
					double hobs  = re[1]/tot;
				//	double hexp = 2.0*p*(1-p);
					//double diff = hexp-hobs;
					//double f = diff==0 ? 0 : diff/hexp;
					if(false){pw.print(p+"\t");
					pw.print(hobs+"\t");
					}else{
						pw.print(String.format("%6.4g",p).trim()+"\t");
						pw.print(String.format("%6.4g",hobs).trim()+"\t");
					}
				}
				pw.println(String.format("%5.3g", this.AFS.get(inds).prob));
				if(!it.incr()) break;
			}
			pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	
	
	public double expectation(){
		it.reset();
		double sum=0;
//		double[] prob = new double[this.pop.length];
		while(true){
			int[] inds = it.inds;
			double lh = 1.0;
			for(int k=0; k<inds.length; k++){
				lh *=this.priors[k].lhood.getEntry(inds[k]);
			}
			sum+=this.AFS.getBF(inds,lh);
			if(!it.incr()) break;
		}
		it.reset();
		while(true){
			int[] inds = it.inds;
			this.AFS.addCounts(inds, sum);
			if(!it.incr()) break;
		}
		this.scale = sum;
		double sum1 = Math.log(sum);
		for(int k=0; k<priors.length; k++){
				sum1 = sum1  + this.priors[k].logScale;
		}
	
		return sum1;
	}
	
	
	
	public void train( int numit){
		double logPrev = Double.NEGATIVE_INFINITY;
		for(int i=0; i<numit; i++){
			double lh1 = this.expectation();
			System.err.println(lh1+" "+logPrev);
			if(lh1< logPrev) System.err.println(" warning ");
			logPrev = lh1;
			transfer();
		}
	}
	
	public void transfer(){
		this.AFS.transfer();
		
		if(CHECK) this.AFS.validate();
	}
	
	double scale =1.0; //
	
	public void setData(RealMatrix[]lh, boolean[][]incl, int[] inds){
		for(int k1=0; k1<inds.length; k1++){
			int k = inds[k1];
			int[] rowsToInclude = getIncl(incl[k]);
			if(rowsToInclude!=null) priors[k1].setData(lh[k], this.hweprobs[k1],rowsToInclude);
		}
	}
	
	double[] ci = new double[] {};//0.05,0.5,0.95};
	
	
	public RealMatrix calculate(){
		it.reset();
		double[] lh = it.lh;
		while(true){
			int[] inds = it.inds;
			AFS.calcLHWeightedFuncts(inds, lh, this.pmap);
			if(!it.incr()) break;
		}
		RealMatrix res = new Array2DRowRealMatrix(ci.length+1,pmap.length);
		
		double[] v = new double[ci.length];
		res.setRow(0, it.lh);
		if(ci.length>0){
		for(int j=0; j<pmap.length; j++){
			double[] quants = new double[ci.length];
			pmap[j].getQuantiles(ci, quants, v);
			res.setRow(j+1, quants);;
		}
		}
		//RealVector re =   new ArrayRealVector(it.lh);
		if(scale!=1){
			for(int j=0; j<res.getColumnDimension(); j++){
				res.setColumnVector(j, res.getColumnVector(j).mapDivideToSelf(scale));
			}
		}
		//re.mapDivideToSelf(scale);
		return res;
	}
	
final public String name;

	private int[] getIncl(boolean[]incl) {
	int cnt=0;
	for(int k=0; k<incl.length; k++){
		if(incl[k])cnt++;
	}
	if(cnt==0){
		return null;
	}
	int[] rowsToInclude = new int[cnt];
	cnt=0;
	for(int k=0; k<incl.length; k++){
		if(incl[k]){
			rowsToInclude[cnt]=k;
			cnt++;
		}
	}
	
	return rowsToInclude;
	}

	

	

	private void normalise(RealMatrix lh1) {
		int rd = lh1.getRowDimension();
		//int cd = lh1.getColumnDimension();
		for(int i=0; i<rd; i++){
			RealVector rv = lh1.getRowVector(i);
			
			rv.mapDivideToSelf(rv.getL1Norm());
			lh1.setRowVector(i, rv);
		}
		
	}
	
	public static Object subIndex(Object array, int[] incl){
		if(incl==null) return array;
		Class clazz = Array.get(array, 0).getClass();
		Object res = Array.newInstance(clazz, incl.length);
		for(int k=0; k<incl.length; k++){
			Array.set(res, k, Array.get(array, incl[k]));
		}
		return res;
	}
	public static Object subIndexI(Object array, int[] incl){
		if(incl==null) return array;
		Object res = Array.newInstance(int.class, incl.length);
		for(int k=0; k<incl.length; k++){
			Array.set(res, k, Array.get(array, incl[k]));
		}
		return res;
	}
	
	public  static String getName(int[] pops) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<pops.length; i++){
			sb.append("_"+pops[i]);
		}
		return sb.toString();
	}

	// int[] inds_incl; /// which indices to include if multiple presented
	public SFS(int[] pops, int ploidy, int noFuncts, String name) {
		this.name = name;
		this.v_het = new double[pops.length];
		this.noFuncts = noFuncts;
		this.v = new double[pops.length];
		this.pop = new int[pops.length];
		this.lim = new int[pops.length];
		hweprobs = new RealMatrix[pops.length];
		this.priors = new Likelihood[v.length];
	//	this.prob = new double[v.length];
		it = new ArrayIt(lim, noFuncts);
		this.sfs_out = new File(name+"_sfs_out");
		sfs_out.mkdir();
		this.pmap = new ProbMap[noFuncts];
		for(int j=0; j<pmap.length; j++){
			pmap[j] = new ProbMap();
		}
	}
final int noFuncts;

	public int getDimension() {
		return this.noFuncts;
	}
	

}