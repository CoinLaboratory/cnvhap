package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import pal.math.UnivariateFunction;
import cern.colt.matrix.DoubleMatrix2D;

public class TrainableNegBinomial implements ProbabilityDistribution, UnivariateFunction, Serializable{
public static void main(String[] args){
	try{
		TrainableNegBinomial tnb = new TrainableNegBinomial(10,0.5);
		double  b = 9;
		for(int i=0; i<10; i++){
			System.err.println(i+" "+tnb.probability((double)i/10.0 + b));
		}
	}catch(Exception exc){
		
	}
}
public double  probability(double x, int mixComponent){
	return this.probability(x);
}
public void setMinMax(double min, double max){
	//throw new RuntimeException("!!");
}
public void addCount(Double b, double val,
		SimpleExtendedDistribution1 mixe1, lc1.stats.ProbabilityDistribution disty){
	disty.addCount(b, val);
	mixe1.addCount(0, val);
	//throw new RuntimeException("!!");
}
List<Double> obsx = new ArrayList<Double>();
List<Double> obsv = new ArrayList<Double>();
double sum=0;


public lc1.stats.ProbabilityDistribution clone(double u,
		SimpleExtendedDistribution1 dist1){
	return this.clone(u);
}
	//final NegativeBinomial nb;
	double r;
	double p;
//	List<Double> obsx  = new ArrayList<Double>();
//	List<Double> obsv = new ArrayList<Double>();
	// double sum=0;
	// double sumln=0;
	 double N =0;
	public void initialise(){
		obsx.clear();
		obsv.clear();
		sum=0;
		//sumln=0;
		N=0;
	}
/*	public TrainableGammaDistribution(double shape, double scale){
		//gamma = new Gamma(shape, scale, null);
		this.k = shape;
		this.theta = scale;
		
	}*/
	
	public TrainableNegBinomial(double shape, double scale){
		//gamma = new Gamma(shape, scale, null);
		this.r = shape;
		this.p = scale;
	//	this.r = d*((p)/(1-p))+1;
	//	nb = new NegativeBinomial(r,p, null);
		//if(shape>0){
			//theta = mode/(k-1);
		//
		//}
			//double mde = this.inverse(0.5);
			//double cdf = this.cdf(mde);
			//double mod = (k-1)*theta;
		//Logger.global.info(k+" "+theta);
	//	
	}
	
	
	public TrainableNegBinomial(String string, double mode_i_r, double p_i_r,
			double round, double d) {
		this.p = p_i_r;
		this.r = d*((p)/(1-p))+1;
		this.meanPrior = location();
		this.stddevPrior = variance();
		
		//this(mean_i_r+1, var_i_r);
	}
	public double location(){
		return r*(p/(1-p));
	}
	public double variance(){
		return r*(p/Math.pow((1-p),2));
	}
	

	private double cdf(double mde) {
		throw new RuntimeException("!!");
	}
	
	public double inverse(double y){
		 throw new RuntimeException("!!");
	}
	public double probability(double x){
	//	cern.jet.random.NegativeBinomial nb = new NegativeBinomial(r,p,null);
	//	return nb.pdf((int) Math.round(k));
//		double bin = Arithmetic.binomial(k+r-1, r-1);
//		return bin *;
		return Math.exp(GammaFunction.lnGamma(r+x) - GammaFunction.lnGamma(r) -GammaFunction.lnGamma(x+1) )*Math.pow(p,r)
		*Math.pow(1-p,x) ;
		//return GammaDistribution.pdf(x, k, theta);
	}
 public synchronized void addCount(double d, double w) {
	this.obsx.add(d);
	 this.obsv.add(w);
	 sum+=d*w;
//	 sumln+= Math.log(d*w);
	 N+=w;
 }
 public void maximise(){
	 /*double x = this.initv();
	 UnivariateMinimum uvm = new UnivariateMinimum();
	 this.k = uvm.findMinimum(x, this, 5);
	 this.theta = theta();
	 */
	 
 }
 public String toString(){
	 return this.r+":"+this.p;
 }

@Override
public double evaluate(double argument) {
//	double N = obsx.size();
	//double res = (k-1)*sumln - N * k - N*k*Math.log(sum/(k*N))
	// - N * GammaFunction.lnGamma(k);
	//return -1 * res;
	throw new RuntimeException("!!");
}


@Override
public double getLowerBound() {
	return 0;
}
@Override
public double getUpperBound() {
	return 1e3;
}
public void print(PrintWriter pw) {
	pw.println("--gammaRate  "+String.format("%5.3g", this.r).trim()+":"
			+String.format("%5.3g",this.p).trim());
	
	
}

@Override
public void addCounts(ProbabilityDistribution probabilityDistribution) {
	// TODO Auto-generated method stub
	
}

@Override
public ProbabilityDistribution clone(double u) {
	return new TrainableNegBinomial(this.r, this.p);
}
public ProbabilityDistribution clone(){
	return new TrainableNegBinomial(this.r, this.p);
}




@Override
public int fillVariance( DoubleMatrix2D y, int numObs, double pseudo) {
	int i=0;
double avg = this.location();// this.average(0.0);
	 for( i=0; i<this.obsx.size(); i++){
		 int k = numObs+i;
		 double v = this.obsv.get(i);
		
		 y.setQuick(k,0, v*( Constants.transformVariance(this.obsx.get(i) -avg )));
	 }
	 if(Math.abs(pseudo)>1e-5){
		 int k = numObs+i;
		 double v = Math.abs(pseudo);
		double vv = Constants.transformVariance(pseudo<0 ? this.stddevPrior : this.variance());
			if(Constants.CHECK && (Double.isNaN(vv) || Double.isNaN(vv))){
	   			throw new RuntimeException("!!");
	   		}
		 y.setQuick(k,0, v*( vv));
		 i++;
	 }
	 return i;
}

@Override
public double[] getCount(double[] angle) {
	// TODO Auto-generated method stub
	return null;
}


@Override
public void getInterval(double[] in, double[] res_inR) {
	// TODO Auto-generated method stub
	
}

@Override
public double getMean() {
	return this.r*((1-p)/(p));
}

@Override
public int getParamIndex() {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public double getParamValue(int n1) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public String id() {
	// TODO Auto-generated method stub
	return null;
}

@Override
public void maximise(double d, double e, double f) {
	// TODO Auto-generated method stub
	
}

@Override
public String name() {
	// TODO Auto-generated method stub
	return null;
}

@Override
public int numObs() {
	return this.obsx.size();
}

@Override
public double plotObservations(String string, boolean b, XYSeries obs,
		boolean swtch) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public void plotTheoretical(String string, boolean b, XYSeries theor) {
	// TODO Auto-generated method stub
	
}

@Override
public double prior() {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public void recalcName() {
	// TODO Auto-generated method stub
	
}

@Override
public double sample() {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public double scale() {
	return Math.sqrt(r *(1-p))/p;
}

/*@Override
public void setParam(int type, double rho) {
	// TODO Auto-generated method stub
	
}*/
double meanPrior,stddevPrior;
public void setParam(int type,  double d){
	//if(true) throw new RuntimeException("!!");
	if(type==0){
		this.r = d;
		this.meanPrior = d;
//		this.location = d;
	}
	else if(type==1){
		double mean = Math.max(r,0.001);
		double var = Math.max(d,0.001);
		this.stddevPrior = var;
		p= Math.min(0.99, Math.max(0.001,mean/var +1));  
		this.r = mean*((p-1)/(p));
//		double d1 =Math.sqrt(Math.max(1e-10,d));
	//	if(Constants.CHECK && (Double.isNaN(d1) || d1==0)) {
		//	throw new RuntimeException("!!");
	//	}
///		this.stddevPrior = d1;
	//	this.scale = d1;
	}
	else{
		throw new RuntimeException("!!");
//		this.skewPrior[0] = d;
//		this.shape = d;
	}
//	this.recalcName();
	}

@Override
public void setParamValue(int n1, double val) {
	// TODO Auto-generated method stub
	
}

@Override
public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
	// TODO Auto-generated method stub
	
}

@Override
public void setPriors(ProbabilityDistribution distx, int type) {
	// TODO Auto-generated method stub
	
}

@Override
public double sum() {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public void transfer(double pseudoC) {
	// TODO Auto-generated method stub
	
}

@Override
public void transfercounts(EmissionState innerState, int phen_index, int i) {
	// TODO Auto-generated method stub
	
}

@Override
public void updateParamIndex() {
	// TODO Auto-generated method stub
	
}

@Override
public void variance(double[] sum) {
	// TODO Auto-generated method stub
	
}

@Override
public int compareTo(Object o) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public double evaluate(double[] argument) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public double getLowerBound(int n) {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public int getNumArguments() {
	// TODO Auto-generated method stub
	return 0;
}

@Override
public OrthogonalHints getOrthogonalHints() {
	// TODO Auto-generated method stub
	return null;
}

@Override
public double getUpperBound(int n) {
	// TODO Auto-generated method stub
	return 0;
}
public void setCoverage(double d) {
	//(r-1)*p/(1-p) = d
	this.r = d*((p)/(1-p))+1;
	
}
@Override
public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[]  noCop,
		 double pseudo) {
		    	int i=0;
		    	int nocols = noCop.length;
		    	 for( i=0; i<this.obsx.size(); i++){
		    		 int k = numObs+i;
		    		 double v = this.obsv.get(i);
		    		 for(int kk=0;kk<nocols; kk++){
		    			 x.setQuick(k, kk, v * (double)noCop[kk]);
		    		 }
		    		
		    		
		    		
		    		 y.setQuick(k,0, v*( this.obsx.get(i) ));
		    	 }
		    	if(Math.abs(pseudo)>1e-5){
		    		double v = Math.abs(pseudo);
		    		int k = numObs+i;
			   		 for(int kk=0;kk<nocols; kk++){
			   			 x.setQuick(k, kk, v * (double)noCop[kk]);
			   		 }
			   		
			   		
			   		
			   		 y.setQuick(k,0, v*(pseudo < 0  ? this.meanPrior : this.location() ));
			   		i++;
		    	}
		    	 return i;
			}
}
