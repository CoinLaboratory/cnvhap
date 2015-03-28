package lc1.stats;

import java.io.PrintWriter;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import JSci.maths.statistics.ProbabilityDistribution;
import cern.colt.matrix.DoubleMatrix2D;


/** warning - we only maximise over the first distribution */
public class Mixture  extends ProbabilityDistribution implements lc1.stats.ProbabilityDistribution, MultivariateFunction, Comparable{
	public final lc1.stats.ProbabilityDistribution[] dist;
	final double[] distribution;
	//int[] paramToDist;
	// int[] paramToPos;
	//int totParams;
	String id;
	public void variance( double[] sum){
	  	dist[0].variance( sum);
	        
	  }
	
	public void setCoverage(double d){
    	for(int k=0; k<dist.length; k++){
    		dist[k].setCoverage(d);
    	}
    }
	public void setPriors(lc1.stats.ProbabilityDistribution distx, int type){
		this.dist[0].setPriors(((Mixture)distx).dist[0], type);
	}

	public double  probability(double x, int mixComponent){
		if(true) throw new RuntimeException("!!");
		return this.dist[mixComponent].probability(x)*this.mix.probs[mixComponent];
	}
	public void setParam(int i, double e){
		dist[0].setParam(i,e);
	}
	
	public int fillVariance(DoubleMatrix2D y,  int numObs, double pseudo){
		return dist[0].fillVariance(y, numObs, pseudo);
	}
	
	@Override
	public void addCount(Double r2,double val,
			SimpleExtendedDistribution1 mixe1,
			lc1.stats.ProbabilityDistribution probDistG) {
				for(int i1=0; i1<distribution.length; i1++){
					double v = val*distribution[i1];
						mixe1.addCount(i1, v);
						probDistG.addCount(r2,  v);
					
				}
		
	}
	public Mixture(lc1.stats.ProbabilityDistribution[] dist, SimpleExtendedDistribution1 mix){
		this.dist = dist;
		this.distribution = new double[dist.length];
		this.mix = mix;
		make();
	}
	public String toString(){
		return dist[0].toString();
	}
	public double scale(){
		return dist[0].scale();
	}
	 @Override
		public int numObs() {
			return dist[0].numObs();
		}
	public void make(){
		//totParams =0;
		StringBuffer id = new StringBuffer();
		for(int i=0; i<dist.length; i++){
		//	totParams+=dist[i].getNumArguments();
			id.append(dist[i].id());
		}
	//	this.paramToDist = new int[totParams];
		//this.paramToPos = new int[totParams];
		/*totParams =0;
		for(int i=0; i<dist.length; i++){
			for(int k=0; k<dist[i].getNumArguments(); k++){
			//	paramToDist[totParams+k] = i;
			//	paramToPos[totParams+k] = k;
			}
			totParams+=dist[i].getNumArguments();
		}*/
		this.id = id.toString();
	}
	public Mixture(lc1.stats.ProbabilityDistribution dist, double min, double max, 
			SimpleExtendedDistribution1 mix
			){
		this(new lc1.stats.ProbabilityDistribution[] {dist, new UniformDistribution(min, max)}, 
				mix);
				
				//new double[] {1 - mixCoeff, mixCoeff});
		
	}
	public Mixture(Mixture mixt, double u,boolean keepMix){
		
		this.dist = new lc1.stats.ProbabilityDistribution[mixt.dist.length];
		this.distribution = new double[dist.length];
		for(int i=0; i<dist.length; i++){
			dist[i] = mixt.dist[i].clone(u);
		}
		this.mix = keepMix ? mixt.mix : new SimpleExtendedDistribution1(mixt.mix.probs, u);
		make();
	}
public Mixture(Mixture mixt, double u,SimpleExtendedDistribution1 mix1){
		
		this.dist = new lc1.stats.ProbabilityDistribution[mixt.dist.length];
		this.distribution = new double[dist.length];
		for(int i=0; i<dist.length; i++){
			dist[i] = mixt.dist[i].clone(u);
		}
		this.mix = mix1;
		make();
	}
	SimpleExtendedDistribution1 mix;
	 @Override
	    public double cumulative(double arg0) {
	       throw new RuntimeException("!!");
	    }

	    @Override
	    public double inverse(double arg0) {
	       throw new RuntimeException("!!");
	    }
	    
	    public void calcDist(double arg0){
	    	for(int i=0; i<dist.length; i++){
				distribution[i] = dist[i].probability(arg0)*mix.probs[i];
			}
	    }
	@Override
	public double probability(double arg0) {
		this.calcDist(arg0);
		return Constants.sum(distribution);
	
	}
	public void addCount(double d, double w) {
		//if(Constants.CHECK && Double.isNaN(d)){
			//throw new RuntimeException("!!");
	//	}
		{
		this.calcDist(d);
		Constants.normalise(distribution);
		for(int i=0; i<distribution.length; i++){
			double frac = distribution[i];
			mix.addCount(i, w*frac);
			/*if(i==1 && frac > 0.5){
				System.err.println("h");
			}*/
			dist[i].addCount(d, w*frac);
		}
		}
		
	}
	@Override
	public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution) {
		throw new RuntimeException("!!");
		
	}
	public lc1.stats.ProbabilityDistribution clone(double u) {
		return new Mixture(this, u,false);
	}
	public lc1.stats.ProbabilityDistribution clone() {
		return new Mixture(this, Double.POSITIVE_INFINITY,false);
	}
	
	public lc1.stats.ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1){
		return new Mixture(this, u, dist1);
	}
	public double[] getCount(double[] angle) {
		throw new RuntimeException("!!");
	}
	public double getMean() {
		return this.dist[0].getMean();
	}
	public int getParamIndex() {
		int sum = this.mix.getParamIndex();
		for(int i=0; i<dist.length; i++){
		sum +=this.dist[0].getParamIndex();
		}
		return sum;
	}
	
	
	public double getParamValue(int n1) {
		return this.dist[0].getParamValue(n1);
	//	return this.dist[paramToDist[n1]].getParamValue(paramToPos[n1]);
	}
	public String id() {
		return id;
	}
	public void initialise() {
		for(int i=0; i<dist.length; i++){
			dist[i].initialise();
		}
		this.mix.initialise();
		
	}
	public void maximise(double d, double e, double f) {
		dist[0].maximise(d, e, f);
		
	}
	public String name() {
		return id;
	}
	public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch) {
		return dist[0].plotObservations(string, b, obs, swtch);
	}
	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		dist[0].plotTheoretical(string, b, theor);
		
	}
	public double prior() {
		return dist[0].prior();
	}
	public void recalcName() {
		dist[0].recalcName();
		
	}
	public double sample() {
		throw new RuntimeException("!!");
	}
	public void setParamValue(int n1, double val) {
		this.dist[0].setParamValue(n1, val);
		
	}
	public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution[] tmp) {
		dist[0].setParamsAsAverageOf(tmp);
		
	}
	public double sum() {
		return dist[0].sum();
	}
	public void transfer(double pseudoC) {
		for(int i=0; i<dist.length; i++){
			dist[i].transfer(pseudoC);
		}
		this.mix.transfer(pseudoC);
		
	}
	public void transfercounts(EmissionState innerState, int phen_index, int i) {
		throw new RuntimeException("!!");
		
	}
	public void updateParamIndex() {
		for(int i=0; i<dist.length; i++){
			dist[i].updateParamIndex();
		}
		
	}
	public int compareTo(Object o) {
		return dist[0].compareTo(((Mixture)o).dist[0]);
	}
	public double evaluate(double[] argument) {
		return this.dist[0].evaluate(argument);
	}
	public double getLowerBound(int n) {
		return this.dist[0].getLowerBound(n);
	}
	public int getNumArguments() {
		return this.dist[0].getNumArguments();
	}
	public OrthogonalHints getOrthogonalHints() {
		return this.dist[0].getOrthogonalHints();
	}
	public double getUpperBound(int n) {
		return this.dist[0].getUpperBound(n);
	}
	@Override
	public void getInterval(double[] in, double[] res_inR) {
		 dist[0].getInterval(in, res_inR);
		
	}
	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[]  noCop,
			 double pseudo) {
		return dist[0].fill(x, y, numObs, noCop, pseudo);
	}
	@Override
	public void print(PrintWriter pw) {
		this.dist[0].print(pw);
		
	}
	
	public void setMinMax(double min, double max){
		this.dist[1].setMinMax(min, max);
	}
}
