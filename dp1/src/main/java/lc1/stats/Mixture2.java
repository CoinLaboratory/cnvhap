package lc1.stats;

import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;
/** warning - we only maximise over the first distribution */
public class Mixture2   implements lc1.stats.ProbabilityDistribution2,  Comparable{
	public final lc1.stats.ProbabilityDistribution2[] dist;
	final public double[] distribution;
	//int[] paramToDist;
	// int[] paramToPos;
	//int totParams;
	String id;
	public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
	){
	 this.dist[0].setPriors(((Mixture2)pr2).dist[0], type, x);	
	}
	public Mixture2(lc1.stats.ProbabilityDistribution2[] dist, SimpleExtendedDistribution1 mix){
		this.dist = dist;
		this.distribution = new double[dist.length];
		this.mix = mix;
		make();
	}
	 public void variance( int type, double[] sum){
    	dist[0].variance( type, sum);
          
    }
	 
	 public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
				DoubleMatrix2D covar, int numObs,  double pseudo){
		return  dist[0].fillVariance(y, yb, covar, numObs, pseudo);
	 }
	
	 public void setParam(int type,int i,double rho){
		 dist[0].setParam(type, i,rho);
	 }
	public String toString(){
		return dist[0].toString();
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
	public Mixture2(lc1.stats.ProbabilityDistribution2 dist, double minx, double maxx, 
			double miny, double maxy,
			SimpleExtendedDistribution1 mix
			){
		this(new lc1.stats.ProbabilityDistribution2[] {dist, new UniformDistribution2(minx, maxx, miny, maxy)}, 
				mix);
				
				//new double[] {1 - mixCoeff, mixCoeff});
		
	}
	public Mixture2(Mixture2 mixt, double u,boolean keepMix){
		
		this.dist = new lc1.stats.ProbabilityDistribution2[mixt.dist.length];
		this.distribution = new double[dist.length];
		for(int i=0; i<dist.length; i++){
			dist[i] = mixt.dist[i].clone(u);
		}
		this.mix = keepMix ? mixt.mix : new SimpleExtendedDistribution1(mixt.mix.probs, u);
		make();
	}
public Mixture2(Mixture2 mixt, double u,SimpleExtendedDistribution1 mix1){
		
		this.dist = new lc1.stats.ProbabilityDistribution2[mixt.dist.length];
		this.distribution = new double[dist.length];
		for(int i=0; i<dist.length; i++){
			dist[i] = mixt.dist[i].clone(u);
		}
		this.mix = mix1;
		make();
	}
	public SimpleExtendedDistribution1 mix;
	
	    public void calcDist(double arg0, double argy){
	    	for(int i=0; i<dist.length; i++){
				distribution[i] = dist[i].probability(arg0, argy)*mix.probs[i];
			}
	    }
	@Override
	public double probability(double arg0, double argy) {
		this.calcDist(arg0, argy);
		double res =  Constants.sum(distribution);
		return res;
	}
	
	
	public double probability(double arg0, double argy, int mixComp) {
		this.calcDist(arg0, argy);
		double res =  distribution[mixComp];
		return res;
	}
	public void addCount(double dx, double dy, double w) {
		if(w>Constants.countThresh()){
		this.calcDist(dx, dy);
		Constants.normalise(distribution);
		for(int i=0; i<distribution.length; i++){
			mix.addCount(i, w*distribution[i]);
			dist[i].addCount(dx,dy, w*distribution[i]);
		}
		}
		
	}
	public void addCounts(lc1.stats.ProbabilityDistribution2 probabilityDistribution) {
		throw new RuntimeException("!!");
		
	}
	public lc1.stats.ProbabilityDistribution2 clone(double u) {
		return new Mixture2(this, u,false);
	}
	public lc1.stats.ProbabilityDistribution2 clone() {
		return new Mixture2(this, Double.POSITIVE_INFINITY,false);
	}
	public double[] getCount(double[] angle) {
		throw new RuntimeException("!!");
	}
	/*public double[] getMean() {
		return this.dist[0].getMean();
	}*/
	public int getParamIndex() {
		int sum = this.mix.getParamIndex();
		for(int i=0; i<dist.length; i++){
		sum +=this.dist[0].getParamIndex();
		}
		return sum;
	}
	
	
	/*public double getParamValue(int n1) {
		return this.dist[0].getParamValue(n1);
	//	return this.dist[paramToDist[n1]].getParamValue(paramToPos[n1]);
	}*/
	public String id() {
		return id;
	}
	public void initialise() {
		for(int i=0; i<dist.length; i++){
			dist[i].initialise();
		}
		this.mix.initialise();
		
	}
	public void maximise(double d, double e, double f, double d1, double e1, double f1, double g) {
		dist[0].maximise(d, e, f, d1, e1, f1, g);
		
	}
	public String name() {
		return id;
	}
	/*public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch) {
		return dist[0].plotObservations(string, b, obs, swtch);
	}
	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		dist[0].plotTheoretical(string, b, theor);
		
	}*/
	/*public double prior() {
		return dist[0].prior();
	}*/
	public void recalcName() {
		dist[0].recalcName();
		
	}
	public double sample() {
		throw new RuntimeException("!!");
	}
	/*public void setParamValue(int n1, double val) {
		this.dist[0].setParamValue(n1, val);
		
	}
	public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution2[] tmp) {
		dist[0].setParamsAsAverageOf(tmp);
		
	}*/
/*	public double sum() {
		return dist[0].sum();
	}*/
	public void transfer(double pseudoC) {
		for(int i=0; i<dist.length; i++){
			dist[i].transfer(pseudoC);
		}
		//this.mix.transfer(pseudoC);
		
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
		return dist[0].compareTo(((Mixture2)o).dist[0]);
	}
	/*public double evaluate(double[] argument) {
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
	public void getIntervalX(double[] in, double[] res_inR) {
		 dist[0].getIntervalX(in, res_inR);
		
	}
	@Override
	public void getIntervalY(double[] in, double[] res_inR) {
		 dist[0].getIntervalY(in, res_inR);
		
	}*/
	@Override
	public void getInterval(double[] input, DoubleMatrix2D res, double[] mean) {
		this.dist[0].getInterval(input, res, mean);
		
	}
	@Override
	public int numObs() {
		return this.dist[0].numObs();
	}
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y,DoubleMatrix2D yB, int numObs,double[] noCop
			,  double pseudo){
		return dist[0].fill(x, y,yB, numObs,noCop, pseudo);
	}
	@Override
	public void setToExclude() {
	  Arrays.fill(this.mix.probs,0.0);
	  Arrays.fill(this.mix.counts,0.0);
	  Arrays.fill(this.mix.pseudo,0.0);
	mix.probs[1] = 1.0;
	mix.counts[1] = 1.0;
	mix.pseudo[1] = 1.0;
	}
	public void print(PrintWriter pw) {
	this.dist[0].print(pw);
		
	}
	@Override
	public ProbabilityDistribution2 clone(double u,
			SimpleExtendedDistribution1 s1) {
	return new Mixture2(this,u, s1);
	}
	@Override
	public void setMinMax(double minR, double maxR, double min, double max) {
		 ((UniformDistribution2)dist[1]).setMinMax(minR,maxR, min, max);
		
	}
	@Override
	public void addCount(Double r2, Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG) {
	//	double[] distribution = ((Mixture2)dist).distribution;
		//		Mixture2 distG =  this.probeOnly[i] ? : this.getDistributionGlobal(noCop, noB);
				for(int i1=0; i1<distribution.length; i1++){
					
					double v = val*distribution[i1];
					
						mixe1.addCount(i1, v);
						probDistG.addCount(r2, b, v);
					
				}
		
	}
	
	
	
}
