package lc1.stats;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;

public class PseudoMixture extends PseudoDistribution {

	public final lc1.stats.PseudoDistribution[] dist;
	final double[] distribution;
	SimpleExtendedDistribution1 mix;
	public PseudoMixture(PseudoDistribution[] pseudoDistributions,
			SimpleExtendedDistribution1 dist2) {
		this.dist = pseudoDistributions;
		this.mix = dist2;
		this.distribution = new double[dist.length];;
	}
	public double weight(){
		return this.mix.probs[0];
	}
	public EmissionStateSpace getEmissionStateSpace(){
		return this.dist[0].getEmissionStateSpace();
	}
	@Override
	public void addCount(int d, double w) {
		if(true) throw new RuntimeException("!!");
		
	}
	
	@Override
	public double[] calcDistribution(double[] distribution,
			EmissionStateSpace emStSp, int pos) {
		// TODO Auto-generated method stub
		if(true) throw new RuntimeException("!!");
		return null;
	}
	@Override
	public PseudoDistribution clone() {
		// TODO Auto-generated method stub
		if(true) throw new RuntimeException("!!");
		return null;
	}
	@Override
	public PseudoDistribution clone(double swtch) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public double[] counts() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public Integer fixedInteger() {
	
		return null;
	}
	@Override
	public int getMax() {
		if(true) throw new RuntimeException("!!");		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public String getPrintString() {
		// TODO Auto-generated method stub
		return this.dist[0].getPrintString();
	}
	@Override
	public void initialise() {
		for(int k=0; k<dist.length; k++){
			dist[k].initialise();
		}
		// TODO Auto-generated method stub
		
	}
	@Override
	public boolean isMeasured() {
		// TODO Auto-generated method stub
		return false;
	}
	@Override
	public double logProb() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public double[] probs() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return null;
	}
	@Override
protected Double getBEst(EmissionStateSpace emstsp) {
		return this.dist[0].getBEst(emstsp);
	}
	@Override
	public double probs(int obj_i) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public double scoreB(int j, int i) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	
	 
	public void calcDist(EmissionStateSpace emStSp, int j, int i){
    	for(int k=0; k<dist.length; k++){
			distribution[k] = dist[k].scoreBR(emStSp, j, i)*mix.probs[k];
		}
    }
	
	@Override
	public void addRBCount(EmissionStateSpace emStSp, int j, double w, int i) {
		this.calcDist(emStSp, j, i);
		Constants.normalise(distribution);
		for(int k=0; k<distribution.length; k++){
			double frac = distribution[k];
			mix.addCount(k, w*frac);
			
			dist[k].addRBCount(emStSp,j, w*frac,i);
		}
	}
	public double scoreBR(EmissionStateSpace emstsp, int j, int i, int mixComponent) {
		this.calcDist(emstsp, j,i);
		//if(distribution[0]/Constants.sum(distribution)<0.1){
		//	System.err.println("h");
		//}
		return this.distribution[mixComponent];
//		return this.dist[mixComponent].scoreBR(emstsp,j, i);//*this.mix.probs[mixComponent];
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emStSp, int j, int i) {
		this.calcDist(emStSp,j,i);
		double sum = Constants.sum(distribution);
		
		return sum;
		
	}
	@Override
	public void setCounts(int i1, double cnt) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		
	}
	@Override
	public void setFixedIndex(int k) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		
	}
	@Override
	public void setProb(double[] prob) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		
	}
	@Override
	public void setProbs(int to, double d) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		
	}
	@Override
	public double sum() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public PseudoDistribution swtchAlleles() {
		// TODO Auto-generated method stub
		for(int k=0; k<this.dist.length; k++){
			dist[k].swtchAlleles();
		}
		return this;
	}
	@Override
	public double totalCount() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public void transfer(double pseudoC1) {
		// TODO Auto-generated method stub
		for(int k=0; k<this.dist.length; k++){
			this.dist[k].transfer(pseudoC1);
		}
	}
	@Override
	public void validate() {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void addCounts(ProbabilityDistribution probabilityDistribution) {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		
	}
	@Override
	public double sample() {
		if(true) throw new RuntimeException("!!");
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
		// TODO Auto-generated method stub
		if(true) throw new RuntimeException("!!");		
	}
	
}
