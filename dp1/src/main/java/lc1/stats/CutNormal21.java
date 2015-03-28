package lc1.stats;

import java.io.PrintWriter;

import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;

public class CutNormal21  implements ProbabilityDistribution2, Comparable{
	
	//double lessThanZeroProb; //in y direction
	 //double greaterThanOneProb;
	final TrainableNormal2 dist;
	 double mlscale;
	
	final double min;
	final double max;
	public String toString(){
		return dist.toString();
	}
	public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
	){
	 this.dist.setPriors(((CutNormal21)pr2).dist, type, x);	
	}
	 public CutNormal21(TrainableNormal2 dist, double min, double max, boolean extrapolate) {
		this.dist = dist;
		this.min = min;
		this.max = max;
		this.topbottom = dist.meany==min || dist.meany==max;
	//	this.update();
		this.extrapolate = extrapolate;
		this.mlscale = dist.sigma_y/5.0;
	}
	 final boolean extrapolate;
	
	public CutNormal21(CutNormal21 cutNormal, double u) {
	//	this.lessThanZeroProb = cutNormal.lessThanZeroProb;
	//	this.greaterThanOneProb = cutNormal.greaterThanOneProb;
		this.dist = cutNormal.dist.clone(u);
		this.min = cutNormal.min;
		this.max = cutNormal.max;
		this.topbottom = cutNormal.topbottom;
		this.extrapolate = cutNormal.extrapolate;
		this.mlscale = cutNormal.mlscale;
	}

	public CutNormal21(CutNormal21 cutNormal) {
		//this.lessThanZeroProb = cutNormal.lessThanZeroProb;
		//this.greaterThanOneProb = cutNormal.greaterThanOneProb;
		this.dist = cutNormal.dist.clone();
		this.min = cutNormal.min;
		this.max = cutNormal.max;
		this.topbottom = cutNormal.topbottom;
		this.extrapolate = cutNormal.extrapolate;
		this.mlscale = cutNormal.mlscale;
	}


	public CutNormal21 clone(){
	        return new CutNormal21(this);
	    }
	
	public CutNormal21 clone(double u){
        return new CutNormal21(this, u);
    }
//boolean upper  = false; // 

	
	
	@Override
	public void addCount(double x, double y, double value) {
		
		
	
		
		
			if(y==this.min){
				double quantileL = probLessThan(x,min, true);
				if(quantileL>1e-10){
					double v = sampleLessThan(cdf, muy1, sigmay1);
					this.dist.addCount(x,v, value);
				}
				else{
					this.dist.addCount(x,y, value);
				}
			}
			else if(y==this.max){
			
				double quantileG = probLessThan(x,max,false);
				if(quantileG>1e-10){
					double v =  this.sampleGreaterThan(1-cdf, muy1, sigmay1);
					dist.addCount(x,v, value);
				}
				else{
					this.dist.addCount(x,y, value);
				}
			}
			else{
				this.dist.addCount(x,y, value);
			}
		
		
		
		/*double sigma_x = dist.sigma_x;
		double sigma_y = dist.sigma_y;
		if(extrapolate && (y==this.min || y==max)){
				double muy1 = 
					dist.meany+dist.rho()*(sigma_y/sigma_x)*(x - sigma_x);
				double sigmay1 = Math.sqrt(1 - Math.pow(dist.rho(), 2))*sigma_y;
				double cdf = normal.cdf(y, muy1,sigmay1);
				if(y==min){
				     dist.addCount(x,sampleLessThan(cdf, muy1, sigmay1), value);
				}
				else{
					dist.addCount(x, this.sampleGreaterThan(1-cdf, muy1, sigmay1), value);
				}
		}
		else{
			this.dist.addCount(x,y, value);
		}*/
		
	}
	
	public double sampleLessThan(double quantile, double meany, double sigmay) {
		if(quantile==0) return min;
		double res =   normal.quantile(Constants.rand.nextDouble()*quantile,
				meany,sigmay);
		//if(Constants.CHECK){
   		 if(Double.isNaN(res)){
   			 return min;
//   			 throw new RuntimeException("!!");
   		 }
   	 //}
		return res;
	}
	
	public double sampleGreaterThan(double topQuantile, double meany, double sigmay) {
		if(topQuantile==0) return max;
		double res = normal.quantile(1.0 - Constants.rand.nextDouble()*topQuantile,
				meany, sigmay);
		if(Double.isNaN(res)) return max;
		return res;
	}

	

	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, DoubleMatrix2D yb, int numObs,
			double[] noCop,  double pseudo) {
		// TODO Auto-generated method stub
		return dist.fill(x, y, yb, numObs, noCop, pseudo);
	}

	
	@Override
	public void getInterval(double[] input, DoubleMatrix2D res, double[] mean) {
		dist.getInterval(input, res, mean);
		
	}

	@Override
	public int getParamIndex() {
		return dist.getParamIndex();
	}

	@Override
	public String id() {
		// TODO Auto-generated method stub
		return dist.id();
	}

	@Override
	public void initialise() {
		dist.initialise();
		
	}

	@Override
	public void maximise(double d, double e, double f, double d1, double e1,
			double f1, double g) {
		dist.maximise(d, e, f, d1, e1, f1, g);
		if(this.topbottom){
			this.dist.setRho(0);
		}
		this.mlscale = Math.sqrt(dist.variance(dist.meanx,dist.meany,e,e1,g)[1]);
		//this.update();
		
	}
	final boolean topbottom;

	@Override
	public String name() {
		return dist.name();
	}

	@Override
	public int numObs() {
		return dist.numObs();
	}

	@Override
	public double probability(double x, double y) {
		double res = dist.probability(x, y);
		if(y==min){
			res += probLessThan(x,min, true);
		}
		else if(y==max){
			res+= probLessThan(x,max,false);
		}
		return res;/// (1-(probLessThan(x,min, true)+ probLessThan(x,max,false)));
	}
	
	double muy1, sigmay1,mux1,sigmax1,cdf;
	
	public double probLessThan(double x, double y, boolean lt){
		 muy1 = 
			dist.meany+dist.rho()*(dist.sigma_y/dist.sigma_x)*(x - dist.sigma_x);
		 sigmay1 = Math.sqrt(1 - Math.pow(dist.rho(), 2))*this.mlscale;
		mux1 = 
			dist.meanx+dist.rho()*(dist.sigma_x/dist.sigma_y)*(y - dist.sigma_y);
		sigmax1 = Math.sqrt(1 - Math.pow(dist.rho(), 2))*dist.sigma_x;
		 cdf = normal.cdf(y, muy1,sigmay1);
		return lt ? cdf*normal.pdf(x, mux1, sigmax1) : 
			(1-cdf) * normal.pdf(x, mux1, sigmax1);
	}
	
	static pal.statistics.NormalDistribution normal;

	@Override
	public void recalcName() {
		dist.recalcName();
	}

	@Override
	public void setParam(int type, int i, double d) {
		dist.setParam(type, i, d);
	//	this.update();
		
	}

	@Override
	public void transfer(double pseudoC) {
		dist.transfer(pseudoC);
	//	this.update();
		
	}

	@Override
	public void updateParamIndex() {
		dist.updateParamIndex();
		
	}

	@Override
	public void variance(int type, double[] sum) {
		dist.variance(type, sum);
		
	}

	@Override
	public int compareTo(Object arg0) {
		return dist.compareTo(((CutNormal21)arg0).dist);
	}

	@Override
	public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
			DoubleMatrix2D covar, int numObs,  double pseudo) {
		return dist.fillVariance(y, yb, covar, numObs, pseudo);
	}

	@Override
	public void print(PrintWriter pw) {
	this.dist.print(pw);
		
	}
	

	
	@Override
	public ProbabilityDistribution2 clone(double u,
			SimpleExtendedDistribution1 dist1) {
		throw new RuntimeException("!!");
	}
	@Override
	public double probability(double r, double b, int mixComponent) {
		// TODO Auto-generated method stub
		throw new RuntimeException("!!");
	}
	@Override
	public void setMinMax(double minR, double maxR, double min, double max) {
		throw new RuntimeException("!!");
		
	}
	@Override
	public void setToExclude() {
		throw new RuntimeException("!!");
		
	}
	public void addCount(Double r2,  Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG){
		throw new RuntimeException("!!");
	}
	
}
