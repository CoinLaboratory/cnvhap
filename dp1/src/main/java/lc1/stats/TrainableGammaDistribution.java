package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.logging.Logger;

import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;

public class TrainableGammaDistribution implements UnivariateFunction, Serializable{

	
	double k, theta;
//	List<Double> obsx  = new ArrayList<Double>();
//	List<Double> obsv = new ArrayList<Double>();
	 double sum=0;
	 double sumln=0;
	 double N =0;
	public void initialise(){
	//	obsx.clear();
		//obsv.clear();
		sum=0;
		sumln=0;
		N=0;
	}
/*	public TrainableGammaDistribution(double shape, double scale){
		//gamma = new Gamma(shape, scale, null);
		this.k = shape;
		this.theta = scale;
		
	}*/
	
	public TrainableGammaDistribution(double shape, double scale, double mode){
		//gamma = new Gamma(shape, scale, null);
		this.k = shape;
		this.theta = scale;
		//if(shape>0){
			//theta = mode/(k-1);
			k = (mode/theta) +1.0;
		//}
			//double mde = this.inverse(0.5);
			//double cdf = this.cdf(mde);
			//double mod = (k-1)*theta;
		Logger.global.info(k+" "+theta);
		
	}
	
	public TrainableGammaDistribution(double mode, double factor){
		this.theta = mode/factor;
		k = (mode/theta) +1.0;
		Logger.global.info("new Gamma "+k+" "+theta);
	}
	private double cdf(double mde) {
		return GammaDistribution.cdf(mde, k, theta);
	}
	public TrainableGammaDistribution(double[] gammaRate, double r) {
		this(gammaRate[0], gammaRate[1], r);
	}
	public TrainableGammaDistribution(double[] gammaRate) {
		this(gammaRate[0], gammaRate[1]);
	}
	public double inverse(double y){
		return GammaDistribution.quantile(y, k, theta);
	}
	public double probability(double x){
		
		return GammaDistribution.pdf(x, k, theta);
	}
 public synchronized void addCount(double d, double w) {
	// this.obsx.add(d);
	// this.obsv.add(w);
	 sum+=d*w;
	 sumln+= Math.log(d*w);
	 N+=w;
 }
 public void maximise(){
	 double x = this.initv();
	 UnivariateMinimum uvm = new UnivariateMinimum();
	 this.k = uvm.findMinimum(x, this, 5);
	 this.theta = theta();
	 
	 
 }
 public String toString(){
	 return this.k+":"+this.theta;
 }

@Override
public double evaluate(double argument) {
//	double N = obsx.size();
	double res = (k-1)*sumln - N * k - N*k*Math.log(sum/(k*N))
	 - N * GammaFunction.lnGamma(k);
	return -1 * res;
}
public double initv(){
	//double N = obsx.size();
	double s = Math.log(sum/N) - sumln/N;
	return (3 -s + Math.sqrt(Math.pow(s-3, 2)+24*s))/12*s;
}
public double theta(){
	return sum / (k*N);
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
	pw.println("--gammaRate  "+String.format("%5.3g", this.k).trim()+":"
			+String.format("%5.3g",this.theta).trim());
	pw.println("--gammaRateG  "+String.format("%5.3g", this.k).trim()+":"
			+String.format("%5.3g",this.theta).trim());
	
}
}
