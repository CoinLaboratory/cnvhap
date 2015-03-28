package lc1.stats;

import java.io.PrintWriter;
import java.util.logging.Logger;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix2D;

public class CutNormal  implements ProbabilityDistribution, Comparable{
	
	//double lessThanZeroProb;
//	 double greaterThanOneProb;
	double mlscale;
	final TrainableNormal dist;
	final double min;
	final double max;
	public void setPriors(ProbabilityDistribution distx, int type){
		dist.setPriors(((CutNormal)distx).dist, type);
	}
	public double  probability(double x, int mixComponent){
		return this.probability(x);
	}
	public void setCoverage(double d){
    	throw new RuntimeException("!!");
    }

	public CutNormal(TrainableNormal dist){
		this(dist, 0.0, 1.0);
		this.quantileG = this.greaterThanMaxProb();
		
		this.quantileL = this.lessThanMinProb();
		this.quantileRem = 1 - (quantileL+quantileG);
	}
	
	 public CutNormal(TrainableNormal dist, double min, double max) {
		 if(dist.scale==0) {
			 throw new RuntimeException("!!");
		 }
		this.dist = dist;
		this.min = min;
		this.max = max;
		this.mlscale = dist.scale;///5.0;
	this.quantileG = this.greaterThanMaxProb();
		
		this.quantileL = this.lessThanMinProb();
		this.quantileRem = 1 - (quantileL+quantileG);
		//this.update();
	}
	
	public CutNormal(CutNormal cutNormal, double u) {
	//	this.lessThanZeroProb = cutNormal.lessThanZeroProb;
	//	this.greaterThanOneProb = cutNormal.greaterThanOneProb;
		this.dist = cutNormal.dist.clone(u);
		this.min = cutNormal.min;
		this.max = cutNormal.max;
		this.mlscale = cutNormal.mlscale;
		this.quantileG = cutNormal.quantileG;
		this.quantileRem = cutNormal.quantileRem;
		this.quantileL = cutNormal.quantileL;
	}

	public CutNormal(CutNormal cutNormal) {
		//this.lessThanZeroProb = cutNormal.lessThanZeroProb;
		//this.greaterThanOneProb = cutNormal.greaterThanOneProb;
		this.dist = cutNormal.dist.clone();
		this.min = cutNormal.min;
		this.max = cutNormal.max;
		this.mlscale = cutNormal.mlscale;
		this.quantileG = cutNormal.quantileG;
		this.quantileRem = cutNormal.quantileRem;
		this.quantileL = cutNormal.quantileL;
	}
//public void update(){
	//this.lessThanZeroProb = lessThanMinProb();
	//this.greaterThanOneProb = greaterThanMaxProb();
	//}
public double lessThanMinProb(){
	double res =  dist.normal.cdf(min, dist.location, this.mlscale);
	return res;
}
public double greaterThanMaxProb(){
	double res =  1.0-  dist.normal.cdf(max, dist.location, this.mlscale);
	return res;
}
	public CutNormal clone(){
	        return new CutNormal(this);
	    }
	
	public CutNormal clone(double u){
        return new CutNormal(this, u);
    }
//boolean upper  = false; // 

	
	double quantileL;
	double quantileG;
	double quantileRem;
	
/*	public void addCount1(double obj_index, double value) {
		quantileL = this.lessThanMinProb();
		quantileG = this.greaterThanMaxProb();
		quantileRem = 1-(quantileL+quantileG);
	   if(quantileL>1e-10){
				double v = obj_index==min ? obj_index : sampleLessThan(quantileL);
				this.dist.addCount(v, value*quantileL);
		}
		if(quantileG>1e-10){
				double v = obj_index==max ? obj_index : sampleGreaterThan(quantileG);
				dist.addCount(v, value*quantileG);
		}
		this.dist.addCount(obj_index, value*(quantileRem));
		
		
	}*/
	
	
	//List<Double> maxToAdd =new ArrayList<Double>();
	//List<Double> minToAdd =new ArrayList<Double>();;
	public static boolean asCut =true;
	@Override
	public void addCount(double obj_index, double value) {
		//this.sumMin+=value;
		/*this.quantileL = this.lessThanMinProb();
		this.quantileG =  this.greaterThanMaxProb();
		this.quantileRem = 1 - (quantileL + quantileG);*/
		/*if(rem==0){
			throw new RuntimeException("!!");
		}*/
		/*
		if(asCut && obj_index<=min){
			this.minToAdd.add(value);
			sumMin+=value;
		}
		else if(asCut && obj_index>=max){
			this.maxToAdd.add(value);
			sumMax+=value;
			
		}
	
	
		else{*/
	
		
			//double qv = quantile
			if(quantileL*value>Constants.countThresh2()){
			double lowerq = this.sampleLessThan(this.quantileL);
			this.dist.addCount(lowerq, quantileL*value);
			}
			if(quantileG*value > Constants.countThresh2()){
			 double upperq = this.sampleGreaterThan(this.quantileG);
		      this.dist.addCount(upperq, quantileG*value );
			}
			 if(quantileRem*value > Constants.countThresh2()){
			this.dist.addCount(obj_index,quantileRem*value);
			}
		
		//}
		
	}
	
	public double sampleLessThan(double quantile) {
		double rand = Constants.rand.nextDouble();
		
		double res =  dist.normal.quantile(rand*quantile,
				dist.location, mlscale);
		if(Constants.CHECK && Double.isNaN(res)){
			Logger.global.info("IS NAN "+res);
		}
		return res;
	}
	
	public double sampleGreaterThan(double topQuantile) {
		double res =  dist.normal.quantile(1.0 - Constants.rand.nextDouble()*topQuantile,
				dist.location, mlscale);
		return res;
	}

	@Override
	public void addCounts(ProbabilityDistribution probabilityDistribution) {
		throw new RuntimeException("!!");
		
	}

	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs,
			double []noCop,  double pseudo) {
		// TODO Auto-generated method stub
	//	this.addRemainder();
		return dist.fill(x, y, numObs, noCop, pseudo);
	}

	@Override
	public int fillVariance(DoubleMatrix2D y, int numObs, double pseudo) {
		// TODO Auto-generated method stub
		//if(this.addRemainder()) throw new RuntimeException("!!");
		return dist.fillVariance(y, numObs, pseudo);
	}

	@Override
	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void getInterval(double[] in, double[] res_inR) {
		this.dist.getInterval(in, res_inR);
		
	}

	@Override
	public double getMean() {
		// TODO Auto-generated method stub
		return this.dist.getMean();
	}

	@Override
	public int getParamIndex() {
		return dist.getParamIndex();
	}

	@Override
	public double getParamValue(int n1) {
	return dist.getParamValue(n1);
	}

	@Override
	public String id() {
		return dist.id();
	}

	@Override
	public void initialise() {
		dist.initialise();
	//	this.maxToAdd.clear();
	//	this.minToAdd.clear();
		//sumMin =0;
		//sumMax =0;
	}

	/*public boolean addRemainder(){
		if(maxToAdd.size()>0 || minToAdd.size()>0){
			// double sumMin = sum(minToAdd);
		  //      double sumMax = sum(maxToAdd);
		double avg = average(1e-3, sumMin, sumMax);
		this.mlscale = Math.sqrt(variance(avg, 1e-3, sumMin, sumMax));
		for(int i=0; i<minToAdd.size(); i++){
			double obj_index = min;
			double value = minToAdd.get(i);
			quantileL = this.lessThanMinProb();
			if(quantileL>1e-10){
				for(int i1=0; i1<Constants.numToSample(); i1++){
				double v = sampleLessThan(quantileL);
				this.dist.addCount(v, value/(double)Constants.numToSample());
				}
			}
			else{
				dist.addCount(obj_index, value);
			}
		}
		this.minToAdd.clear();
		for(int i=0; i<maxToAdd.size(); i++){
			double obj_index = max;
			double value = maxToAdd.get(i);
			quantileG = this.greaterThanMaxProb();
			if(quantileG>1e-10){
				for(int i1=0; i1<Constants.numToSample(); i1++){
					double v = sampleGreaterThan(quantileG);
					dist.addCount(v, value/(double)Constants.numToSample());
				}
			}
			else{
				dist.addCount(obj_index, value);
			}
		}
		
		this.maxToAdd.clear();
		
		sumMin=0; 
		sumMax=0;
		return true;
		}
		return false;
	}*/
	
	 public double average(double pseudo, double sumMin, double sumMax){
	        double tot =pseudo +sumMin+sumMax;
	        double exp =pseudo * dist.meanPrior+sumMin*min + sumMax*max;
	       
	      
	        for(int i=0; i<dist.obsx.size(); i++){
	           
	            double sc = dist.obsx.get(i);
	            double w = dist.obsv.get(i);
	           
	            exp+=sc*w;
	            tot+=w;
	        }
	       // if(Math.abs(mean- Constants.r_mean()[3])<0.01){
	        //   System.err.println(this.meanPrior.getMean()+" rvalues "+observations);
	        //   System.err.println(this.meanPrior.getMean()+" rvalues "+weights);
	       // }
	        double res =  exp/tot;
	        if(Constants.CHECK && Double.isNaN(res)){
	        	throw new RuntimeException("!!");
	        }
	        return res;
	    }
	 
	 public double variance(double mu, double pseudo, double sumMin, double sumMax){
	    	
	        double tot =pseudo +sumMin+sumMax;
	      double exp =Math.pow(dist.stddevPrior,2)*pseudo + sumMin * Math.pow(min - mu, 2) +
	      sumMax * Math.pow(max - mu, 2);
	      for(int i=0; i<dist.obsx.size(); i++){
	      
	          double sc =Math.pow(dist.obsx.get(i)-mu, 2);
	          double w =dist.obsv.get(i);
	            exp+=sc*w;
	            tot+=w;
	          //  if(exp < 0 || tot < 0) throw new RuntimeException("!! "+w+" "+sc);
	        }
	        if(Constants.CHECK && Double.isInfinite(exp/tot) || tot==0){
	            throw new RuntimeException(" ");//+this.observations+"\n"+this.observations);
	        }
	       
	        double res =  exp/tot;
	        if(Constants.CHECK && (res==0  || Double.isNaN(res))){
	        	throw new RuntimeException("!!");
	        }
	        return res;
	    }
	
	/*private double sum(List<Double> minToAdd2) {
		double sum=0;
		for(int i=0; i<minToAdd2.size(); i++){
			sum+=minToAdd2.get(i);
		}
		return sum;
	}*/
	// double sumMin;//, sumMax;


	@Override
	public void maximise(double d, double e, double f) {
	/*	this.addRemainder();
		if(this.dist.location<0.05 || this.dist.location>0.95){
		//	System.err.println("skipping");
			String bef = dist.location+" "+dist.scale+ " "+dist.meanPrior+" "+dist.stddevPrior;
			dist.maximise(d, e, f);
		String aft = dist.location+" "+dist.scale+ " "+dist.meanPrior+" "+dist.stddevPrior;
			if(!bef.equals(aft)){
				System.err.println("BEF/AFT"+dist.id+"\n"+bef+"\n"+aft);
		}
//			System.err.println("after "+dist.location+" "+dist.scale+ " "+dist.meanPrior[0]+" "+dist.stddevPrior[0]);
			
		}
		else{*/
			dist.maximise(d, e, f);
	//	}
		this.mlscale =dist.scale;//Math.sqrt( dist.variance(dist.location, 1e-5));
	this.quantileG = this.greaterThanMaxProb();
		
		this.quantileL = this.lessThanMinProb();
		this.quantileRem = 1 - (quantileL+quantileG);
		//	if(this.greaterThanMaxProb()+this.lessThanMinProb()>=1.0){
		//	throw new RuntimeException("!!");
	//	}
	}

	@Override
	public String name() {
		return dist.name();
	}

	@Override
	public int numObs() {
	return dist.numObs();//+this.maxToAdd.size()+this.minToAdd.size();
	}

	@Override
	public double plotObservations(String string, boolean b, XYSeries obs,
			boolean swtch) {
		return dist.plotObservations(string, b, obs, swtch);
	}

	@Override
	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		dist.plotTheoretical(string,b, theor);
		
	}

	@Override
	public double prior() {
		return dist.prior();
	}

	
	public double probability1(double x) {
		double res = this.dist.probability(x);
		if(res==0) return 0;
		double lt = this.lessThanMinProb() + this.greaterThanMaxProb();;
		if(lt>=1.0){
			throw new RuntimeException("!!");
		}
		double res1 =  res/(1-lt);
		if(Constants.CHECK && Double.isNaN(res1)){
			throw new RuntimeException("!!");
		}
		return res1;
	}
	
	
	@Override
	public double probability(double x) {
		
		double res = 
			 this.dist.probability(x)/this.quantileRem;
		if(x < min || x> max || Double.isNaN(res)){
			System.err.println("WARNING no remainder");
			res = 0.0;
		}
		/*if(Math.abs(x-0.0)<0.01 && this.dist.location>=0.99 ){
			System.err.println("res" + res);
		}*/
		//if(true) return res;
	//	if(asCut && x<=min) return res+this.lessThanMinProb();
	//	else if(asCut && x>=max) return res+this.greaterThanMaxProb();
//		else
			return res;
		/*if(res==0) return 0;
		double lt = this.lessThanMinProb() + this.greaterThanMaxProb();;
		if(lt>=1.0){
			throw new RuntimeException("!!");
		}
		double res1 =  res/(1-lt);
		if(Constants.CHECK && Double.isNaN(res1)){
			throw new RuntimeException("!!");
		}
		return res1;*/
	}

	@Override
	public void recalcName() {
		dist.recalcName();
		
	}

	@Override
	public double sample() {
		double d = dist.sample();
		if(d<=this.min) return min;
		else if(d>=this.max) return max;
		else return d;
	}

	@Override
	public double scale() {
		return dist.scale();
	}

	@Override
	public void setParam(int type, double rho) {
		//if(type==0 && this.dist.location >0.95){
		//	System.err.println("h");
		//}
		//if(this.dist.location<0.05 || this.dist.location>0.95) return;
		dist.setParam(type, rho);
	this.quantileG = this.greaterThanMaxProb();
		
		this.quantileL = this.lessThanMinProb();
		this.quantileRem = 1 - (quantileL+quantileG);
		//this.update();
	}

	
	
	@Override
	public void setParamValue(int n1, double val) {
		throw new RuntimeException("!!");
//		dist.setParamValue(n1, val);
	//	this.update();
	}

	@Override
	public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
		throw new RuntimeException("!!");
		
	}

	@Override
	public double sum() {
		return dist.sum;
	}

	@Override
	public void transfer(double pseudoC) {
		throw new RuntimeException("!!");
//		dist.transfer(pseudoC);
		//this.update();
		
	}

	@Override
	public void transfercounts(EmissionState innerState, int phen_index, int i) {
	throw new RuntimeException("!!");
		
	}

	@Override
	public void updateParamIndex() {
		dist.updateParamIndex();
		
	}

	@Override
	public void variance(double[] sum) {
		dist.variance(sum);
		
	}
	
	public String toString(){
		return this.dist.toString();
	}

	@Override
	public int compareTo(Object arg0) {
		return dist.compareTo(((CutNormal)arg0).dist);
	}
	@Override
	public double evaluate(double[] argument) {
	return dist.evaluate(argument);
	}

	@Override
	public double getLowerBound(int n) {
		// TODO Auto-generated method stub
		return dist.getLowerBound(n);
	}

	@Override
	public int getNumArguments() {
		return dist.getNumArguments();
	}

	@Override
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getUpperBound(int n) {
		return dist.getUpperBound(n);
	}
	@Override
	public void print(PrintWriter pw) {
		this.dist.print(pw);
	}
	public void setMinMax(double min, double max){
		this.dist.setMinMax(min, max);
	this.quantileG = this.greaterThanMaxProb();
		
		this.quantileL = this.lessThanMinProb();
		this.quantileRem = 1 - (quantileL+quantileG);
	}
	public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, ProbabilityDistribution disty){
		disty.addCount(b, val);
    	mixe1.addCount(0, val);
    	throw new RuntimeException("!!");
	}
	@Override
	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1) {
		// TODO Auto-generated method stub
		return this.clone(u);
	}
}
