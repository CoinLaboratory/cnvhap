/**
 * 
 */
package lc1.dp.transition;

import java.io.Serializable;
import java.util.logging.Logger;

import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.TrainableGammaDistribution;
import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;


public class FreeRateTransitionProbs  extends FreeTransitionProbs1 implements Serializable, UnivariateFunction{
final TrainableGammaDistribution sn;
	public double logrelativeRate=0.0;
	double distance;
	MatrixExp global;
//	double pseudo;
	//final double sd = 10.0;
	//final double min;
	//final double max;
	public FreeRateTransitionProbs(FreeExpTransitionProbs probs, double distance, TrainableGammaDistribution sn ) {
		super(probs);
		this.sn = sn;
		this.distance = distance;
		this.global = probs.mat;
	
	}
	
	double pseudo =0; 
	@Override
	public double transferQ(double[] ds,double pseudoAlpha, double pseudoRate,  MatrixExp initial1, int i1, double distance, int it) {
		//  if(ds[i1]<1e10){
		 // if(true) return 0;
		try{
			  UnivariateMinimum uvm= new UnivariateMinimum();
			 pseudo = ds[i1];
			  logrelativeRate = uvm.findMinimum(this.logrelativeRate, this,5);
			  if(logrelativeRate > 5){
				  Logger.global.info("rr "+logrelativeRate);
			  }
			 // 
			  double rate = Math.exp(logrelativeRate)*this.global.currentRate;
			  this.sn.addCount(rate, 1.0);
		//  }
		global.setDistance(distance*Math.exp(logrelativeRate));
		//System.err.println("pi "+initial.pi);
		  double sum=0;
		  for(int i=0; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
				  if(i>0){
					  int k = stateToGroup==null ? i-1 : stateToGroup[i]-1;
			       sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k),ds[i1]);
				  }
				  else{
					  sum+= ((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(ds[i1]);
				  }
			 
			  }
		  }
		  return sum;
		}catch(Exception exc){
			exc.printStackTrace();
			return Double.NEGATIVE_INFINITY;
		}
		}
	
	public double transfer1(double r, double pseudo){
		
		double sum=0;
		global.setDistance(distance * Math.exp(logrelativeRate+r));
		 for(int i=0; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
				  if(i>0){
					  int k = stateToGroup==null ? i-1 : stateToGroup[i]-1;
			       sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k),pseudo);
				  }
				  else{
					  sum+= ((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(pseudo);
				  }
			 
			  }
		  }
		  return sum;
		
	}



	@Override
	public double evaluate(double argument) {
		try{
		this.logrelativeRate = argument;
		global.setDistance(distance*Math.exp(logrelativeRate));
		  double sum=0;
		  for(int i=1; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
					  int k = stateToGroup==null ? i-1 : stateToGroup[i]-1;
			       sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k), pseudo);
			  }
		  }
		  double prior = Math.log(sn.probability(Math.exp(argument)*this.global.currentRate));//normal.logpdf(argument, 0, sd);
		  //System.err.println(argument+" "+(-1*sum)+" "+(-1*prior)+" "+(-1*(sum+prior)));
			
		  return -1*(sum +prior);
		}catch(Exception exc){
			Logger.global.warning("problem with distance "+logrelativeRate+ " "+argument+" "+exc.getMessage());
			return Double.POSITIVE_INFINITY;
		}
	}
	@Override
	public double getRate(int i) {
		// TODO Auto-generated method stub
		return Math.exp(this.logrelativeRate)*this.global.currentRate;
	}



	@Override
	public double getLowerBound() {
	return -20;
	//Math.log(sn.inverse(1e-5)/ this.global.currentRate);
	}



	@Override
	public double getUpperBound() {
		return 20;
//		return Math.log(sn.inverse(1-1e-5)/ this.global.currentRate);
	}

	

	
	  
	
   
       
      
}