/**
 * 
 */
package lc1.dp.transition;

import java.io.Serializable;
import java.util.Arrays;
import java.util.logging.Logger;

import lc1.stats.Dirichlet;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.TrainableGammaDistribution;
import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;


public class FreeRateTransitionProbs1  extends FreeTransitionProbs1 implements Serializable, UnivariateFunction{
final TrainableGammaDistribution sn,sn1;
	public double[] logrelativeRate;
	//public double logrelativeRateBase;
	double distance;
	MatrixExp global;
	int start = 1;
//	double pseudo;
	//final double sd = 10.0;
	//final double min;
	//final double max;
	public FreeRateTransitionProbs1(FreeExpTransitionProbs probs, double distance, TrainableGammaDistribution sn,TrainableGammaDistribution sn1 ) {
		super(probs);
		this.sn = sn;
		
		this.sn1 = sn1;
		logrelativeRate = new double[probs.mat.len];
		Arrays.fill(logrelativeRate,0);
		this.distance = distance;
		this.global = probs.mat;
	
	}
	


	public FreeRateTransitionProbs1(double rate, int[] is, int i, Dirichlet dirichlet,
			TrainableGammaDistribution sn, TrainableGammaDistribution sn1,
			double dist, MatrixExp global, double[] logrelativeRate) {
		super(is.length, dirichlet);
		//this.logrelativeRateBase = rate;
		this.stateToGroup = new int[is.length];
		Arrays.fill(stateToGroup,i-1);
		this.sn = sn;
		this.sn1 = sn1;
		this.logrelativeRate = logrelativeRate;
		
		this.distance = dist;
		this.start=0;
		this.global =global;
		// TODO Auto-generated constructor stub
	}



	double pseudo =0; 
	public double transferQ(double[] ds,double pseudoAlpha, double pseudoRate,  MatrixExp initial1, int i1, double distance) {
		 double sum=0;
		  for(int i=start; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
						global.setDistance(distance*Math.exp(logrelativeRate[i-start]));//+this.logrelativeRateBase));
					  int k = stateToGroup==null ? i-start : stateToGroup[i];
			          sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k),ds[i1]);
			  }
		  }
		  for(int i=0; i<start; i++){
			  if(transitionsOut[i]!=null){
				  sum+= ((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(ds[i1]);
			  }
		  }
		  return sum;
		//throw new RuntimeException("!!");
	}
	public static void transferQ(final FreeRateTransitionProbs1[] probs, double[] ds,
			double pseudoAlpha, double pseudoRate,  MatrixExp initial1, int i1, double distance) {
		try{
			final double[] logrelativeRate = probs[0].logrelativeRate;
			//final boolean hasswtch = logrelativeRate[0] > -4 && probs.length>1;
			/*if(probs.length>1){
			for(int i=0; i<logrelativeRate.length; i++){
				if(probs[1].transitionsOut[i].probs()[2] < 0.99){
					hasswtch = true;
				}
			}
			}*/
			 final UnivariateMinimum uvm= new UnivariateMinimum();//OrthogonalSearch();
			 for(int k=0; k<probs.length; k++){
				 probs[k].pseudo = ds[i1];
			 }
			
			final  TrainableGammaDistribution sn1 = probs[0].sn1;
			final  MatrixExp global = probs[0].global;
			if(true){
				for(int k=0; k<logrelativeRate.length; k++){
			
				final int k1 = k;
				final UnivariateFunction uvf = new UnivariateFunction(){

					@Override
					public double evaluate(double argument) {
					//System.arraycopy(argument, 	0,logrelativeRate,0,argument.length);// = argument;
						double sum=0;
						for(int kk=0; kk<probs.length; kk++){
							probs[kk].currIndex = k1;
						  sum+= probs[kk].evaluate(argument);
						}
						 double prior =0;// Math.log(sn.probability(Math.exp(logrelativeRateBase)*this.global.currentRate));
						//for(int k=0; k<argument.length; k++){
							prior+= Math.log(sn1.probability(Math.exp(argument)*global.currentRate));//normal.logpdf(argument, 0, sd);
						//}
						
						
					/*if(hasswtch){
							System.err.println(k1+" "+argument+" -> "+sum+"\t"+prior+"\t"+(sum-prior));
								}*/
						return sum - prior;
					}

					@Override
					public double getLowerBound() {
						return probs[0].getLowerBound();
					}

					@Override
					public double getUpperBound() {
						return probs[0].getUpperBound();
					}

/*					@Override
					public int getNumArguments() {
						return logrelativeRate.length;
					}

					@Override
					public OrthogonalHints getOrthogonalHints() {
						// TODO Auto-generated method stub
						return null;
					}*/
					
				};
				 Runnable run = new Runnable(){
		                public void run(){
		                	try{
		                		logrelativeRate[k1] = uvm.findMinimum(logrelativeRate[k1],uvf,1);
		                	}catch(Exception exc){
		                		exc.printStackTrace();
		                		
		                	}
		                }
		                };
		                Thread th = new Thread(run);
		                th.run();
		                try{
		                for(int i=0; i<100; i++){
		                    Thread.sleep(100);//this.wait(100);
		                    if(!th.isAlive()) break;
		                }
		                }catch(Exception exc){
		                	}
		               
				
			/*	 double rate = Math.exp(logrelativeRate[k])*probs[0].global.currentRate;
				 probs[0].sn1.addCount(rate, 1.0);*/
			}
				
			}
			 for(int k=0; k<probs.length; k++){
					probs[k].transferQ(ds, pseudoAlpha, pseudoRate, initial1, i1, distance);
				 }
		 
		}catch(Exception exc){
			exc.printStackTrace();
			//return Double.NEGATIVE_INFINITY;
		}
		}
	
	public double transfer1(double r, double pseudo){
		
		double sum=0;
		
		 for(int i=start; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
					  global.setDistance(distance * Math.exp(logrelativeRate[i-start]+r));//+this.logrelativeRateBase));
					  int k = stateToGroup==null ? i-start
							 : stateToGroup[i];
			       sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k),pseudo);
			  }
		  }
		 for(int i=0; i<start; i++){
			  if(transitionsOut[i]!=null){
				  sum+= ((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(pseudo);
			  }
		 }
		  return sum;
		
	}



	
	
	 int currIndex=0;
	 @Override
	public double evaluate(double argument) {
		try{
	//	this.logrelativeRate = argument;
	
		  double sum=0;
		  int i = currIndex;
		//  for(int i=start; i<this.transitionsOut.length; i++){
			  if(transitionsOut[i]!=null){
					global.setDistance(distance*Math.exp(argument));//+this.logrelativeRateBase));
					
					  int k = stateToGroup==null ? i-start : stateToGroup[i];
			       sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(global.M.viewRow(k), pseudo);
			  }
		  //}
		 
		  return -1*(sum);// +prior;
		}catch(Exception exc){
			Logger.global.warning("problem with distance "+logrelativeRate+ " "+argument+" "+exc.getMessage());
			return Double.POSITIVE_INFINITY;
		}
	}
	
	@Override
	/*public double evaluate(double argument) {
		this.logrelativeRateBase = argument;
		return this.evaluate(this.logrelativeRate);
	}*/
	
	
	public double getRate(int i) {
		// TODO Auto-generated method stub
		return Math.exp(this.logrelativeRate[i])//+this.logrelativeRateBase)
		*this.global.currentRate;
	}



	/*@Override
	public double getLowerBound() {
	return -10;
	//Math.log(sn.inverse(1e-5)/ this.global.currentRate);
	}



	@Override
	public double getUpperBound() {
		return 10;
//		return Math.log(sn.inverse(1-1e-5)/ this.global.currentRate);
	}

	@Override
	public int getNumArguments() {
		return this.logrelativeRate.length;
	}

	@Override
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}
*/
	

	//@Override
	public double getLowerBound() {
		return -10;
	}

	//@Override
	public double getUpperBound() {
		return 10;
	}

	

	
	  
	
   
       
      
}