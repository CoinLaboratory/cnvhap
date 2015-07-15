/**
 * 
 */
package lc1.dp.model;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.states.DotState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.BetweenWithinTransitionProbs1;
import lc1.dp.transition.ExpTransProb;
import lc1.dp.transition.ExponentialTransitionProbs;
import lc1.dp.transition.FreeExpTransitionProbs;
import lc1.dp.transition.FreeRateTransitionProbs;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.dp.transition.MatrixExp;
import lc1.dp.transition.MultiExpProbs;
import lc1.stats.Dirichlet;
import lc1.stats.NormalDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.TrainableGammaDistribution;
import lc1.util.Constants;
import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;
import cern.colt.matrix.DoubleMatrix2D;

public abstract class SiteTransitions implements Serializable,UnivariateFunction{
  
	
	FreeExpTransitionProbs globalTrans;
	MatrixExp initial;
	
	private void makeRateDists(double rate, double gamma1, double gamma2){
		//(Constants.probCrossOverBetweenBP * Constants.expModelIntHotSpot[0])
		
		rateDistribution = new TrainableGammaDistribution(
				
				rate,
				gamma1 );
		rateDistributionG = new TrainableGammaDistribution(
				rate,
				gamma2 );
	}
	 TrainableGammaDistribution rateDistribution, rateDistributionG;
	/*= new TrainableGammaDistribution(
		
			(Constants.probCrossOverBetweenBP * Constants.expModelIntHotSpot[0])
			                          ,
			(Constants.gammaRate(false) ));
	
	TrainableGammaDistribution rateDistributionG = new TrainableGammaDistribution(
		
			(Constants.probCrossOverBetweenBP * Constants.expModelIntHotSpot[0])
			                          ,
			(Constants.gammaRate(true) ));*/

	//new CopyOfSkewNormal("rates",0.0,10.0,0.0, -15, 15, 1000, 1);
	
	
    final int len;
    final double minRate = -15;
	
    final double maxRate = 15;
   // double[] pseudoCExp;
    @Override
	public double getLowerBound() {
	return minRate;
	}
	@Override
	public double getUpperBound() {
	return maxRate;
	}
	
	public double evaluate(double logr){
		  try{
			
		double sum=0;
		//globalTrans.mat.
		for(int i=1; i<this.transProbs.length; i++){
			sum+=((FreeRateTransitionProbs)this.transProbs[i]).transfer1(logr, 1e10);
		}
		//
		double rate = Math.exp(logr)*this.globalTrans.mat.currentRate;
		 double prior = Math.log(this.rateDistributionG.probability(rate));//normal.logpdf(argument, 0, sd);
		// System.err.println("eval g "+logr+" "+(-sum)+" "+(-prior));
		return -(sum + prior);
		  }catch(Exception exc){
				Logger.global.warning("problem with distance "+logr+" "+exc.getMessage());
				return Double.POSITIVE_INFINITY;
			}
	}
	public void maximiseGlobalRate(){
		UnivariateMinimum uvm = new UnivariateMinimum();
		double min  = uvm.findMinimum(0, this, 5);
		double newRate = Math.exp(min)*globalTrans.mat.currentRate;
		Logger.global.info("new rate1 "+min+" "+newRate);
		globalTrans.mat.resetRate(newRate);
		globalTrans.mat.update();
	}
    
    public String info() {
       if(transProbs.length>1) return this.transProbs[1].info();
       else return transProbs[0].info();
      }
     public AbstractTransitionProbs[] transProbs; //note offset: transProbs[i] is for transition from snp i-1 to snp i;
     final List<Integer> loc;
    // final List<State> states;
     final int numFounders;
     final Double[] exp_p1;
    protected Double[][] r; //should call this theta
   // protected Gamma[][] dist; // this is prior distributions on r
   //  Double[] r_prev;//=new double[] {0};
    // double[] r
   //  final Class transProbType;
     public Dirichlet getExp_p(int i){
    	// if(true) throw new RuntimeException("!!");
         if(loc==null && exp_p1==null) return null;
         double exp_p;
         if(loc==null || loc.size()==0){
        	 exp_p = exp_p1==null ? null : exp_p1[0];
        	 
         }
         else{
        	 double dist =  loc.get(i)-loc.get(i-1);
        	 exp_p = 
             
                 (r==null ? null :Math.exp(-1*r[0][0]*(dist))
                         );
         }
        
        
              return new Dirichlet(new Double[] {exp_p, 1-exp_p}, Constants.u_global[2]);
      }
     public double getDist(int i){
    	 return this.loc.get(i) - this.loc.get(i-1);
     }
      protected double[][] alpha_overall;
    final Integer[] cn;
      public SiteTransitions(int len, int numStates, Integer[] cn) throws Exception{
         this.len = len;
         this.makeRateDists((Constants.probCrossOverBetweenBP * Constants.expModelIntHotSpot[0]), 
        		 Constants.gammaRate(0), Constants.gammaRate(1));
         this.cn = cn;
         this.transLast = makeTransLast(cn, numStates);
         this.cached = false;
         this.transProbs = new AbstractTransitionProbs[len];
         int numF = numStates-1;
         double[] u = new double[numF+1];
         Arrays.fill(u, 1.0/ (double) numF);
         u[0] = 0;
         Dirichlet dir = new Dirichlet(u, Double.POSITIVE_INFINITY);
         transProbs[0] =new FreeTransitionProbs1(true, dir, numF+1);
         for(int i=1; i<len; i++){
             transProbs[i] =new FreeTransitionProbs1(false, dir, numF+1);
                 ((FreeTransitionProbs1)transProbs[i]).transitionsOut[0] = null;
         }
         this.loc =null;
         this.numFounders = numF;
    //     this.states = null;
         this.exp_p1 = null;
         this.r = null;
        // this.dist = null;
         alpha_overall = new double[1][numStates];
     }
     //transProbs[0] is from magic state, and transProbs[noSnps] is transition to magic state
      public static Integer[] getCN(List<State> states){
    	 Integer[] cn = new Integer[states.size()];
         for(int i=0; i<states.size(); i++){
        	 if(states.get(i) instanceof DotState){
        		 cn[i] = -1;
        	 }
        	 else{
        		 cn[i] = ((EmissionState)states.get(i)).noCop();
        	 }
         }
         return cn;
     }
      
      public SiteTransitions(List<Integer>loc, int numF ,Double[] exp_p1,Double[] r,  int length, int index, Integer[] cn, boolean makeGlobal){
        if(r!=null && exp_p1!=null) throw new RuntimeException("should not both be non null");
         this.transProbs = new AbstractTransitionProbs[length];
         this.makeRateDists((Constants.probCrossOverBetweenBP * Constants.expModelIntHotSpot[0]), 
        		 Constants.gammaRate(0), Constants.gammaRate(1));
         this.transLast = makeTransLast(cn, numF);
         this.loc = loc;
         this.len  = numF+1;
         this.cn = cn;
         this.cached = false;
      //   this.states = states;
         this.numFounders =numF;
         this.exp_p1 = exp_p1;
         int[] transMode0 = Constants.transMode0;
         if(r!=null){
         this.r = new Double[transMode0.length==1 ? 1 : 2][];
       //  this.dist = new Gamma[transMode0.length==1 ? 1 : 2][];
         for(int i=0; i<this.r.length; i++){
        	 this.r[i] = new Double[(transMode0[i]==3 || transMode0[i]==4) ? Constants.expModelIntHotSpot1(i).length:1];
        	// this.dist[i] = new Gamma[this.r[i].length];
        	 this.r[i][0] = r[i]; // out of start?
        	// this.dist[i][0] = makeDist(this.r[i][0],index);
        	 if(r!=null && this.r[i].length>1){
        		 double[] hotspot1 = Constants.expModelIntHotSpot1(i);
     	     	for(int k=0; k<this.r[i].length; k++){
     	     		
     	     	
     	     			this.r[i][k] = hotspot1[k]*r[i];
     	     		
     	     	}
          	
            }
         }	
         if(makeGlobal){
         double[][] trans = Constants.transitionMatrix();
         this.globalTrans = 
        	 trans!=null ? 
        	 new FreeExpTransitionProbs(0, trans, this.r[0][0]):
        	 
        	 new FreeExpTransitionProbs(Constants.modifyFrac(0),  Constants.expModelIntHotSpot1(0),this.r[0][0]);
         this.initial = (MatrixExp) globalTrans.mat.clone();
         }
         
        		// if((transMode0[i]==3 || transMode0[i]==4)){
                	// if(Constants.expModelIntHotSpot1(i).length!=(Constants.modify0.length+1)) throw new RuntimeException("mismatch between expModelIntHotSpot1 and modify0");
                // }
         }
        
       
        		
     
     	
  
     }
   
	
	public boolean converged(){
         return true;
     }
     
 
   
     public abstract SiteTransitions clone(boolean swtch);
     
     /** if pseudocount=0 we do not use this as pseudo */
     public SiteTransitions(SiteTransitions trans_init, boolean swtch){
         this(trans_init.loc,trans_init.numFounders, trans_init.exp_p1, getFirst(trans_init.r),  trans_init.transProbs.length,0, trans_init.cn, true);
         transProbs[0] = trans_init.transProbs[0].clone(swtch);
         for(int i=1; i<transProbs.length; i++){
             transProbs[i] =new FreeTransitionProbs1( trans_init.transProbs[i]);//.clone(swtch);
         }
       //  this.MAGIC = states.get(0);
     }
      private static Double[] getFirst(Double[][] r2) {
		Double[] res = new Double[r2.length];
		for(int i=0; i<res.length; i++){
			res[i] = r2[i][0];
		}
		return res;
	}
	public void setTransitionScore(int from,  int to, int indexOfToEmission, double d, double d_ps) {
         FreeTransitionProbs1 transProbs =  (FreeTransitionProbs1) this.transProbs[indexOfToEmission];
         transProbs.setTransitionScore(from, to, d, transProbs.length());
     }
     
  
     public void initialiseTransitionCounts() {
    	 if(globalTrans!=null){
    		 this.rateDistribution.initialise();
    		 this.globalTrans.initialiseCounts(false, false);
    	 }
         for(int i=0; i< this.transProbs.length; i++){
            if(transProbs[i]!=null) transProbs[i].initialiseCounts(i==0, i==transProbs.length-1);
         }
     }

    
     
   
    public Collection<PseudoDistribution> getDistributions(int i){
        return this.transProbs[i].getDistributions();
    }
  // double[] trans;
   static double[] init_trans = new double[] {100,100};
   double r_prev_prev =0;
  
//private   double logLAll =0;
// double[] pseudoTransFirst = new double[2];  
 boolean allFree = false; //if all of transition probs are FreeTransitionProb1 or are BetweenWithinTransProbs1 with between groups as free.
 //double[] pseudoTransFirst = new double[2];
 
 public final static boolean hierarchy = false;
 NormalDistribution normal;
 
 
 public void transferEquilToStart(){
	 if(!Constants.transferEquilToStart(transProbs.length)) return;
	 DoubleMatrix2D nullspace = this.globalTrans.mat.pi;
		for(int i=0; i<nullspace.rows(); i++){
			((FreeTransitionProbs1)this.transProbs[0]).transitionsOut[0].setProbs(i+1, nullspace.getQuick(i, 0));
		}
 }
 

 public double transferTransitions(final double[] pseudoTrans,final double[] pseudoTransFirst,   final int index){
	double logLAll=0;
	
	//  logLAll+=transProbs[0].transferAlpha(pseudoTransFirst, nu, 0); 
	if(Constants.measureGlobal() && globalTrans!=null ){
    	  this.globalTrans.transferQ(pseudoTrans,0,0, initial, 3,1000,index);
    	   if(false) {
    		   this.maximiseGlobalRate();
    	   }
    	   this.makeRateDists(globalTrans.mat.currentRate, Constants.gammaRate(0), Constants.gammaRate(1));
    	   this.rateDistribution = new TrainableGammaDistribution(
    				//Constants.gammaRate(false), 
    				globalTrans.mat.currentRate,
    				(Constants.trainGlobal ?Constants.gammaRate(0) : Constants.gammaRate(1)));
	}
	if(!this.allFree) {
		if(Constants.transferEquilToStart(this.transProbs.length))
		this.transferEquilToStart();
		else{
			SimpleExtendedDistribution dist1 = ( (SimpleExtendedDistribution)
					 ((FreeTransitionProbs1)transProbs[0]).transitionsOut[0]);
			logLAll+=dist1.evaluate(pseudoTransFirst[0]);
		}
	}
	else{
		SimpleExtendedDistribution dist1 = ( (SimpleExtendedDistribution)
				 ((FreeTransitionProbs1)transProbs[0]).transitionsOut[0]);
		if(Constants.transferEquilToStart(transProbs.length)){
			DoubleMatrix2D nullspace = this.globalTrans.mat.pi;
			this.transferState(dist1,nullspace, pseudoTransFirst[0]);
		}
			else{
				
				logLAll+= dist1.evaluate(pseudoTransFirst[0]);
//					this.transferState(dist1,null, pseudoTransFirst[0]);
			}
	}
	
	
	
	
    	    		double[] pseudoCExp2  = pseudoTrans;
    	    		/*if(allFree){
    	    			pseudoCExp2 = new double[pseudoTrans.length];;
    	    			for(int i1=0; i1< pseudoCExp2.length; i1++){
    	    			pseudoCExp2[i1] = pseudoTrans[i1]/(double) Constants.expandCN('1');	
    	    			}
    	    		}*/
    	    		
    		//	 logLAll = 0;
    			  double[] start = new double[transProbs[0].noStates()];
    		         start[0] = 1.0;
    		      double[] probs = new double[start.length];
    		      fillProbs(transProbs[0], probs, start);
    			   for(int i=1; i<transProbs.length; i++){
    				   for(int update_index=0; update_index<r.length; update_index++){
    				   double d = -1*(loc.get(i)-loc.get(i-1));
    				  // final  boolean singleExpModel  = r[update_index].length==1;
       	              //  final   int len = singleExpModel ? 1: arg1[update_index].length-1;
	                   
    				   }
    		               AbstractTransitionProbs probs_ = transProbs[i];
    	              	  probs_.setHP(probs);
    		                  double d = loc.get(i)-loc.get(i-1);
    		             //   Logger.global.info("updating "+i+" of "+loc.size());
    		                 //  logLAll+=  probs_.transferAlpha(pseudoTrans, alpha_overall, allFree ? 2 : 0);
    		                  if(globalTrans==null){
    		                	if(Constants.updateAlpha())  logLAll+=  probs_.transferAlpha(pseudoTrans, alpha_overall, allFree ? 2 : 0);
    		                	  logLAll+= probs_.transfer(pseudoCExp2, alpha_overall, allFree ? 2 : 0);
    		                  }else{
    		                  
    		                   logLAll+= probs_.transferQ(pseudoCExp2, 1e10, 0, this.globalTrans.mat, allFree ? 2: 0,d,index);
    		                  }
    		           //  if(true)  {
    		           // 	   Logger.global.info(i+" RATE "+((FreeRateTransitionProbs)probs_).logrelativeRate);
    		            //   }
    		                   //probs_.transfer(pseudoCExp2, arg1, allFree ? 2 : 0);
    		                   System.arraycopy(probs,0, start,0,start.length);
    	                         fillProbs(transProbs[i], probs, start);
    		                     //  transProbs[i].logProb();
    		                }
    			   
    			   if(Constants.measureGlobal() ){
    			// rateDistribution.maximise();
    				   
    				   Logger.global.info("new rate distribution is "+this.rateDistribution.toString());
    			   }
    			   if(globalTrans!=null){
    			Logger.global.info("RATES "+this.globalTrans.mat.currentRate);
    			//Logger.global.info(getRates(transProbs));
    			   }
    			   return logLAll;
    		  
    	     
    }
 
 
 private String getRates(AbstractTransitionProbs[] transProbs2) {
	StringBuffer sb = new StringBuffer();
	for(int i=1; i<transProbs2.length; i++){
		sb.append(String.format("%5.3g ",(transProbs2[i]).getRate(1)).trim()+"\t");
	}
	return sb.toString();
}
protected double transferState(SimpleExtendedDistribution dist1,
			DoubleMatrix2D nullspace, double pseudo) {
		// TODO Auto-generated method stub
	 return  dist1.evaluate(pseudo);
	}
public static void fillProbs(AbstractTransitionProbs abstractTransitionProbs,
			double[] probs, double[] start) {
		for(int i=0; i<probs.length; i++){
			probs[i] = 0;
			for(int j=0; j<start.length; j++){
				probs[i]+=start[j]*abstractTransitionProbs.getTransition(j, i);
			}
		}
		
	}
    	         
    private int[] getOrder(double[] ds) {
    	Comparable[][] obj = new Comparable[ds.length-1][];
    	for(int i=1; i<ds.length; i++){
    		obj[i-1] = new Comparable[] {ds[i], i};
    	}
    	Arrays.sort(obj, new Comparator<Comparable[]>(){

			public int compare(Comparable[] o1, Comparable[] o2) {
				return -1*o1[0].compareTo(o2[0]);
			}

			
    		
    	});
    	int[] order = new int[obj.length];
    	for(int i=0; i<order.length; i++){
    		order[i] = (Integer)obj[i][1];
    	}
    	return order;
	}
	/*
        ( (FreeTransitionProbs1) this.transProbs[0]).transfer(pseudoTrans, pseudoCExp);
             for(int i=1; i<this.transProbs.length; i++){
                     ((AbstractTransitionProbs)this.transProbs[i]).transfer(pseudoTrans, pseudoCExp);
             }
          }*/
    
    //final Set<Integer> freeTransitions;
   /* public void transferTransitions(double pseudoTrans, double pseudoCExp){
      if(true) throw new RuntimeException("!!");
   //    Logger.global.info("transitions at 0 "+getOutString( (FreeTransitionProbs1) this.transProbs[0]));
       ( (FreeTransitionProbs1) this.transProbs[0]).transfer(pseudoTrans, pseudoCExp);
        for(int i=1; i<this.transProbs.length; i++){
         //   if(freeTransitions.contains(i)){
                
         //   }
          //  else{
                ((AbstractTransitionProbs)this.transProbs[i]).transfer(pseudoTrans, pseudoCExp);
            //}
        }
     }*/
       
      /*   this.transProbs[0].transfer(0, pseudoCExp);
        for(int i=1; i<this.transProbs.length-1; i++){
                 this.transProbs[i].transfer(pseudoTrans, pseudoCExp);
        }*/
  
     public Double[] getTrans(int from, int to){
         Double[] d = new Double[this.transProbs.length];
         d[0] = transProbs[0].getTransition(0, to);
         for(int i=1; i<d.length-1; i++){
             d[i] = this.getTransitionScore(from, to, i);
                // transProbs[i].getTransition(from, to);
         }
         d[d.length-1] = transProbs[this.transProbs.length-1].getTransition(from, 0);
         return d;
     }
     
   
     
     public void validateTransAt(int indexOfToEmission) {
         if(indexOfToEmission>=this.transProbs.length) return ;//throw new RuntimeException("!!");
         this.transProbs[indexOfToEmission].validate();
         
     }
     
   final double[] transLast;
     
     double[] makeTransLast(Integer[] cn, int numF){
    	 double[] tl = Constants.modifyFracStart();	 
    	/* if(cn[1]>0) tl = new double[numF+1];
    	 tl[1] = 1;
    	 tl[0] = 0;
    	*///double[] tl = new double[2];
    	if(!Constants.penaliseEnd()){
    		Arrays.fill(tl, 1.0);
    	}
    	return tl;
    }
     
    public boolean cached = false;
     public double getTransitionScore(int from, int to, int indexOfToEmission) {
    	
    	 
         if(indexOfToEmission>=this.transProbs.length) return 0;//throw new RuntimeException("!!");
          if(from==0 && 
              indexOfToEmission!=0) return 0;
              //return alpha[0].probs[to];
          else if(to==0){
        	  
              if(indexOfToEmission==this.transProbs.length-1 ){
            	  if(!cached && from>0 ){
            		  Integer ind = cn[from];
            		  if(ind==null || ! Constants.penaliseEnd()) return 1.0;
            		  else return transLast[ ind];
            	  }
            	  else{
            		 return  this.transProbs[indexOfToEmission].getTransition(from, to);
            	  }
              }
              else{
            	  return 0;
              }
          }
          else if(indexOfToEmission==0 && from!=0) return 0;
          else{
              AbstractTransitionProbs tp =  this.transProbs[indexOfToEmission];
              double d = tp==null ? 0 :  tp.getTransition(from, to);
            
              return  d;
          }
     }
     
public double getTransitionScoreToPaint(int from, int to, int indexOfToEmission) {
    	if(true) return this.getTransitionScore(from, to, indexOfToEmission);
    	 
         if(indexOfToEmission>=this.transProbs.length) return 0;//throw new RuntimeException("!!");
          if(from==0 && 
              indexOfToEmission!=0) return 0;
              //return alpha[0].probs[to];
          else if(to==0){
        	  
              if(indexOfToEmission==this.transProbs.length-1 ){
            	  if(!cached && from>0 ){
            		  return transLast[cn[from]];
            	  }
            	  else{
            		 return  this.transProbs[indexOfToEmission].getTransition(from, to);
            	  }
              }
              else{
            	  return 0;
              }
          }
          else if(indexOfToEmission==0 && from!=0) return 0;
          else{
              AbstractTransitionProbs tp =  this.transProbs[indexOfToEmission];
              double d = tp==null ? 0 :  tp.getTransitionToPaint(from, to);
            
              return  d;
          }
     }
     
    
     public void validate(){
         for(int i=0; i<this.transProbs.length-1; i++){
        	 if(i>0) validate(transProbs[i], transProbs[i].noStates(), i);
             transProbs[i].validate();
         }
     }
  
    public abstract void initialise(double[] dist, double permute, double u_glob) throws Exception;
    
     public void initialise( List<State>states, Double[] exp_p, double u_g){
      this.initialise(states,  null, null, exp_p,  null, null, u_g);   
     }

     
     public void initialise( List<State>states, List<Integer> loc, Double[] r, Double[] exp_p1, 
              double[] rel, double[] relst, double u_g){
         try{
             
          //   this.trans_init = new double[states.size()];
             boolean hasStart = states.get(0) instanceof DotState;
             int numFounders = states.size()-1;
          //   Dirichlet  dir1 =  null;
        //     this.MAGIC = states.get(0);
             double[] d2 = new double[numFounders+1];
             double[] d3 = new double[numFounders+1];
             if(rel!=null){
                 for(int i=0; i<rel.length;i ++){
                     d2[i+1] = rel[i];
                     d3[i+1] = relst[i];
                 }
//                 System.arraycopy(rel, 0, d2, 1, rel.length);
             }
             else
             {
                 double inv_col = 1.0/(double)numFounders;
                 Arrays.fill(d2, inv_col);
             }
             d2[0] = 0.0;
             d3[0] = 0.0;
             //Constants.normalise(d3);
            // dir1= new Dirichlet(d2, Constants.u_global(0)[1]);
           
             transProbs[0] =new FreeTransitionProbs1(true, new Dirichlet(d3, 1e100), states.size());
             if(Constants.CHECK ){
                 transProbs[0].validate();
             }
             initialise(  d2,0, u_g);// Constants.samplePermute());
             int[] transMode0 = Constants.transMode(0);
         
             if(transMode0.length==1){
            	    int numF = numFounders;
           	  this.alpha_overall =new double[(transMode0[0]==2 || transMode0[0]==4) ?numF+1:1][numF+1];
           	for(int i=0; i<alpha_overall.length; i++){
          	  alpha_overall[i] = d2;
            }
            }
            else{
                int numF = ((BetweenWithinTransitionProbs1)this.transProbs[1]).transBetweenGroups.noStates();
            this.alpha_overall =new double[(transMode0[0]==2 || transMode0[0]==4) ?numF:1][numF+1];
            for(int i=0; i<alpha_overall.length; i++){
          	  alpha_overall[i] = this.transProbs[1].getAlphaPrior().clone();
            }
            }
             if(transProbs.length>1 && this.transProbs[1] instanceof ExponentialTransitionProbs){
            	ExpTransProb alpha =  ((ExponentialTransitionProbs)transProbs[1]).alpha;
            	if(alpha instanceof MultiExpProbs){
	             for(int k=0; k<alpha_overall.length; k++){
	                // for(int i=0; i<alpha_overall[k].length; i++){
	               	  alpha_overall[k] =((MultiExpProbs)alpha).getProbs(k).clone();
	                 //}
	              }
	           //  if(alpha_overall[0]==alpha_overall[1]) throw new RuntimeException("!!");
            	}
             }
          //   this.allFree = transMode0[0]==1;  
             // transProbs[i] =new FreeTransitionProbs(states, false, dir);
           //  transProbs[transProbs.length-1] =new FreeTransitionProbs( false, dir1, states.size());
         }catch(Exception exc){exc.printStackTrace();}
     }
  




   protected void validate(AbstractTransitionProbs probs, int nostates, int index) {
       for(int i=1; i<nostates; i++){
           double sum=0;
           double[] dist = new double[nostates];
           Arrays.fill(dist,0.0);
           for(int k=1; k<nostates; k++){
               dist[k] = probs.getTransition(i, k);
               sum+= dist[k];
                 
           }
           if(Math.abs(sum-1.0)>SimpleDistribution.tolerance) {
               probs.validate();
               throw new RuntimeException("!! "+sum+" "+index+" "+i+" "+probs.getTransition(i, i));
           }
        //   else{
          //     Logger.global.info("right");
         //  }
           
       }
        
    }



    



    public void print(final PrintWriter sb, final String sbS, final List<State> states,
            double[][]hittingProb) {
    //	sb.println("r prior is "+Arrays.asList(r));
    	
        Double[] top = new Double[hittingProb[0].length];
        for(int i=0; i<this.transProbs.length; i++){
            sb.print(i+"  ");
            if(transProbs[i]==null) sb.print("null");
            else{
                sb.println(transProbs[i].getClass());
                if(i!=0){
                for(int k=0; k<top.length; k++){
                    top[k] = hittingProb[i-1][k];
                }
                }
                    transProbs[i].print(sb, i==0 ? null : top, i==0  || loc==null? 1 : loc.get(i)-loc.get(i-1));
               
            }
        }
    }
	public void printR(File fileout) {
		try{
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(fileout, "transitionModel.txt"))));
		//pw.print("rates: ");
		//pw.println("");
		//rateDistribution.print(pw);
		// TODO Auto-generated method stub
		//if(globalTrans!=null)this.globalTrans.mat.print(pw);
		//pw.println();
		//pw.println("relative rates");
		for(int i=1; i<this.transProbs.length; i++){
			if(transProbs[i] instanceof FreeRateTransitionProbs){
				
				pw.println(i+"\t"+
						String.format("%5.2f",((FreeRateTransitionProbs)(transProbs[i])).logrelativeRate));
			}
			pw.print(DataCollection.datC.snpid.get(i));
			pw.print("\t");
			transProbs[i].print(pw, null, (double) this.loc.get(i) - (double) this.loc.get(i-1));
		}
		for(int ii=0; ii<r.length; ii++){
			pw.println("r priors");
			for(int i=0; i<r[ii].length; i++){
				pw.println(i+":" +r[ii][i]);
			}
		pw.println("r prior is for "+ii+" : " +r[ii][0]/Constants.probCrossOverBetweenBP);//+" "+r[ii][1]/Constants.probCrossOverBetweenBP);
		pw.print("--expModelIntHotSpot1");
		for(int k=0; k<r[ii].length; k++){
			//pw.println(r[ii][k]);
			pw.print((k==0 ? "\t":":")+String.format("%7.4f",r[ii][k]/r[ii][0]).replaceAll(" ",""));
		}
		pw.println();
		pw.print("--modifyFrac0");
		if(alpha_overall!=null){
		for(int ind =0; ind<alpha_overall.length; ind++){
			pw.print((ind==0 ? "\t" : ":")+alpha_overall[ind][1]);
			for(int i=2; i<alpha_overall[ind].length; i++){
				pw.print(";"+alpha_overall[ind][i]);
			}
		
		}
		
		pw.print("\n--modifyFrac0");
		for(int ind =0; ind<alpha_overall.length; ind++){
			pw.print("\n");//+alpha_overall[ind][1]);
			for(int i=1; i<alpha_overall[ind].length; i++){
				pw.print(";"+String.format("%5.3g",alpha_overall[ind][i]));
			}
		
		}
		}
		pw.println();
		}
		pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public  double updateRates(double[][] mat, double[] pi) {
		// TODO Auto-generated method stub
		this.globalTrans.updateRates(mat, pi);
		return this.globalTrans.mat.currentRate;
	}
	
	 
	
	
  
    
    
    
 
    
   
 }