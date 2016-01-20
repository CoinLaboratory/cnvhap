/**
 * 
 */
package lc1.dp.model;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Logger;

//import lc1.dp.data.classification.ROC;
import lc1.dp.states.DotState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.BetweenWithinTransitionProbs1;
import lc1.dp.transition.ExpTransProb;
import lc1.dp.transition.ExponentialTransitionProbs;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.dp.transition.MultiExpProbs;
import lc1.stats.Dirichlet;
import lc1.stats.NormalDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;
import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;
import cern.jet.random.Gamma;

public abstract class CopyOfSiteTransitions implements Serializable{
     //State MAGIC;
   // Set<Integer> special;
    final int len;
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
    protected Gamma[][] dist; // this is prior distributions on r
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
    final int[] cn;
      public CopyOfSiteTransitions(int len, int numStates, int[] cn) throws Exception{
         this.len = len;
         this.cn = cn;
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
         this.dist = null;
         alpha_overall = new double[1][numStates];
     }
     //transProbs[0] is from magic state, and transProbs[noSnps] is transition to magic state
     public static int[] getCN(List<State> states){
    	 int[] cn = new int[states.size()];
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
      
      public CopyOfSiteTransitions(List<Integer>loc, int numF ,Double[] exp_p1,Double[] r,  int length, int index, int[] cn){
        if(r!=null && exp_p1!=null) throw new RuntimeException("should not both be non null");
         this.transProbs = new AbstractTransitionProbs[length];
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
         this.dist = new Gamma[transMode0.length==1 ? 1 : 2][];
         for(int i=0; i<this.r.length; i++){
        	 this.r[i] = new Double[(transMode0[i]==3 || transMode0[i]==4) ? Constants.expModelIntHotSpot1(i).length:1];
        	 this.dist[i] = new Gamma[this.r[i].length];
        	 this.r[i][0] = r[i]; // out of start?
        	 this.dist[i][0] = makeDist(this.r[i][0],index);
        	 if(r!=null && this.r[i].length>1){
        		 double[] hotspot1 = Constants.expModelIntHotSpot1(i);
     	     	for(int k=0; k<this.r[i].length; k++){
     	     		
     	     		//for(int j=0; j<this.r[i][k].length; j++){
     	     			this.r[i][k] = hotspot1[k]*r[i];
     	     			this.dist[i][k] = makeDist(this.r[i][k],index);
     	     		//}
     	     	}
          	
            }
         }	
        		// if((transMode0[i]==3 || transMode0[i]==4)){
                	// if(Constants.expModelIntHotSpot1(i).length!=(Constants.modify0.length+1)) throw new RuntimeException("mismatch between expModelIntHotSpot1 and modify0");
                // }
         }
        
       
        		
     
     	
        // this.r_prev = new Doub
     //   this.transProbType =transProbType;
     }
     private Gamma makeDist(double mode, int i) {
		//double var = Math.pow(Constants.expSd(i),2);
	
	//	double alpha =1.0/Constants.expSd(i);// (mean*mean) /var;
		double lambda =1.0/(1-Constants.expSd(i));
		double alpha = 1.0/lambda;
		//alpha = mode /lambda +1;
		
		System.err.println("shape scale "+alpha+" "+lambda);
		
		Gamma dist =  new Gamma(1.0, 100, null);
		if(false){
double[] d = new double[1000];
double[] x = new double[1000];
		double len = d.length;
		for(int i1=0; i1<len; i1++){
			
			
			x[i1] = Math.pow(10,10* (((double) i1)/((double)len)) - 10.0);
			d[i1] = dist.pdf(x[i1]);
		}
//		ROC.plot(x,d,true, false);
	//	Logger.global.info("h");
		}
		return dist;
	}
	
	public boolean converged(){
         return true;
     }
     
 
   
     public abstract CopyOfSiteTransitions clone(boolean swtch);
     
     /** if pseudocount=0 we do not use this as pseudo */
     public CopyOfSiteTransitions(CopyOfSiteTransitions trans_init, boolean swtch){
         this(trans_init.loc,trans_init.numFounders, trans_init.exp_p1, getFirst(trans_init.r),  trans_init.transProbs.length,0, trans_init.cn);
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
         for(int i=0; i< this.transProbs.length; i++){
             transProbs[i].initialiseCounts(i==0, i==transProbs.length-1);
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
 
 public double transferTransitions(final double[] pseudoTrans,final double[] pseudoTransFirst,  final double[] pseudoCExp, final int index){
	  double  logLAll =0;
	//  pseudoTransFirst[0] = pseudoTrans[0];
	 final  double[] pseudoCExp1 = pseudoCExp;//new double[pseudoCExp.length];
	
//	  Arrays.fill(pseudoCExp1, 1e10);
	 // pseudoTransFirst[1] = 1000;
	  //NOTE WE ARE NOT CURRENTLY TRAINNIG START!!!
	  transProbs[0].transferAlpha(pseudoTransFirst, null, 0); 
	   final    double[][] arg1 = new double[r.length][];
		for(int i=0; i<arg1.length; i++){
			arg1[i] = new double[r[i].length];
		}
	 if(!allFree || Constants.newTrans()){
	   final  int update_index =0;
	  
	   if(hierarchy){
	   if(true){
    	
		 
    	//  r_prev_prev = r_prev;
    	       //  r_prev = r[0];
    	    //    ( (FreeTransitionProbs1) this.transProbs[0]).transfer(pseudoTrans, pseudoCExp);
    	      
    	     if(pseudoCExp[0]<1e-3 ){
    	    	 for(int i=1; i<transProbs.length; i++){
    	    		 transProbs[i].transfer(pseudoCExp, null, 0);
    	    	 }
    	     }
    	     else{
    	        
    	        	final  boolean singleExpModel  = r[update_index].length==1;
    	            final   int len = singleExpModel ? 1: arg1[update_index].length-1;
    
    	
    	 System.err.println(" before r "+Arrays.asList(r[0]));
    
    	final double[] upperBound =new double[len]; 
    	
    	final  double[] lowerBound = new double[len];
    	final Gamma[] dist1 = new Gamma[len];
    	   //  final double[] arg2 = new double[this.numFounders];
    	      //  final double[][] pseudoC = new double[2][2];
    	        MultivariateFunction uvf = new MultivariateFunction(){
    	        	double[] exp = new double[len];
    	            public double evaluate(double[] arg) {
    	                double logProb=0;
    	                for(int k=0; k<len; k++){
    	                	exp[k] = Math.exp(arg[k]);
    	                }
    	                for(int i=1; i<transProbs.length; i++){
    	                	
    	                    double d = -1*(loc.get(i)-loc.get(i-1));
    	                    if(singleExpModel){
    	                    	arg1[update_index][0] = exp[0]*d;
    	                    }
    	                    else{
	    	                    arg1[update_index][0] =  r[update_index][0]*d;
	    	                    for(int k=0; k<len;k++){
	    	                    	arg1[update_index][k+1] = exp[k]*d;
	    	                    }
    	                    }
    	                    logProb+= transProbs[i].transfer(pseudoCExp1, arg1,update_index);
    	                     //  transProbs[i].logProb();
    	                }
    	       
    	           
    	           
    	          
    	                double prior =0;
    	                for(int i=0; i<arg.length; i++){
    	                	prior += Math.log(dist1[i].pdf(exp[i]));
    	                		//normal.logpdf(arg[i], init[i], Constants.expSd());
    	                }
    	          //   System.err.println("evaluate "+Math.exp(arg[0])+
    	    	  //		 " "+Math.exp(arg[1])+" "+Math.exp(arg[2])+
    	    	    //		  " "+logProb+" "+prior+" "+(logProb+prior));
    	                double res = -logProb;
    	                 res=-prior;
    	                return res;
    	             //  return -1*(logProb+prior);
    	            }
    	           
					public double getLowerBound(int n) {
						// TODO Auto-generated method stub
						return lowerBound[n];
					}

					public int getNumArguments() {
						// TODO Auto-generated method stub
						return len;
					}

					public OrthogonalHints getOrthogonalHints() {
						// TODO Auto-generated method stub
						return null;
					}

					public double getUpperBound(int n) {
						// TODO Auto-generated method stub
						return upperBound[n];
					}
    	        };
    	        double[] xvec = new double[len];
    	        if(singleExpModel){
    	        	xvec[0] = Math.log(r[update_index][0]);
    	        	lowerBound[0] = Math.log(r[update_index][0]*1e-6);
    	        	upperBound[0] = Math.log(r[update_index][0]*1e6);
    	        	dist1[0] =this.dist[update_index][0];
    	        }
    	        else{
	    	        for(int k=0;k<len; k++){
	    	        	xvec[k] = Math.log(r[update_index][k+1]);
	    	        	lowerBound[k] = Math.log(r[update_index][k+1]*1e-6);
	    	        	upperBound[k] = Math.log(r[update_index][k+1]*1e6);
	    	        	dist1[k] =this.dist[update_index][k+1];
	    	        }
    	        }
    	   
    	      //  uvf.evaluate(xvec);
    	      //  System.exit(0);
    	        MultivariateMinimum uvm =new OrthogonalSearch();// new pal.math.ConjugateDirectionSearch();
    	     //   double[] x_init = new double[xvec.length];
    	    //    System.arraycopy(xvec, 0, x_init, 0, xvec.length);
    	       if(hierarchy) uvm.optimize(uvf, xvec, 0.01,0.01);
    	       // else uvf.evaluate(xvec);
    	        logLAll-=uvf.evaluate(xvec);//uvf.evaluate(xvec);
    	        if(Constants.CHECK && Double.isInfinite(logLAll)){
    	        	throw new RuntimeException("is infinite");
    	        }
    	        if(singleExpModel){
    	        	r[update_index][0] = Math.exp(xvec[0]);
    	        }
    	        else{
	    	        for(int k=0;k<len; k++){
	    	        	r[update_index][k+1] = Math.exp(xvec[k]) ;
	    	        }
    	        }
    	       System.err.println("r update_index "+Constants.print(r[update_index]));
    	        for(int i=1; i<transProbs.length; i++){
    	            double d = -1*(loc.get(i)-loc.get(i-1));
    	            for(int k=0; k<r[update_index].length;k++){
                    	arg1[update_index][k] = r[update_index][k]*d;
                    }
    	                // pseudoC[0] = pseudoCExp*Math.exp(d*r[0]);
    	                // pseudoC[1] = pseudoCExp -pseudoC[0];
    	                 transProbs[i].transfer(pseudoCExp, arg1,update_index);
    	             //   ((ExponentialTransitionProbs)transProbs[i]).transfer(pseudoTrans, pseudoC);
    	        }
    	        }
	   }
    	     if(pseudoTrans[0]<1e-3 ){
    	    	 for(int i=1; i<transProbs.length; i++){
    	    		 transProbs[i].transferAlpha(pseudoTrans, null, 0);
    	    	 }
    	     }
    	     else{
    	    //  final   double[] arg2 = new double[alpha_overall.length-2];
    	        //arg2[0] = 0;
    	        //now do alpha
    	     // final double[][] alpha_overall1 = new double[alpha_overall.length][alpha_overall[0].length];//new double[alpha_overall.length][alpha_overall[0].length];
    	      final	int[] order = getOrder(((double[][])alpha_overall)[0]);// {1,2,3,4};//new int[] {2,3,1,4};
     	    int start = alpha_overall.length==1  ? 0 : 1; //first one is non-informative if more than one
    	       for(int ind1=start; ind1<
    	       alpha_overall.length; 
    	       ind1++) {
    	    	   final int ind=ind1;
    	       	 final   int len = alpha_overall[ind].length-2;
    	      // 	 alpha_overall[0] = 0;
    	       final	 double[] xvec = new double[len];
    	      
    	       if(Math.abs(1.0 - Constants.sum(alpha_overall[ind]))>0.00001){
    	    	   throw new RuntimeException("!!");
    	       }
    	       int max_ind = Constants.getMax(alpha_overall[ind]);
    	       if(alpha_overall[ind][max_ind]>1.0 - 1e-10) {
    	    	   continue;
    	       }
    	       //	 System.err.println(" before r "+Arrays.asList(alpha_overall));
    	       	   //  final double[] arg2 = new double[this.numFounders];
    	       	      //  final double[][] pseudoC = new double[2][2];
    	       	        MultivariateFunction uvf = new MultivariateFunction(){
    	       	        	public boolean transfer(double[] arg){
    	       	        	  double rem = 1.0;
	       	                    for(int k=0; k<len;k++){
	       	                    	if(arg[k]>1.0 || arg[k]<0.0){
	       	                    		return false;
	       	                    	}
	       	                    	if(Constants.CHECK && Double.isNaN(arg[k])) {
	       	    	        			throw new RuntimeException("!!");
	       	    	        		}
	       	                    	alpha_overall[ind][order[k]] = arg[k]*(rem);
	       	                    	rem-=alpha_overall[ind][order[k]];
	       	                    	if(rem<0){
	       	                    		throw new RuntimeException("!!");
	       	                    	}
	       	                    }
	       	             	if(rem<0){
 	                    		throw new RuntimeException("!!");
 	                    	}
	       	                    alpha_overall[ind][order[len]] = rem;
	       	                    //Constants.normalise(alpha_overall[ind]);
		       	              /*  for(int i=0; i<alpha_overall1.length; i++){
		       	    	        	
		       	    	        	//for(int j=0; j<alpha_overall1[i].length; j++){
		       	    	        		for(int k=0; k<alpha_overall1[i].length; k++){
		       	    	        		alpha_overall1[i][k] = alpha_overall[i][k]*pseudoTrans[update_index]; 
		       	    	        		}
		       	    	        		
		       	    	        	//}
		       	    	       }*/
		       	                return true;
    	       	        	}
    	       	            public double evaluate(double[] arg) {
    	       	                double logProb=0;
    	       	                transfer(arg);
    	       	           //  logProb+= transProbs[0].transferAlpha(pseudoTrans, alpha_overall1, 0);
    	       	                for(int i=1; i<transProbs.length; i++){
    	       	            
    	       	               
    	       	                	
    	       	                 //   double d = -1*(loc.get(i)-loc.get(i-1));
    	       	                	
    	       	                    logProb+= transProbs[i].transferAlpha(pseudoTrans, alpha_overall, 0);
    	       	                   
    	       	                     //  transProbs[i].logProb();
    	       	                }
    	       	        //   System.err.println("evaluate alpha "+arg[0]+" "+arg[1]+" "+logProb);
    	       	          /*   if(Double.isNaN(logProb)){
	       	                    	throw new RuntimeException("!!");
	       	                    }*/
    	       	               return -1*(logProb);
    	       	            }

    	       	           
    	   					public double getLowerBound(int n) {
    	   						//if((ind==1 || ind==3) && n==0) return getUpperBound(n);
    	   						// TODO Auto-generated method stub
    	   						return 1e-5;//Math.max(0.0, xvec[n] - 0.05);
    	   					}

    	   					public int getNumArguments() {
    	   						// TODO Auto-generated method stub
    	   						return len;
    	   					}

    	   					public OrthogonalHints getOrthogonalHints() {
    	   						// TODO Auto-generated method stub
    	   						return null;
    	   					}

    	   					public double getUpperBound(int n) {
    	   						// TODO Auto-generated method stub
    	   					//	if(ind==2)
    	   						return 1.0-1e-5;//Math.min(1.0,xvec[n]+0.05);
    	   					}
    	       	        };
    	       	       
    	       	        double rem = 1.0;
    	       	        for(int k=0;k<len; k++){
    	       	        	xvec[k] =  alpha_overall[ind][order[k]]/ rem;
    	       	        	rem-=alpha_overall[ind][order[k]];
    	       	        }
    	       	        
    	       	      //  uvf.evaluate(xvec);
    	       	      //  System.exit(0);
    	       	        MultivariateMinimum uvm = new OrthogonalSearch();
    	       	      if(hierarchy)  uvm.optimize(uvf, xvec, 0.01,0.0001);
    	       	        rem  = 1.0;
    	       	  double log1 = -  uvf.evaluate(xvec);
    	       	        System.err.println("alpha "+ind+Constants.print(alpha_overall[ind])+" "+log1);
    	       	     logLAll-= log1;
    	       	    if(Constants.CHECK && Double.isInfinite(logLAll)){
        	        	throw new RuntimeException("is infinite");
        	        }
    	       	        }
    	     }
	 }
	 }
    	     if((allFree || ! hierarchy )){
    	    		double[] pseudoCExp2 = pseudoCExp;
    	    		if(allFree){
    	    			for(int i1=0; i1< pseudoCExp2.length; i1++){
    	    			pseudoCExp2[i1] = pseudoCExp[i1]/(double) Constants.expandCN("1");	
    	    			}
    	    		}
    		//	 logLAll = 0;
    			  double[] start = new double[transProbs[0].noStates()];
    		         start[0] = 1.0;
    		      double[] probs = new double[start.length];
    		      fillProbs(transProbs[0], probs, start);
    			   for(int i=1; i<transProbs.length; i++){
    				   for(int update_index=0; update_index<r.length; update_index++){
    				   double d = -1*(loc.get(i)-loc.get(i-1));
    				   final  boolean singleExpModel  = r[update_index].length==1;
       	                final   int len = singleExpModel ? 1: arg1[update_index].length-1;
	                    if(singleExpModel){
	                    	throw new RuntimeException("!!");
	                    	
	                    }
	                    else{
	   	                    arg1[update_index][0] =  r[update_index][0]*d;
	   	                    for(int k=0; k<len;k++){
	   	                    	arg1[update_index][k+1] = r[update_index][k+1]*d;
	   	                    }
	                    }
    				   }
    		               AbstractTransitionProbs probs_ = transProbs[i];
    	              	  probs_.setHP(probs);
    		                 //   double d = -1*(loc.get(i)-loc.get(i-1));
    		                	
    		                   logLAll+=  probs_.transferAlpha(pseudoTrans, alpha_overall, allFree ? 2 : 0);
    		                   logLAll+= probs_.transfer(pseudoCExp2, arg1, allFree ? 2 : 0);
    		                   System.arraycopy(probs,0, start,0,start.length);
    	                         fillProbs(transProbs[i], probs, start);
    		                     //  transProbs[i].logProb();
    		                }
    			   return logLAll;
    		   }
    	     System.err.println("transition logL "+logLAll);
    	     return logLAll;
    	     //   System.err.println(" setting r "+Arrays.asList(r));
    	      // System.exit(0);
    	     
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
     
   
     public static double[] transLast = Constants.modifyFracStart();
    static{
    	if(!Constants.penaliseEnd()){
    		Arrays.fill(transLast, 1.0);
    	}
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
              double d = tp==null ? 0 :  tp.getTransition(from, to);
            
              return  d;
          }
     }
     
    
     public void validate(){
         for(int i=0; i<this.transProbs.length-1; i++){
        	 if(i>0) validate(transProbs[i], transProbs[i].noStates(), i);
             transProbs[i].validate();
         }
     }
  
    public abstract void initialise(double[] dist, double permute, int level) throws Exception;
    
     public void initialise( List<State>states, Double[] exp_p, int level){
      this.initialise(states,  null, null, exp_p,  null, null, level);   
     }

     
     public void initialise( List<State>states, List<Integer> loc, Double[] r, Double[] exp_p1, 
              double[] rel, double[] relst, int level){
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
           
             transProbs[0] =new FreeTransitionProbs1(true, new Dirichlet(d3, 1e10), states.size());
             if(Constants.CHECK ){
                 transProbs[0].validate();
             }
             initialise(  d2,0, level);// Constants.samplePermute());
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
             this.allFree = transMode0[0]==1;  
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
	public void printR(PrintWriter pw) {
		// TODO Auto-generated method stub
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
	}
  
    
    
    
 
    
   
 }