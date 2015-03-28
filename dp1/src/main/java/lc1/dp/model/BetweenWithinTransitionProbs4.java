package lc1.dp.model;

import java.util.Arrays;
import java.util.logging.Logger;

import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.BetweenWithinTransitionProbs3;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.util.Constants;

public class BetweenWithinTransitionProbs4 extends FreeTransitionProbs1 {

	private final AbstractTransitionProbs transBetweenGroups;
	private final int[] stateToGroup;
	 private final int[][] groupToState;
	private int[] stateToIndexWithinGroup;
	

	


	public BetweenWithinTransitionProbs4(
			BetweenWithinTransitionProbs3 betwith,
			AbstractTransitionProbs betweenGroups) {
		super(betwith);
		this.transBetweenGroups = betweenGroups;
		this.stateToGroup = betwith.stateToGroup;
		this.groupToState = betwith.groupToState;
		this.stateToIndexWithinGroup= betwith.stateToIndexWithinGroup;
		this.pseudo = new double[groupToState.length];
	}

	
	
	  
	  @Override
	    public void addCount(int from, int to, double val){
	       super.addCount(from, to, val);
	        this.transBetweenGroups.addCount(this.stateToGroup[from], this.stateToGroup[to], val);
	    }
	  
	  /* (non-Javadoc)
	     * @see lc1.dp.AbstractTransitionProbs#initialiseCounts(boolean, boolean)
	     */
	     public void initialiseCounts(boolean start, boolean end){
	         super.initialiseCounts(start, end);
	         this.transBetweenGroups.initialiseCounts(start, end);
	     }
	     
	     
	     @Override
	 	public double transfer(double[] pseudoCExp, double[][] d, int pos_index) {
	    	 if(pos_index==2){
	    		return  super.transfer( pseudoCExp,d==null ? null : expand1(d),pos_index);
	    	 }
	    	 else{
	    		 return this.transBetweenGroups.transfer(pseudoCExp, d, pos_index);
	    	 }
	 	}
	     @Override
	     public double  transferAlpha(double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
	    	 if(pos_index==2){
	    		return  super.transferAlpha( pseudoTrans, expand(alpha_overall), pos_index, this.groupToState);
	    	 }
	    	 else{
	    		return this.transBetweenGroups.transferAlpha(pseudoTrans, alpha_overall, pos_index);
	    	 }
	     }
	   public static final boolean print = true;
	   final double[] pseudo;
	   
	   @Override
	     public void harmonise(int[] gToS, double[] d, double[] pseudo1, int index, int group){
	 		//  System.err.println("harmonising");
		   getTrans(group);
	 		  Arrays.fill(d,0);
	 			 double total =0;
	 			 for(int i1=0; i1<gToS.length; i1++){
	 				 int state = gToS[i1];
	 				 double cnt = Constants.sum(transitionsOut[state].counts());
	 				 double[] prob = transitionsOut[state].probs();
	 			if(print)	System.err.println("before "+i1+"\t\t\t\t"+Constants.print(shorten(prob))+"\t"+cnt);
	 				 for(int j=0; j<prob.length; j++){
	 					 d[stateToGroup[j]]+=cnt*prob[j];
	 				 }
	 				 total+=cnt;
	 			 }
	 			 if(total==0) return;
	 			 for(int j=0; j<d.length; j++){
	 				 d[j] = d[j]/total;
	 			 }
	 		if(print){
	 			System.err.println(index +" target "+Constants.print(shorten(pseudo))+"\n actual "+Constants.print(d));
	 			 Logger.global.info("h");
	 		}
	 			
	 			 for(int i1=0; i1<gToS.length; i1++){
	 				 int state = gToS[i1];
	 				 double[] prob = transitionsOut[state].probs();
	 				  double sum = 0;
	 				  for(int j=0; j<prob.length; j++){
	 					
	 					 if(d[stateToGroup[j]]!=0)  prob[j] = prob[j] *(pseudo[stateToGroup[j]]/d[stateToGroup[j]]);
	 					  sum+=prob[j];
	 				  }
	 				Constants.normalise(prob);
	 				// System.err.println("after "+i1+"\t\t\t\t"+Constants.print(prob));
	 			 }
	 	  }
	     
	     
	     private void getTrans(int group) {
		for(int i=0; i<this.pseudo.length; i++){
			pseudo[i] = this.transBetweenGroups.getTransition(group, i);
		}
		
	}




		/*  @Override
		   public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
		   	//flip?
		   	//flip(pseudoTrans);flip( pseudoCExp);//flip(d);
		   	double logL= this.transBetweenGroups.transferAlpha( pseudoTrans, expand(alpha_overall), pos_index, this.groupToState);
		  return logL;
				
				  @Override
		public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
			double logL = 
			}
			}*/
	  
	     private double[] shorten(double[] pseudo) {
			double[] res = new double[this.groupToState.length];
			Arrays.fill(res,0.0);
			for(int i=0; i<pseudo.length; i++){
				res[this.stateToGroup[i]]+=pseudo[i];
			}
			return res;
		}




		public double[] expand(double[] d, int[] stateToGroup){
	  	   double[] results = new double[stateToGroup.length];
	  	   for(int i=0; i<stateToGroup.length; i++){
	  		   results[i] = d[stateToGroup[i]];
	  	   }
	  	   return results;
	     }
	    
	     public double[][] expand(double[][] d){
	  	   if(d==null) return null;
	  	   double[][] res = new double[stateToGroup.length][];
	  	   for(int i=0; i<res.length; i++){
	  		   res[i] = d[stateToGroup[i]];
	  	   }
	  	   return res;
	     }
	     public double[][] expand1(double[][] d){
	  	   double[][] res = new double[d.length][];
	  	   res[0] = expand(d[0], stateToGroup);
	  	   for(int i=1; i<res.length; i++){
	  		   res[i] = d[i];
	  	   }
	  	   return res;
	     }
	
	

}
