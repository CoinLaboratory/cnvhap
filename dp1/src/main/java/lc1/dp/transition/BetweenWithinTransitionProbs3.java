package lc1.dp.transition;

import lc1.stats.Dirichlet;
import lc1.stats.Sampler;
import lc1.util.Constants;


public class BetweenWithinTransitionProbs3 extends BetweenWithinTransitionProbs1 {
    public final int[][] groupToState;


	@Override
    public double getBetweeenGroupTrans(int from, int groupFrom, int groupTo){
        return  transBetweenGroups.getTransition(from,
                groupFrom
                ,  groupTo);
     }
   @Override
    public void addGroupTransCount(int from, int groupFrom, int groupTo, double val, double dist){
        this.transBetweenGroups.addCount(
                from,
                groupFrom
                , groupTo, val);
    }
   @Override
   public  AbstractTransitionProbs makeTransBetweenGroups(Object clazz, Sampler samplerFirst, Sampler exp_p, int no_states, int[] stateToGroup) throws Exception{
       AbstractTransitionProbs res =   ExponentialTransitionProbs.get(clazz, samplerFirst, exp_p,
              //samplerFirst.dist.length
                 no_states, expand(Constants.expModelIntHotSpot1(0), stateToGroup)
               );
       if(Constants.CHECK){
           res.validate();
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
   
   @Override
	public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
		if(pos_index==2){
		double logL = this.transBetweenGroups.transfer( pseudoExp,expand1(d),0);
		
	        for(int i=0; i<transWithinGroups.length; i++){
	        	//should be pos_index+1
	                logL+= this.transWithinGroups[i][i].transfer( pseudoExp, null, 1);
	        }
	    	return logL;
		}
		else{
			 return this.transBetweenGroups.transfer( pseudoExp,expand1(d),pos_index);
		}
	}
   
   
   
  
   
   public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
   	if(pos_index==2){
   		double logL = this.transBetweenGroups.transferAlpha( pseudoTrans, expand(alpha_overall), 0);
   		 for(int i=0; i<transWithinGroups.length; i++){
   			 for(int i1=0; i1<transWithinGroups.length; i1++){
   			 ///note should be pos_index + 1, but this would lead to 3rd index
                logL+=this.transWithinGroups[i][i1].transferAlpha(pseudoTrans, null, i==i1 ? 1 : 2);//, alpha_overall, pos_index+1);
   			 }
   		 }
   		//logL+=this.optimiseBack(false);
   		 return logL;
   	}
   	else{
   		return this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, pos_index);
   	}
		
	}
   
   @Override
   public double  transferQ(double[] ds, double pseudoAlpha, double pseudoRate, MatrixExp initial, int pos_index, double distance, int it) {
		if(pos_index==2){
			double logL = this.transBetweenGroups.transferQ( ds,pseudoAlpha, pseudoRate, initial,0,distance,it);
			
		        for(int i=0; i<transWithinGroups.length; i++){
		        	//should be pos_index+1
		                logL+= this.transWithinGroups[i][i].transfer( ds, null, 1);
		                for(int i1=0; i1<transWithinGroups.length; i1++){
		          			 ///note should be pos_index + 1, but this would lead to 3rd index
		                       logL+=this.transWithinGroups[i][i1].transferAlpha(ds, null, i==i1 ? 1 : 2);//, alpha_overall, pos_index+1);
		          	     }
		        }
		    	return logL;
			}
			else{
				 return this.transBetweenGroups.transferQ( ds,pseudoAlpha, pseudoRate, initial,0,distance,it);
			}
		
	}
   public double[][] expand1(double[][] d){
	   double[][] res = new double[d.length][];
	   res[0] = expand(d[0], stateToGroup);
	   for(int i=1; i<res.length; i++){
		   res[i] = d[i];
	   }
	   return res;
   }
   
   
   
   
   @Override
   public AbstractTransitionProbs cloneBetweenGroupsTrans(AbstractTransitionProbs transBetweenGroups
		   ,int[] stateToGroup,
		   //int[][] statesToGroupTrans, 
		   double[] u){
   	AbstractTransitionProbs probs =  transBetweenGroups.clone(stateToGroup, u);
   //	((FreeTransitionProbs1)probs).modify(statesToGroupTrans, stateToGroup);
 //modify((FreeTransitionProbs1) probs, stateToGroup);
   	return probs;
   }
    
  
@Override
   public AbstractTransitionProbs clone(boolean swtch){
       
       if(swtch) return new FreeTransitionProbs1(this);
       else return new BetweenWithinTransitionProbs3(this, swtch);
   }
   
   @Override
   public AbstractTransitionProbs makeTransBetweenGroups(AbstractTransitionProbs transBetweenGroups, boolean swtch, int[] statesToGroup){
       return transBetweenGroups.clone(statesToGroup, null);
   }
   
   public BetweenWithinTransitionProbs3(BetweenWithinTransitionProbs1 probs1) {
       super(probs1, false);
       this.groupToState = null;
   }

  
    
    public BetweenWithinTransitionProbs3( 
          
    //        AbstractTransitionProbs transBetweenGroups,
            Sampler samplerFirst, 
            Sampler[] samplers,
            int[] stateToGroup,
            int[] stateToIndexWithinGroup,
            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState) throws Exception{
       super(samplerFirst, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, clazz, stateToGroup.length, groupToState);
       this.groupToState = groupToState;
    }
    
    public BetweenWithinTransitionProbs3( 
            
    	    //        AbstractTransitionProbs transBetweenGroups,
    	            Sampler samplerFirst, 
    	            Sampler[] samplers,
    	            AbstractTransitionProbs[] transWithinGroups, 
    	            int[] stateToGroup,
    	            int[] stateToIndexWithinGroup,
    	            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState) throws Exception{
    	       super(samplerFirst, samplers,transWithinGroups, stateToGroup, stateToIndexWithinGroup, exp_p, clazz, stateToGroup.length, groupToState);
    	       this.groupToState = groupToState;
    	    }
    public BetweenWithinTransitionProbs3( 
            double[] hittingProbs,
    	            AbstractTransitionProbs transBetweenGroups,double[] u,
    	          //  Sampler samplerFirst, 
    	            Sampler[] samplers,
    	            int[] stateToGroup,
    	            int[] stateToIndexWithinGroup,
    	            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState,
    	           // int[][] stateToGroupTrans, 
    	          //  int[][] stateToWithinGroupTrans,
    	            int i) throws Exception{
    	       super( transBetweenGroups, u,samplers, stateToGroup, stateToIndexWithinGroup, exp_p, clazz, groupToState
    	    		  // stateToGroupTrans, stateToWithinGroupTrans
    	    		   );
    	       this.groupToState = groupToState;
    	      //  if(this.transBetweenGroups instanceof FreeTransitionProbs1){
    	    	//   ((FreeTransitionProbs1)this.transBetweenGroups).harmonise(this.groupToState, hittingProbs);
    	       // }
    	    }
    
    public BetweenWithinTransitionProbs3( 
            double[] hittingProbs,
    	            AbstractTransitionProbs transBetweenGroups,double[] u,
    	          //  Sampler samplerFirst, 
    	            Sampler[] samplers,
    	            AbstractTransitionProbs[] transWithin,
    	            int[] stateToGroup,
    	            int[] stateToIndexWithinGroup,
    	            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState,
    	            int[][] stateToGroupTrans, 
    	            int[][] stateToWithinGroupTrans,
    	            int i) throws Exception{
    	       super( transBetweenGroups, u,samplers, transWithin, stateToGroup, stateToIndexWithinGroup, exp_p, clazz, groupToState
//    	    		   stateToGroupTrans, stateToWithinGroupTrans
    	    		   );
    	       this.groupToState = groupToState;
    	     //   if(this.transBetweenGroups instanceof FreeTransitionProbs1){
    	    //	   ((FreeTransitionProbs1)this.transBetweenGroups).harmonise(this.groupToState, hittingProbs);
    	     //   }
    	    }
   

    public BetweenWithinTransitionProbs3(BetweenWithinTransitionProbs3 probs, boolean swtch){
        super(probs, swtch);
        this.groupToState = ((BetweenWithinTransitionProbs3)probs).groupToState;
    }
   
  


    
}
