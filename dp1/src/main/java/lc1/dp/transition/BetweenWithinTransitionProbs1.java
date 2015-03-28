package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.Dirichlet;
import lc1.stats.Sampler;
import lc1.util.Constants;
public class BetweenWithinTransitionProbs1 extends AbstractTransitionProbs {

//   final private SimpleExtendedDistribution1 exp_rd;
//   final  SimpleExtendedDistribution alpha;
   
   final public AbstractTransitionProbs transBetweenGroups;
   public double[] getAlphaPrior(){
	   return this.transBetweenGroups.getAlphaPrior();
   }
 // final   ExpTransProb[] alphaWithinGroup;
  final   AbstractTransitionProbs[][] transWithinGroups;
  public String info() {
      StringBuffer info = new StringBuffer(super.info());
      info.append("between "+transBetweenGroups.info()+"  within "+transWithinGroups[1][1].info());//+" alpha wihtin "+alphaWithinGroup[1].info());
      return info.toString();
  }
  @Override
  public AbstractTransitionProbs clone(boolean swtch){
      
      if(swtch) return new BetweenWithinTransitionProbs3(this);
      else return new BetweenWithinTransitionProbs1(this, swtch);
  }
  /*  public BetweenWithinTransitionProbs1(BetweenWithinTransitionProbs probs, boolean swtch){
        this.transLength   = probs.transLength;
        this.transBetweenGroups = probs.transBetweenGroups.clone(swtch);
       
        alphaWithinGroup = new ExpTransProb[probs.alphaWithinGroup.length];
        for(int i=0; i<alphaWithinGroup.length; i++){
            alphaWithinGroup[i] = alphaWithinGroup[i].clone(swtch, probs.noStates());
        }
        stateToGroup  = new int[probs.stateToGroup.length];
        System.arraycopy(probs.stateToGroup.length, 0, stateToGroup, 0, stateToGroup.length);
        stateToIndexWithinGroup  = new int[probs.stateToIndexWithinGroup.length];
        System.arraycopy(probs.stateToIndexWithinGroup.length, 0, stateToIndexWithinGroup, 0, stateToIndexWithinGroup.length);
        transWithinGroups = new AbstractTransitionProbs[probs.transWithinGroups.length];
        for(int i=0; i<transWithinGroups.length; i++){
            transWithinGroups[i] = probs.transWithinGroups[i].clone(swtch);
        }
    }*/
   // int[][] indices; //i.e. indices[i] are the states which are part of free_trans;
    //boolean[][] membership;
    @Override
    public int noStates() {
       return stateToGroup.length;
    }
   public final  int[] stateToGroup;
   public final int[] stateToIndexWithinGroup;
  //  public SimpleExtendedDistribution getExp(int groupFrom){
   //     return exp_rd;
   // }
   
   public double getBetweeenGroupTrans(int from, int groupFrom, int groupTo){
      return  transBetweenGroups.getTransition(groupFrom, groupTo);
   }
   
    public final double getTransition(int from, int to) {
        int groupFrom = stateToGroup[from];
        int indexFrom = stateToIndexWithinGroup[from];
        int groupTo  = stateToGroup[to];
        int indexTo = stateToIndexWithinGroup[to];
        double groupSc = getBetweeenGroupTrans(from, groupFrom, groupTo);
       // double exp = getExp(groupFrom).probs[0];
       // double toProb = alpha.probs[groupTo];
        
            return groupSc * transWithinGroups[groupFrom][groupTo].getTransition(indexFrom, indexTo);
          //  return exp * (transWithinGroups[groupFrom].getTransition(indexFrom, indexTo))  //within group prob
           // +(1-exp)*toProb * alphaWithinGroup[groupTo].probs[indexTo];  //between group prob
      
    }
    
    
    
  public void addGroupTransCount(int from, int groupFrom, int groupTo, double val, double dist){
      this.transBetweenGroups.addCount(groupFrom, groupTo, val,dist);
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
  
  @Override
  public final void addCount(int from, int to, double val) {
	  this.addCount(from, to, val,0);
  }
    @Override
    public final void addCount(int from, int to, double val, double dist) {
        int groupFrom = stateToGroup[from];
        int indexFrom = stateToIndexWithinGroup[from];
        int groupTo  = stateToGroup[to];
        int indexTo = stateToIndexWithinGroup[to];
        this.addGroupTransCount(from, groupFrom, groupTo, val,  dist);
//        this.transBetweenGroups.addCount(groupFrom, groupTo, val);
     //   SimpleExtendedDistribution exprd = this.getExp(groupFrom);
       
            this.transWithinGroups[groupFrom][groupTo].addCount(indexFrom, indexTo, val,dist);
      //     double exp =exprd.probs[0];
      //     double non_jump_prob = exp;
      //     double jump_prob = (1 - exp)*alpha.probs[groupTo];
      //     double alloc = val*(jump_prob/(jump_prob+non_jump_prob));
      //     alpha.counts[groupTo]+= alloc;
      //     exprd.counts[1]+=alloc;
        //   alphaWithinGroup[groupTo].counts[indexTo]+=alloc;
        //   exprd.counts[0]+=(val-alloc);
        //   transWithinGroups[groupFrom].addCount(indexFrom, indexTo, val-alloc);
          // transWithinGroups[groupFrom].transitionsOut[indexFrom].counts[indexTo]+=val-alloc;
        
    }

 // final int transLength;
  public  AbstractTransitionProbs makeTransBetweenGroups(Object clazz,Sampler samplerFirst, Sampler exp_p, int no_states, int[] stateToGroup) throws Exception{
      return ExponentialTransitionProbs.get(clazz, samplerFirst, exp_p,
              samplerFirst.dist.length, Constants.expModelIntHotSpot1(0));
  }
  public AbstractTransitionProbs makeTransBetweenGroups(AbstractTransitionProbs transBetweenGroups, boolean swtch, int[] statesToGroup){
      return transBetweenGroups.clone(swtch);
  }
    public BetweenWithinTransitionProbs1( 
          
    //        AbstractTransitionProbs transBetweenGroups,
            Sampler samplerFirst, 
            Sampler[] samplers,
            int[] stateToGroup,
            int[] stateToIndexWithinGroup,
            Dirichlet[] exp_p, Object[] clazz, int no_st, int[][] groupToState
            ) throws Exception{
       
       
       //   this.transLength = transLength;
          this.transBetweenGroups =makeTransBetweenGroups(clazz[0], samplerFirst, exp_p[0], no_st, stateToGroup);
//        this.alpha = new SimpleExtendedDistribution(samplerFirst);
  //      exp_rd = new SimpleExtendedDistribution1(new double[] {exp_p, 1-exp_p}, Double.POSITIVE_INFINITY);
        this.stateToGroup =stateToGroup;
        this.stateToIndexWithinGroup = stateToIndexWithinGroup;
       // this.alphaWithinGroup = new ExpTransProb[samplers.length];
        this.transWithinGroups = new AbstractTransitionProbs[samplers.length][samplers.length];
        for(int i=0; i<samplers.length; i++){
        	int grpsize = groupToState[i].length;
        	for(int j=0; j<samplers.length; j++){
        		if(i==j){
            transWithinGroups[i][i] =ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
        		}
        		else{
        			int grpsize1 = groupToState[j].length;
					this.transWithinGroups[i][j] = new FreeTransitionProbs1(grpsize, BetweenWithinTransitionProbs6.getSampler(grpsize1, 0.5,1.0, false));
        		}
        		}
        	}
     
       
        
        if(Constants.CHECK) this.validate();
    }
    
    public BetweenWithinTransitionProbs1( 
            
    	    //        AbstractTransitionProbs transBetweenGroups,
    	            Sampler samplerFirst, 
    	            Sampler[] samplers,
    	            AbstractTransitionProbs[] transWithinGroups,
    	            int[] stateToGroup,
    	            int[] stateToIndexWithinGroup,
    	            Dirichlet[] exp_p, Object[] clazz, int no_st, int[][] groupToState
    	            ) throws Exception{
    	       
    	       
    	       //   this.transLength = transLength;
    	          this.transBetweenGroups =makeTransBetweenGroups(clazz[0], samplerFirst, exp_p[0], no_st, stateToGroup);
//    	        this.alpha = new SimpleExtendedDistribution(samplerFirst);
    	  //      exp_rd = new SimpleExtendedDistribution1(new double[] {exp_p, 1-exp_p}, Double.POSITIVE_INFINITY);
    	        this.stateToGroup =stateToGroup;
    	        this.stateToIndexWithinGroup = stateToIndexWithinGroup;
    	       // this.alphaWithinGroup = new ExpTransProb[samplers.length];
    	        this.transWithinGroups = new AbstractTransitionProbs[samplers.length][samplers.length];
    	        for(int i=0; i<samplers.length; i++){
    	        	int grpsize = groupToState[i].length;
    	        	for(int j=0; j<samplers.length; j++){
    	        		if(i==j){
    	            this.transWithinGroups[i][i] =ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
    	        		}
    	        		else{
    	        			int grpsize1 = groupToState[j].length;
    						this.transWithinGroups[i][j] = new FreeTransitionProbs1(grpsize, BetweenWithinTransitionProbs6.getSampler(grpsize1, 0.5,1.0,false));
    	        		}
    	        		}
    	        	}
    	      
    	       
    	        
    	        if(Constants.CHECK) this.validate();
    	    }
    
    public AbstractTransitionProbs cloneBetweenGroupsTrans(AbstractTransitionProbs transBetweenGroups
    		, int[] stateToGroup, 
    		//int[][] groupToState,
    		double[] u){
    	return transBetweenGroups.clone(false);
    }
   
    /*statesToWithinGroupTrans is, for each to group, the 'most likely' state within that group to transition into 
     * statesToGroupTrans is, for each state, which groups it is encourage to transition into 
     * */
    
    public BetweenWithinTransitionProbs1( 
            	
    	           AbstractTransitionProbs transBetweenGroups, double[] u,
    	         //   Sampler samplerFirst, 
    	            Sampler[] samplers,
    	            AbstractTransitionProbs[] transWithin,
    	            int[] stateToGroup,
    	            int[] stateToIndexWithinGroup,
    	           
    	            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState
    	           // int[][] statesToGroupTrans, int[][] statesToWithinGroupTrans
    	            ) throws Exception{
    	       
    	       
    	       //   this.transLength = transLength;
    	          this.transBetweenGroups = cloneBetweenGroupsTrans(transBetweenGroups, stateToGroup,
    	        		  //statesToGroupTrans, 
    	        		  u);
//    	        this.alpha = new SimpleExtendedDistribution(samplerFirst);
    	  //      exp_rd = new SimpleExtendedDistribution1(new double[] {exp_p, 1-exp_p}, Double.POSITIVE_INFINITY);
    	        this.stateToGroup =stateToGroup;
    	        this.stateToIndexWithinGroup = stateToIndexWithinGroup;
    	     //   this.alphaWithinGroup = new ExpTransProb[samplers.length];
    	        this.transWithinGroups = new AbstractTransitionProbs[samplers.length][samplers.length];
    	      //  AbstractTransitionProbs[] transWithinGroups = new AbstractTransitionProbs[samplers.length];
    	        for(int i=0; i<samplers.length; i++){
    	        	int grpsize = groupToState[i].length;
    	        	for(int j=0; j<samplers.length; j++){
    	        		if(i==j){
    	            this.transWithinGroups[i][i] =ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
    	        		}
    	        		else{
    	        			int grpsize1 = groupToState[j].length;
    						this.transWithinGroups[i][j] = new FreeTransitionProbs1(grpsize, BetweenWithinTransitionProbs6.getSampler(grpsize1, 0.5,1.0,false));
    	        		}
    	        		}
    	        	}
    	       
    	//       System.arraycopy(transWithinGroups, 0, this.transWithinGroups, 0, transWithinGroups.length);
    	        
    	        if(Constants.CHECK) {
    	        	/*Set<Integer> groups = new HashSet<Integer>();
    	        	for(int i=0; i<stateToGroup.length; i++){
    	        		groups.add(stateToGroup[i]);
    	        	}
    	        	if(groups.size()!=this.transBetweenGroups.noStates()){
    	        		throw new RuntimeException("!!");
    	        	}*/
    	        	//if(this.transBetweenGroups.noStates()!=samplers.length) throw new RuntimeException("!!");
    	        	this.validate();
    	        }
    	    }
    
    public BetweenWithinTransitionProbs1( 
        	
	           AbstractTransitionProbs transBetweenGroups, double[] u,
	         //   Sampler samplerFirst, 
	            Sampler[] samplers,
	            
	            int[] stateToGroup,
	            int[] stateToIndexWithinGroup,
	           
	            Dirichlet[] exp_p, Object[] clazz, int[][] groupToState
	     //       int[][] statesToGroupTrans, int[][] statesToWithinGroupTrans
	            ) throws Exception{
	       
	       
	       //   this.transLength = transLength;
	          this.transBetweenGroups = cloneBetweenGroupsTrans(transBetweenGroups,  stateToGroup,
	        		 //statesToGroupTrans, 
	        		  u);
//	        this.alpha = new SimpleExtendedDistribution(samplerFirst);
	  //      exp_rd = new SimpleExtendedDistribution1(new double[] {exp_p, 1-exp_p}, Double.POSITIVE_INFINITY);
	        this.stateToGroup =stateToGroup;
	        this.stateToIndexWithinGroup = stateToIndexWithinGroup;
	       // this.alphaWithinGroup = new ExpTransProb[samplers.length];
	        this.transWithinGroups = new AbstractTransitionProbs[samplers.length][samplers.length];
	        for(int i=0; i<samplers.length; i++){
	        	int grpsize = groupToState[i].length;
	        	for(int j=0; j<samplers.length; j++){
	        		if(i==j){
	            this.transWithinGroups[i][i] =ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
	        		}
	        		else{
	        			int grpsize1 = groupToState[j].length;
						this.transWithinGroups[i][j] = new FreeTransitionProbs1(grpsize, BetweenWithinTransitionProbs6.getSampler(grpsize1, 0.5,1.0,false));
	        		}
	        		}
	        	}
	       
	        if(Constants.CHECK) {
	        	/*Set<Integer> groups = new HashSet<Integer>();
	        	for(int i=0; i<stateToGroup.length; i++){
	        		groups.add(stateToGroup[i]);
	        	}
	        	if(groups.size()!=this.transBetweenGroups.noStates()){
	        		throw new RuntimeException("!!");
	        	}*/
	        	//if(this.transBetweenGroups.noStates()!=samplers.length) throw new RuntimeException("!!");
	        	this.validate();
	        }
	    }

    public BetweenWithinTransitionProbs1(BetweenWithinTransitionProbs1 probs, boolean swtch){
      
        this.transBetweenGroups =  makeTransBetweenGroups(probs.transBetweenGroups, swtch, probs.stateToGroup);
//       this.exp_rd = new SimpleExtendedDistribution1(probs.exp_rd);
 //      this.alpha = new SimpleExtendedDistribution(probs.alpha);
     //  this.transWithinGroups = new AbstractTransitionProbs[probs.transWithinGroups.length][probs.transWithinGroups.length];
       this.stateToIndexWithinGroup = probs.stateToIndexWithinGroup;
       this.stateToGroup = probs.stateToGroup;
      this.transWithinGroups = probs.transWithinGroups;
       this.validate();
    }
  

  
    public final Collection getDistributions(){
    	throw new RuntimeException("!");
        //List l = new ArrayList();
      /*  l.addAll(transBetweenGroups.getDistributions()); 
        //l.addAll(getExpRdColl());
        l.addAll(Arrays.asList(alphaWithinGroup));
        for(int i=0; i<transWithinGroups.length; i++){
            l.addAll(transWithinGroups[i].getDistributions());
        }
        return l;*/
     }
     
    

  
  

    public final void initialiseCounts(boolean start, boolean end) {
     this.transBetweenGroups.initialiseCounts(start, end);
        // initialiseExpRd();
     //  this.alpha.initialise();
      
       for(int i=0; i<transWithinGroups.length; i++){
    	   for(int i1=0; i1<transWithinGroups.length; i1++){
           this.transWithinGroups[i][i1].initialiseCounts(start, end);
    	   }
       }
    }
 
    /*@Override
    public void transfer(double[] pseudoTrans, double[] pseudoExp) {
        this.transBetweenGroups.transfer(pseudoTrans, pseudoExp[0]);
        for(int i=0; i<alphaWithinGroup.length; i++)
        {
                alphaWithinGroup[i].transferExp(pseudoTrans[1]);
        }
        for(int i=0; i<transWithinGroups.length; i++){
                this.transWithinGroups[i].transfer(pseudoTrans[1], pseudoExp[1]);
        }
        
    }*/
    public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
    	if(pos_index==2){
    		double logL = this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, 0);
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
	public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
		if(pos_index==2){
		double logL = this.transBetweenGroups.transfer( pseudoExp,d,0);
		
	        for(int i=0; i<transWithinGroups.length; i++){
	        	//should be pos_index+1
	                logL+= this.transWithinGroups[i][i].transfer( pseudoExp, null, 1);
	        }
	    	return logL;
		}
		else{
			 return this.transBetweenGroups.transfer( pseudoExp,d,pos_index);
		}
	}
   
	

   /* private void flip(double[] pseudoTrans) {
		double tmp = pseudoTrans[1];
		pseudoTrans[1] = pseudoTrans[0];
		pseudoTrans[0] = tmp;
		
	}
	public final void transfer(double pseudoTrans, double pseudoExp) {
        this.transBetweenGroups.transfer(pseudoTrans, pseudoExp);
        for(int i=0; i<alphaWithinGroup.length; i++)
        {
                alphaWithinGroup[i].transferExp(pseudoTrans);
        }
        for(int i=0; i<transWithinGroups.length; i++){
                this.transWithinGroups[i].transfer(pseudoTrans, pseudoExp);
        }
    }*/

   

    public double transitionDistance(AbstractTransitionProbs probs) {
        throw new RuntimeException("!!");
//        ExponentialFreeTransitionProbs probs1 = (ExponentialFreeTransitionProbs)probs;
  //      return alpha.KLDistance(probs1.alpha) + exp_rd.KLDistance(probs1.exp_rd);
    }
   
    public final void validate() {
        this.transBetweenGroups.validate();
//       alpha.validate();
 //      validateExp();
       
       for(int i=0; i<transWithinGroups.length; i++){
    	   for(int i1=0; i1<transWithinGroups.length; i1++){
           this.transWithinGroups[i][i1].validate();
    	   }
       }
    }
    
   // public Collection getExpRdColl(){
  //      return Arrays.asList(new SimpleExtendedDistribution[] {exp_rd});
  //  }
//public void initialiseExpRd(){
//    this.exp_rd.initialise();
//}
  // public void transferExp(double pseudoExp){
  //     exp_rd.transfer(pseudoExp);
  // }
  // void validateExp(){
  //     exp_rd.validate();
  // }
 
   

    void printExp(PrintWriter pw, double dist){
      //  if(Constants.annotate()) pw.print("exp_probs_0_"+this.getClass());
        this.transBetweenGroups.print(pw, null, dist);
      //  pw.print(transform(exp_rd.probs[0],dist));
        pw.print("; ");
    }
    public final void print(PrintWriter pw, Double[] hittingProb, double dist) {
        if(Constants.annotate()){
           
           
            pw.println("state to group "+print(stateToGroup));
            pw.println("state to index in group "+print(this.stateToIndexWithinGroup));
        }
      
     // printExp(pw, dist);
      //if(Constants.annotate()) 
          pw.print("trans between groups: ");
      this.transBetweenGroups.print(pw, hittingProb, dist);
    //  alpha.print(pw, true, alpha.getPrintString(), "; ");
      if(Constants.annotate()) pw.print("alpha within");
    
      pw.println();
        
    }
    private String print(int[] s) {
        Integer[] res = new Integer[s.length];
        for(int i=0; i<res.length; i++){
            res[i] = s[i];
        }
        return Arrays.asList(res).toString();
    }
    @Override
    public AbstractTransitionProbs clone( int[] statesToGroup, double[] u) {
      throw new RuntimeException("!!");
    }
	
   

    
}
