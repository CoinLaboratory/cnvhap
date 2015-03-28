package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.PermutationSampler;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalSearch;
public class BetweenWithinTransitionProbs5 extends AbstractTransitionProbs implements Between {

 int[][] groupToState;
//   final private SimpleExtendedDistribution1 exp_rd;
//   final  SimpleExtendedDistribution alpha;
   
final SimpleExtendedDistribution[][] within; //[from states, to groups, prob over states withing groups] 

	final SimpleExtendedDistribution[] frac;
	
	public BetweenWithinTransitionProbs5(
			int[] stateToGroup,
			int[][] groupToState,
			int[] stateToIndexWithinGroup,
			AbstractTransitionProbs betweenGroups, AbstractTransitionProbs[] transWithin, double[] hp,
			Object[] clazz, Sampler[] samplers, Sampler[] exp_p) {
		this.transBetweenGroups = betweenGroups;
		this.transWithinGroups = new AbstractTransitionProbs[transWithin.length];
		
		this.stateToGroup = stateToGroup;
		this.groupToState = groupToState;
		this.stateToIndexWithinGroup= stateToIndexWithinGroup;
		//this.pseudo = new double[groupToState.length];
		//this.pseudo1 = new double[groupToState.length];
		int no_states = stateToGroup.length;
		int no_groups = groupToState.length;
		this.within = new SimpleExtendedDistribution[no_states][];
		this.backward = new SimpleExtendedDistribution[no_groups][];
		this.frac = new SimpleExtendedDistribution[no_groups];
		for(int i=0; i<transWithin.length; i++){
			
		//	Arrays.fill(d,1.0/(double)d.length);
			frac[i] = new SimpleExtendedDistribution(groupToState[i].length);
			int grpsize = this.groupToState[i].length;
			if(transWithin[i]==null){
				try{
				this.transWithinGroups[i] = ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
				}catch(Exception exc){
					exc.printStackTrace();
				}
				//this.transWithinGroups[i] = new FreeTransitionProbs1(grpsize, Double.POSITIVE_INFINITY);
			
			}
			else{
				this.transWithinGroups[i] = new StartRemovedTransitionProbs(transWithin[i]);
			}
			backward[i] = new SimpleExtendedDistribution[no_groups];
			for(int k=0; k<no_groups; k++){
			if(i!=k){
					backward[i][k] = new SimpleExtendedDistribution(getSampler(groupToState[i].length,0.5));
				}
			}
		}
		
		for(int i=0; i<within.length; i++){
			within[i] = new SimpleExtendedDistribution[no_groups];
		
			int grp = stateToGroup[i];
			for(int k=0; k<no_groups; k++){
				if(grp!=k){
					within[i][k] = new SimpleExtendedDistribution(getSampler(groupToState[k].length,0.5));
				}
			}
		}
		this.setHP(hp);
		this.optimiseBack();
	/*	 for(int j=0; j<no_states; j++){
	            double[] probs1 = new double[no_states];
	            for(int j1=0; j1<no_states; j1++){
	         //    Arrays.fill(probs, 0.0);
	            }
	            between.transitionsOut[j] =  new SimpleExtendedDistribution1(probs1, Double.POSITIVE_INFINITY);
	        }
		 for(int k=0; k<this.groupToState.length; k++){
			 this.getTrans(k);
		 }*/
		//this.harmonise();
	//	this.pseudo = new double[groupToState.length];
	}
	
	 private Sampler getSampler(int length, double v) {
		double[] prob = new double[length];
		for(int i=0; i<length; i++){
			prob[i] = Math.pow(v, i);
		}
		Constants.normalise(prob);
		return new PermutationSampler(prob, Double.POSITIVE_INFINITY);
	}
	 


	public BetweenWithinTransitionProbs5(
				BetweenWithinTransitionProbs5 betwith,
				boolean swtch) {
		this.frac = betwith.frac;
		 this.transBetweenGroups = betwith.transBetweenGroups;
		 this.transWithinGroups = betwith.transWithinGroups;
			this.backward = betwith.backward;
			this.stateToGroup = betwith.stateToGroup;
			this.groupToState = betwith.groupToState;
			this.stateToIndexWithinGroup= betwith.stateToIndexWithinGroup;
		//	this.pseudo = new double[groupToState.length];
		//	this.pseudo1 = new double[groupToState.length];
			
			 this.within = betwith.within;
			// TODO Auto-generated constructor stub
		}
	
  

final public AbstractTransitionProbs transBetweenGroups;
   final public SimpleExtendedDistribution[][] backward;  //[to groups, from groups, prob over states within each group]
   public double[] getAlphaPrior(){
	   return this.transBetweenGroups.getAlphaPrior();
   }
  final   AbstractTransitionProbs[] transWithinGroups;
  public String info() {
      StringBuffer info = new StringBuffer(super.info());
      info.append("between "+transBetweenGroups.info()+"  within "+transWithinGroups[1].info());
      return info.toString();
  }
  @Override
  public AbstractTransitionProbs clone(boolean swtch){
      
     return new BetweenWithinTransitionProbs5(this, swtch);
  }
    @Override
    public int noStates() {
       return stateToGroup.length;
    }
   public final  int[] stateToGroup;
   public final int[] stateToIndexWithinGroup;
  //  public SimpleExtendedDistribution getExp(int groupFrom){
   //     return exp_rd;
   // }
   
 /*  public double getBetweeenGroupTrans(int from, int groupFrom, int groupTo){
      return  transBetweenGroups.getTransition(groupFrom, groupTo);
   }*/
  

   public final synchronized double getTransition(int from, int to) {
    	
        int groupFrom = stateToGroup[from];
        int indexFrom = stateToIndexWithinGroup[from];
        int groupTo  = stateToGroup[to];
        int indexTo = stateToIndexWithinGroup[to];
        
        double groupSc =this.transBetweenGroups.getTransition(groupFrom, groupTo);
  //    if(Double.isNaN(groupSc)){
  //  	  throw new RuntimeException("!!");
   //   }
        if(groupTo==groupFrom){
        	
           
            double gt = transWithinGroups[groupFrom].getTransition(indexFrom, indexTo);;
           double res = groupSc * gt;
           return res;
          //  return exp * (transWithinGroups[groupFrom].getTransition(indexFrom, indexTo))  //within group prob
           // +(1-exp)*toProb * alphaWithinGroup[groupTo].probs[indexTo];  //between group prob
        }
        else{
        	double back = this.backward[groupFrom][groupTo].probs(indexFrom);
        	double within = this.within[from][groupTo].probs(indexTo);
            double res =( groupSc * back*within) / frac[groupFrom].probs(indexFrom);
            return res;
//            return (1-exp)*toProb*alphaWithinGroup[groupTo].probs[indexTo];  //between group prob
        }
        
    }
  
   
   
   
    
 
    
    @Override
    public final void addCount(int from, int to, double val) {
        int groupFrom = stateToGroup[from];
        int indexFrom = stateToIndexWithinGroup[from];
        int groupTo  = stateToGroup[to];
        int indexTo = stateToIndexWithinGroup[to];
        this.transBetweenGroups.addCount(groupFrom, groupTo, val);
        this.frac[groupFrom].addCount(indexFrom, val);
        
        if(groupTo!=groupFrom){//!state.equals(st)){
        	this.within[from][groupTo].addCount(indexTo, val);
        	  this.backward[groupFrom][groupTo].addCount(indexFrom, val);
        	
        }
        else{
            this.transWithinGroups[groupTo].addCount(indexFrom, indexTo, val);
        }
        
    }

 // final int transLength;
  public  AbstractTransitionProbs makeTransBetweenGroups(Object clazz,Sampler samplerFirst, Sampler exp_p, int no_states, int[] stateToGroup) throws Exception{
      return ExponentialTransitionProbs.get(clazz, samplerFirst, exp_p,
              samplerFirst.dist.length, Constants.expModelIntHotSpot1(0));
  }
  public AbstractTransitionProbs makeTransBetweenGroups(AbstractTransitionProbs transBetweenGroups, boolean swtch, int[] statesToGroup){
      return transBetweenGroups.clone(swtch);
  }
   
    
    public AbstractTransitionProbs cloneBetweenGroupsTrans(AbstractTransitionProbs transBetweenGroups, int[] stateToGroup, int[][] groupToState, double[] u){
    	return transBetweenGroups.clone(false);
    }
   
    /*statesToWithinGroupTrans is, for each to group, the 'most likely' state within that group to transition into 
     * statesToGroupTrans is, for each state, which groups it is encourage to transition into 
     * */
    
   
  
    public final Collection getDistributions(){
       throw new RuntimeException("!!");
     }
     
    

  
  

    public final void initialiseCounts(boolean start, boolean end) {
     this.transBetweenGroups.initialiseCounts(start, end);
        // initialiseExpRd();
     //  this.alpha.initialise();
     
       for(int i=0; i<transWithinGroups.length; i++){
           this.transWithinGroups[i].initialiseCounts(start, end);
       }
       for(int i=0; i<backward.length; i++){
    	   for(int j=0; j<backward[i].length; j++){
    		   if(backward[i][j]!=null) backward[i][j].initialise();
    	   }
       }
       for(int i=0; i<frac.length; i++){
    	   frac[i].initialise();
       }
       for(int i=0; i<this.within.length; i++){
    	   for(int k=0; k<within[i].length; k++){
    		   if(within[i][k] !=null) within[i][k].initialise();
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
    
   // public final double[] pseudo, pseudo1;
   public  double[]  getTrans(int group) {
    	double[] pseudo = new double[this.groupToState.length];
		for(int i=0; i<pseudo.length; i++){
		  pseudo[i] = this.transBetweenGroups.getTransition(group, i);
		}
		return pseudo;
	}
  
    
  /*  public void getTrans1(int state){
   	 //Arrays.fill(pseudo1,0);
   		 for(int j=0; j<this.groupToState.length; j++){
   			 pseudo1[j]=this.between.getTransition(state, j);
   		 }
   }
  
    
    public void getTrans1(int state, double weight){
      		 for(int j=0; j<this.groupToState.length; j++){
      			 pseudo1[j]+=weight*this.backward.getTransition(state, j);
      		 }
      }

    private double[] shorten(double[] pseudo) {
		double[] res = new double[this.groupToState.length];
		Arrays.fill(res,0.0);
		for(int i=0; i<pseudo.length; i++){
			res[this.stateToGroup[i]]+=pseudo[i];
		}
		return res;
	}
   // public final static boolean print = false;
 
  
   
    
   public double getDist(int group){
	  
	   int[] gToS = this.groupToState[group];
		 Arrays.fill(pseudo1, 0.0);
		 double total =0;
			 for(int i1=0; i1<gToS.length; i1++){
				 int state = gToS[i1];
				 double cnt = Constants.sum(backward.transitionsOut[state].counts());
				 this.getTrans1(state, cnt);
				 total+=cnt;
			 }
			 if(total==0) return 0;
			 for(int j=0; j<pseudo1.length; j++){
				 pseudo1[j] = pseudo1[j]/total;
			 }
			 getTrans(group);
			 double dist =  dist(pseudo1, pseudo);;
			 if(dist>0.01){
				 Logger.global.info("h");
			 }
			 return dist;
   }
    
    public void harmonise(int[] gToS,  int index, int group){
		//  System.err.println("harmonising");
    	
		
			
			 for(int i1=0; i1<gToS.length; i1++){
				 int state = gToS[i1];
				 double[] prob = backward.transitionsOut[state].probs();
				  double sum = 0;
			//		if(print) System.err.println("before "+i1+"\t\t\t\t"+Constants.print(shorten(prob)));
				  for(int j=0; j<prob.length; j++){
					 
					  if(pseudo1[j]!=0)  prob[j] = prob[j] *(pseudo[j]/pseudo1[j]);
					  else if(pseudo[j]==0){
						  prob[j] = 0;
					  }
					  sum+=prob[j];
				  }
				Constants.normalise(prob);
			//	if(print) System.err.println("after "+i1+"\t\t\t\t"+Constants.print(shorten(prob)));
			 }
	  }
    private double dist(double[] d, double[] pseudo2) {
		double sum =0;
		for(int i=0; i<d.length; i++){
			if(pseudo2[i]>0) sum+= pseudo2[i] *Math.log(pseudo2[i]/d[i]);
		}
		return sum;
	}

	
   
    
    public void harmonise(){
		  double[] d = new double[groupToState.length];
			 for(int i=0; i<groupToState.length; i++){
					 int[] gToS = groupToState[i];
					 double dist = getDist(i);
					 if(dist>0.1){
							System.err.println(" target "+Constants.print(pseudo)+"\n actual "+Constants.print(pseudo1));
							 Logger.global.info("h");
						}
					harmonise(gToS,  0,i);
					
					 if(dist>0.1){

						 double distA = getDist(i);
						 System.err.println(" target "+Constants.print(pseudo)+"\n actual "+Constants.print(pseudo1));
						 Logger.global.info("h");
						 Logger.global.info("h");
					 }
					 	//if(print) System.err.println("dist is "+dist+ " "+index);
					 //}
				//	 if(pseudo[2]<0.6){
					//	 System.err.println();
				//	 }
					 
				 }
	  }
     */
    @Override
    public void setHP(double[] probs) {
		for(int i=0; i<frac.length; i++){
			frac[i].initialise();
		}
		for(int i=0; i<probs.length; i++){
			int grp = this.stateToGroup[i];
			frac[grp].addCount(this.stateToIndexWithinGroup[i], probs[i]);
		}
		for(int i=0; i<frac.length; i++){
			frac[i].transfer(1e-5);
		}
		
	}
	
 
    public double evaluteBack() {
    	double logL =0;
    	for(int i=0; i<backward.length; i++){
			for(int j=0; j<backward[i].length; j++){
				if(backward[i][j]!=null) logL+=backward[i][j].evaluate(1e-5);
			}
		}
    	return logL;
	}
   public double optimiseBack(){
	   double logL=0;
	   for(int i1=0; i1<groupToState.length; i1++){
			 if(this.frac[i1].probs.length==1){
				 
			 }
			 else{
				 double[] trans =getTrans(i1);
				 int num=0;
				 for(int i=1; i<trans.length; i++){
					 if(trans[i]>1e-3){
						 num++;
					 }
				 }
				 int[] grps = new int[num];
				 {
					 int k=0;
					 for(int i=1; i<trans.length; i++){
						 if(trans[i]>1e-3){
							 grps[k] = i;
							 k++;
						 }
					 }
				 }
				 BackCalc mvf  = null;//new BackCalc(this, i1, trans, grps, grps1);
				 MultivariateMinimum mvm = new OrthogonalSearch();
				 double[] xvec = mvf.x_init;
				 mvm.optimize(mvf, xvec, 0.01, 0.01);
				 double ll = -1.0*mvf.evaluate(xvec);
				 logL+=ll;
			 }
		 }
	   return logL;
   }
    
    public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
    	if(pos_index==2){
    		double logL = this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, pos_index);
    		
    		
    	//	this.harmonise();
    		 for(int i=0; i<transWithinGroups.length; i++){
    			 ///note should be pos_index + 1, but this would lead to 3rd index
                 logL+=this.transWithinGroups[i].transferAlpha(pseudoTrans, null, 1);//, alpha_overall, pos_index+1);
         }
    		 for(int i=0; i<this.within.length; i++){
    			
    			 for(int j=0; j<within[i].length; j++){
    				 if(within[i][j]!=null) logL+=within[i][j].evaluate(1e-5);
    			 }
    		 }
    		// for(int i=0; i<frac.length; i++){
    		//	logL+= frac[i].evaluate(1e-5);
    		// }
    		logL+=optimiseBack();
    		 return logL;
    	}
    	else{
    		return this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, pos_index);
    	}
		
	}
	@Override
	public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
		if(pos_index==2){
		double logL = this.transBetweenGroups.transfer( pseudoExp,d,pos_index);
		
	        for(int i=0; i<transWithinGroups.length; i++){
	        	//should be pos_index+1
	                logL+= this.transWithinGroups[i].transfer( pseudoExp, null, 1);
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
     for(int i=0; i<backward.length; i++){
    	 for(int j=0; j<backward[i].length; j++){
    		if(backward[i][j]!=null)  backward[i][j].validate();
    	 }
     }
     for(int i=0; i<frac.length; i++){
    	frac[i].validate();
     }
     for(int i=0; i<within.length; i++){
    	 for(int j=0; j<within[i].length; j++){
    		if(within[i][j]!=null) within[i][j].validate();
    	 }
     }
       for(int i=0; i<transWithinGroups.length; i++){
           this.transWithinGroups[i].validate();
       }
       /* for(int i=0; i<this.groupToState.length; i++){
    	 
    	  this.getTrans(i);
    	   this.getTrans1(i);
    	   for(int j=0; j<pseudo.length; j++){
    		   if(Math.abs(pseudo[j]-pseudo1[j])>1e-5){
    			   System.err.println(Constants.print(pseudo));
    			   System.err.println(Constants.print(pseudo1));
    			   throw new RuntimeException("!!");
    			   
    		   }
    	   }
       }*/
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
        for(int i=0; i<transWithinGroups.length; i++){
            if(this.transWithinGroups[i].noStates()>1){
            pw.println("trans Within Groups "+i);
            this.transWithinGroups[i].print(pw, null, dist);
            }
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

	public void setProb(int groupFrom, int groupTo, int indexFrom, double prob) {
		// TODO Auto-generated method stub
		 backward[groupFrom][groupTo].setProbs(indexFrom,(prob*this.frac[groupFrom].probs(indexFrom) )/
				 this.transBetweenGroups.getTransition(groupFrom, groupTo));
	}

	@Override
	public SimpleExtendedDistribution frac(int groupFrom) {
		return this.frac[groupFrom];
	}

	@Override
	public int[][] groupToState() {
		return this.groupToState;
	}

	@Override
	public void validate1() {
		throw new RuntimeException("!!");
		
	}

	@Override
	public double getGroupProb(int groupFrom, int groupTo) {
		throw new RuntimeException("");
	}

	@Override
	public double getProb(int groupFrom, int groupTo, int indexFrom) {
		// TODO Auto-generated method stub
		throw new RuntimeException("");
	}
	
   

    
}
