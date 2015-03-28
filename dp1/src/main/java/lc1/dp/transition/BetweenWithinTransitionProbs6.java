package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.Dirichlet;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.TrainableGammaDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
public class BetweenWithinTransitionProbs6 extends AbstractTransitionProbs implements Between {

int[][] groupToState;
//   final private SimpleExtendedDistribution1 exp_rd;
//   final  SimpleExtendedDistribution alpha;
   

//public  double[] hp=null;
final SimpleExtendedDistribution[] frac;
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



/*public void updateBetween(){
	for(int i=0; i<frac.length; i++){
		double t = this.transBetweenGroups.getTransition(i, i);
		for(int indexFrom=0; indexFrom<frac[i].probs.length; indexFrom++){
			this.setProb(i, i, indexFrom, t*frac[i].probs(indexFrom));
		}
	}
}*/
	//final SimpleExtendedDistribution[] frac;
public FreeTransitionProbs1[] freeTransBetweenGroups; //[groupFrom, probs are then from within to fromGroup
MatrixExp[] global;


public double getRate(int j1){
	int j = j1+1;
	if(Constants.svnTransM()==0){
	
	return this.freeTransBetweenGroups[this.stateToGroup[j]].getRate(this.stateToIndexWithinGroup[j]);
	}
	else{
		return this.transBetweenGroups.getRate(this.stateToGroup[j]-1);
	}
}





final FreeRateTransitionProbs1[][]sameSize;

public BetweenWithinTransitionProbs6(
			int[] stateToGroup,
			int[][] groupToState,
			int[] stateToIndexWithinGroup,
			AbstractTransitionProbs betweenGroups, AbstractTransitionProbs[] transWithin,Object[] clazz, Sampler[] samplers, 
			Sampler[] exp_p,  double[] hp, double dist, MatrixExp[] global, MatrixExp global1,
			TrainableGammaDistribution sn, TrainableGammaDistribution sn1, int[][]sameSize) {
		this.transBetweenGroups = betweenGroups;
		this.global = global;
		
		this.transWithinGroups = new AbstractTransitionProbs[transWithin.length][transWithin.length];
		this.freeTransBetweenGroups = new FreeTransitionProbs1[groupToState.length];
	
		this.pseudo = new double[groupToState.length];
		this.pseudoC = new double[groupToState.length];
		if(Constants.svnTransM()==0 && ! Constants.diffRatesPerState()){
			this.sameSize = new FreeRateTransitionProbs1[sameSize.length][];
			for(int i2=0; i2<sameSize.length; i2++){
				this.sameSize[i2] = new FreeRateTransitionProbs1[sameSize[i2].length];
				double[] logrelativeRate = new double[groupToState[sameSize[i2][0]].length];
				double rate = Math.log(transBetweenGroups.getRate(0)/global1.currentRate);
				for(int i1=0; i1<logrelativeRate.length; i1++){
					logrelativeRate[i1] = lc1.stats.NormalDistribution.quantile(Constants.rand.nextDouble(), rate, 0.1);
				}
			   for(int i1=0; i1<sameSize[i2].length; i1++){
				   int i = sameSize[i2][i1];
					this.getTrans(i);
				   this.sameSize[i2][i1] = new FreeRateTransitionProbs1(rate, groupToState[i], i,new Dirichlet(pseudo.clone(),Constants.expand_init_prior(3)) , 
							sn, sn1,dist, global1,logrelativeRate);
				   this.freeTransBetweenGroups[i] = this.sameSize[i2][i1];
			   }
			}
		}
		else{
			this.sameSize=null;
		for(int i=1; i<freeTransBetweenGroups.length; i++){
			this.getTrans(i);
			
			double rate = Math.log(transBetweenGroups.getRate(i-1)/global1.currentRate);
			double[] logrelativeRate = new double[groupToState[i].length];
			for(int i1=0; i1<logrelativeRate.length; i1++){
				logrelativeRate[i1] = lc1.stats.NormalDistribution.quantile(Constants.rand.nextDouble(), rate, 0.1);
			}
			freeTransBetweenGroups[i] =
				Constants.svnTransM()==0  ? 
					new FreeRateTransitionProbs1(rate, groupToState[i], i,new Dirichlet(pseudo.clone(),Constants.expand_init_prior(3)) , 
							sn, sn1,dist, global1,logrelativeRate) : 
				
				new FreeTransitionProbs1(groupToState[i].length,new Dirichlet(pseudo.clone(), 
					Double.POSITIVE_INFINITY));
		}
		}
		this.stateToGroup = stateToGroup;
		this.groupToState = groupToState;
		
		this.stateToIndexWithinGroup= stateToIndexWithinGroup;
		this.frac = new SimpleExtendedDistribution[groupToState.length];
		
	
		for(int i=0; i<transWithinGroups.length; i++){
			int grpsize = this.groupToState[i].length;
			frac[i] = new SimpleExtendedDistribution(groupToState[i].length);
			for(int i1=0; i1<transWithinGroups.length; i1++){
				if(i==i1 || (Constants.sameSizeWithin() && groupToState[i].length==groupToState[i1].length)){
					
					if(transWithin[i]==null){
						try{
					    if((i>0 && Constants.modify(0)[i-1]=='0' && Constants.svnTransM()!=0) || grpsize==1){
					    	int grpsize1 = this.groupToState[i1].length;
					    	this.transWithinGroups[i][i1] = new FreeTransitionProbs1(grpsize, this.getSampler(grpsize1, 0.5,
					    			Constants.expand_init_prior(4),false));
					    	{
					    		((FreeTransitionProbs1) this.transWithinGroups[i][i1]).setToIdentity();
					    	}
					    		
					    }
					    else{
					    	double [] pi = new double[grpsize];
					    	Arrays.fill(pi, 1.0/(double)pi.length);
					    	this.transWithinGroups[i][i1] =
					    		ExponentialTransitionProbs.get( clazz[1],samplers[i], exp_p[1], samplers[i].dist.length, Constants.expModelIntHotSpot1(1));
			    	           
							//new FreeExpTransitionProbs1(dist, Constants.transform(samplers[i].dist),Constants.expModelIntHotSpot[1] * Constants.probCrossOverBetweenBP );
					    }
						//	
					//	 = new FreeTransitionProbs1(grpsize, Double.POSITIVE_INFINITY);
						}catch(Exception exc)
						{exc.printStackTrace();}
					}
					else{
						this.transWithinGroups[i][i1] = new StartRemovedTransitionProbs(transWithin[i]);
					}
				}
				else{
					int grpsize1 = this.groupToState[i1].length;
					Dirichlet dir;
				//	double[] d = new double[grpsize1];
					int primary_index=Math.min(i-1,grpsize1-1);
					/*if(i==1 && i1==2){
					//    d[0] = 1.0;
						primary_index = 0;
						dir = new Dirichlet(d, Double.POSITIVE_INFINITY);
					}
					else if(i==3 && i1 ==2){
						d[d.length-1] =1.0;
						dir = new Dirichlet(d, Double.POSITIVE_INFINITY);
					}*/
				//	else{*/
					//if(i==2 ){
						dir = getSampler(grpsize1, 100,1,primary_index);
					/*}
					else{
					dir = this.getSampler(grpsize1, 1,0.1);
					}*/
				//	}
			 		this.transWithinGroups[i][i1] = new FreeTransitionProbs1(grpsize, dir);
				}
			}
		}
		
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
		this.setHP(hp);
	if(Constants.svnTransM()>0)	{
		this.optimiseBack(true,Constants.svnTransM(),1);
	}
	//	for(int i=0; i<freeTransBetweenGroups.length; i++){
	//		freeTransBetweenGroups[i].initialiseCounts(false,false);
	//	}
	//	this.pseudo = new double[groupToState.length];
	}


	
	
	public  static Dirichlet getSampler(int length, double v, double v1, boolean rev) {
		double[] prob = new double[length];
		for(int i=0; i<length; i++){
			prob[i] = Math.pow(v, rev ? length-1-i : i);
		}
		Constants.normalise(prob);
		return new Dirichlet(prob, v1);
	}
	
	
	public  static Dirichlet getSampler(int length, double v, double v1, int primary) {
		double[] prob = new double[length];
		for(int i=0; i<length; i++){
			prob[i] = i==primary ? v : 1;
		}
		Constants.normalise(prob);
		return new Dirichlet(prob, v1);
	}

	public BetweenWithinTransitionProbs6(
				BetweenWithinTransitionProbs6 betwith,
				boolean swtch) {
		this.pseudo = betwith.pseudo;
		this.pseudoC = betwith.pseudoC;
		this.frac = betwith.frac;
	//	this.frac = betwith.frac;
		this.sameSize = betwith.sameSize;
		 this.transBetweenGroups = betwith.transBetweenGroups;
		 this.transWithinGroups = betwith.transWithinGroups;
		//	this.backward = betwith.backward;
			this.stateToGroup = betwith.stateToGroup;
			this.groupToState = betwith.groupToState;
			this.stateToIndexWithinGroup= betwith.stateToIndexWithinGroup;
		//	this.pseudo = new double[groupToState.length];
		//	this.pseudo1 = new double[groupToState.length];
			
			// this.within = betwith.within;
			// TODO Auto-generated constructor stub
		}
	
  

final public AbstractTransitionProbs transBetweenGroups;
final public AbstractTransitionProbs[][] transWithinGroups;
 //  final public SimpleExtendedDistribution[][] backward;  //[to groups, from groups, prob over states within each group]
   public double[] getAlphaPrior(){
	   return this.transBetweenGroups.getAlphaPrior();
   }
  
  public String info() {
      StringBuffer info = new StringBuffer(super.info());
      info.append("between "+transBetweenGroups.info()+"  within "+transWithinGroups[1][1].info());
      return info.toString();
  }
  @Override
  public AbstractTransitionProbs clone(boolean swtch){
      
     return new BetweenWithinTransitionProbs6(this, swtch);
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
      
        double groupSc =groupFrom==0 ? transBetweenGroups.getTransition(groupFrom, groupTo) : this.freeTransBetweenGroups[groupFrom].getTransition(indexFrom, groupTo);
        //	this.transBetweenGroups.getTransition(groupFrom, groupTo);
  //    if(Double.isNaN(groupSc)){
  //  	  throw new RuntimeException("!!");
   //   }
      
        	
           
            double gt = transWithinGroups[groupFrom][groupTo].getTransition(indexFrom, indexTo);;
           double res = groupSc * gt;
           if(Constants.CHECK && res < 0){
        	   throw new RuntimeException("!!");
           }
          /* if(groupFrom==1 && groupTo==1 && from!=to){
        	   System.err.println("h");
           }*/
           return res;
          //  return exp * (transWithinGroups[groupFrom].getTransition(indexFrom, indexTo))  //within group prob
       
        
    }
   
    
 
    
    @Override
    public final void addCount(int from, int to, double val, double dist) {
        int groupFrom = stateToGroup[from];
        int indexFrom = stateToIndexWithinGroup[from];
        int groupTo  = stateToGroup[to];
        int indexTo = stateToIndexWithinGroup[to];
        this.transBetweenGroups.addCount(groupFrom, groupTo, val, dist);
      // if(val>0.5 && groupFrom ==2 && groupTo==2 && indexFrom!=indexTo){
    //	   Logger.global.info("h");
     //  }
        this.freeTransBetweenGroups[groupFrom].addCount(indexFrom, groupTo, val,dist);
            this.transWithinGroups[groupFrom][groupTo].addCount(indexFrom, indexTo, val,dist);
       
        
    }
    
    @Override
    public final void addCount(int from, int to, double val) {
      this.addCount(from, to, val,0);
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
    	   if(i>0)this.freeTransBetweenGroups[i].initialiseCounts(start, end);
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
    
   // public final double[] pseudo, pseudo1;
    public final double[] pseudo, pseudoC;
   public  double[]  getTrans(int group) {
    	//double[] pseudo = new double[this.groupToState.length];
		for(int i=0; i<pseudo.length; i++){
		  pseudo[i] = this.transBetweenGroups.getTransition(group, i);
		}
	//	if(group==2 && pseudo[1]>0.01){
	//		System.err.println("h");
	//	}
		if(group==0 && Constants.sum(pseudo)==0){
			Arrays.fill(pseudo,1.0/(double)pseudo.length);
		}
		if(Constants.CHECK && pseudo[Constants.getMin(pseudo)]<0){
			throw new RuntimeException("!!");
		}
		return pseudo;
	}
   
   public  double[]  getCounts(int group) {
   	//double[] pseudoC = new double[this.groupToState.length];
	   double sum=0;
		for(int i=0; i<pseudoC.length; i++){
		  pseudoC[i] = this.transBetweenGroups.getTransitionCount(group, i);
		  sum+=pseudoC[i];
		}
	//	if(group==2 && pseudoC[1]>0.01){
	//		System.err.println("h");
	//	}
		if(group==0 && sum==0){
			Arrays.fill(pseudoC,1.0/(double)pseudoC.length);
		}
		if(sum>0){
			for(int i=0; i<pseudoC.length; i++){
				pseudoC[i] = pseudoC[i]/sum;
			}
		}
		if(Constants.CHECK && pseudoC[Constants.getMin(pseudoC)]<0){
			throw new RuntimeException("!!");
		}
		return pseudoC;
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
						}
					harmonise(gToS,  0,i);
					
					 if(dist>0.1){

						 double distA = getDist(i);
						 System.err.println(" target "+Constants.print(pseudo)+"\n actual "+Constants.print(pseudo1));
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
	public SimpleExtendedDistribution frac(int groupFrom) {
		return this.frac[groupFrom];
	}

	@Override
	public int[][] groupToState() {
		return this.groupToState;
	}
  
	BackCalc[] backCalc;
	public double optimiseBack(boolean initial, double ent, double pseudo){
 	   double logL=0;
 	   if(initial){
 		  backCalc = new BackCalc[groupToState.length];
 	   }
 	   for(int i=0; i<groupToState.length; i++){
 			 if(this.frac[i].probs.length==1 ){
 				 this.setUniform(i);
 			 }
 			 else{
 				 double thresh = 1e-4;
 				 double numGM =0;
 				double[] trans = this.getTrans(i);
 				double[] counts = initial ? null : this.getCounts(i);
 				 int num=0;
 				 for(int i1=1; i1<trans.length; i1++){
 					 if(trans[i1]>thresh){
 						 num++;
 					 }
 					 if(trans[i1]>1e-4){
 						 numGM++;
 					 }
 				 }
 				 if(initial){
 				 if(num>1 && numGM>1){
 				int[]  grps = new int[num];
 				int[] grps1 = new int[trans.length-num];
 				 {
 					 int k=0; int k1=0;
 					 for(int i1=1; i1<trans.length; i1++){
 						 if(trans[i1]>thresh){
 							 grps[k] = i1;
 							 k++;
 						 }
 						 else{
 							 grps1[k1] = i1;
 							 this.setUniformTo(i,i1, trans[i1]);
 							 k1++;
 						 }
 					 }
 				 }
 				
 				backCalc[i]  = new BackCalc(this,i, trans, grps, grps1);
 				BackCalc mvf = backCalc[i];
 				{
	 				
	 				 double[] xvec = new double[mvf.x_init.length];
	 				 System.arraycopy(mvf.x_init, 0, xvec, 0, xvec.length);
	 				//if(initial ){
	 					System.err.println("before "+ mvf.n+" "+mvf.evaluate(mvf.x_init));
	 					mvf.initialise(xvec);
	 					System.err.println("after "+ mvf.n+" ");
	 				//}
	 				
	 				
 				}
 				 }else{
 				 this.setUniform(i);
 				 }
 				 
 				 }else if(backCalc[i]!=null){
 					BackCalc mvf = backCalc[i];
 					 double[] xvec = new double[mvf.x_init.length];
 					 for(int k=0; k<xvec.length; k++){
 						 xvec[k] = mvf.x.get(k);
 					 }
 					 logL+=mvf.update(trans, counts,ent, pseudo);
	 				/*	 MultivariateMinimum mvm = new pal.math.OrthogonalSearch();
	 					boolean nan = false;
	 					try{
	 				 mvm.optimize(mvf, xvec, 0.01, 0.01);
	 				 
	 					}catch(Exception exc){
	 						Logger.global.warning(exc.getMessage());
	 						nan = true;
	 					}
	 				 double ll = nan ? -1.0*mvf.evaluate(mvf.x_init) :  -1.0*mvf.evaluate(xvec);
	 				 logL+=ll;*/
	 				
 				 }
 			 }
 	   }
 	   return logL;
    }
    
   
    
   private void getExtreme(DoubleMatrix1D xvec, int numStates) {
	  int offset =0;
	  int index = 0;
	  for(int i=0; i<numStates; i++){
		  if(i==index) continue;
		  double d = xvec.get(offset+i);
		  xvec.set(offset+index,xvec.get(offset+index)+d);
		  xvec.set(offset+i,xvec.get(offset+i)-d);
		  xvec.set(numStates+index, xvec.get(numStates)-d);
		  xvec.set(numStates+i,xvec.get(numStates+i)+d);
	  }
	  offset = 2*numStates;
	  index = numStates -1;
	  for(int i=0; i<numStates; i++){
		  if(i==index) continue;
		  double d = xvec.get(offset+i);
		  xvec.set(offset+index,xvec.get(offset+index)+d);
		  xvec.set(offset+i,xvec.get(offset+i)-d);
		  xvec.set(numStates+index, xvec.get(numStates)-d);
		  xvec.set(numStates+i,xvec.get(numStates+i)+d);
	  }
	  System.err.println(xvec);
		
	}


private double[] avg(double[] xvec, double[] xmin, double[] ds) {
		double[] res = new double[xvec.length];
		for(int i=0; i<res.length; i++){
			res[i] = (xvec[i]*ds[0]+xmin[i]*ds[1])/2.0;
		}
		return res;
	}




public  double evaluteBack() {
    	double logL =0;
    	double entropy = 0;
    	for(int i=1; i<this.freeTransBetweenGroups.length; i++){
    		
    		
			
			 logL+=freeTransBetweenGroups[i].logprob();
		
			 if(Constants.CHECK && Double.isNaN(logL)){
				 throw new RuntimeException("!!");
			 }
		}
    	return logL;
	}
    public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
    	if(pos_index==2){
    		double logL = this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, 0);
    		 for(int i=0; i<transWithinGroups.length; i++){
    			 for(int i1=0; i1<transWithinGroups.length; i1++){
    				 boolean same =  i==i1 || (Constants.sameSizeWithin() && this.groupToState[i].length==groupToState[i1].length);
    			 ///note should be pos_index + 1, but this would lead to 3rd index
                 logL+=this.transWithinGroups[i][i1].transferAlpha(pseudoTrans, null,same ? 1 : 2);//, alpha_overall, pos_index+1);
    			 }
    		 }
    		logL+=this.optimiseBack(false,1.0, pseudoTrans[2]);
    		 return logL;
    	}
    	else{
    		return this.transBetweenGroups.transferAlpha( pseudoTrans, alpha_overall, pos_index);
    	}
		
	}
    
 
    
    @Override
    public double  transferQ(double[] ds, double pseudoAlpha, double pseudoRate, MatrixExp initial, int pos_index, double distance, int it) {
		if(pos_index==2){
		
			double logL = this.transBetweenGroups.transferQ( ds, pseudoAlpha,  pseudoRate,initial,0,distance,it);
			
		        for(int i=0; i<transWithinGroups.length; i++){
		        	//should be pos_index+1
		               // logL+= this.transWithinGroups[i][i].transferQ( ds, 0,0, this.global[i],2,distance);
		                for(int i1=0; i1<transWithinGroups.length; i1++){
		                	boolean same = i==i1 || (Constants.sameSizeWithin() && this.groupToState[i].length==this.groupToState[i1].length);
		                	if(this.groupToState[i1].length>1){
		                	if( same){
		          			 ///note should be pos_index + 1, but this would lead to 3rd index
		                       logL+=this.transWithinGroups[i][i1].transferQ( ds, 0,0,
		                    		   same ? this.global[i] : null,
		                    		   2,distance,it);//, alpha_overall, pos_index+1);
		                	}
		                	else{
		                		this.transWithinGroups[i][i1].transferAlpha(ds, null, 2);
		                	}
		                	}
		          	     }
		        }
		        if(this.sameSize!=null){
		        	for(int i=0; i<this.sameSize.length; i++){
		        		FreeRateTransitionProbs1.transferQ( sameSize[i], ds, 0,0, initial, 1, distance);
		        	}
		        }
		        else if(Constants.svnTransM()==0){
		        
		        	
		        	for(int i=1; i<this.freeTransBetweenGroups.length; i++){
		        		freeTransBetweenGroups[i].transferQ(ds, 0.0, 0.0, initial, 1, distance,it);
		        		
		        	}
		        
		        	
		      
		        	
		        }
		        else{
		       logL+=this.optimiseBack(false,Math.pow(it -Constants.numIt[0],1.0),ds[1]);
		        }
		    	return logL;
			}
			else{
				 return this.transBetweenGroups.transferQ( ds,pseudoAlpha, pseudoRate, initial,0,distance,it);
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
   
    
    public final void validate1(boolean check){
    	for(int i=0; i<this.freeTransBetweenGroups.length; i++){
    		if(i>0){
    		this.freeTransBetweenGroups[i].validate(check);
    		}
    	}
    }
    public final void validate() {
    	validate1(true);
        this.transBetweenGroups.validate();
//       alpha.validate();
 //      validateExp();
    
       for(int i=0; i<transWithinGroups.length; i++){
    	   for(int i1=0; i1<transWithinGroups.length;i1++){
           this.transWithinGroups[i][i1].validate();
    	   }
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
            if(this.transWithinGroups[i][i].noStates()>1){
            pw.println("trans Within Groups "+i);
            this.transWithinGroups[i][i].print(pw, null, dist);
            }
        }
     // printExp(pw, dist);
      //if(Constants.annotate()) 
          pw.print("trans between groups: ");
      //this.transBetweenGroups.print(pw, hittingProb, dist);
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
		if(groupFrom!=0){
	this.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom].setProbs(groupTo, prob);
		}
	}
	
	public double  getProb(int groupFrom, int groupTo, int indexFrom) {
		return this.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom].probs(groupTo);
		
		}
	
	public void setUniform(int groupFrom){
		double[] trans = this.getTrans(groupFrom);
		for(int j=0; j<this.groupToState.length; j++){
			this.setUniformTo(groupFrom,j, trans[j]);
		}
	}
	
	public void setUniformTo(int groupFrom, int groupTo, double transj){
		for(int i=0; i<this.frac[groupFrom].probs.length; i++){
			this.setProb(groupFrom, groupTo, i, 
				 transj);
			}
		}




	@Override
	public double getGroupProb(int groupFrom, int groupTo) {
		return this.transBetweenGroups.getTransition(groupFrom, groupTo);
	}


	@Override
	public void validate1() {
		this.validate1(true);
		
	}
	
   

    
}
