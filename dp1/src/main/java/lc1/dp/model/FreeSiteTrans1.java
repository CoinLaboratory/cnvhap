package lc1.dp.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import lc1.dp.states.State;
import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.BetweenWithinTransitionProbs1;
import lc1.dp.transition.BetweenWithinTransitionProbs3;
import lc1.dp.transition.BetweenWithinTransitionProbs6;
import lc1.dp.transition.ExponentialTransitionProbs;
import lc1.dp.transition.FreeExpTransitionProbs;
import lc1.dp.transition.FreeRateTransitionProbs;
import lc1.dp.transition.FreeRateTransitionProbs1;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.dp.transition.MatrixExp;
import lc1.dp.transition.MultipliedProbs;
import lc1.stats.Dirichlet;
import lc1.stats.PermutationSampler;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;

public class FreeSiteTrans1 extends SiteTransitions{
List<State> states;
    final Object clazz;
    public FreeSiteTrans1(List<Integer> loc, List<State> states, Double[] exp_p1, Double[] r, int length, Object clazz, int index) {
        super(loc, states.size()-1, exp_p1, r, length, index,
        		SiteTransitions.getCN(states), Constants.globalTrans());
        this.states = states;
        this.clazz = clazz;
        // TODO Auto-generated constructor stub
    }

    /** if pseudocount=0 we do not use this as pseudo */
    public FreeSiteTrans1(FreeSiteTrans1 trans_init, boolean swtch){
        super(trans_init, swtch);
        clazz = trans_init.clazz;
//        this(trans_init.loc,trans_init.states, trans_init.exp_p1, trans_init.r,  trans_init.transProbs.length, trans_init.clazz);
    }
    public FreeSiteTrans1(List<State> states, Integer noSnps, Object clazz) {
       this(null, states, null, null, noSnps,clazz,0);
    }

    public FreeSiteTrans1(Integer noSnps, int i, Integer[] cn)throws Exception {
       super(noSnps, i, cn);
       this.clazz = null;
    }
    int[][] groupToState;
    
   
	
    
    @Override
    protected double transferState(SimpleExtendedDistribution dist1, DoubleMatrix2D nullspace, double pseudo) {
    	// TODO Auto-generated method stub
    	if(groupToState!=null){
    	return dist1.evaluate(1e-5, nullspace, this.groupToState);
    	}
    	else{
    		return super.transferState(dist1,nullspace, pseudo);
    	}
    	
     }

    public FreeSiteTrans1(FreeSiteTrans1 trans, int numFounders, int noSnps) throws  Exception{
    	super(noSnps, numFounders+1, trans.cn);
		this.states = new ArrayList(trans.states.subList(0, numFounders+1));
		this.clazz = trans.clazz;
		for(int i=0; i<trans.transProbs.length; i++){
			transProbs[i] = new TruncatedTransitionProbs(trans.transProbs[i], numFounders+1);
		}
	}

	@Override
    public SiteTransitions clone(boolean swtche) {
        return new FreeSiteTrans1(this, swtche);
    }

    public int[] initialise1(int[][] m) throws Exception{
           //super.initialise1();
              Object[] transProbType=
   HaplotypeHMMIterator.getMode(Constants.transMode(2));
            int[] stateToGroup = new int[states.size()];
         
            int[] groupSizes =  new int[m.length];
          int[] transMode2 = Constants.transMode(2);
            int[] stateToIndexWithinGroup = new int[states.size()];
            Double[][] r_new = new Double[Constants.transMode(2).length==1 ? 1 : 2][];
            r_new[0] = r[0];
            for(int i=1; i<r_new.length; i++){
           	 r_new[i] = new Double[]{r[0][0]};
           	 
            }	
            this.r = r_new;
            double[]d2 = new double[this.states.size()];
            Arrays.fill(d2,0);
            /* for(int i=1; i<d2.length; i++){
                  d2[i] = 1.0/(double)(d2.length-1);
             }*/
            for(int i=0; i<m.length; i++){
                for(int j=0; j<m[i].length; j++){
                   int k = m[i][j];
                   stateToGroup[k] = i;
                   stateToIndexWithinGroup[k] = j;
                   double totProb = this.transProbs[0].getTransition(0, i);
                       d2[k] = totProb/(double) m[i].length;
                }
                groupSizes[i] = m[i].length;
           
            }
             groupToState = new int[groupSizes.length][];
            for(int i=0; i<groupSizes.length; i++){
            	groupToState[i] = new int[groupSizes[i]];
            	
            }
            for(int i=0; i<stateToGroup.length; i++){
            	groupToState[stateToGroup[i]][stateToIndexWithinGroup[i]] = i;
            }
        //    double[] u = new double[stateToGroup.length];
           /* for(int i=0; i<u.length; i++){
            	if(groupSizes[stateToGroup[i]]>1) u[i] = Constants.expand_init_prior(1);
            	else u[i] = Double.POSITIVE_INFINITY;;
            }*/
          
            int[][] statesToGroupTrans = new int[stateToGroup.length][];
            int[][] statesToWithinGroupTrans = new int[groupToState.length][];
            makeStateToGroupTrans(stateToGroup,stateToIndexWithinGroup, 
            		groupToState, Constants.transitions(1), 
            		statesToGroupTrans, statesToWithinGroupTrans);
            
            double[] start = new double[states.size()];
            start[0] = 1.0;
         double[] probs = new double[states.size()];
            
            {
               
                           
                         //   this.trans_init = new double[states.size()];
                           
                      
                         
                              
                           // dir1= new Dirichlet(d2,Constants.u_global(0)[1]);
                       
                            transProbs[0] =new FreeTransitionProbs1(true,
   new Dirichlet(d2, 1e10), states.size());
                            if(Constants.CHECK ){
                                transProbs[0].validate();
                            }
                         
                        fillProbs(transProbs[0], probs, start);
                         //   initialise(  d2, Constants.samplePermute());
                          
                                // transProbs[i] =new FreeTransitionProbs(states, false, dir);
                          //  transProbs[transProbs.length-1] =newFreeTransitionProbs( false, dir1, states.size());
                  
            }
         
                Double[] d = new Double[groupSizes.length];
                Arrays.fill(d, 0.0);
                for(int i=1; i<d.length; i++){
                   d[i] = 1.0/(double)(d.length-1);
                }
               
                //for(int k=1; k<d.length; k++){
                 //   d[k] = (double) groupSizes[k] / (double) numFounders;
               // }
              //  double[] u0 = Constants.u_global(0);
            Dirichlet   dir= new Dirichlet(d, Constants.u_global(0)[1]);
            Sampler[] samplers = new Sampler[groupSizes.length];  // within group alpha for transitions in
            for(int k=0; k<groupSizes.length; k++){
                Double[] d1 = new Double[groupSizes[k]];
                double incr = 1.0/(double)groupSizes[k];
                double samplePermute = Constants.samplePermute();
                if(samplePermute>0){
                    int j=0;
                    for(j=0; j<d1.length; j++){
                        d1[j] = //(double) j+1; //
                        Math.pow((double)Constants.samplePermute(), j);
                    }
                //   d1[j-1] = d1[j-1]*Constants.samplePermute();
                    SimpleExtendedDistribution.normalise(d1);
                    samplers[k] = new PermutationSampler(d1,
   Constants.u_global(1)[1]);
                }
                else{
                    Arrays.fill(d1, incr);
                    samplers[k] = new Dirichlet(d1, Constants.u_global(1)[1]);
                }
            }
            //Double[] r_0 = new Double[] {1e-60,r[1]};
            for(int i=1; i<transProbs.length; i++){
              //  Double[] r_ = special==null || this.special.contains(i)?  r : r_0;
                double[] expp =
                   
                  
                    new double[] {
                        loc==null || loc.size()==0 ? exp_p1[0] :
   Math.exp(-1*r[0][0]*(loc.get(i)-loc.get(i-1))),
                        loc==null || loc.size()==0 ? exp_p1[1] :
   Math.exp(-1*r[1][0]*(loc.get(i)-loc.get(i-1)))
                    };
              /* if(special!=null && !special.contains(i)){
                   expp[0] = 1.0;  // 0 is between groups
               }*/
            
                Dirichlet[] exp_p = new Dirichlet[] {
                        new Dirichlet(new double[] {expp[0], 1-expp[0]},
   Constants.u_global(0)[2]),
                        new Dirichlet(new double[] {expp[1], 1-expp[1]},
   Constants.u_global(1)[2]),
                }; //exp_p[0] is between groups
     
              if(transProbs[i]==null)  {transProbs[i] =
                    Constants.trans1() ?
                    new BetweenWithinTransitionProbs3(dir, samplers,
   stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState):
                        new BetweenWithinTransitionProbs1(dir, samplers,
   stateToGroup, stateToIndexWithinGroup, exp_p, transProbType,samplers.length, groupToState);
              }
                        else{
                        	throw new RuntimeException("!!");
                            
                        }
              System.arraycopy(probs,0, start,0,start.length);
              fillProbs(transProbs[i], probs, start);
              transProbs[i].validate();
                if(Constants.CHECK){
                    validate(transProbs[i], states.size(), i);
                }
            }
            return stateToGroup;
    }
    
    public int[] initialise2(SiteTransitions orig, int[][] m, SiteTransitions withinGroup, double initialU) throws Exception{
    	this.allFree = true;
         int[] stateToGroup = new int[states.size()];
      
         int[] groupSizes =  new int[m.length];
         int[] stateToIndexWithinGroup = new int[states.size()];
        
         double[]d2 = new double[this.states.size()];
         Arrays.fill(d2,0);
         /* for(int i=1; i<d2.length; i++){
               d2[i] = 1.0/(double)(d2.length-1);
          }*/
         for(int i=0; i<m.length ; i++){
             for(int j=0; j<m[i].length; j++){
                int k = m[i][j];
                stateToGroup[k] = i;
                stateToIndexWithinGroup[k] = j;
                double totProb = orig.transProbs[0].getTransition(0, i);
                    d2[k] = totProb/(double) m[i].length;
             }
             groupSizes[i] = m[i].length;
        
         }
          groupToState = new int[groupSizes.length][];
         for(int i=0; i<groupSizes.length; i++){
         	groupToState[i] = new int[groupSizes[i]];
         	
         }
         for(int i=0; i<stateToGroup.length; i++){
         	groupToState[stateToGroup[i]][stateToIndexWithinGroup[i]] = i;
         }
      
         
         double[] start = new double[states.size()];
         start[0] = 1.0;
      double[] probs = new double[states.size()];
         
         {
          
                         transProbs[0] =new FreeTransitionProbs1(true,
new Dirichlet(d2, initialU), states.size());
                         if(Constants.CHECK ){
                             transProbs[0].validate();
                         }
                         fillProbs(transProbs[0], probs, start);
         }
             Double[] d = new Double[groupSizes.length];
             Arrays.fill(d, 0.0);
             for(int i=1; i<d.length; i++){
                d[i] = 1.0/(double)(d.length-1);
             }
         
         //Double[] r_0 = new Double[] {1e-60,r[1]};
        // int[][]sameSize = this.getSameSize(groupSizes);
         for(int i=1; i<transProbs.length; i++){
                    	 AbstractTransitionProbs between = orig.transProbs[i];
                    	 
                         transProbs[i] = new MultipliedProbs(between, withinGroup.transProbs[i]);
                        	
           transProbs[i].validate();
             if(Constants.CHECK){
                 validate(transProbs[i], states.size(), i);
             }
         }
         return stateToGroup;
     	
    }
    public int[] initialise1(SiteTransitions orig, int[][] m, SiteTransitions[] withinGroups, double initialU) throws Exception{
    	//Arrays.fill(withinGroups,null);
    	this.allFree = true;
    	AbstractTransitionProbs[] transWithin = new AbstractTransitionProbs[withinGroups.length];
        //super.initialise1();
    	//if(true) throw new RuntimeException("!!");
           Object[] transProbType=
HaplotypeHMMIterator.getMode(Constants.transMode(2));
         int[] stateToGroup = new int[states.size()];
      
         int[] groupSizes =  new int[m.length];
       int[] transMode2 = Constants.transMode(2);
         int[] stateToIndexWithinGroup = new int[states.size()];
         Double[][] r_new = new Double[Constants.transMode(2).length==1 ? 1 : 2][];
         r_new[0] = r[0];
         for(int i=1; i<r_new.length; i++){
        	 r_new[i] = new Double[]{r[0][0]};
        	 
         }	
         this.r = r_new;
         double[]d2 = new double[this.states.size()];
         Arrays.fill(d2,0);
         /* for(int i=1; i<d2.length; i++){
               d2[i] = 1.0/(double)(d2.length-1);
          }*/
         for(int i=0; i<m.length ; i++){
             for(int j=0; j<m[i].length; j++){
                int k = m[i][j];
                stateToGroup[k] = i;
                stateToIndexWithinGroup[k] = j;
                double totProb = orig.transProbs[0].getTransition(0, i);
                    d2[k] = totProb/(double) m[i].length;
             }
             groupSizes[i] = m[i].length;
        
         }
          groupToState = new int[groupSizes.length][];
         for(int i=0; i<groupSizes.length; i++){
         	groupToState[i] = new int[groupSizes[i]];
         	
         }
         for(int i=0; i<stateToGroup.length; i++){
         	groupToState[stateToGroup[i]][stateToIndexWithinGroup[i]] = i;
         }
       //  double[] u = new double[stateToGroup.length];
         MatrixExp[] global = new MatrixExp[groupToState.length];
         for(int i=0; i<global.length; i++){
        	 if(groupToState[i].length>1)
        	 global[i] = make(groupToState[i].length);
         }
         /*for(int i=0; i<u.length; i++){
         	if(groupSizes[stateToGroup[i]]>1) u[i] = Constants.expand_init_prior(1);
         	else u[i] = Double.POSITIVE_INFINITY;;
         }*/
       
        /* int[][] statesToGroupTrans = new int[stateToGroup.length][];
         int[][] statesToWithinGroupTrans = new int[groupToState.length][];
         makeStateToGroupTrans(stateToGroup,stateToIndexWithinGroup, 
         		groupToState, Constants.transitions(1), 
         		statesToGroupTrans, statesToWithinGroupTrans);*/
         
         double[] start = new double[states.size()];
         start[0] = 1.0;
      double[] probs = new double[states.size()];
         
         {
            
                        
                      //   this.trans_init = new double[states.size()];
                        
                   
                      
                           
                        // dir1= new Dirichlet(d2,Constants.u_global(0)[1]);
                    
                         transProbs[0] =new FreeTransitionProbs1(true,
new Dirichlet(d2, initialU), states.size());
                         if(Constants.CHECK ){
                             transProbs[0].validate();
                         }
                         fillProbs(transProbs[0], probs, start);
         }
      
             Double[] d = new Double[groupSizes.length];
             Arrays.fill(d, 0.0);
             for(int i=1; i<d.length; i++){
                d[i] = 1.0/(double)(d.length-1);
             }
            
             //for(int k=1; k<d.length; k++){
              //   d[k] = (double) groupSizes[k] / (double) numFounders;
            // }
           //  double[] u0 = Constants.u_global(0);
         Dirichlet   dir= new Dirichlet(d, Constants.u_global(0)[1]);
         Sampler[] samplers = new Sampler[groupSizes.length];  // within group alpha for transitions in
         for(int k=0; k<groupSizes.length; k++){
             Double[] d1 = new Double[groupSizes[k]];
             if(withinGroups[k]!=null && withinGroups[k].numFounders!=  groupSizes[k]){
            	 throw new RuntimeException("!!");
             }
             double incr = 1.0/(double)groupSizes[k];
             double samplePermute = Constants.samplePermute();
             if(samplePermute>0){
                 int j=0;
                 for(j=0; j<d1.length; j++){
                     d1[j] = //(double) j+1; //
                     Math.pow((double)Constants.samplePermute(), j);
                 }
             //   d1[j-1] = d1[j-1]*Constants.samplePermute();
                 SimpleExtendedDistribution.normalise(d1);
                 samplers[k] = new PermutationSampler(d1,
Constants.u_global(1)[1]);
             }
             else{
                 Arrays.fill(d1, incr);
                 samplers[k] = new Dirichlet(d1, Constants.u_global(1)[1]);
             }
         }
         //Double[] r_0 = new Double[] {1e-60,r[1]};
         int[][]sameSize = this.getSameSize(groupSizes);
         for(int i=1; i<transProbs.length; i++){
        	 for(int k=0; k<transWithin.length; k++){
        		 transWithin[k] = withinGroups[k]==null ?  null : withinGroups[k].transProbs[i];
        	 }
           //  Double[] r_ = special==null || this.special.contains(i)?  r : r_0;
             double[] expp =
                
               
                 new double[] {
                     loc==null || loc.size()==0 ? exp_p1[0] :
Math.exp(-1*r[0][0]*(loc.get(i)-loc.get(i-1))),
                     loc==null || loc.size()==0 ? exp_p1[1] :
Math.exp(-1*r[1][0]*(loc.get(i)-loc.get(i-1)))
                 };
           /* if(special!=null && !special.contains(i)){
                expp[0] = 1.0;  // 0 is between groups
            }*/
         
             Dirichlet[] exp_p = new Dirichlet[] {
                     new Dirichlet(new double[] {expp[0], 1-expp[0]},
Constants.u_global(0)[2]),
                     new Dirichlet(new double[] {expp[1], 1-expp[1]},
Constants.u_global(1)[2]),
             }; //exp_p[0] is between groups
  
           if(orig.transProbs[i]==null)  {transProbs[i] =
                 Constants.trans1() ?
                 new BetweenWithinTransitionProbs3(dir, samplers,transWithin,
stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState):
                     new BetweenWithinTransitionProbs1(dir, samplers,transWithin,
stateToGroup, stateToIndexWithinGroup, exp_p, transProbType,samplers.length, groupToState);
           }
                     else{
                    	 AbstractTransitionProbs between = orig.transProbs[i];
                    	 if(!(between  instanceof FreeTransitionProbs1)){
                    		 between = new FreeTransitionProbs1(between);
                    	 }
                    	 
                 
                    	 double dist = this.getDist(i);
                         transProbs[i] =
                        	 new BetweenWithinTransitionProbs6(stateToGroup, groupToState, 
                        		 stateToIndexWithinGroup, orig.transProbs[i], transWithin, 
                        		 transProbType, samplers, exp_p, probs,dist, global, this.globalTrans.mat, 
                        		 this.rateDistributionG, this.rateDistribution, sameSize);
                    //	new FreeTransitionsProbs2( 
                        			//  new BetweenWithinTransitionProbs3(probs, orig.transProbs[i], u,  samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState,
                                          	//	 statesToGroupTrans, statesToWithinGroupTrans,
                                     //     		 i);
                        			 
                        			
                    	//		)
                        	 ;
                         System.arraycopy(probs,0, start,0,start.length);
                         fillProbs(transProbs[i], probs, start);
                      
                    /*         Constants.trans1() ?
                             new
BetweenWithinTransitionProbs3(probs, between, u, samplers,transWithin,  stateToGroup,
stateToIndexWithinGroup, exp_p, transProbType, groupToState,statesToGroupTrans, statesToWithinGroupTrans, i):
                                 new
BetweenWithinTransitionProbs1(between,u, samplers, transWithin, stateToGroup,
stateToIndexWithinGroup, exp_p, transProbType, groupToState, statesToGroupTrans, statesToWithinGroupTrans);
                     }
          
          */
                     }
           transProbs[i].validate();
        //  if(Constants.newTrans()) transProbs[i] = new BetweenWithinTransitionProbs5((BetweenWithinTransitionProbs3)transProbs[i], orig.transProbs[i],u);
             if(Constants.CHECK){
                 validate(transProbs[i], states.size(), i);
             }
         }
         return stateToGroup;
 }
    
    private int[][] getSameSize(int[] groupSizes) {
		TreeMap<Integer, List<Integer>> m = new TreeMap<Integer, List<Integer>>();
		for(int i=1; i<groupSizes.length; i++){
			List<Integer> l = m.get(groupSizes[i]);
			if(l==null){
				m.put(groupSizes[i], l = new ArrayList<Integer>());
			}
			l.add(i);
		}
		int[][]res = new int[m.size()][];
		int k=0;
		for(Iterator<Integer> it = m.keySet().iterator();it.hasNext();k++){
			Integer v = it.next();
			List<Integer> l = m.get(v);
			res[k] = new int[l.size()];
			for(int i=0; i<l.size(); i++){
				res[k][i] = l.get(i);
			}
		}
		return res;
	}

	public static MatrixExp make(int numStates){
    	double[] pi = new double[numStates];
    	Arrays.fill(pi, 1.0/(double)numStates);
    	return new MatrixExp(pi,Constants.expModelIntHotSpot[1] * Constants.probCrossOverBetweenBP);
    }



public static void makeStateToGroupTrans(int[] stateToGroup, int[] stateToWithinGroup,
    		int[][] groupToState, int[][]trans1, int[][] statesToGroupTrans, int[][] statesToWithinGroupTrans){
    	// makeStateToGroupTrans(stateToGroup, groupToState, trans, 
         //		statesToGroupTrans, statesToWithinGroupTrans);//
    	Integer[][] trans=  new Integer[stateToGroup.length][];
    	Set<Integer>[] transS = new Set[stateToGroup.length];
    	for(int i=0; i<transS.length; i++){
    		transS[i] = new HashSet<Integer>();
    	}
    	for(int i=0; i<trans1.length; i++){
    		for(int j=0; j<trans1[i].length; j++){
    			Set<Integer> transS_j = transS[trans1[i][j]];
    			for(int j1=0; j1<trans1[i].length; j1++){
    				if(j1!=j){
    					transS_j.add(trans1[i][j1]);
    				}
    			}
    		}
    	}
    	for(int i=0; i<transS.length; i++){
    		trans[i] = (new ArrayList<Integer>(transS[i])).toArray(new Integer[0]);
    		
    	}
    	for(int i=1; i<groupToState.length; i++){
        	statesToWithinGroupTrans[i] = new int[stateToGroup.length];
        	Arrays.fill(statesToWithinGroupTrans[i], -1);
        	
        }
        for(int i=1; i<stateToGroup.length; i++){
        	Integer[] trans_i = trans[i];
        	statesToGroupTrans[i] = new int[trans_i.length];
        	for(int j=0; j<trans_i.length; j++){
        		int group= stateToGroup[trans_i[j]];
        		int withinGroup = stateToWithinGroup[trans_i[j]];
        		statesToGroupTrans[i][j] = group;
        		statesToWithinGroupTrans[group][i] = withinGroup;
        	}
        }
//      Logger.global.info("done");
     //   return statesToGroupTrans;
	}

	public void initialise11(int[][] m) throws Exception{
 	   //super.initialise1();
    	   Object[] transProbType= HaplotypeHMMIterator.getMode(Constants.transMode(2));
         int[] stateToGroup = new int[states.size()];
       
         int[] groupSizes =  new int[m.length];
        
         int[] stateToIndexWithinGroup = new int[states.size()];
        
         for(int i=0; i<m.length; i++){
        	 for(int j=0; j<m[i].length; j++){
        		int k = m[i][j];
        		stateToGroup[k] = i;
        		stateToIndexWithinGroup[k] = j;
        	 }
        	 groupSizes[i] = m[i].length;
         
         }
        groupToState = new int[groupSizes.length][];
         for(int i=0; i<groupSizes.length; i++){
         	groupToState[i] = new int[groupSizes[i]];
         	
         }
         for(int i=0; i<stateToGroup.length; i++){
         	groupToState[stateToGroup[i]][stateToIndexWithinGroup[i]] = i;
         }
         double[] start = new double[states.size()];
         start[0] = 1.0;
      double[] probs = new double[states.size()];
         {
         	 
         	            
         	         //   this.trans_init = new double[states.size()];
         	            
         	       
         	          
         	               
         	           // dir1= new Dirichlet(d2, Constants.u_global(0)[1]);
         	          double[]d2 = new double[this.states.size()];
         	          for(int i=1; i<d2.length; i++){
         	        	   d2[i] = 1.0/(double)(d2.length-1);
         	          }
         	            transProbs[0] =new FreeTransitionProbs1(true, new Dirichlet(d2, 1e10), states.size());
         	           
                        fillProbs(transProbs[0], probs, start);
         	            if(Constants.CHECK ){
         	                transProbs[0].validate();
         	            }
         	         //   initialise(  d2, Constants.samplePermute());
         	           
         	                // transProbs[i] =new FreeTransitionProbs(states, false, dir);
         	          //  transProbs[transProbs.length-1] =new FreeTransitionProbs( false, dir1, states.size());
         	   
         }
         int[][] statesToGroupTrans = new int[stateToGroup.length][];
         int[][]  statesToWithinGroupTrans = new int[groupToState.length][];
         FreeSiteTrans1.makeStateToGroupTrans(stateToGroup,stateToIndexWithinGroup, 
          		groupToState, Constants.transitions(0), 
          		statesToGroupTrans, statesToWithinGroupTrans);  
             Double[] d = new Double[groupSizes.length];
             Arrays.fill(d, 0.0);
             for(int i=1; i<d.length; i++){
                d[i] = 1.0/(double)(d.length-1);
             }
             
             //for(int k=1; k<d.length; k++){
              //   d[k] = (double) groupSizes[k] / (double) numFounders;
            // }
           //  double[] u0 = Constants.u_global(0);
         Dirichlet   dir= new Dirichlet(d, Constants.u_global(0)[1]);
         Sampler[] samplers = new Sampler[groupSizes.length];  // within group alpha for transitions in
         for(int k=0; k<groupSizes.length; k++){
             Double[] d1 = new Double[groupSizes[k]];
             double incr = 1.0/(double)groupSizes[k];
             if(Constants.samplePermute()>0){
                 int j=0;
                 for(j=0; j<d1.length; j++){
                     d1[j] = //(double) j+1; //
                     Math.pow((double)Constants.samplePermute(), j);
                 }
             //   d1[j-1] = d1[j-1]*Constants.samplePermute();
                 SimpleExtendedDistribution.normalise(d1);
                 samplers[k] = new PermutationSampler(d1, Constants.u_global(1)[1]);
             }
             else{
                 Arrays.fill(d1, incr);
                 samplers[k] = new Dirichlet(d1, Constants.u_global(1)[1]);
             }
         }
         //Double[] r_0 = new Double[] {1e-60,r[1]};
         for(int i=1; i<transProbs.length; i++){
           //  Double[] r_ = special==null || this.special.contains(i) ?  r : r_0;
             double[] expp = 
                 
                
                 new double[] {
                     loc==null || loc.size()==0 ? exp_p1[0] : Math.exp(-1*r[0][0]*(loc.get(i)-loc.get(i-1))),
                     loc==null || loc.size()==0 ? exp_p1[1] : Math.exp(-1*r[0][1]*(loc.get(i)-loc.get(i-1)))
                 };
           /* if(special!=null && !special.contains(i)){
                expp[0] = 1.0;  // 0 is between groups
            }*/
          
             Dirichlet[] exp_p = new Dirichlet[] {
                     new Dirichlet(new double[] {expp[0], 1-expp[0]}, Constants.u_global(0)[2]),
                     new Dirichlet(new double[] {expp[1], 1-expp[1]}, Constants.u_global(1)[2]),
             }; //exp_p[0] is between groups
   
           if(transProbs[i]==null)  {transProbs[i] = 
                 Constants.trans1() ? 
                 new BetweenWithinTransitionProbs3(dir, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState):
                     new BetweenWithinTransitionProbs1(dir, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType,samplers.length, groupToState);
           }
                     else{
                     	transProbs[i] = 
                             Constants.trans1() ? 
                             new BetweenWithinTransitionProbs3(probs, transProbs[i], null,  samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState,
                            	//	 statesToGroupTrans, statesToWithinGroupTrans,
                            		 i):
                                 new BetweenWithinTransitionProbs1(transProbs[i], null, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, 
                                		 transProbType, groupToState
                                		 //statesToGroupTrans, statesToWithinGroupTrans
                                		 );
                     }
           System.arraycopy(probs,0, start,0,start.length);
           fillProbs(transProbs[i], probs, start);
           transProbs[i].validate();
             if(Constants.CHECK){
                 validate(transProbs[i], states.size(), i);
             }
         }
    
 }
    
    
	
	@Override
    public void initialise(double[] rel1, double permute, double u_glob) throws Exception{
        lc1.stats.Sampler dir1;
      //  double[][] trans = Constants.transitionMatrix();
        double[] rel2 = new double[rel1.length-1];
        System.arraycopy(rel1, 1,rel2, 0, rel2.length);
        if(permute>0.1){
            double[] rel = new double[rel1.length];
           // rel[0] = 0;
            int j;
          //  if(rel.length>1)rel[1]  = permute*rel[1]; 
            for( j=0; j<rel.length; j++){
                double mult = 1.0;//(j+1);
                   rel[j] = rel1[j] * //mult;;//
                     Math.pow(permute, j);
            }
//           rel[j-1]  = permute*rel[j-1]; 
            SimpleExtendedDistribution.normalise(rel);
            dir1 =  new PermutationSampler(rel1, u_glob);
        }
        else dir1 =  new Dirichlet(rel1, u_glob);
      //  if(dir1!=null && dir1.dist.length!=no_states) throw new RuntimeException("!!");
        
            for(int i=1; i<transProbs.length; i++){
            	if(globalTrans!=null){
            	this.globalTrans.mat.setDistance(this.getDist(i));
            	if(clazz.equals(FreeRateTransitionProbs1.class) 
                    	|| (!(clazz instanceof Class) &&  ((Class[])clazz)[0].equals(FreeRateTransitionProbs1.class) )
                            	){
            		  transProbs[i] = 	
            			  Constants.diffRatesPerState() ? 
            					  new FreeRateTransitionProbs1((FreeExpTransitionProbs)this.globalTrans, this.getDist(i), this.rateDistributionG
            							  ,this.rateDistribution):
            			  new FreeRateTransitionProbs((FreeExpTransitionProbs)this.globalTrans, this.getDist(i), this.rateDistributionG);
            			 // new FreeTransitionProbs1(this.globalTrans);
            	}else if(clazz.equals(FreeTransitionProbs1.class) 
                    	|| (!(clazz instanceof Class) &&  ((Class[])clazz)[0].equals(FreeTransitionProbs1.class))){
            		 transProbs[i] = new FreeTransitionProbs1(this.globalTrans);
            	}else{
            		transProbs[i] = new FreeExpTransitionProbs(globalTrans, this.getDist(i));
            			
            			//(
                     	//	trans==null ? 
                     	//new FreeExpTransitionProbs(rel1, this.getDist(i), this.len, Constants.expModelIntHotSpot1(0),this.r[0][0]) :
                     		//new FreeExpTransitionProbs(this.getDist(i), trans, this.r[0][0])
                     //	) 
                     	//)
                     	   ;
            	}
            	}
            	else{
            		 if(clazz.equals(FreeTransitionProbs1.class) || clazz instanceof Class[] && ((Class[])clazz)[0] ==FreeTransitionProbs1.class){
            	            transProbs[i] =  new FreeTransitionProbs1(dir1);
            	        }
            		 else{
            		transProbs[i] = 
            			 ExponentialTransitionProbs.get( clazz, dir1, this.getExp_p(i), this.len, Constants.expModelIntHotSpot1(0));
            		//new FreeExpTransitionProbs(this.getDist(i), rel2,this.r[0][0] );
            		 }
                	 // 
                //	new FreeTransitionProbs1(
                			
            	}	
                		  
                      
                   
                
                                
                                
//                   new FreeTransitionProbs1(dir1,  states.size());
               if(transProbs[i] instanceof FreeTransitionProbs1){
                   ((FreeTransitionProbs1)transProbs[i]).transitionsOut[0] = null;
               }
               //transProbs[i].validate();          
             //  if(Constants.CHECK && dir1!=null) validate(transProbs[i], rel1.length, i);
            }
        }

}
