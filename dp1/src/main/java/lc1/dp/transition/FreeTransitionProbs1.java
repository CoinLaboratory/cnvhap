/**
 * 
 */
package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.Dirichlet;
import lc1.stats.IntegerDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;


public class FreeTransitionProbs1  extends AbstractTransitionProbs implements Serializable{
    public  PseudoDistribution[] transitionsOut ;
    @Override
    public double[] countsFrom(int i){
    	return transitionsOut[i].counts();
    }
    @Override
    public int noStates() {
       return this.no_states;
    }
    @Override
    public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
     FreeTransitionProbs1 res =  new FreeTransitionProbs1(this, statesToGroup,u);
     return res;
    }
   
    public double[] getAlphaPrior() {
		double[] d = ((SimpleExtendedDistribution)this.transitionsOut[0]).probs;;//.get;
		double[] res = new double[d.length];
		System.arraycopy(d, 0, res, 0, d.length);
	return d;
	}
  
    
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#transferProbsToPseudo()
    
    public void transferProbsToPseudo() {
       for(int i=0; i<transitionsOut.length; i++){
           if(transitionsOut[i]!=null){
               transitionsOut[i].transferProbToPseudo();
           }
       }
    } */
    
    public FreeTransitionProbs1(AbstractTransitionProbs probs){
        no_states = probs.noStates();
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
            double[] probs1 = new double[no_states];
            for(int j1=0; j1<no_states; j1++){
                probs1[j1] = probs.getTransition(j, j1);
            }
            transitionsOut[j] =  new SimpleExtendedDistribution1(probs1, Constants.switchU());
        }
    }
   
    
    public FreeTransitionProbs1(AbstractTransitionProbs probs, double[] u){
      no_states = probs.noStates();
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
            double[] probs1 = new double[no_states];
            for(int j1=0; j1<no_states; j1++){
                probs1[j1] = probs.getTransition(j, j1);
            }
            transitionsOut[j] =  new SimpleExtendedDistribution1(probs1, u[j]);
        }
    }
    
    public FreeTransitionProbs1(FreeTransitionProbs1 probs){
        no_states = probs.noStates();
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
            if(probs.transitionsOut[j]==null) continue;
            else{
                transitionsOut[j] =  probs.transitionsOut[j].clone(Constants.switchU());
            }
        }
    }
    
  final int no_states;
    public FreeTransitionProbs1(FreeTransitionProbs1 probs1, int[] statesToGroup, double[] u) {
     int no_states = statesToGroup.length;
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
        	int j1 = statesToGroup[j];
            if(probs1.transitionsOut[j1]!=null){
                transitionsOut[j] = probs1.transitionsOut[j1].clone(u[j]);
            }
        }
        this.no_states = probs1.no_states;
        this.stateToGroup = statesToGroup;
    }
    
   /* public FreeTransitionProbs1(int noStart, Sampler sampler) {
        int no_states = noStart;
           this.transitionsOut = new PseudoDistribution[no_states];
           for(int j=0; j<no_states; j++){
           	//int j1 = statesToGroup[j];
                   transitionsOut[j] = new SimpleExtendedDistribution(sampler);
           }
           this.no_states = transitionsOut[0].probs().length;
       }*/
    
  /** statesToGroups is which groups, each states should have as a target
   * value of v close to 1 has no effect on skewing distirubiotn
   *  */
    public void modify(int[][] statesToGroups, int[] groups){
    	int no_states = transitionsOut.length;
    	double v = Constants.skewTransitions(0);
    	double v2 = Math.pow(v, 2);
     //   this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=1; j<no_states; j++){
        	int[] states = statesToGroups[j];
        	if(states.length>0){
        	SimpleExtendedDistribution trans = (SimpleExtendedDistribution)this.transitionsOut[j];
			if(trans!=null){
				double[] probs = (trans).probs;
				for(int i=1; i<probs.length; i++){
					probs[i] =v2;
				}
//				probs[groups[j]] = 0.1;  //same group transition
				for(int i=0; i<states.length; i++){
					probs[states[i]] =v;
				}
				Constants.normalise(probs);
			}
			}
		
        	
        }
    }
  public FreeTransitionProbs1 clone(boolean swtch){
     return new FreeTransitionProbs1(this);
  }
   
   /* public FreeTransitionProbs1(FreeTransitionProbs1 tp_init){
   //     this.pseudoCountTrans = tp_init.pseudoCountTrans;
        this.transitionsOut = new SimpleExtendedDistribution1[tp_init.transitionsOut.length];
        for(int i=0; i<transitionsOut.length; i++){
            transitionsOut[i] = 
                tp_init.transitionsOut[i]==null ? null : tp_init.transitionsOut[i].clone();
              
        }
    }
    
   
   /* public FreeTransitionProbs1( Dirichlet samplerFirst, int no_states) throws Exception{
        this( false, samplerFirst, no_states);
    }*/
   public FreeTransitionProbs1(int no_states){
	   this.no_states = no_states;
       this.transitionsOut = new PseudoDistribution[no_states];
   }
    public FreeTransitionProbs1(Sampler samplerFirst){
      no_states = samplerFirst.dist.length;
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
            transitionsOut[j] =  new SimpleExtendedDistribution1(samplerFirst);
        }
    }
    
    public FreeTransitionProbs1(int no_states, double u){
    	double[] d1 = new double[no_states];
    	this.no_states  = no_states;
    	Arrays.fill(d1,1.0/(double)no_states);
        Dirichlet dir = new Dirichlet(d1, u);
          this.transitionsOut = new PseudoDistribution[no_states];
          for(int j=0; j<no_states; j++){
              transitionsOut[j] =  new SimpleExtendedDistribution1(dir);
          }
      }
    
    public FreeTransitionProbs1(int no_states,int no_states_out, double u){
    	double[] d1 = new double[no_states_out];
    	
    	Arrays.fill(d1,1.0/(double)no_states_out);
        Dirichlet dir = new Dirichlet(d1, u);
          this.transitionsOut = new PseudoDistribution[no_states];
          for(int j=0; j<no_states; j++){
              transitionsOut[j] =  new SimpleExtendedDistribution1(dir);
          }
          this.no_states  = no_states_out;
      }
    
    public FreeTransitionProbs1(int no_states,Dirichlet dir){
    
          this.transitionsOut = new PseudoDistribution[no_states];
          for(int j=0; j<no_states; j++){
              transitionsOut[j] =  new SimpleExtendedDistribution1(dir);
          }
          this.no_states  =dir.sample().length;
      }
    
    
    public FreeTransitionProbs1(Sampler samplerFirst,int length){
         no_states = length;
        this.transitionsOut = new PseudoDistribution[no_states];
        for(int j=0; j<no_states; j++){
            transitionsOut[j] =  
                samplerFirst==null ?
                        new SimpleExtendedDistribution(length):       
                new SimpleExtendedDistribution1(samplerFirst);
        }
    }
   
    public FreeTransitionProbs1(   boolean first, Dirichlet samplerFirst, int no_states) throws Exception{
     //   this.states = states;
    	this.no_states = no_states;
        this.transitionsOut = new PseudoDistribution[no_states];
        //this.pseudoCountTrans = new double[no_states];
        int numF = no_states-1;
        Double conc = Constants.initialConcentration(); //0.5;
        if(samplerFirst!=null){
          
            if(!first){
             //   int[] list =  getList(samplerFirst.sample());
                for(int j=1; j<no_states; j++){
                    double[] d1 = new double[samplerFirst.dist.length];
                    Arrays.fill(d1, (1-conc)/ (double)  numF);
                    d1[j]  += conc;
                    d1[0] = 0.0;
                  // double sum = Constants.sum(d1);
                        transitionsOut[j] = 
                           new SimpleExtendedDistribution1(new Dirichlet(d1,samplerFirst.u()));
                }
            }
            else  {
                transitionsOut[0] = new SimpleExtendedDistribution1(samplerFirst);
            }
        }
    }
   
    @Override
    public void addCount(int from, int to, double val){
        transitionsOut[from].addCount(to,val);
     
    }
    
    @Override
    public void addCount(int indexFrom,int from, int to, double val, double dist){
        transitionsOut[indexFrom].addCount(to,val);
     
    }
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#addCounts(lc1.dp.genotype.StateDistribution[], java.util.List)
    
    public void addCounts(StateDistribution[] observed) {
        for(int j=0; j<observed.length; j++){
            if(observed[j]==null) continue;
           // State k = states.get(j);
            if(transitionsOut[j]!=null){
                transitionsOut[j].addCounts(observed[j]);
            }
            else{
                if(Constants.CHECK1 && observed[j].sum()>0) throw new RuntimeException("!!");
            }
        }
        
    } */

  
    
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#getTransition(int, int)
     */
    public double getTransition(int from, int to){
        PseudoDistribution dist = this.transitionsOut[from];
        if(dist!=null){
          /* if(to>=dist.probs.length){
               System.err.println("out of bounds");
           }*/
                return dist.probs(to);
        }
        return 0;
        
    }
    
    public double getTransitionCount(int from, int to){
        PseudoDistribution dist = this.transitionsOut[from];
        if(dist!=null){
          /* if(to>=dist.probs.length){
               System.err.println("out of bounds");
           }*/
                return dist.counts()[to];
        }
        return 0;
        
    }
    
    public double getTransition(int index, int from, int to){
        PseudoDistribution dist = this.transitionsOut[index];
        if(dist!=null){
          /* if(to>=dist.probs.length){
               System.err.println("out of bounds");
           }*/
                return dist.probs(to);
        }
        return 0;
        
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#getTransitionPseudo(int, int)
   
    public double getTransitionPseudo(int from, int to){
        SimpleExtendedDistribution dist = this.transitionsOut[from];
        if(dist!=null){
                return dist.getPseudo(to);
        }
        return 0;
        
    }  */
    
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#initialiseCounts(boolean, boolean)
     */
     public void initialiseCounts(boolean start, boolean end){
          for(int j=0; j<this.transitionsOut.length; j++){
                 if(transitionsOut[j]!=null) transitionsOut[j].initialise();
                 
          }
     }
    
    
 
  
    
   
   
  
     
  
    public String toString(){
         return this.transitionsOut.toString();
     }
  //  final double[] pseudoCountTrans;
     /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#transfer(double)
   
    
    public void transfer(double[] pseudoC, double[] pseudo_exp){
        transfer(pseudoC[0], pseudo_exp[0]);
    }  */
    public void transfer(double pseudoC, double pseudo_exp){
      
         for(int j=0; j<this.transitionsOut.length; j++){
           
             if(transitionsOut[j]!=null && transitionsOut[j].fixedInteger()==null){
           
                 transitionsOut[j].transfer(pseudoC);
               //  PseudoDistribution fixed = ((SimpleExtendedDistribution)transitionsOut[j]).makeFixed();
                // if(fixed!=null) transitionsOut[j] = fixed;
                 if(Constants.CHECK) transitionsOut[j].validate();
             }
         }
     }
     
     /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#getDistributions()
     */
    public Collection getDistributions(){
           return Arrays.asList(transitionsOut);
     }
     
  
    
   

    
  

    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#validate()
     */
    public void validate(){
        for(int j=0; j<this.transitionsOut.length; j++){
            PseudoDistribution dist = transitionsOut[j];
            if(dist==null) continue;
            dist.validate();
          
        }
      //  if(transitionsOut.length>1){
      //  super.validate(this.transitionsOut[1]==null);
       // }
    }
    
    public void validate(boolean check){
        for(int j=0; j<this.transitionsOut.length; j++){
            PseudoDistribution dist = transitionsOut[j];
            if(dist==null) continue;
           
            dist.validate(check);
           
        }
      //  if(transitionsOut.length>1){
      //  super.validate(this.transitionsOut[1]==null);
       // }
    }
    public void validate(int start){
        for(int j=start; j<this.transitionsOut.length; j++){
            PseudoDistribution dist = transitionsOut[j];
            if(dist==null) continue;
            dist.validate();
          
        }
        if(transitionsOut.length>1){
        super.validate(this.transitionsOut[1]==null);
        }
    }




    public int length() {
        return this.transitionsOut.length;
    }



    public void setTransitionScore(int from, int to, double d,  int length) {
        PseudoDistribution dist = transitionsOut[from];//[from];
        if(dist==null){
            dist = new SimpleExtendedDistribution(length);
           // dist.fillPseudo(0);
            
           transitionsOut[from] = dist;
        }
        dist.setProbs(to, d);
       // dist.setPseudo(to,d_ps);
        
        
    }

    public void print(PrintWriter pw, Double[] hittingProb, double dist) {
        String st = transitionsOut[0]==null ? transitionsOut[1].getPrintString() : 
            transitionsOut[0].getPrintString() ;
        pw.println("out");
      
       //
        for(int i=0; i<transitionsOut.length; i++){
          //  Integer fixed = transitionsOut[i].fixedInteger();
          //  if(fixed!=null){
          //      prob
          //  }
            if(transitionsOut[i]==null){
               pw.print("{"+i+"};");
            }
            
            else{
              //  double thresh =  1.0/(transitionsOut.length-1)-0.01;
                transitionsOut[i].printSimple(pw, i+" ", ";" ,0.0);
//                pw.print(transitionsOut[i].probs[i]+"; ");//;print(pw, true, st, ";  ");
            }
        }
        pw.println("\nin");
        Integer fixed = transitionsOut[0].fixedInteger();
        if(fixed!=null){
        	pw.print(fixed);
        }else{
        SimpleExtendedDistribution distr = new SimpleExtendedDistribution
        
        (
        		transitionsOut[0]==null ? transitionsOut[1].counts().length:
        		transitionsOut[0].counts().length);//double[] prob = new double[transitionsOut.length];
        for(int i=0; i<distr.counts.length; i++){
            distr.initialise();
            for(int j=0; j<transitionsOut.length; j++){
               // if(transitionsOut[j]!=null && hittingProb!=null){
              //     distr.counts[j]+= transitionsOut[j].probs(i)*hittingProb[j];
              //  }
               if(transitionsOut[j]!=null){
                    distr.counts[i]+= transitionsOut[j].probs(i);
                }
            }
            if(Constants.sum(distr.counts)!=0){
             distr.transfer(0);
            distr.printSimple(pw, i+" ", ";", 0.0);
            }
        }
        }
        pw.println();
      /*  pw.print(exp_rd.probs[0]);
        pw.print("; ");
        alpha.print(pw, true, alpha.getPrintString(), "; ");
        */
          
      }
	@Override
	public double transfer(double[] pseudoCExp, double[][] d, int pos_index) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double transfer(double d){
		double sum=0;
		 for(int i=0; i<this.transitionsOut.length; i++){
				
			  if(transitionsOut[i]!=null){
		sum+=	((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(d);
			  }
			
		}
		 return sum;
	}
	
	public double logprob(){
		double sum=0;
		 for(int i=0; i<this.transitionsOut.length; i++){
				
			  if(transitionsOut[i]!=null){
		sum+=	((SimpleExtendedDistribution)this.transitionsOut[i]).logProb();
			  }
			
		}
		 return sum;
	}
	
	
	
	int[] stateToGroup;
	
	    public double transferQ(double[] ds,double pseudoAlpha, double pseudoRate,  MatrixExp initial, int i1, double distance, int it) {
		 
	    	initial.setDistance(distance);
		  double sum=0;
		
		 
		  for(int i=0; i<this.transitionsOut.length; i++){
			
			  if(transitionsOut[i]!=null){
				  if( i>0 ){
					  int k = stateToGroup==null ? i-1 : stateToGroup[i]-1;
			       sum+=(this.transitionsOut[i]).evaluate(initial.M.viewRow(k), ds[i1]);
				  }
				  else{
					  sum+= (this.transitionsOut[i]).evaluate(ds[i1]);
				  }
			 
			  }
		  }
		  return sum;
			
		}
	
	public double  transferAlpha(double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
		 // double[][]v = new double[][] {alpha_overall};
		  double sum=0;
		  double[][] alpha = alpha_overall;
		  if(alpha==null){
			 sum+=transfer(pseudoTrans[pos_index]);
		  }
		  else{
		  for(int i=0; i<this.transitionsOut.length; i++){
			
			  if(transitionsOut[i]!=null){
			sum+=((SimpleExtendedDistribution)this.transitionsOut[i]).evaluate(alpha.length==1 ? alpha[0] : alpha[i]);
			  }
			
		}
		  }
		  return sum;
	  }
	//  static double[] zero = new double[] {1e-10,1e-10};
	  public double transferAlpha(double[] pseudoTrans, double[][] ds, int pos_index, int[][] groupToState) {
			if(ds!=null) {
				return this.transferAlpha(pseudoTrans, ds, pos_index);
			}
			else{
				//previous.initialiseCounts(false, false);
				 for(int i=0; i<this.transitionsOut.length; i++){
					 
					  if(transitionsOut[i]!=null){
						  
					  ((SimpleExtendedDistribution)this.transitionsOut[i]).transfer(pseudoTrans[pos_index]);//pseudoTrans[pos_index]);
					  }
					
				}
				//harmonise(groupToState);
				// previous.transferAlpha(pseudoTrans, null, pos_index);
				 
				 
				 
			}
			if(Constants.CHECK) this.validate();
			return 0;
		}
	  /** makes sure the weighted sum of probabilities between super states remains same as pseudo 
	  public void harmonise(int[][] groupToState){
		  double[] d = new double[groupToState.length];
			 for(int i=0; i<groupToState.length; i++){
				 double[] pseudo = transitionsOut[groupToState[i][0]].pseudo();
				 if(groupToState[i].length==1){
					 double[] prob = ((SimpleExtendedDistribution1)transitionsOut[groupToState[i][0]]).probs;
					 if(i==1 && pseudo[2]>0.1){
					 }
					 System.arraycopy(pseudo, 0, prob,0,prob.length);
				 }
				 else{
				
					 int[] gToS = groupToState[i];
					// for(int index =0; index < 2; index++){
					 harmonise(gToS,d, pseudo, 0,i);
					// }
				//	 if(pseudo[2]<0.6){
					//	 System.err.println();
				//	 }
					 
				 }
			 }
	  }*/
	  
	  public void harmonise(int[][] groupToState, double[] hittingProbs){
		  double[] d = new double[groupToState.length];
			 for(int i=0; i<groupToState.length; i++){
				 double[] pseudo = transitionsOut[groupToState[i][0]].pseudo();
				 if(groupToState[i].length==1){
					 double[] prob = ((SimpleExtendedDistribution)transitionsOut[groupToState[i][0]]).probs;
					 System.arraycopy(pseudo, 0, prob,0,prob.length);
				 }
				 else{
				
					 int[] gToS = groupToState[i];
					// for(int index =0; index < 2; index++){
					 harmonise(gToS,d, pseudo, 0, hittingProbs);
					// }
				//	 if(pseudo[2]<0.6){
					//	 System.err.println();
				//	 }
					 
				 }
			 }
	  }
	  public void harmonise(int[] gToS, double[] d, double[] pseudo, int index, int group){
		//  System.err.println("harmonising");
		  Arrays.fill(d,0);
			 double total =0;
			 for(int i1=0; i1<gToS.length; i1++){
				 int state = gToS[i1];
				 double cnt = Constants.sum(transitionsOut[state].counts());
				 double[] prob = transitionsOut[state].probs();
				//System.err.println("before "+i1+"\t\t\t\t"+Constants.print(prob)+"\t"+cnt);
				 for(int j=0; j<prob.length; j++){
					 d[j]+=cnt*prob[j];
				 }
				 total+=cnt;
			 }
			 if(total==0) return;
			 for(int j=0; j<d.length; j++){
				 d[j] = d[j]/total;
			 }
			 //System.err.println(index +" target "+Constants.print(pseudo)+"\n actual "+Constants.print(d));
			 for(int i1=0; i1<gToS.length; i1++){
				 int state = gToS[i1];
				 double[] prob = transitionsOut[state].probs();
				  double sum = 0;
				  for(int j=0; j<prob.length; j++){
					 if(d[j]!=0)  prob[j] = prob[j] *(pseudo[j]/d[j]);
					  sum+=prob[j];
				  }
				Constants.normalise(prob);
				// System.err.println("after "+i1+"\t\t\t\t"+Constants.print(prob));
			 }
	  }
	  
	  public void harmonise(int[] gToS, double[] d, double[] pseudo, int index, double[] hittingProbs){
			//  System.err.println("harmonising");
			  Arrays.fill(d,0);
				 double total =0;
				 for(int i1=0; i1<gToS.length; i1++){
					 int state = gToS[i1];
					 double cnt = hittingProbs[state];//Constants.sum(transitionsOut[state].counts());
					 double[] prob = transitionsOut[state].probs();
					// System.err.println("before "+i1+"\t\t\t\t"+Constants.print(prob)+"\t"+cnt);
					 for(int j=0; j<prob.length; j++){
						 d[j]+=cnt*prob[j];
					 }
					 total+=cnt;
				 }
				 if(total==0) return;
				 for(int j=0; j<d.length; j++){
					 d[j] = d[j]/total;
				 }
				 //System.err.println(index +" target "+Constants.print(pseudo)+"\n actual "+Constants.print(d));
				 for(int i1=0; i1<gToS.length; i1++){
					 int state = gToS[i1];
					 double[] prob = transitionsOut[state].probs();
					  double sum = 0;
					  for(int j=0; j<prob.length; j++){
						 if(d[j]!=0)  prob[j] = prob[j] *(pseudo[j]/d[j]);
						  sum+=prob[j];
					  }
					Constants.normalise(prob);
					// System.err.println("after "+i1+"\t\t\t\t"+Constants.print(prob));
				 }
		  }
	public void setCounts(double d, Sampler sampler) {
		for(int i=0; i<this.transitionsOut.length; i++){
			Double[] d1 = sampler.sample();
			d1[0] = 0.0;
			((SimpleExtendedDistribution)this.transitionsOut[i]).setCounts(d,  d1);
		}
		
	}
	public void setToIdentity() {
		for(int i=0; i<this.transitionsOut.length; i++){
			transitionsOut[i] = new IntegerDistribution(i,null,null);
		}
		
		
	}
  
    
    
  
   
   
       
      
}