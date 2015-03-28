package lc1.dp.core;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.logging.Logger;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CachedHMM;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.MergedEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.IntegerDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.StateDistribution;
import lc1.util.Constants;
import lc1.util.Executor;
import lc1.util.PseudoIterator;

import org.apache.commons.pool.PoolableObjectFactory;

/** need to 1) initialise 2) calculateExpectation 3) do maximisation (ie. reassing parameters ) 
 * for each iteration
 * */
public class BaumWelchTrainer {
    public static ExecutorService es = Constants.numThreads()== 0 ? null:  
       Executor.getEs(BaumWelchTrainer.class, Constants.numThreads());;
        public static void involeTasks(List l, boolean seq) throws Exception{
           if(!seq && Constants.numThreads()>1) es.invokeAll(l);
            
           else{
           for(Iterator it = l.iterator(); it.hasNext();){
               try{
                ((Callable)it.next()).call();
               }catch(Exception exc){
                   exc.printStackTrace();
                   System.exit(0);
               }
            }
           }
        }
      //  
   public MarkovModel hmm;
  // final int no_copies;
   final double[] logprob;
    public EmissionState[] data;
    int firstNonZeroIndex =-1;
    //DP[] dp;
  //  DP dp1;
    final MyObjectPool dp;
    final MyObjectPool stateDist;
   // PseudoIterator weights;
  final  double[] weight; //
    
    final int seqLength;
   public  boolean trainDists = false;
   
   private double weight(int index){
     return 1.0;//Constants.pseudoMod1(index);
     /*  if(true) return 
       if(index <0){
         
         // System.err.println("WARNING index less than 0 "+data[j].getName());
          
       }
      // System.err.println("weight is "+weight[index]+" for "+data[j].getName());
       return 1.0;//weight[index];*/
   }
   final String[] dataNames;
    public BaumWelchTrainer(final MarkovModel hmm, EmissionState[] obj, PseudoIterator weights, String[] dataNames, double[] weight){
        this.weight = weight;
    	this.hmm = hmm;
    	//this.probRDists = (IlluminaProbR[]) probDists[0];
    	//this.probBDists = (IlluminaProbB[]) probDists[1];
        this.dataNames = dataNames;
        this.pcs = new PropertyChangeSupport(this);
        this.data = obj;
        this.dp = new MyObjectPool(new PoolableObjectFactory(){

            public void activateObject(Object arg0) throws Exception {}
            public void destroyObject(Object arg0) throws Exception {}
            public Object makeObject() throws Exception {
                System.err.println("MAKING NEW DP OBJECT ");
                return new DP(hmm,"",Constants.isLogProbs(), hmm.noSnps , false);
            }
            public void passivateObject(Object arg0) throws Exception {}
            public boolean validateObject(Object arg0) {return true;}}   ,Math.max(2, Constants.numThreads()), Constants.numThreads());
        this.stateDist = new MyObjectPool(new PoolableObjectFactory(){

            public void activateObject(Object arg0) throws Exception {}
            public void destroyObject(Object arg0) throws Exception {}
            public Object makeObject() throws Exception {
              //  System.err.println("MAKING NEW DP OBJECT ");
                return  new StateDistribution(hmm.modelLength());
            }
            public void passivateObject(Object arg0) throws Exception {}
            public boolean validateObject(Object arg0) {return true;}}   ,Math.max(2, Constants.numThreads()), Constants.numThreads());
      //  this.emissionCount = new StateDistribution[obj.length];
       
        this.logprob = new double[obj.length];
     
        seqLength = obj[0].length();
     
        for(int j=0; j<obj.length; j++){
            if(!trainDists && ((HaplotypeEmissionState)obj[j]).hasIlluminaDist()){
                trainDists =true;
            }
            int index = ((HaplotypeEmissionState)this.data[j]).dataIndex();
            double weightj = this.weight(index);
            if(weightj>1e-10 || Constants.getMinValue(Constants.r_train(0))<1e7 || Constants.getMinValue(Constants.b_train(0))<1e7){
       
                if(firstNonZeroIndex<0) firstNonZeroIndex = j;
            }
        }
       
           
              /*  for(int i=0; i<seqLength; i++){
                   
                    PseudoDistribution[] transD =((FreeTransitionProbs1)((CachedHMM)hmm).trans.transProbs[i]).transitionsOut;
                    transProbs[i] = transD;
                 
               
                }*/
         //   this.probDists = new ProbDists[Constants.format().length];
            //for(int k=0; k<hmm2.size(); k++){
               /* CompoundMarkovModel cmm =  (CompoundMarkovModel)hmm;
                for(int i=0; i<probDists.length; i++){
                	if(!Constants.format()[i].startsWith("geno")  || (Constants.plot()>0 && Constants.plot(i))){
                    probDists[i] = new ProbDists();
                         //   if(k!=1) throw new RuntimeException("! need to check  for this case - we could be adding the same probR multiple times to the list");
                            cmm.probB(i).probDists(probDists[i].s2, probDists[i].probBIndices);
                            cmm.probR(probDists[i].stateIndices,probDists[i].stateIndexCNV, probDists[i].s1, probDists[i].s1_memb, i);
                    }
                	//}
                  //   m2 = transform(s2);
                    //ss2= new ArrayList<ProbabilityDistribution>();
                   // fill(s2, ss2, )
                   // if(Constants.trainEnsemble()>0){
                            if((Constants.plot(i) && Constants.plot()>0) || !Constants.format()[i].startsWith("geno")){
                        probDists[i].mvf =
                        	new ProbMultivariate(probDists[i].s1, probDists[i].s1_memb);
                        	probDists[i].initColors();
                            }
                            else{
                            	//probDists[i]=null;
                            }
                    //}
                   
                }*/
    }
   
    



    private String print(double[] weight2) {
       Double[] d = new Double[weight2.length];
       for(int i=0; i<d.length; i++){
           d[i] = weight2[i];
       }
       return Arrays.asList(d).toString();
    }
    public double validate(Collection<StateDistribution> dist, int no){
        double sum=0;
       // double[] sum_i;
        for(Iterator<StateDistribution> it = dist.iterator(); it.hasNext();){
            StateDistribution d = it.next();
            if(d==null) continue;
            double s = d.sum();
            sum+=s;
            //logger.info(hmm.getState(i)+" "+s);
        }
      if(Math.abs(sum/(double)no-1.0)>0.01){
          if(Math.abs(sum / (double) no- 1.0)>0.3){
              throw new RuntimeException("sum not right "+sum+" "+no);
          }
          else{
              Logger.global.warning("sum was not 1 "+sum*no+".  Normalising !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
              for(Iterator<StateDistribution> it = dist.iterator(); it.hasNext();){
                  StateDistribution d = it.next();
                  if(d==null) continue;
                  d.multiplyValues((double)no/sum);
                  //logger.info(hmm.getState(i)+" "+s);
              }
          }
      }
      
       return sum;
    }
    
    // [position, from_state_index] 
   //final PseudoDistribution[][] transProbs;
    public PseudoDistribution[] transProbs(int i){
        return ((FreeTransitionProbs1)((CachedHMM)hmm).trans.transProbs[i]).transitionsOut;
    }
    public DP dp(int j) throws Exception{
        return (DP) this.dp.getObj(j);
    }
    public void addTransitionProbSum(int j) throws Exception{
   	 int index = ((HaplotypeEmissionState)this.data[j]).dataIndex();
     double weight = this.weight(index);
    	for(int i=-1; i<seqLength-1; i++){
        
            dp(j).addTransitionPosterior(i, 
            		  transProbs(i+1), weight);
        }
    }
    
    
    public void validateCounts(int total){
        /*for(int i=-1; i<seqLength; i++){
        
            try{ 
            double sum = validate(Arrays.asList(transProbs[i+1]), this.data.length);
            if(Math.abs(sum/total-1.0) > 0.01) throw new RuntimeException("!! "+sum+" "+total);
            }catch(Exception exc){
                exc.printStackTrace();
                System.exit(0);
            }
        }*/
    }
    
   
    
   public void compare(Map<State, SimpleDistribution> dist1, Map<State, SimpleDistribution> dist2){
       for(Iterator<State> it = dist1.keySet().iterator(); it.hasNext();){
           State k= it.next();
           SimpleDistribution d1 = dist1.get(k);
           SimpleDistribution d2 = dist2.get(k);
           if(d1.different(d2)){
               Logger.getAnonymousLogger().info("changed "+k+"-->"+"\n"+d1.toString()+" to "+d2.toString());
           }
       }
   }
  
 // DP backgroundDP=null;
 
  
  
  public Callable expectationStep(final int j1, final double[] pseudo, final boolean last, final boolean updateEmissions){
       Callable run = new Callable(){
           public Object call(){
               try{
                
             //  double time=  System.currentTimeMillis();
         DP  dp1;
       //  if(backgroundDP!=null && j1==0){
       // 	 dp1 =backgroundDP;
      //   }
        // else{
        	 dp1= dp(j1);
        	 dp1.setData(data[j1]);
               dp1.reset(true);
       //  }
      //   System.err.println("doing "+j1);
       logprob[j1] =  dp1.search(true, false);
            //      System.err.println("done "+j1);
       Object[] tmp = new Object[]{j1, dp1};
       firePropertyChange("expec_i", null, tmp);
               if(Double.isInfinite(logprob[j1])) throw new RuntimeException("is infinite!");
               try{
                   addTransitionProbSum(j1);
                 
                 //  if(pseudo[emiss_col]<1000 || pseudo[data_col]<1000){ 
                       addEmissionProbSum(j1, last, updateEmissions);
                //   }
             
               }catch(Exception exc){
                   exc.printStackTrace();
                   System.exit(0);
               }
               dp.returnObj(j1);
               }catch(Exception exc){
                   exc.printStackTrace();
               }
             
               return logprob[j1];
           }
       };
       return run;
     //  run.run();
//       this.es.execute(run);
       //return weight[j1]*dp[j1].logprob;
   }
   
  public static int emiss_col = 0;
  public static int data_col = 3;
   
 // public static Long[] time = new Long[11];
  public static Long[] t = new Long[12];
  
  static String tPrint = "%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %5i";
  public static long t1;
    /** calculates the counts for the hmm 
     *  */
    public double expectationStep(final double[] pseudo, final boolean last, final boolean updateEmissions){
      /*for(int i=0; i<seqLength; i++){
            for(int j=0; j<hmm.modelLength(); j++){
                if(transProbs[i][j]!=null)
                transProbs[i][j].initialise();
            }
      }*/
      
        double logprob =0;//((FreeHaplotypeHMM)((CompoundMarkovModel)this.hmm).getMemberModels()[0]).trans.logLAll;
        Arrays.fill(t, (long)0);//Arrays.fill(time, (long)0);
        List tasks = new ArrayList<Callable>();
       /* if(Constants.modelbg()){
        	try{
        	// this.backgroundDP =  dp(0); //assumes first is bg
            //  backgroundDP.setData(data[0]);
               //   dp[j1] = dp1; 
             //       backgroundDP.reset(true);
        	}catch(Exception exc){
        		exc.printStackTrace();
        	}
        }*/
        for(int j=0 ;j<data.length; j++){
        	  int index = ((HaplotypeEmissionState)this.data[j]).dataIndex();
              if(this.include(index) || updateEmissions){
            tasks.add(this.expectationStep(j,pseudo, last, updateEmissions));
              }
            t1 = System.currentTimeMillis();
       
        }
        try{    
        /*  for(int i=0; i<tasks.size(); i++){
                ((Callable)tasks.get(i)).call();
            }*/
         // es.
            this.involeTasks(tasks,true);
            /*if(this.backgroundDP!=null){
            	backgroundDP.exponentiateEmissions();
            	this.expectationStep(0, pseudo, last, updateEmissions).call();
            }*/
//           List res =  this.es.invokeAll(tasks);
//          this.es.awaitTermination(60, TimeUnit.SECONDS);
            //logprob = this.hmm
            
          for(int j=0; j<data.length; j++){
              logprob+=weight(((HaplotypeEmissionState)data[j]).dataIndex())*this.logprob[j];
          }
          this.firePropertyChange("expectation1", null, this);
          this.firePropertyChange("expectation2", null, this);
        }catch(Exception exc){
            exc.printStackTrace();
        }
       t[3] -=t[2];
       t[2]-=t[1];
       t[1]-=t[0];
        t1 = System.currentTimeMillis();
        if(firstNonZeroIndex<0) return logprob;
        for(int k=0; k<hmm.modelLength(); k++){
            for(int i=-1; i<seqLength-1; i++){
                PseudoDistribution[] dist_i = transProbs(i+1);
              
                PseudoDistribution dist_ik = dist_i[k];
                if(dist_ik==null) continue;
                double[] counst = dist_ik.counts();
                for(int j=0; j<hmm.modelLength(); j++){
                     counst[j]*=hmm.getTransitionScore(k, j, i+hmm.getState(j).adv);
                }
               }
        }
        t[4] = System.currentTimeMillis()-t1;
        if(Constants.CHECK){
            this.validateCounts(data.length );
        }
       // hmm.setStatesChanged(false);
        tasks.clear();
        
            for(int i=-1; i<seqLength-1; i++){
                hmm.addCounts(transProbs(i+1), i, this.data.length); //this adds directly to underlying model if it is a pair markov model
            }
            t[5] = System.currentTimeMillis()-t1-t[4];
            if(hmm instanceof CachedHMM){
                ((CachedHMM)hmm).transferEmissionCountsToMemberStates();
                
            }
           
            t[6] = System.currentTimeMillis()-t1-t[4]-t[5];
           // weight = weights.next();
           // System.err.println("weight is "+print(weight));
        return logprob;
    }
    private boolean include(int index) {
    	//if(true) return true;
    	   return (this.weight[index] > 1e-9 || !Constants.format[index].startsWith("geno"));
    			   //(pseudo[3]>1e3 || (Constants.getMinValue(Constants.r_train(0))>1e7 && Constants.getMinValue(Constants.b_train(0))>1e7))) continue;
           
	}





	public static  boolean maximisationStep(final MarkovModel hmm, final int index, List tasks){
    	Logger.global.info("maximising hmm params");
    	tasks.add(new Callable(){
           public Object call(){
               try{
               hmm.transferCountsToProbs(index);
               }catch(Exception exc){
                   exc.printStackTrace();
               }
               return null;
           }
       });
        t[t.length-1] = System.currentTimeMillis();
      //  System.err.println(Format.sprintf(tPrint, t));
        Logger.global.info("done maximising hmm params");
       return true;
    }
    
   
    final static boolean ML = false;
  
    public void addEmissionProbSum(Integer ll, boolean last, boolean updateEmissions) throws Exception{
        Object[] tmp = new Object[7];
        Arrays.fill(tmp, null);
        tmp[1] = ll;
        HaplotypeEmissionState hes = (HaplotypeEmissionState)data[ll];
        tmp[3] = hes.dataIndex();
        StateDistribution emissionCount = (StateDistribution) this.stateDist.getObj(ll) ;
        this.firePropertyChange("init", null, tmp);
   //   String name = this.data[ll].name;
     /* for(int lk=0; lk<this.data.length; lk++){
    	  if(data[lk].name.equals("471111")) {
    		  System.err.println(lk);
    		  System.exit(0);
    	  }
      }*/
   ///   if(name.equals("471111")){
    //	  System.err.println("h "+ll);
     // }
        for(int i=0; i<seqLength; i++){
            emissionCount.reset();
            tmp[2] = i;
     // if(ML)    dp[l].getPosteriorML(i, emissionCount);
     // else 
           dp(ll).getPosterior(i, emissionCount);
           /*if(ll>0 && this.backgroundDP!=null){
        	   
        	   this.backgroundDP.updateEmissions((HaplotypeEmissionState) data[ll], emissionCount, i);
           }*/
           if(updateEmissions && hes.emissions[i] instanceof IlluminaRDistribution){
        	   EmissionStateSpace emstsp =  hes.getEmissionStateSpace();
        	   double[] prob =  PairEmissionState.pool.getObj(emstsp.genoListSize());
        	 //  if(true) throw new RuntimeException("!!");
        	   
        	   Sampler.getProbOverStates(emissionCount, hmm, hes, i, prob,Constants.isLogProbs());
        	   int m = Constants.getMax(prob);
        	   if(prob[m]>0.99){
        		 
        		  hes.emissions[i] = new IntegerDistribution(m, emstsp);
        	   }
        	   else{
        		  hes.emissions[i]  = new SimpleExtendedDistribution(prob, Double.POSITIVE_INFINITY,emstsp);
        	   }
        	   PairEmissionState.pool.returnObj(prob);
//        	   ((HaplotypeEmissionState)tmp[3]).emissions[i];
           }
          
          tmp[0] = emissionCount;
          tmp[4] = data[ll].getName();
          tmp[5] = data[ll].noCop();
          tmp[6] = data[ll].name;
          this.firePropertyChange("emiss", null, tmp);
          //this.dataIndex(i);
         // for(int ik=0; ik<this.hmmplots.size(); ik++){
          //    hmmplots.get(ik).update(emissionCount[l], 1, i);
          //}
        if(Constants.CHECK) try{
            emissionCount.validate();
        }catch(Exception exc){
            System.err.println("prob at "+ll+" "+i+" of "+this.seqLength);
       //     System.err.println(Arrays.asList(this.data[l].getBestIndices()));
            System.err.println(emissionCount);
            exc.printStackTrace();
            System.exit(0);
        }
  
          for(int j=0; j<hmm.modelLength(); j++){
             // State st = hmm.getState(j);
              if(dp(ll).emiss[j]){
                  EmissionState k = (EmissionState)hmm.getState(j);
                  double val = emissionCount.get(j);
                  if(val>Constants.bwThresh()){
                	  int index = ((HaplotypeEmissionState)this.data[ll]).dataIndex();
                	 // Constants.pseudoMod1(index);
                      k.addCount((HaplotypeEmissionState) data[ll], val,i, weight(index), Constants.isLogProbs());
                 }
              }
          }
         
        }
        stateDist.returnObj(ll);
        this.firePropertyChange("finished", null, tmp);
           
    }
    private double sum(Iterator<SimpleDistribution>it){
        double sum=0;
        while(it.hasNext()){
            sum+= it.next().sum();
        }
        return sum;
    }

  
 //  final List<Integer> copyNumber = new ArrayList<Integer>();
 //  public final  ProbDists[] probDists;
  //  public final IlluminaProbR[] probRDists;
 //   public final IlluminaProbB[] probBDists;
    
  
  //  final List<ProbabilityDistribution> ss2,ss1 ;
  // final Map<String, List<ProbabilityDistribution>> m1,m2;
 
    public static final List<ProbabilityDistribution> extract(Map<String, List<ProbabilityDistribution>> m){
        List<ProbabilityDistribution> res = new ArrayList<ProbabilityDistribution>();
        for(Iterator<List<ProbabilityDistribution>> it = m.values().iterator(); it.hasNext();){
            res.add(it.next().get(0));
        }
        return res;
    }
  /*  public static final Map<String, List<ProbabilityDistribution>> transform(Collection<ProbabilityDistribution> s2){
        Map<String, List<ProbabilityDistribution>> m = new HashMap<String, List<ProbabilityDistribution>>();
        
        for(Iterator<ProbabilityDistribution> it = s2.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
            String st = pdist.toString();
            List<ProbabilityDistribution > l = m.get(st);
            if(l==null){
                m.put(pdist.toString(), l = new ArrayList<ProbabilityDistribution>());
            }
            l.add(pdist);
        }
        return m;
    }
    public final static void train(ProbabilityDistribution[][] m, final double pseudo, List tasks){
        final double[] meanvarskew = Constants.meanvarskewprior();
        for(int i=0; i<m.length; i++){
       //  for(Iterator<ProbabilityDistribution> it = m.iterator(); it.hasNext();){
             final ProbabilityDistribution[] pdist0 =m[i];
            
             tasks.add(new Callable(){
            public Object call(){
                try{
                	for(int j=0; j<pdist0.length; j++){ 
                		if(pdist0[j].sum()>Constants.trainThresh()){
             	//	Logger.global.info("minimizing skew normal "+pdist0);
                 ((ProbabilityDistribution)pdist0[j])
                 .maximise(pseudo*meanvarskew[0], pseudo*meanvarskew[1], pseudo*meanvarskew[2]);
             	Logger.global.info("done minimizing skew normal "+pdist0);
                		}
                	}
                }catch(Exception exc){
                    exc.printStackTrace();
                }
                 return null;
            }
            
            // }
             });
           //  }
         }
        
         
     }*/
/*public final static void train(List<ProbabilityDistribution> m, final double pseudo, List tasks){
       final double[] meanvarskew = Constants.meanvarskewprior();
        for(Iterator<ProbabilityDistribution> it = m.iterator(); it.hasNext();){
            final ProbabilityDistribution pdist0 =it.next();
            if(pdist0.sum()>Constants.trainThresh()){
            tasks.add(new Callable(){
           public Object call(){
               try{
            		Logger.global.info("minimizing skew normal "+pdist0);
                ((ProbabilityDistribution)pdist0)
                .maximise(pseudo*meanvarskew[0], pseudo*meanvarskew[1], pseudo*meanvarskew[2]);
            	Logger.global.info("done minimizing skew normal "+pdist0);
               }catch(Exception exc){
                   exc.printStackTrace();
               }
                return null;
           }
           
           // }
            });
            }
        }
        
    }*/

public IlluminaNoBg getState(int i, int index){
    EmissionState state =data[i];
    if(state instanceof MergedEmissionState){
        return (IlluminaNoBg) ((MergedEmissionState)state).em[index];
     }
     else{
         return  (IlluminaNoBg)state;
     }
}


/* k is individual, i is position 
 * prob_b is prob emission state space indices
 * prob_r is prob over different copy indices
 * */
//double[] condb;

final public PropertyChangeSupport pcs;
//List<HMMPanel> hmmplots = new ArrayList<HMMPanel>();

public void addPropertyChangeListener(PropertyChangeListener arg0){
    this.pcs.addPropertyChangeListener(arg0);
}

public void firePropertyChange(String propertyName, Object oldValue, Object newValue){
    pcs.firePropertyChange(propertyName, oldValue, newValue);
}





public Set<Short> getDataIndices(boolean all) {
   Set<Short> s = new HashSet<Short>();
   for(int i=0; i<this.data.length; i++){
	   if(all){
		   ((HaplotypeEmissionState)data[i]).dataIndices(s);
	   }
	   else s.add((short)((HaplotypeEmissionState)data[i]).dataIndex());
   }
   return s;
}
String[] indexString = null;
public String[] getIndexString() {
	if(indexString!=null) return indexString;
	Set<String>[] res = new Set[this.hmm.noSnps];
	indexString = new String[hmm.noSnps];
	for(int i=0; i<res.length; i++){
		res[i] = new HashSet<String>();
	}
	for(int j=0; j<this.data.length; j++){
		for(int i=0; i<res.length; i++){
			short di = ((HaplotypeEmissionState)data[j]).emissions[i].getDataIndex();
		
			if(di>=0){
				//din.add(di);
				String dn = dataNames[(int)di];
				res[i].add(dn.substring(0, Math.min(dn.length()-1, Constants.indexLength())));
			}
		}
	}
	for(int i=0; i<res.length; i++){
		indexString[i] = res[i].toString().replaceAll("[\\[\\],]", "");
		//indexString[i] = indexString[i];
	}
	return indexString;
}



short[][] data_indices = null;
public Set<Set<Short>> fillDataIndices(){
	//if(data_indices!=null) return;
	data_indices = new short[data.length][];
	Map<Set<Short>, short[]> m = new HashMap<Set<Short>, short[]>();
	for(int i=0; i<data.length; i++){
		Set<Short> res  = new HashSet<Short>();
		((HaplotypeEmissionState)this.data[i]).getDataIndices(res);
		if(m.containsKey(res)){
			data_indices[i] = m.get(res);
		}
		else{
			short[] res1 = new short[res.size()];
			int ii=0;
			for(Iterator<Short> it = res.iterator(); it.hasNext();){
				res1[ii] = it.next();
				ii++;
			}
			m.put(res, res1);
			data_indices[i] = res1;
		}
	}
	return m.keySet();
}
public short[] dataIndex(int l) {
	return data_indices[l];
	
}

/*public void addHMMPlot(HMMPanel hmmP){
    this.hmmplots.add(hmmP);
}*/
     
}

    

