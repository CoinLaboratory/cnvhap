package lc1.dp.states;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.PseudoMixture;
import lc1.stats.SimpleDistribution;
import lc1.stats.SkewNormal;
import lc1.util.Constants;

public abstract class EmissionState extends  State {
    static double round = Constants.round();
   // protected boolean changedParams = true;
    final protected SimpleDistribution lengthDistrib;

    public abstract void append(EmissionState emissionState);

   
   /* public Integer noCop(int data_index){
        return null;
    }*/
    public Integer noCop(){
        return null;
    }
    
  
 


   /* public static IlluminaRArray getProbRGroup(String name, int noCop, double bg, double[][] meanvarskew){
    	 lc1.stats.ProbabilityDistribution[] res = new lc1.stats.ProbabilityDistribution[meanvarskew.length];
    	 for(int i=0; i<meanvarskew.length; i++){
    		 if(meanvarskew[i]!=null){
    			 res[i] = make(name, meanvarskew[i], i);
    		 }
    	 }
    	 return new IlluminaRArray(res,noCop);
     
    }*/
    
   /* (non-Javadoc)
 * @see lc1.dp.states.EmissState#noSnps()
 */
public abstract int noSnps();
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#getEmissionStateSpace()
     */
    public abstract EmissionStateSpace getEmissionStateSpace();
    public  EmissionState(String name, int adv){
        super(name, adv);
        this.lengthDistrib = new SimpleDistribution(new int[] {(int) adv}, new double[] {1});
    }
    
    public synchronized void addCountSynch(int obj_index,  double value, int i){
       this.addCount(obj_index,  value, i);
   }
    
    /*obj_index is the index in the state space (which is stored in the hmm)*/
 /*   public  void addCount(int obj_index,  double value, int i){
        addCount(obj_index, value, i);
      //  if(this.e!=null){
            for(int k=0; k<phenotype.length; k++){
                addCountDT(phenotype[k], k, value, i);
          //  }
        }
    }*/
    abstract public  void addCount(int obj_index,  double value, int i);
    abstract public  void addCountDT(double phenotype, int phen_index,  double value, int i);
    
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#adv(lc1.dp.states.State)
     */
    public SimpleDistribution adv(State st){
        return this.lengthDistrib;
    }
   // public abstract void soften();
    public EmissionState(EmissionState st1){
        this(st1.name, st1.adv);
    }
  
   
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#print(java.io.PrintWriter, java.lang.String, java.util.List)
     */
    public abstract void print(PrintWriter pw, String st, List<Integer>columns);
    
   
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#initialiseCounts()
     */
    public void initialiseCounts(){
        //if(probR() instanceof TrainableSkewNormal) ((TrainableSkewNormal) probR()).initialiseCounts();
        
    }
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#KLDistance(lc1.dp.states.EmissState)
     */
    abstract public int sample(int i);
    /* obj_index is the object index in the hmm state space */
    abstract public  double score(int obj_index,int i);
  abstract public double[] getEmiss(int i);
    
 //   abstract public void setRandom(double emiss, boolean restart);
    public String toString(int i){
        return "";
    }
    public  String toString(){
      return this.getName();
        /*   StringWriter stw = new StringWriter();
            PrintWriter pw  = new PrintWriter(stw);
            this.print(pw, this.getName(), null);
            pw.close();
            return stw.getBuffer().toString();*/
    }

   abstract public int mostLikely(int pos);
   
   /* (non-Javadoc)
 * @see lc1.dp.states.EmissState#transferCountsToProbs(double)
 */
       public  boolean transferCountsToProbs( double pseudo) {
    	   return true;
         //  if(Constants.getMin(Constants.b_train(0))<1e7 || 
        	//	   Constants.getMin(Constants.r_train(0))<1e7) return true;
          // else return false;
       }
   
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#getFixedInteger(int)
 */
public abstract  Integer getFixedInteger(int i);
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#isFixed()
 */
//public abstract boolean isFixed();
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#score(lc1.dp.states.EmissionState, int, int)

public  double score(HaplotypeEmissionState data_state,int i_data) {
    throw new RuntimeException("!!");
    } */

/* (non-Javadoc)
 * @see lc1.dp.states.EmissionState#addCount(lc1.dp.states.EmissionState, java.lang.Double, int, double)
 */

 public   void addCount(HaplotypeEmissionState data_state, Double double1, int i, double weight, boolean log) {
	 if(double1 < Constants.countThresh()) return;
 	EmissionState hmm_state = this;
 	EmissionStateSpace emstsp = hmm_state.getEmissionStateSpace();
 //	int di = data_state.dataIndex(i);
 //	int nocop = hmm_state.noCop(di);
	
 	/*if(!Constants.joint){
 		int nocop = hmm_state.noCop();
 	 
	   if(bg!=null){
		   if(!Constants.joint)
		  data_state.emissions[i].addRCount(bg,nocop,  double1* weight,i);
	   }
	   else{
		  double[] d = DataCollection.datC.calcDist(nocop,i, data_state.emissions[i]);
		  
		  
		   for(int j=0; j<d.length; j++){
			  data_state.emissions[i].addRCount(j,nocop,  double1* weight*d[j],i);
		   }
		  
	   }
 	}*/
   
    
      Integer index = hmm_state.getFixedInteger(i);
       
      if(index!=null){
     	//  if(Constants.joint ){
         	  if(double1 > Constants.countThresh1()){
         	
         		  data_state.emissions[i].addRBCount(emstsp,index,double1,i);   
         	  }
          /*  }
     	  else{
          data_state.emissions[i].addBCount( index, double1,i);
     	  }*/
          hmm_state.addCountSynch(index, weight*double1, i);
      }
      else{
    	   PseudoDistribution distr = data_state.emissions[i];
    	 // double mixP = 1.0;//distr.weight();
    	  if(distr.weight()>0){
//    	  double sumMix = hmm_state.calcDistribution(data_state, i,0);
            double sum = hmm_state.calcDistribution(data_state, i,log);
   //         mixP = sumMix/sum;
           
          
            for(int j=0; j<hmm_state.distribution.length; j++){
                double prob_j =(log ? Math.exp(hmm_state.distribution[j]-sum) : hmm_state.distribution[j] / sum);
                if(prob_j>0){  
                	
                    double val =  prob_j*double1;
                    hmm_state.addCountSynch(j,val*weight, i);
               
                    if(Constants.joint ){
                 	 if(val > Constants.countThresh1()){
                 	
                 		 data_state.emissions[i].addRBCount(emstsp,j,val,i);   
                 	 }
                    }
                    else{
                     data_state.emissions[i].addBCount( j, val,i);
                    }
                }
            }
    	  }
      }
 }

 /**Note, if log, then log(pr) returned */
 public  synchronized double calcDistribution(HaplotypeEmissionState data_state, int i,boolean log){
    EmissionState hmm_state = this;
 	EmissionStateSpace emStSp = 	hmm_state.getEmissionStateSpace();
 	
 	double sum=0;
     
     if(hmm_state.distribution==null){
         hmm_state.distribution = new double[emStSp.size()];
         hmm_state.cn = new int[hmm_state.distribution.length];
         for(int j=0; j<hmm_state.distribution.length; j++){
         	hmm_state.cn[j] = emStSp.getCN(j);
         }
     }
     Arrays.fill(hmm_state.distribution, log ? Double.NEGATIVE_INFINITY: 0.0);
     for(int j=0; j<hmm_state.distribution.length; j++){
             double prob_j =hmm_state.score(j,i);
             if(prob_j>0){
               double ems =  
             	Constants.joint ?
             			 data_state.emissions[i].scoreBR(emStSp,j, i):
             	  data_state.emissions[i].scoreB(j,i);
              double weight = emStSp.getWeight(j);
                hmm_state.distribution[j] =  log ? Math.log(prob_j*weight) + ems : prob_j   * ems *weight;  
                sum+=log ? Math.exp(hmm_state.distribution[j]) :hmm_state.distribution[j] ;
             
             }
             
     }
     if(log && sum==0 ){
    	int maxind =  Constants.getMax(hmm_state.distribution);
    	
    	 for(int j=0; j<hmm_state.distribution.length; j++){
    		 if(j!=maxind) hmm_state.distribution[j]=Double.NEGATIVE_INFINITY;
    	 }
    	 return hmm_state.distribution[maxind];
     }
   
     return log ? Math.log(sum) : sum;
 }

 public  synchronized double calcDistribution(HaplotypeEmissionState data_state, int i, int mixComponent){
	    EmissionState hmm_state = this;
	 	EmissionStateSpace emStSp = //Emiss.getSpaceForNoCopies(data_state.noCop())
	 		hmm_state.getEmissionStateSpace();
	 	double sum=0;
	     
	     if(hmm_state.distribution==null){
	         hmm_state.distribution = new double[emStSp.size()];
	         hmm_state.cn = new int[hmm_state.distribution.length];
	         for(int j=0; j<hmm_state.distribution.length; j++){
	         	hmm_state.cn[j] = emStSp.getCN(j);
	         }
	     }
	   //  int di = data_state.dataIndex(i);
	     Arrays.fill(hmm_state.distribution, 0.0);
	 //    IlluminaProbB[] probBState = hmm_state.probB();
	    // SkewNormal[] probRState =  hmm_state.probR();
	 //    double probr = this.emissions[i_data].score( probRState);
	     for(int j=0; j<hmm_state.distribution.length; j++){
	    //    if(hmm_state.cn[j]==hmm_state.noCop().intValue()){
	             double prob_j =hmm_state.score(j,i);
	             //int j1 = hmm_state.mod(j,di);
	          //   if(j==1) scores.println("in "+prob_j);
	             if(prob_j>0){
	            //		int j1 = hmm_state.mod(j,di);
	               double ems =  
	             	Constants.joint ?
	             			 data_state.emissions[i].scoreBR(emStSp,j, i,mixComponent):
	             	  data_state.emissions[i].scoreB(j,i);
	            //   double scR =  ( Constants.suppressR() ? 1.0 : this.emissions[i].score(hmm_state.probR()));
	              double weight = emStSp.getWeight(j);
	                hmm_state.distribution[j] =  prob_j 
	                    * ems
	                   // *this.emissions[i_data].score(j, probRState)
	                    *weight;  //do we include the weight term???;
	                sum+=hmm_state.distribution[j];
	              //  if(Constants.CHECK && Double.isNaN(sum)) {
	                //	throw new RuntimeException("!!");
	               // }
	             }
	       //      if(j==1) scores.println("in "+distribution[j]);
	         //}
	     }
	   
	     return sum;
	 }






public   int getBestIndex(HaplotypeEmissionState data_state,  int i, boolean sample, boolean log){
     EmissionState hmm_state = this;
 	double sum = hmm_state.calcDistribution(data_state, i,log);
   //  if(state_indices==null) return this.best_index[i];
    
         int res =  (sample & !log) ? Constants.sample(hmm_state.distribution, sum) :
             Constants.getMax(hmm_state.distribution);
       //  if(res==0){
       //      System.err.println("res ==0");
       //  }
         return res;
    
 }
public   double score(HaplotypeEmissionState data_state , int i_hmm,boolean log) {
	  EmissionState hmm_state = this;
	   double sc =hmm_state.calcDistribution(data_state, i_hmm,log);
	 
    return sc;
}
  



public final  int getBestIndex(int i) {
    Integer index1 = getFixedInteger(i);
    if(index1!=null) return index1;
    else return Constants.getMax(this.getEmiss(i));
}

/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#setFixedIndex(int, int)
 */
public void setFixedIndex(int i, int k) {
   throw new RuntimeException("wrong class "+this.getClass());
 }


//taken from LikelihoodData


public  double[] distribution ;
private int[] cn;

/*{
    if(true) 
   //  if(state_indices==null) return this.best_index[i];
  Integer index = data_state.getFixedInteger(i_hmm);
  if(index!=null) return index;
  index = this.getFixedInteger(i_hmm);
  if(index!=null) return index;
  
  else{
    double sum = calcDistribution(data_state, i_hmm);
    try{ 
        int res =  sample ? Constants.sample(distribution, sum) :
            Constants.getMax(distribution);
        return res;
    }catch(Exception exc){
        exc.printStackTrace();
        return 0;
    }
  }
}

protected double calcDistribution(EmissionState hmm_state, int i_data){
    double sum=0;
    if(true) throw new RuntimeException("!!");
        if(this.distribution==null){
            distribution = new double[this.getEmissionStateSpace().size()];
            cn = new int[distribution.length];
            for(int j=0; j<distribution.length; j++){
            	cn[j] = this.getEmissionStateSpace().getCN(j);
            }
        }
        Arrays.fill(distribution,  0.0);
        
        for(int j=0; j<distribution.length; j++){
            double prob_j =this.score(j,i_data);
            if(prob_j>0){
               distribution[j] =  ( prob_j * hmm_state.score(j,i_data));
               sum+=distribution[j];
            }
        }
    return sum;
}
*/



/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#length()
 */
public final int length() {
  return this.noSnps();//this.getEmissionStateSpace().size();
}

//Assume this emission state is the 'data' and state is the state from the hmm
// both may be fixed, i.e. concentrated on one point
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#score(lc1.dp.states.EmissionState, boolean)
 */


public final  double[] score(HaplotypeEmissionState data_state, boolean logspace, boolean isLog){
    double[] score = new double[data_state.noSnps()];
    EmissionState hmm_state = this;
    if(data_state.getName().equals("107665_19") && hmm_state.getName().indexOf("1")>=0 || hmm_state.getName().indexOf("2")>=0){
		Arrays.fill(score, logspace ? Double.NEGATIVE_INFINITY : 0);  
	  }
    else if(data_state.getName().equals("107665_9") && hmm_state.getName().indexOf("3")>=0 || hmm_state.getName().indexOf("4")>=0){
		Arrays.fill(score, logspace ? Double.NEGATIVE_INFINITY : 0);  
	  }else{
   // int len = this.getEmissionStateSpace().size();
   // EmissionState indi = this; 
    for(int i=0; i<score.length; i++){
           double  sc =hmm_state.score(data_state, i,isLog); //switch
        score[i] = logspace ? Math.log(sc) : sc ;
      
        if(Constants.CHECK && Double.isNaN(score[i])){
            throw new RuntimeException("!!" +sc+" "+logspace);
     
        }
    }
	  }
   
    return score;
}
public Double[] phenValue(){
    throw new RuntimeException("!!");
}
public int dataIndex(int i){
    return -1;
}
// state is hmm state



public void addCountSynchDT(Double[] phenValue, double d, int i) {
  for(int k=0; k<phenValue.length; k++){
      if(phenValue[k]!=null){
      this.addCountDT(phenValue[k], k, d, i);
      } 
  }
    
}
/*public static HaplotypeEmissionState getEmissionState(PhasedDataState obj, boolean replaceNWithUnknown, double pseudo,int[] se){
    EmissionStateSpace emStSp = Emiss.getEmissionStateSpace(obj.noCopies()-1);
    if(emStSp.copyNumber.size()==1) {
        replaceNWithUnknown = true;
    }
    HaplotypeEmissionState hes =  getEmissionState(obj, emStSp, replaceNWithUnknown, se);
   
    if(pseudo >0){
        Map<Integer, List<Integer>> diffCNV = new  HashMap<Integer, List<Integer>>();
        for(int j=0; j<emStSp.defaultList.size(); j++){
            int cn = emStSp.getCN(j);
            if(cn==2) continue;
            List<Integer> l = diffCNV.get(cn);
            if(l==null){
                diffCNV.put(cn, l = new ArrayList<Integer>());
                
            }
            l.add(j);
        }
        double[] incr = new double[emStSp.defaultList.size()];
        for(int j=0; j<incr.length; j++){
            int cn = emStSp.getCN(j);
            if(cn==2) continue;
            List<Integer> l = diffCNV.get(cn);
            incr[j] = (pseudo / (double) diffCNV.size())/(double) l.size();
        }
        if(se!=null){
        for(int i=se[0]; i<se[1]; i++){
            SimpleExtendedDistribution dist =( SimpleExtendedDistribution) hes.emissions[i];
            double[] probs = dist.probs;
          
           // double incr = pseudo / (double)probs.length;
            for(int j=0; j<probs.length; j++){
                probs[j] =  (probs[j]*(1-pseudo)+ incr[j]);
//                int cn = hes.getEmissionStateSpace().getCN(j);
            }
        }
        }
       
    }
    return hes;
}
public static HaplotypeEmissionState getEmissionState(PhasedDataState obj, EmissionStateSpace emStSp, boolean replaceNWithUnknown, int[] se){
   // if(!Constants.trainWithGenotypes()){
    boolean fixed =false;
    Integer n_index = emStSp.getFromString("");
  //  Integer[] fixedI = new Integer[obj.length()];
        HaplotypeEmissionState hes = new HaplotypeEmissionState(obj.getName(), obj.length(), emStSp.size(), emStSp,null, null, obj.data_index);
        for(int i=0; i<obj.length(); i++){
            fixed = se==null || !(i>=se[0] && i<se[1]);
            double[] prob = new double[emStSp.size()];
            Arrays.fill(prob,  0.0);
            Comparable comp = obj.getElement(i);
           Integer genoIndex = emStSp.getGenotype(comp);
            int[] indices;
            double[] weights;
            if(replaceNWithUnknown && (genoIndex==null ||(n_index!=null && genoIndex ==n_index))){
                indices = emStSp.getGenoForCopyNo(obj.noCopies());
                weights = new double[indices.length];
                Arrays.fill(weights,  1.0 / (double)indices.length);
            }
            else{
                 indices = emStSp.getGenotypeConsistent( genoIndex);
                weights = emStSp.getWeights(genoIndex);
            }
                if(fixed && indices.length>1) fixed = false;
                for(int k=0; k<indices.length; k++){
                     prob[indices[k]] =  weights[k];
                }
                if(indices.length>1 || !fixed){
                    hes.setTheta(prob, i);
                }
                else{
                    hes.setTheta(indices[0], i);
                }
                hes.emissions[i].setDataIndex(obj.emissions[i].getDataIndex());
        }
            return hes;
 }*/

/** converts to a bigger state space */
public static EmissionState getEmissionState(EmissionStateSpace emStSp, 
        EmissionState st, EmissionStateSpaceTranslation trans
        ) {
    return new WrappedEmissionState(st, emStSp, trans);
}
/*public static HaplotypeEmissionState getEmissionState(String name, int len, 
        EmissionStateSpace emStSp, EmissionStateSpace subSp, double cn_ratio){
    
    List<String> genoList = new ArrayList<String>();
    for(int k=0; k<subSp.getGenotypeList().size(); k++){
        char[] ch =  ((ComparableArray)subSp.getGenotype(k)).toStringShort().replaceAll("_", "").toCharArray();
        Arrays.sort(ch);
         genoList.add(new String(ch));
    }
        HaplotypeEmissionState hes = new HaplotypeEmissionState(name, len, emStSp.size(), emStSp, null, null);
     
        for(int i=0; i<len; i++){
            double[] prob = new  double[emStSp.size()];
            Arrays.fill(prob,  0.0);
            for(int k=0; k<genoList.size(); k++){
                String compa = genoList.get(k);
              
                int genoIndex = emStSp.getFromString(compa);
                int cn= compa.length();
                int[] indices = emStSp.getGenotypeConsistent( genoIndex);
                double[] weights = emStSp.getWeights(genoIndex);
                double diff = cn-((CompoundEmissionStateSpace)emStSp).noCopies();
                double val =  diff==0 ?cn_ratio :1;
                        //10 : (cn==1 ? 1 : (cn==0 ? 0.1 : 0.01));
                for(int k1=0; k1<indices.length; k1++){
                     prob[indices[k1]] =  (weights[k1]*val);
                }
            }
            SimpleExtendedDistribution.normalise(prob);
            hes.setTheta(prob, i);
        }
       // hes.train_j = false;
       
            return hes;
   
 }*/
public static HaplotypeEmissionState getEmissionState(SimpleScorableObject obj, EmissionStateSpace emStSp, double soften){
    // if(!Constants.trainWithGenotypes()){
  //   boolean fixed = true;
   //  Integer[] fixedI = new Integer[obj.length()];
         HaplotypeEmissionState hes = new HaplotypeEmissionState(obj.getName(), obj.length(), emStSp.size(), emStSp, null, null);
         double incr = soften / (double) emStSp.defaultList.size();
         for(int i=0; i<obj.length(); i++){
             double[] prob = emStSp.getArray((String)obj.getElement(i));
             if(soften>0){
                 double[] cp = new double[prob.length];
                 System.arraycopy(prob, 0, cp, 0, prob.length);
                 for(int k=0; k<prob.length; k++){
                     prob[k] = (cp[k] *(1-soften) +incr);
                 }
             }
            // if(fixed) fixed = prob[Constants.getMax(prob)] > 0.999;
           //  SimpleExtendedDistribution.normalise(prob);
             hes.setTheta(prob, i);
         }
       //  if(fixed){
        //     return new FixedHaplotypeEmissionState(hes);
       //  }
       //  else {
             return hes;
       //  }
     //}
     //else  {
      //   return new FixedHaplotypeEmissionState(obj, emStSp);
     //}
// new AffyEmissionState("", obj, hmm.size(), hmm));
  }

public static EmissionState getEmissionState(EmissionState  mother, EmissionState  child, EmissionStateSpace emStSp){
   // super(mother.name+";"+child.name,emStSp.size());
   return
        new CachedEmissionState(
        new PairEmissionState(Arrays.asList(new EmissionState[] { mother, child}), 
           (CompoundEmissionStateSpace) emStSp, true), emStSp.size());
}


public static EmissionState getEmissionState(EmissionState father, EmissionState mother, EmissionState  child, EmissionStateSpace emStSp){
   // super(father.name+";"+mother.name+";"+child.name,emStSp.size());
    return
        new CachedEmissionState(
        new PairEmissionState(Arrays.asList(new EmissionState[] {father, mother, child}), 
            (CompoundEmissionStateSpace)emStSp, true), emStSp.size());
}


public static HaplotypeEmissionState getEmissionState( EmissionState  mother, EmissionState  child){
    return (HaplotypeEmissionState) getEmissionState( mother, child,Emiss.getEmissionStateSpace(2) );
}
public static EmissionState getEmissionState(EmissionState father,EmissionState  mother, EmissionState child){
    return getEmissionState(father, mother, child,Emiss.getEmissionStateSpace(5) );
}

public EmissionState[] split() {
   return ((CompoundState)this).getMemberStates(true);
    /*String[] name = this.name.split(";");
    LikelihoodData[] res = new LikelihoodData[name.length];
    for(int i=0; i<res.length; i++){
        res[i] = new LikelihoodData(name[i], states[i]);
    }
    return res;
*/  }

public static EmissionState getEmissionState(String name, EmissionStateSpace emStSp, int noSnps){
  return null;
}

public PIGData getGenotypeData() {
    EmissionStateSpace stSp = this.getEmissionStateSpace();
    PIGData res = SimpleScorableObject.make(name, this.noSnps(), stSp, (short) this.dataIndex(0));
    for(int i=0; i<this.noSnps(); i++){
       Comparable em = stSp.get(this.getBestIndex(i));
        res.addPoint(i, em);
    }
    return res;
}
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#noCopies()
 
public int noCopies() {
	return Constants.backgroundCount();
//   return this.getEmissionStateSpace().noCopies();
}*/
/*public EmissionStateSpace getEmissionStateSpaceAtPos(int i) {
   return getEmissionStateSpace();
}*/

public Collection<Integer> getConstantPos() {
   Set<Integer>s = new TreeSet<Integer>();
   for(int i=0; i<this.noSnps(); i++){
       if(this.score(this.getBestIndex(i), i)>0.9999) {
    	   s.add(i);
       }
   }
   return s;
}
public Collection<? extends Integer> getMultPos() {
    Set<Integer>s = new HashSet<Integer>();
    for(int i=0; i<this.noSnps(); i++){
        if(this.getFixedInteger(i)!=null) continue;
        int[] ind = Constants.getMax2(this.getEmiss(i));
        double sum = this.score(ind[0], i) + this.score(ind[1], i);
        if(sum < 0.99999) s.add(i);
//        if(this.score(this.getBestIndex(i), i)>0.9999) s.add(i);
    }
    return s;
}
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#reverse()
 */
public abstract void reverse();
/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#getUnderlyingData(int)
 */
public String getUnderlyingData(int i) {
    return this.getBestIndex(i)+"";
}

/* (non-Javadoc)
 * @see lc1.dp.states.EmissState#getParamIndex()
 */
public abstract int getParamIndex();
public  void switchAlleles(int i) {
    throw new RuntimeException("!! "+this.getClass());
}
public void fixLikelihood(boolean X, int i){
    throw new RuntimeException("!!");
}

/*
public  void fillLikelihood(Locreader mid, List<Integer> loc) {
  
      //boolean X = Constants.chrom[0].equals("X");
    for(int i=0; i<this.length(); i++){
        Location li = new Location(Constants.chrom0(), loc.get(i), loc.get(i));
        Integer fixed = getFixedInteger(i);
        if(fixed!=null) {
            //emissions[i] =new IntegerDistribution(fixed);
        }
        else{
          
         
            Location overl = mid.overlaps(li, 0);
            if(overl!=null){
            }
            else{
              this.fixLikelihood(false, i);
            }
        }
    }
    
}*/
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#removeAll(java.util.List)
     */
    public void removeAll(List<Integer> toDrop) {
        throw new RuntimeException("!!");
      }
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#setAsMissing(java.util.List, double)
     */
    public void setAsMissing(List<Integer> toDrop, double cn_ratio) {
       throw new RuntimeException("!!");
        
    }
    /* (non-Javadoc)
     * @see lc1.dp.states.EmissState#applyAlias(int[])
     */
    public void applyAlias(int[] alias) {
       throw new RuntimeException("!!");
        
    }
    public void setProbR(SkewNormal probR1) {
      throw new RuntimeException("!!");
        
    }
    public double scoreEmiss(Double[] phenValue, int i) {
       throw new RuntimeException("!!");
    }





	public PseudoDistribution emissions(int i) {
		// TODO Auto-generated method stub
		return null;
	}





	public ProbabilityDistribution[][] emissionsDatatype() {
		// TODO Auto-generated method stub
		return null;
	}





	public void refreshSiteEmissions() {
		// TODO Auto-generated method stub
		
	}





	




	





	
	/*public final  int mod(int j, int di){
		return j;
	}*/
	
	
	
    
  
 


}
