package lc1.dp.core;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.logging.StreamHandler;

import lc1.dp.core.StatePath.StatePosEmission;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.BackgroundEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.State;
import lc1.stats.ChiSq;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;
import lc1.stats.StateDistribution;
import lc1.util.Constants;




public class DP{

    
   
    public static boolean  verbose = true;
	
    protected TraceMatrix backwardTrace;
    public final double[][] emissions; //just the emissions in the current position, in order to save space;
    protected TraceMatrix forwardTrace;
    protected MarkovModel hmm;
    final protected boolean logspace;
    final protected int modelLength;
    EmissionState obj;  //this is the data
     public String protName;
    final protected int seqLength;
	
    final double nullEms;//
  // final Double[][] hittingProb;
      
  
        
        double threshold = -200;
  double thresholdMax = 200;
  double thresholdMax1 = -10;
  static Logger logger = Logger.getLogger("lc1.dp.core.DP");
    {
            logger.setLevel(Level.OFF);
            File out = new File("tmp_out");
            try{
           OutputStream os = new FileOutputStream(out);
           logger.addHandler(new StreamHandler(os, new SimpleFormatter()));
         //  logger.removeHandler(handler)
            }catch(Exception exc){exc.printStackTrace();} 
            /*logger.getParent().getHandlers()[0].setFormatter( new Formatter(){
                public String format(LogRecord record){
                    return "DP/"+record.getSourceMethodName()+":\n"+record.getMessage()+"\n";
                }
            });*/
    }
    
    public DP(MarkovModel hmm, String protName, boolean logspace, int seqLength, boolean complete){
        this.logspace = logspace;
        this.modelLength = hmm.modelLength();
        this.seqLength = seqLength;
        this.protName = protName;
        this.hmm = hmm;
        this.nullEms = logspace ? 0 : 1;
        emissions = new double[modelLength][];
        this.emiss = new boolean[modelLength];
        this.adv = new int[modelLength];
    //    this.nonZero = new boolean[seqLength];
        for(int i=0; i<modelLength; i++){
           State st =  hmm.getState(i);
            emiss[i] =  st instanceof EmissionState;
            adv[i] = st.adv;
        }
        Class clazz = complete ? ComplexTerm.class : Term.class;
        try{
            forwardTrace = new TraceMatrix(modelLength, seqLength, true, logspace, clazz);
            backwardTrace = new TraceMatrix(modelLength, seqLength, false, logspace, clazz);
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
        }
        this.stateIndex = new int[hmm.modelLength()];
        Arrays.fill(stateIndex, -1);
    }
   
   public  DP(MarkovModel hmm, String protName, EmissionState obj, boolean logspace, boolean complete){
       this(hmm, protName,  logspace, obj.noSnps(), complete);
       if(Constants.CHECK){
       //    for(int i=0; i<obj.length(); i++){
           int s1 = ((EmissionState)hmm.getState(1)).getEmissionStateSpace().size();
           int s2 = obj.getEmissionStateSpace().size();
           if(s1!=s2) throw new RuntimeException("!! "+s1+" "+s2);
          // }
       }
    this.obj = obj;
    
    //this.hittingProb = hmm.getHittingProb(this.seqLength);
}
   
   
   public void setData(EmissionState emissionState) {
       this.obj = emissionState;
       this.protName = obj.getName();
      this.hmm.allowTransitions(!Constants.parentObjContains(emissionState.getName()));
    if(DataCollection.datC.dc!=null) DataCollection.datC.dc.setIndiv(protName);
    }
   int objIndex = -1; //index keeps track of which param set we are using
   final int[] stateIndex; //index keeps track of param set
 //final boolean[] nonZero;
   double[] cnByState;
   double[][][] genoCount; //for hwe calc
   public final void updateEmissions(HaplotypeEmissionState data, StateDistribution emissionCount, int i){
	   if(cnByState==null){
		   cnByState = new double[data.getEmissionStateSpace().copyNumber.size()];
		   genoCount = new double[this.seqLength][modelLength][cnByState.length];
		   for(int ii=0; ii<seqLength;ii++){
			   for(int j=0; j<modelLength; j++){
			   Arrays.fill(genoCount[ii][j], 0.0);
			   }
		   }
	   }
	   BackgroundEmissionState bg = (BackgroundEmissionState) obj;
	   Arrays.fill(cnByState, 0.0);
	  // int data_index = data.dataIndex();
	   for(int k=1; k<emissionCount.dist.length; k++){
		   int nocop = ((EmissionState)hmm.getState(k)).noCop();
		 //  int nocop = ((EmissionState)hmm.getState(k)).noCop(data_index);
		   cnByState[nocop]+=emissionCount.dist[k];
	   }
	   for(Iterator<int[]> it = hmm.equivalenceClasses(); it.hasNext();){
           int[] equiv = it.next();
           EmissionState state = (EmissionState) hmm.getState(equiv[0]);
          // if(){
               //Logger.global.info("updating "+state.getName()+" "+obj.getName());
               double d  =   (bg).score(state, true,data,cnByState,i);
               for(int k=0; k<equiv.length; k++){
                   //if(logspace) 
                	   emissions[equiv[k]][i]  += d;
                	   for(int kk=0; kk<bg.cn_distribution.length; kk++){
                		   genoCount[i][equiv[k]][kk]+=bg.cn_distribution[kk];
                	   }
                   //else emissions[equiv[k]][i]  *= d;
               }
              
           //}
         //  stateIndex[equiv[0]] = state.getParamIndex();
       }
	  // if(logspace)
     }
   public void exponentiateEmissions(){
	   int bgtype = Constants.bgtype();
	   if(bgtype >0){
		   for(int i=0; i<seqLength; i++){
			   for(int j=0; j<emissions.length; j++){
				   if(emiss[j]){
				   emissions[j][i] =
					   bgtype== 1?
					   //emissions[j][i]+
					  hwep(this.genoCount[i][j],j) : emissions[j][i]+hwep(this.genoCount[i][j],j);
				   }
			   }
		   }
	   }
	   if(!logspace){
		  
		   for(int i=0; i<this.seqLength; i++){
			 //  double min = Double.POSITIVE_INFINITY;
			   double max = Double.NEGATIVE_INFINITY;
			   for(int j=0; j<emissions.length; j++){
				   if(this.emiss[j]){
				//	   if(emissions[j][i]<min) min = emissions[j][i];
					   if(emissions[j][i]>max) max = emissions[j][i];
				   }
			   }
			   for(int j=0; j<emissions.length; j++){
				   if(this.emiss[j]){
					   emissions[j][i] = Math.exp(emissions[j][i] - max);
				   }
			   }
		   }
	   }
	   
   }
   private double hwep(double[] ds,int j) {
	   double sum = Constants.sum(ds);
	   if(sum==0) return Double.NEGATIVE_INFINITY;
      double p = (ds[0]+ds[1]/2.0)/sum;
      double q = (ds[2]+ds[1]/2.0+ds[3]/2.0)/sum;
      double r = (ds[3]/2.0+ds[4])/sum;
     double[] exp = {p*p, 2*p*q,q*q+2*p*r, 2*q*r, r*r}; //0,1,2,3,4 
     
      double obs_exp = 0;
      for(int i=0; i<exp.length; i++){
    	  if(exp[i]>0){
    		  obs_exp+=Math.pow(ds[i] - exp[i],2)/exp[i];
    	  }
      }
      double hwep = ChiSq.chi2prob1(2, obs_exp);
    //  System.err.println(hwep+" "+j);
      if(hwep==0) return Double.NEGATIVE_INFINITY;
    //  if(hwep== 0 || Double.isNaN(hwep)) throw new RuntimeException("is NAN "+hwep+" "+obs_exp);
      else return Math.log(hwep);
}

public final void  fillEmissions(boolean first ) {
      // if(hmm.trainEmissions() || first){
       boolean objUpdated = obj.getParamIndex()!=this.objIndex;
       //double max=Double.NEGATIVE_INFINITY;
       for(Iterator<int[]> it = hmm.equivalenceClasses(); it.hasNext();){
           int[] equiv = it.next();
           EmissionState state = (EmissionState) hmm.getState(equiv[0]);
           
           if(objUpdated || state.getParamIndex()!=stateIndex[equiv[0]] || first){
               //Logger.global.info("updating "+state.getName()+" "+obj.getName());
        	   boolean isLog = Constants.isLogProbs();
               double[] d  =   state.score((HaplotypeEmissionState)obj, logspace &!isLog,isLog); ///swap!!!
               for(int k=0; k<equiv.length; k++){
                   emissions[equiv[k]]  = d;
               }
              
           }
           stateIndex[equiv[0]] = state.getParamIndex();
       }
      
      /* if(((HaplotypeEmissionState)obj).dataIndex()==1){
    	   System.err.println("h");
       }*/
       this.objIndex = obj.getParamIndex();
      // double avg = 1.0/(double) modelLength;
      if(first){
    	//   Arrays.fill(nonZero, true);
    	   outer: for(int i=0; i<this.seqLength; i++){
    		 //  if(emiss[i]){
    		   	double sum=0;
    		   	double cnt=0;
	    		   for(int j=0; j<this.modelLength; j++){
	    			   if(emiss[j]){
	    				  sum+=this.emissions[j][i];
	    				  cnt++;
	    			   }
	    		   }
	    		     if(sum==0){ 
	    		    	 int loc =DataCollection.datC.loc.get(i); 
	    		    	 if(true) throw new RuntimeException("!! "+obj.getName()+" "+i+" "+loc);
			    		    
	    		    	((EmissionState) hmm.getState(2)).calcDistribution((HaplotypeEmissionState)obj, 2, Constants.isLogProbs());
	    		       	   CompoundEmissionStateSpace emStSp =  
		    		    		   Emiss.getSpaceForNoCopies(obj.noCop());
//		    		    		   (CompoundEmissionStateSpace) ((HaplotypeEmissionState)this.obj).getEmissionStateSpace();
		    		      //     double[] init = new double[emStSp.size()];
		    		          
		    		        //   Arrays.fill(init, 1.0/(double) init.length);
		    		   
		    		       //                                          SimpleExtendedDistribution simex = new SimpleExtendedDistribution( 
		    		    //		init, Double.POSITIVE_INFINITY);
		    		//  simex.emstsp = emStSp;
		    		  ((HaplotypeEmissionState)this.obj).emissions[i] = emStSp.getHWEDist1(null);
		//    		    		((HaplotypeEmissionState)this.obj).emissions[i]);
		    		    	// System.err.println("IS ALL ZERO FOR "+this.obj.getName()+" "+i);
		    		   
		    		  for(int j=0; j<this.modelLength; j++){
		    			   if(emiss[j]){
			    			 
			    			   double sc = ((EmissionState) hmm.getState(j)).score((HaplotypeEmissionState) obj, i, Constants.isLogProbs()); ///switch!!
			    			   emissions[j][i] = logspace ? Math.log(sc) : sc;
	;	    			   }
	    		   }
	    		     }
	    		  
	    	//	   nonZero[i] = false;
	    		// }
    	   }
       }
   }
   
 
   
   public void reset(boolean first){
     
     forwardTrace.clear();
     backwardTrace.clear();
     fillEmissions( first);
   this.statePath = null;
    this.logprob = null;
     
   }
   
   
   
  
   
  
 
  
  
  /** l1 is the l1'th state which goes from k*/
  
  
   final int[] adv;
   final boolean[] emiss; //is this emission state
   
  /*iPlus1 
   * k is index of from state */
double backward(int k, int i){
	int max_j = -1;
    int max_i = 0;
	double max = 0;
	double sum = 0;
    double summ=0;
    int[] out = hmm.statesOut(k, i);
    for(int k1=0; k1<out.length; k1++){
        int j1 = out[k1];
        int i1 = i+adv[j1];
        double val = hmm.getTransitionScore1(k, j1, i1);
        summ+=val;
		
		double score;
        AbstractTerm AbstractTerm = backwardTrace.getTrace(j1,i1);
        if(AbstractTerm==null) score=0;
        else{
            double d = emiss[j1] ? emissions[j1][i1] : nullEms;
            score =  val * AbstractTerm.score()*d;
            if(Constants.CHECK && (score > 1e210 || Double.isNaN(score) || Double.isInfinite(score))){
            	throw new RuntimeException(" is nan "+score+" "+k+" "+i+" "+val+" "+d+" "+AbstractTerm.score());
            }
        }
		if(score>max){
			max = score;
			max_j = j1;
            max_i = i1;
		}
      
		sum+=score;
		
	}
    if(Constants.CHECK && i+1 != seqLength && Math.abs(summ-1.0)>SimpleDistribution.tolerance){// && this.hittingProb[k.getIndex()][i]>0.001){
      try{
          hmm.validate(hmm.noSnps);
      }catch(Exception exc){
          exc.printStackTrace();
      }
        throw new RuntimeException("sum not right "+summ+" "+i+" "+this.obj.length()+" "+k+" ");//+this.hittingProb[k.getIndex()][i]);
    }
	if(i>=0){
    	backwardTrace.setTrace(k,i,max_j, max_i, sum);
	}
    if(Constants.CHECK && Double.isNaN(sum) || Double.isInfinite(sum)) throw new RuntimeException(" is nan "+k+" "+i);
	return sum;
}
protected void calcScoresBackward(){
    if(logspace) throw new RuntimeException("only run backward in prob space");
	//initialise
    backwardTrace.setTrace(hmm.MAGIC.getIndex(),seqLength-1, 0,seqLength-1,1);
	//backward
	for(int i=seqLength-1; i>=0; i--){
        double logscale = backwardTrace.getLogScale(i+1);
        for(int j=this.modelLength-1; j>=1; j--){
           // State j = it.next();
            //if(j==hmm.MAGIC
                   // || this.hittingProb[j.getIndex()][i]==0  
              //      ) continue;
			this.backward(j, i);
		}
		double[] min = backwardTrace.minScore(i);
    //    System.err.println("minmaxb "+min[0]+" "+min[1]+" "+logscale);
        if(Constants.CHECK && Double.isInfinite(min[1])){
            
            throw new RuntimeException("max should be greater than zero at "+i+" "+this.obj.getName());
         //           Arrays.asList(d) +"\n"+);
        }
		if(min[0] < threshold){
			double scale = Math.min(-min[0], thresholdMax-min[1]);
			backwardTrace.scale(scale, i);
            logscale+=scale;
			
		}
        else if(min[1]> thresholdMax){
            double scale = 
                Math.min((thresholdMax-10)-min[1], Math.max(-min[1], threshold-min[0]));
            backwardTrace.scale(scale, i);
            logscale+=scale;
        }
	   
        backwardTrace.setLogscale(i, logscale);
	}
	//AbstractTermination
	backwardTrace.overall = this.backward(0, -1);
   /* HmmerProfileHMM hm = (HmmerProfileHMM)hmm;
    AbstractTerm[] t_n = this.backwardTrace.trace.get(hm.n);
    AbstractTerm[] t_b = this.backwardTrace.trace.get(hm.begin);
    AbstractTerm[] t_m = this.backwardTrace.trace.get(hm.MAGIC);
    AbstractTerm[] t_c = this.backwardTrace.trace.get(hm.c);
    AbstractTerm[] t_e = this.backwardTrace.trace.get(hm.end);*/
    logger.fine("bsc "+backwardTrace.overall);
}
   
  public void calcScoresForward(boolean full) {
	if(logspace) throw new RuntimeException("only run forward in prob space");
	double[] res = new double[modelLength];
	for(int i=0; i<seqLength; i++){
        double logscale = forwardTrace.getLogScale(i-1);
        for(int j=1; j<modelLength; j++){
			if( full)  res[j] = this.forwardComplete(j, i,  emiss[j]  ? emissions[j][i] : nullEms);
            else{
                try{
                	res[j] =   this.forward(j, i,  this.emiss[j]   ? emissions[j][i] : nullEms);
                
                }catch(Exception exc){
                	System.err.println("problem at "+i+" with state "+j+" and "+this.obj.getName());
                    exc.printStackTrace();
                    System.exit(0);
                }
            }
		}
      if(Constants.CHECK && res[Constants.getMax(res)]<=0){
    	  double[] ems = new double[res.length];
    	  for(int j=1; j<modelLength; j++){
    		  ems[j] = emiss[j]  ? emissions[j][i] : nullEms;
    	  }
    		Double[][]trans = getTrans(i);
    	  Double[] tr = forwardTrace.getTr(i-1);
    	  throw new RuntimeException("!!" +Constants.print(res) +"\n"+Constants.print(ems)+"\n"+Constants.print(tr));
      }
      // if(Constants.CHECK) forwardTrace.allZero(i);
		double[] min = forwardTrace.minScore(i);
         if(Double.isInfinite(min[1])){
        	 Double[] res1 = forwardTrace.getTr(i);
        	 System.err.println(Arrays.asList(res1));
        //	 DistributionCollection.dc.print(i);
        	Double[][]trans = getTrans(i);
        	
        	//	System.err.println("trans from "+k+" "+Arrays.asList(trans[k]));
        	
             EmissionState es = obj;
             //double[] pr;
            // pr= es.getEmiss(i);
            /* if(es instanceof CachedEmissionState) {
                 pr= ((CachedEmissionState) es).emissions[i].probs;
             }
             else{
               
                   //  ((HaplotypeEmissionState) es).emissions[i].probs;
             }*/
             EmissionStateSpace ss = es.getEmissionStateSpace();
            List<Double> d = new ArrayList<Double>();
            List<String> str = new ArrayList<String>();
           int size =0;
           StringBuffer sb = new StringBuffer();  StringBuffer sb1 = new StringBuffer();
          /* for(int k=0; k<pr.length; k++){
               if(pr[k]>0){
                   d.add(pr[k]);
                   String string = ss.get(k).toString();
                   str.add(string);
                   if(string.length()>size) size = string.length();
               }
              
              
           }*/
           for(int k=0; k<d.size(); k++){
               sb.append("%"+(size+2)+"s ");
               sb1.append("%"+(size+2)+"."+(size-1)+"g ");
           }
           System.err.println(this.obj.getName());
           System.err.println(String.format(sb.toString(), str.toArray()));
           System.err.println(String.format(sb1.toString(), d.toArray()));
          EmissionStateSpace emstsp =  Emiss.getSpaceForNoCopies(obj.noCop());
           ProbabilityDistribution2 dist[] = new ProbabilityDistribution2[emstsp.size()]; 
           for(int j=0; j<dist.length; j++){
        	 int cn =   emstsp.getCN(j);
        	 int nob = emstsp.getBCount(j);
        	 if(DistributionCollection.dc!=null &&(short)obj.emissions(0).getDataIndex()>=0 )
        	   dist[j] = DistributionCollection.dc.getDistributionRB((short)obj.emissions(0).getDataIndex(), cn, nob, i);
           }
          PseudoDistribution dist1  = obj.emissions(i);
//           try{
                throw new RuntimeException("max should be greater than zero at "+i+" "+this.obj.getName()+" "+
                		DataCollection.datC.loc.get(i)+" exiting");
  //         }catch(Exception exc){
    //    	   exc.printStackTrace();
        	//   System.exit(0);
      //     }
             //           Arrays.asList(d) +"\n"+);
            }
         if(min[1]> thresholdMax){
             double scale = 
                 Math.min((thresholdMax-10)-min[1], Math.max(-min[1], threshold-min[0]));
             forwardTrace.scale(scale, i);
             logscale+=scale;
             if(Constants.CHECK && min[1]+scale>thresholdMax){
            	 throw new RuntimeException("!!");
             }
         }
         else if(min[0] < threshold){
			double scale = Math.min(-min[0], thresholdMax-min[1]); //multiply by exp of this
			//double sc1 = Math.exp(scale);
			//if(Double.isInfinite(sc1)) throw new RuntimeException("!! "+scale+" "+min[0]+" "+min[1]);
			forwardTrace.scale(scale, i);
			logscale+=scale;
		    if(Constants.CHECK && min[1]+scale>thresholdMax+10){
		    	
		    	throw new RuntimeException(min[0]+" "+min[1]+" "+scale+" "+thresholdMax);
		    }
		}
         else if(min[1] < thresholdMax1){
        	  double scale = 50-min[1];
                  Math.min(20-min[1], Math.max(-min[1], threshold-min[0]));
              forwardTrace.scale(scale, i);
              logscale+=scale;
            //  if(min[1]+scale>thresholdMax) throw new RuntimeException("!!");
         }
     
        forwardTrace.setLogscale(i, logscale);
	}
	//AbstractTermination automatic - magic index at state seqLength -1 
	forwardTrace.overall = full ? this.forwardComplete(0,seqLength-1 , nullEms ) : this.forward(0,seqLength-1 , nullEms);
	
 //  System.exit(0);

}
   
  private Double[][] getTrans(int i) {
		Double[][] trans = new Double[this.modelLength][this.modelLength];
		for(int k=0; k<trans.length; k++){
			for(int k1=0; k1<trans.length; k1++){
	    		trans[k][k1] = this.hmm.getTransitionScore1(k, k1, i);
	    	}
		}
			return trans;
	}
   
   /*assumes a log emission matrix but probability for transitions*/ 
  private double calcScoresViterbi() {
          //initialise matrix
      
          for(int i=0; i<seqLength; i++){
              for(int j=1; j<modelLength;j++){
//                  State j = it.next();
  //                if(j==hmm.MAGIC) continue;
               //   logger.info(j);
                  this.viterbi(j, i, this.emiss[j]  ? emissions[j][i] : this.nullEms);
              }
          }
          //AbstractTermination automatic - magic index at state seqLength -1 
          forwardTrace.overall = this.viterbi(0,seqLength-1 , this.nullEms);
         // logger.info(forwardTrace.overall);
          return forwardTrace.overall;
      }
   

  public void clear(){
    forwardTrace.clear();
    backwardTrace.clear();
}
  
 
     
  
  public double fillPosteriorProbabilities(int i,double[][] posterior){
    	double total =0;
      //  if(posterior.size()!=this.modelLength ) throw new RuntimeException("doesn't match");
       if(logspace){
    		
    		StatePosEmission spe = this.statePath.getStatePosForSeqIndex(i);
    		posterior[spe.j][i]=1.0;
    		total+= 1.0;
    	}else{
    		double scale_i = -forwardTrace.getLogScale(i)-backwardTrace.getLogScale(i) + forwardTrace.getLogScale(seqLength-1);
    	    
       		for(int j=0; j<modelLength; j++){
              //  State j = it.next();
       		//	if(hmm.statesL.get(j) instanceof DotState) continue;
                AbstractTerm forward = forwardTrace.getTrace(j,i);
                AbstractTerm backward = backwardTrace.getTrace(j,i);
                if(forward==null || backward==null) posterior[j][i] =0;// logger.info("forward null for "+hmm.getState(j));
                else{
                    double unscaled = (forward.score()*backward.score())/ forwardTrace.overall;
                    posterior[j][i] =  Math.exp(scale_i)*unscaled;
                    total+=posterior[j][i];
                }
       		}
    	}
       		return total;
    }
public void fillPosteriorProbabilities(double[][] posterior){
    for(int i=0; i<seqLength; i++){
        this.fillPosteriorProbabilities(i, posterior);
    }
}
double forward(int j, int i, double emissionScore)  {
	int max_j = -1;
    double max = 0;
	double sum = 0;
    int adv = this.adv[j];
   
    int[] toStates = hmm.statesIn(j, i);
	for(int k1=0; k1<toStates.length; k1++){
        int j1 = toStates[k1];
		double score;
        AbstractTerm AbstractTerm = forwardTrace.getTrace(j1,i-adv);
        if(AbstractTerm==null) score=0;
        else score =hmm.getTransitionScore1(j1, j, i)*AbstractTerm.score();
		if(score>max){// || (score==max && Math.random()>0.5)){
			max = score;
			max_j = j1;
		}
		
		
			
		if(Constants.CHECK){
			if(score<0){
				double tr = hmm.getTransitionScore1(j1, j, i);
				System.err.println("h"+tr);
				
			}
		/*	if(j==5 && j1==5  ){
				if(score<=0){
				//	Logger.global.info("h");
				}
				else if(score*emissionScore<=0){
				//	Logger.global.info("h");
				}
			}*/
       if( (Double.isNaN(score))){
           throw new RuntimeException("is nan "+AbstractTerm.score()+" "+j1+" "+emissionScore+" "+k1+" "+hmm.getTransitionScore1(j1, j, i));
    
       }
		}
		sum+=score;
	}
   if(Constants.CHECK &&( Double.isNaN(sum) || Double.isNaN(emissionScore))) throw new RuntimeException("is nan "+sum+" "+emissionScore);
   double res = sum* emissionScore;
   
	forwardTrace.setTrace(j,i, max_j,i-adv, res);
	return res;
}


double forwardComplete(int j, int i, double emissionScore)  {
    int max_j = -1;
    double max = 0;
    double sum = 0;
    int adv = this.adv[j];
  
    int[] toStates = hmm.statesIn(j, i);
    double[] score = forwardTrace.getDoubleArray(j, i);
    Arrays.fill(score, 0.0);
    for(int k1=0; k1<toStates.length; k1++){
        int j1 = toStates[k1];
        AbstractTerm AbstractTerm = forwardTrace.getTrace(j1,i-adv);
        if(AbstractTerm==null) score[j1]=0;
        else score[j1] =hmm.getTransitionScore1(j1, j, i)*AbstractTerm.score()*emissionScore;
        if(score[j1]>max){
            max = score[j1];
            max_j = j1;
        }
       if(Constants.CHECK && (Double.isNaN(score[j1]))){
           throw new RuntimeException("is nan "+AbstractTerm.score()+" "+j1+" "+emissionScore);
       }
        sum+=score[j1]; 
    }
   // if(Constants.CHECK &&( Double.isNaN(newTerm.score()))) throw new RuntimeException("is nan "+sum+" "+emissionScore);
    forwardTrace.setTrace(j,i, max_j,i-adv, sum);
    return sum;
}

/*double forwardComplete(int j, int i, double emissionScore)  {
    double sum = 0;
    int adv = this.adv[j];
    double[] score= forwardTrace.getDoubleArray(j,i);
    for(int j1 =0; j1<modelLength; j1++){
        AbstractTerm term = forwardTrace.getTrace(j1,i-adv);
        if(term==null) score[j1]=0;
        else score[j1] =hmm.getTransitionScore1(j1, j, i)*term.score()*emissionScore;
       if(Constants.CHECK && Double.isNaN(score[j1])) throw new RuntimeException("is nan "+term.score()+" "+j1+" "+emissionScore);
        
        sum+=score[j1]; 
    }
    forwardTrace.setTrace(j,i,   i-adv, sum);
    return sum;
}*/

private State  getBestEmission(int i){
    State max_j =null;
    double max = Double.NEGATIVE_INFINITY;
    for(Iterator<State> it = hmm.states(); it.hasNext();){
        State j = it.next();
        if(j instanceof EmissionState){
            double sc = emissions[j.getIndex()][i];
             if(     sc>max){
                max_j = j;
                max = sc;
            }
        }
    }
    return max_j;
}


   private  double[] getOverallScore(){
    	 double scoreForward = Math.log(forwardTrace.overall)-forwardTrace.getLogScale(seqLength-1);
    	 double scoreBackward = Math.log( backwardTrace.overall)-backwardTrace.getLogScale(0);
        if(Double.isNaN(scoreForward) || Double.isNaN(scoreBackward) || Double.isInfinite(scoreBackward)) throw new RuntimeException("is nan "+forwardTrace.overall+" "+forwardTrace.getLogScale(seqLength-1)+" "+backwardTrace.overall+" "+backwardTrace.getLogScale(0)+" "+this.protName);
    	 return new double[]{scoreForward,scoreBackward};
    }
    
   
   public StateDistribution getPosteriorML( int i, StateDistribution id){
         for(int j=1; j<hmm.modelLength(); j++){
            if(!emiss[j]) continue;
  //             AbstractTerm forward = forwardTrace.getTrace(j,i);
    //           AbstractTerm backward = backwardTrace.getTrace(j,i);
      //         if(forward!=null && backward!=null) // logger.info("forward null for "+hmm.getState(j));
        //      {
          //         double unscaled = (forward.score()*backward.score())/ forwardTrace.overall;
              //     if(unscaled>0) id.put(j,  scale_i*unscaled);
            id.put(j, this.emissions[j][i]);
            //   }
              
               
           }
       // }
           id.normalise();
           return id;
       }
    
    public StateDistribution getPosterior( int i, StateDistribution id){
    	   Arrays.fill(id.dist, 0.0);
    	if(logspace){
    		
    		StatePosEmission spe = this.statePath.getStatePosForSeqIndex(i);
    		id.dist[spe.j]=1;
    	}else{
    	
   
           double logscale_i = -forwardTrace.getLogScale(i)-backwardTrace.getLogScale(i) + forwardTrace.getLogScale(seqLength-1);
           for(int j=1; j<hmm.modelLength(); j++){
            if(!emiss[j]) continue;
               AbstractTerm forward = forwardTrace.getTrace(j,i);
               AbstractTerm backward = backwardTrace.getTrace(j,i);
               if(forward!=null && backward!=null) // logger.info("forward null for "+hmm.getState(j));
              {
                   double unscaled = (forward.score()*backward.score())/ forwardTrace.overall;
                   if(unscaled>0) id.put(j,  Math.exp(logscale_i+Math.log(unscaled)));
               }
              
               
           }
    	}
           if(Constants.CHECK) id.validate();
           return id;
       }
    
    public double[][] getPosteriorMatch(){
         //warning - used to only give result for match states
              double[][] posterior = new double[modelLength][this.seqLength];
              this.fillPosteriorProbabilities(posterior);
              return posterior;
      }
    
    public double[] getPosteriorMatchSum(State[] states){
      double[][] d = getPosteriorMatch();
      double[] d1 = new double[d[0].length];
      Arrays.fill(d1, 0);
      for(int j=0; j<states.length; j++){
          for(int i=0; i<d1.length; i++){
              d1[i]+=d[states[j].getIndex()][i];
          }
      }
      return d1;
     }


    



    public StatePath getStatePath(boolean sample){
    	if(!sample && this.statePath!=null) return statePath;
        StatePath path = new StatePath(protName, this.forwardTrace.overall, this.seqLength);
        getStatePath(path, new StoppingCriterion(){
            public boolean stop(Object emission){
                return false;
            }
        }, sample);
        return path;
    }
    
   
    
    public int getStatePath(StatePath path, StoppingCriterion stop, boolean sample)  { 
        State firstState = path.getFirstState();
        int j= (firstState==null) ? 0 : firstState.getIndex();
        Integer st = path.getFirstEmissionPos();
        int i  = st==null ? seqLength-1: st;
       // logger.info("starting with "+j+" at "+i);
       
        AbstractTerm tr = forwardTrace.getTrace(j,i);
        if(tr ==null) return i ;
        i = tr.i;
        j = tr.getBestPath();
       // int a = tr.a;
         while(j!=0 ){
        
             boolean emiss = this.emiss[j];
             try{   
                 
             Object obj_i =null;
             if(emiss){ 
                 EmissionState state_j = (EmissionState)hmm.getState(j);
               //  System.err.println("best "+state_j);
                 obj_i = 
                      //   Constants.fillGaps() ? 
                                 state_j.getBestIndex((HaplotypeEmissionState) obj,i, sample, Constants.isLogProbs()) ;
                           //      : obj.getBestIndices()[i]
                        
             }
                
                 path.add(hmm.getState(j), obj_i,  i ,j);
                 if(stop.stop(obj_i)) return i;
             }catch(Exception exc){
                exc.printStackTrace();
                logger.info(i+"");
                logger.info(j+"");
                logger.info(path+"");
                //logger.info(obj.getElement(i)+"");
                System.exit(0);
                 path.add(hmm.getState(j), null, i,j);
             }
            
             tr =  forwardTrace.getTrace(j,i);
             i = tr.i;
             j = tr.getBestPath();
           //  a = tr.a;
         }
         return i;
        // path.add(j, null, i);
    }
    
    /** need to multiply by a_kl  and scale to get real transition prob*/
    public  void addTransitionPosterior( int i, PseudoDistribution[] transitionPosterior, double weight){
    	
      if(logspace){
    		
    		StatePosEmission spe1 = this.statePath.getStatePosForSeqIndex(i+1);
    		 int j = i==-1 ? 0 : this.statePath.getStatePosForSeqIndex(i).j;
    		  PseudoDistribution d1 =transitionPosterior[j];
    		  d1.addCount(j, 1.0);
    	}
          double scale_i =  forwardTrace.getLogScale(seqLength-1)-forwardTrace.getLogScale(i)  -backwardTrace.getLogScale(i+1);
       /*  if( Double.isInfinite(scale_i)){
        	 double d1 = forwardTrace.getLogScale(seqLength-1);
        	 double d2 = forwardTrace.getLogScale(i);
        	 double d3 = backwardTrace.getLogScale(i+1);
        	 throw new RuntimeException("!! scale is infinite "+d1+" "+d2+" "+d3);
         }*/
          for(int k=0; k<modelLength; k++){
              PseudoDistribution d1 = transitionPosterior[k];
            //  d1.reset();
              for(int j=0; j<modelLength; j++){
                  double sc = Math.exp(Math.log(getTransitionProb(k,j, i)) + scale_i) *weight ;
                  if(Double.isInfinite(sc)){
                	  throw new RuntimeException("!!");
                  }
                  if(sc>0) d1.addCount(j, sc);
              }
          }
      }

    /** need to multiply by a_kl  and scale to get real transition prob*/
    private double getTransitionProb(int k, int l,  int i){
      int adv = this.adv[l];
      if(i+adv<0) return 0;
      double f_ki = forwardTrace.getScore(k, i, false);
      double b_li =backwardTrace.getScore(l, i+adv, false);
      double p_x = forwardTrace.overall  ;
      double e_li =  
          emiss[l]  ? 
              (i+adv==seqLength ? 0 :emissions[l][i+adv]) : 1.0;
      return  (f_ki*e_li*b_li) / p_x;
      //return unscaled;
     
  }

    
    
    
    public double searchML() {
     double sc = 0;
     double adj = emiss.length;
       for(int i=0; i<this.seqLength; i++){
           double sc_i =0;
           for   (int j=0; j<this.emiss.length; j++){
               if(emiss[j]){
                   sc_i+=this.emissions[j][i]*adj;
               }
           }
           sc+=Math.log(sc_i);
       }
        return sc;
     }
    
    Double logprob = null;
    StatePath statePath;
    public double search(boolean incBackward, boolean complete) {
    	if(logspace){ ///do viterbi training if emissions are in log-space
    		double sc =  this.calcScoresViterbi();
    		statePath = this.getStatePath(false);
    		logprob = sc;
    		return sc;
    	}else{
		 	   calcScoresForward(complete);
		 	  if(incBackward) calcScoresBackward();
		 	   double[] sc = getOverallScore();
		    
		      // if(Constants.CHECK){
		       if(Math.abs(sc[0]-sc[1])>0.001 || Double.isNaN(sc[0]) || Double.isNaN(sc[1])){
		           if(incBackward) throw new RuntimeException("these should match "+sc[0]+" "+sc[1]+"\n");
		     //  }
		       }
		    logprob = sc[0];
		       return sc[0];
    	}
    }
    public StatePath searchViterbi(){
       this.calcScoresViterbi();
       this.statePath = this.getStatePath(false);
        return statePath;	   
    }
    
   
    public double validate(StateDistribution dist, double sum){
          double sum1 = dist.sum();
          if(Math.abs(sum-sum1) > 0.01) throw new RuntimeException("sum is not right "+sum+" "+sum1);
          return sum1;
      }
/* works in log scale */
double viterbi(int j, int i, double emiss)  {
    int max_j =-1;
    int max_i = 0;
   
    double max = Double.NEGATIVE_INFINITY;
    int[] toStates = hmm.statesIn(j, i);
   for(int j3=0; j3<toStates.length; j3++){
        int j2 = toStates[j3];
        double score;
        int adv = this.adv[j];
        AbstractTerm term = forwardTrace.getTrace(j2,i-adv);
        if(term ==null) continue;
                  score =Math.log(hmm.getTransitionScore1(j2, j, i))+
                        +term.score();
            
                if(score>max){
                    max = score;
                    max_j = j2;
                    max_i = i-adv;
                }
    }
   
    double result = max;
    forwardTrace.setTrace(j,i, max_j, max_i, result+emiss);
    return result;
}






    
   
    
 


}

/*  public List<Integer> getPath(int j, int i){
List<Integer> l = new ArrayList<Integer>();
while(j>=0 && i>=0){
    l.add(0,j);
    State state = hmm.statesL.get(j);
    int j1 = j;
    j = this.forwardTrace.getTrace(j,i);
    i-= this.forwardTrace.adv[j1][i];
}
return l;
}*/


/* protected class TraceMatrix{
final public double[][] score;
final public int[][] adv; 
final public int[][] trace;


final public double[] logNullScore;



TraceMatrix(int modelLength, int seqLength)
{
    score = new double[modelLength][seqLength];
    trace = new int[modelLength][seqLength];
  adv = new int[modelLength][seqLength];
    
    logNullScore = new double[seqLength];

    for(int j=0; j<score.length; j++){
        Arrays.fill(score[j], 0);
        Arrays.fill(trace[j],-1);
     Arrays.fill(adv[j], 0);
    }
   
}



 /*public StatePath getStatePath()  { 
        int i = this.seqLength-1;
        int j = hmm.magicIndex;
        while(i>=0){
            int best_j=0;
            for(j=0; j<hmm.statesL.size(); j++){
                if(this.forwardTrace.score()[j][i] > this.forwardTrace.score()[best_j][i]){
                    best_j = j;//this.forwardTrace.score[j][i];
                }
            }
            if(this.forwardTrace.score[best_j][i] > Double.NEGATIVE_INFINITY){
                j =best_j;
                break;
            }
            i--;
        }
        logger.info("starting with "+j+" at "+i);
        StatePath path = new StatePath(hmm.MAGIC, hmm.MAGIC, this.protName);
         while(j>=0 && j<modelLength && i>=0){
             char st = hmm.statesL.get(j).name;
             logger.info("added "+st);
             path.add(hmm.statesL.get(j),st, i);
             State state = hmm.statesL.get(j);
             int j1 = j;
             j = this.forwardTrace.getTrace(j,i);
             if(state instanceof EmissionState) i-= forwardTrace.adv[j1][i];
           
         }
         return path;
  
         
    }



}*/




	
	
