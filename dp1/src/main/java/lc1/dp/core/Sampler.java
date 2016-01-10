  package lc1.dp.core;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import lc1.dp.data.collection.AssociationCalculator;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.HWECalculator;
import lc1.dp.data.representation.CSOData;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.genotype.io.trio.HalfTrioComparableArray;
import lc1.dp.genotype.io.trio.TrComparableArray;
import lc1.dp.genotype.io.trio.TrioComparableArray;
import lc1.dp.model.CachedHMM;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.model.FreeHaplotypeHMM;
import lc1.dp.model.MarkovModel;
import lc1.dp.model.WrappedModel;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.dp.states.State;
import lc1.dp.states.WrappedEmissionState1;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.apache.commons.pool.PoolableObjectFactory;
/*
 * @Todo problem sampling 'state'  genotypes.
 */
public class Sampler {
    //this is what will be phased.  Should be a copy of the data rather than the original data
    public DataCollection data;
    final String sb_i1;
  
 //   static final boolean buffer1 = false;
 //   static final boolean buffer =false;
    public File directory ;
    
    Map<String, List<CSOData>[] > spList = null;//new HashMap<String,List< CSOData>[]>();
  
    public void removeSamples(boolean delete){
        File[] sampleDir = directory.listFiles();
        for(int i=0; i<sampleDir.length; i++){
           if(!Constants.resample()){
               deleteDir(sampleDir[i]);
               if(sampleDir[i].exists()) throw new RuntimeException("exists!");
           }
        }
    }
   private void setDirectory(File f, boolean delete){
        this.directory = f;
        directory.mkdir();
        if( !directory.exists()) directory.mkdir();
       removeSamples(delete);
        //directory.deleteOnExit();
    }
   
    public List<CSOData>[] getBuffer(String key){
        if(this.spList.containsKey(key)){
            return spList.get(key);
        }
        else{
            List<CSOData>[] l = new List [] {
                    new ArrayList<CSOData>(),
                    new ArrayList<CSOData>()
            };
            spList.put(key, l);
            return l;
        }
    }
 public static PrintWriter[][] getOutputStream(File directory, String key) throws Exception{
           PrintWriter[][] res = new PrintWriter[2][2];
           for(int j1=0; j1<res.length; j1++){
               for(int i1=0; i1<res.length; i1++){
                   File f =  new File(directory, (i1==0 ? "hap_" :"states_")+(j1==0 ?"h_" : "g_") +key);
                   res[j1][i1] = new PrintWriter(new BufferedOutputStream(new FileOutputStream(f)));
               }
           }
           return res;
   }
   public static File[][] getInputStream(File directory, String key) throws Exception{
       File[][] res = new File[2][2];
        for(int j1=0; j1<res.length; j1++){
           for(int i1=0; i1<res.length; i1++){
               File f =  new File(directory, (i1==0 ? "hap_" :"states_")+(j1==0 ?"h_" : "g_") +key);
               res[j1][i1] = f;
           }
        }
     //  }
       return res;
   }
   
   List<CSOData> [][] get( String key, List<Integer> rep, EmissionStateSpace[] emStSp, int noSnps) throws Exception{
       List[][] res = new ArrayList[][] {{new ArrayList(), new ArrayList()},{new ArrayList(), new ArrayList()}};
       //or(int j1=0; j1<res.length; j1++){
       for(int i=0; i<rep.size(); i++){
          File[][] os =getInputStream(new File(directory, "hmm_"+rep.get(i)), key);
          for(int j=0; j<os.length; j++){
           res[j][0].addAll( getIterator(os[j][0], Emiss.class, emStSp[0], noSnps)); 
           res[j][1].addAll(getIterator(os[j][1], Integer.class, emStSp[1], noSnps));     
          }
           
       }
       //}
           return res;
        
   }
    
    public Sampler(DataCollection data, File dirname, boolean delete){
        this.data = data;
        setDirectory(dirname, delete);
        
       
       
       // if(!directory.exists()) directory.mkdir();
        StringBuffer sb_is1 = new StringBuffer();
        for(int i=0; i<data.length(); i++){
            sb_is1.append("%1i");
        }
        this.sb_i1  = sb_is1.toString();
        if(Constants.calcAssoc()){
        	 ac = data.getArmitageCalculator();
        	
        }else ac = null;
        if(Constants.hwe()){
        	 hwe_calc = data.getHWECalculator();
        }else hwe_calc = null;
    }
    static Comparator entryComp  = new Comparator<Map.Entry<List<Comparable>, Integer>>(){

        public int compare(Entry<List<Comparable>, Integer> o1, Entry<List<Comparable>, Integer> o2) {
           return o1.getValue().compareTo(o2.getValue());
        }
        
    };
    
  
    public String getPrintString(){
        return sb_i1;
    }
    
    
   /* private static boolean compare(ComparableArray orig, ComparableArray orig_mod){
    
       if(!Arrays.asList(orig.compString()).equals(Arrays.asList(orig_mod.compString()))) return false;
       else return true;
         
     }*/
    
    
    /** samples a new configuration of the original, based on the statepath sp*  */
    public static  void  sample(StatePath sp,CompoundMarkovModel hmm, boolean sample, PhasedDataState[] states){
       
      
        PhasedDataState result_phased =states[0];
        PhasedDataState result_switches =states[1];
       State[] prev_true = null;
        State[] prev_within = null;
         for(int i=0; i<sp.size(); i++){
        	// boolean plot = hmm.modelLength()>70 && states[0].name.equals("196811") && i==20;
             StatePath.StatePosEmission spe = sp.getSPE(i);
             CompoundState state_i = (CompoundState) spe.state;
             if(spe.pos!=i) throw new RuntimeException();
      //       EmissionState[] states_i = hmm.disambiguate(state_i.getMemberStates(false), prev, i,sample);
           State[] states_i_true = hmm.disambiguate(state_i.getMemberStates(true), prev_true, i, sample);
          State[] within = new State[states_i_true.length];
           for(int j=0; j<within.length; j++){
        	   //within phasing
        	   if(prev_within==null ){
        		   within[j] = states_i_true[j];
        		   if(within[j] instanceof WrappedEmissionState1){
        			   within[j] = ((WrappedEmissionState1)within[j]).inner;
        		   }
        	   }
        	   else{
        		   within[j] =hmm.disambiguate((EmissionState)states_i_true[j], (EmissionState)prev_within[j], i, sample,j);
        	   }
           }
            if(result_switches!=null) {
            	result_switches.addPoint(i, makeArray(states_i_true));
            }
         
          
        	   Integer selected = 
                   (Integer) sp.getEmission(i) ;// orig_i;
         CompoundEmissionStateSpace emStSp = (CompoundEmissionStateSpace)state_i.getEmissionStateSpace();  
         
             int[] possibilities = emStSp.getHaps(selected.intValue()); //indices of hapls
        //   
             int max_ind = 0;
             if(Constants.CHECK && possibilities.length==2 && (emStSp.getHapl(possibilities[0]).toString().equals(emStSp.getHapl(possibilities[1]).toString())))
                 throw new RuntimeException("!! "+possibilities[0]+" "+possibilities[1]+ " "+emStSp.getHapl(possibilities[0])+" "+emStSp.getHapl(possibilities[1]));
             if(possibilities.length>1){
                 double[] prob = PairEmissionState.pool.getObj(possibilities.length);
                     //new double[possibilities.length];
                 double sum=0;
                 inner: for(int k=0; k<prob.length; k++){
                	  int[] type = emStSp.getMemberTypeIndices(possibilities[k]);
                     int[] me = emStSp.getMemberIndices(possibilities[k]);
                   //  int[] types = emStSp.getMemberTypeIndices(possibilities[k]);
                     prob[k] = 1.0;//states_i_true[0].score(me[0], i);
                     for(int k1 =0; k1<me.length; k1++){
                    	 EmissionStateSpace emstsp1 = emStSp.getMembers()[type[k1]];
                    	 EmissionStateSpace emstsp2 = ((EmissionState)within[k1]).getEmissionStateSpace();
                    	boolean equals = emstsp1.defaultList.size()==emstsp2.defaultList.size();
                    	 if(equals){
                             prob[k] *= ((EmissionState)within[k1]).score(me[k1], i);
                           
                    	 }
                    	 else{
                    		 
                    		 prob[k] =Double.NaN;
                    		 continue inner;
                    	 }
                     }
//                         state_i.score(possibilities[k], i, true, false);//*w[k];
                     sum+=prob[k];
                 }
                 max_ind =  sample? Constants.sample(prob, sum) : Constants.getMax(prob);
                
                  PairEmissionState.pool.returnObj(prob);                                    
                                                                  
                                                                  
             }
             
             selected =possibilities[  max_ind   ]   ;
            // int[] types = emStSp.getMemberTypeIndices(selected);
             int[] type = emStSp.getMemberTypeIndices(selected);
             ComparableArray hapl= (ComparableArray) emStSp.getHapl(selected);
             
             //now phase within
            for(int j=0; j<hapl.size(); j++){
            	 if(hapl.get(j) instanceof ComparableArray){
            		//int[] ind =  ((CompoundEmissionStateSpace) emStSp).getMemberTypeIndices(selected);
            		
            		 CompoundEmissionStateSpace inner =(CompoundEmissionStateSpace) ( emStSp).getMembers()[type[j]];
            		 int[] possibilities1 = inner.getHaps(inner.getGenotype(hapl.get(j)));
            		 if(possibilities1.length>1){
            			  double[] prob = PairEmissionState.pool.getObj(possibilities1.length);
                          //new double[possibilities.length];
                      double sum=0;
                    
                      EmissionState[] states_i_true1 = ((CompoundState)within[j]).getMemberStates(true);
                      for(int k=0; k<prob.length; k++){
                          int[] me = inner.getMemberIndices(possibilities1[k]);
                         int[] types1 = inner.getMemberTypeIndices(possibilities1[k]);
                          prob[k] = 1.0;//states_i_true[0].score(me[0], i);
                          for(int k1 =0; k1<me.length; k1++){
                        	  EmissionStateSpace emstsp1 = inner.getMembers()[types1[k1]];
                         	 EmissionStateSpace emstsp2 = states_i_true1[k1].getEmissionStateSpace();
                         	boolean equals = emstsp1.defaultList.size()==emstsp2.defaultList.size();
                         
                         	 if(equals){
                                  prob[k] *= states_i_true1[k1].score(me[k1], i);
                         	 }
                         	 else prob[k] =0.0;
                          }
//                              state_i.score(possibilities[k], i, true, false);//*w[k];
                          sum+=prob[k];
                      }
                       int selected_ =possibilities1[    sample? Constants.sample(prob, sum) : Constants.getMax(prob)  ]   ;
                       PairEmissionState.pool.returnObj(prob);    
                       hapl.set(j, inner.get(selected_));
            		 }
            	 }
             }
             result_phased.addPoint(i, hapl);//selected.copy());
           //  prev = states_i;
             prev_true = states_i_true;
             prev_within = within;
         }
         
        // return new PhasedDataState[] {result_phased,result_switches};
        }
    
   
    
   
   boolean check = false;
   
   public void addObject(Object[] obj, CSOData[] objects,   int k) throws Exception{
      /* if(buffer){
         ( ( List<CSOData>) obj[0]).add(objects[0]);
         ( ( List<CSOData>) obj[1]).add(objects[1]);
       }
       else{*/
          objects[0].print((PrintWriter)obj[0], true, false, false, null, null, null,null);
          ((PrintWriter)obj[0]).flush();

          if(objects[1]!=null){
          objects[1].print((PrintWriter)obj[1], true, false, false, null, null, null, null);
        
          ((PrintWriter)obj[1]).flush();
          }
       //}
     
   }
   
   public static boolean deleteDir(File dir) {
       if (dir.isDirectory()) {
           String[] children = dir.list();
           for (int i=0; i<children.length; i++) {
               boolean success = deleteDir(new File(dir, children[i]));
               if (!success) {
                   return false;
               }
           }
       }
   
       // The directory is now empty so delete it
       return dir.delete();
   }
   
  public static  double getProbOverStates(StateDistribution emissionC,
           MarkovModel hmm, HaplotypeEmissionState obj, int i, double[] res, boolean log, double[] dist) {
	   EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
       Arrays.fill(res,  0.0);
       double sum = 0;
       for(int j=1; j<emissionC.dist.length; j++){
           double p = emissionC.dist[j];
           EmissionState state_j = (EmissionState) hmm.getState(j);
           if(p>0){
          	 // double[] dist = state_j.distribution();
               sum+=HaplotypeEmissionState.calcDistribution(obj, (EmissionState) state_j, i,  p,log, dist);
               EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
               for(int k=0; k<dist.length; k++){
            	   int haploPair = emstsp.get(emstsp1.get(k)).intValue();
            	   int geno = emstsp.getGenoForHaplopair(haploPair);
             	  res[geno]+=dist[k];
               }
           }
       }
       for(int k=0; k<res.length; k++){
           res[k] = res[k] / sum;
       }
//int cn =        emstsp.getCN(Constants.getMax(res));
	//System.err.println("cn is "+cn);
         return sum;
   }
  
  public static  double getProbOverStatesHaplo(StateDistribution emissionC,
          MarkovModel hmm, HaplotypeEmissionState obj, int i, double[] res) {
	   EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
      Arrays.fill(res,  0.0);
      double sum = 0;
      for(int j=1; j<emissionC.dist.length; j++){
          double p = emissionC.dist[j];
          EmissionState state_j = (EmissionState) hmm.getState(j);
          if(p>0){
        	  double[] dist = state_j.distribution();
              sum+=HaplotypeEmissionState.calcDistribution(obj, (EmissionState) state_j, i,  p,Constants.isLogProbs(), dist);
              EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
        	 
              for(int k=0; k<dist.length; k++){
           	   int haploPair = emstsp.get(emstsp1.get(k)).intValue();
           	 
            	  res[haploPair]+=dist[k];
              }
          }
      }
      for(int k=0; k<res.length; k++){
          res[k] = res[k] / sum;
      }
//int cn =        emstsp.getCN(Constants.getMax(res));
	//System.err.println("cn is "+cn);
        return sum;
  }
  
  public static  double getProbOverStates(StateDistribution emissionC,
          MarkovModel hmm, HaplotypeEmissionState obj, int i, double[] res, int mixComp) {
	   EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
      Arrays.fill(res,  0.0);
      double sum = 0;
      for(int j=1; j<emissionC.dist.length; j++){
          double p = emissionC.dist[j];
          EmissionState state_j = (EmissionState) hmm.getState(j);
          if(p>Constants.countThresh3()){
        	  double[] dist = state_j.distribution();
              sum+=HaplotypeEmissionState.calcDistribution(obj, (EmissionState) state_j, i,  p,Constants.isLogProbs(),mixComp, dist);
              EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
        	
              for(int k=0; k<dist.length; k++){
           	   int haploPair = emstsp.get(emstsp1.get(k)).intValue();
           	   int geno = emstsp.getGenoForHaplopair(haploPair);
            	  res[geno]+=dist[k];
              }
          }
      }
      for(int k=0; k<res.length; k++){
          res[k] = res[k] / sum;
      }
        return sum;
  }
  
  
   List<MyObjectPool> dp_pool = new ArrayList<MyObjectPool>();
   List< MyObjectPool> stateDist = new ArrayList<MyObjectPool>();
   //int modelCount=0;;
   boolean unwrapToSample = Constants.unwrapForSampling();
  
   final List<AssociationCalculator>[][] ac;
   final HWECalculator[] hwe_calc;
List<EmissionStateSpace> stateEm1 = new ArrayList<EmissionStateSpace>();
   public Callable sampleFromHMM(final int j1, final EmissionState dat1, final int noSamples, final List models, final File sampleDir) throws Exception{
	  return new Callable(){
		   
	   public Object call(){
		   
	    try{
	   boolean sample = noSamples >1;
       
    
      
       
      // EmissionState dat1 = it.next();//bwt[ik2].get(ik1);
       //  no_samples++;
       String key = dat1.getName();
       int nocop = dat1.noCop();
         CompoundMarkovModel hmm =(CompoundMarkovModel)models.get(nocop-1);
       if(  unwrapToSample){
           while(hmm instanceof WrappedModel){
            hmm = ((WrappedModel)hmm).getHMM();
           } 
       }
       else if(hmm instanceof CachedHMM ){
           ((CachedHMM)hmm).refresh();
       }
       
       if(ac!=null && ac[nocop-1]!=null){
    	   for(int ii=0; ii<ac[nocop-1].length; ii++){
    		   for(int jj=0; jj<ac[nocop-1][ii].size(); jj++){
    			   ac[nocop-1][ii].get(jj).setModel(hmm);
    		   }
    	   }
       }
      // System.err.println("getting "+j1);
       EmissionStateSpace stateEm=  null;
      
       if(Constants.transMode(1)==null && Constants.saveStates()){
    	   
    	   stateEm = hmm.getStateSpace();;//
    	  for(int i=stateEm1.size(); i<=nocop; i++){
    		  stateEm1.add(null);
    	  }
    	   stateEm1.set(nocop, stateEm);
       }
         //   new DP(hmm, dat1.getName(),dat1, false);
            //bwt.dp[ik1];
       DP dp = (DP) dp_pool.get(data.noCopies(key)-1).getObj(j1);
       dp.setData(dat1,j1);
       dp.reset(true);  
       
        if(!sample) dp.searchViterbi();
        else  dp.search(true, sample);
     //   int len=1;
        Object[][] os1 =  !sample && Constants.numRep()<=1 ? null : getOutputStream(sampleDir, key) ;
        CompoundEmissionStateSpace emStSp =Emiss.getSpaceForNoCopies(nocop);
        	//(CompoundEmissionStateSpace) ((EmissionState)hmm.getState(1)).getEmissionStateSpace(); 
      
                
        PhasedDataState[] sam = new PhasedDataState[] {
                SimpleScorableObject.make(dat1.getName(), dat1.noSnps(), Emiss.getSpaceForNoCopies(nocop),
                		(short) ((HaplotypeEmissionState)dat1).dataIndex()),
                Constants.saveStates() ? SimpleScorableObject.make(dat1.getName(), 
                		dat1.noSnps(), stateEm,
                		(short) ((HaplotypeEmissionState)dat1).dataIndex()) : null
               
        };
        PhasedDataState[] samL = null; 
        double[] certainty = null; //can I get rid of this??
        boolean updateDirect = !sample && Constants.numRep()<=1;
        if(updateDirect) { //update data directly and don't write out anything
        	samL = new PhasedDataState[] {
           		  SimpleScorableObject.make(dat1.getName(), dat1.noSnps(), emStSp
           				  , (short) ((HaplotypeEmissionState)dat1).dataIndex()),
                     Constants.saveStates() ? SimpleScorableObject.make(dat1.getName(), dat1.noSnps(), 
                    		 stateEm,(short) ((HaplotypeEmissionState)dat1).dataIndex()) : null
           };
        	
        }
        if(Constants.calcLD()){
        	certainty = new double[dat1.noSnps()];//
        }
        StateDistribution emissionC = (StateDistribution) stateDist.get(data.noCopies(key)-1).getObj(j1);
      
           for(int i=0; i<hmm.noSnps; i++){
	            dp.getPosterior(i, emissionC);
	         //   int i1 = Constants.getMax(emissionC.dist);
	            EmissionStateSpace emiss = Emiss.getSpaceForNoCopies(nocop);
	            if(sample){
	          	  double[] prob =  PairEmissionState.pool.getObj(emiss.genoListSize());
	          	 double sumNonMix = 
	          		 Constants.allowComponent() ? 	getProbOverStates(emissionC, hmm,(HaplotypeEmissionState) dat1, i,prob, 0) : 0;
	          	
	          	  
	          	 double sum = getProbOverStates(emissionC, hmm,(HaplotypeEmissionState) dat1, i,prob,Constants.isLogProbs(), dp.distribution);
	          	if(!Constants.allowComponent()) sumNonMix = sum;
	          	 
		          double maxV = prob[Constants.getMax(prob)];
		         
	          //	 if(Double.isNaN(prob[0]) || sum==0 ){
	          	//	 Logger.global.info("h");
	          	 //}   
	          	 if(sumNonMix/sum>= Constants.imputedThresh(1)){
	            if(Constants.calcAssoc() && ac!=null && ac[nocop-1]!=null){
	            	List<AssociationCalculator>[] l = ac[nocop-1];
	            	for(int jj=0; jj<l.length; jj++){
		            	for(int k=0; k<l[jj].size(); k++){
		            		if(maxV>Constants.imputedThresh(0)){
		            		l[jj].get(k).scoreChi1(prob,  i, true,  ((HaplotypeEmissionState)dat1).getName());
		            		}
		            		l[jj].get(k).scoreChi1(emissionC,  i, true,  ((HaplotypeEmissionState)dat1).getName());
		            	}
	            	}
		         }
	            if(Constants.hwe() && hwe_calc!=null && hwe_calc[nocop-1]!=null ){
	            	hwe_calc[nocop-1].scoreChi1(prob, ((HaplotypeEmissionState)dat1).dataIndex(),i);
	            }
	           
	          	 }
	          else{
	        	  Arrays.fill(prob, Double.NaN);
	          }
	            
	            if(Constants.saveStates()){ 
	            	   if(os1!=null){
	            		   PrintWriter states = (PrintWriter)os1[1][1];
	            		   //warning - we assume hmm states are same order as emissionstate space
	            		   print(states, emissionC.dist, 1);
	            	   }
	            	    double[] state =PairEmissionState.pool.getObj(stateEm.defaultList.size());
	            	    Arrays.fill(state, 0.0);
	            	    for(int k=0; k<state.length; k++){
	            	    	int[] hap = stateEm.getHaploFromHaploPair(k);
	            	    	for(int j=0; j<hap.length; j++){
	            	    		if(hap[j]+1<emissionC.dist.length){
	            	    		state[k]+=emissionC.dist[hap[j]+1];
	            	    		}
	            	    	}
	            	    }
	            	  
	            	    if(Constants.CHECK){
	            			double sum2 = Constants.sum(state);
	            			if(Math.abs(sum2-1.0)>0.01){
	            				emissionC.validate();
	            				throw new RuntimeException("!!");
	            			}
	            			}	
	            	  if(samL!=null ) samL[1].emissions[i] = new SimpleExtendedDistribution(state, Double.POSITIVE_INFINITY);
	            }
	            if(samL!=null){
	            	
	            	samL[0].emissions[i] = new SimpleExtendedDistribution(prob, Double.POSITIVE_INFINITY);
	            	if(certainty!=null){
	            		certainty[i] = prob[Constants.getMax(prob)];
	            	}
	            }
	            else if(os1!=null){
	            	 PrintWriter geno = (PrintWriter)os1[1][0];
	            	  double[] prob1 =  PairEmissionState.pool.getObj(Emiss.getSpaceForNoCopies(nocop).defaultList.size());
	            	 getProbOverStatesHaplo(emissionC, hmm,(HaplotypeEmissionState) dat1, i,prob1);
	            	print(geno, prob1,0);
	            	   PairEmissionState.pool.returnObj(prob1);
	            }
	            PairEmissionState.pool.returnObj(prob);
	            }
        }
        for(int k=0; k<noSamples; k++){
            StatePath sp = dp.getStatePath(sample);
            sample(sp, hmm, sample, sam);
            sam[0].setName(noSamples==1 ? key : key+"_"+k);
            if(sam[1]!=null) sam[1].setName(noSamples==1 ? key : key+"_"+k);
            if(os1!=null) addObject(os1[0],sam ,   k);
        }
        dp_pool.get(data.noCopies(key)-1).returnObj(j1);
        stateDist.get(data.noCopies(key)-1).returnObj(j1);
       // if(!buffer){
        if(os1!=null){
            for(int i=0; i<os1.length; i++){
                for(int j=0; j<os1[i].length; j++){
                    ((PrintWriter)os1[i][j]).close();
                }
            }
        }
       // }
       // System.err.println("returning "+j1);
      
        if(updateDirect){
        	sam[0].setNoCop(nocop);
        	samL[0].setNoCop(nocop);
        if(sam[1]!=null){	sam[1].setNoCop(nocop);
        	samL[1].setNoCop(nocop);
        }
        	//if(Constants.saveInference()){
      // FreeHaplotypeHMM hmmin =  (FreeHaplotypeHMM)hmm.getMemberModels()[0];
   //   EmissionStateSpace stsp =  hmmin.getStateSpace();
        if(!sample){
        	((DataCollection)data).setData(key, sam[0], sam[0], sam[1], sam[1], certainty);
        }else{
          	((DataCollection)data).setData(key, sam[0], samL[0], sam[1], samL[1], certainty);
        }
        	//else{
        		//((DataCollection)data)
        	//}
        }
        }catch(Exception exc){
        	System.err.println("problem with "+j1);
        	exc.printStackTrace();
        }
	    return null;
	   }
	  };
   }
   private void print(PrintWriter states, double[] dist, int start) {
       Object[] obj = new Object[1];
       for(int i=start; i<dist.length; i++){
           obj[0] = dist[i]; 
           states.print(String.format("%5.3f", obj));
           if(i<dist.length-1){
               states.print("\t");
           }
           else{
               states.println();
           }
   }
    
}
public void sampleFromHMM(final List models, final int noSamples, int modelCount,List<EmissionState>[] bwt){
       //if(bwt.data.length!=this.data.size()) throw new RuntimeException("!!");
	   if(ac!=null){
		   for(int i=0; i<ac.length; i++){
			   if(ac[i]!=null){
				   for(int jj=0; jj<ac[i].length; jj++){
				   for(int k=0; k<ac[i][jj].size(); k++){
					   ac[i][jj].get(k).initialise(); 
				   }
				   }
			   }
		   }
	   }
      for(int i=0; i<models.size(); i++){
          final int i1 = i;
          DPPool pool = null;
          MyObjectPool poolDist = null;
          if(models.get(i)!=null){
             pool =  new DPPool((MarkovModel)models.get(i1), noSamples, unwrapToSample,Math.max(2, Constants.numThreads()), Constants.numThreads());
             
               poolDist =  new MyObjectPool(new PoolableObjectFactory(){

                 public void activateObject(Object arg0) throws Exception {}
                 public void destroyObject(Object arg0) throws Exception {}
                 public Object makeObject() throws Exception {
                   //  System.err.println("MAKING NEW DP OBJECT ");
                     CompoundMarkovModel hmm =(CompoundMarkovModel)models.get(i1);
                     if(unwrapToSample){
                         while(hmm instanceof WrappedModel){
                          hmm = ((WrappedModel)hmm).getHMM();
                         } 
                     }
                     return  new StateDistribution(hmm.modelLength());
                 }
                 public void passivateObject(Object arg0) throws Exception {}
                 public boolean validateObject(Object arg0) {return true;}}   ,Math.max(2, Constants.numThreads()), Constants.numThreads());
          }
          this.dp_pool.add(pool);
          this.stateDist.add(poolDist);
      }
       if( !directory.exists()){ 
         directory.mkdir();
         if(!directory.exists()) throw new RuntimeException("could not make dir!!");
       }
       File sampleDir = new File(directory, "hmm_"+modelCount);
           sampleDir.mkdir();
      //     Logger.global.info("sampling "+modelCount+" "+sampleDir);
       modelCount++;
       boolean sample = noSamples > 1;
  //     PrintWriter[]pw_rep=null;
    //   if(buffer1)pw_rep= new PrintWriter[noSamples];
      
       int index=0;
       try{
         /*  if(pw_rep!=null){
               pw_rep = new PrintWriter[noSamples];
               for(int k=0; k<pw_rep.length; k++){
                   pw_rep[k] = new PrintWriter(new BufferedWriter(new FileWriter(new File(sampleDir,"rep_"+k))));
               }
           }*/
       // outer2: for(int ik2 =0; ik2<bwt.length; ik2++){
          //  if(bwt[ik2]==null) continue;
    	   List<String> keys =data.indiv();
    	   List tasks = new ArrayList<Callable>();
           outer: for(//int ik1 =0; ik1<bwt[ik2].size(); ik1++){
                   Iterator<String> it = keys.iterator(); it.hasNext();){
        	   String key = it.next();
        	 //  System.err.println("sampling "+key);
             EmissionState emst =   data.dataL.get(key);
        //       Logger.global.info("sampling "+modelCount+" "+emst.getName());
             try{
           tasks.add(sampleFromHMM(index, emst, noSamples, models, sampleDir));// pw_rep);
             }catch(Exception exc){
            	 exc.printStackTrace();
             }
           index+=1;
          }
    	   BaumWelchTrainer.involeTasks(tasks,true);
       if(Constants.bufferCompress()){
    	   data.printcompressed(keys);
       }
          /* if(pw_rep!=null){
               for(int k=0; k<pw_rep.length; k++){
                   pw_rep[k].close();
               }
           }*/
       // }
       }catch(Exception exc){
           exc.printStackTrace();
          // System.exit(0);
       }
    }
   
   


 /*public PIGData readNext(BufferedReader br, List<String> l, String name) throws Exception{
     String st = br.readLine();
     if(st==null) return null;
     l.clear();
     while(st!=null && !st.startsWith("#")){
         l.add(st);
         st = br.readLine();
     }
     if(l.size()==2){
         return SimpleScorableObject.make(name,l, Emiss.class);
     }
     else if(l.size()==6){
         PIGData[] pid = new PIGData[3];
         for(int i=0; i<pid.length; i++){
             pid[i] = SimpleScorableObject.make(name, l.subList(2*i, 2*(i+1)), Emiss.class);
         }
         PIGData res =SimpleScorableObject.make(pid, false, ";");
         res.setName(name);
         return res;
     }
     else throw new Exception("!!");
 }*/
   
 
/** updates dat_em and dat with best path */
 public void calcBestPathSampling(HaplotypeEmissionState dat_em, 
         EmissionStateSpace[] ststSp,
         List<Integer> rep) {
     try{
        String key = dat_em.getName();
        List<CSOData>[][] spList =  get(key, rep, ststSp, dat_em.noSnps());
         double[] uncertaintyPhase = data.uncertaintyPhase(key);
         double[] uncertaintyVitPhase = data.uncertaintyVitPhase(key);
     
         Arrays.fill(uncertaintyPhase, 1.0);
         Arrays.fill(uncertaintyVitPhase, 1.0);
        PhasedDataState datvit, datvitL;
        datvit = datvitL = null;
        if(Constants.saveStates()){
	        datvit=  ((PhasedDataState) spList[0][1].get(0)).clone();
	        datvitL =  ((PhasedDataState) spList[1][1].get(0)).clone();
        }
        PhasedDataState dat = ((PhasedDataState) spList[0][0].get(0)).clone();
         HaplotypeEmissionState datL
          = ((PhasedDataState) spList[1][0].get(0)).clone();
        if( dat_em.dataIndex()>=0)  datL.setDataIndex(dat_em.dataIndex());
         datL.setPhenotype(dat_em.phenValue());
         datL.setNoCop(dat.noCopies());
       
      //   ((DataCollection)data).setData(key, sam[0], samL[0], sam[1], samL[1], certainty);
         
         if(spList[0][0].size()==1){
            
            
         }
         else{
             datL.sampleGenotype(spList[1][0]);
             transferMax(datL, dat);
             dat.samplePhase(spList[0][0], uncertaintyPhase,null);
             if(datvitL!=null){
                 datvitL.sampleGenotype( spList[1][1]); ///NO LONGER SAMPLING STATE!!!!!
                 transferMax(datvitL, datvit);
                 
             }
             datvit.samplePhase(spList[0][1], uncertaintyVitPhase, this.stateEm1.get(dat.noCopies()));
        }
         
         ((DataCollection)
                 data).setData(key, dat, datL, datvit, datvitL, uncertaintyPhase);
         if(Constants.saveStates()) data.setRecSites(dat.getName(), sampleRecSites(spList[0][1]));
       }catch(Exception exc){
           exc.printStackTrace();
       }
     
 }
 
 
 
 private void transferMax(HaplotypeEmissionState datL, PhasedDataState dat) {
  for(int i=0; i<datL.noSnps(); i++){
      dat.emissions[i].setFixedIndex(datL.emissions[i].getMax());
  }
    
}
private int noShared(ComparableArray comp_train, ComparableArray comp__no_train) {
     List<Comparable> l = new ArrayList(comp_train.elements());
     List<Comparable> l1 = new ArrayList(comp__no_train.elements());
     int cnt =0;
     for(int i=0; i<l.size(); i++){
         if(l1.remove(l.get(i))) cnt++;
     }
     return cnt;
  }
 /*
public  void calcPairing(String key1, String key2, Double[] uncertainty,Integer[] shared, int modelCount){
     try{
         List<PIGData>statesL1 = new ArrayList<PIGData>();
         List<double[]> cert1 = new ArrayList<double[]>();
         read(key1, statesL1, cert1, modelCount);
         List<PIGData>statesL2 = new ArrayList<PIGData>();
         List<double[]> cert2 = new ArrayList<double[]>();
         read(key2, statesL2, cert2, modelCount);
         double[][]count = new  double[uncertainty.length][];
         for(int i=0; i<count.length; i++){
             count[i] = new  double[] {0,0,0};
         }
         for(int j=0; j<statesL1.size(); j++){
             for(int i=0; i<count.length; i++){
                 count[i][noShared((ComparableArray) statesL1.get(j).getElement(i), (ComparableArray) statesL2.get(j).getElement(i))]
                          +=cert1.get(j)[i]*cert2.get(j)[i];
             }
         }
         for(int i=0; i<count.length; i++){
             int max =0;
             for(int j=1; j<count[i].length; j++ ){
                 if(count[i][j] > count[i][max]){
                     max = j;
                 }
             }
             shared[i] = max;
             uncertainty[i] = (double) count[i][max] / (double) statesL1.size() ;
         }
     }catch(Exception exc){
         exc.printStackTrace();
     }
 }*/

 
    public void  calcBestPathSampling(List<Integer> rep) {
       Logger.global.info("sampling from hmm");
     //   this.data.uncertainty =new HashMap<String,double[]>();
        int index=0;
       // boolean viterbi = Constants.viterbi();
     
         
        for(Iterator<String> it = this.data.dataL.keySet().iterator(); it.hasNext();index++){
        	String key  = it.next();
        	HaplotypeEmissionState emSt= (HaplotypeEmissionState) data.dataL.get(key);
         //  PhasedDataState dat= (PhasedDataState) it.next();
           // HaplotypeEmissionState emSt = (HaplotypeEmissionState) this.data.getL(dat.getName());
            if(!(emSt instanceof HaplotypeEmissionState)){
            	HaplotypeEmissionState em1 = new HaplotypeEmissionState(emSt.getName(), emSt.length(),
            			emSt.getEmissionStateSpace().size(), emSt.getEmissionStateSpace(),  null, null);
                em1.setDataIndex(emSt.dataIndex());
                ((DataCollection)data).dataL.put(emSt.getName(), emSt = em1);   }
            int ploidy = emSt.noCop();
            EmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
            calcBestPathSampling(emSt, 
                    new EmissionStateSpace[] {stsp, 
                    stateEm1.get(ploidy)},  rep);
            }
        Logger.global.info("... done sampling");
        if(!Constants.keepSamples())this.removeSamples(true);
       // this.data.calculateMaf(true);
        }
    

// int no_samples =0;
   public List<PIGData> getIterator(File f, Class clazz, EmissionStateSpace emStsp,  int noSnps) throws Exception{
	   if(!f.exists() || f.length()==0) {
		   System.err.println("warning file did not exist "+f.getAbsolutePath());
		   return new ArrayList();
	   }
	   final BufferedReader stream = new BufferedReader(new FileReader(f));
       boolean haps = f.getName().indexOf("_h_")>=0;
    
      List<PIGData> res =   haps ? new ArrayList(
               DataCollection.readFastPhaseOutput(stream, clazz, emStsp).data.values()) : 
                   Arrays.asList(new PIGData[] {DataCollection.readBasic(stream, clazz, emStsp, noSnps)});
               stream.close();
               return res;
        }
     
    private SortedMap<Integer, Integer>[] sampleRecSites(List<CSOData> spList) {
        ComparableArray obj = (ComparableArray)spList.get(0).getElement(0);
        int gap = 0;
       if(!(obj instanceof TrComparableArray)) return null;
        SortedMap<Integer, Integer>[] m = new SortedMap[((TrComparableArray)obj).third.size()];
        for(int j=0; j<m.length; j++){
            m[j] = new TreeMap<Integer, Integer>();
        }
        for(Iterator<CSOData> it = spList.iterator(); it.hasNext();){
            CSOData nxt = it.next();
             ComparableArray prev =((TrComparableArray)nxt.getElement(0)).third; 
            for(int i=1; i<nxt.length(); i++){
               ComparableArray third = ((TrComparableArray)nxt.getElement(i)).third;
               for(int j=0; j<m.length; j++){
                   if(!third.get(j).equals(prev.get(j))){
                       int i1 = i;
                   //    for(int i1 = Math.max(0, i-gap); i1 < Math.min(nxt.length(), i+gap); i1++){
                           Integer cnt = m[j].get(i1);
                           m[j].put(i1,cnt==null ? 1 : cnt+1 );
                       //}
                   }
               }
               prev = third;
            }
        }
       for(int j=0; j<m.length; j++){
           for(Iterator<Entry<Integer, Integer>> it = m[j].entrySet().iterator(); it.hasNext();){
               if(it.next().getValue() < Constants.sampleThresh() *(double) spList.size()){
                   it.remove();
               }
           }
       }
        return m;
    }


  /*  public void  calcRecSitesSampling()  throws Exception{
           int index=0;
           for(Iterator<String> it = this.data.getKeys().iterator(); it.hasNext();index++){
              // Arrays.fill(data.uncertainty[index], 1.0);
                //note we need to change for trisomic
                   try{
                    String key = it.next();
              
                     List<CSOData> spList;
                         spList = this.getSpList(key, 1);
                             //spLists[0].get(index);
                     double[] avgNo =new double[] {0,0};
                     for(int i=0; i<spList.size(); i++){
                         Set[] switches = ((PIGData)spList.get(i)).getSwitches();
                         for(int k=0; k<avgNo.length; k++){
                             avgNo[k]+=switches[k].size();
                         }
                     }
                     for(int k=0; k<avgNo.length; k++){
                         avgNo[k]=avgNo[k]/(double) spList.size();
                     }
                     data.noSwitches.put(key, avgNo);
                   //  Logger.global.info("done sampling "+index);
                   }catch(Exception exc){
                       exc.printStackTrace();
                   }
               }
          // Logger.global.info("  done calculating best path sampling");
           }
   */    
    
   
   /* private double calculateAverage(List<StatePath> spList, ComparableArray result) {
        Map<List<Comparable>, Integer> m = new HashMap<List<Comparable>, Integer>();
        List<List<Comparable>> samples = new ArrayList<List<Comparable>>();
       for(int k=0; k<spList.size(); k++){
           StatePath sp = spList.get(k);
           int i = sp.getFirstEmissionPos();
           List<Comparable> obj = Arrays.asList((ComparableArray)((PairEmissionState)((GenotypeEmissionState)sp.getFirstState()).innerState).getSampleDecomposition(sp.getFirstEmission(), i));
               samples.add(obj);
               Integer count = m.get(obj);
               m.put(obj, count==null ? 1 : count+1);
       }
       SortedSet<Map.Entry<List<Comparable>, Integer>> set = new TreeSet<Map.Entry<List<Comparable>, Integer>>(entryComp);
       set.addAll(m.entrySet());
       List<Comparable> res =  set.last().getKey();
       for(int k=spList.size()-1; k>=0; k--){
           if(!res.equals(samples.get(k))){
               spList.remove(k);
               samples.remove(k);
            //   logger.info(rem+" "+res);
           }
       }
       int size = spList.size();
       while(spList.size() < 200){
           spList.add(new StatePath(spList.get(
               Constants.nextInt(size))));
       }
      res.toArray(result);
      return ((double) set.last().getValue())/ sum(set);
    }*/
   
    int validate(double[] d){
        double sum=0;
        int maxIndex =0;
        for(int i=0; i<d.length; i++){
            sum+=d[i];
            if(d[i]>d[maxIndex] ) maxIndex = i;
        }
        //if(Math.abs(sum-1.0) > 0.001) throw new RuntimeException("!!");
        return maxIndex;
    }
    private double sum(SortedSet<Entry<List<Comparable>, Integer>> set) {
        double sum=0;
        for(Iterator<Map.Entry<List<Comparable>, Integer>> it = set.iterator(); it.hasNext();){
            sum+=it.next().getValue();
        }
        return sum;
     }
     
    
    private static String toString(List<ComparableArray> poss) {
        StringBuffer sb = new StringBuffer();
        sb.append("{");
        for(int i=0; i<poss.size(); i++){
            sb.append(Arrays.asList(poss).toString());
        }
        sb.append("}");
        return sb.toString();
    }

   
    
   private static ComparableArray makeArray(State[] st){
      // boolean order_known = st[0] instanceof CompoundState;
       List<Comparable>res = new ArrayList<Comparable>();
       for(int i=0; i<st.length; i++){
           if(false && st[i] instanceof CompoundState){
               res.add(makeArray(((CompoundState)st[i]).getMemberStates(true)));
           }
           else if(st[i] instanceof State){
               res.add(st[i].getIndex());
           }
           else{
               res.add(st[i]);
           }
       }
       if(false && st[0] instanceof CompoundState && st.length==3) {
           return new TrioComparableArray(new ComparableArray(res));
       }
       else if(false && st[0] instanceof CompoundState && st.length==2){ //redo this!!!
           return new HalfTrioComparableArray(new ComparableArray(res));
       }
       else return new ComparableArray(res);
   }
    
    public void  calcBestPathViterbi(List models)  throws Exception{
        data.clearViterbi();
        int cntSwitch=0;
        EmissionStateSpace[] stateEm=  new EmissionStateSpace[2];//null;
        if(Constants.transMode(1)==null && Constants.saveStates()){
            stateEm[0] =   Emiss.getStateEmissionStateSpace(new int[] {Constants.numF(0)});
            stateEm[1] =   Emiss.getStateEmissionStateSpace(new int[] {Constants.numF(0), Constants.numF(0)});
           
        
        }
     //   data.reinitialiseData();
    for(Iterator<EmissionState> it =data.dataLvalues(); it.hasNext();){
            try{
               EmissionState dat1 =  it.next();
                String key = dat1.getName();
                int nocop = dat1.noCop();
             //  PIGData dat = (PIGData) data.get(key);
               CompoundMarkovModel hmm =(CompoundMarkovModel) models.get(nocop-1);
               while(hmm instanceof WrappedModel){
                   hmm = ((WrappedModel)hmm).getHMM();
               }
            //   StateIndices dat1 =  new StateScorableObject(dat, hmm);
               DP dp  =(DP) dp_pool.get(data.noCopies(key)-1).getObj(dat1.getName());
               dp.reset(true);
               StatePath sp = dp.searchViterbi();
              
             //  if(sp.size()!=dat.length()) throw new RuntimeException("!!  "+dat.noCopies());
              State[] previous = null;
               PIGData viterbi = SimpleScorableObject.make(dat1.getName(), dat1.length(), stateEm[1], 
            		   (short)((HaplotypeEmissionState)dat1).dataIndex());
              
            for(int i=0; i<sp.size();i++){
                CompoundState st = (CompoundState) sp.getState(i);
               State[] states = hmm.disambiguate(st.getMemberStates(true), previous, i, false);//
                if(i>1){for(int k=0; k<states.length; k++){
                    if(states[k]!=previous[k]) cntSwitch++;
                }
                }
                ComparableArray compa = makeArray(states);
               viterbi.addPoint(i, compa);
                
                 previous = states;
             }
            System.err.println("viterbi "+viterbi);
            data.setViterbi(key,  viterbi);
             dp_pool.get(data.noCopies(key)-1).returnObj(dat1.getName());
          //  pw.println();
            }catch(Exception exc){
                exc.printStackTrace();
            }
    }
       // data.getRecSitesFromViterbi();
        System.err.println("switch count is "+cntSwitch);
       
    }

  /*  public void printStates(int j) throws Exception{
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(directory, "hmm_"+j+"/"+"phased_states.txt"))));
        data.writeFastphase(data.viterbi, data.uncertaintyVit, data.uncertaintyVitPhase,pw,  true , false, null);
               pw.close();
        
    }*/
  


  
    
    /** returns a sublist which match 
    public List<PIGData> filter(ArrayList<Integer[][]> vit_cases, Integer[] sigSites, Integer[] sigStates) {
        List<PIGData> list = new ArrayList<PIGData>();
        outer: for(int k=0; k<this.data.size(); k++){
           Integer[][] vit = vit_cases.get(k);
           boolean match = false;
           inner: for(int k1= 0; k1<vit.length; k1++){
               boolean allMatch = true;
               for(int i=0; i<sigSites.length; i++){
                   int s = sigSites[i];
                   int j = sigStates[i];
                   int j1 = vit[k1][s];
                   allMatch = allMatch && j1==j;
               }
               match = match || allMatch;
           }
           if(match)  list.add((PIGData) data.get(k));
       }
        for(Iterator<PIGData> it = list.iterator(); it.hasNext();){
            it.next().restrictTo(sigSites[0], sigSites[sigStates.length-1]);
        }
       return list;
    }*/
    
  
   
}
