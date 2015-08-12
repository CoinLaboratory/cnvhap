
package lc1.dp.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.State;
import lc1.dp.states.WrappedEmissionState1;
import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.Dirichlet;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;



/**
 * @author Lachlan Coin
 */
 public class ExpandedHMM  extends FreeHaplotypeHMM{
    
    // final public SiteTransitions trans;
   
  static final long serialVersionUID = 1;
  static double epsilon = 0.00000001;
  Dirichlet dir;
  public String info(){
      return this.trans.info();
  }
  @Override
  public int rateLength(){
	  return this.m.length;
//	   return states.size();
  }
  public boolean converged(){
	    return trans.converged();
	}
  public void initialiseTransitionCounts() {
	  for(int i=0; i<cop.length; i++){
		  ((CachedHMM)this.cop[i]).refresh();
	  }
	  for(int i=0; i<core.length; i++){
		  this.core[i].initialiseTransitionCounts();
	  }
	   trans.initialiseTransitionCounts();
	}
  public double getTransitionScore(int from, int to, int indexOfToEmission){
		
	    return this.trans.getTransitionScore(from, to, indexOfToEmission);
	}

  
  public void addCounts(PseudoDistribution[] observed, int i, int numIndiv) {
	   
	    if(i+1>=trans.transProbs.length ) return;
	    double dist = i>=0 ? this.trans.loc.get(i+1)-trans.loc.get(i) : 0;
	    trans.transProbs[i+1].addCounts(observed, dist);
	    if(trans.globalTrans!=null && i>0) (trans.globalTrans).addCounts(observed, dist, stateToGroup);
  }
int[] stateToGroup;
  public  void transferTransitionCountsToProbs(int index){
	/*  for(int i=0; i<cop.length; i++){
		 ((CachedHMM)cop[i]).transferEmissionCountsToMemberStates();
	  }*/
	 // this.core.transferTransitionCountsToProbs(index);
	    double[] p = new double[pseudocountWeights.length];
	    double[] p1 = new double[pseudocountWeights.length];
	   // double[] e = new double[pseudocountWeights.length];
	    for(int k=0; k<p.length; k++){
	    	p[k] = this.pseudocountWeights[k][1];
	    	//e[k] = this.pseudocountWeights[k][2];
	    	p1[k] =  this.pseudocountWeights[k][5];
	    }
	   // this.core.transferCountsToProbs(index);
	    trans.transferTransitions(p,p1,  index);
	    if(Constants.CHECK){
	    	trans.validate();
	    }
	    for(int i=0; i<this.core1.length; i++){
			  if(this.num_states.get(i)!=coreLength){
			  AbstractTransitionProbs[] probs =core1[i].trans.transProbs;
			  for(int k=0; k<probs.length; k++){
				  ((TruncatedTransitionProbs)probs[k]).reinitialise();
			  }
			  }
		  }
	    for(int i=0; i<cop.length; i++){
			 ((CachedHMM)cop[i]).refresh();
		  }
	    
	   
	}
  
  public void transferCountsToProbs(int index){
	  
	  this.original.transferEmissionCountsToProbs(index);
	  super.transferCountsToProbs(index);
	 
	  for(int k=0; k<core1.length; k++){
		  this.core1[k].transferEmissionCountsToProbs(index);
	  }
	  for (int j = this.modelLength() - 1; j > 0; j--) { // backwards so
			// we transfer
			// from paired
			// states first
		EmissionState hes = (EmissionState) getState(j);
		hes.refreshSiteEmissions();
		}
	  
  }
  public ExpandedHMM(final ExpandedHMM hmm, boolean swtch) throws Exception{
	  super(hmm);
	  this.coreLength = hmm.coreLength;
	  //  this.special = hmm.special;
	 //   modifyWithData = hmm.modifyWithData;
	 //   this.trans = hmm.trans.clone(swtch);
	    if(Constants.CHECK)try{
	        this.validate(this.noSnps);
	    }
	    catch(Exception exc){
	        exc.printStackTrace();
	    }
	        this.stateToGroup = hmm.stateToGroup;
  }
  public MarkovModel clone(boolean swtch){
	  try{
	    return new ExpandedHMM(this, swtch);
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	  return null;
	}
  FreeHaplotypeHMM[] core;
  List<Integer> num_states = new ArrayList<Integer>(); //how many core states to use for modelling cnv states
  int[] aliasToNumStates;
  FreeHaplotypeHMM[] core1;
  CompoundMarkovModel[] cop;
  int[][] m;
  
  SiteTransitions[] transWithin;
 
 final int coreLength;
 // Map<EmissionState, EmissionState> wrappedStates = new HashMap<EmissionState, EmissionState>();
  // FreeHaplotypeHMM core1;//this is the part of core to be multiplied up
 
 FreeHaplotypeHMM original;
  
 @Override
	public void setPseudoCountWeights(double[][] d) {
		super.setPseudoCountWeights(d);
		this.original.setPseudoCountWeights(d);
		for(int i=0; i<this.core1.length; i++){
			core1[i].setPseudoCountWeights(d);
		}
	}
 
 /* hmm_core is hmm with model between alleles */
  public ExpandedHMM(FreeHaplotypeHMM hmm_original, 
		int[] expand, DataCollection datac){
	 
	  super(hmm_original.getName()+"_expanded", datac.loc.size());
	  boolean allsamesize = true;
	  for(int k=1; k<expand.length; k++){
		  if(expand[k]!=expand[0]) allsamesize = false;
	  }
	//  Constants.gammaRate = Constants.gammaRate1;
	  transWithin = new FreeSiteTrans1[hmm_original.states.size()];
	  List<Integer> no_cop = new ArrayList<Integer>();
	  for(int i=1; i<hmm_original.states.size(); i++){
		  no_cop.add(((EmissionState)hmm_original.states.get(i)).noCop());
	  }
	  this.original = hmm_original;
	 
	  this.cop = new CompoundMarkovModel[hmm_original.states.size()-3];
	 
	  this.aliasToNumStates = new int[hmm_original.states.size()-3];
	  for(int i=0; i<aliasToNumStates.length; i++){
		   
		  int ind = num_states.indexOf(expand[no_cop.indexOf(i+2)]);
		  if(ind<0){
			  ind = num_states.size();
			  num_states.add(expand[no_cop.indexOf(i+2)]);
		  }
		  aliasToNumStates[i] = ind;
	  }
	  this.core = new FreeHaplotypeHMM[1+num_states.size()];
	  this.core[0] =getCoreModel(hmm_original, expand[no_cop.indexOf(1)], datac);
	 this.coreLength = core[0].modelLength()-1;
	  this.core1  = new FreeHaplotypeHMM[num_states.size()];
	 boolean useSame = Constants.useSameModelForAllCN();
	  for(int i=0; i<core1.length; i++){
		  try{
			  this.core[i+1] = !useSame ? getCoreModel(original, num_states.get(i), datac) : core[0];
			  this.core1[i] = 
				  num_states.get(i)==coreLength ? 
						  core[i+1] : 
							   
									
				  extractFromCore( core[i+1] , num_states.get(i));
			//  core1[i].validate(noSnps);
			  }catch(Exception exc){
				  exc.printStackTrace();
			  }
	  }
	  for(int i=0; i<cop.length; i++){
		  boolean onlyPairs = false;
		  if(expand[i+2]== this.coreLength){
				onlyPairs=true;
			  }
		  cop[i] = makeHMM(i+2, core1[aliasToNumStates[i]], onlyPairs);
	  }
		 m = new int[hmm_original.states.size()][];
			m[0]  = new int[] {0};
	  for(int kk=1; kk<hmm_original.states.size(); kk++){
		  int noCop =((EmissionState) hmm_original.states.get(kk)).noCop();
		  if(noCop==0)
	  {
		  	
		  	//add deletion states
			  EmissionStateSpace emStSp= //Emiss.mergedSpace;
					Emiss.spaceByCN[0];
			  int zero_ind = no_cop.indexOf(0);
			  double[][] probs = new double[Constants.inputDir.length][];
			  for(int k=0; k<probs.length; k++){
				 
			  
				probs[k] =   emStSp.getArray("0");
			  }
			  m[kk] = new int[expand[zero_ind]];
			  for(int i=0; i<expand[zero_ind]; i++){
				  EmissionState sta=null;
				
				  //WARNING NEED TO THINK ABOUT THIS!!!
		          sta = (EmissionState) makeState(states.size()+"", noSnps, probs, true,   emStSp);
		          sta.setIndex(states.size());
		          m[kk][i] = states.size();
		         this.addState(sta);
		         
			  }
			  if(Constants.svnTransM()==0 && expand[zero_ind]== this.coreLength){
					 transWithin[kk] = (core[0]).trans;
				  }
			 
	  }
		  
	 
		  else if(noCop==1)
	  {
		  //add normal states
		  transWithin[kk] = (core[0]).trans;
		  m[kk] = new int[core[0].states.size()-1];
		  for(int i=1; i<core[0].states.size(); i++){
			 
			  m[kk][i-1] = states.size();
			 WrappedEmissionState1 state = new WrappedEmissionState1((EmissionState)core[0].states.get(i));
			//  this.wrappedStates.put((EmissionState)core.states.get(i), state);
			  state.name = states.size()+"";
			  state.setIndex(states.size());
			  this.addState(state);
		  }
	  }
	
		  else if(noCop>1){
			  int i = noCop-2;
			 if(Constants.svnTransM()==0 && expand[noCop]== this.coreLength){
				 transWithin[kk] = (core[0]).trans;
			  }
			 else{
			  transWithin[kk] = ((CachedHMM)cop[i]).trans;
			 }
			  m[kk] = new int[cop[i].states.size()-1];
			  for(int k=1; k<cop[i].states.size(); k++){
				  m[kk][k-1] = states.size();
				  WrappedEmissionState1 state = new WrappedEmissionState1((EmissionState)cop[i].states.get(k));
				 // this.wrappedStates.put((EmissionState)cop[i].states.get(k), state);
				  state.setIndex(states.size());
				  this.addState(state);
			  }
		 
	}
	  }
	  in = new int[this.states.size()-1];
      for(int jk=1; jk<states.size(); jk++){
          in[jk-1] = jk;
      }
      Double[] r = 
          new Double[] {   Constants.expModelIntHotSpot(0)*Constants.probCrossOverBetweenBP,
          Constants.expModelIntHotSpot(1)*Constants.probCrossOverBetweenBP};
      trans = new FreeSiteTrans1(datac.loc, this.states, null, r, noSnps, null,1);
     // trans.globalTrans = 
      trans.globalTrans = hmm_original.trans.globalTrans;
      Logger.global.info("state space size is "+this.states.size()+" for "+this.getName());
      try{
     this.stateToGroup =  allsamesize ?  ((FreeSiteTrans1)trans).initialise2(hmm_original.trans, m, core[0].trans, Constants.expand_init_prior(2)) : ((FreeSiteTrans1)trans).initialise1(hmm_original.trans, m, transWithin, Constants.expand_init_prior(2));
     SimpleExtendedDistribution dist1 = ((SimpleExtendedDistribution)((FreeTransitionProbs1)trans.transProbs[0]).transitionsOut[0]);
     
     int[][] gToS = ((FreeSiteTrans1)trans).groupToState;
 	double[] tmp = ((SimpleExtendedDistribution)((FreeTransitionProbs1)hmm_original.trans.transProbs[0]).transitionsOut[0]).probs;
	//dist1.fill(gToS, tmp);
     dist1.normalise1(gToS,tmp);
      if(Constants.newTrans()) {   
    	  trans.alpha_overall = hmm_original.trans.alpha_overall;
      trans.r = hmm_original.trans.r;
      
      
  }
      }catch(Exception exc){
    	  //this.stateToGroup=null;
    	  exc.printStackTrace();
      }
  }
  
  private FreeHaplotypeHMM extractFromCore(FreeHaplotypeHMM core2, int numFounders)  throws Exception{
	FreeHaplotypeHMM res = new FreeHaplotypeHMM(core2.getName(), 
			new ArrayList<State>(core2.states.subList(1, numFounders+1)),core2.noSnps);
	 
	 
       res.trans = new FreeSiteTrans1((FreeSiteTrans1)core2.trans, numFounders, core2.noSnps);
      
	return res;
}
private FreeHaplotypeHMM getCoreModel(HaplotypeHMM hmm_original, int expand1, DataCollection datac) {
	   Double[] r = 
           new Double[] {   Constants.expModelIntHotSpot[1]*Constants.probCrossOverBetweenBP,
           Constants.expModelIntHotSpot[1]*Constants.probCrossOverBetweenBP};
	   double[] modifyFrac0 = new double[expand1];
	   Arrays.fill(modifyFrac0, 1.0/(double) expand1);
	   int ind =1;
	   for(;ind<hmm_original.states.size();ind++){
		   if(((EmissionState)hmm_original.states.get(ind)).noCop()==1) break;
	   }
	   HaplotypeEmissionState orig = (HaplotypeEmissionState) hmm_original.states.get(ind);
	  FreeHaplotypeHMM res =  new FreeHaplotypeHMM(hmm_original.name+"_expanded_core", expand1, datac.loc.size(), 
			 orig, 
	    		
	    	(Object[])	HaplotypeHMMIterator.getMode(new int[] {0}), datac.loc, r, 
	            null,
	           modifyFrac0,modifyFrac0,
                false,
	            ((DataCollection)datac).probeOnly, datac.numLevels()) ;
	  if(Constants.transMode2[0]=='1'){
		for(int i=1; i<res.trans.transProbs.length; i++){
			res.trans.transProbs[i] = new FreeTransitionProbs1(res.trans.transProbs[i]);
		}
	  }
	  res.trans.globalTrans = null;
	  if(Constants.CHECK){
	  try{
	  res.validate(this.noSnps);
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	  }
	  return res;
}
  
  //ublic FreeHaplotypeHMM(String name, int numFounders, int noSnps,  HaplotypeEmissionState orig,   Object []clazz, List<Integer> locs, 
 // Double[] r, Double[] exp_p1,  double[] rel, double[]rel_start,
  //boolean correlateR, Boolean[] probeOnly, ProbabilityDistribution[]  numLevels)  {

private static CompoundMarkovModel makeHMM(int no_cop, FreeHaplotypeHMM hapHMM, boolean onlyPairs){
      CompoundMarkovModel result = null;
    
          result =  makeBasicModel(no_cop, hapHMM, onlyPairs);
      if(Constants.fast() && result.noCopies()>1 && hapHMM.modelLength()>2 && !result.allOneLength()) {
    	  result = new CollapsedHMM(result);
      }
     // return result;
      Logger.global.info(" mem is "+Runtime.getRuntime().freeMemory());
      if(Constants.cache()){
          
          CompoundMarkovModel res1 =  new CachedHMM(result);
          Logger.global.info(" mem is "+Runtime.getRuntime().freeMemory());
          return res1;
      }
      else return result;
  }
  private static CompoundMarkovModel makeBasicModel(int no_cop, MarkovModel toCopy, boolean onlyPairs){
      int[] no_copies = new int[no_cop];
      MarkovModel[] m = new MarkovModel[] {toCopy};
      for(int i=0; i<no_copies.length; i++){
        //  m[i] = toCopy;
          no_copies[i] = 0;
      }
  
        return  new PairMarkovModel(m, no_copies, PairEmissionState.class,   true, 
        		onlyPairs );
      
     
  }
  public EmissionState makeState(String st, int noSnps, double[][] init, boolean fixed,  
			
			 EmissionStateSpace emissionStateSpace){
	  
	  int[] index = new int[init.length];
		Integer[] cn = new Integer[init.length];
		for(int j=0; j<init.length; j++){
			index[j] = Constants.getMax(init[j]);
		    cn[j] = emissionStateSpace.getCN(index[j]);
		  inner: for(int k=0; k<init[j].length; k++){
			if(init[j][k]>0 && emissionStateSpace.getCN(k)!=cn[j].intValue()){
				cn[j] = null;
				break inner;
			}
		  }
		}
		
	  
			//int index = Constants.getMax(init);
	    	return new  HaplotypeEmissionState(st, noSnps, Constants.u_global(0)[0], init,  emissionStateSpace,  
	            		cn, null);
	   
	  }
public EmissionState getCompoundState(EmissionState[] inner) {
	return this.cop[inner.length-2].getCompoundState(inner);
}
public EmissionState disambiguate(WrappedEmissionState1 state,
		EmissionState prev,
		int index, boolean sample) {
	// TODO Auto-generated method stub
	
	if(state.inner instanceof CompoundState){
		
		EmissionState[] memberStates = ((CompoundState)state.inner).getMemberStates(true);
		CompoundMarkovModel model = cop[memberStates.length-2];
		State[] prevStates ;
		 if(prev instanceof CompoundState) {	
			
			prevStates = ((CompoundState) prev).getMemberStates(true);
			
		 }
		 else{
			prevStates = new State[]{prev};
		/*	prevStates[0] = memberStates[0];
			for(int i=1; i<prevStates.length; i++){
				prevStates[i] = this.core.MAGIC;
			}*/
		 }
		
		 if(prevStates.length>memberStates.length){
				State[] prevStates1 = prevStates;
				prevStates = new State[memberStates.length];
				Arrays.fill(prevStates,core[0].MAGIC);
				int len = Math.min(prevStates.length, prevStates1.length);
				System.arraycopy(prevStates1, 0, prevStates, 0, len);
				////NOTE : when it goes from more to fewer states, we are only phasing relative to first n (i.e. assuming
				/// that we go to the original chroms
	     }
	 
	 if(prevStates.length < memberStates.length){
		CompoundMarkovModel cmm =  ((CachedHMM)model).innerModel;
		if(cmm instanceof PairMarkovModel) return state.inner;
		 CollapsedHMM hmm1 = //model.
			 (CollapsedHMM)cmm;
			
		 int[] possibilities = hmm1.expand.get(hmm1.hmm.getCompoundState(memberStates).getIndex());
		 if(possibilities.length==1){
			 return state.inner;
		 
		 }
		 else{
			  double[] prob = PairEmissionState.pool.getObj(possibilities.length);
			 if(prevStates.length==1){
				 for(int i=0; i<possibilities.length; i++){
					 EmissionState[] res =  ((CompoundState)hmm1.hmm.getState(possibilities[i])).getMemberStates(false);
					 prob[i] = core[0].getTransitionScore(prevStates[0].getIndex(), res[0].getIndex(), index);
				 }
			 }
			 else{
				 CollapsedHMM hmm2 = ((CollapsedHMM)((CachedHMM)cop[prevStates.length-2]).innerModel);
				 State[] memb1 = new State[prevStates.length];
				 for(int i=0; i<possibilities.length; i++){
					 EmissionState[] res =  ((CompoundState)hmm1.hmm.getState(possibilities[i])).getMemberStates(false);
					 System.arraycopy(res, 0, memb1, 0, memb1.length);
					 int index1 = hmm2.hmm.getCompoundState(prevStates).getIndex();
					 int index2 = hmm2.hmm.getCompoundState(memb1).getIndex();
					 prob[i] =  hmm2.hmm.getTransitionScore(index,index2,index);
				 }
			 }
			 double sum = Constants.sum(prob);
			 int chosen = sample ? Constants.sample(prob, sum) : Constants.getMax(prob);
		        EmissionState[] order =  ((CompoundState)hmm1.hmm.getState(possibilities[chosen])).getMemberStates(false);
		     //   if(!res[0].equals(state[0]))Logger.global.info("sw "+Arrays.asList(state)+" -> "+Arrays.asList(res));
		       // else System.err.println("no sw");
		        //return res;
		        PairEmissionState.pool.returnObj(prob);
		        CompoundState res = model.getCompoundState(order);
		        return res;
		 }
	 }
	
		State[] order = model.disambiguate( memberStates, prevStates, index, sample);
		//EmissionState tmp = order[0];
		//order[0] = order[1];
		//order[1] = tmp;
		                
		CompoundState res = model.getCompoundState(order);
		/*if(((CachedEmissionState)state.inner).innerState.equals(res)){
			return state;
		}
		//EmissionState res1 = this.wrappedStates.get(res);
		else {*/
			return res;
		 }
		 else return state.inner;
		//}
	
 }
 

 /*public EmissionState[] disambiguate(CollapsedHMM hmm, EmissionState[] state, EmissionState[] previous, int positionOfToEmiss, boolean sample) {
     if(previous==null) return state;
     int[] possibilities = hmm.expand.get(hmm.getCompoundState(state).getIndex());
     if(possibilities.length==1) return state;
    // int index = this.hmm.getCompoundState(previous).getIndex();
     double[] prob = PairEmissionState.pool.getObj(possibilities.length);
     double sum=0;
     for(int i=0; i<prob.length; i++){
    	  int[] frome =hmm.expand.get(from);
          int[] toe = hmm.expand.get(possibilities[i]);
          
        
          double sc =0;
         // for(int i=0; i<frome.length; i++){
              for(int j=0; j<toe.length; j++){
                  sc+=hmm.getTransitionScore(frome[0], toe[j], positionOfToEmiss);
               
              }
         // }
          prob[i] = sc;
          sum+=prob[i];
     }
     int chosen = sample ? Constants.sample(prob, sum) : Constants.getMax(prob);
     EmissionState[] res =  ((CompoundState)hmm.getState(possibilities[chosen])).getMemberStates(false);
  //   if(!res[0].equals(state[0]))Logger.global.info("sw "+Arrays.asList(state)+" -> "+Arrays.asList(res));
    // else System.err.println("no sw");
     //return res;
     PairEmissionState.pool.returnObj(prob);
    return res;//this.hmm.disambiguate(res, previous, index, sample);
 }*/
 }
