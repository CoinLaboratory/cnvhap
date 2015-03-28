/**
 * 
 */
package lc1.dp.model;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.Callable;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.states.AbstractCachedEmissionState;
import lc1.dp.states.CachedEmissionState;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.Dirichlet;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;
import lc1.util.CopyEnumerator;

public class CachedHMM extends CompoundMarkovModel implements WrappedModel{
    public CompoundMarkovModel innerModel;
    final public SiteTransitions trans;
    
  /*  public  SimpleExtendedDistribution[] probHomoIsHemi(){
        return innerModel.probHomoIsHemi();
    }*/
    public CompoundMarkovModel getHMM(){
        return innerModel;
    }
    public PairMarkovModel unwrapModel(){
        CompoundMarkovModel hmm = innerModel;
        while(hmm instanceof WrappedModel){
            hmm = ((WrappedModel)hmm).getHMM();
        }
        return (PairMarkovModel) hmm;
    }
   @Override 
   public boolean trainEmissions(){
	   if(Constants.sum(Constants.numIt)==0) return false;
       return innerModel.trainEmissions();
   }
    transient List tasks;
    public void transferEmissionCountsToMemberStates(){
        if(trainEmissions()){
            
            for(int j=1; j<this.modelLength(); j++){
                final int j1 = j;
                tasks.set(j-1,  new Callable(){
                    public Object call(){
                        try{
                        	// if(j1==5){
                            // 	System.err.println("p1");
                             //}
                        ((AbstractCachedEmissionState)getState(j1)).transferCountsToMemberStates();
                       
                        }catch(Exception exc){
                            exc.printStackTrace();
                        }
                        return null;
                    }
                });
            }
            try{
            BaumWelchTrainer.involeTasks(tasks, false);
            }catch(Exception exc){
                exc.printStackTrace();
            }
        }
       // if(Constants.modelInaccurateData()){
     /*  for(int i=0; i<this.getEmissionStateSpace().size(); i++){
                CachedEmissionState st = ((CachedEmissionState)this.getEmissionStateSpaceDist(i));
             if(st!=null){
                st.transferCountsToMemberStates();
             }
        }*/
        //}
    }
    
 //  boolean fixed;
 
    
    public CachedHMM(final CompoundMarkovModel innerModel){
        super(innerModel.getName()+"x"+innerModel.getName(),  innerModel.noSnps);
       
        this.innerModel = innerModel;
       
        
    //   fixed =
        for(int j=1; j<innerModel.modelLength(); j++){
           EmissionState newState;
           CompoundState state_j =   (CompoundState)innerModel.getState(j);
          /*  if( state_j.isFixed() ){
                newState =  new FixedCachedEmissionState(
                        (CompoundState)innerModel.getState(j));
            
            }
            else{*/
                newState = new CachedEmissionState(
                      state_j, 
                              state_j.getEmissionStateSpace().size());
            //}
            this.addState(newState);            
          
                                           
                 
        }
        trans = new FreeSiteTrans1(this.states, noSnps, FreeTransitionProbs1.class);
        trans.cached = true;
       // double[] d=  new double[] {1e6, 1e6, 1e6};
        double[] rel = new double[this.states.size()];
        Arrays.fill(rel, 1.0 / (double)(states.size()-1));
        rel[0] = 0;
        try{
            trans.transProbs[0] =new FreeTransitionProbs1(true, new Dirichlet(rel, 1e10), states.size());
            
        this.trans.initialise(rel, 0,Constants.u_global(0)[1]);    ///CHECK FOLLOWING LINE - BE CAREFUL!!!
       
        }catch(Exception exc){
            exc.printStackTrace();
        }
       /* for(int i=0; i<this.trans.transProbs.length; i++){
           PseudoDistribution[] dist = ((FreeTransitionProbs1) this.trans.transProbs[i]).transitionsOut;
           for(int i1=0; i1<dist.length; i1++){
               if(dist[i1]==null) continue;
               ((SimpleExtendedDistribution) dist[i1]).counts = null;
           }
        }*/
        this.tasks = Arrays.asList(new Callable[this.modelLength()-1]);
        
    }
    
    private void readObject(ObjectInputStream ois) throws IOException, ClassNotFoundException{
    	ois.defaultReadObject();
    	this.tasks = Arrays.asList(new Callable[this.modelLength()-1]);
    }
   public boolean needInit(){
       return this.trans.transProbs[0]==null;
   }
    
    private void refreshSiteTransitions(){
         CopyEnumerator dblIterator = new CopyEnumerator(2){
            public  Iterator getPossibilities(int depth) {
                 return states();
             }
             public void doInner(Comparable[] list) {
                 State from = (State) list[0];
                 State to = (State) list[1];
                 for(int i=0; i<noSnps; i++){
                     double sc = innerModel.getTransitionScore(from.getIndex(), to.getIndex(), i);
                    // if(sc>0){
                       //  if(i==noSnps-1 || i==0) Logger.global.info(from+"->"+to+" "+sc+" "+i);
                   /*  if(Constants.CHECK && Double.isNaN(sc)){
                    	 throw new RuntimeException("!!");
                     }*/
                         trans.setTransitionScore(from.getIndex(), to.getIndex(), i,sc, sc);
                    // }
                 }
             }
             public boolean exclude(Comparable[] list){
            	 return false;
             }
             public boolean exclude(Object obj, Object previous){
                 return false;
//                 if(!(obj instanceof EmissionState)) return true;
//                      else return false;
             }
             
         };
         dblIterator.run();
         if(Constants.CHECK){
         try{
         this.validate(this.noSnps);
         }catch(Exception exc){
             try{
             innerModel.validate(this.noSnps);
             }catch(Exception exc1){
                 exc1.printStackTrace();
             }
             this.trans.transProbs[0].validate();
             exc.printStackTrace();
             System.exit(0);
         }
         }
    }
    
    public Object clone(boolean swtch) {
      //  if(m1!=m2) throw new RuntimeException("!!");
      CompoundMarkovModel m = (CompoundMarkovModel)this.innerModel.clone(swtch);
      return new CachedHMM(m);
    }

   

@Override
public  double getTransitionScore(int from, int to, int indexOfToEmission){
    return this.trans.getTransitionScore(from, to, indexOfToEmission);
}

/*@Override
public  double getTransitionScorePseudo(int from, int to, int indexOfToEmission){
    return this.trans.getTransitionScorePseudo(from, to, indexOfToEmission);
}*/
@Override
public void transferCountsToProbs(int index) {
    innerModel.transferCountsToProbs( index);
  // this.probHomoIsHemi.transfer(this.getPseudoCountWeights()[1]);
}
@Override
public void initialiseEmissionCounts(){
    super.initialiseEmissionCounts();
   /* for(int i=0; i<Constants.format().length; i++){
      //  this.probB(i).initialiseBCounts();
    
    
    List<ProbabilityDistribution> s = new ArrayList<ProbabilityDistribution>();
 //   this.probR(null, null, s, null, i);
    for(Iterator<ProbabilityDistribution> it = s.iterator(); it.hasNext();){
        ProbabilityDistribution dist = it.next();
      if(dist!=null)
           dist.initialise();
       
    }
    }*/
}
public void refresh(){
    this.initialiseEmissionCounts();
    this.innerModel.refresh();
   // BaumWelchTrainer.t[8]=System.currentTimeMillis();
    this.refreshSiteTransitions();
   
    //BaumWelchTrainer.t[9]=System.currentTimeMillis();
    if(trainEmissions() ){
        for(int j=1; j<this.modelLength(); j++){
           final  State st = this.getState(j);
            this.tasks.set(j-1,
                    new Callable(){
               
                    public Object call(){
                        ((AbstractCachedEmissionState)st).refreshSiteEmissions();
                        return null;
                    }
            });
            
        }
        try{
        BaumWelchTrainer.involeTasks(tasks, true); //warning this should always be true until we sought out why this causes a problem if done sequentially
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
  //  BaumWelchTrainer.t[10]=System.currentTimeMillis();
}

@Override
public void addCounts(PseudoDistribution[] transProbs, int i, int numIndiv) {
  
   innerModel.addCounts(transProbs, i, numIndiv);
    
}
/*
@Override
public void validateTrans(int length){
 this.innerModel.validateTrans(length);   
}*/


@Override
public void initialiseTransitionCounts() {
   this.trans.initialiseTransitionCounts();
    innerModel.initialiseTransitionCounts();
    
}

@Override
public void setPseudoCountWeights(double[][] d) {
  innerModel.setPseudoCountWeights(d);
    
}

/*@Override
public void setRandomTransitions(double u, boolean restart, boolean lastOnly) {
 innerModel.setRandomTransitions(u, restart, lastOnly);
    
}*/

@Override
public int[] statesIn(int j, int i) {
    return innerModel.statesIn(j, i);
}

@Override
public int[] statesOut(int j, int i) {
  return innerModel.statesOut(j, i);
}



   
/*@Override
public EmissionStateSpace getEmissionStateSpace(){
  return this.innerModel.getEmissionStateSpace();
}
@Override
/** returns the index of the output object in the state space 
public int getEmissionStateSpaceIndex(Object element) {
return this.innerModel.getEmissionStateSpaceIndex(element);
}*/



/*@Override
public void updateEmissionStateSpaceDist(int stage){
    throw new RuntimeException("!!");
   // this.innerModel.updateEmissionStateSpaceDist();
}*/

/*@Override
public  EmissionState getEmissionStateSpaceDistribution(int index) {
    return this.emissionStateSpaceDist[index];
}*/

public MarkovModel getMarkovModel(int i) {
   
    return innerModel.getMarkovModel(i);
}

public int noCopies() {
   return innerModel.noCopies();
}


@Override
public State[] disambiguate(State[] memberStates, State[] prev, int i, boolean sample) {
    return this.innerModel.disambiguate(memberStates, prev, i, sample);
}

public  EmissionState disambiguate(EmissionState state, EmissionState previous, int i, boolean sample, int j){
	return this.innerModel.disambiguate(state, previous, i, sample,j);
}

@Override
public CompoundState getCompoundState(State[] res) {
    return this.innerModel.getCompoundState(res);
 //throw new RuntimeException("!!");
}


@Override
public MarkovModel[] getMemberModels() {
  return innerModel.getMemberModels();
}
/*@Override
public EmissionStateSpace getStateEmissionStateSpace() {
   return innerModel.getStateEmissionStateSpace();
}*/
/*@Override
public IlluminaProbB probB(int i) {
   return innerModel.probB(i);
}*/


   
}