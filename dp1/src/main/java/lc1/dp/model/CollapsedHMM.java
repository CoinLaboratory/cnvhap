package lc1.dp.model;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.core.DoubleDoublePool;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.CachedEmissionState;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.State;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

/** the point of this class is to collapse equivalent states from a markov model */
public class CollapsedHMM extends CompoundMarkovModel implements WrappedModel {
    public CompoundMarkovModel hmm;
   /* @Override
    public IlluminaProbB probB(int i) {
       return hmm.probB(i);
    }*/
    public final int[] collapse; //maps from states in hmm to states in this
    public final List<int[]>expand = new ArrayList<int[]>(); //maps from states in this to states in hmm
    //List expand = new ArrayList<int[]>();
    
    public CompoundMarkovModel getHMM(){
        return hmm;
    }
   /* @Override
    public EmissionStateSpace getStateEmissionStateSpace() {
       return hmm.getStateEmissionStateSpace();
    }*/
    public boolean trainEmissions(){
        return hmm.trainEmissions();
    }
  /*  public  SimpleExtendedDistribution[] probHomoIsHemi(){
        return hmm.probHomoIsHemi();
    }*/
    
    public PairMarkovModel unwrapModel(){
        CompoundMarkovModel hmm = this.getHMM();
        while(hmm instanceof WrappedModel){
            hmm = ((WrappedModel)hmm).getHMM();
        }
        return (PairMarkovModel) hmm;
    }
    
    
    public void splitCounts() {
       for(int j=1; j<expand.size(); j++){
           CachedEmissionState ges = (CachedEmissionState) this.getState(j);
           PseudoDistribution[] emissions = ges.emissions;
           int[] ex = expand.get(j);
           PseudoDistribution[][] l = new SimpleExtendedDistribution[ex.length][];
           for(int k=0; k<ex.length; k++){
               l[k] =  ((CachedEmissionState) this.getState(ex[k])).emissions;
           }
           for(int i=0; i<emissions.length; i++){
               double[] count = emissions[i].counts();
               for(int i1=0; i1<count.length; i1++){
                   double cnt = count[i1]/(double)ex.length;
                   for(int k=0; k<ex.length; k++){
                       l[k][i].setCounts(i1, cnt);
                   }
               }
           }
       }
    }
    
    
   int max;
    public CollapsedHMM(CompoundMarkovModel hmm){
        super(hmm.getName()+"f",  hmm.noSnps);
        this.hmm = hmm;
        EmissionStateSpace emstsp = hmm.getStateSpace();
        this.collapse = new int[hmm.states.size()];
        collapse[0] = 0;
        expand.add(new int[]{0});
        boolean allOneLength = true;
       max = 0;
        for(int i1=0; i1<emstsp.defaultList.size(); i1++){
        	int[] hapl = emstsp.getHaploFromHaploPair(i1);
        	int[] hapl1 = new int[hapl.length];
        	for(int k=0; k<hapl1.length; k++){
        		hapl1[k] = hapl[k]+1;
        	}
        	if(hapl1.length>1) allOneLength = false;
        	State st = this.addState(hmm.getState(hapl1[0]));
        	 for(int i=0; i<hapl1.length; i++){
                 collapse[hapl1[i]]= st.getIndex();
             }
        	 if(st.getIndex()!=expand.size()) throw new RuntimeException("!!");
        	 expand.add(hapl1);
        	
        	 if(hapl1.length>max) max = hapl1.length;
        }
      /*  for(Iterator<int[]> it = hmm.equivalenceClasses.iterator(); it.hasNext();){
            int[] equiv = it.next();
            if(equiv.length>1) allOneLength = false;
           State st =  this.addState(hmm.getState(equiv[0])); //note
           for(int i=0; i<equiv.length; i++){
               collapse[equiv[i]]= st.getIndex();
           }
           if(st.getIndex()!=expand.size()) throw new RuntimeException("!!");
           expand.add(equiv);
           if(equiv.length>max) max = equiv.length;
        }*/
        in = new int[this.states.size()-1];
        for(int jk=1; jk<states.size(); jk++){
            in[jk-1] = jk;
        }
      if(allOneLength) throw new RuntimeException("do not need collapsedHMM");
      if(Constants.CHECK){
          try{
              this.validate(this.noSnps);
          }catch(Exception exc){
              exc.printStackTrace();
          }
      }
        Logger.global.info("state space size for collapsed "+this.states.size()+" for "+this.getName());
        this.pool = new DoubleDoublePool(max);
    }

    private void readObject(ObjectInputStream ois) throws IOException, ClassNotFoundException{
    	ois.defaultReadObject();
    	pool = new DoubleDoublePool(max);
    }
    
    PseudoDistribution[] transProbs1;
     transient DoubleDoublePool pool;
    @Override
    public synchronized void addCounts(PseudoDistribution[] transProbs, int i,int numIndiv) {
        if(transProbs1==null){
            transProbs1 = new PseudoDistribution[hmm.states.size()];
           for(int j=0; j<hmm.states.size(); j++){
               transProbs1[j] = new SimpleExtendedDistribution(hmm.states.size());
           }
        }
        else{
            for(int j=0; j<hmm.states.size(); j++){
                transProbs1[j].initialise();
            }
        }
        for(int j=0; j<transProbs.length; j++){
            if(transProbs[j]==null) continue;
            double[] counts = transProbs[j].counts();
            for(int k=0; k<counts.length; k++){
                int[] from = expand.get(j);
              
                double sum=0;
                double value = counts[k];///(double) from.length;
                if(value==0) continue;
                int[] to = expand.get(k);
                double[][] prob = pool.getObj(from.length, to.length);//new double[from.length][];
                for(int m=0; m<from.length; m++){
                     prob[m] = new double[to.length];
                    for(int n=0; n<to.length; n++){
                        prob[m][n] = hmm.getTransitionScore(from[m], to[n], i+this.hmm.getState(k).adv);  //is this index right?
                        sum+=prob[m][n];
                    }
                    
                }
           //     System.err.println(j+" "+k+" "+sum+":"+toString(from)+"->"+toString(to));
                for(int m=0; m<from.length; m++){
                    for(int n=0; n<to.length; n++){
                        transProbs1[from[m]].addCount(to[n], value*(sum==0 ? 1.0 : (prob[m][n]/sum)));
                    }
                }
                pool.returnObj(prob);
            }
        }
      //  System.err.println(sum(transProbs));
      /*  if(Constants.CHECK){
        double diff =sum(transProbs) - sum(transProbs1);
        if(Math.abs(diff) > 0.01){
            throw new RuntimeException("!! "+diff);
        }
        }*/
        
      //  if(transProbs1[0]!=null && transProbs1[0].dist[0]!=transProbs[0].dist[0]) throw new RuntimeException("!!");
        hmm.addCounts(transProbs1, i, numIndiv);
        
    }
    
    public String toString(int[] d){
        Integer[] d1 = new Integer[d.length];
        for(int i=0; i<d.length; i++){
            d1[i] = d[i];
        }
        return Arrays.asList(d1).toString();
    }
    public void refresh(){
        this.hmm.refresh();
    }
    private double sum(StateDistribution[] transProbs){
        double sum=0;
        for(int i=0; i<transProbs.length; i++){
            sum+=transProbs[i].sum();
        }
        return sum;
    }
    private double sum(PseudoDistribution[] transProbs){
        double sum=0;
        for(int i=0; i<transProbs.length; i++){
            if(transProbs[i]!=null)
            sum+=transProbs[i].sum();
        }
        return sum;
    }

    @Override
    public Object clone(boolean swtch) {
        return null;
    }
    
    @Override
    public double getTransitionScore(int from, int to, int positionOfToEmission) {
        int[] frome =expand.get(from);
        int[] toe = expand.get(to);
        
      
        double sc =0;
       // for(int i=0; i<frome.length; i++){
            for(int j=0; j<toe.length; j++){
                sc+=((PairMarkovModel)hmm).getTransitionScore(frome[0], toe[j], positionOfToEmission);
             
            }
       // }
        return sc;///(double) frome.length;
    }
    
  /*  @Override
    public double getTransitionScorePseudo(int from, int to, int positionOfToEmission) {
        int[] frome =expand.get(from);
        int[] toe = expand.get(to);
        double sc =0;
        for(int i=0; i<frome.length; i++){
            for(int j=0; j<toe.length; j++){
                sc+=hmm.getTransitionScorePseudo(frome[i], toe[j], positionOfToEmission);
            }
        }
        return sc/(double) frome.length;
    }*/

  
    @Override
    public void initialiseTransitionCounts() {
       hmm.initialiseTransitionCounts();
        
    }

    @Override
    public void setPseudoCountWeights(double[][] d){
       hmm.setPseudoCountWeights(d);
        
    }

   /* @Override
    public void setRandomTransitions(double u, boolean restart, boolean lastOnly) {
        hmm.setRandomTransitions(u, restart, lastOnly);
        
    }*/



    @Override
    public void transferCountsToProbs(int index) {
        hmm.transferCountsToProbs(index);
        
    }

   

   
    
    final int[] in;// = new int[0];
    final int[] in0 = new int[] {0};

    public int[] statesIn(int j, int i) {
        if(i==0 && j!=0) return in0;
        else  return in;
      
    }
    public int[] statesOut(int j, int i) {
        if(i==noSnps-1 && j!=0) return in0;
        else return in;
    }
    
    /*@Override
    public EmissionStateSpace getEmissionStateSpace(){
      return this.hmm.getEmissionStateSpace();
    }
    @Override
    /** returns the index of the output object in the state space
    public int getEmissionStateSpaceIndex(Object element) {
        return this.hmm.getEmissionStateSpaceIndex(element);
    } */
  /*  @Override
    public  EmissionState getEmissionStateSpaceDistribution(int index) {
        return hmm.getEmissionStateSpaceDistribution(index);
    }
   
    @Override
    public void updateEmissionStateSpaceDist(int stage){
        throw new RuntimeException("!!");
    }

    public EmissionState getEmissionStateSpaceDist(int i){
        return hmm.getEmissionStateSpaceDist(i);
    }*/
    
    public MarkovModel getMarkovModel(int i) {
       return this.hmm.getMarkovModel(i);
    }


    public int noCopies() {
        return this.hmm.noCopies();
    }


    @Override
    public State[] disambiguate(State[] state, State[] previous, int positionOfToEmiss, boolean sample) {
        if(previous==null) return state;
        int[] possibilities = this.expand.get(this.hmm.getCompoundState(state).getIndex());
        if(possibilities.length==1) return state;
        int index = this.hmm.getCompoundState(previous).getIndex();
        double[] prob = PairEmissionState.pool.getObj(possibilities.length);
        double sum=0;
        for(int i=0; i<prob.length; i++){
            prob[i] = this.hmm.getTransitionScore(index,possibilities[i], positionOfToEmiss);
            sum+=prob[i];
        }
        int chosen = sample ? Constants.sample(prob, sum) : Constants.getMax(prob);
        EmissionState[] res =  ((CompoundState)hmm.getState(possibilities[chosen])).getMemberStates(false);
     //   if(!res[0].equals(state[0]))Logger.global.info("sw "+Arrays.asList(state)+" -> "+Arrays.asList(res));
       // else System.err.println("no sw");
        //return res;
        PairEmissionState.pool.returnObj(prob);
       return res;//this.hmm.disambiguate(res, previous, index, sample);
    }
    
    public  EmissionState disambiguate(EmissionState state, EmissionState previous, int i, boolean sample, int j){
    	return this.hmm.disambiguate(state, previous, i, sample,j);
    }

    
    @Override
    public CompoundState getCompoundState(State[] res) {
     return this.hmm.getCompoundState(res);
    }


    @Override
    public MarkovModel[] getMemberModels() {
        return hmm.getMemberModels();
    }

   /* @Override
    public Set<Integer> getHemizygous(int emissionIndex) {
      return hmm.getHemizygous(emissionIndex);
    }*/
    

}
