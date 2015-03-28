/*
 * Created on 17-Aug-2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package lc1.dp.model;

import java.io.PrintWriter;
import java.util.List;

import lc1.dp.states.State;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.StateDistribution;



/**
 * @author lc
 * puts things in a matrix to avoid looking up maps for transitions
 * 
 */
public  abstract class FastMarkovModel extends MarkovModel{
    

    FreeTransitionProbs1 transProbs; 
    
    
  
    
   /* public Object clone(){
        return new FastMarkovModel(this);
    }*/
  
    public  FastMarkovModel(String name){
            super(name, 0);
    }
    public void initialiseTransitions(){
        try{
        this.transProbs = new FreeTransitionProbs1(false, null, this.modelLength());
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    public FastMarkovModel(FastMarkovModel fmm){
        super(fmm);
        throw new RuntimeException("!!");
      //  this.transProb = new TransitionProbs(fmm.transProb, this.states);
    }
    
    
     public void addCounts( StateDistribution[] observed, int i) {
        throw new RuntimeException("!!");
    }
     
    public State addState(State st){
        super.addState(st);
        return st;
    }
     
    
    public void print(PrintWriter pw){
        super.print(pw, null,0);
        transProbs.print(pw, null, 1);
        pw.print("\nTrans: ");
    }
      
      
     /* public void fix(){
          this.initialiseStateSpace();
          for(int i=0; i<transProbs.length; i++){
              ExtendedDistribution outDist = this.transProbs[i];
              for(Iterator<Entry<Object, Double>> it = outDist.probs.iterator(); it.hasNext();){
                  Entry<Object, Double> nxt = it.next();
                statesIn[((Integer)nxt.getKey())].put(i,  nxt.getValue());
              }
          }
      }*/
	
      
	 @Override
    public double getTransitionScore(int from, int to, int fromIndex) {
        return this.transProbs.getTransition(from,to);
    }
     
   
        
    public void validateTransitions(){
       this.transProbs.validate();
    }
     
	public boolean modelChanged = false;
  
  /*  public void initialiseTransitions(){
       this.transProbs = new ExtendedDistribution[states.size()];
       this.statesIn = new ExtendedDistribution[states.size()];
       for(int i=0; i<transProbs.length; i++){
           transProbs[i] = new ExtendedDistribution();
           statesIn[i] = new ExtendedDistribution();
       }
    }*/
    
    public void setTransition(State from, State to, double sc){
            this.transProbs.setTransitionScore(from.getIndex(),to.getIndex(), sc, transProbs.length());
    }
    
	 




@Override
public List<State> statesOut(State k, int i){
throw new RuntimeException("!!");
}


@Override
public List<State> statesIn(State k, int i){
 throw new RuntimeException("!!");
}




	    


/*public void setTransitions(Distribution[] counts){
     // this.statesOut = counts;
      for(int k=0; k<this.statesOut.length; k++){
         if(!getState(k).trainTransition) continue;
         statesOut[k] =counts[k];
          statesOut[k].validate();
          for(int l1=0; l1<statesOut[k].length(); l1++){
              IntegerDouble id_l = statesOut[k].dist[l1];
              setTransitionIn(k,id_l);
          }
      }
  }
     
  private void setTransitionIn(int k, IntegerDouble id_l){
      int l = id_l.stateId;
      for(int k1 = 0; k1<statesIn[l].length();k1++){
          if(statesIn[l].dist[k1].stateId==k){
              statesIn[l].dist[k1].d = id_l.d;
          }
      }
  }*/




/*transient static final Comparator<State>COMP = new Comparator<State>(){
public int compare(State st1, State st2){
    if(st1.getClassId()==(st2.getClassId())){
        if( st1.getName()==st2.getName()){
           return 0;
        }
        else{
            return st1.getName()< st2.getName() ? -1 :1;
        }
    }
    else{
        return st1.getClassId() < st2.getClassId() ? -1 :1;
    }
}
};*/



}
/** maps state to index */


