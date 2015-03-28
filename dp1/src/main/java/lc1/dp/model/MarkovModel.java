package lc1.dp.model;

import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.IntegerEmiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.SimpleEmissionStateSpace;
import lc1.dp.states.CompoundState;
import lc1.dp.states.DotState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.stats.IntegerDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;


public abstract class MarkovModel implements Serializable {
    public static final DotState MAGIC = new DotState(new DotState("!", SimpleDistribution.noOffset,SimpleDistribution.noOffset));
    protected String name;
    protected List<State> states = new ArrayList<State>();
   
 public abstract boolean converged();
 //  public abstract double[] getPseudoCountWeights();
  // public abstract EmissionStateSpace getStateEmissionStateSpace();
   public State getState(int j){
       return states.get(j);
   }
   public abstract Object clone(boolean swtch);
  
   int length=0;
// protected 
  protected EmissionStateSpace emstsp = null;
   
   public EmissionStateSpace getStateSpace(){
 	  if(emstsp==null){
 		  
 		  List<Comparable> l= new ArrayList<Comparable>();
 		  List<String>nme = new ArrayList<String>();
 		  for(int j=1; j<this.states.size(); j++){
 			  State state =  states.get(j);
 			  l.add(getEl((EmissionState)state, false));
 			  nme.add(state.getName());
 		  }
 		  emstsp = new SimpleEmissionStateSpace(l, nme);
 		
 		  
 	  }
 	  return emstsp;
   }
   
   public static Comparable getEl(EmissionState st, boolean expand){
	   if(expand && st instanceof CompoundState){
		   EmissionState[] membs =  ((CompoundState)st).getMemberStates(true);
			 Comparable[] compa = new Comparable[membs.length];
			for(int k=0; k<membs.length; k++){
				compa[k] = getEl(membs[k], expand);
			}
			return new ComparableArray(compa);
	   }
	   else{
		   return new IntegerEmiss(st.getIndex());
	   }
   }
  
  
  
   
  
   public final double[][] getCounts(int length, double[][]hittingProb){
	   for(int i=0; i<length; i++){
		   double sum=0;
		   for(int j=1; j<this.modelLength(); j++){
			   PseudoDistribution dist = ((EmissionState)this.states.get(j)).emissions(i);
			   if(dist instanceof IntegerDistribution){
				  double cnt =  ((IntegerDistribution)dist).cnt;
				   hittingProb[i][j] = ((IntegerDistribution)dist).cnt;
			   }
			   else{
				   hittingProb[i][j] = Constants.sum(((SimpleExtendedDistribution)dist).counts);
			   }
			   sum+=hittingProb[i][j];
		   }
		   //System.err.println(" i "+sum);
		   for(int j=1; j<this.modelLength(); j++){
			   hittingProb[i][j] = hittingProb[i][j]/sum;
		   }
	   }
	   return hittingProb;
   }
	   
   public final double[][] getHittingProb(int length, double[][]hittingProb){
       
       
       if(hittingProb==null || length!=this.length || true){
          this.length = length;
        
           for(int j=1; j<this.modelLength(); j++){
               //hittingProb[j] = new Double[length];
               hittingProb[0][j] = this.getTransitionScore(0, j, 0);
           }
           for(int i=1; i<length; i++){
               hittingProb[i][0] = 0.0;
               for(int j=1; j<this.modelLength(); j++){
                   double sum=0;
                   double[] trans = new double[this.modelLength()];
                   double[] hp = new double[this.modelLength()]; 
                   for(int k=1; k<this.modelLength(); k++){
                       trans[k] = this.getTransitionScore(k, j, i);
                       hp[k] = hittingProb[i-1][k];
                       sum+=hp[k]*trans[k];
                   }
                   hittingProb[i][j] = sum;
               }
           }
       }
       if(Constants.CHECK){
           for(int i=0; i<length; i++){
               double sum = 0;
              
               for(int j=0; j<modelLength(); j++){
                  // if(hittingProb[j]==null)continue;
                   sum+=hittingProb[i][j];
               }
               if(Math.abs(1.0-sum)>0.01) throw new RuntimeException("!!! "+sum+" "+i);
           }
       }
       return hittingProb;
   }
   
   public abstract void setPseudoCountWeights(double[][] d);
   
   public MarkovModel(MarkovModel m){
       this(m.getName() , m.noSnps);
    
       try{
        //   Iterator<State> it_ps = pseudo==null ? null : pseudo.states();
       for(Iterator<State> it = m.states(); it.hasNext();){
           State st = (State)it.next();
         //  State st_ps = it_ps==null ? null : (State) it_ps.next();
          // workingDist.put(st, 0.0);
           if(st==m.MAGIC) continue;
           State clone = (State) st.clone();
           clone.setIndex(states.size());
           this.addState(clone);//st.getClass().getConstructor(new Class[] {st.getClass()}).newInstance(new Object[] {st}));
       }
     //  this.emissionStateSpaceDist = new EmissionState[this.emissionStateSpace.size()];
     //  boolean allNull = true;
     
       //if(allNull) emissionStateSpaceDist=null;
       }catch(Exception exc){
           exc.printStackTrace();
           System.exit(0);
       }
   }
   
   public MarkovModel(String name, int noSnps){
      this.name = name;
      states.add(MAGIC);
      MAGIC.setIndex( 0);
      this.noSnps = noSnps;
  //    this.emissionStateSpace = emissionStateSpace;
   }
   
   public MarkovModel(String name, int noSnps, boolean includeStart){
	      this.name = name;
	      if(includeStart){
	    	  states.add(MAGIC);
	    	  MAGIC.setIndex( 0);
	      }
	    
	      this.noSnps = noSnps;
	  //    this.emissionStateSpace = emissionStateSpace;
	   }
   /** list of indices of emission states which have the same emission distribution */
     List<int[]> equivalenceClasses = new ArrayList<int[]>();
      public abstract void addCounts( PseudoDistribution[] transProbs, int i, int numIndiv);
      public State addState(State st){
        //    workingDist.put(st, 0.0);
           if(states.contains(st)) return st;
           else {
               st.setIndex(states.size());
               if(st instanceof EmissionState){
                   equivalenceClasses.add(new int[] {st.getIndex()});
               }
               states.add(st);
               return st;
           }
       }
      
  /* public void append(StringBuffer sb){
       StringWriter sw = new StringWriter();
       this.print(new PrintWriter(sw));
       
   }*/
      
    
      
    /** note - to start in begin state startState =0 
           * only adds emission states
           */
          public State[] emitStatePath(State startState, Class clazz){
              List<State> l = new ArrayList<State>();
              if(clazz.isInstance(startState)) l.add(startState);
              State state = startState;
              int cumLength = -1;
              while(cumLength<0 || state!=MAGIC){
                  try{
                        if(state instanceof EmissionState) cumLength++;
                        if(clazz.isInstance(state))l.add(state);
                        state = sampleState(state, cumLength);
                  }catch(Exception exc){
                      Exception exc1 = new Exception("problem with sampling from state "+state+" at "+state+" "+cumLength);
                      exc1.initCause(exc);
                     exc1.printStackTrace();
                     System.exit(0);
                  }
              }
              return (State[]) l.toArray(new State[0]);
          }
   
    public String getName(){
            return name;
       }
   
   public abstract double getTransitionScore(int from, int to, int positionOfToEmission);
   public void initialiseEmissionCounts(){
       for(Iterator<State> states = this.states(); states.hasNext();){
           State st = states.next();
           if(st instanceof EmissionState){
               ((EmissionState) st).initialiseCounts();
           }
       }
   }
   public void initialiseCounts(){
     initialiseTransitionCounts();
      initialiseEmissionCounts();
     //  if(this.emissionStateSpaceDist!=null){
     //  for(int i=0; i<this.getEmissionStateSpace().size(); i++){
     //      EmissionState st = this.getEmissionStateSpaceDist(i);
     //      if(st!=null)st.initialiseCounts();
      // }
       //}
      //ps = ps*0.8;
   }
   
   static  public double ps = 0.0;
   public abstract void initialiseTransitionCounts();
   
   public int maxAdv(){
       int adv=0;
       for(Iterator<State> it = states.iterator(); it.hasNext();){
           State state = it.next();
           if(state.adv > adv) adv = state.adv;
       }
       return adv;
   }
   
   
   
   public  int modelLength(){
       return states.size();
   }
   public int rateLength(){
	   return states.size();
   }
   
   public  State sampleState(State state, int i ) {
       Iterator<State> it = this.statesOut(state, i).iterator();
       double pr = Constants.rand.nextDouble();
       double sum = 0;
      while( it.hasNext()){
           State entry = it.next();
           sum+=this.getTransitionScore(state.getIndex(), entry.getIndex(), i+entry.adv);
           if(sum>=pr){
               return entry;
           }
       }
       throw new RuntimeException("did not find sample "+sum+" "+pr);
   }
   
   public static double sum(Map<Object, Double> m){
       double sum =0;
       for(Iterator<Double > d = m.values().iterator(); d.hasNext();){
           sum+=d.next();
       }
       return sum;
   }
 protected static Map getMap(Iterator<Entry<Object, Double>> it){
     Map m = new HashMap();
     while(it.hasNext()){
         Entry<Object, Double> entry = it.next();
         m.put(entry.getKey(), entry.getValue());
     }
     return m;
 }
   /** sets all distributions to uniform distribution 
   public void setRandom(double trans, double emiss, boolean restart, boolean lastOnly){
       if(emiss<Double.POSITIVE_INFINITY){
           for(Iterator<State> states = this.states(); states.hasNext();){
       
           State st = states.next();
           if(st instanceof EmissionState && (!lastOnly || !states.hasNext() )){
               ((EmissionState)st).setRandom(emiss, restart);
           }
       }
       }
       if(trans<Double.POSITIVE_INFINITY) this.setRandomTransitions(trans, restart, lastOnly);
   }*/
   
 
   
   
  // public abstract void setRandomTransitions(double u, boolean restart, boolean lastOnly);
   
   public  Iterator<State> states(){
       return states.iterator();
   }
   public  List<State> statesIn(final State to, final int positionOfToEmission) {
       return this.states;
   }
   //end at state to in position i1
  /* public Iterator<State> statesInDefault(final State to, final int positionOfToEmission) {
       return this.states();
       //final int i1 = to.isMagic() ? i+1 :i;
      /* return new Iterator<Entry<Object, Double>>(){
           Iterator<Entry<Object, Double>> it =  workingDist.iterator();
           public boolean hasNext(){
               return it.hasNext();
           }
           public Entry<Object, Double> next(){
               Entry<Object, Double> entry = it.next();
               entry.setValue(getTransitionScore( (State) entry.getKey(),to, positionOfToEmission));
               return entry;
           }
           public void remove(){}
       };*/
   
   
   public   List<State> statesOut(final State from,  final int beforeToEmission){
       return this.states;
   }
   //start at state from in position i1
   /*public Iterator<State> statesOutDefault(final State from,  final int beforeToEmission) {
       //final int i1 = from.isMagic() ? i :i+1;
      return this.states();
     /*  return new Iterator<Entry<Object, Double>>(){
           Iterator<Entry<Object, Double>> it =  workingDist.iterator();
           public boolean hasNext(){
               return it.hasNext();
           }
           public Entry<Object, Double> next(){
               Entry<Object, Double> entry = it.next();
               State to =  (State) entry.getKey();
               entry.setValue(getTransitionScore(from,to, beforeToEmission+to.adv));
               return entry;
           }
           public void remove(){}
       };*/
   
   public void print(PrintWriter pw, List<Integer>cols, int popsize){
       pw.println(name);
       for(Iterator<State>it = states(); it.hasNext();){
             State st = it.next();
             if(st instanceof EmissionState){
                 ((EmissionState)st).print(pw, "State_"+st.getName()+"    ", cols);
                 pw.print("\n");
             }
       }
     //  pw.println("misclassification");
       
       /*  int len = this.getEmissionStateSpace().size();
       for(int i=0; i<len; i++){
         State st = this.getEmissionStateSpaceDist(i);
           if(st==null) continue;
           if(st instanceof EmissionState){
               ((EmissionState)st).print(pw, "State_"+st.getName()+" ", cols);
               pw.print("\n");
           }
       }*/
   }
   public String toString(){
        StringWriter sb = new StringWriter();
        PrintWriter pw = new PrintWriter(sb);
        this.print(pw, null, 100);
        return sb.getBuffer().toString();
  }
   /*public final double totalTransitionEmissionDist(MarkovModel m1){
   //    if(m1==this) throw new RuntimeException("!!");
       double sum=0;
       if(m1.modelLength()!= this.modelLength()) throw new RuntimeException("cannot compare");
       for(int i=0; i<this.states.size(); i++){
           State st1 = states.get(i);
           State st2 = m1.states.get(i);
           if(st1 instanceof EmissionState){
               sum+=((EmissionState)st1).KLDistance((EmissionState)st2);
             
           }
       }
         sum+= this.totalTransitionDistance(m1);
       return sum;
   }*/
   
  
   public abstract void transferCountsToProbs(int index);
  
   
public Iterator<int[]> equivalenceClasses(){
    return this.equivalenceClasses.iterator();
}

   

    
    public void validate(int length) throws Exception{
        for(int j=0; j<this.states.size(); j++){
            State st = this.states.get(j);
            if(st instanceof EmissionState){
                ((EmissionState)st).validate();
            }
        }
        //if(!(this instanceof CachedHMM) ){
            this.validateTrans(length);
        //}
      
    }
    
    
   public void validateTrans(int length) throws Exception{
        int start = 0; int end = 1;
        for(int i=-1; i<length-1; i++){
            if(i>=0){
                start = 1;
                end = this.modelLength();
            }
           
            for(int k=start; k<end; k++){
                double sum=0;
                Double[] d = new Double[this.modelLength()];
                for(int j=0; j< modelLength(); j++){
                    int adv = states.get(j).adv;
                    d[j] = this.getTransitionScore(k, j, adv+i) ;
                    sum+=d[j];
                }
               if(Math.abs(1.0-sum)>SimpleDistribution.tolerance){
                   this.validateTransAt(i);
                   int adv = states.get(1).adv;
                   this.getTransitionScore(k, 1, adv+i);
                   throw new Exception(sum+" at  "+i+" "+k+" "+this.getClass());
               }
            }
        }
    }
  // public abstract double getTransitionScorePseudo(int k, int j, int i);

    protected void validateTransAt(int i) {
    // TODO Auto-generated method stub
    
}
  
    
    public final Integer noSnps;
    abstract public int[] statesIn(int j, int i) ;
        
    abstract public int[] statesOut(int j, int i) ;
   
   // protected final EmissionStateSpace  emissionStateSpace;
   
 
    
    
    /*stage is training stage 
    public abstract void updateEmissionStateSpaceDist(int stage);*/
    
   /* public EmissionStateSpace getEmissionStateSpace(){
        if(emissionStateSpace==null) throw new RuntimeException("must initialise state space");
        return emissionStateSpace;
    }
    /** returns the index of the output object in the state space
    public int getEmissionStateSpaceIndex(Object element) {
    //    System.err.println(element+" "+stateSpaceToIndex.get(element));
        try{
           
      return this.emissionStateSpace.get(element);
        }catch(Exception exc){
            System.err.println("prob with "+this.getName());
            exc.printStackTrace();
            return -1;
        }
    } */
   /* public  EmissionState getEmissionStateSpaceDistribution(int index) {
        if(emissionStateSpaceDist==null) return null;
       return emissionStateSpaceDist[index];
    }*/
    public abstract boolean trainEmissions() ;
    public String info() {
       return "";
    }
	public State getState(Comparable comparable) {
		return this.getState(((IntegerEmiss)comparable).v.intValue());
	}
    
  

      
}
