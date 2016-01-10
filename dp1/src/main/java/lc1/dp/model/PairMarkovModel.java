package lc1.dp.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.IntegerEmiss;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.State;
import lc1.dp.states.WrappedEmissionState1;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SkewNormal;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

/** A class which consists of two underlying markov models chained together 
 *  Each of the emission states is paired, and each of the non-emitting states (eg start/end) are unpaired
 * @author Lachlan Coin
 *
 */
public    class PairMarkovModel extends CompoundMarkovModel {
   public  MarkovModel[] m1;
   public final int[] num_copies;  //the pattern of models i.e. [0,1,1] means first model 0 then model 1 then model 2 
   public final int[] count_copies; 
   final Class emissionClass;
   @Override
   
   public State getState(Comparable comparable) {
	   ComparableArray compa = (ComparableArray) comparable;
	  outer: for(int i=1; i<this.states.size();i++){
		   PairEmissionState st = (PairEmissionState) states.get(i);
		   EmissionState[] membs = st.getMemberStates(true);
		   for(int j=0; j<compa.size(); j++){
			   Comparable c = compa.get(j);
			   int ind = membs[j].getIndex();
			   int ind1 = ((IntegerEmiss)c).v.intValue();
			   if(ind!=ind1){
				   continue outer;
			   }
		   }
		   return st;
	   }
	   throw new RuntimeException("!!");
//		return this.getState(((IntegerEmiss)comparable).v.intValue());
	}
   
   
 // final  EmissionStateSpace stateEmissionStateSpace;
    @Override
    public void setPseudoCountWeights(double[][] d){
        for(int i=0; i<m1.length; i++){
          m1[i].setPseudoCountWeights(d);
      }
    }
    public boolean trainEmissions(){
        for(int i=0; i<m1.length; i++){
            if(m1[i].trainEmissions())return true;
        }
        return false;
    }
   
   
    
    public void refresh(){
        for(int i=0; i<m1.length; i++){
            if(m1[i] instanceof CompoundMarkovModel){
                ((CompoundMarkovModel)m1[i]).refresh();
            }
        }
    }
    
    /*public  SimpleExtendedDistribution[] probHomoIsHemi(){
        return this.probHomoIsHemi;
    }*/
    
  
    public final State[] disambiguate(State[] state, State[] previous, int indexOfToEmiss, boolean sample) {
      // if(previous==null)
            return state;
     /*  EmissionState[] res = new EmissionState[state.length];
       System.arraycopy(state, 0, res, 0, state.length);
       for(int j=0; j<res.length; j++){
         
          
         
       }
       return res;*/
    }
    
    public final EmissionState disambiguate(EmissionState state, EmissionState previous, int indexOfToEmiss, boolean sample, int j){
    	 MarkovModel hmm_inner = this.getMarkovModel(j);
         if(hmm_inner instanceof ExpandedHMM && state instanceof WrappedEmissionState1){
       	
       	 // if(innerSt instanceof CompoundState && innerPr instanceof CompoundState){
	              return 
	                  ((ExpandedHMM)hmm_inner).disambiguate(
	                    (WrappedEmissionState1) state, previous,
	               indexOfToEmiss, sample);
	           //   res[j] = ((ExpandedHMM)hmm_inner).getCompoundState(inner);
       	 // }
         }
         else return state;
         
    }
    
    
    public void transferCountsToProbs(int index) {
     //   for(int i=0; i<probB.length; i++){
      //      this.probB[i].bMaximistation(this.pseudo1[3], index);
     //   }
        for(int i=0; i<m1.length; i++){
            m1[i].transferCountsToProbs(index);
        }
    }
    
  

 
  public Object clone(boolean swtch) {
      //  if(m1!=m2) throw new RuntimeException("!!");
      MarkovModel[] m = new MarkovModel[this.m1.length];
      for(int i=0; i<m.length; i++){
          m1[i] = (MarkovModel)this.m1[i].clone(swtch);
      }
      return new PairMarkovModel(m, num_copies, this.emissionClass,true,this.onlyRepeats);
    }
  final boolean onlyRepeats;
  

   
    
   
   
  /*  public void updateEmissionStateSpaceDist(int stage){
        
        if(this.m1[0].emissionStateSpaceDist!=null && this.m1[0].emissionStateSpaceDist[0]!=null){
            emissionStateSpaceDist = new EmissionState[this.emissionStateSpace.size()];
            for(int i=0; i<emissionStateSpace.size(); i++){
                ComparableArray compa = (ComparableArray)this.emissionStateSpace.get(i);
                List l1 = new ArrayList<EmissionState>();
           //     boolean allOne = true;
                for(int k=0; k<compa.size(); k++){
                    EmissionState st = 
                        this.m1[num_copies[k]].getEmissionStateSpaceDistribution(m1[num_copies[k]].getEmissionStateSpaceIndex(compa.get(k)));
                //    if(((HaplotypeEmissionState)st).emissions[0].probs[k]!=1.0) allOne = false;
                    if(st==null) {
                        throw new RuntimeException("!! is null "+k);
                    }
                    l1.add(st);
                }
             //   if(!allOne){
            //    *** need to fix this!!!!
                try{
                   CompoundState res =       (PairEmissionState) emissionClass.getConstructor(new Class[] {List.class, PairMarkovModel.class, boolean.class}).newInstance(new Object[] {l1, this, true});
                    //   this.constructPair(l1)
                    //   new PairEmissionState(l1, this, true); ///*** rethink this for TrioEmissionState!!!
                 //   if(this.probHomoIsHemi!=null) res = new DistributedEmissionState(res, probHomoIsHemi);
                    emissionStateSpaceDist[i] =res; 
                }catch(Exception exc){
                    exc.printStackTrace();
                }
               // }
            }
        }
    }*/
    
  /*List<Set<Integer>> homoToHemizygous = new ArrayList<Set<Integer>>();
   
   public Set<Integer> getHemizygous(int emissionIndex){
       if(homoToHemizygous==null) return null;
       else return homoToHemizygous.get(emissionIndex);
   }*/
    
 // public IlluminaProbB probB(int i){
 //     return probB[i];
 // }
    
  //  final private 
   final private Map<List<Comparable>, CompoundState> compoundStates = new HashMap<List<Comparable>, CompoundState>();
   
   public CompoundState getCompoundState(State[] res){
       CompoundState state =  compoundStates.get(Arrays.asList(res));
       if(state==null) throw new RuntimeException("!!");
       return state;
   }
   
  
   
  
  
   public final Comparator<CompoundState> getComparator(){
       return new Comparator<CompoundState>(){

        public int compare(CompoundState l1, CompoundState l2) {
                    ComparableArray c1 = l1.getMemberStatesRecursive(false);
                    ComparableArray c2 =  l2.getMemberStatesRecursive(false);
                    List<Comparable> co1 = new ArrayList<Comparable>(c1.elements());
                    List<Comparable> co2 = new ArrayList<Comparable>(c2.elements());
                    Collections.sort(co1);
                    Collections.sort(co2);
                    for(int i=0; i<co1.size(); i++){
                        int res = co1.get(i).compareTo(co2.get(i));
                        if(res!=0) return res;
                    }
                    return 0;
                //    return c1.getComparator().compare(c1, c2);
                 //Comparable[] compStr1 = c1.compString();
                //  Comparable[] compStr2 = c2.compString();
                   
                    //int res = c1.compareTo(c2);
                    //return res; 
//                  return comparison(c1,c2);
        }
           
       };
//       return ((ComparableArray)this.emissionStateSpace.get(0)).getComparator();
   }
   
   static String getName(MarkovModel[] m1, int[] no_copies){
       StringBuffer sb = new StringBuffer("(");
       for(int i=0; i<no_copies.length; i++){
           sb.append(m1[no_copies[i]].getName()+(i<no_copies.length-1 ? "X":""));
       }
       sb.append(")");
       return sb.toString();
   }

   public PairMarkovModel(final MarkovModel[] m1, int[] no_copies, Class emClass, final boolean decomp,
		   boolean onlyRepeats
		   ){
        super(getName(m1, no_copies),  m1[0].noSnps);
        this.onlyRepeats = onlyRepeats;
     //   this.probHomoIsHemi =no_copies.length==1 ? null :  probHemiIsHomo;
        int noSources = Constants.format.length;
        
        this.emissionClass = emClass;
        this.num_copies = no_copies;
        this.count_copies = new int[m1.length];
        Arrays.fill(count_copies, 0);
        for(int i=0; i<num_copies.length; i++){
            count_copies[num_copies[i]]++;
        }
        this.from1 = new int[no_copies.length];
        this.to1 = new int[no_copies.length];
        this.m1 = m1;
      //  final  Comparator<CompoundState> comp = getComparator();
       // final SortedMap<CompoundState, List<CompoundState>> m = new TreeMap<CompoundState, List<CompoundState>>(comp);
        this.emstsp = getStateSpace(m1, count_copies, onlyRepeats);
        {
        	
        
      
     
        
      
            //this just gets all of the pair states ready
        	
        	for(int i=0; i<emstsp.haploSize(); i++){
            	ComparableArray compa = (ComparableArray) emstsp.get(i);
            	Comparable[] states = new Comparable[compa.size()];
            	for(int k=0; k<states.length; k++){
            		states[k] = m1[no_copies[k]].getState(compa.get(k));
            	}
            	List<Comparable> l = Arrays.asList(states);
            	   CompoundState st = constructPair(l, decomp);
            	 //  List<CompoundState> equiv = m.get(st);
                 //  if(equiv==null){
                 //      m.put(st, equiv = new ArrayList<CompoundState>());
                 //  }
                 //  equiv.add(st);
                   this.addState(st);
                   compoundStates.put(l, st);
            }
             
            /*CopyEnumerator dblIterator = new CopyEnumerator(this.num_copies.length){
                public Iterator getPossibilities(int depth) {
                    return PairMarkovModel.this.m1[num_copies[depth]].states();
                }
                public void doInner(Comparable[] list) {
                 //   if(excl(list)) return;
                   ComparableArray l =new ComparableArray(Arrays.asList(list));
                   CompoundState st = constructPair(l.elements(), decomp);
                  //comp.compare(l, m.firstKey());
                   List<CompoundState> equiv = m.get(st);
                   if(equiv==null){
                       m.put(st, equiv = new ArrayList<CompoundState>());
                   }
                   equiv.add(st);
                   compoundStates.put(l.elements(), st);
                }
                public boolean exclude(Object obj, Object prv){
                    if(obj==MAGIC) return true;
                    else if( ((State)obj)instanceof DotState) throw new RuntimeException("!! "+obj);
                  //  else if(Constants.exclude()){
                   //     if(ob)
                 //   }
                    else return false;
                }
            };
            dblIterator.run();*/
          singleLine =   new PseudoDistribution[m1.length][];
            for(int i1=0; i1<singleLine.length; i1++){
                singleLine[i1] = new PseudoDistribution[m1[i1].modelLength()];
                for(int i2 =0; i2<singleLine[i1].length; i2++){
                    singleLine[i1][i2] = new SimpleExtendedDistribution(m1[i1].modelLength());
                }
                
            }
        }
    
        /*this is to get things in the right order (i.e. duplicates are at the end )
        for(Iterator<List<CompoundState>> it = m.values().iterator(); it.hasNext();){
            List<CompoundState> nxt = it.next();
            this.addState(nxt.get(0));
        }
        for(Iterator<List<CompoundState>> it = m.values().iterator(); it.hasNext();){
            List<CompoundState> nxt = it.next();
            for(int j=1; j<nxt.size(); j++){
                this.addState(nxt.get(j));
            }
        }*/
          //  this.makeEquivalenceClasses(m.values().iterator());
            System.err.println(this.emissionClass);
        in = new int[this.states.size()-1];
        for(int jk=1; jk<states.size(); jk++){
            in[jk-1] = jk;
        }
        Logger.global.info("state space size is "+this.states.size()+" for "+this.getName());
        if(Constants.CHECK){
            try{ 
           
                validate(noSnps);
            }catch(Exception exc){
                exc.printStackTrace();
            };
        }
      //Logger.getAnonymousLogger().info("equiv classes are "+this.equivalenceClasses);
     
      
    }
  
	  
  
 /*  protected void makeEquivalenceClasses(Iterator<List<CompoundState>> it){
       this.equivalenceClasses.clear();
       while( it.hasNext()){
           List<CompoundState> l = it.next();
           int[] val = new int[l.size()];
           for(int i=0; i<l.size(); i++){
               val[i] = l.get(i).getIndex();
           }
           
           this.equivalenceClasses.add(val);
       }
   }*/
  
   
   
    public String toString(){
        StringBuffer sb = new StringBuffer(m1[0].toString());
        for(int i=1; i<m1.length; i++){
            sb.append(m1[i].toString());
        }
        return sb.toString();
    }
    
   public PairMarkovModel(PairMarkovModel m1){
       this(m1.m1, m1.num_copies, m1.emissionClass,  true, m1.onlyRepeats);
   }
    
  // final IlluminaProbB[] probB; //need one for each data type
                                 
   SkewNormal probR2;
   /*public int noCop(List st1){
       int noCop =0;
       for(int i=0; i<st1.size(); i++){
           noCop+=((EmissionState)st1.get(i)).noCop();
       }
       return noCop;
   }*/
   private  PairEmissionState constructPair(List st1, boolean decompose) {
	   
	   
       PairEmissionState res = null;
       try{
        res =   (PairEmissionState) this.emissionClass.getConstructor(new Class[] {List.class, boolean.class}).newInstance(new Object[] {st1,decompose});
       }catch(Exception exc){
    	   exc.printStackTrace();
       }
    return res;
   }
    
  /*  public void setRandomTransitions(double u, boolean restart, boolean lastOnly) {
        for(int i=0; i<m1.length; i++){
            this.m1[i].setRandomTransitions(u, restart, lastOnly);
        }
//       if(m1!=m1)throw new RuntimeException("expecting these to be same");
    }*/

    public void initialiseTransitionCounts() {
        for(int i=0; i<m1.length; i++){
      this.m1[i].initialiseTransitionCounts();
        }
  //    if(m1!=m1) this.m1.initialiseTransitionCounts();
    }
    @Override
    public void initialiseEmissionCounts(){
        super.initialiseEmissionCounts();
       /* for(int i=0 ; i<probB.length; i++){
            this.probB[i].initialiseBCounts();
       
        List<ProbabilityDistribution> l = new ArrayList<ProbabilityDistribution>();
      //  this.probR(null, null,l, null,i);
        for(Iterator<ProbabilityDistribution> it = l.iterator(); it.hasNext();){
            it.next().initialise();
           
            
           
        }
        }*/
    }

   
final int[] from1, to1;
/** index is the index of sequence in which transition ends */
    public  double getTransitionScore(int from, int to, int endIndex){
       updateFrom(from);
       updateTo(to);
        double sc = 1.0;
        for(int j=0; j<from1.length; j++){
            sc*=m1[num_copies[j]].getTransitionScore(from1[j], to1[j], endIndex);
            
        }
      //  if(Constants.CHECK && sc > 1.0){
       // 	throw new RuntimeException("!!");
        //}
       /* if(this.correction!=null){
        	if(to!=0){
        	int ind = 1+endIndex-this.getState(to).adv;
        	double corr = this.correction[ind][from];
        	sc = sc*corr;
        	}
        }*/
        return sc;
    }
    
  //  double[][] correction;
    /*
    @Override
   protected void validateTransAt(int i, int k, double sum) {
    	if(correction==null) correction = new double[this.noSnps][this.modelLength()];
	    // TODO Auto-generated method stub
    	correction[i+1][k] = 1.0/sum;
    	
	}*/
    
 /*   protected  List getStateSpaceEquivalents(Comparable em){
        return ((ComparableArray)em).decompose(true);
    }
    
    public  double getTransitionScorePseudo(int from, int to, int endIndex){
        updateFrom(from);
        updateTo(to);
         double sc = 1.0;
         for(int j=0; j<from1.length; j++){
             sc*=m1[num_copies[j]].getTransitionScorePseudo(from1[j], to1[j], endIndex);
         }
         return sc;
     }*/
    
   boolean updateFrom(int from){
        if(from ==0){
            Arrays.fill(from1, from);
            return true;
        }
        else{
            PairEmissionState pes1 = (PairEmissionState) getState(from);
            for(int i=0; i<this.num_copies.length; i++){
                from1[i] = pes1.getInnerState(i, true).getIndex();
            }
            return false;
        }
    }
   boolean updateTo(int to){
       if(to==0){
           Arrays.fill(to1, to);
          return true;
       }
       else{
           PairEmissionState pes1 = (PairEmissionState)getState(to);
           for(int i=0; i<this.num_copies.length; i++){
               to1[i] = pes1.getInnerState(i, true).getIndex();
           }
           return false;
       }
   }
 //  final StateDistribution[][] map1; //used for adding counts
 //  final StateDistribution[] dist1; //used for adding counts;
   @Override
  public void validate(int length) throws Exception{
       for(int i=0; i<m1.length; i++){
           m1[i].validate(length);
       }
      super.validate(length);
       
  }
   
   public synchronized void transform(PseudoDistribution[] transProbs, int i){
       for(int j=0; j<transProbs.length; j++){
           updateFrom(j);
           double [] counts = transProbs[j].counts();
           for(int j1=0; j1<this.modelLength(); j1++){
              Double val = counts[j1];
           // if(val>Constants.bwThresh()){
              updateTo(j1);
              for(int k=0;k<this.num_copies.length; k++){
                  int ind = this.from1[k];
                  PseudoDistribution distr = singleLine[num_copies[k]][ind];
                //  if(distr==null) {
                //      distr = new StateDistribution(m1[num_copies[k]].modelLength());
                //      singleLine[num_copies[k]][ind] = distr;
                //  }
                  distr.addCount(this.to1[k], val);
              }
            }
          //}
       }
      
   }
   final PseudoDistribution[][]singleLine;
   @Override
    public void addCounts( PseudoDistribution[] transProbs, int i, int numIndiv) {
       // if(Constants.CHECK) validate(Arrays.asList(transProbs), numIndiv);
    for(int i1=0; i1<singleLine.length; i1++){
          PseudoDistribution[] ss = singleLine[i1];
          for(int i2 = 0; i2 < ss.length; i2++){
              ss[i2].initialise();
          }
    }
       
       this.transform(transProbs,  i);
       for(int i1=0; i1<singleLine.length; i1++){
          // if(Constants.CHECK) validate(Arrays.asList(singleLine[i1]), numIndiv*count_copies[i1]);
           m1[i1].addCounts(singleLine[i1],i, numIndiv*count_copies[i1]);
       }
    }
    
    void validate(List<StateDistribution> obs1, double total){
        double sum=0;
        for(Iterator<StateDistribution> obs = obs1.iterator(); obs.hasNext();){
            StateDistribution nxt = obs.next();
            if(nxt==null) continue;
            sum+=nxt.sum();
        }
        if(Math.abs(sum/total-1)>0.001){
            if(Math.abs(sum-total)>0.1) throw new RuntimeException("!! "+sum+" "+total+" "+this.emissionClass);
            for(Iterator<StateDistribution> obs = obs1.iterator(); obs.hasNext();){
                StateDistribution nxt = obs.next();
                if(nxt==null) continue;
                nxt.multiplyValues(sum/total);
            }
            
        }
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

    @Override
    public MarkovModel getMarkovModel(int i) {
       return this.m1[num_copies[i]];
    }

    @Override
    public int noCopies() {
       return num_copies.length;
    }
    /*  @Override
    public EmissionState[] disambiguate(EmissionState[] state1, EmissionState[] prev) {
       EmissionState[] state = new EmissionState[state1.length];
        for(int i=0; i<state.length; i++){
            if(state1[i] instanceof CompoundState){
               
                ((CompoundMarkovModel)this.m1[this.num_copies[i]])
                .disambiguate(((CompoundState)state[i]),
                              ((CompoundState)prev[i]) );
            }
        }
        return state;
        
    }*/
    
    @Override
    public MarkovModel[] getMemberModels() {
        return this.m1;
    }
   /* @Override
    public EmissionStateSpace getStateEmissionStateSpace() {
       return this.stateEmissionStateSpace;
    }*/
    
  
    
}
