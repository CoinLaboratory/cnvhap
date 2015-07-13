
package lc1.dp.model;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.Dirichlet;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import lc1.util.CopyEnumerator;



/**
 * @author Lachlan Coin
 */
 public class FreeHaplotypeHMM  extends HaplotypeHMM{
    
     public SiteTransitions trans;
   
  static final long serialVersionUID = 1;
  static double epsilon = 0.00000001;
  Dirichlet dir;
  public String info(){
      return this.trans.info();
  }
 
  
  public void swtch(){
	  for(int i=1; i<this.trans.transProbs.length; i++){
			trans.transProbs[i] = new FreeTransitionProbs1(trans.transProbs[i]);
		}
		this.trans.allFree = true;
  }
  public void expand(){
	
  	int[][] m =  super.expand(Constants.expand());
 
  
  	if(trans instanceof FreeSiteTrans1){
  		try{
  		 ((FreeSiteTrans1)trans).initialise1(m);
  		/* double[] u = new double[stateToGroup.length];
  		u[0] = Constants.switchU();
  		for(int i=0; i<u.length; i++){
  			int num1 = m[stateToGroup[i]].length;
  			if(num1<=1){
  				u[i] = Constants.switchU();
  				
  			}
  			else u[i] = 10;
  		}*/
  		/*if(swtch){
  	  		for(int i=1; i<this.trans.transProbs.length; i++){
  	  			trans.transProbs[i] = new FreeTransitionProbs1(trans.transProbs[i], u);
  	  		}
  	  	}*/
  		}catch(Exception exc){
  			exc.printStackTrace();
  		}
  	}
  	
  }
  
 /* public FreeHaplotypeHMM(DataCollection dc,List<String> names, Object []clazz,
		  Double[] r, Double[] exp_p1,  double[] rel, double[]rel_start
	      ){
	  super(dc,names);
	  List<Integer> locs = dc.loc;
	  if(clazz.length>1){
          trans = new FreeExpSiteTrans(locs, this.states, exp_p1, r, noSnps, clazz);
      }
      else{
          trans = new FreeSiteTrans1(locs, this.states, exp_p1, r, noSnps, clazz[0],0);
      }
    
      this.trans.initialise( this.states,  locs, r, exp_p1, rel, rel_start,Constants.u_global(0)[1]);
   this.trans.transferEquilToStart();
     if(Constants.CHECK){
         try{this.validate(noSnps);
         
         }catch(Exception exc){
             exc.printStackTrace();
             System.exit(0);
         }
     }
  }*/
  /** u[0] is for emissions, u[1] for within group transitions and u[2] for between group transitions */
public FreeHaplotypeHMM(String name, int numFounders, int noSnps,   double[] init,  EmissionStateSpace emStSp, Object []clazz, List<Integer> locs, 
        Double[] r, Double[] exp_p1, char[][] mod, double[] rel, double[]rel_start,
        boolean correlateR, Boolean[] probeOnly, ProbabilityDistribution[]  numLevels)  {
      super(name, numFounders, noSnps,   init, emStSp, mod, locs,  
    		  probeOnly,numLevels);
  
      if(clazz.length>1){
          trans = new FreeExpSiteTrans(locs, this.states, exp_p1, r, noSnps, clazz);
      }
      else{
          trans = new FreeSiteTrans1(locs, this.states, exp_p1, r, noSnps, clazz[0],0);
      }
    
      this.trans.initialise( this.states,  locs, r, exp_p1, rel, rel_start,Constants.u_global(0)[1]);
   this.trans.transferEquilToStart();
     if(Constants.CHECK){
         try{this.validate(noSnps);
         
         }catch(Exception exc){
             exc.printStackTrace();
             System.exit(0);
         }
     }
}

public FreeHaplotypeHMM(String name, int numFounders, int noSnps,  HaplotypeEmissionState orig,   Object []clazz, List<Integer> locs, 
        Double[] r, Double[] exp_p1,  double[] rel, double[]rel_start,
        boolean correlateR, Boolean[] probeOnly, ProbabilityDistribution[]  numLevels)  {
      super(name, numFounders, noSnps,   orig);
  
      if(clazz.length>1){
          trans = new FreeExpSiteTrans(locs, this.states, exp_p1, r, noSnps, clazz);
      }
      else{
          trans = new FreeSiteTrans1(locs, this.states, exp_p1, r, noSnps, clazz[0],0);
      }
    trans.globalTrans = null;
      this.trans.initialise( this.states,  locs, r, exp_p1, rel, rel_start,Constants.expand_init_prior(1));
   
     if(Constants.CHECK){
         try{this.validate(noSnps);
         
         }catch(Exception exc){
             exc.printStackTrace();
             System.exit(0);
         }
     }
}
public FreeHaplotypeHMM(final FreeHaplotypeHMM hmm, boolean swtch){
    super(hmm);
  //  this.special = hmm.special;
 //   modifyWithData = hmm.modifyWithData;
    this.trans = hmm.trans.clone(swtch);
    if(Constants.CHECK)try{
        this.validate(this.noSnps);
    }
    catch(Exception exc){
        exc.printStackTrace();
    }
        
    }
public FreeHaplotypeHMM(final MarkovModel hmm) throws Exception{
    super(hmm);
    /*if(hmm instanceof FreeHaplotypeHMM){
        special = ((FreeHaplotypeHMM)hmm).special;
    }
    else{
        special = new HashSet<Integer>();
    }*/
  //  modifyWithData = false;
//    modifyWithData = hmm.modifyWithData;
  //  this.trans = hmm.trans.clone(swtch);
    this.trans = new FreeSiteTrans1(hmm.noSnps, hmm.modelLength(),
    		SiteTransitions.getCN(hmm.states));//, hmm.trans.exp_p1, hmm.trans.r, noSnps, FreeTransitionProbs1.class);
        
            CopyEnumerator dblIterator = new CopyEnumerator(2){
                public  Iterator getPossibilities(int depth) {
                     return states();
                 }
                 public void doInner(Comparable[] list) {
                     State from = (State) list[0];
                     State to = (State) list[1];
                     for(int i=0; i<noSnps; i++){
                         double sc =  hmm.getTransitionScore(from.getIndex(), to.getIndex(), i);
                         if(i==noSnps-1 && to.getIndex()==0) sc =0;
                         if(sc>0){
                             trans.setTransitionScore(from.getIndex(), to.getIndex(), i,sc, sc);
                         }
                            
                     }
                 }
                 public boolean exclude(Comparable[] list){
                	 return false;
                 }
                 public boolean exclude(Object obj, Object previous){
                     return false;
                 }
                 
             };
             dblIterator.run();
             if(Constants.CHECK){
                 try{
                 this.validate(this.noSnps);
                 }catch(Exception exc){
                     exc.printStackTrace();
                 }
             }
        
    }

/*public FreeHaplotypeHMM(final CachedHMM hmm){
    super(hmm);
    modifyWithData = false;
    this.trans = new FreeSiteTrans1(hmm.trans.loc, this.states, hmm.trans.exp_p1, hmm.trans.r, noSnps, FreeTransitionProbs1.class);
        
    this.trans.initialise(this.states, Constants.u_global(), hmm.trans.loc, hmm.trans.exp_p1, hmm.trans.r, false, null);
            CopyEnumerator dblIterator = new CopyEnumerator(2){
                public  Iterator getPossibilities(int depth) {
                     return states();
                 }
                 public void doInner(Comparable[] list) {
                     State from = (State) list[0];
                     State to = (State) list[1];
                     for(int i=0; i<noSnps; i++){
                         double sc =  hmm.getTransitionScore(from.getIndex(), to.getIndex(), i);
                         if(i==noSnps-1 && to.getIndex()==0) sc =0;
                         if(sc>0){
                             trans.setTransitionScore(from.getIndex(), to.getIndex(), i,sc, sc);
                         }
                            
                     }
                 }
                 public boolean exclude(Object obj, Object previous){
                     return false;
                 }
                 
             };
             dblIterator.run();
             if(Constants.CHECK){
                 try{
                 this.validate(this.noSnps);
                 }catch(Exception exc){
                     exc.printStackTrace();
                 }
             }
        
    }*/


public FreeHaplotypeHMM(String string, int size) {
	super(string,size);
}


public FreeHaplotypeHMM(String name, List<State> subList, Integer noSnps) {
	// TODO Auto-generated constructor stub
super(name, subList, noSnps);
}


public MarkovModel clone(boolean swtch){
    return new FreeHaplotypeHMM(this, swtch);
} 



public boolean converged(){
    return trans.converged();
}


@Override
public void validateTransAt(int to){
    this.trans.validateTransAt(to);
}

public double getTransitionScore(int from, int to, int indexOfToEmission){
	
    return this.trans.getTransitionScore(from, to, indexOfToEmission);
}

public double getTransitionScoreToPaint(int from, int to, int indexOfToEmission){
	
    return this.trans.getTransitionScoreToPaint(from, to, indexOfToEmission);
}
@Override
public void transferCountsToProbs(int index){
    super.transferCountsToProbs(index);
}


@Override
public void validate(int length) throws Exception{
    try{
        super.validate(length);
    }
    catch(Exception exc){
       exc.printStackTrace();
       System.err.println("validating transitions ");
     this.trans.validate();
     System.err.println("validating transitions done");
     System.exit(0);
    }
    
}
/** endindex is sequence index for _end_ state */
  
 
  
  private int[] getRandomTransformation(){
      int len = this.modelLength();
      List<Integer> l = new ArrayList<Integer>(len);
      for(int i=1; i<len; i++){
          l.add(i);
      }
      int[] res = new int[len];
      res[0] = 0;
      for(int i=1; i<len; i++){
         res[i] = l.remove(Constants.nextInt(l.size()));
      }
     
      return res;
  }
  
  /*public void reorderStates(boolean random, boolean emissOnly){
   //   if(true)return;
      if(!((EmissionState)this.getState(1)).isFixed()) return;
     // Double[][] hittingProbs =random ? null:  this.getHittingProb(this.noSnps);
      for(int i=1; i<this.noSnps; i++){
          int[] transformation = 
              random ? getRandomTransformation() : 
             getTransformation((FreeTransitionProbs1)this.trans.transProbs[i]);
          int[] cp = new int[transformation.length];
          for(int j=1; j<transformation.length; j++){
              cp[j] = ((FixedHaplotypeEmissionState)this.getState(j)).getFixedInteger(i);
          }
          for(int j=1; j<transformation.length; j++){
              ((FixedHaplotypeEmissionState)this.getState(transformation[j])).setFixedIndex(i, cp[j]);
          }
          if(!emissOnly){
              {
                  FreeTransitionProbs tp =(FreeTransitionProbs) this.trans.transProbs[i];
                  if(i==0){
                      SimpleExtendedDistribution dist = tp.transitionsOut[0];
                      dist.apply(transformation);
                  }
                  else{
                      for(int j=1; j<tp.transitionsOut.length; j++){
                          tp.transitionsOut[j].apply(transformation);
                      }
                  }
              }
                  
                  if(i<this.noSnps-1){
                      FreeTransitionProbs tp =(FreeTransitionProbs) this.trans.transProbs[i+1];
                      SimpleExtendedDistribution[] tmp = new SimpleExtendedDistribution[tp.transitionsOut.length];
                      System.arraycopy(tp.transitionsOut, 0, tmp, 0, tmp.length);
                      for(int j=0; j<tmp.length; j++){
                          tp.transitionsOut[transformation[j]] = tmp[j];
                      }
                  }
          }
      }
  }*/
  
  class DoubleInteger implements Comparable{
      double d;
      int i;
      DoubleInteger(int i, double d){
          this.i = i;
          this.d = d;
      }
      public String toString(){
          return i+":"+d+" ";
      }
    public int compareTo(Object o) {
      DoubleInteger d1 = (DoubleInteger)o;
      if(d== d1.d) return 0;
      else if(d < d1.d) return +1;
      else return -1;
    }
  }
  private List<DoubleInteger> getList(SimpleExtendedDistribution[] probs, int j){
      List<DoubleInteger> l = new ArrayList<DoubleInteger>();
      for(int i=1; i<probs.length; i++){
          l.add(new DoubleInteger(i, probs[i].probs[j]));
      }
      return l;
  }
 /*private int[] getTransformation(FreeTransitionProbs probs) {
     int len  =probs.transitionsOut.length;
     int[] res = new int[len];
     Set<Integer> done = new HashSet<Integer>();
    
     for(int i=1; i<len; i++){
         List<DoubleInteger> l = getList(probs.transitionsOut, i);
         Collections.sort(l);
        // probs.transitionsOut[i].probs;
         int j=0;
         while(true){
             if(done.contains(l.get(j).i)){
                 j++;
                 continue;
             }
             else  break;
         }
         done.add(l.get(j).i);
         res[i] = l.get(j).i;// = i;
     }
    return res;
}*/
private static int[] getTransformation(Double[] doubles) {
     int len  = doubles.length;
    Double[] copy = new Double[len];
   
    int[] res = new int[len];
    System.arraycopy(doubles, 0, copy, 0, len);
    Arrays.sort(copy);
    Set<Integer> done = new HashSet<Integer>();
    for(int i=0; i<len; i++){
        Double d = doubles[i];
        int j=0;
        while(true){
            if(done.contains(j)){
                j++;
                continue;
            }
            else if(copy[j]==d) break;
            j++;
        }
        done.add(j);
        res[i] = j;
    }
    if(done.size()!=doubles.length) throw new RuntimeException("!!");
    if(res[0]!=0) throw new RuntimeException("!!");
    for(int i=1; i<len; i++){
        res[i]  = (len-1) - res[i]+1;
    }
    return res;
}
 public void print(final PrintWriter sb, List<Integer>cols, int popsize){
      super.print(sb, cols, popsize);
   //   this.trans.print(sb, sbS, states, hittingProb)
     if(false){ StringBuffer sb1= new StringBuffer();
      Double[] u= new Double[noSnps+1];
      for(int i=0; i<this.noSnps; i++){
       
          sb1.append("%8.3g ");
          u[i] = pseudocountWeights[0][0];
      }
      final String sbS = sb1.toString();
      sb.println("transitions");
     double[][] hittingProb = new double[length][this.modelLength()];
     
     this.trans.print(sb,sbS, states, this.getHittingProb(this.noSnps, hittingProb));
     }
     
  }
  
static Double zero = new Double(0);

public void addCounts(PseudoDistribution[] observed, int i, int numIndiv) {
   
    if(i+1>=trans.transProbs.length ) return;
    double dist = i>=0 ? this.trans.loc.get(i+1)-trans.loc.get(i) : 0;
    if(dist==0 && i>=0){
    	Logger.global.info("warning dist is zero for "+i);
    }
    trans.transProbs[i+1].addCounts(observed, dist);
   if(trans.globalTrans!=null && i>0) (trans.globalTrans).addCounts(observed, dist);
}


public void initialiseTransitionCounts() {
   trans.initialiseTransitionCounts();
}

public  void transferTransitionCountsToProbs(int index){
	 double[] p = new double[pseudocountWeights.length];
	 double[] p1 = new double[pseudocountWeights.length];
	  //  double[] e = new double[pseudocountWeights.length];
	    for(int k=0; k<p.length; k++){
	    	p[k] = this.pseudocountWeights[k][1];
	    	p1[k] = this.pseudocountWeights[k][5];
	    //	e[k] = this.pseudocountWeights[k][2];
	    }
    trans.transferTransitions(p, p1,  index);
    if(Constants.CHECK){
    	trans.validate();
    }
   //trans.transferTransitions(this.pseudocountWeights[1], this.pseudocountWeights[2]);
}
//final boolean modifyWithData;



  
}
