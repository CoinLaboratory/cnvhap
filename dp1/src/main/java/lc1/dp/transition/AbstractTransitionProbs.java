package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Collection;

import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleDistribution;

public abstract class AbstractTransitionProbs implements Serializable{

  //  public abstract void transferProbsToPseudo();
public void validate(boolean start){
	if(start){
		 double sum=0;
         for(int j=0; j<this.noStates(); j++){
             sum+=this.getTransition(0, j);
         }
         if(Math.abs(sum-1.0) >SimpleDistribution.tolerance) throw new RuntimeException("!!");
         
	}
	else{
	  for(int i=1; i<this.noStates(); i++){
          double sum=0;
          for(int j=0; j<this.noStates(); j++){
              sum+=this.getTransition(i, j);
          }
          if(Math.abs(sum-1.0) >SimpleDistribution.tolerance) throw new RuntimeException("!!");
          
      }
	}
}

public double  transferQ(double[] ds,double pseudoAlpha, double pseudoRate,  MatrixExp initial, int i, double distance, int it) {
	double l = this.transferAlpha(ds, null, i);
	l+= this.transfer(ds, null, i);
	return l;
}

public abstract void validate();

    //public abstract void addCounts(StateDistribution[] observed
         //   );
   public  static String transform(double prob, double dist){
        return String.format("%5.3g ", new Object[] {-Math.log(prob)/dist});
    }
    public abstract double getTransition(int from, int to);
    
    public double getTransitionCount(int group, int i) {
		throw new RuntimeException("!!");
	}
    
    
    public double getTransitionToPaint(int from, int to) {
		return getTransition(from, to);
	}
    /** index allows flexibility of having multiple models */
    public  double getTransition(int index, int from, int to){
        throw new RuntimeException("!!");
        //return this.getTransition(from, to); 
    }

  //  public abstract double getTransitionPseudo(int from, int to);

    /** initialises counts as pseudo-counts distribution multiplied by pseudocount */
    public abstract void initialiseCounts(boolean start, boolean end);

   
   // public abstract void transfer(double pseudoTrans, double pseudoExp);
   // public abstract void transfer(double[] pseudoTrans, double[] pseudoExp);
    public abstract Collection getDistributions();


   // public abstract void validate();


    public final void addCounts(PseudoDistribution[] observed) {
        int no_states = observed.length;
        for(int j=0; j<no_states; j++){
            PseudoDistribution dist1 =observed[j];
            int st = j;
            double[] counts = dist1.counts();
            if(dist1==null) continue;
            for(int j1 = 0; j1<no_states; j1++){
                int state = j1;
                Double val =counts[state];
                if(val==0) continue;
               //  if(st!=0 && state==5 && val>30){
                 //   Logger.global.info("h");
               // }
                addCount(st, state, val);
            }
            
        }
    }
    
    public final void addCounts(PseudoDistribution[] observed, double dist) {
        int no_states = observed.length;
        for(int j=0; j<no_states; j++){
            PseudoDistribution dist1 =observed[j];
            int st = j;
            double[] counts = dist1.counts();
            if(dist1==null) continue;
            for(int j1 = 0; j1<no_states; j1++){
                int state = j1;
                Double val =counts[state];
                if(val==0) continue;
               
                addCount(st, state, val, dist);
            }
            
        }
    }


     void addCount(int st, int state, double val, double dist) {
		this.addCount(st, state,val);
		
	}

	public abstract AbstractTransitionProbs clone(boolean swtch);
    public abstract void print(PrintWriter pw, Double[] hittingProb, double dist);

    public abstract void addCount(int indexFrom, int indexTo, double d);
    public abstract int noStates();
    public abstract AbstractTransitionProbs clone( int[] statesToGroup, double[] u) ;
    public  void addCount(int from, int groupFrom, int groupTo, double val, double dist){
        throw new RuntimeException("!!");
    }
    public String info() {
        return this.getClass().toString();
    }
	/** pos is the index for taking the information on pseudo counts from 
	 * d is multiplier of pseudo count for each different 'category' 
	 * */
    public abstract double transfer( double[] pseudoCExp, double[][] d, int pos_index);
	public abstract double transferAlpha(double[] pseudoTrans, double[][] alpha_overall, int pos_index) ;
	//	// TODO Auto-generated method stub
		//if(true) throw new RuntimeException("!!");
		//return 0;
	//}
	public abstract double[] getAlphaPrior();
	public double transferAlpha(double[] pseudoTrans, double[][] ds, int pos_index, int[][] groupToState) {
		return this.transferAlpha(pseudoTrans, ds, pos_index);
	}
	public void setHP(double[] probs) {
		// TODO Auto-generated method stub
		
	}
	public double[] countsFrom(int i){
    	throw new RuntimeException("!!");
    }

	public Object mat() {
		// TODO Auto-generated method stub
		return null;
	}

	public double getRate(int i) {
		// TODO Auto-generated method stub
		return 0;
	}

	

	

	
	
	

}