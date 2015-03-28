package lc1.dp.transition;

import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

public class FreeTransitionsProbs2 extends FreeTransitionProbs1 {
public AbstractTransitionProbs prob;
	
	public FreeTransitionsProbs2(AbstractTransitionProbs prob){
		super(prob.noStates());
		this.initialise();
		this.prob = prob;
	}
	@Override
	public void addCount(int from, int to, double val){
		prob.addCount(from, to, val);
	}
	@Override
	public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
		return prob.transfer(pseudoExp, d, pos_index);
	}
	 public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
		   double res =  prob.transferAlpha(pseudoTrans, alpha_overall, pos_index);
		   if(pos_index==2) {
			   initialise();
		   }
		   return res;
	 }
	 
	 public void initialise(){
		 for(int j=0; j<no_states; j++){
	            double[] probs1 = new double[no_states];
	            for(int j1=0; j1<no_states; j1++){
	                probs1[j1] = getTransition(j, j1);
	            }
	         double v =   Constants.sum(probs1);
	         if(Math.abs(v-1.0)>1e-12){
	        	 try{
	        	throw new RuntimeException("problem ");
	        	 }catch(Exception exc){
	        	//	 exc.printStackTrace();
	        	 }
	        	 Constants.normalise(probs1);
	         }
	            transitionsOut[j] =  new SimpleExtendedDistribution1(probs1, Double.POSITIVE_INFINITY);
	        }
	 }
	 
	
}
