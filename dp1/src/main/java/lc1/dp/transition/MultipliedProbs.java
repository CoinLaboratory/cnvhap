package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Collection;

public class MultipliedProbs extends AbstractTransitionProbs{

	AbstractTransitionProbs[] orig;
	AbstractTransitionProbs within;
	final int numwithin;
   public MultipliedProbs(MultipliedProbs multipliedProbs, boolean swtch) {
		this.numwithin = multipliedProbs.numwithin;
		this.orig = new AbstractTransitionProbs[numwithin];
		for(int k=0; k<numwithin; k++){
			this.orig[k] = multipliedProbs.orig[k].clone(swtch);
		}
		this.within = multipliedProbs.within.clone(swtch);
	
	}

//	int[][] groupToState;
	
	public MultipliedProbs(AbstractTransitionProbs orig2, AbstractTransitionProbs siteTransitions) {
		this.numwithin = siteTransitions.noStates()-1;
	this.orig = new AbstractTransitionProbs[numwithin];
	for(int i=0; i<numwithin; i++){
		orig[i] = new FreeTransitionProbs1(orig2);
	}
	this.within = siteTransitions;
	
}

	@Override
	public void validate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getTransition(int from, int to) {
		
		 int groupFrom = (int) Math.floor((double)(from-1)/(double)numwithin);
		 int groupTo = (int) Math.floor((double)(to-1)/(double) numwithin);
	     int indexFrom = from - groupFrom*numwithin;
	     int indexTo = to - groupTo*numwithin;
		// TODO Auto-generated method stub
		return orig[indexFrom-1].getTransition(groupFrom+1, groupTo+1)* within.getTransition(indexFrom, indexTo);
	}

	@Override
	public void initialiseCounts(boolean start, boolean end) {
		for(int k=0; k<orig.length; k++){
			orig[k].initialiseCounts(start, end);
		}
		within.initialiseCounts(start, end);
		
	}

	@Override
	public Collection getDistributions() {
		// TODO Auto-generated method stub
		  throw new RuntimeException("!!");
	}

	@Override
	public AbstractTransitionProbs clone(boolean swtch) {
		return new MultipliedProbs(this,swtch);
	}

	@Override
	public void print(PrintWriter pw, Double[] hittingProb, double dist) {
		for(int k=0; k<orig.length; k++){
		this.orig[k].print(pw, hittingProb, dist);
		
		}
		this.within.print(pw, hittingProb, dist);
		
	}

	@Override
	public void addCount(int from, int to, double d) {
		int groupFrom = (int) Math.floor((double)(from-1)/(double)numwithin);
		 int groupTo = (int) Math.floor((double)(to-1)/(double) numwithin);
	     int indexFrom = from - groupFrom*numwithin;
	     int indexTo = to - groupTo*numwithin;
		// TODO Auto-generated method stub
		this.orig[indexFrom-1].addCount(groupFrom+1, groupTo+1,d);
		this.within.addCount(indexFrom, indexTo, d);
	}

	@Override
	public int noStates() {
		// TODO Auto-generated method stub
		return this.numwithin * this.orig[0].noStates();
	}

	@Override
	public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
		// TODO Auto-generated method stub
	    throw new RuntimeException("!!");
	}

	@Override
	public double transfer(double[] pseudoExp, double[][] d, int pos_index) {
		if(pos_index==2){
			double logL = this.orig[0].transfer( pseudoExp,d,0);
			for(int k=1; k<orig.length; k++){
				logL+=orig[k].transfer( pseudoExp,d,0);
			}
		        	//should be pos_index+1
		                logL+= this.within.transfer( pseudoExp, null, 1);
		    	return logL;
			}
			else{
				 double logL = this.orig[0].transfer( pseudoExp,d,0);
				 for(int k=1; k<orig.length; k++){
						logL+=orig[k].transfer( pseudoExp,d,0);
					}
				 return logL;
			}
	}

	@Override
	public double transferAlpha(double[] pseudoTrans, double[][] alpha_overall,
			int pos_index) {
		if(pos_index==2){
    		double logL = this.orig[0].transferAlpha( pseudoTrans, alpha_overall, 0);
    		 for(int k=1; k<orig.length; k++){
    			 logL += this.orig[k].transferAlpha( pseudoTrans, alpha_overall, 0);
    		 }
    			 ///note should be pos_index + 1, but this would lead to 3rd index
                 logL+=this.within.transferAlpha(pseudoTrans, null, 1 );//, alpha_overall, pos_index+1);
                 return logL;
    	}
    	else{
    		double logL = this.orig[0].transferAlpha( pseudoTrans, alpha_overall, pos_index);
    		 for(int k=1; k<orig.length; k++){
    			 logL += this.orig[k].transferAlpha( pseudoTrans, alpha_overall, 0);
    		 }
    		return logL;
    	}
	}

	@Override
	public double[] getAlphaPrior() {
		// TODO Auto-generated method stub
		return this.orig[0].getAlphaPrior();
	}

}
