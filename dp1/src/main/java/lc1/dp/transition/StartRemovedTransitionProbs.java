package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Collection;

public class StartRemovedTransitionProbs extends AbstractTransitionProbs {
AbstractTransitionProbs inner;
	public StartRemovedTransitionProbs(
		StartRemovedTransitionProbs startRemovedTransitionProbs, boolean swtch) {
			this.inner = startRemovedTransitionProbs.inner.clone(swtch);
}

	public StartRemovedTransitionProbs(
			StartRemovedTransitionProbs startRemovedTransitionProbs,
			int[] statesToGroup, double[] u) {
		this.inner = startRemovedTransitionProbs.inner.clone(statesToGroup, u);
	}

	public StartRemovedTransitionProbs(
			AbstractTransitionProbs abstractTransitionProbs) {
		this.inner = abstractTransitionProbs;
	}

	@Override
	public void addCount(int indexFrom, int indexTo, double d) {
		inner.addCount(indexFrom+1, indexTo+1, d);

	}
	
	@Override
	public void addCount(int indexFrom, int indexTo, double d, double dist) {
		inner.addCount(indexFrom+1, indexTo+1, d, dist);

	}

	@Override
	public AbstractTransitionProbs clone(boolean swtch) {
		return new StartRemovedTransitionProbs(this, swtch);
	}

	@Override
	public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
		return new StartRemovedTransitionProbs(this,statesToGroup, u);
	}

	@Override
	public double[] getAlphaPrior() {
		return inner.getAlphaPrior();
	}

	@Override
	public Collection getDistributions() {
		return inner.getDistributions();
	}

	@Override
	public double getTransition(int from, int to) {
	return inner.getTransition(from+1, to+1);
	}

	@Override
	public void initialiseCounts(boolean start, boolean end) {
		inner.initialiseCounts(start, end);

	}

	@Override
	public int noStates() {
		return inner.noStates();
	}

	@Override
	public void print(PrintWriter pw, Double[] hittingProb, double dist) {
		inner.print(pw, hittingProb, dist);

	}

	@Override
	public double transfer(double[] pseudoCExp, double[][] d, int pos_index) {
		 double res = inner.transfer(pseudoCExp, d, pos_index);
	//	 inner.validate();
		 return res;
		 
	}
	public double  transferQ(double[] ds,double pseudoAlpha, double pseudoRate,  MatrixExp initial, int i, double distance, int it) {
		return this.inner.transferQ(ds, pseudoAlpha, pseudoRate, initial, i, distance,it);
		
	}

	@Override
	public double transferAlpha(double[] pseudoTrans, double[][] alpha_overall,
			int pos_index) {
		return inner.transferAlpha(pseudoTrans, alpha_overall, pos_index);
	}

	@Override
	public void validate() {
		if(inner instanceof FreeTransitionProbs1) ((FreeTransitionProbs1)inner).validate(1);
		else inner.validate();

	}

}
