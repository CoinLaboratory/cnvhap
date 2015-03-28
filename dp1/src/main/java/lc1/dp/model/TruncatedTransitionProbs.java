package lc1.dp.model;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.dp.transition.AbstractTransitionProbs;
import lc1.dp.transition.MatrixExp;

public class TruncatedTransitionProbs extends AbstractTransitionProbs {
final AbstractTransitionProbs prob;
final double[] totalSum;
final int num;

public TruncatedTransitionProbs(AbstractTransitionProbs p, int len){
	this.prob = p;
	this.totalSum = new double[len];
	this.num = len;
	this.reinitialise();
	//System.err.println("done");
	
}
public TruncatedTransitionProbs(
		TruncatedTransitionProbs prob, boolean swtch) {
	this.prob = prob.clone(swtch);
	this.totalSum = new double[prob.num];
	this.reinitialise();
	this.num = prob.num;
}

public TruncatedTransitionProbs(
		TruncatedTransitionProbs prob) throws Exception{
	this.prob = (AbstractTransitionProbs)prob.clone();
	this.totalSum = new double[prob.num];
	this.reinitialise();
	this.num = prob.num;
}
public TruncatedTransitionProbs(TruncatedTransitionProbs prob, int[] statesToGroup, double[] u) {
	this.prob = prob.clone(statesToGroup, u);
	this.totalSum = new double[prob.num];
	this.reinitialise();
	this.num = prob.num;
}
public void reinitialise(){
	Arrays.fill(totalSum,0.0);
	for(int i=0; i<num; i++){
		for(int j=0; j<num; j++){
			totalSum[i]+=prob.getTransition(i, j);
		}
	}
}
	@Override
	public void addCount(int indexFrom, int indexTo, double d) {
		prob.addCount(indexFrom, indexTo, d);

	}

	@Override
	public AbstractTransitionProbs clone(boolean swtch) {
	return new TruncatedTransitionProbs(this, swtch);
	}

	@Override
	public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
		return new TruncatedTransitionProbs(this, statesToGroup, u);
	}

	@Override
	public double[] getAlphaPrior() {
		return prob.getAlphaPrior();
	}

	@Override
	public Collection getDistributions() {
		throw new RuntimeException("!!");
	}

	@Override
	public double getTransition(int from, int to) {
		if(totalSum[from]==0) return 0;
		return prob.getTransition(from, to)/this.totalSum[from];
	}

	@Override
	public void initialiseCounts(boolean start, boolean end) {
		//this.reinitialise();
		prob.initialiseCounts(start, end);

	}

	@Override
	public int noStates() {
	return this.num;
	}

	@Override
	public void print(PrintWriter pw, Double[] hittingProb, double dist) {
	prob.print(pw, hittingProb, dist);

	}

	@Override
	public double transfer(double[] pseudoCExp, double[][] d, int pos_index) {
	return prob.transfer(pseudoCExp, d, pos_index);
	}

	@Override
	public double transferAlpha(double[] pseudoTrans, double[][] alpha_overall,
			int pos_index) {
		return prob.transferAlpha(pseudoTrans, alpha_overall, pos_index);
	}
	
	 public double  transferQ(double[] ds, double pseudoAlpha, double pseudoRate, MatrixExp initial, int pos_index, double distance, int i) {
			
		 throw new RuntimeException("!!"); 
	 }

	@Override
	public void validate() {
		prob.validate();
		   super.validate(this.totalSum[1]==0);
	}

}
