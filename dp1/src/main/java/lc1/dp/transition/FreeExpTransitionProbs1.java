package lc1.dp.transition;

public class FreeExpTransitionProbs1 extends FreeExpTransitionProbs {

	
	 public FreeExpTransitionProbs1(double dist, double[] transform, double d) {
		super(dist, transform, d);
	}

	public double getTransition(int from, int to){
	    	return super.getTransition(from+1, to+1);
	    }
	
	 public double getTransitionToPaint(int from, int to) {
		 return super.getTransitionToPaint(from, to);
	 }
	 @Override
	    void addCount(int from, int to, double val, double dist) {
	    	super.addCount(from+1, to+1, val);
		}
}
