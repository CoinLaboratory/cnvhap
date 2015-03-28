/**
 * 
 */
package lc1.sfs.freqfunct;


public class FstFunctThresh extends Fstfunct{
	final double thresh;
	final double upper_thresh;
	public FstFunctThresh(Double thresh0, Double thresh1){
		this.thresh = thresh0;
		this.upper_thresh = thresh1;
	}
	public double calc(double[] v, int[] pop1) {
		double v1 = super.calc(v, pop1);
		return v1<=upper_thresh && v1>=thresh ? 1:0;
	}
	public String getName() {
		return thresh+"<=fst <= "+upper_thresh;
	}
}