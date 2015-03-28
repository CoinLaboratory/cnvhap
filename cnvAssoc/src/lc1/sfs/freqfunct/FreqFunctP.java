/**
 * 
 */
package lc1.sfs.freqfunct;

public class FreqFunctP extends Funct{
	final double lower, upper;
	public static Short numArgs = 1;
	public FreqFunctP(Double thresh, Double thresh1){
		this.lower = thresh;
		this.upper = (thresh1);

	}
	
	public double calc(double[] v, int[] pops) {
		return v[0]>=lower && v[0] <=upper ? 1.0 : 0;
	}
    public	String getName() {
		return lower+"<=frequency <="+upper;
	}
	
	
}