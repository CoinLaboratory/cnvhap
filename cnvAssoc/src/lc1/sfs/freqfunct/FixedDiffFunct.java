/**
 * 
 */
package lc1.sfs.freqfunct;

public class FixedDiffFunct extends Funct{
	final double fixedThresh;
	final double fixedThreshUpper;
	public static Short numArgs = 2;
	public FixedDiffFunct(Double fixed0,Double fixed1){
		this.fixedThresh = fixed0;
		this.fixedThreshUpper = fixed1;
	}
    public 	String getName() {
		return fixedThresh+ "<= fixed_diff <= "+fixedThreshUpper;
	}
	public double calc(double[]v, int[] pops) {
		double res = Math.abs(v[0] -v[1]);
		if(res>=fixedThresh && res < fixedThreshUpper) return 1.0;
		else return 0.0;
	}
}