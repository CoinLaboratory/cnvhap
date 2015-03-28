/**
 * 
 */
package lc1.sfs.freqfunct;

//NOTE - THIS IS FREQUENCY OF REFERENCE ALLELE!!!
public class FreqFunct extends Funct{
	
	
	public static Short numArgs = 1;
	
	
	public double calc(double[] v, int[] pops) {
		if(v[0]<0) {
			throw new RuntimeException("!!");
		}
		return v[0];
	}
	public String getName() {
		return "Frequency";
	}
	
}