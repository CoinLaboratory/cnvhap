/**
 * 
 */
package lc1.sfs.freqfunct;

public class AlleleDifffunct extends Funct{
public static Short numArgs = 2;
	
	public String toString(){
		return this.getName();
	}
	public double calc(double[] v, int[] pops) {
		return Math.abs(v[1] - v[0]);
	}
	public String getName() {
		return "AlleleDiff";
	}
	
}