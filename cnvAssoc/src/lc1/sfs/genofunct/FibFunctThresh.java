package lc1.sfs.genofunct;

public class FibFunctThresh extends FibFunct {
final double min,max;
public static Short numArgs = 1;
	//pops is number of individuals, ploidy is 2, v is allele frequency
	@Override
	public double calc(double[] v, double[] vhet, int[] pops) {
		double r = super.calc(v, vhet, pops);
		return r>min && r<=max ? 1 : 0;
	}

	public FibFunctThresh(Double l,Double u){
		this.max = u;this.min = l;
	}
	public FibFunctThresh(Double l){
		this.max = Double.POSITIVE_INFINITY;this.min = l;
	}
	
	
	@Override
	public String getName() {
		return min+"<inbreeding_coeff<="+max;
	} 

}
