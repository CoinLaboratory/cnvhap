package lc1.sfs.genofunct;

import lc1.stats.ChiSq;

public class HWEFunct extends FibFunct {
public 	static Short numArgs = 1;
	@Override
	public double calc(double[] v, double[] vhet, int[] pops) {
		double res = super.calc(v, vhet, pops);
		return ChiSq.chi2prob(1, Math.pow(res,2)*pops[0]);
	}

	@Override
	public String getName() {
		return "HWE";
	}

}
