package lc1.sfs.genofunct;

public class FibFunct extends GFunct {

	public static Short numArgs = 1;
	
	//pops is number of individuals, ploidy is 2, v is allele frequency
	@Override
	public double calc(double[] v, double[] vhet, int[] pops) {
		// TODO Auto-generated method stub
		double H0 = 2*v[0]*(1-v[0]);
		double H =vhet[0];
		double diff = H0-H;
		double res = Math.abs(diff) < 1e-5 ? 0 : diff/H0;
		//if(res < -1000) {
		//	throw new RuntimeException("!!");
		//}
//		return v[0];
//		return v[0];
		return res;
	}

	@Override
	public String getName() {
		return "inbreeding_coeff";
	} 

}
