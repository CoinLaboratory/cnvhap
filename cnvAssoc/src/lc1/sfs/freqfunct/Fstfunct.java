/**
 * 
 */
package lc1.sfs.freqfunct;


public class Fstfunct extends Funct{
	public static Short numArgs = 2;
	public String getName() {
		return "Fst";
	}
	@Override
	public double calc(double[] v, int[] pop1) {
			double isum=0;
			double popsum=0;
			double pqsum=0;
			for(int k=0; k<2; k++){
				isum+=v[k]*(double) pop1[k];
				popsum+=pop1[k];
				double p = v[k];
				pqsum +=2*p*(1-p)*pop1[k];
			}
			double avg2pq = pqsum/(double)popsum;
			double p = isum/popsum;
			double q = 1-p;
			double res =  1 - (avg2pq==0 ? 0 : avg2pq/(2*p*q));
			if(res<0){
				System.err.println("warning should not be negative fst:"+res+"\n v=c("+v[0]+","+v[1]+")\npop1=c("+pop1[0]+","+pop1[1]+")");
			}
			if(Double.isNaN(res)){
				throw new RuntimeException("!! NA" +v[0]+v[1]);
			}
			//if(res<0){
			//	throw new RuntimeException("!!");
			//}
			return res;
	}
	
}