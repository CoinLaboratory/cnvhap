package lc1.stats;

import pal.statistics.ChiSquareDistribution;

public class ArmitageTrendTest {
	ChiSq chisq = new ChiSq();
	public static void main(String[] args){
		double[] cse = new double[] {89, 369, 342}; //[AA, Aa, aa]
		double[] cntrl = new double[] {56, 250, 266};//[AA, Aa, aa]
		ArmitageTrendTest att = new ArmitageTrendTest();
		att.set(cse, cntrl);
		System.err.println(att.chisq()+" "+att.getSig());
		double[] cse1 = new double[] {0.000,60.518,576.480,0.002,0.000};
		double[] cntrl1 = new double[] {0.000,173.169,321.831,1.000,0.000};//};//[AA, Aa, aa]
		att.set(cse1, cntrl1);
		System.err.println(att.chisq()+" "+att.getSig());
		System.err.println(ChiSquareDistribution.quantile(0.95, 1.0));
		System.err.println(  ChiSquareDistribution.cdf(3.84, 1));
		double[] d = new double[] {1e-3, 1e-5, 1e-6};
		for(int i=0; i<d.length; i++){
			double y = ChiSquareDistribution.quantile(1.0-d[i], 1);
			System.err.println(d[i]+" "+y);
		}
	}
	double[] cse;
	double[] cntrl;
	double[] total;
	double sum_nx=0;
	double sum_Nx=0;
	double sum_Nx2=0;
	double T=0;
	double t=0;
	/* cse is count of cases : [AA, Aa, aa]  similarly for cntrl */
	public void set(double[] cse, double[]cntrl){
		this.cse = cse;
		this.cntrl = cntrl;
		this.total = new double[cse.length];
		sum_nx =0;
		sum_Nx =0;
		sum_Nx2 = 0;
		T=0; 
		t=0;
		for(int i=0; i<total.length; i++){
			total[i]= cse[i]+cntrl[i];
			T+=total[i];
			t+=cse[i];
			sum_nx+=i*cse[i];
			sum_Nx+=i*total[i];
			sum_Nx2+=Math.pow(i,2)*total[i];
		}
		
	}
	public double getSig(){
		double res =  chisq.chi2prob(1, chisq());
		return res;
	}
	public double chisq(){
		double num = T *Math.pow(T*sum_nx - t*sum_Nx,2);
		double denom = t*(T-t)*(T*sum_Nx2 - Math.pow(sum_Nx, 2));
		if(num==0 || (denom==0 && num<1e-5) ) return 0;
	//	double p =  ChiSq.chi2prob1(1, num/denom);
		return num/denom;
	}
	
	
}
