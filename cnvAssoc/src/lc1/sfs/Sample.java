/**
 * 
 */
package lc1.sfs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import lc1.util.Constants;

import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import data.pca.ProcessVCF;

class Sample{
	
	
	
	public static void main(String[] args){
		try{
			boolean[][] incl = new boolean[2][];
			RealMatrix[] rm = new RealMatrix[2];
			int[] pop = new int[] { 60,50};
			
			double[] p = new double[] {0.15,0.5};
			double[] f = new double[] {0.0,0.0};
			Sample[] samp = new Sample[2];
			int noreps = 1;
			List<RealMatrix> res = new ArrayList<RealMatrix>();
			List<String> nme = new ArrayList<String>();
			for(int rep=0; rep<noreps; rep++){
			for(int k=0; k<samp.length; k++){
				samp[k] = new Sample(p[k], pop[k],f[k]);
				incl[k] = new boolean[pop[k]];
				Arrays.fill(incl[k], true);
				rm[k] = 	new Array2DRowRealMatrix(samp[k].lh);
			//	for(int j=0; j<(int) (Math.round((double)pop1[k]/1.05)); j++) incl[k][j] = false;
				
			}
			double fib = samp[0].Fib();
		System.err.println(fib);
	//	String[] str =   "Fstfunct:AlleleDifffunct:FixedDiffFunct_0.99;0.9;0.8;0.5:FstFunctThresh_0.99;0.9;0.9;0.8;0.5 ".split(":");
	//	Funct[] functs = Funct.getFuncts(str,nme);
		
		int[] inds = new int[] {0,1};
		int[] pop1 = (int[]) SFS.subIndexI(pop, inds);
		//RealMatrix[] rm1 = (RealMatrix[]) SFS.subIndex(rm, inds);
		//boolean[][] incl1 = (boolean[][]) SFS.subIndex(incl, inds);
		String functString = "AlleleDifffunct";
		List<SFS> fst_ = new ArrayList<SFS>();
		List<int[]> inds_ = new ArrayList<int[]>();
		List<String>nme1 = new ArrayList<String>();
	ProcessVCF.getSFS(functString, nme1, pop,fst_, inds_, "");
	SFS[] fst = fst_.toArray(new SFS[0]);
		for(int k=0; k<fst.length; k++){
			
			fst[k].setData(rm, incl,inds);
			fst[k].train(10);
			fst[k].expectation();
			res.add(fst[k].calculate());
			
			
			fst[k].print("sample.txt");
		
		}
			}
			RealMatrix rm0 = res.get(0);
			/*RealMatrix avg = new Array2DRowRealMatrix(rm0.getRowDimension(), rm0.getColumnDimension());
			RealVector var = new ArrayRealVector(rm0.getRowDimension(), rm0.getColumnDimension());
			for(int k=0; k<res.size(); k++){
				avg = avg.add(res.get(k));
			}
			avg.mapDivideToSelf((double)res.size());
			for(int k=0; k<res.size(); k++){
				var = var.add((res.get(k).subtract(avg)).mapPowToSelf(2.0));
			}
			
			var.mapDivideToSelf((double)res.size()).mapPowToSelf(0.5);*/
			for(int i=0; i<rm0.getRowDimension(); i++){
				System.err.println(nme.get(i)+" :  "+rm0.getRowVector(i));
			}		
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public double p,n;
	public int[] geno;
	public double[] counts = new double[3];
	int no_trials = 50;
	public double[][]lh;
	
	public double Fib(){
		double p = (counts[1]+counts[0]*2)/(2.0*(counts[0]+counts[1]+counts[2]));
		double H0 = 2*p*(1-p);
		double H = counts[1]/(counts[0]+counts[1]+counts[2]);
		double diff = H0-H;
		System.err.println(counts[0]+" "+counts[1]+" "+counts[2]);
		System.err.println(H0+" "+H);
		return diff/H;
	}
	public static Random random = new Random(31094);// );
	Sample(double p, int n, double f) throws Exception{
		geno = new int[n];
		 lh = new double[n][3];
		 
		double[] d = new double[] { Math.pow(p, 2)*(1-f)+p*f,(1-f)*2*p*(1-p),Math.pow(1-p, 2)*(1-f)+(1-p)*f};
		BinomialDistribution dist = new BinomialDistributionImpl(no_trials,0);
		for(int k=0; k<geno.length; k++){
			geno[k] =Constants.sample(d);
			counts[geno[k]]++;
			dist.setProbabilityOfSuccess((double)geno[k]/2.0);
			double rand =random.nextDouble();
			/*int breads = 0 ;//dist.inverseCumulativeProbability(Math.random());
			while(dist.cumulativeProbability(breads)<rand){
				breads++;
			}*/
			int breads = dist.inverseCumulativeProbability(rand)+1;
			for(int j=0; j<=2; j++){
				dist.setProbabilityOfSuccess((double)j/2.0);
				lh[k][j] = dist.probability(breads);
				
			}
//			Arrays.fill(lh[k],25);
		}
	}
}