package data.cc;


import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

/* Experimental class for carrying out linear regression with uncertainty!!!!
 * Best not to use at this stage, still incomplete
 *  */
public class LinearReg {

	
	NormalDistribution dist = new NormalDistributionImpl(0.0,1.0);
	RealMatrix x; //genotype matrix, rows are indiv
	RealMatrix y; //phenotype vector (single column matrix)
	
	RealMatrix x1; //expanded genotype matrix, rows are assignments of indiv to different risk cats
	RealMatrix y1; //exppaned phenotype matrix
	
	int numlevels;
	double[][] probs; // probs k is assignment of probabilities of kth individual
	
	RealMatrix beta;
	final double[] levels;
	int max = 10;
	double sd = 1.0;
	final boolean onesk;
	final RealMatrix Q;
	final double[] p;
	LinearReg(RealMatrix x, RealMatrix y, RealMatrix beta_est,RealMatrix Q, int numlev,  boolean onesk){
		this.x = x;
		this.onesk = onesk;
		this.Q= Q;
		this.y = y;
		this.numlevels = numlev;
		
		beta = beta_est.copy();
		for(int k=0; k<beta_est.getRowDimension(); k++){
			beta.setEntry(k, 0, Math.random());
		}
		this.probs = new double[x.getRowDimension()][numlev];
		this.x1 =  new Array2DRowRealMatrix(x.getRowDimension()*numlev,x.getColumnDimension());
		this.y1 =  new Array2DRowRealMatrix(x.getRowDimension()*numlev,1);
		this.levels = new double[numlevels];
		this.p = new double[numlevels];
		double mid = ((double)numlevels -1.0)/2.0;
		//double max = 
			//
		for(double k=0; k<numlevels; k++){
			 p[(int) k] = (k+1)/((double)numlevels+1);
			levels[(int)k] = -Math.log(1/p[(int)k] - 1);
				//((k-mid)/mid);
		}		
	}
	public RealMatrix run(){
		System.err.println("starting");
		double logL = Double.NEGATIVE_INFINITY;
		for(int k=0; k<5; k++){
			RealVector v0 = beta.getColumnVector(0).copy();
			double sd_p = this.sd;
			double logL1 =expectation();
			 maximisation();
			 
			 double norm = v0.subtract(beta.getColumnVector(0)).getNorm();
			 if(norm < 0.005 && Math.abs(sd - sd_p)<0.01){
			 break;
			 }
			 logL = logL1;
			 System.err.println(logL+"\t"+norm+"\t"+sd);
		}
		return this.beta;
	}
	public void maximisation(){
	
		//System.err.println(beta);
		for(int k=0; k<probs.length; k++){
			for(int j=0; j<numlevels; j++){
				int row = k*numlevels+j;
				double prob = probs[k][j];
				for(int l =0; l<x.getColumnDimension(); l++){
					x1.setEntry(row, l, prob*x.getEntry(k, l));
				}
				y1.setEntry(row, 0, levels[j]*prob);
			}
		}
		beta = solve(onesk, x1,  y1, Q);
		this.sd = Math.sqrt(sum/cnt);
		
	}
	
	/*basic method for calculating the regression */
	public  RealMatrix solve(boolean onesk, RealMatrix A, RealMatrix b,
			RealMatrix Q) {
		RealMatrix AT = A.transpose();
		RealMatrix prod = AT.multiply(A).add(Q);
		RealMatrix b1 = b.copy();
		LUDecompositionImpl lu = new LUDecompositionImpl(prod);
		DecompositionSolver solver = lu.getSolver();
		RealMatrix res = solver.solve(AT.multiply(b1));
		return res;
	}

	//List<Double> obs = new ArrayList<Double>();
	
	double sum=0;
	double cnt=0;
	//public double logL =0;
	public double  expectation(){
		sum=0; cnt=0;
		double logL=0;
		RealMatrix risk = x.multiply(beta);
		for(int k=0; k<probs.length; k++){
			double[] pr = probs[k];
			double m = risk.getEntry(k, 0);
			double val = y.getEntry(k, 0);
			boolean zero = Math.abs(val)<0.0001;
			boolean one  = Math.abs(val-1)<0.001;
			if(!zero && !one) {
				throw new RuntimeException("!!");
			}
			double sum1=0;
			for(int k1=0; k1<numlevels; k1++){
				double p_ = p[k1];
				dist.setMean(m);dist.setStandardDeviation(sd);
				double pdf =  dist.density(levels[k1]);
				double pow = Math.pow(p_, val)*Math.pow(1-p_, 1-val);;
				pr[k1] = pdf* pow;
				sum1+=pr[k1];
				logL+=Math.log(pr[k1]);
			}
			for(int k1=0; k1<numlevels; k1++){
				pr[k1] = pr[k1]/sum1;
				sum+=Math.pow(m-levels[k1],2)*pr[k1];
				cnt+=pr[k1];
			}
			
		//	System.err.println(cnt);
		}
		return logL;
	}
	
	
}
