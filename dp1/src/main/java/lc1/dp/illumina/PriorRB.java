package lc1.dp.illumina;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class PriorRB extends Prior {

	public DoubleMatrix1D[]priorBaf,priorBafVar, priorRho;
	public PriorRB(int ik, int[][] indices, EmissionStateSpace emstsp, int len) {
		super(ik, indices,emstsp, len);
		
		
		// TODO Auto-generated constructor stub
	}
	
	 public boolean equals(Object obj){
		 Prior p1 = (Prior) obj;
		 return super.equals(p1) && priorBaf.equals(p1.priorR) && priorBafVar.equals(p1.priorRVar);
	 }
	@Override
	public void makePriors(int length) {
		super.makePriors(length);
		this.priorBaf =new DoubleMatrix1D[length];
		this.priorBafVar =new DoubleMatrix1D[length];
	}
	@Override
	public void makePriorNormal(int len, int i,int ik){
		super.makePriorNormal(len ,i, ik);
		this.priorBaf[i] = new DenseDoubleMatrix1D(len);
		this.priorBafVar[i] = new DenseDoubleMatrix1D(len);
		copy(Constants.b_mean(ik,i),priorBaf[i]);
		copy(Constants.b_var1(ik,i), priorBafVar[i]);
	}
	
	 public PriorRB(PriorRB prior) {
		 super(prior);
			this.priorBaf = copy(prior.priorBaf);
			this.priorBafVar = copy(prior.priorBafVar);
			this.priorRho  =copy(prior.priorRho);
		}

}
