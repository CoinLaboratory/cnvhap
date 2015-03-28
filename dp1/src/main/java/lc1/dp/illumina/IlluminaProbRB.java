package lc1.dp.illumina;

import java.io.File;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.CutNormal;
import lc1.stats.OrthogonalProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.TrainableNormal;
import lc1.stats.UniformDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class IlluminaProbRB extends IlluminaProbR{

	public IlluminaProbRB(EmissionStateSpace emstsp, boolean train, int ik,
			int noSNPS, File clusterFile) {
		super(emstsp, train, ik, noSNPS,  clusterFile);
		// TODO Auto-generated constructor stub
	}
@Override
	 void makePriors(int ik, int[][] all_indices2, EmissionStateSpace emstsp, int len) {
		 this.prior = new PriorRB(ik,all_indices, emstsp,len);
	 	this.initial =new PriorRB((PriorRB)prior);
		
	}
	
@Override
	protected ProbabilityDistribution2[] makeDists() {
		  return  
		  new ProbabilityDistribution2[emstsp.genoListSize()];
	}


@Override
	public void calculatePriorsFromLocs( double[] pseudo , double pseudoFill, int it){
		
		MultipleDataRegressionRB mdr = new MultipleDataRegressionRB(new IlluminaProbR[] {this}, new int[] {-1},pseudo,
				pseudoFill, this.initial);
		mdr.makeRegress(true,true);
		mdr.rst(this.globalRst);
		mdr.bst(this.globalRst);
		prior.priorR = copy(mdr.paramR);
		prior.priorRVar = copy(mdr.paramRVar);
		((PriorRB)prior).priorBaf = copy(mdr.paramBaf);
		((PriorRB)prior).priorBafVar = copy(mdr.paramBafVar);
		((PriorRB)prior).priorRho = copy(mdr.paramCov);

	 }
	
	 public IlluminaProbRB(IlluminaProbRB probR, boolean sameMix) {
		 super(probR, sameMix);
	 }
	 

	 protected void makeBases() {
	 	this.basis_mean = new BasisFunctionRB(emstsp, Constants.backgroundCount(index), Constants.basisNme(index));
	     this.basis_var = new BasisFunctionRB(emstsp, Constants.backgroundCount(index),Constants.basisNme(index));
	 	
	 }

	 
	 
	 protected ProbabilityDistribution2 getInnerDistribution(int ind
	 	  ) {
	 	int ind1 = this.all_alias[ind];
	 	DoubleMatrix1D v_mean = (new DenseDoubleMatrix1D(this.basis_mean.getVals(ind)));
	 	DoubleMatrix1D v_var = (new DenseDoubleMatrix1D(this.basis_var.getVals(ind)));
	 	Comparable compa = this.emstsp.get(ind);
		double mean_i_r = v_mean.zDotProduct(this.prior.priorR[ind1]);
		double var_i_r = Math.sqrt(v_var.zDotProduct(this.prior.priorRVar[ind1]));
	 	double mean_i_b = v_mean.zDotProduct(((PriorRB)prior).priorBaf[ind1]);
	 	double var_i_b =Math.sqrt(v_var.zDotProduct(((PriorRB)prior).priorBafVar[ind1]));
	 	//if(Double.isNaN(var_i_b)){
	 	//	throw new RuntimeException("!!");
	 	//}
	 	lc1.stats.ProbabilityDistribution disty = 
	 	emstsp.getCN(ind)==0 ? new UniformDistribution(Constants.minB(this.index), Constants.maxB(this.index)):
	 	!Constants.beta(this.index) ? new CutNormal(new TrainableNormal(compa+":"+mean_i_b, mean_i_b, var_i_b, var_i_b, round,1.0),0,1):
	 		new lc1.stats.Beta(compa+":"+mean_i_b+"", mean_i_b, var_i_b, var_i_b, round,1.0);
	 	return  // Constants.orthogonal() ?
	       		new OrthogonalProbabilityDistribution(
	       				new TrainableNormal(compa+":"+mean_i_r, mean_i_r, var_i_r, var_i_r, round,1.0),
	       				disty
	       			//.clone()
	       				);
	   //    			   :
	    //   	   new CutNormal21(
	      // 		new TrainableNormal2(compa+":", mean_i_r, var_i_r,var_i_r,
	       //			   mean_i_b, var_i_b, var_i_b),0.0,1.0,true);
	       			
	 }
	 
	 public final ProbabilityDistribution2 getDistributionGlobal(int noCop, int noB
		) {
   int al = this.aliasNB[noCop][noB];// b_alias[j];
ProbabilityDistribution2 dist =  
	  rGlobal[al];
  return dist;

}
}
