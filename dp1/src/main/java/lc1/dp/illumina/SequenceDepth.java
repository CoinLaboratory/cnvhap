package lc1.dp.illumina;

import java.io.File;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.OrthogonalProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.stats.TrainableNegBinomial;
import lc1.stats.UniformDistribution;
import lc1.util.Constants;

public class SequenceDepth extends IlluminaProbR {
    
   /* public static void main(String[] args){
        CompoundEmissionStateSpace stSp = Emiss.getEmissionStateSpace(1, 2);
        IlluminaProbB probB = new IlluminaProbB(stSp, true, 1.0,0);
        List<Double> l = new ArrayList<Double>();
      //  probB.ad
      ///  probB.s
    }*/
    
	
	
    public SequenceDepth(EmissionStateSpace stSp,
			 boolean b,
			int index, int noSnps,  File clusterFile) {
		super(stSp, b,index, noSnps, clusterFile);
	}
    
    public void maximiseGlobal(double[] pseudo,  int it) {
    	
    }
   /* public Prior calculatePriorsFromLocs(boolean initialise, double[] pseudo , double pseudoFill, int it){
    	return null;
    }*/
    
   
     @Override
      protected ProbabilityDistribution2 getInnerDistribution(int index){
    		// TODO Auto-generated method stub
    	 int ind1 = this.all_alias[index];
    	 DoubleMatrix1D v_mean = (new DenseDoubleMatrix1D(this.basis_mean.getVals(index)));
    		DoubleMatrix1D v_var = (new DenseDoubleMatrix1D(this.basis_var.getVals(index)));
    		Comparable compa = this.emstsp.get(index);
    		double mean_i_r = v_mean.zDotProduct(this.prior.priorR[ind1]);
    		double var_i_r = v_var.zDotProduct(this.prior.priorRVar[ind1]);
    		return   new OrthogonalProbabilityDistribution(
    			      new TrainableNegBinomial(compa+"", mean_i_r, var_i_r, round,1.0) ,
    	      		   new UniformDistribution(Constants.minB(index),Constants.maxB(index)));
    	}
     
		
}
