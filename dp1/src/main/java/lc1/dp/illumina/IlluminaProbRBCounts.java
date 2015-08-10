package lc1.dp.illumina;

import java.io.File;
import java.io.PrintWriter;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.OrthogonalProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.stats.TrainableBinomialDistr;
import lc1.stats.UniformDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class IlluminaProbRBCounts extends IlluminaProbRB{

	public double[][] errors;

	
	public void addCollection(IlluminaProbR probRB) {
		super.addCollection(probRB);
		errors = (double[][]) Constants.join(errors, ((IlluminaProbRBCounts)probRB).errors).toArray(new double[0][]);
	}
	
	
	public void reverse(){
		super.reverse();
		Constants.reverse(errors);
		
	}
	
	public IlluminaProbRBCounts(EmissionStateSpace emstsp, boolean train,
			int ik, int noSNPS, File clusterFile) {
		super(emstsp, train, ik, noSNPS, clusterFile);
		errors = new double[noSNPS][2];
		for(int k=0; k<errors.length; k++){
			errors[k][0] =A2B;
			errors[k][1] = B2A;
		}
		
		// TODO Auto-generated constructor stub
	}
	@Override
	 public void print(PrintWriter pw) {
		for(int i=0; i<this.errors.length; i++){
 		   pw.println(errors[i][0]+" "+errors[i][1]);
 	   }
	 }
static double A2B = Constants.A2B();
static double B2A = Constants.A2B();

	@Override
	protected ProbabilityDistribution2 getInnerDistribution(int ind
		 	  ) {
		int noB = this.emstsp.getBCount(ind);
		int noC = this.emstsp.getCN(ind);
		int noA = noC -noB;
		double p =  noC ==0 ?  Math.random() : (noA*A2B + noB*(1-B2A))/noC;
//			noB==0 ?  A2B : (noB==2 ? 1-B2A : 0.5+(A2B-B2A)/2.0);
		 	int ind1 = this.all_alias[ind];
		 	DoubleMatrix1D v_mean = (new DenseDoubleMatrix1D(this.basis_mean.getVals(ind)));
		 	DoubleMatrix1D v_var = (new DenseDoubleMatrix1D(this.basis_var.getVals(ind)));
		 	Comparable compa = this.emstsp.get(ind);
			double mean_i_r = v_mean.zDotProduct(this.prior.priorR[ind1]);
			double var_i_r = Math.sqrt(v_var.zDotProduct(this.prior.priorRVar[ind1]));
		 	double mean_i_b = p;//v_mean.zDotProduct(((PriorRB)prior).priorBaf[ind1]);
		 	double var_i_b =Math.sqrt(v_var.zDotProduct(((PriorRB)prior).priorBafVar[ind1]));
		 	//if(Double.isNaN(var_i_b)){
		 	//	throw new RuntimeException("!!");
		 	//}
		 	//if(mean_i_b<0){
		 	//	return 
		 	//}
		 	
		 	return
		 	
		 	new TrainableBinomialDistr(compa+":"+mean_i_r, mean_i_r, var_i_r, var_i_r, mean_i_b, var_i_b,var_i_b);
		 		
		 //	return  // Constants.orthogonal() ?
		  //     		new OrthogonalProbabilityDistribution(
		   //    				new UniformDistribution(Constants.minR(index), Constants.maxR(index)),
		    //   				disty
		       			//.clone()
		      // 				);
		   //    			   :
		    //   	   new CutNormal21(
		      // 		new TrainableNormal2(compa+":", mean_i_r, var_i_r,var_i_r,
		       //			   mean_i_b, var_i_b, var_i_b),0.0,1.0,true);
		       			
		 }
	
}
