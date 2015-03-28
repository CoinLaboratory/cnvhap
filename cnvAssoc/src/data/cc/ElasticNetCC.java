package data.cc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public class ElasticNetCC extends CC{
	static double thresh_ = 0.00;

	public ElasticNetCC(AbstractData[] data, File dir, Double[] q,  int[] band) throws Exception {
		super(data, dir,  q,  band);
		super.startEng();
		re.eval("library(glmnet)");
	
	}
	
	@Override
	public RealMatrix solve(boolean binom, boolean onesProbe, RealMatrix realMatrix,
			RealMatrix b,  int k) {
		RealVector yj = b.getColumnVector(0);
		String family =binom ? "binomial" : "gaussian";
		if (re == null)
			startEng();
		StringBuffer assocline = new StringBuffer("y~1");
		List<String> names_ = new ArrayList<String>();
		StringBuffer matr = new StringBuffer("x = cbind(");
		re.assign("y", yj.getData());// new double[] {0,1,2,3});//pheno);
		for (int j = 0; j < realMatrix.getColumnDimension(); j++) {
			if (j > 0 || !onesProbe) {
				re.assign("x_" + j, realMatrix.getColumnVector(j).getData());// new
																				// double[]
																				// {0,0,1,1});//geno);
				assocline.append("+x_" + j);
				names_.add("x_" + j);
				matr.append("x_"+j);
				matr.append(j< realMatrix.getColumnDimension()-1 ? "," : ")");
			} else if (j == 0 && onesProbe) {
				names_.add("(Intercept)");
			}
		}
		RealMatrix beta = new Array2DRowRealMatrix(names_.size(), 1);
		StringBuffer form =new StringBuffer( "y ~ 1");
	
		    	re.eval(matr.toString());
		    	re.eval("fit = glmnet(x,y,nlambda=500,family=\""+family+"\")");
		    	//re.eval("print(dim(fit$beta))");
		    	double[] be = this.getRow(re.eval("fit$beta"),re.eval("dim(fit$beta)").asIntArray(),q[k], thresh_);
		    	if(be.length>=names_.size()-1){
		    	for (int j = 0; j < names_.size(); j++) {
					beta.setEntry(j, 0, j==0  || Math.abs(be[j-1]) <thresh_? 0 : be[j-1]);
		    		if(j>0 && beta.getEntry(j, 0)>=thresh_){
						form.append("+"+names_.get(j));
					}
				}
		    	}
		return beta;
	}
	

}
