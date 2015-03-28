package data.pca;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import lc1.util.CompressDir;
import lc1.util.Constants;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class PCAExact1 extends ProcessSNP {

	double[][] covariance;
	public void setSize(int[] sze) {
		
		super.setSize(sze);
		int totsze = Constants.sum(sze);
		this.covariance = new double[totsze][totsze];
		
	}
	


	Jama.Matrix U;
	
   
	
	//List<Double> norm= new ArrayList<Double>();
	
	public PCAExact1(File[] dir, String nme, boolean str, GramSchmidt gs) {
		super(dir, nme, str,gs);
		
		// TODO Auto-generated constructor stub
	}

	
	@Override
	public void print(int j, int[][] aliasRevs, CompressDir[] comp, String nme) {
		try{
	
		int start =0;
		for(int kk=0; kk<this.sze.length; kk++){
			OutputStreamWriter pw = comp[kk].getWriter(nme, true);
			int[] alias = aliasRevs[kk];
			int len = alias.length;
			  for(int k=0; k<len; k++){
				  int k1 = alias[k];
			//	for(int j=0; j<gs.pcs.length; j++){
					pw.write(k1<0 ? "NaN":String.format("%5.3g", U.get(k1+start,j)*1000.0).trim());
					pw.write(k<len-1 ? "\t":"\n");
				//}
			}
			comp[kk].closeWriter(pw);
			start+=sze[kk];
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	@Override
	public double getHistogram() throws Exception {
		if(true){
	EigenvalueDecomposition ev = new EigenvalueDecomposition(new Matrix(this.covariance));
	 U = ev.getV();
		}
/*		else{	Jama.SingularValueDecomposition svd =     new Jama.SingularValueDecomposition(matrix);
//				  new Jama.Matrix(data.transpose().getData()));//alg.transpose(X));
			// double[] d =  svd.getSingularValues();
			// Jama.Matrix U1 = svd.getU();
		//Matrix	
		U = svd.getV();//.transpose();
	}*/
		
		//	SingularValueDecompositionImpl svd = new SingularValueDecompositionImpl();
		//	U = svd.getU();
			
		return 0;
	}
	
	

	@Override
	public RealVector process(double[] res) {
		double mean =0;
		double cnt=0;
		this.cnt++;
		double v;
		int dim = res.length;
		for(int k=0; k<dim; k++){
			v = res[k];
			
//			d.setEntry(k, v);
			if(!Double.isNaN(v)){
				mean+=v;
				cnt++;
			}
		}
		mean = mean/cnt;
		double var =0;
		for(int k=0; k<dim; k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				var+=Math.pow(v-mean,2);
			}
			
		}
		if(var>0){
			var = Math.sqrt(var/cnt);
		for(int k=0; k<dim; k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				res[k] =   (v-mean)/var;
			}else{
			res[k] = 0;
			}
		}
		
		for(int k=0; k<dim; k++){
			for(int k1 = 0; k1< dim; k1++){
				covariance[k][k1]+=res[k]*res[k1];
			}
		}
		}
		return null;
		
	}
	private double count(char[] st, char c) {
		int b =0;
		for(int j=0; j<st.length; j++){
			if(st[j]==c) b++;
		}
		return b;
	}

	
	

	

}
