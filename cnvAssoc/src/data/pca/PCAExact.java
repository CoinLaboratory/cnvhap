package data.pca;

import java.io.File;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import lc1.util.CompressDir;

import org.apache.commons.math.linear.RealVector;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class PCAExact extends ProcessSNP {

	public void refresh(){
		this.cols.clear();
	}
	
	List<RealVector> cols = new ArrayList<RealVector>();
	public void add(RealVector rv){
		cols.add(rv);
	}
	Jama.Matrix getRM()
	{
	//	RealMatrix rm =   new Array2DRowRealMatrix(cols.get(0).getDimension(),cols.size());
		Jama.Matrix rm = new Jama.Matrix(cols.get(0).getDimension(),cols.size());
		int i=0;
		for(Iterator<RealVector> it = cols.iterator(); it.hasNext();i++){
			RealVector nxt = it.next();
			for(int j=0; j<nxt.getDimension(); j++){
				rm.set(j, i, nxt.getEntry(j));
			}
			it.remove();
		}
		return rm;
	}


	Jama.Matrix U;
	
   
	
	//List<Double> norm= new ArrayList<Double>();
	
	public PCAExact(File[] dir, String nme, boolean str, GramSchmidt gs) {
		super(dir, nme, str,gs);
		
		// TODO Auto-generated constructor stub
	}

	
	@Override
	public void print(int j, int[][] aliasRevs, CompressDir[] comp, String nme) {
		try{
		int start =0;
		 int n = !AbstractProcessZip.singleEntry ? 1: Math.min(ProcessZip.no_reps, U.getRowDimension());
		for(int kk=0; kk<this.sze.length; kk++){
			OutputStreamWriter pw = comp[kk].getWriter(nme, true);
			int[] alias = aliasRevs[kk];
			int len = alias.length;
			  for(int k=0; k<len; k++){
				  int k1 = alias[k];
			//	for(int j=0; j<gs.pcs.length; j++){
				  if(j<0){
					 
						 for(int j1=0; j1<n; j1++){
							 pw.write(k1<0 ? "NaN":String.format("%5.3g", U.get(k1+start,j1)*1000.0).trim());
							pw.write(j1==n-1 ? "\n":"\t");
						 }
					  }else{
				  pw.write(k1<0 ? "NaN":String.format("%5.3g", U.get(k1+start,j)*1000.0).trim());
				  pw.write(j<n-1 ? "\t":"\n");
					  }
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
		Matrix matrix = getRM();
		if(false){
	EigenvalueDecomposition ev = new EigenvalueDecomposition(	matrix.transpose().times(matrix));
	 U = ev.getV();
		}
		else{	
			System.err.println("calculating SVD "+matrix.getColumnDimension()+" "+matrix.getRowDimension());
			Jama.SingularValueDecomposition svd =     new Jama.SingularValueDecomposition(matrix);
//				  new Jama.Matrix(data.transpose().getData()));//alg.transpose(X));
			// double[] d =  svd.getSingularValues();
			// Jama.Matrix U1 = svd.getU();
		//Matrix	
		U = svd.getU();
	}
		
		//	SingularValueDecompositionImpl svd = new SingularValueDecompositionImpl();
		//	U = svd.getU();
			
		return 0;
	}
	
	

	@Override
	public RealVector process(double[] res) {
		double mean =0;
		double cnt=0;
		
		double v;
		for(int k=0; k<res.length; k++){
			v = res[k];
			if(!Double.isNaN(v)){
				mean+=v;
				cnt++;
			}
		}
		mean = mean/cnt;
		double var =0;
		for(int k=0; k<res.length; k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				var+=Math.pow(v-mean,2);
			}
			
		}
		if(var>0){
			var = Math.sqrt(var/cnt);
		for(int k=0; k<d.getDimension(); k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				d.setEntry(k,  (v-mean)/var);
			}else{
				d.setEntry(k, 0);
			}
		}
		if(this.pcsToAdjust!=null){
			d = pcsToAdjust.removeProj(d);
		}
		this.cnt++;
		this.cols.add(d.copy());
		}else{
			System.err.println("var is zero");
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
