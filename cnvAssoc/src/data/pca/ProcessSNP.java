package data.pca;

import java.io.File;

import lc1.util.CompressDir;
import lc1.util.Constants;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

abstract class ProcessSNP {
	
	GramSchmidt pcsToAdjust;
	public ProcessSNP(File[] dir, String nme, boolean str, GramSchmidt pcsToAdjust){
		this.outfile = new File[dir.length];
		this.pcsToAdjust = pcsToAdjust;
		for(int k=0; k<outfile.length; k++){
			
		outfile[k] = new File(dir[k], nme);
		}
		this.str = str;
	}
	final boolean str;
	int cnt=0;
	final File[] outfile;
	int[] sze;
	RealVector d = null;
	/* (non-Javadoc)
	 * @see data.ProcessSNP#setSize(int)
	 */
	public void setSize(int[] sze) {
		this.sze = sze;
		d = new ArrayRealVector(Constants.sum(sze));
		
	}


	public abstract RealVector process(double[] res);

	public abstract double  getHistogram() throws Exception;


	public void print(int k, int[][] aliasRevs, CompressDir[] comp, String nme) {
		// TODO Auto-generated method stub
		
	}


	public double process(RealVector realVector) {
		return 0;
		
	}


	public void refresh() {
		// TODO Auto-generated method stub
		
	}


	public int count() {
		// TODO Auto-generated method stub
		return this.cnt;
	}

}