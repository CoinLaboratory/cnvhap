/**
 * 
 */
package data.pca;

import java.util.Arrays;

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public class GramSchmidt{
	 RealVector[] pcs;
	    RealVector[] proj;
	  //  int[] naInds;
	    int noRow;
	    boolean[] na;
		public GramSchmidt(RealMatrix realMatrix, int[] rowInds, int[] naInds) {
			int[] colInds = new int[realMatrix.getColumnDimension()];
			for(int k=0; k<colInds.length; k++) colInds[k] = k;
			RealMatrix rm = realMatrix.getSubMatrix(rowInds, colInds);
			this.pcs = new RealVector[realMatrix.getColumnDimension()];
			this.proj = new RealVector[realMatrix.getColumnDimension()];
		//	this.naInds = naInds;
			this.noRow = realMatrix.getRowDimension();
			
			na = new boolean[noRow];
			Arrays.fill(na, false);
			for(int k=0; k<pcs.length; k++){
				RealVector t = realMatrix.getColumnVector(k);
				t =  t.mapDivide(rm.getColumnVector(k).getNorm());
				pcs[k] =t;
				for(int j=0; j<naInds.length; j++){
					pcs[k].setEntry(naInds[j], 0);
				}
			}
		}
		public GramSchmidt() {
			  pcs= new RealVector[0];
			 proj = new RealVector[0];
		}
		public RealVector removeProj(RealVector d) {
			System.err.println("");
			for(int j=0; j<noRow; j++){
				if(Double.isNaN(d.getEntry(j))){
					d.setEntry(j, 0);
					na[j] = true;
				}
			}
			double n = d.getNorm();
			System.err.println(n);
			for(int k=0; k<this.pcs.length; k++){
				this.proj[k] = d.projection(pcs[k]);
			}
			for(int k=0; k<this.pcs.length; k++){
				d = d.subtract(proj[k]);
			}
			for(int j=0; j<noRow; j++){
				if(na[j]) {d.setEntry(j, Double.NaN);
					na[j] = false;
				}
			}
			return d;
		}
		public RealVector removeProjLast(RealVector d) {
			int k = pcs.length-1;
			if(k>=0){
				this.proj[k] = d.projection(pcs[k]);
				d = d.subtract(proj[k]);
			}
			return d;
		}
		
		public boolean checkOrth(RealVector p){
			for(int k=0; k<pcs.length; k++){
				double dp = Math.abs(p.dotProduct(pcs[k]));
				System.err.println("dot product "+dp);
				if(dp>1e-3){
					System.err.println("not orthogonal "+k);
					return false;
				}
			}
			return true;
		}
		public void update(RealVector p) {
		
			RealVector[] pcn = new RealVector[pcs.length+1];
			System.arraycopy(pcs, 0, pcn, 0, pcs.length);
			pcn[pcs.length] = p;
		    this.pcs = pcn;
		    this.proj = new RealVector[pcs.length];
		  
		}
}