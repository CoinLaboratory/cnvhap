/**
 * 
 */
package lc1.dp.data.collection;

import java.util.Arrays;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

//import assoc.Constants;

public class GramSchmidt{
	 RealVector[] pcs;
	    RealVector[] proj;
	    
	    boolean[] na;
	    int noRow;
	    
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
	    
		public GramSchmidt(RealMatrix realMatrix) {
			this.pcs = new RealVector[realMatrix.getColumnDimension()];
			this.proj = new RealVector[realMatrix.getColumnDimension()];
			for(int k=0; k<pcs.length; k++){
				RealVector t = realMatrix.getColumnVector(k);
				t =  t.mapDivide(t.getNorm());
				pcs[k] =t;
				
			}
		}
		public GramSchmidt() {
			  pcs= new RealVector[0];
			 proj = new RealVector[0];
		}
		/*public RealVector removeProj(RealVector d) {
			for(int k=0; k<this.pcs.length; k++){
				this.proj[k] = d.projection(pcs[k]);
			}
			for(int k=0; k<this.pcs.length; k++){
				d = d.subtract(proj[k]);
			}
			return d;
		}*/
	
		
		static boolean standardise = true;
		public RealVector removeProj(RealVector d) {
			
			double mean =0;
			double cnt=0;
			double v;
			for(int k=0; k<this.noRow; k++){
				v = d.getEntry(k);
				if(!Double.isNaN(v)){
					mean+=v;
					cnt++;
				}
			}
			mean = mean/cnt;
			double var =0;
			for(int k=0; k<noRow; k++){
				 v = d.getEntry(k);
				if(!Double.isNaN(v)){
					var+=Math.pow(v-mean,2);
				}
				
			}
			if(var==0 || ! standardise) var = 1;
			else{
				var = Math.sqrt(var/cnt);
			}
			for(int j=0; j<noRow; j++){
				v = d.getEntry(j);
				if(Double.isNaN(v)){
					d.setEntry(j, 0);
					na[j] = true;
				}else{
					d.setEntry(j,  (v-mean)/var);
				}
			}
//			double n = d.getNorm();
		//	System.err.println(n);
			for(int k=0; k<this.pcs.length; k++){
				this.proj[k] = d.projection(pcs[k]);
			}
			for(int k=0; k<this.pcs.length; k++){
				d = d.subtract(proj[k]);
			}
			for(int j=0; j<noRow; j++){
				if(na[j]) {
					d.setEntry(j, Double.NaN);
					na[j] = false;
				}else{
					v = d.getEntry(j);
					d.setEntry(j,  (v*var+mean));
				}
			}
			return d;
		}
	/*	public RealVector removeProjLast(RealVector d) {
			int k = pcs.length-1;
			if(k>=0){
				this.proj[k] = d.projection(pcs[k]);
				d = d.subtract(proj[k]);
			}
			return d;
		}*/
		
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