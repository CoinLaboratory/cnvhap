/**
 * 
 */
package lc1.sfs;

import java.util.Arrays;

import org.apache.commons.math.linear.RealVector;

public class ArrayIt{
	
	int[] inds;
	int[] lim;
//	double[] v;
	double[] lh; //this accumulates likelihood
	
//	double sum =0;
	
	boolean hasNext = true;
	int len;
	ArrayIt(int[] sze, int noFuncts){
		inds = new int[sze.length];
		lim = sze;
			//new int[sze.length];
		//for(int k=0; k<lim.length; k++){
		//	lim[k]= sze[k];
		//}
		len = inds.length-1;
		//this.v = new double[sze.length];
		this.lh = new double[noFuncts];
	}
	public boolean incr(){
		return this.incr(len);
	}
	
	
	
	
	
	
	
	public boolean incr(int i){
		inds[i]++;
		//v[i] = (double) inds[i] / (double)lim[i];
		boolean res = true;
		if(inds[i]==lim[i]){
			if(i==0) return false;
			else{
				res = incr(i-1);
				inds[i] =0;
				//v[i] = 0;
			}
		}
		return res;
	}
	public void reset(){
		Arrays.fill(inds, 0);
	//	Arrays.fill(v,0);
	//	sum=0;
		Arrays.fill(lh,0);
	}
}