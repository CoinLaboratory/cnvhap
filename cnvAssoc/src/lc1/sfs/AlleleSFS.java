package lc1.sfs;


import java.util.Arrays;
import java.util.List;

import lc1.sfs.freqfunct.Funct;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;


public class AlleleSFS extends SFS {

	

	
	/*pop1, pop2 are the  total number of indiv */
	public AlleleSFS(int[] pops, int ploidy, Funct[] functs, String nme){
		super(pops,ploidy, functs.length,nme+ "_allele"+getName(pops));
		System.err.println(Arrays.asList(functs));
		System.err.println("using functs ");
		
		for(int k=0; k<pop.length; k++){
			pop[k] = pops[k]*ploidy;
			lim[k] = pop[k]+1;
			hweprobs[k] = getHWEMatrix(lim[k]);
			priors[k] = new Likelihood(lim[k]);
		}
		AFS = AbstractMultiDM.getMultiDM(lim, lim.length,1.0);
		//for(int k=0; k<fst.length; k++){
		//	fst[k] = 
		//}
		while(true){
			//	for(int k=0; k<fst.length; k++){
			double[] v = calc(functs,  pop);
			for(int j=0; j<this.pmap.length; j++){
				pmap[j].put(v[j], 0.0);
			}
					AFS.set(it.inds,v);
				//}
				if(!it.incr()) break;
		}
		
	}
	
	protected double[] calc(Funct[] functs, int[] pop2) {
		double[] res = new double[functs.length];
		int[] inds = this.it.inds;
		for(int j=0; j<inds.length; j++){
			v[j] = (double) inds[j]/ (double)pop2[j];
		}
		//System.err.println(v[0]+" "+pop2[0]);
		for(int k=0; k<functs.length; k++){
			res[k] = functs[k].calc(v, pop2);
		/*	if(res[k]<0){
				throw new RuntimeException("exc");
			}*/
		}
		return res;
	}

public AlleleSFS(int[] sze, int ploidy, String string, List<String >nme, String name) {
		this(sze, ploidy, Funct.getFuncts(string.split(":"),nme),name);
	}
	
protected RealMatrix getHWEMatrix(int pop) {
		RealMatrix re = new Array2DRowRealMatrix(3,pop);
		for(int i=0; i<pop; i++){
			double p = (double)i/(double)pop;
			double p2 = Math.pow(p, 2);
			double q2 = Math.pow(1-p, 2);
			double twopq = 2*(1-p)*p;
			re.setEntry(0, i, p2);
			re.setEntry(1, i, twopq);
			re.setEntry(2, i, q2);
		}
		return re;
}

}
