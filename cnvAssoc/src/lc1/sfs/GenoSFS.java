package lc1.sfs;


import java.util.Arrays;
import java.util.List;

import lc1.sfs.genofunct.GFunct;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

/*note genotypes are arranged in order first of number of count(AB), then count(AA) */
public class GenoSFS extends SFS {

//	int[][] het_offset;  //this is an index which says which position for starting the list of genotypes with hets = i
//	int[][][] allele_counts; //[pop_ind, cell_index, AA:AB:BB] count

	
	//private int getIndex(int nohets, int no_aa){
	//	return het_offset[nohets]+no_aa;
	//}
	
/*pop1, pop2 are the  total number of indiv */
	public GenoSFS(int[] pops, int ploidy, GFunct[] functs, String name){
		super(pops,ploidy, functs.length, name+"_geno"+getName(pops));
		
		double[] d = new double[3];
		for(int k=0; k<pops.length; k++){
			pop[k] = (pops[k]+2)*(pops[k]+1)/2;
			lim[k] = pop[k]; //CHECK !!!
			this.priors[k] = new Likelihood(pop[k]);
			RealMatrix hwek =new Array2DRowRealMatrix(3,pop[k]);
			hweprobs[k]= hwek;
			int ind=0;
			for(int i=0; i<=pops[k]; i++){ //i is number of hets AB
				d[1] = ((double) i) /  (double)pops[k];
				for(int j=0; j<=pops[k]-i; j++){ //j is number of BB
					d[2] = ((double) j) /  (double)pops[k];
					d[0] = 1 - d[1] - d[2];
					hwek.setColumn(ind, d.clone());
					ind++;
				}
			}
			if(ind!=pop[k]) {
				throw new RuntimeException("!! "+ind+" "+pop[k]+" "+pops[k]);
			}
		}
		System.err.println("using functs ");
		System.err.println(Arrays.asList(functs));
		AFS = AbstractMultiDM.getMultiDM(lim, lim.length,1.0);
		while(true){
			double[] v = calc(functs,  pop);
			for(int j=0; j<this.pmap.length; j++){
				pmap[j].put(v[j], 0.0);
			}
					AFS.set(it.inds,v);
					
					if(!it.incr()) break;
		}
		
	}
	
	
	
	



	protected double[] calc(GFunct[] functs,  int[] pop2) {
		double[] res = new double[functs.length];
		int[] inds = this.it.inds;
		for(int j=0; j<inds.length; j++){
			double[] re = this.hweprobs[j].getColumn(inds[j]);
			double tot = re[0]+re[1]+re[2];
			v[j] =(re[1] + 2*re[0]) / (2*tot);
			this.v_het[j] = re[1]/tot;
		}
		
		for(int k=0; k<functs.length; k++){
			res[k] = functs[k].calc(v, v_het,pop2);
		/*	if(res[k]<0){
				throw new RuntimeException("exc");
			}*/
		}
		return res;
	}
public GenoSFS(int[] sze, int ploidy, String string, List<String >nme, String name) {
		this(sze, ploidy, GFunct.getFuncts(string.split(":"),nme),name);
	}
	


}
