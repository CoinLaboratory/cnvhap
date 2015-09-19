package lc1.stats;


import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.dp.data.collection.MatchedDistributionCollection.BackgroundDistribution;
import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;
import pal.math.UnivariateMinimum;

public class Multidimmax extends AbstractMultiDimMax implements MultivariateFunction {
	
	public Multidimmax(MatchedDistributionCollection mdc){
		 super(mdc);
	 }
	 @Override
	public double evaluate(double[] arg0) {
		int cum=0;
		double lp =0;
		
		for(int i=0; i<this.l.size(); i++){
			int len = l.get(i).length;
		//	for(int k=0; k<len; k++){
		//	lp+=-Math.abs(arg0[k+cum] - initial.get(i)[k])/1e10;
		//	}
			System.arraycopy(arg0,cum, l.get(i), 0, len);
			cum +=len;
		}
		double res =  lp+mdc.evaluate();
		//System.err.println(arg0.length+" "+arg0[0]+" "+res);
		return res;
	}

	@Override
	public double getLowerBound(int arg0) {
		int[] pos = alias.get(arg0);
		return this.lower.get(pos[0])[pos[1]];
	}

	@Override
	public int getNumArguments() {
		return this.alias.size();
	}

	@Override
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getUpperBound(int arg0) {
		int[] pos = alias.get(arg0);
		return this.upper.get(pos[0])[pos[1]];
	}

	public void minimize() {
		 MultivariateMinimum mvm =new OrthogonalSearch(); // new ConjugateGradientSearch();//n
		   mvm.findMinimum(this, vals(), 1, 3);
		
		
	}

}
