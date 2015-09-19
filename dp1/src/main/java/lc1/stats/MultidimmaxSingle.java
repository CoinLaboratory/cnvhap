package lc1.stats;


import java.util.ArrayList;
import java.util.List;

import lc1.dp.data.collection.MatchedDistributionCollection;
import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;

public class MultidimmaxSingle extends AbstractMultiDimMax {
	
	public MultidimmaxSingle(MatchedDistributionCollection mdc){
		 super(mdc);
	 }
	
	public double evaluate1(double arg0, int k) {
		int cum=0;
		double lp =0;
		inner: for(int i=0; i<this.l.size(); i++){
			int len = l.get(i).length;
			if( k < cum+len){
				 l.get(i)[k-cum] = arg0;
				 break inner;
			}
		//	for(int k=0; k<len; k++){
		//	lp+=-Math.abs(arg0[k+cum] - initial.get(i)[k])/1e10;
		//	}
			//System.arraycopy(arg0,cum, l.get(i), 0, len);
			cum +=len;
		}
		double res =  lp+mdc.evaluate(k);
		//System.err.println(arg0.length+" "+arg0[0]+" "+res);
		return res;
	}
	
	@Override
	public void add(String name, double[] d, double[] lower, double[] upper) {
		super.add(name,d,lower,upper);
		for(int k=0; k<d.length; k++){
			inn.add(new Inner(inn.size()));
		}
	}
	
	final List<Inner> inn = new ArrayList<Inner>();
	
class Inner implements UnivariateFunction{
	final int k;
	Inner(int arg0){
		this.k = arg0;
	}
	@Override
	public double getLowerBound() {
		int[] pos = alias.get(k);
		return lower.get(pos[0])[pos[1]];
	}

	@Override
	public double getUpperBound() {
		int[] pos = alias.get(k);
		return upper.get(pos[0])[pos[1]];
	}

	@Override
	public double evaluate(double arg0) {
		// TODO Auto-generated method stub
		return evaluate1(arg0,k);
	}
}
@Override
	public void minimize() {
		UnivariateMinimum mvm =new UnivariateMinimum(); // new ConjugateGradientSearch();//n
		super.vals();
		for(int i=0; i<this.inn.size(); i++){
			//System.err.println(i);
		   mvm.findMinimum(vals[i], (UnivariateFunction) inn.get(i), 2);
		}
		
		
	}


}
