package lc1.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.dp.data.collection.MatchedDistributionCollection.BackgroundDistribution;
import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;
import pal.math.UnivariateMinimum;

public class Multidimmax implements MultivariateFunction {
	
	List<double[]> l = new ArrayList<double[]>();
	List<double[]> initial = new ArrayList<double[]>();
	List<double[]> lower = new ArrayList<double[]>();
	List<double[]> upper = new ArrayList<double[]>();
	//double[] initial;
	List<String> name = new ArrayList<String>();
	List<int[]> alias = new ArrayList<int[]>(); //converts args to i
	//List<Integer> alias1=new ArrayList<Integer>(); //converts args to j  l.get(i)[j] 
	int ratios_ind = -1;
	
	final MatchedDistributionCollection mdc;
	public void add(String name, double[] d, double[] lower, double[] upper){
		for(int k=0; k<d.length; k++){
			alias.add(new int[] {l.size(),k});
		}
		if(name.equals("ratios")) ratios_ind = l.size();
		this.name.add(name);
		this.l.add(d);
		initial.add(d.clone());
		this.lower.add(lower);
		this.upper.add(upper);
	}
	
	public double[]  add(String name, double d, double lower, double upper, int len, boolean add){
		double[] toadd = new double[len];
		double[] toaddlower = new double[len];
		double[] toaddupper = new double[len];
		Arrays.fill(toadd, d);
		Arrays.fill(toaddlower, lower);
		Arrays.fill(toaddupper, upper);
		if(add) this.add(name, toadd, toaddlower, toaddupper);
		return toadd;
		
	}
	
	
	 public Multidimmax(MatchedDistributionCollection mdc){
		 this.mdc = mdc;
	 }
	
	
	public void updateBounds(){
		if(ratios_ind<0) return;
		double[] ratios = l.get(ratios_ind);
		double[] lower = this.lower.get(ratios_ind);
		double[] upper = this.upper.get(ratios_ind);
			for(int k=0; k<ratios.length; k++){
				if(k==Constants.maxPloidy1()){
					lower[k]=1.0;
					upper[k]=1.0;
				}else{
					  lower[k] = k==0 ? 0.001: ratios[k-1]+0.0001;
					  upper[k] = k<ratios.length-1 ? ratios[k+1]-0.0001 : ratios[k]+0.5;
					  if(k<Constants.maxPloidy1()){
						  lower[k] = Math.min(lower[k], 0.99999999);
						  upper[k] = Math.min(upper[k], 0.99999999);
					  }
					  else if(k>Constants.maxPloidy1()){
						  lower[k] = Math.max(lower[k], 1.0000001);
						  upper[k] = Math.max(upper[k], 1.0000001);
					  }
				  if(Math.signum(lower[k]-1) != Math.signum(upper[k]-1)){
					  throw new RuntimeException("!!");
				  }
				}
			}
			System.err.println("h");
		}
	 
	
	public void print(){
		System.err.println("new vals");
		for(int k=0; k<this.l.size(); k++){
		
			System.err.print(this.name.get(k));
			double[] v = l.get(k);
			for(int j=0; j<v.length; j++){
				System.err.print(String.format("%5.3g ", v[j]));
			}
		}
	}

	double[] vals;
	public double[] vals() {
		int cum=0;
		if(vals==null) vals = new double[this.alias.size()];
		for(int i=0; i<this.l.size(); i++){
			int len = l.get(i).length;
			System.arraycopy(l.get(i),0, vals, cum, len);
			cum +=len;
		}
		// TODO Auto-generated method stub
		return vals;
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
		//System.err.println(arg0[0]+" "+res);
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

}
