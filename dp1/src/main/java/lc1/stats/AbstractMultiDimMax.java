package lc1.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.data.collection.MatchedDistributionCollection;
import pal.math.MultivariateFunction;

public abstract class AbstractMultiDimMax {

	protected List<double[]> l = new ArrayList<double[]>();
	List<double[]> initial = new ArrayList<double[]>();
	protected List<double[]> lower = new ArrayList<double[]>();
	protected List<double[]> upper = new ArrayList<double[]>();
	List<String> name = new ArrayList<String>();
	protected List<int[]> alias = new ArrayList<int[]>();
	int ratiosL_ind;
	int ratiosR_ind = -1;
	protected final MatchedDistributionCollection mdc;
	double[] vals;
	
	public AbstractMultiDimMax(MatchedDistributionCollection mdc){
		 this.mdc = mdc;
	 }
public abstract void minimize();	
	public void add(String name, double[] d, double[] lower, double[] upper) {
		for(int k=0; k<d.length; k++){
			alias.add(new int[] {l.size(),k});
		}
		if(name.equals("ratiosL")) ratiosL_ind = l.size();
		else if(name.equals("ratiosR")) ratiosR_ind = l.size();
		this.name.add(name);
		this.l.add(d);
		initial.add(d.clone());
		this.lower.add(lower);
		this.upper.add(upper);
	}

	public double[] add(String name, double d, double lower, double upper,
			int len, boolean add) {
				double[] toadd = new double[len];
				double[] toaddlower = new double[len];
				double[] toaddupper = new double[len];
				Arrays.fill(toadd, d);
				Arrays.fill(toaddlower, lower);
				Arrays.fill(toaddupper, upper);
				if(add) this.add(name, toadd, toaddlower, toaddupper);
				return toadd;
				
			}

	public void updateBounds() {
		 this.updateBounds(ratiosL_ind, 0);
		 this.updateBounds(ratiosR_ind, 1);
	 }

	public void updateBounds(int ratiosL_ind, double min) {
	if(ratiosL_ind<0) return;
	{
	double[] ratios = l.get(ratiosL_ind);
	double[] lower = this.lower.get(ratiosL_ind);
	double[] upper = this.upper.get(ratiosL_ind);
		for(int k=0; k<ratios.length; k++){
			{
				  lower[k] = k==0 ? min+0.001: ratios[k-1]+0.0001;
				  upper[k] = k<ratios.length-1 ? ratios[k+1]-0.0001 : ratios[k]+0.5;
			}
		}
	}
	
		System.err.println("h");
	}

	public void print() {
		System.err.println("new vals");
		for(int k=0; k<this.l.size(); k++){
		
			System.err.print(this.name.get(k));
			double[] v = l.get(k);
			for(int j=0; j<v.length; j++){
				System.err.print(String.format("%5.3g ", v[j]));
			}
		}
	}

	

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

}