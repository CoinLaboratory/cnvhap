/**
 * 
 */
package lc1.dp.illumina;

import java.io.Serializable;
import java.util.Arrays;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;

class BasisFunction implements Serializable{
	//DoubleMatrix1D matr;
	double[] d;
	final EmissionStateSpace emstsp;
	final double bg;
	String[] name;
	
	public boolean equals(Object o){
		BasisFunction bf1 = (BasisFunction)o;
		return  bf1.name.length==name.length && Arrays.asList(name).equals(Arrays.asList(bf1.name));
	}
	
	public Double CN(Double cn, Double cb, Double bg){
		
		return cn;
	}
	public Double lCN(Double cn, Double cb, Double bg){
	
		return cn==0 ? Constants.lrr0() :log2((double)cn/bg);
	}
	public Double lCN2(Double cn, Double cb, Double bg){
		return cn==0 ? Constants.lrr0(): Math.pow(log2((double)cn/bg),2);
	}
	
	public Double lCN3(Double cn, Double cb, Double bg){
		return cn==0 ? Constants.lrr0() :Math.pow(log2((double)cn/bg),3);
	}

	public Double noA(Double cn, Double cb, Double bg){
		return cn - cb;
	}

	public Double noB(Double cn, Double cb, Double bg){
		return cb;
	}
	public Double zero(Double cn, Double cb, Double bg){
		if(cn==0) return 1.0; 
		else return 0.0;
	}
	
	BasisFunction(EmissionStateSpace emstsp, double bg, String[] nme){
		 d = new double[nme.length];
		d[0] = 1;
		this.emstsp = emstsp;
		this.bg = bg;
		this.name = nme;
		//this.matr = new DenseDoubleMatrix1D(d);
		//matr.setQuick(0, 1);
	}
	public double[] getVals(int index){
		int cn =emstsp.getCN(index);
		int cb = emstsp.getBCount(index);
		return getVals(cn,cb);
	}
	public double[] getVals(double cn, double cb) {
		//double baf = (double) cb / (double) cn;
		try{
			Object[] args = new Object[] {cn,cb,bg};
		double cnv = cn==0 ? -1e10: log2((double)cn/bg);
		for(int k=1; k<this.name.length; k++){	
		  d[k] = (Double)this.getClass().getMethod(name[k], new Class[] {Double.class, Double.class, Double.class}).invoke(this, args);
		}
		}catch(Exception exc){
			Arrays.fill(d,0);
			d[1] = 1;
			exc.printStackTrace();
		}
		return d;
	}
	static double log2 = Math.log(2);
	private double log2(double e) {
		return Math.log(e)/log2;
	}
	

	public int length() {
		return d.length;
	}
	public String toString(){
		return Arrays.asList(name).toString();
	}
	public String getName(int k) {
		return name[k];
	}
}