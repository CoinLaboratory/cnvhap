package lc1.dp.illumina;

import lc1.dp.emissionspace.EmissionStateSpace;

public class BasisFunctionRB extends BasisFunction {

	public Double BAF(Double cn, Double cb, Double bg){
		return  cn==0 ? 0 : (double) cb / (double) cn;
	}
	public Double BAF1(Double cn, Double cb, Double bg){
		double baf = BAF(cn,cb,bg);
		return (baf*(1-baf))/norm1;
	}
	static double norm1 =Math.pow(0.5,2); 
	static double norm2 = (Math.pow(0.85-0.5,2)*(0.85*(1-0.85)));
	public Double BAF2(Double cn, Double cb, Double bg){
		double baf = BAF(cn,cb,bg);
		return (Math.pow(baf-0.5,2)*(baf*(1-baf)))/ norm2;
	}
	
	
	
	BasisFunctionRB(EmissionStateSpace emstsp, double bg, String[] nme) {
		super(emstsp, bg,nme);
		this.d = new double[nme.length];
		d[0] = 1;
		this.name =nme; 
			//"C:CN:CN^2:BAF:(BAF-1)*(BAF-0)".split(":");
		// TODO Auto-generated constructor stub
	}

}
