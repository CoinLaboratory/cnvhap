/**
 * 
 */
package lc1.dp.illumina;

import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.reflect.Field;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;

public class Prior implements Serializable{
	 public Prior(Prior prior) {
		this.priorR = copy(prior.priorR);
		this.priorRVar = copy(prior.priorRVar);
	}
	 
	 public boolean equals(Object obj){
		 Prior p1 = (Prior) obj;
		 return priorR.equals(p1.priorR) && priorRVar.equals(p1.priorRVar);
	 }
	public Prior(int ik, int[][] indices, EmissionStateSpace emstsp, int len) {
		makePriors(indices.length);
		//double[] rmean = Constants.r_mean(ik);
		//int len = rmean.length;
		for(int i=0; i<indices.length; i++){
			int cn = emstsp.getCN(indices[i][0]);
			makePriorNormal(len, i,ik);
		}
	}
	public void makePriors(int length) {
		this.priorR =new DoubleMatrix1D[length];
		this.priorRVar =new DoubleMatrix1D[length];
		
	}
	public void makePriorNormal(int len, int i,int ik){
		priorR[i] = new DenseDoubleMatrix1D(len);
		priorRVar[i] = new DenseDoubleMatrix1D(len);
		copy(Constants.r_mean(ik,i),priorR[i]);
		copy(Constants.r_var(ik,i), priorRVar[i]);
	}
	public void copy(double[] r_var, DoubleMatrix1D doubleMatrix1D) {
		for(int k=0; k<r_var.length && k<doubleMatrix1D.size(); k++){
			doubleMatrix1D.set(k, r_var[k]);
		}
		
	}
	
	
	
	protected DoubleMatrix1D[] copy(DoubleMatrix1D[] priorRho2) {
		 if(priorRho2==null) return null;
		else return priorRho2.clone();
	}
	public DoubleMatrix1D[] priorR, priorRVar;
	
	//  public DoubleMatrix1D 
	public void print(PrintWriter pw) {
		try{
	Field[] f = 	Prior.class.getFields();
		for(int i=0; i<f.length; i++){
			if(f[i]==null) pw.println("null");
			else pw.println(f[i].getName()+"\t"+getString((DoubleMatrix1D[])f[i].get(this)));
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public void print() {
		try{
	Field[] f = 	Prior.class.getFields();
		for(int i=0; i<f.length; i++){
			if(f[i]==null) System.err.println("null");
			else  System.err.println(f[i].getName()+"\t"+getString((DoubleMatrix1D[])f[i].get(this)));
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	private String getString(DoubleMatrix1D[] d) {
		StringBuffer sb = new StringBuffer();
		if(d!=null){
		for(int i=0; i<d.length; i++){
			if(d[i]==null){
				sb.append("null");
			}
			else{
				for(int j=0; j<d[i].size(); j++){
					sb.append(String.format("%5.3g", d[i].getQuick(j)));
					if(j<d[i].size()-1) sb.append(";");
				}
			}
			if(i<d.length-1) sb.append(":");
		}
		}
		return sb.toString();
	}
	public void parse(String st) {
		if(!st.startsWith("prior")) return ;
		try{
		String[] str = st.split("\t");
		String[] str1 = str[1].split(":");
		DoubleMatrix1D[] mat = new DoubleMatrix1D[str1.length];
		Prior.class.getField(str[0]).set(this, mat);
	
		for(int i=0; i<mat.length; i++){
			if(str1[i].equals("null")) continue;
			else{
				mat[i] = getMatr(str1[i].split(";"));
			}
		}
		}catch(Exception exc){
			
		}
	}
	private DoubleMatrix1D getMatr(String[] split) {
		DoubleMatrix1D res = new DenseDoubleMatrix1D(split.length);
		for(int i=0; i<split.length; i++){
			res.setQuick(i, Double.parseDouble(split[i]));
		}
		return res;
	}

		
 }