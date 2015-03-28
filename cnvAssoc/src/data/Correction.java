/**
 * 
 */
package data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealVector;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;

class Correction{

final boolean centralised = true;
final boolean standardised = true;
	final double[] prod;
	final boolean[] na;
	private double[] adjust(RealVector v, double mean, double stddev){
		RealVector v1 = v.mapAdd(-mean).mapDivide(stddev);
		for(int k=0; k<vec.length; k++){
		//for(int k=0; k<mats.length; k++){
			RealVector basis = this.vec[k];
		
		
			double denom =  basis.dotProduct(basis);
			double num = v.dotProduct(basis);
			prod[k] =num==0 ? 0 : num/denom;
		}
		for(int k=0; k<prod.length; k++){
				RealVector basis = this.vec[k];
				v1 = v1.subtract(basis.mapMultiply(prod[k]));
		}
		
		double[] d = v1.mapMultiply(stddev).mapAdd(mean).getData();
		for(int k=0; k<d.length; k++){
			if(na[k]){
				d[k] = Double.NaN;
			}
		}
		return d;
	}
/*
	public static List<String> read(ZipFile zf,  String probe) throws Exception {
		ZipEntry ent = zf.getEntry(probe);
		if(ent==null) return null;
		List<String> res = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st);
		}
		br.close();
		return res;
		
	}*/
	private RealVector getVector(Double[] split) {
		double[] d = new double[split.length];
		for(int k=0; k<d.length; k++){
			d[k] = split[k];
		}
		return new ArrayRealVector(d);
	}
String[][]input;
	public Correction(File corrections2, List<String> indiv,int len1, int max) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(corrections2));
		this.input = new String[indiv.size()][len1];
		List<Double>[] vecs=null;
		String st;
		int k=0;
		for(k=0; (st=br.readLine())!=null; k++){
		
			String [] str = st.split("\\s+");
		
		if(k==0){
				vecs = new List[Math.min(str.length-1,max)];
				for(int j=0; j<vecs.length; j++){
					vecs[j] = new ArrayList<Double>();
				}
			}
		if(!indiv.get(k).split("\t")[0].equals(str[0])) throw new RuntimeException("!!");
		//	indv.add(str[0]);
			for(int j=0; j<vecs.length; j++){
				vecs[j].add(Double.parseDouble(str[j+1]));
			}
		}
		
		int len = k;
		this.prod = new double[vecs.length];
		this.na = new boolean[len];
		br.close();
		this.vec = new RealVector[vecs.length];
		boolean orthog = true;
		outer:for(int j=0; j<vec.length; j++){
			vec[j] = this.getVector(vecs[j].toArray(new Double[0]));
			for(int kk=0; kk<j; kk++){
				double dp = vec[j].dotProduct(vec[kk]);
				if(Math.abs(dp)>1e-5) {
					orthog = false;
					//break outer;
				}
			}
		}
		if(!orthog){
			System.err.println("WARNING WAS NOT ORTHOGONAL");
			this.orthogonalise();
		}
		double[]d = new double[len];
		this.inp = new ArrayRealVector(d);
	}
	private void orthogonalise() {
		// private DoubleMatrix2D getGramSchmidt(List<DoubleMatrix1D> x2) {
			 //int len = this.vec[0].getDimension();
			
//			DoubleMatrix2D inp1 = this.getNewNullspace(inp, x2); //first column is target
			for(int i=1; i<vec.length; i++){
				double[] weight = new double[i];
				for(int j=0; j<i; j++){
					weight[j] = this.vec[j].dotProduct(this.vec[i]) / this.vec[j].dotProduct(this.vec[j]);
					
				}
				for(int j=0; j<i; j++){
					vec[i] = vec[i].subtract(this.vec[j].mapMultiply(weight[j]));
				}
			}
	}
	File corrections;
	RealVector[] vec;
	RealVector inp;
	public void correct(int lrr_id) {
		Arrays.fill(na,false);
		double cnt=0;
		double sum=0;
		for(int k1=0; k1<input.length; k1++){
		double v = Double.parseDouble(input[k1][lrr_id]);
		if(Double.isNaN(v) && centralised){
			this.inp.setEntry(k1,0 );
			na[k1] = true;
		}
		else{
			this.inp.setEntry(k1,  v);
			cnt++;
			sum+=v;
		}
		}
		double mean = sum/cnt;
		double stddev = 1;
		if(this.standardised){
			 sum=0;
			for(int k1=0; k1<input.length; k1++){
				double v = inp.getEntry(k1);
				if(!Double.isNaN(v) ){
						sum+=Math.pow(v-mean,2);
				}
			}
			stddev = Math.sqrt(sum/cnt);
		}
		if(!centralised){
			for(int k=0; k<na.length; k++){
				if(na[k]) inp.setEntry(k, mean);
			}
			mean=0;
		}
	
		double[] d = adjust(inp,mean, stddev);
		 
		for(int k=0; k<input.length; k++){
			String[] inp = input[k];
			inp[lrr_id] = String.format("%5.3g", d[k]);
		
		}
	}
	
	
	public String[][] correct(BufferedReader br, int lrr_id) throws Exception{
		String st = "";
		for(int k=0;(st =br.readLine() )!=null;k++){
			int k1 = k;
				input[k1]=st.split("\t");
		}
		
          correct( lrr_id);
		 br.close();
			return this.input;
	}
	
	public void correct1(List<String[]> l, int lrr_id, boolean centralised) throws Exception{
		String st = "";
		for(int k=0;k<l.size();k++){
				input[k]=l.get(k);
		}
		
          correct( lrr_id);
		for(int k=0; k<l.size(); k++){
			l.set(k, input[k]);
		}
		//System.err.println("corrected");
	}
	
	public void correct(List<String> l, int lrr_id) throws Exception{
		for(int k=0;(k<l.size() );k++){
			int k1 = k;
				input[k1]=l.get(k).split("\t");
		}
		
          correct( lrr_id);
          for(int k=0;(k<l.size() );k++){
				l.set(k,getString(input[k]));
			}
	}
	private String getString(String[] strings) {
		StringBuffer sb  = new StringBuffer(strings[0]);
		for(int k=1; k<strings.length; k++){
			sb.append("\t"+strings[k]);
		}
		return sb.toString();
	}

	
}