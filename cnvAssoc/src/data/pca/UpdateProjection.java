package data.pca;

import java.io.File;
import java.io.OutputStreamWriter;

import lc1.util.CompressDir;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public class UpdateProjection extends ProcessSNP {
double var;
	@Override
	public void refresh(){
		this.gs = new GramSchmidt();
	}
	
  public  GramSchmidt gs = new GramSchmidt();
    
	RealVector p = null;
	RealVector t= null;
	public void setP(RealVector p){
		this.p = p;
		
		this.t = new ArrayRealVector(p.getDimension());
		t.set(0);
	}
	//List<Double> norm= new ArrayList<Double>();
	
	public UpdateProjection(File[] dir, String nme, boolean str, GramSchmidt pcsToAdjust) {
		super(dir, nme, str, pcsToAdjust);
		
		// TODO Auto-generated constructor stub
	}

	@Override
	//j<1 indicates that we should write all columns
	public void print(int j, int[][] aliasRevs, CompressDir[] comp, String nme) {
		try{
		int start =0;
		for(int kk=0; kk<this.sze.length; kk++){
			OutputStreamWriter pw = comp[kk].getWriter(nme, true);
			int[] alias = aliasRevs[kk];
			int len = alias.length;
			  for(int k=0; k<len; k++){
			//	for(int j=0; j<gs.pcs.length; j++){
				  int k1 = alias[k];
				 if(j<0){
					 for(int j1=0; j1<gs.pcs.length; j1++){
					    pw.write(k1<0 ? "NaN": String.format("%5.3g", gs.pcs[j1].getEntry(k1+start)*1000.0).trim());
						pw.write(j1==gs.pcs.length-1 ? "\n":"\t");
					 }
				  }else{
					pw.write(k1<0 ? "NaN": String.format("%5.3g", gs.pcs[j].getEntry(k1+start)*1000.0).trim());
					pw.write("\n");
				}
			}
			comp[kk].closeWriter(pw);
			start+=sze[kk];
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	@Override
	public double getHistogram() throws Exception {
		t = t.mapDivide(this.t.getNorm());
		double dist = t.getDistance(p);
		p = t.copy();
		t.set(0);
		return dist;
	}
	public boolean addPC() {
		this.cnt=0;
		if(pcsToAdjust!=null && !this.pcsToAdjust.checkOrth(p) ) {
			return false;
		}
		if(!this.gs.checkOrth(p) ) return false;
		this.gs.update(this.p);
		return true;
		}

	

	@Override
	public RealVector process(double[] res) {
		double mean =0;
		double cnt=0;
		this.cnt++;
		double v;
		for(int k=0; k<res.length; k++){
			 v = res[k];
//			d.setEntry(k, v);
			if(!Double.isNaN(v)){
				mean+=v;
				cnt++;
			}
		}
		mean = mean/cnt;
		 var =0;
		for(int k=0; k<res.length; k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				var+=Math.pow(v-mean,2);
			}
			
		}
		if(var>0){
			var = Math.sqrt(var/cnt);
		for(int k=0; k<d.getDimension(); k++){
			 v = res[k];
			if(!Double.isNaN(v)){
				d.setEntry(k,  (v-mean)/var);
			}else{
				d.setEntry(k, 0);
			}
		}
		if(this.pcsToAdjust!=null){
			d = pcsToAdjust.removeProj(d);
		}
		d = gs.removeProj(d);
		
		//double dp = p.dotProduct(p);
		//double dp1 = d.dotProduct(d);
		double dp = p.dotProduct(d);
		//d.mapMultiplyToSelf(dp);
		RealVector proj = d.mapMultiply(dp);//p.projection(d);
		this.t = this.t.add(proj);
		}
		return null;
		
	}
	private double count(char[] st, char c) {
		int b =0;
		for(int j=0; j<st.length; j++){
			if(st[j]==c) b++;
		}
		return b;
	}

	public double process(RealVector realVector) {
		d = realVector;
		d = gs.removeProjLast(d);
		double dp = p.dotProduct(d);
		//d.mapMultiplyToSelf(dp);
		RealVector proj = d.mapMultiply(dp);//p.projection(d);
		this.t = this.t.add(proj);
		return dp;
		
	}

	
	

	

}
