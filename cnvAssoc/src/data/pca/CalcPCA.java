package data.pca;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.ZipFile;

import lc1.util.Compressor;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.linear.SingularValueDecompositionImpl;
import org.apache.commons.math.stat.correlation.Covariance;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import assoc.Constants;

public class CalcPCA extends ProcessZip {

static boolean check = false;


//int max_cnt = 100;//1*1000*1000;//Integer.MAX_VALUE;



	public CalcPCA(File dir, String[] nmes, String[] lrrst,String[] strtype, String[] chroms, String[] thin, String thresh, boolean merge,Integer no_pcs) {
		super(dir, nmes, lrrst, strtype, chroms, thin, thresh, merge,no_pcs);
		
		// TODO Auto-generated constructor stub
	}
	
	
	
	public void makeHists(String nme, String[] lrrst){
		int len = lrrst.length;
		 this.l1 = new UpdateProjection[len];
		 for(int k=0; k<len; k++){
			 this.l1[k] =new UpdateProjection(dir, nme+".hist.txt"+lrrst[k], this.type[k], this.gs[k]);
		 }
	 }
	
	//List<RealVector> pcs_auto = new ArrayList<RealVector>(0);
	//List<RealVector> pcs_sex = new ArrayList<RealVector>(0);
	RealVector p_auto = null;
	
	
	
	
	public void run() throws Exception{
		this.refresh();
		for(int k=0; k<no_reps; k++){
		
			if(!this.run1(k)) {
				break;
			}
		}
		for(int j=0; j<l1.length; j++){
			if(AbstractProcessZip.singleEntry){
				l1[j].print(-1,aliasRev, this.output, this.currentRegion.toString()+"PC");
				this.snps[j].write(currentRegion.line()+"PC"+"\n");
			}else{
				for(int k=0; k<no_reps; k++){
					l1[j].print(k,aliasRev, this.output, this.currentRegion.toString()+"PC"+k);
					this.snps[j].write(currentRegion.line()+"PC"+k+"\n");
				}
			}
		}
		
	}
	
	public void refresh(){
		this.p_auto = null;
		
		this.makeP();
		for(int i=0; i<this.l1.length; i++){
			this.cache[i] = new ArrayList<RealVector>();
			this.l1[i].refresh();
			 ((UpdateProjection)l1[i]).setP(p_auto.copy());
		}
		
	}
	
	
	protected void makeP() {
		if(p_auto==null){
			try{
			
				double[] d = new double[Constants.sum(sze)];
				double sum=0;
				for(int k=0; k<d.length; k++){
					double v = Math.random();
					d[k] = v;
					sum+=v;
				}
				
				//Arrays.fill(d,1.0);
				p_auto = new ArrayRealVector(d);
				p_auto.mapDivide(p_auto.getNorm());
			//	p_sex = new ArrayRealVector((new Dirichlet(d,u)).sample());
				for(int k=0; k<l1.length; k++){
				 ((UpdateProjection)l1[k]).setP(p_auto.copy());
				}
				// ((UpdateProjection)l1_sex).setP(p_sex, pcs_sex);
				}catch(Exception exc){
					exc.printStackTrace();
				}
		}
		
	}



	public boolean run1(int index) throws Exception{
	System.err.println("rep "+index);
	boolean updated = true;
	double[] t_norm = new double[l1.length];
		for(int ii=0; ;ii++){
			if(tocache && (ii>0 || index > 0)){
				this.runInnerQuick();
			}else{
				this.runInner(index, ii);
			}
			double maxdiff = 0;
			//if(l1[0].count()==0) return false;
		 for(int k=0; k<l1.length; k++){
			t_norm[k] = ((UpdateProjection)l1[k]).t.getNorm();
			 double diff =  l1[k].getHistogram( );
			 if(diff> maxdiff) maxdiff = diff;
			 System.err.print(String.format("%5.3g",diff).trim()+"\t");
		 }
	 System.err.println();
	
	 if(maxdiff <thresh) {
		 for(int k=0; k<l1.length; k++){
		 updated = updated && ((UpdateProjection)l1[k]).addPC();
		 }
		/* if(this.weights!=null){
			 for(int k=0; k<weights.length; k++){
				 for(int j=0; j<weights[k].length; j++){
					 weights[k][j] = Math.abs(weights[k][j]/t_norm[k]);
					 this.snp_list.set(j,snp_list.get(j)+"\t"+weights[k][j]);
				 }
			 }
			 
		 }*/
		 
		 break;
	 }
	}
		System.err.println();
		return updated;
	}



	private File[] getFile(File[] dir, String string) {
		// TODO Auto-generated method stub
		return null;
	}

	
}
