package data.pca;

import java.io.File;

import lc1.util.Compressor;

import org.apache.commons.math.linear.RealVector;

public class CalcCovariance extends ProcessZip {

	public CalcCovariance(File dir, String[] lrrst, String[] strtype,String[] chroms) {
		super(dir, lrrst,strtype, chroms);
		// TODO Auto-generated constructor stub
	}
	public void makeHists(String nme, String[] lrrst){
		 int len = lrrst.length;
		 this.l1 = new Hist[len];
		 for(int k=0; k<len; k++){
			 this.l1[k] = new Hist(dir, nme+".hist.txt."+lrrst[k], this.type[k]);
		 }
		// this.l1_sex = new Hist(dir, nme+"sex.hist.txt");
	 }
	
	public void run() throws Exception{
		this.runInner(0, 0);
	}

}
