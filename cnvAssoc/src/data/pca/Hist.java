/**
 * 
 */
package data.pca;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealVector;


class Hist extends ProcessSNP{
	public Hist(File[] dir, String nme, boolean str) {
		super(dir, nme, str);
		// TODO Auto-generated constructor stub
	}
	List<Double> l1 = new ArrayList<Double>();
	
	RealVector rv = new ArrayRealVector(2);
	
	/* (non-Javadoc)
	 * @see data.ProcessSNP#process(java.lang.String[])
	 */
	public RealVector process(String[] res) {
		double mean =0;
		double cnt=0;
		for(int k=0; k<res.length; k++){
			String st = res[k];
			if(st.indexOf('N')>=0) st = "NaN";
			double v = Double.parseDouble(st);
			d.setEntry(k, v);
			if(!Double.isNaN(v)){
				mean+=v;
				cnt++;
			}
		}
		mean = mean/cnt;
		double var =0;
		for(int k=0; k<res.length; k++){
			double v = d.getEntry(k);
			if(!Double.isNaN(v)){
				var+=Math.pow(v-mean,2);
			}
		}
		var = var/cnt;
		rv.setEntry(0, mean);
		rv.setEntry(1,var);
			l1.add(var);
		return rv;
	}
	/* (non-Javadoc)
	 * @see data.ProcessSNP#getHistogram()
	 */
	public double getHistogram() throws Exception {
		PrintWriter[] pw = new PrintWriter[this.sze.length];
		for(int k=0; k<sze.length; k++){
			pw[k] = new PrintWriter(new FileWriter(outfile[k]));
		}
		// TODO Auto-generated method stub
		 Collections.sort(l1);
		 double step = 0.01;
		 double nextp=0;
		 nextp+=step;
		for(int k=0; k<l1.size(); k++){
			double perc = (double) k / (double) l1.size();
			if(perc>=nextp){
				for(int kk=0; kk<sze.length; kk++){
				pw[kk].println(String.format("%5.3g",perc )+"\t"+String.format("%5.3g",l1.get(k)));
				}
				nextp+=step;
			}
		}
		for(int k=0; k<sze.length; k++){
			pw[k].close();
		}
		
		return 0;
	}
	
	
}