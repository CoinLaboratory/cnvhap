package data.cc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.linear.RealVector;


/* finds correlated vectors to plate */
public class PlateCorrection {

	static int start = 0;
	static int end = 1*1000*1000; //Integer.MAX_VALUE;
	static boolean plot = false;
	
	public static void main(String[] args) {
		try {
			File base = new File(args[0]);
			String[] type = args[1].split(":");
			File[][] out = new File[type.length][2];
			List<String> chr = Arrays.asList(args[3].split(":"));
			String logr_ = args[4];
			String[] thi = args[5].split(":");
			double[] thin_ = new double[2];
			if(thi.length==1){
				double v = Double.parseDouble(thi[0]);
				thin_[0] = Math.max(0, 0.5 -v);
				thin_[1] = Math.min(1, 0.5+v);
			}
			else{
				thin_[0] = Double.parseDouble(thi[0]);
				thin_[1] = Double.parseDouble(thi[1]);
			}
			String[] reg = args[6].split(":");
			File[] dir = new File[type.length];
			for (int k = 0; k < dir.length; k++) {
				dir[k] = new File(base, type[k]);
				out[k] =new File[2];
				for(int k1 =0; k1<2; k1++){
					out[k][k1] = new File(base, type[k]+"/"+(k1==0 ? args[2]+".txt":"plate_vecs.txt"));
				}
			}
			PlateCorrection pl = new PlateCorrection(dir,  chr, logr_.split(":"), thin_);
			if(pl.plate!=null){
			pl.run(Double.parseDouble(reg[0]), Double.parseDouble(reg[1]));
			for(int k=0; k<out.length; k++){
				pl.print(out[k], pl.indivs[k]);
			}
			}else{
				pl.runPCA();
				
			}
			
		} catch (Exception exc) {
			exc.printStackTrace();
			System.exit(0);
		}
		System.exit(0);
	}

	private void runPCA() {
	this.lrr.getPCs();
		
	}

	Data plate;
	Data lrr;
	//File out;
CC d1;
	public PlateCorrection(File[] dir, List< String> chr, String[] logr_,
			double[] var_bound) throws Exception {
		Data.perc = var_bound;
		boolean pca =true;
		File[] plate_f = new File[dir.length];
		//this.out = out;
		File[] lrr_f = new File[dir.length];
		for(int k=0; k<logr_.length; k++){
			logr_[k] = logr_[k].replaceAll("_", " ");
		}
		int[] start = new int[] { PlateCorrection.start, PlateCorrection.end };
	
		for (int k = 0; k < plate_f.length; k++) {
			plate_f[k] = new File(dir[k], "plate.zip");
			lrr_f[k] = new File(dir[k], chr + ".zip");
		}
		indivs = new List[plate_f.length];
		for(int k=0; k<plate_f.length;k++){
			Data lr = new Data(lrr_f[k],logr_[k].split(":"), chr, start, true, false, null,true,
					"lrr",false,false);
			Data pl =  pca? null : new Data(plate_f[k], null, "plate", start, false,false,
					("PLATE".split(":")),true, "plate",false,false);
							
			
			
			//Data.period = thin;
			
			indivs[k] = new ArrayList<String>(lr.indiv);
			if(lrr==null){
				lrr = lr;
			}
			else{
				lrr.append(lr);
			}
			if(plate==null){
				plate = pl;
			}
			else plate.append(pl);
			
		}
		//System.err.println(indivs[1].size());
	}
	List<String>[] indivs;

	List<RealVector>[] correlatedVectors;
	List<String> indiv;

	public void run(double q, double thresh) throws Exception {
		CC d1 = new CC(new AbstractData[] { lrr, plate },
				new File("."), new Double[] {q,q}, new int[] {0,0});
		d1.centralise(true, true,true);
		d1.updateY();
		d1.decompose( false, thresh);
		System.err.println("h");
		correlatedVectors = d1.correlatedVectors;
		
		indiv = d1.indiv;
		if(plot){
		for(int k=0; k<d1.correlatedVectors.length; k++){
			d1.plotCC1("correlation_"+d1.data[k].type_name+".pdf",k,-1);
		}
		}
		normalise();
	}

	
	public void normalise(){
		for(int j=0; j<correlatedVectors.length; j++){
		for (int k = 0; k < correlatedVectors[j].size(); k++) {
			RealVector v1 = correlatedVectors[j].get(k);
			double norm = v1.getNorm();
			correlatedVectors[j].set(k,v1.mapDivide(norm / (double) indiv.size()));
		
		}
		}
	}
	public void print(File[] out, List<String> indiv_) throws Exception {
		
		for(int j1=0; j1<correlatedVectors.length; j1++){
		PrintWriter pw = new PrintWriter(
				new BufferedWriter(new FileWriter(out[j1])));
		// pw.close();
		for(int k=0; k<correlatedVectors[j1].size(); k++){
		//	double mean = sum(correlatedVectors.get(k));
			for(int j=0; j<k; j++){
				double dp = correlatedVectors[j1].get(k).dotProduct(correlatedVectors[j1].get(j));
				if(Math.abs(dp)>1e-5) {
					throw new RuntimeException("!! "+dp);
				}
			//	System.err.println(dp);
			}
		}
		for(int k=0; k<indiv_.size(); k++){
			pw.print(indiv_.get(k));
			int k1 = this.indiv.indexOf(indiv_.get(k));
			if(k1!=k) throw new RuntimeException("!!");
			for(int jk=0; jk<this.correlatedVectors[j1].size(); jk++){
				pw.print("\t"+
//						String.format("%7.5f", 
						correlatedVectors[j1].get(jk).getEntry(k1)
	//					).trim()
						);
			}
			pw.println();
		}
		pw.close();
		}
	}

	private double sum(RealVector realVector) {
		double sum=0;
		for(int k=0; k<realVector.getDimension(); k++){
			sum+=realVector.getEntry(k);
		}
		return sum;
	}


	private static Data getData(File[] files, String[] object, String chr,
			int[] start, boolean modify, boolean c, String[] asList,
			String type_name) throws Exception {
		boolean remZero = files.length == 1;
		Data d = null;
		for (int k = 0; k < files.length; k++) {
			try {
				Data d1 = new Data(files[k], object, chr, start, modify, c,
						asList, true,type_name, remZero,false);
				if (d == null)
					d = d1;
				else
					d.append(d1);
			} catch (Exception exc) {
				System.err.println(exc.getMessage());
			}
		}
		return d;
	}
}
