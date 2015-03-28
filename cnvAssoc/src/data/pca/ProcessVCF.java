package data.pca;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import lc1.sfs.AlleleSFS;
import lc1.sfs.GenoSFS;
import lc1.sfs.SFS;
import lc1.util.ApacheCompressor;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

public class ProcessVCF extends AbstractProcessZip {
	
	static int no_reps;
	
	 public static void main(String[] args){
		 try{
			 
		 
			String[] n = args[10].split(":");
			regionsF = args[11];
			buildF = args[12];
			for(int k=0; k<n.length; k++){
			ProcessVCF cv;
			no_reps = Integer.parseInt(args[7]);
			tocache  = Boolean.parseBoolean(args[8]);
			functString = args[6];
			String[] limits = args[4].split(":");
			limitDone = limits.length<=3 ? Integer.MAX_VALUE-1 : Integer.parseInt(limits[3]);
			boolean mergeRegions = Boolean.parseBoolean(args[9]);
			Integer no_pcs = n[k].equals("null") ? null : Integer.parseInt(n[k]);
			
				cv= new ProcessVCF(new File("."), args[0].split(":"), args[1].replace('_',' ').split(":"), args[2].split(":"),
						args[3].split(":"),args[4].split(":"), args[5], mergeRegions,no_pcs);
			
			if(cv.regions==null){	
				//cv.keepZipOpen = false;
				cv.currentRegion.chrom_out = cv.chr_out;
				cv.run();
			  
			
			}
			else{
			//	cv.keepZipOpen = true;
			  for(Iterator<Location> it = cv.regions.iterator();it.hasNext();){
				  Location loc = it.next();
				  cv.currentRegion = loc;
				  cv.currentRegion.chrom_out = cv.chr_out;
				  loc.chrom_out = cv.chr_out;
				  cv.run();
			  }
			}
			  cv.finish();
			
			}
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
		}
	public ProcessVCF(File pdir, String[] dir, String[] lrrst,
			String[] strtype, String[] chroms1_, String[] thin, String thresh,
			boolean mergeRegions, Integer no_pcs) {
		super(pdir, dir, lrrst, strtype, chroms1_, thin, thresh, mergeRegions, no_pcs);
		this.toinclude = new boolean[dir.length][];
		for(int k=0; k<toinclude.length; k++){
			toinclude[k] =new boolean[sze[k]];
		}
		
		// TODO Auto-generated constructor stub
	}

	static String functString="FreqFunct";
	@Override
	public void makeHists(String nme, String[] lrrst) {
		// TODO Auto-generated method stub
		
		try{
			File resultsDir = new File("results");
			resultsDir.mkdir();
			FileOutputStream fos = new FileOutputStream(new File(resultsDir, this.chroms.keySet().iterator().next()+"_results.txt"+
					(compress ? ".gz" :"")));
			this.respw = new PrintWriter(new BufferedWriter(
					new OutputStreamWriter(
							compress ? new GZIPOutputStream(fos) : fos)));
			
			
			//fst1 = new SFS[1];
			List<SFS> fst = new ArrayList<SFS>();
			List<int[]> inds = new ArrayList<int[]>();
			List<String>nme1 = new ArrayList<String>();
		getSFS(functString, nme1, sze, fst, inds, chroms.firstKey()+"-"+chroms.lastKey());
		this.fst1 = fst.toArray(new SFS[0]);
		this.inds = inds.toArray(new int[0][]);
		
		this.loglh1  = new double[fst1.length];
		
		this.respw.print("snpid\tpos");
		for(int k=0;k<nme1.size(); k++){
			respw.print("\t"+nme1.get(k));
		}
		respw.println();
		}catch(Exception exc){
			exc.printStackTrace();
		}	
		
	}

	public static void getSFS(String functString2, List<String> nme1, int[] sze, List<SFS> fst, List<int[]> inds, String name) throws Exception{
		String[] functs = functString2.split(":");
		List<String> singleFunctsAllele = new ArrayList<String>(); //applied to each cohort sep
		List<String> jointFunctsAllele = new ArrayList<String>();
		List<String> singleFunctsGeno = new ArrayList<String>(); //applied to each cohort sep
		List<String> jointFunctsGeno = new ArrayList<String>();
		List<String>[] all = new List[] {singleFunctsAllele, jointFunctsAllele, singleFunctsGeno, jointFunctsGeno};
		boolean[] joint = new boolean[] {false, true, false, true};
		boolean[] allele = new boolean[] {true, true, false, false};
		for(int i=0; i<functs.length; i++){
			String[] str1 = functs[i].split("_");
			try{
			Class clazz = Class.forName("lc1.sfs.freqfunct."+str1[0]);
			if((Short)clazz.getField("numArgs").get(null)>1) jointFunctsAllele.add(functs[i]);
			else singleFunctsAllele.add(functs[i]);
			}catch(ClassNotFoundException exc){
				Class clazz = Class.forName("lc1.sfs.genofunct."+str1[0]);
				if((Short)clazz.getField("numArgs").get(null)>1) jointFunctsGeno.add(functs[i]);
				else singleFunctsGeno.add(functs[i]);
			}
		}
		
		
		getSFS(singleFunctsAllele,sze, true, AlleleSFS.class, nme1 ,fst, inds,name);
		getSFS(jointFunctsAllele,sze, false, AlleleSFS.class, nme1 ,fst, inds,name);
		getSFS(singleFunctsGeno,sze, true, GenoSFS.class, nme1 ,fst, inds,name);
		getSFS(jointFunctsGeno,sze, false, GenoSFS.class, nme1 ,fst, inds,name);
		//"Fstfunct:AlleleDifffunct:FixedDiffFunct_0.99;0.9;0.8;0.5:FstFunctThresh_0.99;0.9;0.9;0.8;0.5 ");
		
	/*if(sze.length==2){
		fst1[0]  = new AlleleSFS(sze,2, "Fstfunct:AlleleDifffunct:FixedDiffFunct_0.99;0.98;0.97;0.96;0.95:FstFunctThresh_0.99;0.98;0.97;0.96;0.95",nme1);
	}else{
		fst1[0]  = new AlleleSFS(sze,2,"FreqFunct",nme1);
		
	}*/
	
		
		// TODO Auto-generated method stub
	}

	private static void getSFS (
			List<String> functs, int[] sze, boolean sep,
			Class class1, List<String> nme1, List<SFS> res, List<int[]> inds, String name) throws Exception{
		if(functs.size()>0){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<functs.size();i++){
				sb.append(functs.get(i));
				sb.append(i<functs.size()-1  ? ":" : "");
			}
			String str = sb.toString();
		Constructor c = class1.getConstructor(new Class[] {sze.getClass(), int.class, String.class, List.class, String.class});
		if(sep){
			for(int i=0; i<sze.length; i++){
				inds.add(new int[] {i});
				res.add((SFS)c.newInstance(new Object[] {new int[] {sze[i]},2,str,nme1,name}));
			}
		}else{
			int[] inds_ = new int[sze.length];
			for(int i=0; i<sze.length; i++) inds_[i] = i;
			inds.add(inds_);
			res.add((SFS)c.newInstance(new Object[] {sze,2,str,nme1,name}));
		}
		}
	}

	RealMatrix[]res;  //pop_index x individual_index x genotype_index 
	
	@Override
	void makeRes() {
		res = new RealMatrix[sze.length];
		for(int k=0; k<sze.length; k++){
			res[k] = new Array2DRowRealMatrix(sze[k],lrr_id[k].length);
		}
	}
	
	SFS[] fst1; //does inbreeding coeffs
	int[][] inds; //which data corresponds to which fst
	double[] loglh1;
	 PrintWriter respw, respw1;
	 
	boolean calculatePerSite = false;// whether to calculate per site
	 
	
	public void process(int k){
		
	}
	
	@Override
	public void process(String id, int pos) {
		// TODO Auto-generated method stub
		if(calculatePerSite) respw.print(id+"\t"+pos);
		for(int kk=0; kk<fst1.length; kk++){
			if(check(readFile, inds[kk])){
			fst1[kk].setData(res,this.toinclude, inds[kk]);
				loglh1[kk]+=fst1[kk].expectation();
				if(calculatePerSite){
					RealMatrix resu = fst1[kk].calculate();
					for(int k=0; k<resu.getRowDimension(); k++){
						respw.print("\t");
						for(int k1=0; k1<resu.getColumnDimension(); k1++){
							if(k1>0)respw.print(";");
							respw.print(String.format("%6.4g",resu.getEntry(k,k1)).trim());
						}
					}
				}
			}else{
				if(calculatePerSite){
					for(int k=0; k<fst1[kk].getDimension(); k++){
						respw.print("\tNC");
					}
				}
			}
		}
		respw.println();
		respw.flush();
	}
	
	

	private boolean check(boolean[] readFile, int[] is) {
	for(int k=0; k<is.length; k++){
		if(!readFile[is[k]]) {
			return false;
		}
	}
	return true;
	}

	final boolean[][] toinclude;
	String alleles ="";
	String rsid = "";
	@Override
	boolean readData(String id, boolean complete, int k) {
		Arrays.fill(toinclude[k], true);
		if(k==0) {
			alleles ="" ;
			rsid = "";
		}
		try{
			String[] comms 

			
		//	 = ApacheCompressor.readZip(zf[k], ze[k],id, res[k],toinclude[k], lrr_id[k], alleles, rsid);

			 = ApacheCompressor.readZip(zf[k], ze[k],id, res[k],toinclude[k], lrr_id[k], alleles, rsid, gl[k][0]);
//>>>>>>> .r258
			if(comms!=null && comms.length>2){alleles = comms[0];
			rsid = comms[1];
			}
			complete = complete & alleles!=null;
		}catch(Exception exc){
			exc.printStackTrace();
		}
	return complete;	
	}

	
	@Override
	protected void finish(){
		System.err.println("finishing");
		super.finish();
		respw.close();
		//this.missingpw.close();
	}
	
	
	@Override
	public void run() throws Exception{
		for(int k=0; k<this.no_reps; k++){
			Arrays.fill(this.loglh1, 0);
			if(k==no_reps-1) this.calculatePerSite = true;
			this.runInner(k, 0);
			for(int j=0; j<fst1.length; j++){
				if(fst1[j]!=null){
					System.err.println("log lh_"+fst1[j].name+" "+this.loglh1[j]);
					fst1[j].transfer();
				}
			}
			print(k+".txt");
		
		}
		print("final.txt");
		}
	
	public void print(String outnme){
		for(int kk=0; kk<fst1.length; kk++){
			fst1[kk].print(outnme);
		}
	}
	@Override
	void saveSampleInfo(List<String> sample_info, int k) {
		// TODO Auto-generated method stub

	}

}
