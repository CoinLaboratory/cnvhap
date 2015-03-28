package data.pca;

import java.io.File;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import lc1.util.ApacheCompressor;
import lc1.util.CompressDir;
import lc1.util.Constants;

import org.apache.commons.compress.archivers.zip.ZipFile;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import data.util.ApacheUtil;


public abstract class ProcessZip extends AbstractProcessZip {
	
	static int no_reps = 4;
	List<RealVector>[] cache = null;
	 final CompressDir[] output;
	 //	BufferedReader br;
	
	
		protected OutputStreamWriter[] snps;
	 double[][] res;
	 public static void main(String[] args){
	String[] n = args[10].split(":");
	regionsF = args[11];
	buildF = args[12];
	for(int k=0; k<n.length; k++){
	ProcessZip cv;
	no_reps = Integer.parseInt(args[7]);
	tocache  = Boolean.parseBoolean(args[8]);
	
	boolean mergeRegions = Boolean.parseBoolean(args[9]);
	Integer no_pcs = n[k].equals("null") ? null : Integer.parseInt(n[k]);
	if(args[6].equals("exact")){
		tocache = false;
	cv= new CalcPCAExact(new File("."), args[0].split(":"), args[1].replace('_',' ').split(":"), args[2].split(":"),
			args[3].split(":"),args[4].split(":"), args[5], mergeRegions,no_pcs);
	}
	else{
		cv= new CalcPCA(new File("."), args[0].split(":"), args[1].replace('_',' ').split(":"), args[2].split(":"),
				args[3].split(":"),args[4].split(":"), args[5], mergeRegions,no_pcs);
	}
	try{
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
	}catch(Exception exc){
		exc.printStackTrace();
	}
	}
}


RealMatrix[] pcs_to_adjust; //genomewide pcs to project out
GramSchmidt[] gs;





@Override
void makeRes() {
	 int totsze = assoc.Constants.sum(sze);
	res = new double[lrr_id[0].length][totsze];
}
@Override
void saveSampleInfo(List<String> sample_info, int k) {
	try{
		
	
	OutputStreamWriter ow = this.output[k].getWriter("Samples", true);
	for(int ii=0; ii<sample_info.size(); ii++){
		ow.write(sample_id[k].get(ii));
	}
	 output[k].closeWriter(ow);
	 ow = this.output[k].getWriter("Name", true);
	  if(ProcessZip.singleEntry){//assumes lrrst has length 1
		  for(int kk=0; kk<no_reps; kk++){
			  ow.write("PC"+kk);
			  ow.write(kk<no_reps-1 ? "\t" :"\n");
		  }
	  }else{
	  for(int kk=0; kk<lrrst.length; kk++){
		  ow.write(lrrst[kk]);
		  ow.write(kk<lrrst.length-1 ? "\t" :"\n");
	  }
	  }
	  ow.write("chr\tbuild36_start\tbuild36_end\tid\n");
	  ow.write("id\tvariance\n");
	  output[k].closeWriter(ow);
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

public ProcessZip(File pdir, String[] dir, String[] lrrst, String[] strtype, String[] chroms1_, 
		String[] thin, String thresh, boolean mergeRegions,  Integer no_pcs) {
	super(pdir, dir, lrrst, strtype, chroms1_, thin, thresh, mergeRegions, no_pcs);
	if(tocache){
		this.cache = new List[dir.length];
		for(int k=0; k<cache.length; k++) cache[k]= new ArrayList<RealVector>();
		
	}
	subdir_ind = 4;
	output = new CompressDir[dir.length];
	this.snps = new OutputStreamWriter[dir.length];
	for(int k=0; k<dir.length; k++){
		String nme__ = this.chroms.size()==1 ? "pcs_"+this.chroms.keySet().iterator().next()+"_out" : "pcs_out";
		chr_out = nme__+(no_pcs==null ? "" : no_pcs+"");
		File outf = new File(dir[k],chr_out);
		if(!overwrite){
		for(int i=0; outf.exists(); i++){
			chr_out = nme__+(no_pcs==null ? "" : no_pcs+"")+"_"+i;
			outf = new File(dir[k],chr_out);
		}
		}
		try{
		output[k] = new CompressDir(outf);
		output[k].delete = true;
		this.snps[k] = output[k].getWriter("SNPS", false);
		}catch(Exception exc){ exc.printStackTrace();}
	}
	boolean pcFilesExist = true;
	ZipFile[] pcs_in = new ZipFile[dir.length];
	File[] pcs_inf = new File[pcs_in.length];
	for(int k=0; k<pcs_in.length; k++){
		File pcs_ink = new File(dir[k],"pcs_in.zip");
		if(pcs_ink.exists()) pcs_inf[k] = pcs_ink;
		else pcFilesExist = false;
	}
	if(pcFilesExist) pcs_in = this.getZip("pcs_in",pcs_inf);
	this.pcs_to_adjust = new RealMatrix[lrrst.length];
	this.gs = new GramSchmidt[lrrst.length];
	if(pcFilesExist && no_pcs!=null && no_pcs>0){
		
			//
		List<String> SNPS = ApacheUtil.read(pcs_in[0], "SNPS", 3);
		if(no_pcs < SNPS.size()){
			SNPS = SNPS.subList(0, no_pcs);
		}
		for(int kk=0; kk<pcs_to_adjust.length; kk++){
			pcs_to_adjust[kk] = new Array2DRowRealMatrix(Constants.sum(sze), SNPS.size());
		}
		for(int k=0 ;k<pcs_in.length; k++){
			List<String> nme = Arrays.asList(ApacheUtil.read(pcs_in[k], "Name").get(0).split("\t"));
			List<String> samples = ApacheUtil.read(pcs_in[k], "Samples", 0);
			RealMatrix mat = new Array2DRowRealMatrix(samples.size(), lrrst.length);
			for(int i=0; i<SNPS.size(); i++){
				ApacheUtil.read(pcs_in[k], SNPS.get(i),mat);
				for(int j=0; j<lrrst.length; j++){
					int j1 = nme.indexOf(lrrst[j]);
					if(j1<0 && nme.size()==1) j1 =0;
					for(int k1=0; k1<alias[k].length; k1++){
						this.pcs_to_adjust[j].setEntry(k1,i,mat.getEntry(alias[k][k1],j1));
					}
				}
			}
		}
	
		for(int k=0; k<gs.length; k++){
			gs[k]= new GramSchmidt(pcs_to_adjust[k]);
		}
	}
}

@Override
public void reInitialise(String nme1) throws Exception{
	super.reInitialise(nme1);
	if(checkOrder){
		for(int jk=0; jk < dir.length; jk++){
			  List<String> indiv = ApacheCompressor.readZipFrom(zf[jk], "Samples", include[jk]);
				if(!this.sample_id[jk].equals(indiv)) throw new RuntimeException("!!");
		}
	
	
	}
}


@Override
boolean readData(String id, boolean complete, int k) {
complete = complete & ApacheCompressor.readZip(zf[k], id, res, lrr_id[k], offset[k], this.type, threshNA, this.include[k]);
if(complete && depth_index>=0){
	int no_cols = lrr_id[k].length;
	for(int j=0; j<no_cols; j++){
		for(int i=0; i<depth.size(); i++){
			res[j][i] = res[j][i]/this.depth.get(i);
		}
	}
}
return complete;	
}
/*@Override
public void process(){
for(int k=0; k<l1.length; k++){
	// System.err.println("running "+chr+" "+pos);
	  process(k);
	  
  }
}*/
void process(int k) {
	 l1[k].process(res[k]);
	  if(tocache && ((UpdateProjection)l1[k]).var>0){
		 // this.snp_list.add(chr+"\t"+pos+"\t"+id);
		  this.cache[k].add(l1[k].d.copy());
	  }
}
public void runInnerQuick() {
	
	for(int k=0; k<cache.length; k++){
		for(int j=0; j<cache[k].size(); j++){
				(l1[k]).process(cache[k].get(j));
			cache[k].set(j, l1[k].d);
		}
	}
}
@Override
protected void finish(){
	super.finish();
	try{
		for(int i=0; i<this.output.length; i++){
			this.snps[i].close();
			output[i].run();
			output[i].close();
		}
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

 
 
 
 /* returns start/end location */
 

 
 
} // Example