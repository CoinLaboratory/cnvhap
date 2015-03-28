package data.pca;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;

public class CalcPCAExact extends ProcessZip {

static boolean check = false;


//int max_cnt = 100;//1*1000*1000;//Integer.MAX_VALUE;
int no_reps = 2;

	
public void makeHists(String nme, String[] lrrst){
	int len = lrrst.length;
	 this.l1 = new ProcessSNP[len];
	 for(int k=0; k<len; k++){
		 this.l1[k] =new PCAExact(dir, nme+".hist.txt"+lrrst[k], this.type[k], this.gs[k]);
	 }
 }

	public CalcPCAExact(File dir, String[] nmes, String[] lrrst,String[] strtype, String[] chroms, String[] thin, 
			String thresh, boolean merge, Integer no_pcs) {
		super(dir, nmes, lrrst, strtype, chroms, thin, thresh, merge, no_pcs);
		
		// TODO Auto-generated constructor stub
	}
	
	
	
	
	//List<RealVector> pcs_auto = new ArrayList<RealVector>(0);
	//List<RealVector> pcs_sex = new ArrayList<RealVector>(0);
	//RealVector p_auto = null;
	
	
	public void run() throws Exception{
			this.run1();
			
			for(int j=0; j<l1.length; j++){
				if(AbstractProcessZip.singleEntry){
					l1[j].print(-1,aliasRev, this.output, this.currentRegion.toString()+"PC");
					this.snps[j].write(currentRegion.line()+"PC"+"\n");
			}else{
				for(int k=0; k<Math.min(ProcessZip.no_reps, ((PCAExact)l1[j]).U.getColumnDimension()); k++){
					l1[j].print(k,aliasRev, this.output, this.currentRegion.toString()+"PC"+k);
					this.snps[j].write(currentRegion.line()+"PC"+k+"\n");
				}
			}
			}
	}
	
	
	public void run1() throws Exception{
	
		this.runInner(0, 0);
	
	 for(int k=0; k<l1.length; k++){
		 double diff =  l1[k].getHistogram( );
	 }
	}

	

	private File[] getFile(File[] dir, String string) {
		// TODO Auto-generated method stub
		return null;
	}

	private void close(ZipFile[] zf) throws IOException {
		for(int k=0; k<zf.length; k++)zf[k].close();
		
	}
}
