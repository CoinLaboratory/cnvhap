package lc1.dp.data.collection;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.print.attribute.standard.MediaSize.NA;

import org.apache.commons.compress.archivers.zip.ZipFile;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

//import data.pca.GramSchmidt;
import data.util.ApacheUtil;

public class DataProjection {

	private RealMatrix pcs_to_adjust;
	private GramSchmidt gs = null;
	RealVector todivid;
	Set<String> todrop = new HashSet<String>();
	public DataProjection(File pcsInFile, int no_pcs, List<String> l1, List<String> l2, boolean isdepth) throws Exception{
		this.isdepth = isdepth;
		ZipFile pcs_in = new ZipFile(pcsInFile);
		List<String> SNPS = ApacheUtil.read(pcs_in, "SNPS", 3);
		List<String> samples = ApacheUtil.read(pcs_in, "Samples", 0);
		todivid = new ArrayRealVector(l2.size());
		for(int i=0; i<l2.size(); i++){
			todivid.setEntry(i, Double.parseDouble(l2.get(i)));
		}
		if(no_pcs < SNPS.size()){
			SNPS = SNPS.subList(0, no_pcs);
		}
		
		if(SNPS.size()>0){
			pcs_to_adjust = new Array2DRowRealMatrix(l1.size(), SNPS.size());
		}
		Set<Integer> NA = new HashSet<Integer>();
			List<String> nme = Arrays.asList(ApacheUtil.read(pcs_in, "Name").get(0).split("\t"));
			RealMatrix mat = new Array2DRowRealMatrix(samples.size(), 1);
			for(int i=0; i<SNPS.size(); i++){
				ApacheUtil.read(pcs_in, SNPS.get(i),mat);
					int j1 = nme.indexOf("PC0");
					if(j1<0 && nme.size()==1) j1 =0;
					//for(int k1=0; k1<alias[k].length; k1++){
					for(int k=0; k<l1.size(); k++){
						String v1 =l1.get(k);
						int k1 = samples.indexOf(v1);
						double entry = k1<0 ? Double.NaN : mat.getEntry(k1,0);
						if(k1<0){
							
							//throw new RuntimeException(" problem with "+v1);
							System.err.println("could not find "+v1);
						}
						this.pcs_to_adjust.setEntry(k,i,entry);
						if(Double.isNaN(entry)){
							NA.add(k);
							this.todrop.add(l1.get(k));
						}
					}
					
				//	}
			}
			if(SNPS.size()>0){
			int noRows = pcs_to_adjust.getRowDimension();
			int[] nonNaIndex = new int[noRows - NA.size()];
			int[] naIndex = new int[NA.size()];
			int k1=0;
			int k2=0;
			for(int k=0; k<noRows; k++){
				if(!NA.contains(k)){
					nonNaIndex[k1] = k;
					k1++;
				}else{
					naIndex[k2] =k;
					k2++;
				}
			}
			 gs= new GramSchmidt(pcs_to_adjust, nonNaIndex, naIndex);
			}
			d  = new ArrayRealVector(l1.size());
			pcs_in.close();
		
	}

	
	RealVector d;
	final boolean isdepth;
	public void convert(List<String> l, int lrr_index) {
		
		int len = l.size();
		for(int i=0; i<len; i++){
			try{
				String[] str = l.get(i).split("\t");
				double v = Double.parseDouble(str[lrr_index]);
				
				/*if(Double.isNaN(v)){
					
					System.err.println("h");
				}*/
				
			d.setEntry(i, isdepth ? v /  this.todivid.getEntry(i) : v - this.todivid.getEntry(i));
			}catch(Exception exc){
				exc.printStackTrace();
				System.exit(0);
			}
		}
		RealVector d1  = gs!=null ? gs.removeProj(d) : d;
		for(int i=0; i<len; i++){
			String[] str = l.get(i).split("\t");
			double v = d1.getEntry(i) ;
			str[lrr_index] = ( isdepth ? v *  this.todivid.getEntry(i) : v + this.todivid.getEntry(i))+"";
			l.set(i, getStr(str,"\t"));
		}
		
	}
	private RealVector standardise(RealVector d2) {
		// TODO Auto-generated method stub
		return null;
	}
	private String getStr(String[] str, String string) {
		StringBuffer sb = new StringBuffer(str[0]);
		for(int i=1; i<str.length; i++){
			sb.append("\t"+str[i]);
		}
		return sb.toString();
	}
	
}
