package conversion;

import java.io.File;

public class CompressedAgilentCNV extends CompressCNVResults {

	CompressedAgilentCNV(File f, String chr, File rawDir) throws Exception {
		super(f, chr, true, rawDir, null);
	}

	@Override
	public int getCN(String[] str) {
		String chrom  = getChrom(str);
		double d = Double.parseDouble(str[22]);
		if(chrom.indexOf('X')>=0){
			if(d<-0.5) return 0;
			else return 2;
		}
		else{
		
		if(d<-1.0) return 0;
		else if(d<0) return 1;
		else if (d>0.6) return 4;
		else return 3;
		}
	}

	@Override
	public int getEnd(String[] st1) {
		return Integer.parseInt(st1[11]);
	}

	@Override
	public String getId(String[] str) {
		return str[6];
	}

	@Override
	public int getStart(String[] str) {
		return Integer.parseInt(str[10]);
	}
	@Override
	public String getChrom(String[] str){
		if(str.length<=1) return "";
		return "chr"+str[9];
	}

}
