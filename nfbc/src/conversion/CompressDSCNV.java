package conversion;

import java.io.File;

public class CompressDSCNV extends CompressCNVResults {

    public int getStart(String[]str){
	
	  return Integer.parseInt(str[2]);
  }
  public int getEnd(String[]st1){
	  return Integer.parseInt(st1[3]);
  }
  public int getCN(String[] str){
	  double d =  Double.parseDouble(str[5]) ;
	  if(str[6].indexOf("DEL")>=0){

			if(d<-1.0) return 0;
			else return 1;
	  }
	  else   if(str[6].indexOf("CNV")>=0){
		   if (d>0.6) return 4;
			else return 3;
	  }
	  else  throw new RuntimeException("!!");
		
  }
  public String getId(String[] str){
	return   str[0].split("#")[0];
  }
  public String getChrom(String[] str){
		return "chr"+str[1];
	}
    CompressDSCNV(File f, String chr, File samples) throws Exception{
       super(f, chr, true, samples, null);
    }
   // public CompressCNVResults(File f, String chr, boolean header, File rawDir) throws Exception {
    	
}
