package conversion;

import java.io.File;

public class CompressPennCNV extends CompressCNVResults {

    public int getStart(String[]str){
	  String[] st1 = str[0].split(":");
	  return Integer.parseInt(st1[1]);
  }
  public int getEnd(String[]st1){
	  return Integer.parseInt(st1[1]);
  }
  public int getCN(String[] str){
	  return 	Integer.parseInt(str[5].split("=")[1]);
  }
  public String getId(String[] str){
	return   str[6].split("\\.")[1];
  }
  public String getChrom(String[] str){
		return str[0].split(":")[0];
	}
    CompressPennCNV(File f, String chr, File samples) throws Exception{
       super(f, chr, false, samples, null);
    }
   // public CompressCNVResults(File f, String chr, boolean header, File rawDir) throws Exception {
    	
}
