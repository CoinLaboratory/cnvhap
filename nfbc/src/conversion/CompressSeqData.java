package conversion;

import java.io.File;

public class CompressSeqData extends CompressCNVResults {

    public int getStart(String[]str){
	
	  return Integer.parseInt(str[5]);
  }
    
   
  public int getEnd(String[]st1){
	  return Integer.parseInt(st1[6]);
  }
  public int getCN(String[] str){
	  return 	str[1].equals("deletion") ? 1:3;
  }
  public String getId(String[] str){
	return  "yanhuang";
  }
  public String getChrom(String[] str){
		return str[4];
	}
 
 public   CompressSeqData(File f, String chr, File samples) throws Exception{
       super(f, chr, false, samples, null);
    }
   // public CompressCNVResults(File f, String chr, boolean header, File rawDir) throws Exception {
    	
}
