package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;

public class StringAsZipLike implements ZipFileAccess {
 public StringAsZipLike(File in, String first, String last, int locind, int chrind) throws Exception {
		this.f =in.getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in)))): new BufferedReader(new FileReader(in)) ;
		String st = "";
		while((st = f.readLine()).startsWith("#")){
			currentline = st;
		}
		this.curr = currentline.split("\t");
		this.locind = locind;
		this.chrind  =chrind;
		this.formatind = Arrays.asList(curr).indexOf("FORMAT");
		//locind = Arrays.asList(curr).indexOf("START");
		//chrind = Arrays.asList(curr).indexOf("#CHROM");
		List<String> currL = Arrays.asList(curr);
		start = currL.indexOf("FORMAT")>=0 ? currL.indexOf("FORMAT")+1 : currL.indexOf("END") + 1;
		end = curr.length;
		this.sum = new double[end - start];
		//nextLine();
		currentline = st;
		this.curr = currentline.split("\t");
		while(curr!=null && !(curr[chrind]+"_"+curr[locind]).equals(first)){
			nextLine();
		}
		this.last = last;
	}
final int locind, chrind, formatind;
final String last;
final int start, end;
BufferedReader f;
String currentline;
String[] curr;

double[] sum;
boolean trysum = true;
void nextLine() throws Exception{
	currentline = f.readLine();
	if(currentline!=null){
	this.curr = currentline.split("\t");
	if(Constants.allowChrom(curr[chrind])>=0  && trysum){ ///////////NOTE NEED TO CHANGE ThIS BACK
		try{
		for(int k=0; k<sum.length; k++){
			sum[k] += Double.parseDouble(curr[k+start]);
		}
		}catch(Exception exc){
			trysum = false;
		}
	}
	}else{
		curr = null;
	}
	
}

/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer)
 */
@Override
public  List<String> getIndiv(String entryName, Integer column) throws Exception{
    return getIndiv(entryName, column, Constants.splString());
}

public void skip(String entryName) throws Exception{
	while(!entryName.endsWith("_"+curr[locind]) || !entryName.startsWith(curr[chrind]+"_")) {
		nextLine();
	}
}

/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer, java.lang.String)
 */
@Override
public  List<String> getIndiv( String entryName, Integer column, String spl) throws Exception{
	
	skip(entryName);
	List<String> res = Arrays.asList(curr).subList(start, end);
	for(int i=0; i<res.size(); i++){
		res.set(i,res.get(i).replace(':', '\t'));
	}
		
		
	   return res;
}
/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer, java.lang.String[])
 */
@Override
public  boolean getIndiv(String entryName, Integer column, String[] indiv) throws Exception{
	
	skip(entryName);
	System.arraycopy(curr, start, indiv, 0, indiv.length);
	nextLine();
  return true;
//   return indiv;
}
/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getBufferedReader(java.lang.String)
 */
@Override
public  BufferedReader getBufferedReader(String string) throws Exception{
 throw new RuntimeException("not implemented");
     
}
private void read( String string, List<String> indiv,
       int k) throws Exception {
	skip(string);
	
	
	for(int i=0; i<indiv.size(); i++){
		indiv.set(i, curr[i+start]);
	}
	nextLine();
   
}
@Override
public void getAvgDepth(String pref, int avgDepthCol, List<Integer> dToInc,
		File samplesFile, List<Integer> ploidy, List header_sample, List avgDepth) {
	try{
		while((currentline )!=null){
			nextLine();
		}
		f.close();
	}catch(Exception exc){
		exc.printStackTrace();
	}
	for(int k_=0; k_<dToInc.size(); k_++){
		 int k = dToInc.get(k_);
  		avgDepth.add(sum[k]);//.split("\\s+");//split(":");
	 }
	//	String[] str1 = ad.get(k).spl
	
	
}


}
