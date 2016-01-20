package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;

public class StringAsZipLike implements ZipFileAccess {
 public StringAsZipLike(File in, String first, String last) throws Exception {
		this.f =new BufferedReader(new FileReader(in));
		currentline = f.readLine();
		this.curr = currentline.split("\t");
		snpind = Arrays.asList(curr).indexOf("ID");
		chrind = Arrays.asList(curr).indexOf("#CHROM");
		start = Arrays.asList(curr).indexOf("FORMAT")+1;
		end = curr.length;
		this.sum = new double[end - start];
		while(!(curr[chrind]+"_"+curr[snpind]).equals(first)){
			nextLine();
		}
		this.last = last;
	}
final int snpind, chrind;
final String last;
final int start, end;
BufferedReader f;
String currentline;
String[] curr;

double[] sum;

void nextLine() throws Exception{
	currentline = f.readLine();
	if(currentline!=null){
	this.curr = currentline.split("\t");
	if(Constants.allowChrom(curr[chrind])>=0){
	for(int k=0; k<sum.length; k++){
		sum[k] += Double.parseDouble(curr[k+start]);
	}
	}
	}
	
}

/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer)
 */
@Override
public  List<String> getIndiv(String entryName, Integer column) throws Exception{
    return getIndiv(entryName, column, Constants.splString());
}



/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer, java.lang.String)
 */
@Override
public  List<String> getIndiv( String entryName, Integer column, String spl) throws Exception{
	
	while(!entryName.endsWith(curr[snpind]) || !entryName.startsWith(curr[chrind])) {
		nextLine();
	}
   return Arrays.asList(curr).subList(start, end);
}
/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer, java.lang.String[])
 */
@Override
public  boolean getIndiv(String entryName, Integer column, String[] indiv) throws Exception{
	
	while(!entryName.equals(curr[snpind])) {
		System.err.println("skipping "+curr[snpind]);
		nextLine();
	}
	
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
	while(!string.equals(curr[snpind])) {
		System.err.println("skipping "+curr[snpind]);
		nextLine();
	}
	
	
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
	for(int i=0; i<this.sum.length; i++){
		avgDepth.add(sum[i]);
	}
	
}


}
