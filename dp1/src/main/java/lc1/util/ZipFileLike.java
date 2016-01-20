package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipFile;

public class ZipFileLike implements ZipFileAccess {
 public ZipFileLike(org.apache.commons.compress.archivers.zip.ZipFile zf2) {
		this.f = zf2;
	}

ZipFile f;


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
   BufferedReader  nxt = 
       new BufferedReader(new InputStreamReader(
       f.getInputStream(f.getEntry(entryName))));
   List<String> indiv = new ArrayList<String>();
   String st = "";
   while((st = nxt.readLine())!=null){
       String str = st.trim();
       if(column==null){
           indiv.add(str);
       }
       else if (column==-1){
    	   indiv.add("1");
       }
       else{
    	   String[] st_ = str.split(spl);
    	   if(st_.length>column){
    		   String res = st_[column];
    		   if(res.length()>0)
    			   indiv.add(res);
    	   }
       }
   }
  nxt.close();
   return indiv;
}
/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getIndiv(java.lang.String, java.lang.Integer, java.lang.String[])
 */
@Override
public  boolean getIndiv(String entryName, Integer column, String[] indiv) throws Exception{
   ZipArchiveEntry ent = f.getEntry(entryName);
   if(ent==null) {
	   return false;
   }
   BufferedReader  nxt = 
       new BufferedReader(new InputStreamReader(
       f.getInputStream(ent)));
 //  List<String> indiv = new ArrayList<String>();
   String st = "";
   int i=0;
  for(; (st = nxt.readLine())!=null; i++){
       String str = st.trim();
       if(column==null){
           indiv[i] = str;
       }
       else{
           indiv[i] = (str.split("\t")[column]);
       }
   }
  
  nxt.close();
  return i==indiv.length;
//   return indiv;
}
/* (non-Javadoc)
 * @see lc1.util.ZipFileAccess#getBufferedReader(java.lang.String)
 */
@Override
public  BufferedReader getBufferedReader(String string) throws Exception{
  ZipArchiveEntry entry =  f.getEntry(string);
  if(entry==null) return null;
   return   new BufferedReader(new InputStreamReader(
           f.getInputStream(entry)));
     
}
private void read( String string, List<String> indiv,
       int i) throws Exception {
   BufferedReader  nxt = 
       new BufferedReader(new InputStreamReader(
       f.getInputStream(f.getEntry(string))));
   String st = "";
   while((st = nxt.readLine())!=null){
       indiv.add(st.split("\t")[i]);
   }
   nxt.close();
   
}



@Override
public void getAvgDepth(String pref, int avgDepthCol, List<Integer> dToInc,
		File samplesFile,List<Integer> ploidy,List header_sample,  List avgDepth) {
	List<String> ad ;
	try{
	//File samplesFile = new File(f.getParentFile(), "Samples");
	if(samplesFile.exists()){
		ad  = ApacheCompressor.getIndiv(samplesFile);//,pref+"Samples");
	}
	
	else{
		ad = ApacheCompressor.getIndiv(f,pref+"Samples");
	}
	 for(int k_=0; k_<dToInc.size(); k_++){
		 int k = dToInc.get(k_);
   		String[] str1 = ad.get(k).split("\\s+");//split(":");
   		int avgDepthCol1 = avgDepthCol;
   		if(avgDepthCol1>=0){
   			avgDepth.add(Double.parseDouble(str1[avgDepthCol1]));
   		}
   		else if(avgDepthCol1<0 && str1.length>header_sample.size()){ avgDepthCol1 = header_sample.size();
   		 	avgDepth.add(Double.parseDouble(str1[avgDepthCol1]));
   		}else{
   			avgDepth.add(ploidy.get(k).doubleValue());
   		}
   		
	 }
	}catch(Exception exc){
		exc.printStackTrace();
	}
}



}
