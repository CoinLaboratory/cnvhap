package conversion;

import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class CompressCooperSeq extends CompressCNVResults {

    public int getStart(String[]st1){
    	 if(!Double.isNaN(this.agiStartEnd[0])){
   		  return (int) agiStartEnd[0];
   	  }
   	  if(!Double.isNaN(this.agiStartEnd[2])){
   		  return (int) agiStartEnd[2];
   	  }
	  return Integer.parseInt(st1[2]);
  }
  public int getEnd(String[]st1){
	  if(!Double.isNaN(this.agiStartEnd[0])){
		  return (int) agiStartEnd[1];
	  }
	  if(!Double.isNaN(this.agiStartEnd[2])){
		  return (int) agiStartEnd[3];
	  }
	  return Integer.parseInt(st1[3]);
  }
  public int getCN(String[] str){
	  for(int i=0; i<del_gain.length; i++){
		  if(del_gain[i]) {
			  return (i==0 || i==2) ? 1 : 3;
		  }
	  }
	  return 	(str[6].equals("S"))  ?3:  1;
  }
  public String getId(String[] str){
	String id =    m.get(str[0]);
	return id;
  }
  
  public Integer val_index ;
  public Integer agi_start_ind ;
  public Integer agi_end_ind ;
  public Integer nim_start_ind ;
  public Integer nim_end_ind ;
  public Integer nim_del ;
  public Integer nim_gain ;
  public Integer agi_del ;
  public Integer agi_gain;
  static String[] prop = "Nim Deletion	Nim_Gain	Agi Deletion	Agi Gain	Overlap With Validated Locus	NimRefinedStart	NimRefinedEnd	AgiRefinedStart	AgiRefinedEnd".split("\t");
  static String[] propField = "nim_del	nim_gain	agi_del	agi_gain	val_index	agi_start_ind	agi_end_ind	nim_start_ind	nim_end_ind".split("\t");
  
  
  
  @Override
  public boolean exclude(String[] str){
	// boolean excl =  !str[val_index].equals("yes");
	 boolean excl1 = setStartEnd(str);
	
	 return  excl1;
	 //return false; 
	 //return true;
  }
  double[] agiStartEnd = new double[4];
  boolean[] del_gain = new boolean[4];
  public boolean  setStartEnd(String[] str){
	 agiStartEnd[0] = getVal(str,this.agi_start_ind);
	 agiStartEnd[1] = getVal(str,this.agi_end_ind);
	 agiStartEnd[2] = getVal(str,this.nim_start_ind);
	 agiStartEnd[3] = getVal(str,this.nim_end_ind);
	 del_gain[0] = getBoolVal(str,this.agi_del);
	 del_gain[1] = getBoolVal(str,this.agi_gain);
	 del_gain[2] = getBoolVal(str,this.nim_del);
	 del_gain[3] = getBoolVal(str,this.nim_gain);
	 return Double.isNaN(agiStartEnd[0]) && Double.isNaN(agiStartEnd[2]);
  }
  private boolean getBoolVal(String[] str, int nim_gain2) {
	  if(nim_gain2>=str.length) return false;
	  else return str[nim_gain2 ].indexOf("1")>=0;
	
}
private double getVal(String[] string1, int ind) {
	  if(string1.length<=ind) return Double.NaN;
	  String string = string1[ind];
	  if(string.length()>0)
	 return Double.parseDouble(string);
	  else return Double.NaN;
}
@Override
  void setHeader(String[] str){
  	List<String> l = Arrays.asList(str);
  	for(int i=0; i<prop.length; i++){
  		try{
  			String fie = (propField[i]);
  			Class clazz = this.getClass();
  			Field fi = clazz.getField(fie);
  			Integer ind = l.indexOf(prop[i]);
  			
  			fi.set(this,ind);
  			System.err.println("set "+fie+" "+ind);
  		}catch(Exception exc){
  			exc.printStackTrace();
  		}
  	}
  	
  	// TODO Auto-generated method stub
  	
  }
  
  
 
 static String[] st = 
  ("NA18517 ABC7 "+
  "NA18507 ABC8 " +
  "NA18956 ABC9 " +
  "NA19240 ABC10 " +
  "NA18555 ABC11 "+
  "NA12878 ABC12 " +
  "NA19129 ABC13 "+ 
  "NA12156 ABC14").split("\\s+");
 
 static Map<String, String> m = new HashMap<String, String>();
 static List<String> sampleL = new ArrayList<String>();
 static{
	 for(int i=0; i<st.length; i+=2){
		 m.put(st[i+1], st[i]);
		 sampleL.add(st[i]);
	 }
 }

  public String getChrom(String[] str){
		return str[1];
	}
    CompressCooperSeq(File f, String chr, File samples) throws Exception{
       super(f, chr, true, samples, sampleL);
    }
   // public CompressCNVResults(File f, String chr, boolean header, File rawDir) throws Exception {
    	
}
