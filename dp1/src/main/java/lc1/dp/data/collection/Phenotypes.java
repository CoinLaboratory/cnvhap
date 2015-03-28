/**
 * 
 */
package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import lc1.stats.ProbabilityDistribution;

public class Phenotypes{
   public  Map<String, Integer>[] phenVals;
    public ProbabilityDistribution[] phenotypeDistribution;
    public List<String> phen = new ArrayList<String>();
 int[] type;
    public Phenotypes(){
        
    }
   public Phenotypes(File incl) throws Exception{
       DataCollection.readPosInfo(incl, new int[] {0}, false,new List[] {phen},new Class[] {String.class});
       this.phenVals = new Map[phen.size()];
       BufferedReader br = DataCollection.getBufferedReader(incl);
       String st = "";
       type = new int[phenVals.length];
       
       for(int i=0; i<phenVals.length; i++){
    	   String sti = br.readLine();
           String[] str = sti.split("\t");
    	   try{
    		  
           type[i] = Integer.parseInt(str[1]);
           if(str.length>2){
               phenVals[i] = new HashMap<String, Integer>();
               for(int j=2; j<str.length; j++){
                   phenVals[i].put(str[j], j-2);
               }
              
           }
       }catch(Exception exc){
    	   exc.printStackTrace();
    	   System.err.println("problem with line "+sti+"\n"+Arrays.asList(str));
    	   System.exit(0);
       }
       }
   }
   public Phenotypes(List<String> header, String string) {
	Set<String> s = new HashSet<String>();
	for(int i=0; i<header.size(); i++){
		if(header.get(i).indexOf(string)>=0){
			s.add(header.get(i).split("_")[0]);
		}
	}
	this.phen = new ArrayList<String>(s);
}
public int size(){
       return phen.size();
   }
   public Map<Double, String> reverse(int i){
       Map<Double, String> m = new HashMap<Double, String>();
if(phenVals[i]==null) return null;
       for(Iterator<Map.Entry<String, Integer>> it = phenVals[i].entrySet().iterator(); it.hasNext();){
           Map.Entry<String, Integer> nxt = it.next();
           m.put(nxt.getValue().doubleValue(), nxt.getKey());
       }
       return m;
   }
   
  public String print(PrintWriter pw) throws Exception{
      StringBuffer format = new StringBuffer();
      for(int i=0; i <phen.size(); i++){
          format.append("%7s \t");
      }
      String fstr = format.toString();
     pw.println("id      \t"+String.format(fstr, phen.toArray()));
      String[] toP = new String[phen.size()+1];
      for(int i=0; i<phenVals.length; i++){
          toP[i] = phenVals[i]==null ? "null" : phenVals[i].toString();
      }
      pw.println("id     \t"+String.format(fstr, toP));
     return fstr;
  }
public int[] type() {
   return type;
}
}