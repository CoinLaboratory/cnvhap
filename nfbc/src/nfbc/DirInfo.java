/**
 * 
 */
package nfbc;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

class DirInfo{
   String[] fields;
     Map<String, Map<Integer, Integer>> header = new HashMap<String, Map<Integer, Integer>>();
     Map<String, List<String[]>> results = new HashMap<String, List<String[]>>();
     String file;
     String pheno;
     String experiment;
     
     
     public Iterator<String[]> getIterator(String key){
         Iterator<String[]> fR = results.get(key).iterator();
         String[] head =  fR.next();
         Map<Integer, Integer> m = header.get(key);
         for(Iterator<Integer> it = m.keySet().iterator(); it.hasNext();){
             Integer j = it.next();
             Integer i = m.get(j);
             if(!fields[j].equals(head[i])){
                 throw new RuntimeException("!!");
             }
         }
        
         return fR;
     }
     public int[][] getOrder(List<String> keys){
         int[][] order = new int[fields.length][];
         outer: for(int i=0; i<order.length; i++){
             for(int j=0; j<keys.size(); j++){
                 Integer res = header.get(keys.get(j)).get(i);
                 if(res!=null){
                     order[j] = new int[] {j, res};
                     continue outer;
                 }
             }
         }
         return order;
     }
     
     public void check(int k){
         String id = null;
         for(Iterator<List<String[]>> it = results.values().iterator(); it.hasNext();){
             if(id==null) id = it.next().get(k)[0];
             else if(!id.equals(it.next().get(k)[0])) throw new RuntimeException("!!");
         }
     }
     /** assume fields[0] is common */
     public void print(PrintWriter pw){
         List<String> keys = new ArrayList<String>(header.keySet());
         int[][] order = getOrder(keys);
         
         String fK = keys.get(0);
         int size = results.get(fK).size();
         String[] toPrint = new String[fields.length];
         for(int k=0; k<size; k++){
             check(k);
             for(int i=0; i<toPrint.length; i++){
                 int[] ind = order[i];
                 if(ind==null){
                     toPrint[i] = "";
                 }
                 else{
                     toPrint[i] = results.get(ind[0]).get(k)[ind[1]];
                 }
             }
             print(pw, toPrint);
             
             for(int i=1; i<keys.size(); i++){
                 
             }
         }
     }
    private void print(PrintWriter pw, String[] toPrint) {
        int i=0;
       for(; i<toPrint.length-1; i++){
           pw.print(toPrint[i]+"\t");
       }
       pw.print(toPrint[i]+"\t");
        
    }
 }