package lc1.ensj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class Analysis {
    
    
  static String key = "gaps";
  static int pos = 4;
  public static void main(String[]args){
      try{
          File dir = new File(".");
          File[] f =  dir.listFiles(new FileFilter(){
             public boolean accept(File pathname) {
                return pathname.isDirectory();
             }
               
           });
          List<Double> set = new ArrayList<Double>();
      //  Map<Long, Double> m = new TreeMap<Long, Double>();
          for(int i=0; i<f.length; i++){
              File summ = new File(f[i], "logfile.txt");
              if(!summ.exists()) continue;
              BufferedReader br = new BufferedReader(new FileReader(summ));
              String st = "";
              inner: while((st = br.readLine())!=null){
                  if(st.startsWith(key)){
                     // System.err.println(st);
                     String[] str = st.split("\\s+");
                     set.add(Double.parseDouble(str[pos]));
                    // m.put(Long.parseLong(f[i].getName()), Double.parseDouble(str[3]));
                      
                      break inner;
                  }
              }
              br.close();
          }
          Collections.sort(set);
         for(Iterator< Double> it = set.iterator(); it.hasNext();){
              System.err.println(it.next());
          }
         System.err.println("median "+set.get((int) Math.round((double)set.size() / 2.0)));
      }catch(Exception exc){
          exc.printStackTrace();
      }
  }
}
