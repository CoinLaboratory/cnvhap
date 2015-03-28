package lc1.ensj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class PedigreeDataCollection  {
    
 //   public DataCollection dat;
   // final Map<String, Boolean> sex =   new HashMap<String, Boolean>();
    final public Map<String, String> mother;
    final public Map<String, String> father; 
  ///  final Set<String>  males = new HashSet<String>(); 
   // final Set<String> females = new HashSet<String>(); 
    
   
    
    public PedigreeDataCollection( File pedigree){
   //     this.dat = dc;
       mother = new TreeMap<String, String>();
      father = new TreeMap<String, String>();
            try{
        readPedigreeFamilies(pedigree);
            }catch(Exception exc){
                exc.printStackTrace();
            }
          /*  for(Iterator<PhasedIntegerGenotypeData> it = dc.iterator(); it.hasNext();){
                String name = it.next().getName();
               boolean male = sex.get(name);
               if(male){
                  males.add(name);
               }
               else{
                   females.add(name);
               }
            }*/
            
    }
   public  PedigreeDataCollection(){
       mother = new TreeMap<String, String>();
       father = new TreeMap<String, String>();
    }
 
public PedigreeDataCollection(Map<String, String> mother2, Map<String, String> father2) {
  this.mother = mother2;
  this.father = father2;
}

private void readPedigreeFamilies(File f) throws Exception{
    BufferedReader br = new BufferedReader(new FileReader(f));
    String st;
    while((st=br.readLine())!=null){
        String[] str = st.split("\\s+");
        String moth =str[1];
        String fath = str[2];
        String chi = str[0];
        
      mother.put( chi, moth);
      father.put( chi,fath);
    }
    br.close();
}
private void readPedigreeFamiliesHapMap(File f) throws Exception{
       BufferedReader br = new BufferedReader(new FileReader(f));
       String st;
       while((st=br.readLine())!=null){
           String[] father = st.split("\\s+");
           String id_f = father[6].split(":")[4];
           String sex_f = father[4];
          // sex.put(id_f, sex_f.equals("1")); //true ==male
       }
       br.close();
      if(true){
          br = new BufferedReader(new FileReader(f));
          while((st=br.readLine())!=null){
              Map<Integer, String[]>m = new HashMap<Integer, String[]>();
              Integer child_id = null;
              Integer mother_id = null;
              Integer father_id = null;
             for(int i=0; i<3; i++){
                 String[] str = st.split("\\s+");
                 m.put(Integer.parseInt(str[1]), str);
                 if(Integer.parseInt(str[2]) !=0 && Integer.parseInt(str[3])!=0){
                     child_id = Integer.parseInt(str[1]);
                 }
                 else if(str[4].equals("2")){
                     mother_id = Integer.parseInt(str[1]);
                 }
                 else if(str[4].equals("1")){
                     father_id = Integer.parseInt(str[1]);
                 }
                 if(i<2)st = br.readLine();
             }
              String id_m = m.get(mother_id)[6].split(":")[4];
              String id_c = m.get(child_id)[6].split(":")[4];
              String id_f = m.get(father_id)[6].split(":")[4];
              this.mother.put(id_c, id_m);
              this.father.put(id_c, id_f);
          }
          br.close();
      }
  }
public void setTrio(String mother, String father, String child) {
    if(this.mother.containsKey(mother)) throw new RuntimeException("!!");
    if(this.father.containsKey(father)) throw new RuntimeException("!!");
   this.mother.put(child, mother);
   this.father.put(child, father);
   
}
public void print(PrintWriter ped_pw) {
   for(Iterator<String> it = this.mother.keySet().iterator(); it.hasNext();){
       String chi = it.next();
       String moth = mother.get(chi);
       String fath = father.get(chi);
       ped_pw.println(chi+"\t"+moth+"\t"+fath);
   }
   ped_pw.close();
    
}
   
  /* public PedigreeDataCollection(PedigreeDataCollection data1, PedigreeDataCollection data2) throws Exception{
       this.sex.putAll(data1.sex); this.sex.putAll(data2.sex);
       this.mother.putAll(data1.mother); this.mother.putAll(data2.mother);
       this.father.putAll(data1.father); this.father.putAll(data2.father);
       this.males.addAll(data1.males); this.females.addAll(data1.females);
       this.males.addAll(data2.males); this.females.addAll(data2.females);
        dat = new DataCollection(data1.dat, data2.dat);
   }*/

}
