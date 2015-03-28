package lc1.CGH;

import java.io.BufferedReader;
import java.util.Set;
import java.util.logging.Logger;

public class ExcelAberationReader extends AbstractAberatiionReader {
    
    //244k file
    //AberrationNo  Chr Cytoband    Start   Stop    NoOfProbes  Amplification   Deletion    (-)log10pval
    
    //185k file
    //AberrationNo  Chr Cytoband    Start   Stop    NoOfProbes  Amplification   Deletion    (-)log10pval


    public ExcelAberationReader(long lengthLim) {
        super(lengthLim, "");
        // TODO Auto-generated constructor stub
    }



    static String[] cols = new String[] {"aaaa", "Chr", "Start", "Stop", "Mean", "NoOfProbes", "Amplification", "Deletion"};
    
    String name;
    @Override
    public String getName(String[] str){
        return name;
     }
     public String getChr(String[] str){
         return str[col[1]].substring(3);
     }
  
   
     public double getNoCopy(String[] str){
         double ampl =Double.parseDouble( str[col[6]]);
         double del =Double.parseDouble( str[col[7]]);
        return Math.abs(ampl)<0.001 ? del: ampl;
     }
     @Override
     public void initialise( BufferedReader dir,String chromosome, Location region,  int st1, int st2, String name, Set<String> indiv)throws Exception{
         if(name.indexOf('/')>=0){
             this.name = name.substring(name.lastIndexOf('/')+1);
         }
         this.name = this.name.split("_")[0];
         if(indiv!=null && !indiv.contains(this.name)) {
             Logger.global.info("excluding "+this.name);
             return;
         }
         if(this.name.length()<4) throw new RuntimeException(" "+name);
         super.initialise(dir, chromosome, region,st1,st2, this.name,  indiv);
     }
     
    
      public String getStart(String[] str){
          return str[col[2]];
      }
      public String getEnd(String[] str){
          return str[col[3]];
      }
      public int getNoProbes(String[] str){
          return Integer.parseInt(str[col[5]]);
      }  
    
      
     
@Override
public String[] getCols() {
   return cols;
}
    

   

}
