package lc1.CGH;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;



public class SNPAberrationReader extends Locreader {
    
  
  

  

    static int gap = 0;
   // public Integer from=null;
   // public Integer to=null;
   // static boolean merge =false;
    //int loc_int = 3;
    public SNPAberrationReader (File dir, String chromosome, int noOfSnps,  Location startEnd, String file, 
            double thresh, Set<String> indiv, String name)throws Exception{
        super(Long.MAX_VALUE, name);
        BufferedReader br = AberationFinder.getBufferedReader(dir,file);
        String st = br.readLine();
       if(st.startsWith("from")){
           String[] str = st.split("\\s+");
      //     from = Integer.parseInt(str[1]);
      //     to = Integer.parseInt(str[3]);
           st = br.readLine();
       }
       List<Location> list = new ArrayList<Location>();
        String chr = "Chr";
        String[] names = st.split("\\s+");
        while((st = br.readLine())!=null){
            if(st.startsWith("from") || st.startsWith("File")) {
                System.err.println("removed "+st);
                continue;
            }
            String[] str = st.split("\\s+");
            if(str.length==1) str = st.split(",");
            String chrom = str[2];
            if(chrom.startsWith("pbs")) chrom = chromosome;
          //  if(str[1].startsWith("NA")) {
          //      System.err.println("removed "+st);
          //      continue;
         //   }
            if(chromosome!=null && !chrom.equals(chromosome) && ! chrom.equals(".")){
              //  System.err.println("removed "+st);
                continue;
            }
      
            long min =Long.parseLong(str[3]);
           long max = Long.parseLong(str[4]);
           int nocop =str[7].startsWith("DEL") ? 0 :
               str[7].startsWith("CNV") ? 2:
               Integer.parseInt(str[7]);
           double cert = str.length>8 ? Double.parseDouble(str[8]) : 1.0;
           int noSnps = Integer.parseInt(str[5]);
           
           if(noSnps<noOfSnps) continue;
           if(noOfSnps>1 && max==min) throw new RuntimeException("!!");
            Location location = max > min ? new Location(chrom, min-gap,max+gap) :new Location(chrom, max-gap,min+gap) ;
          
               int ind = str[1].indexOf('_');
               String nme  = (ind>=0)  ? str[1].substring(0,ind) : str[1];
            
            location.setNoCop(nocop);
          location.setName(nme);
          location.setNoSnps(noSnps);
          if(noSnps < 1) throw new RuntimeException("!!");
         if(indiv!=null && ! indiv.contains(nme)) continue;
          
               if(cert > thresh && (startEnd==null || startEnd.overlaps(location)>=0)) add(location);
               else{
                   System.err.println("removed "+st);
               }
        }
      sort();
        br.close();
        //thin(cntThresh);
    }


   
  
}
