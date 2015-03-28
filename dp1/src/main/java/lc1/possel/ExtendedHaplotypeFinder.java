package lc1.possel;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import lc1.CGH.Aberation;
import lc1.CGH.AberationFinder;
import lc1.CGH.Location;
import lc1.CGH.Locreader;
import lc1.CGH.SNPLocationReader1;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;

public class ExtendedHaplotypeFinder implements Runnable{
    
   //  static int bef = 50;
   // static int aft = 50;
    final File[] files;
    final int chrom;
  
    
    final File  out;
    
   
   
    int thresh = 3;
  
  
  
   
    public void run() {
        try{
     
       
        inner: for(int i=0; i<files.length; i++){
           PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(out, chrom+"_"+i))));
            //int chrom = Integer.parseInt(files[i].getParentFile().getName());
             Locreader snpLocations = new SNPLocationReader1(files[i], chrom+"", null);
              Logger.global.info("file is "+files[i]);
                 List<Integer> locs = (snpLocations).getLocs();
                SimpleDataCollection sdt = 
                 SimpleDataCollection.readFastPhaseOutput((AberationFinder.getBufferedReader(files[i], "phased.txt")),Emiss.class, Emiss.getEmissionStateSpace(1));
                 sdt.loc = locs;
                 sdt.split();
               EHHFinder ehh = new EHHFinder(sdt);
                List<Aberation> obj = sdt.getDeletedPositions(true);
                for(int i1=0; i1<obj.size(); i1++){
                    Aberation ab  = obj.get(i1);
                    if(ab.end - ab.start < 3) continue;
                    System.err.println("aBERATION " +ab);
                    ehh.extend(ab, pw);
             
                    pw.println("##############################################");
                    pw.flush();
//                    List<String> indiv = sdt.getNames();
                      }
                // PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(
               //          new File(out, "reg_"+region.toStringPrint()+".txt"))));
                // for(int i1=0; i1<sdt.loc.size(); i1++){
                 //    if(sdt.loc.get(i1) < region.min) pw.print("-");
                 //    else if(sdt.loc.get(i1) > region.max) pw.print("+");
                 //    else pw.print("=");
                // }
                // pw.println();
              //   sdt.writeLocation(pw);
              //   sdt.writeFastphase(sdt.data, sdt.uncertainty, pw,false);
              //   sdt.writeLocation(pw);
              //   pw.close();
              //   System.err.println("matched region "+region+" at "+files[i]);
             //  continue outer;
             //}
                pw.close();
        }
        
        // Logger.global.warning("no match for "+region);
         
        }catch(Exception exc){
            exc.printStackTrace();
        }
      }
    public  ExtendedHaplotypeFinder(File dir, File[] user, String[] args, PrintWriter out, 
            Map<String, Number> sum, String chromosome, Location reg) throws Exception{
       
       
        this.out = new File(dir.getParentFile(), "interesting regions");
        if(!this.out.exists()) this.out.mkdir(); 
        this.chrom = Integer.parseInt(chromosome);
        this.files = user;
    
        
    }
  

  /*returns first position greater than */
private static Integer getPos(List<Integer> locs, long min) {
   for(int i=0; i<locs.size();i++){
       if(locs.get(i)>min) return i;
   }
   return 0;
}
}
