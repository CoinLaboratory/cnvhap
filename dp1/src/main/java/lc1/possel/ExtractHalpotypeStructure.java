package lc1.possel;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import lc1.CGH.AberationFinder;
import lc1.CGH.AbstractAberatiionReader;
import lc1.CGH.Location;
import lc1.CGH.Locreader;
import lc1.CGH.MultiProbeAberationReader;
import lc1.CGH.SNPLocationReader1;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;

public class ExtractHalpotypeStructure implements Runnable{
    static int bef = 50;
    static int aft = 50;
    final File[] files;
    final int chrom;
    final AbstractAberatiionReader abR;
    
    final File  out;
    
    public void run() {
        try{
        outer: for(Iterator<Location> it = abR.iterator(); it.hasNext();){
            Location region = it.next();
       
        inner: for(int i=0; i<files.length; i++){
            //int chrom = Integer.parseInt(files[i].getParentFile().getName());
            
           if(chrom!=Integer.parseInt(region.chr))  continue inner;
             Locreader snpLocations = new SNPLocationReader1(files[i], chrom+"", null);
           System.err.println(files[i]+" "+snpLocations.getFirst().min+" "+snpLocations.getLast().max+" "+region);
             if(snpLocations.getFirst().min < region.min &&snpLocations.getLast().max > region.max ){
              Logger.global.info("file is "+files[i]);
                 List<Integer> locs = (snpLocations).getLocs();
                 int posMin = getPos(locs, region.min)-1;
                 int posMax = getPos(locs, region.max);
                 posMin = Math.max(0, posMin-bef);
                 posMax = Math.min(locs.size(), posMax+aft);
                 SimpleDataCollection sdt = 
                 SimpleDataCollection.readFastPhaseOutput((AberationFinder.getBufferedReader(files[i], "phased.txt")),Emiss.class, Emiss.getEmissionStateSpace(1));
                 sdt.loc = locs.subList(posMin, posMax);
                 for(Iterator<PIGData> it1 = sdt.data.values().iterator(); it1.hasNext();){
                     PIGData nxt = it1.next();
                     nxt.restrictSites(posMin, posMax);
                     if(nxt.length()!=sdt.loc.size()) throw new RuntimeException("!!");
                 }
                
                 PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(
                         new File(out, "reg_"+region.toStringPrint()+".txt"))));
                 for(int i1=0; i1<sdt.loc.size(); i1++){
                     if(sdt.loc.get(i1) < region.min) pw.print("-");
                     else if(sdt.loc.get(i1) > region.max) pw.print("+");
                     else pw.print("=");
                 }
                 pw.println();
                 sdt.writeLocation(pw, null);
                 sdt.writeFastphase( new File(out, "reg_"+region.toStringPrint()+".txt"),false);
                 sdt.writeLocation(pw, null);
                 pw.close();
                 System.err.println("matched region "+region+" at "+files[i]);
               continue outer;
             }
        }
       
        // Logger.global.warning("no match for "+region);
         
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
      }
    public  ExtractHalpotypeStructure(File dir, File[] user, String[] args, PrintWriter out, 
            Map<String, Number> sum, String chromosome, Location reg) throws Exception{
       
        abR = new MultiProbeAberationReader(Integer.MAX_VALUE, ""){
            public boolean exclude(String[] str){
                return false;
            }
        };
        abR.initialise(dir, null, null,0,0, "valsum1507_1.txt", null);
        this.out = new File(dir.getParentFile(), "reg1");
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
