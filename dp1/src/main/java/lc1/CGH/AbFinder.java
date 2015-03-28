package lc1.CGH;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.IlluminaRDataCollection;

public class AbFinder implements Runnable{
    boolean ran = true;
    AbstractAberatiionReader agilent =    new MultiProbeAberationReader(Long.MAX_VALUE, "agilent");
    Locreader snp = new Locreader(Long.MAX_VALUE, "snp");
    Locreader french = new Locreader(Long.MAX_VALUE, "french");
    Locreader cgh = new Locreader(Long.MAX_VALUE, "cgh");
    int ind =1;
    
    SortedSet<Integer> snps = new TreeSet<Integer>();
    SortedSet<Integer> probes = new TreeSet<Integer>();
    Set<String> individuals = new HashSet<String>();
    //Map<String, Number> sum;
    PrintWriter out;
    String chr;
    
    File logDir;
  //  PrintWriter log;
//    Set<String> abIndividuals = new HashSet<String>();
//    Set<String> cnpIndividuals = new HashSet<String>();
    
    public static List<Integer> readPos(File user, String file) throws Exception{
        BufferedReader br = AberationFinder.getBufferedReader(user, file);
        List<Integer> res = new ArrayList<Integer>();
        DataCollection.readPosInfo(br, new int[] {1}, true, new List[] {res}, new Class[] {Integer.class});
        br.close();
        return res;
    }
    public  AbFinder(File dir, File[] user, String[] args, PrintWriter out, 
            Map<String, Number> sum, String chromosme, Location reg) throws Exception{
        System.err.println("chromosome is "+chromosme);
        File data = new File(dir.getParentFile(), "../data");
        this.chr ="1";//chromosme;
        logDir  = new File("logDir");
        if(!logDir.exists()){
            logDir.mkdir();
        }
        BufferedReader indv = new BufferedReader(new FileReader(new File(dir, "indiv_shared.txt")));
        String st1 = "";
        while((st1 = indv.readLine())!=null){
            this.individuals.add(st1.trim());
        }
        indv.close();
        for(int i=0; i<user.length; i++){
                this.snps.addAll(readPos(user[i], "phased1.txt_"+ind));
                this.probes.addAll(readPos(user[i], "phased1.txt_"+3));
                SNPAberrationReader snps = 
                    new SNPAberrationReader(user[i], chromosme, 5, 
                        null, "cnv.txt_"+ind, AberationFinder.threshold,AberationFinder.limit() ? individuals : null, "snps");
                SNPAberrationReader cgh = 
                    new SNPAberrationReader(user[i], chromosme, AberationFinder.noOfSnps, 
                        null, "cnv.txt_3", AberationFinder.threshold,AberationFinder.limit() ? individuals : null, "cgh");
                this.snp.addAll(snps);
                this.cgh.addAll(cgh);
         }
        snp.sort();
        snp.merge(0);
        cgh.sort();
        cgh.merge(0);
        this.snp.probes = snps;
        this.cgh.probes = probes;
        snp.setChr(chr);
        cgh.setChr(chr);
        File cghFile = new File(dir.getParentFile(),"data/"+chr+"_cghdata.txt");
        File snpF = new File(dir.getParentFile(),"data/"+chr+"_data.txt");
       DataCollection cghdat = null;// new CGHDataCollection(cghFile);  
        IlluminaRDataCollection illdat = null;//new IlluminaRDataCollection(snpF,(short)0, 2);
        this.agilent.probes = probes;
        cgh.dc = null;//cghdat;
        snp.dc = illdat;
      //  this.chr = chromosome;
        Location region = new Location(chr, Math.min(probes.first(), snps.first()), Math.max(probes.last(), snps.last()));
        System.err.println("region is "+region);
       
       // this.sum = sum;
        this.out = out;
            
        agilent.initialise(dir, chr, region, 0,0,"MultiProbeByIndividual.txt", AberationFinder.limit() ? individuals : null);
        agilent.restrict();
        agilent.sort();
        
        
        french = new SNPAberrationReader(dir, chr, AberationFinder.noOfSnps, null, "FrenchSamples.txt", AberationFinder.threshold, AberationFinder.limit() ? individuals : null, "french");
        this.french.probes = snps;
        if(true){ // do we merge?
        snp.mergeNames();
        cgh.mergeNames();
        
        agilent.mergeNames();
        french.mergeNames();
     
        
        snp.sort();
        cgh.sort();
     
            double mergeThresh = 0.0;
            snp.merge(mergeThresh);
            cgh.merge(mergeThresh);
            agilent.merge(mergeThresh);
            snp.merge(mergeThresh);
            snp.removeWithObsLessThan(2);
        }
        french.restrict();
        snp.restrict();
        agilent.restrict();
        cgh.restrict();
    }
    
   

    public void run() {
     try{
         {
             LocCompare lc1 = new LocCompare(french, cgh, logDir);
             }
        LocCompare lc = new LocCompare(snp, cgh,logDir);
         
         {
             LocCompare lc1 = new LocCompare(agilent, cgh, logDir);
             }
         {
             LocCompare lc1 = new LocCompare(french, snp, logDir);
             }
         {
             LocCompare lc1 = new LocCompare(french, agilent, logDir);
             }
         {
             LocCompare lc1 = new LocCompare(snp, agilent, logDir);
             }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    //  log.close();
        
    }
    
  
    

}
