package lc1.dp.appl;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.Sampler;
import lc1.dp.data.collection.DataC;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.external.Fastphase;
import lc1.dp.model.HaplotypeHMMIterator;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.EmissionState;
import lc1.possel.EHHFinder;
import lc1.util.Constants;
import lc1.util.Executor;
import pal.statistics.ChiSquareTest;

public class RunFastPhase {
    static Logger logger = Logger.getAnonymousLogger();
    static Integer[] num;
    
    static Calendar cal = new GregorianCalendar();
  static Double[] chiSq(Integer[][] vit_control, Integer[][] vit_cases){
    Double[] p = new Double[vit_control[0].length];
    for(int i=0; i<p.length; i++){
        List<Number >expected = new ArrayList<Number>();
        List<Integer> obs = new ArrayList<Integer>();
        double sum=0;
        int len = vit_control.length;
        for(int j=0; j<len; j++){
            if(vit_control[j][i]!=0 || vit_cases[j][i]!=0){
                expected.add(vit_control[j][i]);
                obs.add(vit_cases[j][i]);
                sum+=vit_control[j][i];
            }
        }
        double[] ex = new double[expected.size()];
        int[] o =  new int[expected.size()];
        for(int j=0; j<expected.size(); j++){
            ex[j] =  expected.get(j).doubleValue()/sum;
            o[j] = obs.get(j);
        }
        p[i] = ChiSquareTest.compare(ex, o);
    }
    return p;
  }
    
 
  
static String[] tag_report =   new String[] {"snpid","loc", "index",   "maf"}; //"hwe","regression", "regrP", 
static String[] tag_pheno = new String[] {"chisq_state", "armitage_state","odds_state",
        "cases_state",
        "controls_state"};

    public static Runnable getMonitor(final Sampler sampler,  final DataC original, 
            final DataC affy_reference,
           final File clusterFile,  final File parentfile, final String[] toD){
        return new Runnable(){
            Map<String, boolean[]> unc =null;//((DataCollection)affy_reference).getUncertainPositions();
            public void run(){
                runSample();
             //   if(Constants.summarise()) runSummary();
            }
           // Sampler sampler = new Sampler(copy(obj, true));
            public void runSample() {
                try{
                
                if(Constants.sample()){
                    List<Integer> list = new ArrayList<Integer>();
                    if(Constants.keepBest()){
                        list.add(0);
                    }
                    else{
                    for(int i=0; i<Constants.numRep(); i++){
                        list.add(i);
                    }
                    }
                    logger.info("starting sample");
                    if(Constants.numRep()>1 || Constants.noSamples()>1) sampler.calcBestPathSampling(list);
                    ((DataCollection)sampler.data).dropIndiv(toD);
                 //   if(Constants.dropFixed()) ((DataCollection)sampler.data).dropFixed();
                    logger.info("ending sample");
                  //  PrintWriter pw_hap = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "phased"+".txt"))));
                   // pw_hap.println(cal.getTime());
                   // PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "phased1"+".txt"))));
                  //  pw_hap.println(cal.getTime());
                    File pw_hap2 = new File(parentfile, "phased2");
                    pw_hap2.mkdir();
                  if(Constants.savePhasedConfiguration())  ((DataCollection)sampler.data).writeFastphase(pw_hap2, false);
                    
                    
                    
                  if(Constants.writeRes() && Constants.savePhasedConfiguration()){
                	  ((DataCollection)sampler.data).writeCompressed(new File(parentfile, "res"),true);
                  }
                  String[] writAvg = Constants.writeAverages(false);
                  if(writAvg!=null &&writAvg.length>0 && !writAvg[0].equals("null") && Constants.savePhasedConfiguration()){
                	  File avg = new File(parentfile, "avg");
                	  avg.mkdir();
                	  ((DataCollection)sampler.data).writeAverages(avg, writAvg);
              	
                  }
                  ((DataCollection)sampler.data).finishedPrinting();
                  
                 //   PrintWriter pw_hap3 = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "phased3"+".txt"))));
                  
                  
                 //   Collection<Integer> toD =  Constants.drop() ? sampler.data.calculateMaf1().getConstantPos() : null;
                 //   sampler.data.writeLocation(pw_hap, toD);
               //sampler.data.writeDickFormat1(new File(parentfile, "phased3"+".txt"), true);
                   // sampler.data.printHapMapFormat(pw_hap, toD,false);
               //    if(false)  ((DataCollection)sampler.data).calcHWE(false, false);
                   // writeFastphase(sampler.data.data,sampler.data.uncertainty,sampler.data.uncertaintyPhase, pw_hap,  true, true, toD);
                   // PrintWriter pw_st = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "phased_states"+".txt"))));
                  // ((DataCollection) sampler.data).printHapMapFormat(pw_hap1, toD, true, 
                  //         tag_report,
                   //        tag_pheno,
                          //  : new String[] {"snpid", "loc"}, 
                     //      new String[] {}, "%7s");
                  //  sampler.data.writeLocation(pw_hap2,Constants.drop() ?  sampler.data.calculateMaf1().getConstantPos(): null);
                   //pw_hap1.close();
                  
                 //   sampler.data.printHapMapFormat( 
                  //          pw_st, toD, true, 
                  //          new String[] {"snpid", "loc"}, 
                 //           new String[] {"viterbi", "dataLSt", "uncertaintyVitPhase"}, "%7s");
                  // if(Constants.saveStates()) sampler.data.writeFastphase(((DataCollection)sampler.data).viterbi, ((DataCollection)sampler.data).uncertaintyVitPhase, 
                          //  pw_st, false, true, false,null);
                  //  pw_st.close(); 
                    //pw_hap.close();
                  //  PrintWriter pw_del =  new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "cnv"+".txt"))));
                  //  PrintWriter pw_del1 =  new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "cnv1"+".txt"))));
                    //PrintWriter pw_loc =  new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "loc"+".txt"))));
                  if(Constants.savePhasedConfiguration())  sampler.data.printDeletedPositions(parentfile);
                  //  sampler.data.writeLocation(pw_loc, toD);
                  //  pw_del.close();
                  //  pw_del1.close();
                 //   PrintWriter pw_del =  new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "cnv"+".txt"))));
                 //   PrintWriter pw_loc =  new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "loc"+".txt"))));
                 //   sampler.data.printDeletedPositions(pw_del);
                  //  sampler.data.writeLocation(pw_loc, toD);
                 //   pw_del.close();
                  //  pw_loc.close();
                }
                }catch(Exception exc){
                    exc.printStackTrace();
                }
            }
         

             private boolean contains(Set<Integer> trues, Integer pos) {
                 int gap = 2;
               for(int i=pos-gap; i<=pos+gap;i++ ){
                   if(trues.contains(i)) return true;
               }
               return false;
            }

            private double getMax(Collection <Double> set){
                 double max = 0;
                 for(Iterator <Double > it = set.iterator(); it.hasNext();){
                     Double d = it.next();
                     if(d>max){
                         max = d;
                     }
                 }
                 return max;
             }
           
              
         
        };
    }
     

    static String[] tag1 = new String[] {"snpid","loc", "index", "majorAllele", "minorAllele","maf"}; //  "hwe_exclMissing"
    public static void run( DataCollection originalData,final  File  summar, 
            File clusters,   File parentfile) throws Exception{
    //	originalData.restrictToIds(Constants.reportIds());
    	
        if(originalData==null) throw new RuntimeException("!!");
        
      /* if(((DataCollection)originalData instanceof MergedDataCollection) && Constants.useDataAsModel()>0){
        	
        }*/
        
     //   EmissionState ems = originalData.dataL.get("J1d2");        
    long time = System.currentTimeMillis();
       // PrintWriter summary = null;
   //   EmissionState ems = originalData.dataL.values().iterator().next();
      
         /* if(Constants.format()[0].equals("illumina") && originalData instanceof LikelihoodDataCollection && Constants.allowCloning()){
              originalData.clearData();
              ((LikelihoodDataCollection)originalData).calculateMLGenotypeData(Constants.allowCloning());
           }*/
        
        DataCollection copy =//Constants.allowCloning() ? (DataCollection)originalData.clone() : 
        	(DataCollection) originalData;
     if(Constants.no_hmm()){
    	 if(true) throw new RuntimeException("!!");
    	 copy.calculateMLGenotypeData(true);
    	  copy.writeCompressed(new File(parentfile, "res_ML"),true);
    	  if(Constants.plot()==1) System.exit(0);
    	  Executor.shutdown();
    	  return;
     }
     else if(Constants.run() ){
        if(true){
           /*if(Constants.trainWithPedigree() || Constants.sampleWithPedigree()){
               copy.arrangeDataAccordingToPedigree();
           }*/
        	
           String col = Constants.column();
           String[] writAvg = Constants.writeAverages(true);
           if(writAvg!=null && writAvg.length>0 && !writAvg[0].equals("null")){
        	   File avg1 = new File(parentfile,"avg1");
               avg1.mkdir();
               ((DataCollection)originalData).writeAverages(avg1, writAvg);
               
       	
           }
          
       //    PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter( new File(parentfile, "reference"+".txt"))));
           
         //  PrintWriter pw_hap3 = new PrintWriter(new BufferedWriter(new FileWriter( new File(parentfile, "pheno"+".txt"))));
         //  PrintWriter pw_indiv = new PrintWriter(new BufferedWriter(new FileWriter( new File(parentfile, "indiv"+".txt"))));
        Sampler sampler = new Sampler(copy, (new File(parentfile,"samples")), true);//
    
       // Collection<Integer> toD =   Constants.drop() ?  sampler.data.calculateMaf1().getConstantPos() : null;
     //   if(false) ((DataCollection)originalData).calcHWE(false, true);
    
       originalData.writeSNPFile(new File(parentfile, "snp_"+originalData.name+".txt"), Constants.chrom0(), false, null);
     
     //   pw_indiv.close();
       // ((DataCollection)originalData).printPheno(pw_hap3);
     //   pw_hap3.close();
      //rreference.data, null,null, pw_hap1, false, false, toD);
      //  pw_hap.close();
     //   pw_hap1.close();
       // pw_hap2.close();
      //  monitor.run();
      
      
      
       ((DataCollection)sampler.data).initialisePrinting(parentfile
       );
      if(!Constants.resample() && !Constants.no_hmm()){
    	
    	  Fastphase fp= new Fastphase(copy,sampler, Constants.numRep(), parentfile);
        Iterator<MarkovModel> hit;
     /*   if(Constants.transMode(1)!=null){
           // int[] transMode = ;
            final HaplotypeHMMIterator counts = new HaplotypeHMMIterator("counts",originalData.length(), Constants.numRep(), 
                    (DataCollection)originalData,Constants.modify(0).length,
                    Emiss.countSpace(), Constants.transMode(0), Constants.modify(0), Constants.modifyFrac(0), Constants.modifyFracStart(), originalData.numLevels());
            final HaplotypeHMMIterator alleles = new HaplotypeHMMIterator("alleles",originalData.length(), Constants.numRep(), 
                   (DataCollection) originalData, Constants.modify(1).length,
                    Emiss.alleleSpace(), Constants.transMode(1), Constants.modify(1),  Constants.modifyFrac(1),Constants.modifyFracStart(),originalData.numLevels());
             hit = new Iterator<MarkovModel>(){

            public boolean hasNext() {
                return alleles.hasNext() && counts.hasNext();
            }

            public MarkovModel next() {
                MarkovModel cnt=counts.next();
                MarkovModel allele = alleles.next();
               
               return 
             //  new CachedHMM(
               new PairMarkovModel(new MarkovModel[] {cnt,allele }, new int[] {0,1}, AlleleCopyPairEmissionState.class,
               Emiss.stateSpace() , false       
               //)
               );
            }

            public void remove() {}
               
           };
        }*/
       // else{
       String[] l =  Constants.useDataAsModel;
       if(l!=null){
       List<String> l_new = new ArrayList<String>();
       for(int k=0; k<l.length; k++){
    	   if(l[k].endsWith("*")){
    		   String st =l[k].substring(0,l[k].length()-1);
    		   Iterator<String> it;
    		   if(originalData instanceof MergedDataCollection){
    			   it = ((MergedDataCollection)originalData).ldl[0].getKeyIterator();
    		   }else it = originalData.getKeyIterator(); 
    		   for(; it.hasNext();){
    			   String nxt = it.next();
    			   if(nxt.startsWith(st)) l_new.add(nxt);
    		   }
    	   }else{
    		   l_new.add(l[k]);
    	   }
       }
       Constants.useDataAsModel = l_new.toArray(new String[0]);
       }
       l = Constants.useDataAsModel();
           hit = new HaplotypeHMMIterator("", originalData.length(), Constants.numRep(), 
            		(DataCollection)originalData, Constants.modify(0).length,
               
                    Constants.transMode(0), Constants.modify0, Constants.modifyFrac(0),
                    Constants.modifyFracStart(), originalData.numLevels(), null);
           if(false){
        	   String[][] modify1 = Constants.modify1; 
        	  
        	   double[] modFrac1St = new double[modify1.length];
        	   Arrays.fill(modFrac1St, 1.0/(double)modFrac1St.length);
        	  HaplotypeHMMIterator alleles = new HaplotypeHMMIterator("", originalData.length(), Constants.numRep(), 
               		(DataCollection)originalData, Constants.modify1[0].length,
                  
                       Constants.transMode(1), modify1, modFrac1St,
                      modFrac1St, originalData.numLevels(),null);
        	   hit = HaplotypeHMMIterator.mergeModel(hit, alleles);
           }
        //}
     //   hit.set(init);
        Logger.global.info("mem "+Runtime.getRuntime().freeMemory());
      //  EmissionState ems = originalData.dataL.get("NA19240");
        Callable call = fp.train(fp.getTrainingElements(hit));
        List tasks = Arrays.asList(new Callable[]{call});
        if(Constants.numThreads()==0 || Constants.plot()<=1){
       for(int i=0; i<tasks.size(); i++){
            ((Callable)tasks.get(i)).call();
        }
        }
        else{
            BaumWelchTrainer.involeTasks(tasks, true);
        }
      }
    
     //summary = new PrintWriter(new BufferedWriter(new FileWriter(summar)));
      String[] toDI = 	originalData.restrictToIds(Constants.reportIds());
      Runnable monitor = getMonitor(sampler, originalData, originalData,  clusters,  parentfile, toDI);
    //  summary.println(  cal.getTime().toString().toUpperCase()+"------------------------------");
    //  summary.println(Constants.transMode(0)+" "+Constants.transMode(1)+" "+Constants.numRep()+" ");
       monitor.run();
        if(Constants.plot()>1 && false)EHHFinder.writeEHH((DataCollection)sampler.data, new File(summar.getParentFile(), "sweep"),Constants.chrom0()+"", null, ""
                                                                                                                ); 
       }
        }
      /*  if(Constants.runFastPhase()){
            if(summary==null){
                summary = new PrintWriter(new BufferedWriter(new FileWriter(summar)));
            }
       //     DataCollection result = FastPhaseFormat.runFastphase(copy, originalData, originalData,Constants.numF(0), Constants.numRep(),
           //         Constants.sum(Constants.numIt()), summary);
        //    result.writeCompressed(new File(parentfile, "res"));
         //   logger.info("finished in "+(System.currentTimeMillis()-time));
        }*/
       // sampler.finalise();
      //  if(Constants.inputFile.endsWith(".lhood")){
        //    PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, Constants.inputFile+"_summ1"))));
       //     ((DataCollection)originalData).print(pw1);
       //     pw1.close();
      //  }
        logger.info("finished in "+(System.currentTimeMillis()-time));
        
          //    summary.close();
        deleteEmptyFiles(parentfile);
              
    }
    
   
   
 
  
    private static void deleteEmptyFiles(File parentfile) {
		File[] f = parentfile.listFiles();
		for(int k=0; k<f.length; k++){
			if(f[k].isDirectory()){
				deleteEmptyFiles(f[k]);
				if(f[k].listFiles().length==0) f[k].delete();
			}else{
				if(f[k].length()==0) f[k].delete();
			}
			
		}
	}





	static Boolean get(String st){
        if(st.equals("-")) return null;
       int i =  Integer.parseInt(st);
       if(i==0) return zero;
       else if(i==1) return one;
       else throw new RuntimeException("!!");
    }
    static Boolean one = new Boolean(true);  static Boolean zero = new Boolean(false);
    
}

/* public void runSummary() {
logger.info("starting run");
try{

    double tp,fp,tn, fn;
    tp = fp =tn =fn =0;
//   Integer[] tpfp = new Integer[] {0,0,0,0};  //[truepos, trueneg, falsepos, falseneg];
   // summary.println(msg);
   if(Constants.sample()){
      
      if(Constants.sampleWithPedigree()){
          //sampler.calcRecSitesSampling();
       //   if(sampler.data.recSites!=null && sampler.data.recSites.size()>0){
              
              for(Iterator<String> it = sampler.data.getKeys().iterator(); it.hasNext();){
                  String j = it.next();
                  SortedMap<Integer, Integer>[] rec2 = sampler.data.recSites(j);
                  System.err.println(j);
                  String[] j1 = j.split(";");
                  SortedMap<Integer, Integer>[] rec1 = original.recSites(j1[j1.length-1]);
                  if(rec1==null) continue;
                //  if(rec1!=null) {
                          
                          summary.println("orig vs inferred"+j1[j1.length-1]);
                         // Set<Integer>[] sw1 = rec1.getSwitches();
                         // Set<Integer>[] sw2 = rec2.getSwitches();
                          
                          for(int k=0; k<rec1.length; k++){
                              summary.println("orig: "+rec1[k]);
                              List<Entry<Integer, Integer>> l = new ArrayList<Entry<Integer, Integer>>(rec2[k].entrySet());
                              Collections.sort(l, comp);
                              summary.println("inf : "+l);
                              Set<Integer> trues = rec1[k].keySet();
                              Set<Integer> positives = rec2[k].keySet();
                             Set<Integer> negatives = new HashSet<Integer>();
                             for(int i=0; i<sampler.data.length(); i++){
                                 negatives.add(i);
                             }
                           
                           //  negatives.removeAll(positives);
                              for(Iterator<Integer> it1 = positives.iterator(); it1.hasNext();){
                                  Integer pos = it1.next();
                                  if(contains(trues,pos)) tp++;
                                  else{
                                      System.err.println("false positive "+pos+" cf "+trues);
                                      fp++;
                                  }
                              }
                              for(Iterator<Integer> it1 = negatives.iterator(); it1.hasNext();){
                                  Integer fa = it1.next();
                                  if( contains(positives, fa)){continue;}
                                  else if(trues.contains(fa) ){
                                      System.err.println("false negative "+fa+" cf "+trues);
                                      fn++;
                                  }
                                  else tn++;
                              }
                          }
                         
                    //  }
                  
                 
              }
         // }
      }
      summary.println("tp:fn:tn:fp "+tp+" "+fn+" "+tn+" "+fp);
      if(tp+fn>0)  summary.println("sensitivity "+tp / (tp+fn));
      if(tn+fp>0) summary.println("specifitity "+tn/ (double)(tn+fp));
      if(tp+fp>0)summary.println("ppv "+tp/ (double)(tp+fp));
      if(tn+fn>0)summary.println("npv "+fn/ (double)(tn+fn));
   }

   if(Constants.sampleWithPedigree()){
      sampler.data.extractFromTrioData();
   }
 
 if(Constants.fillGaps() && Constants.sample()){
     Set<Integer>gapSites = new HashSet<Integer>();
     Map<Integer, Integer> polySites = new TreeMap<Integer, Integer>();
     Map<Integer, Map<String, Double>> polyCertainty = new TreeMap<Integer, Map<String, Double>>();
     Map<Integer, Integer> realPolySites = new TreeMap<Integer, Integer>();
     int gap_count = 0;
     int poly_count=0;
    for(int i=0; i<original.length(); i++){
        for(Iterator<String> it = original.getKeys().iterator(); it.hasNext();){
            String key = it.next(); 
            if(((ComparableArray)sampler.data.get(key).getElement(i)).containsNull()){
                Integer val = polySites.get(i);
              polySites.put(i, val==null ? 1 : val+1);
               poly_count++;
             
            }
            if(((ComparableArray)original.get(key).getElement(i)).containsNull()){
                Integer val = realPolySites.get(i);
                realPolySites.put(i, val==null ? 1 : val+1);
            }
            if(((ComparableArray)affy_reference.get(key).getElement(i)).containsNull()){
                gapSites.add(i);
                gap_count++;
            }
        }
    }
    summary.println("frac poly " +((double)polySites.size()/(double)gapSites.size())+" "+ (double) poly_count/ (double)gap_count );
    summary.println(polySites);
    summary.println(gapSites);
    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "gap_poly_sites"+".txt"))));
    pw.println("total "+sampler.data.size()+" "+cal.getTime().toString());
    for(Iterator<Entry<Integer, Integer>> it = polySites.entrySet().iterator(); it.hasNext();){
        Entry<Integer, Integer> ent = it.next();
        pw.print(ent);
        pw.print("\t");
        pw.print(sampler.data.snpid.get(ent.getKey()));
        pw.print("\t");
        pw.println(polyCertainty.get(ent.getKey()));
    }
    pw.close();
    
    if(sampler.data.ped()!=null){
        //sampler.data.arrangeDataAccordingToPedigree();
       double[][] cons =  sampler.data.checkConsistency();
       summary.println("Consistency inferred");
       summary.println("Consistency emissions: "+cons[0][0]+" "+cons[0][1]);
       summary.println("Consistency state:     "+cons[1][0]+" "+cons[1][1]);
       summary.println("Consistency original");
        cons =  original.checkConsistency();
       summary.println("Consistency emissions: "+cons[0][0]+" "+cons[0][1]);
       summary.println("Consistency state:     "+cons[1][0]+" "+cons[1][1]);
    }
 
 }

   int numIt = Constants.numItSum();
   int[][] sources = original.getSourcePositions();
       if(original!=affy_reference ){
           if(numIt>0){
               PrintWriter logCNV = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "logCNV"+".txt"))));
               PrintWriter logGeno = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "logGeno"+".txt"))));
               PrintWriter logHaplo = new PrintWriter(new BufferedWriter(new FileWriter(new File(parentfile, "logHaplo"+".txt"))));
               summary.println("COMPARING INFERENCE WITH AFFY");
               
               compare(sampler.data, affy_reference,
                       this.unc,
                       "sampling avg",   
                       0.0, summary, 
                      
                        sources
                       , new PrintWriter[] {logCNV, logGeno, logHaplo}
                      );   
               logCNV.close();
               logGeno.close();
               logHaplo.close();
           }
          
           summary.println("COMPARING ORIGINAL WITH AFFY");
           compare(original,affy_reference,unc,
                   "sampling avg",   
                    0.0, summary,
                 sources
                   , new PrintWriter[] {null, null, null}
                 );   
          
       
       }
       if(numIt >0){
           summary.println("COMPARING INFERENCE WITH ORIGINAL");
               compare(sampler.data, original,unc,
                       "sampling avg",   
                        0.0, summary,
                      sources
                        , new PrintWriter[] {null, null, null}
                      );
           }
   
  
   }catch(Exception exc){
    exc.printStackTrace();
}
}*/
