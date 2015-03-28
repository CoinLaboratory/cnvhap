package lc1.possel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

import lc1.CGH.AberationFinder;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.EmissionState;
import lc1.stats.TrainableNormal;
/** worry about edges!!!!!! */

public class RunEHH {
    static Logger logger = Logger.getLogger("RunEHH");
    static{
        try{
        ConsoleHandler handler = new ConsoleHandler();
        FileHandler handlerF = new FileHandler("logres.txt", false);
       logger.addHandler(handlerF);
        Formatter formatter = 
        new Formatter(){
            public String format(LogRecord record){
                return record.getMessage()+"\n";
            }
        };
        handler.setFormatter(formatter);
        handlerF.setFormatter(formatter);
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
        }
    }
    
    PrintWriter score;
    final Map<String, Emiss> ancestralMap = new HashMap<String, Emiss>();
    //final File file;
   // final File dataDir;
    final File user;
    final  Map<Double, TrainableNormal> normal;
     String chrom;
  //  Map<String, Location> locsCovered = new HashMap<String, Location>();
    
   final boolean scoring;
    RunEHH(File user,   Map<Double, TrainableNormal> normal, PrintWriter score, boolean scoring, Class clazz) throws Exception{
       this.score = score;
       this.clazz = clazz;
       this.scoring = scoring;
     //  this.dataDir = dataDir;
        this.normal =  normal;
        this.user = user;
    }
    public void setChrom(String chrom) throws Exception{
        ancestralMap.clear();
        this.chrom = chrom.split("_")[0];
        File ancestralFile  = new File(user, "data/"+this.chrom+"_anc.txt");
        System.err.println("ancestral file is "+ancestralFile);
      //  if(ancestralFile.exists()){
                BufferedReader br = read(ancestralFile);
                String st = "";
                ancestralMap.clear();
                while((st = br.readLine())!=null){
                   String[] str = st.trim().split("\\s+");
                  Emiss em;
                /*  if(str[3].equals("A")) em = Emiss.A;
                  else if(str[3].equals("B")) em = Emiss.B;
                  else throw new RuntimeException("");
                   ancestralMap.put(str[0], em );*/
               }
               br.close();
    }
  public static BufferedReader read(File f) throws Exception{
	  if(f.exists()){
		  return new BufferedReader(new FileReader(f));
	  }
	  else{
		  File f1 = new File(f.getParentFile(), f.getName()+".gz");
          if(!f1.exists()) return null;
		  return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f1))));
	  }
  }
   
    
   
   
   
    /** returns [derived, ancestral]*/
    public Double getHaps(EmissionState maf, int i, String snp_id, List<String> haps, boolean needAncestral){
        if(maf.getFixedInteger(i)!=null){
            return null;
        }
        Emiss ancestral  = Emiss.a();
        Emiss derived = Emiss.b();
        if(needAncestral){
           ancestral = this.ancestralMap.get(snp_id);
            if(ancestral==null) return null;
            derived = ancestral==Emiss.a() ? Emiss.b() : Emiss.a();
          
        }
        haps.add(derived.toString());
        haps.add(ancestral.toString());
        int derived_index = maf.getEmissionStateSpace().get(derived);
        int ancestral_index = maf.getEmissionStateSpace().get(ancestral);
        double[] emiss = maf.getEmiss(i);
        if(emiss[derived_index]==0 || emiss[ancestral_index]==0) return null;
        return emiss[derived_index];
    }
    public Double getCNVIndices(EmissionState maf, int i, List<String> haps, int copyNumber){
        if(maf.getFixedInteger(i)!=null) return null;
        String cnv_hap = "";
        EmissionStateSpace emStSp = maf.getEmissionStateSpace();
        double[] emiss = maf.getEmiss(i);
        double max_prob=0;
        for(int k=0;k< emiss.length; k++ ){
            if(emiss[k] >0){
                if(emStSp.getCN(k)==copyNumber && emiss[k] > max_prob){
                    max_prob = emiss[k];
                    cnv_hap = emStSp.get(k).toString();
                    
                }
                else if(emStSp.getCN(k)==1)haps.add(emStSp.get(k).toString());
            }
        }
        if(cnv_hap.length()==0) return null;
        else {
            haps.add(0, cnv_hap);
            return max_prob;
        }
    }
   
    
    
    List<Double> leftAnc = new ArrayList<Double>();
    List<Double> leftDerived = new ArrayList<Double>();
    List<Double> rightDerived = new ArrayList<Double>();
    List<Double> rightAnc = new ArrayList<Double>();
    List<Integer> loc = new ArrayList<Integer>();
    List<String> id = new ArrayList<String>();
    List<String> idR = new ArrayList<String>();
    List<Double> derivedFreq = new ArrayList<Double>();
   static double min_p = 0.1;
   static  double max_sc = 1.0;
  //  int lastLoc=0;
 //   public static plot = true;
    public  void run(File f,  int nextStart, int copyNumber) throws Exception{
//        String[] str = f.getName().split("\\.")[0].split("_");
        SimpleDataCollection sdt = 
            SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader(f, "phased1.txt"),Emiss.class, Emiss.getEmissionStateSpace(1));
        File snpFile  = new File(user, "snp.txt");
                sdt.snpid= new ArrayList<String>();
                sdt.readPosInfo(AberationFinder.getBufferedReader(f, "snp.txt"), new int[] {0}, true, new List[] {sdt.snpid}, new Class[] {String.class});
      //  sdt.calculateMaf(false);
        sdt.split();
        AbstractEHH ehh =(AbstractEHH) this.clazz.getConstructor(new Class[] {sdt.getClass(), String.class})
        .newInstance(new Object[] {sdt, f.getName()} );
        List<String> haps = new ArrayList<String>();
        boolean[] doLR = new boolean[2];
        //boolean graphPrev = false;
        for(int i=0; i<sdt.length(); i++){
           /* if(f.s!=null && !f.contains(sdt.loc, i)){
                if(sdt.loc.get(i)==17927329){
                    Logger.global.info("h");
                }
                continue;
            }*/
            try{
            //System.err.print(i+" ");
            int pos = sdt.loc.get(i);
            haps.clear();
            Double[][] scores = new Double[2][2];
            List<Double>[][] in = new List[2][2];
            boolean doLeft = loc.size()==0 || pos > loc.get(loc.size()-1);
            boolean doRight = pos <nextStart;
            doLR[0] = doLeft; doLR[1] = doRight;
            Double maf = null;
       //         copyNumber==1 ?
          //      getHaps(sdt.maf, i, sdt.snpid.get(i), haps, this.clazz.equals(EHH.class)):
            //        getCNVIndices(sdt.maf, i, haps, copyNumber);
            if(maf!=null){ //ie not fixed at this position
                if(haps.size()>2){
                    System.err.println("haps are "+haps);
                  /*  String[] ancs = haps.subList(1, haps.size()).toArray(new String[0]);
                    Integer index =  ehh.findAncestral(i,i, haps.get(0), ancs, doLR);
                    if(index!=null){
                        ehh.setCore(i,i, haps.get(0), ancs[index], scores, doLR);
                    }*/
                }
                else {
                    doLR[0] = true;
                    doLR[1] =true;
                   if(true) throw new RuntimeException("!!");
                 //   ehh.setCore(i,i, haps.get(0), haps.get(1), scores, doLR, in);
                    if(copyNumber!=1 ||!EHHMultiChromosome.overl)  {
                        try{
                            if(this.clazz==EHH.class && scores[0][0]!=null && scores[1][0]!=null && in[0]!=null){
                                  double sc =   Math.log((scores[0][0]+scores[1][0])/(scores[0][1]+scores[1][1]));
                                  double pval =this.pvalue(sc, maf);
                                  if( (( pval < 1e-3  && sc>1.0 )|| !EHHMultiChromosome.overl)){
                                      //min_p = pval;
                                      //max_sc = sc;
                                     // GraphIHH.plot(in, sc,pval, i, i, sdt.loc, f, maf, 1.0-maf, true);
                                    //  graphPrev = true;
                                      
                                  }
                            }
                        }catch(Exception exc){
                            exc.printStackTrace();
                        }
                    }
                  // else graphPrev = false;
                }
            }
            
                 if(doLeft){
                     leftDerived.add(scores[0][0]);
                     leftAnc.add(scores[0][1]);
                     this.loc.add(pos);
                     this.id.add(f.getName().split("\\.")[0]);
                   
                     derivedFreq.add(maf);
                 }
                 if(doRight){
                     if(pos!=loc.get(rightDerived.size())) throw new RuntimeException("!!");
                     rightDerived.add(scores[1][0]);
                     rightAnc.add(scores[1][1]);
                     this.idR.add(f.getName().split("\\.")[0]);
                 }
                 if(doRight && doLeft){
                     if(rightDerived.size()!= leftDerived.size()) throw new RuntimeException("!!");
                 }
            }catch(Exception exc){
                System.err.println("problem at "+i+" "+f);
                exc.printStackTrace();
            }
        }
    }
    
    final Class clazz;
    public  void run(File f,  int start,  String hapString, String anc) throws Exception{
//      String[] str = f.getName().split("\\.")[0].split("_");
      SimpleDataCollection sdt = 
          SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader(f, "phased1.txt"),Emiss.class, Emiss.getEmissionStateSpace(1));
      File snpFile  = new File(user, "snp.txt");
              sdt.snpid= new ArrayList<String>();
              sdt.readPosInfo(AberationFinder.getBufferedReader(f, "snp.txt"), new int[] {0}, true, new List[] {sdt.snpid}, new Class[] {String.class});
     // sdt.calculateMaf(false);
      sdt.split();
      AbstractEHH ehh = (AbstractEHH) clazz.getConstructor(new Class[] {sdt.getClass(), String.class}).newInstance(new Object[] {sdt, f.getName()}); 
          //new EHH(sdt, f.getName() );
      List<String> haps = new ArrayList<String>();
      boolean[] doLR = new boolean[]{true, true};
      int i1 = sdt.loc.indexOf(start);
      int i2 = i1+hapString.length()-1;
          //System.err.print(i+" ");
       Map<String, List<PIGData>> m = sdt.getAllHaplotypes(i1, i2);
       List l = m.remove(hapString);
       if(l==null) throw new RuntimeException("!!");
    //   String anc;
       /*if(m.keySet().size()>1) {
    	   String[] ancestral = m.keySet().toArray(new String[0]);
    	   anc = ancestral[ehh.findAncestral(i1, i2, hapString, ancestral, doLR)];
       }
       else anc = m.keySet().iterator().next();*/
     // anc = "AAA";
       Double maf = (double)l.size() / ((double)m.get(anc).size()+ (double)l.size() );
          haps.clear();
          Double[][] scores = new Double[2][2];
          List<Double>[][] in = new List[2][2];
              if(haps.size()>2){
                /*  String[] ancs = haps.subList(1, haps.size()).toArray(new String[0]);
                  Integer index =  ehh.findAncestral(i,i, haps.get(0), ancs, doLR);
                  if(index!=null){
                      ehh.setCore(i,i, haps.get(0), ancs[index], scores, doLR);
                  }*/
              }
              else {
                  doLR[0] = true;
                  doLR[1] =true;
                  if(true) throw new RuntimeException("!!");
                  //ehh.setCore(i1,i2, hapString, anc, scores, doLR, in);
                 // if(copyNumber!=1 ||!EHHMultiChromosome.overl)  {
                      try{
                          if(scores[0][0]!=null && scores[1][0]!=null){
                                double sc =   ehh.score(scores);//
                             double pval =this.pvalue(sc, maf);
                           double sc1 = this.scoreNorm(sc, maf);
                             //   if( (( pval < 1e-3 || sc > max_sc  )|| !EHHMultiChromosome.overl)){
                              //      min_p = pval;
                               //     max_sc = sc;
                                  //  GraphIHH.plot(in, sc1,pval, i1, i2, sdt.loc, f, maf, 1.0-maf, true);
                              //  }
                          }
                      }catch(Exception exc){
                          exc.printStackTrace();
                      }
                 // }
              }
          
              
         
  }
    
    public double scoreNorm(double ratio, double maf1){

    	double maf = round(maf1);

    	if(Double.isNaN(ratio) || ratio==Double.POSITIVE_INFINITY || ratio ==Double.NEGATIVE_INFINITY

    	|| Double.isInfinite(ratio)) throw new RuntimeException(" !! "+ratio);

    	TrainableNormal norm = normal.get(maf);

    	double sc = (ratio- norm.getMean())/norm.getStdDev();

    	return sc;

    	}


    	 

    	 
    	public void calcScores(){
	    	if(rightDerived.size()!= leftDerived.size()) throw new RuntimeException("!! "+rightDerived.size()+" "+leftDerived.size());
	    	Object[] toPrint = new Object[8];
	    	String fs = "%-5s\t%15s\t%15s\t%10i\t%5.3g\t%5.3g\t%5.3g\t%5.3g";
	    	toPrint[0] =this.chrom;
	    	for(int i=0; i<this.loc.size(); i++){
	    	try{
	    	toPrint[1] = this.id.get(i);
	    	toPrint[2] = this.idR.get(i);
	    	if(toPrint[2].equals(toPrint[1])) toPrint[2] = "-";
	    	toPrint[3] = this.loc.get(i);
	    	Double derL = this.leftDerived.get(i);
	    	Double derR = this.rightDerived.get(i);
	    	Double ancL = this.leftAnc.get(i);
	    	Double ancR = this.rightAnc.get(i);
	    	if(derL!=null && derR !=null && ancL!=null && ancR!=null){
	    	System.err.print(i+" ");
	    	double sc = this.clazz==EHH.class ? Math.log((derL+derR)/(ancL+ancR)) : Math.max(derL,derR);
	    	this.addCount(sc, derivedFreq.get(i));
	    	double pval = this.pvalue(sc, derivedFreq.get(i));
	    	toPrint[4] = derivedFreq.get(i);
	    	toPrint[5] = sc;
	    	toPrint[6] = Math.min(pval, 1-pval);
	    	toPrint[7] = this.scoreNorm(sc, derivedFreq.get(i));
	    	this.score.println(String.format(fs, toPrint));
            //System.err.println(Format.sprintf(fs, toPrint));
	     }

    	}catch(Exception exc){

    	exc.printStackTrace();

    	}

    	}

    	System.err.println();

    	this.score.flush();

    	}
  
    public double round(double maf1){
//        return maf1;
      return Math.round(maf1*20.0)/20.0;
    }
    public void addCount(double ratio, double maf1){
       double  maf = round(maf1);
        if(!Double.isNaN(ratio) &&  !Double.isInfinite(ratio) 
                && ratio!=Double.NEGATIVE_INFINITY
                && ratio!=Double.POSITIVE_INFINITY ) {
           TrainableNormal norm =  normal.get(maf);
           if(norm==null){
               if(scoring){}// throw new RuntimeException(" nothing found for "+maf);
               normal.put(maf, norm = new TrainableNormal(null,0,1, 100,1.0));
           }
            norm.addCount(ratio, 1.0);
        }
        else{
            throw new RuntimeException("!! "+ratio);
        }
            
        }
    public double pvalue(double ratio, double maf1){
       double  maf = round(maf1);
        if(Double.isNaN(ratio) 
                || ratio==Double.POSITIVE_INFINITY 
                || ratio ==Double.NEGATIVE_INFINITY
                || Double.isInfinite(ratio)) throw new RuntimeException("   !!    "+ratio);
            TrainableNormal norm =  normal.get(maf);
            double sc =  (1.0- norm.cumulative(ratio));
            if(sc < 1e-4){
                System.err.println(maf1+"->"+maf+" low pvalue for "+ratio+" "+norm+" pval : "+sc);
            }
            return sc;
          //  double sc = Math.abs(Math.log(ratio));
         //   score.println(file.getName().toString()+" "+start+" "+end+" "+derived+" "+ancestral+" "+sc+" "+norm+" "
         //           +));
          //  score.flush();
          //  log.flush();
       
      //  System.err.println(start+" ihh "+ratio+" for maf "+coreInner.maf+" "+core+" "+core1);
        
    }
}
