package lc1.dp.data.collection;

import java.util.ArrayList;
import java.util.List;

public abstract class DataC implements DataCInterface {

   public String name = "-";
    public Integer length;
    public List<String> snpid = new ArrayList<String>();
    public List<Integer> loc = new ArrayList<Integer>();
    public List<Character> alleleA = new ArrayList<Character>();
    public List<Double> baf = new ArrayList<Double>();
    public List<Double> minR = new ArrayList<Double>();
    public List<Double> maxR = new ArrayList<Double>();
    public List<Character> alleleB= new ArrayList<Character>();
    public List<Boolean> strand = new ArrayList<Boolean>();
    //  public List<String> names  = new ArrayList<String>();
     // public EmissionState maf;
     // public Map<String, Integer> name_index = new HashMap<String, Integer>();
    
   /* public abstract int noCopies(String i);

  //  public abstract void fillLikelihoodData(double pseudo, int[] indices);

    public abstract void writeSNPFile(File file, String chr, boolean header,
            Collection<Integer> toD) throws Exception;

    /** i is position, k is value */
    public  Double scoreChi(int i){
        return 1.0;
    }

    /** i is position, k is value */
    public  void scoreChi1(int i,  boolean useEmSt, int phenIndex){
        throw new RuntimeException("!!");
        //return new Double[] {1.0, 1.0, 1.0};
    }

 //   public abstract int count(int[] indices, int index);

    //public abstract Double[] scoreChi(double[][] ns, double[][] nc);



//    public abstract Object getInfo(String tag, String key, int i, boolean style)
  //          throws Exception;

    //public abstract void reverse();

    //public abstract int[][] getSourcePositions();

    //public abstract void mix();

    //public abstract EmissionState makeMafState(EmissionStateSpace emStSp1);

    //public abstract DataC clone();

    //public abstract int noAllelles();

    //public abstract void extractFromTrioData();

//    public abstract void calculateMaf(boolean state);

  //  public abstract EmissionState calculateMaf1();

    //public abstract void printTrioData(PrintWriter pw);

//    public abstract PIGData get(String st);

  //  public abstract int countSwitches();

    /*  public List<Clusters> getClusters(int numF){
         List<Integer> bound = getBlockBoundaries();
         List<Clusters> result = new ArrayList<Clusters>();
         for(int ik=0; ik<bound.size(); ik++){
             Clusters res = new Clusters(numF, bound.get(ik), ik==bound.size()-1 ? length() : bound.get(ik+1));
             int i = bound.get(ik);
             for(Iterator<PIGData> it = this.viterbi.values().iterator(); it.hasNext();){
                 PIGData nxt = it.next();
                 ComparableArray dat= (ComparableArray) nxt.getElement(i);
               Integer[] val = new Integer[dat.size()];
                 for(int j=0; j<dat.size(); j++){
                     val[j] = (Integer)dat.get(j);
                  
                  
                 }
                 res.add(val, nxt.getName());
             }
             res.sort();
             result.add(res);
         }
         return result;
     }*/
//    public abstract List<Integer> getBlockBoundaries();

  //  public abstract void getBlockBoundaries(String key, Set<Integer> res);

    //public abstract void restricToAlias(Collection<String> alias);

    /** only for half trio emissions states 
         public void getRecSitesFromViterbi(){
          //   if(!((ComparableArray)this.data.get(0).getElement(0)).isNested()) return;
             recSites.clear();
             Integer zero =0;
             Integer one = 1;
             for(Iterator<String> it = viterbi.keySet().iterator(); it.hasNext();){
                 String key = it.next();
                 PIGData dat_i = viterbi.get(key);
                 if((((ComparableArray)dat_i.getElement(0))).size()==1){
                     recSites.put(dat_i.getName(), null);
                   
                 }
                 else{
                     
                     PIGData recS =  new PIGData(key, this.length());
                     int noC = ((ComparableArray)dat_i.getElement(0)).size();
                         for(int i=0; i<this.length(); i++){
                             TrComparableArray tic = (TrComparableArray)dat_i.getElement(i);
                            recS.addPoint(tic.third);
                         }
                         recSites.put(dat_i.getName(), recS);
                 }
             }
             Logger.global.info(""+this.recSites.size());
         }*/

    /*  public void restrictSites(int i) {
          if(i<this.length()){
              
              if(loc.size()>0) loc = loc.subList(0,i);
              if(names.size()>0) names = names.subList(0,i);
              if(maf!=null) maf.restrict(i);
              if(this.homo_count.size()>0) homo_count = homo_count.subList(0, i);
          }
              
      }*/

   // public abstract EmissionState getL(String st);

    //public abstract int getFirstIndexAbove(long pos);

   // public abstract int getIndex(Integer name);

//    public abstract int length();

  //  public abstract void printHapMapFormat(File f, String chr) throws Exception;

    //public abstract void printHapMapFormat(PrintWriter pw,
      //      Collection<Integer> toD, boolean style) throws Exception;

    //public abstract void printHapMapFormat(PrintWriter pw,List<String> indiv, 
      //      Collection<Integer> toD, boolean style, String[] initTags,String[] phenTags,
        //    String[] tags, String len
    //            Map<String, PIGData> data, 
    //          Map<String, Double[]> uncertainty,
    //        Map<String, Double[]> uncertaintyPhase,
    //      Map<String, EmissionState> emState
    //) throws Exception;

    //public abstract int size();

    //public abstract void writeDickFormat(File f, boolean header,
      //      Collection<String> idv, Collection<Integer> toD) throws Exception;

    //public abstract void getDeletedPositions(double mafThresh,
      //      SortedSet<Integer> ampl, SortedSet<Integer> del);

//    public abstract void printLocations(PrintWriter pw);

  //  public abstract List<Aberation> getDeletedPositions(boolean deletion);

//    public abstract void printDeletedPositions(File parentFile);

  //  public abstract void writeLocation(PrintWriter pw,
    //        Collection<Integer> toDrop);

    //public abstract void writeFastphase(
      //      Map<String, PIGData> data1,
         
        //    Map<String, double[]> uncertaintyPhase, PrintWriter pw,
            //    boolean restrictToPairs, 
          //  boolean printUncertainty, boolean mark, Collection<Integer> toDrop)
            //throws Exception;

//    public abstract void writeFastphase(
  //          Map<String, PIGData> data1,
          
    //        Map<String, double[]> uncertaintyPhase, PrintWriter pw,
            //    boolean restrictToPairs, 
      //    boolean printUncertainty, boolean mark, boolean expand,
       //     Collection<Integer> toDrop) throws Exception;

   // public abstract List<String> getNames();

    //public abstract void trim(int i);

    //public abstract List<Integer> getLocations();

    //public abstract void replace(Map<String, PIGData> newData);

  //  public abstract void arrangeDataAccordingToPedigree();

//    public abstract void replaceL(Map<String, EmissionState> newData);

   // public abstract Map<String, PIGData> arrangeDataAccordingToPedigree(
     //       Map<String, PIGData> datai);

  //  public abstract Iterator<PIGData> iterator();

//    public abstract void set(int ij, PIGData newDat);

  //  public abstract void setPedigree(PedigreeDataCollection pedData);

    //public abstract Set<String> getKeys();

   // public abstract double[][] checkConsistency();

   // public abstract double[] checkConsistency(
     //       Map<String, PIGData> datai);

   // public abstract boolean containsKey(String name2);

   // public abstract void summarise();

//    public abstract double[] calcLDAverage();

//    public abstract double[][] calcLD();

  //  public abstract double[][] calcLD(String st1, String st2);

    /** splits the data into haplotype data 
     * Note does not split likelihood data!!!
     * */
 //   public abstract void split();

   // public abstract void calcLD(LDCalculator upper, LDCalculator lower,
     //       PIGData poss, double[][] res);

    //public abstract List<Double> calcLD(LDCalculator upper,
      //      PIGData poss, boolean left, int pos, int lenthresh);

    /*
     * st1 are the different possible haplotypes at pos1
     * st2 are the different possible haplotypes at pos2
     */
    //public abstract double[][] getLD(int pos1, int pos2, String[] st1,
      //      String[] st2);

   // public abstract int countHaplotypes(int[] st, int[] en, String[] str);

   // public abstract void printIndiv(PrintWriter pw_indiv);



//    public abstract double[] uncertaintyPhase(String key);

  //  public abstract double[] uncertaintyVitPhase(String key);

   

  //  public abstract void setViterbi(String key, PIGData datvit);

   // public abstract void setRecSites(String name,
     //       SortedMap<Integer, Integer>[] sampleRecSites);

//    public abstract void clearViterbi();

  //  public abstract Iterator<EmissionState> dataLvalues();

    //public abstract void clearData();

  //  public abstract SortedMap<Integer, Integer>[] recSites(String j);

  //  public abstract void writeFastphase(File pw, boolean states) throws Exception;

   

   // public abstract EmissionState maf();

  //  public abstract Object ped();
//
  //  public abstract EmissionState dataL(String string) ;

    //public abstract EmissionState  getState(String key, EmissionStateSpace stateStSp) ;

 /*   public  double weight(String key) {
        double[] weights = Constants.weights();
     List<Short> dataIndices = ((HaplotypeEmissionState)dataL(key)).getDataIndices();
     double max = 0;
     for(int i=0; i<dataIndices.size(); i++){
         double wt = weights[dataIndices.get(i)];
         if(wt>max) max = wt;
     }
     return max;
    }*/

  //  public abstract void writeDickFormat1(File file, boolean b) throws Exception;

    public String[] getUnderlyingDataSets() {
       return new String[] {this.name};
    }

   // public abstract ProbabilityDistribution[] numLevels() ;

   // public abstract Phenotypes pheno();

	



}