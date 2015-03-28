/**
 * 
 */
package lc1.possel;

import hep.aida.ref.Histogram1D;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.EmissionState;
import lc1.util.Constants;
import calc.DPrimeCalculator;
import calc.LDCalculator;
import calc.R2Calculator;

class LDInner extends AbstractEHH{
    static int numPerm =0;
    static boolean show = false;
   // SimpleDataCollection sdt;
    final PIGData data;
    final PIGData immutable;
   // int start, end;
   // double[] res;
    final LDCalculator lower;
    double maf;
    final String[] keys;
    List<String> keyL=  new ArrayList<String>();
    Comparable gap = Emiss.N();
    Comparable noGap = Emiss.a();
  //  int[] occ0;
    final int length;
    int exp;
   final  Double[]  randomSc ;
    double pval;
    double eval;
    public void setCore(int start, int end, String derived, String ancestral, Double[][]res, boolean[] doLR, List<Double>[][] in1){
        super.setCore(start, end, new String[] {derived, ancestral}, res, doLR, in1);
        EmissionStateSpace emStSp = maf1.getEmissionStateSpace();
        maf = maf1.score(emStSp.get(Emiss.N()), start);
        exp = (int) Math.round(maf *sdt.getKeys().size());
        this.start = start;
        this.end = end;
       Arrays.fill(randomSc, null);
     //   Arrays.fill(this.res, 0.0);
        for(int i=start; i<=end; i++){
                data.set(i,ComparableArray.make(derived.charAt(i-start), ancestral.charAt(i-start)));
            //data2.addPoint(ComparableArray.make(Emiss.a(), Emiss.b()));
        }
        this.calcLD();
        for(int i=0; i<doLR.length; i++){
            Arrays.fill(res[i], null);
            Arrays.fill(in1[i], null);
            if(!doLR[i]) continue;
            List<Double > l =i==0 ? left : right;
            res[i][0] = l.get(max(l));
            res[i][1] = res[i][0];
            in1[i][0] =l;
            in1[i][0] = l;
        }
        for(int i=start; i<=end; i++){
            data.set(i,immutable.getElement(i));
        //data2.addPoint(ComparableArray.make(Emiss.a(), Emiss.b()));
    }
      //  Logger.global.info("h");
        
    }
    
    static double lenthresh = 50000;
    
   
    private int max(List<Double> left) {
       int j=0;
        for(int i=1; i<left.size(); i++){
            if(left.get(i) > left.get(j)) j=i;
        }
       
      return j;
    }
    public static PIGData makeBase(int len){
        PIGData data = SimpleScorableObject.make("", len,     Emiss.getEmissionStateSpace(1),(short)-1);
        for(int i=0; i<len; i++){
                data.addPoint(i, ComparableArray.make(Emiss.a(), Emiss.b()));
        }
        return data;
    }
    public LDInner(SimpleDataCollection sdt, String name, boolean dprime) throws Exception {
        super(sdt, name);
       // this.sdt = sdt;
        this.randomSc = new Double[numPerm];
        keys = sdt.getKeys().toArray(new String[0]);
        length = keys.length;
         maf1 = null;//sdt.maf;
        EmissionStateSpace emStSp = maf1.getEmissionStateSpace();
        data = makeBase(sdt.length());
      
           maf = maf1.score(emStSp.get(Emiss.N()), start);
           exp = (int) Math.round(maf *sdt.getKeys().size());
        
          lower =dprime ? new DPrimeCalculator():new R2Calculator();
          immutable = SimpleScorableObject.make(data);
    }
    final EmissionState maf1;
    public LDInner(SimpleDataCollection sdt2, boolean b, int numPerm2, boolean dprime) throws Exception{
        super(sdt2, "");
        keys = sdt.getKeys().toArray(new String[0]);
        length = keys.length;
        data = makeBase(sdt.length());
        immutable =(PIGData)  data.clone();
        lower = dprime ? new DPrimeCalculator() : new R2Calculator();
        randomSc = new Double[numPerm2];
    //    sdt.calculateMaf();
         maf1 =null;// sdt.maf;
    
      // this.sdt = sdt;
    }
    private void permute(){
        this.keyL.clear();
        this.keyL.addAll(Arrays.asList(keys));
     
        for(int i=0; i<exp; i++){
            String key = (keyL.remove(Constants.nextInt(keyL.size())));
            PIGData data = sdt.get(key);
            for(int k=start; k<=end; k++){
                data.set(k, gap);
            }
        }
        for(int i=0; i<keyL.size(); i++){
            String key = keyL.get(i);
            PIGData data = sdt.get(key);
            for(int k=start; k<=end; k++){
                data.set(k, noGap);
            }
        }
    }
    public double score(Double[][] scores) {
        return scores[0][0]+scores[1][0];
    } 
    public void getRandomScores(double sc){
        for(int i=0; i<randomSc.length; i++){
            this.permute();
           randomSc[i] = calcLD();
        }
     //   EVD evd = null;
        Histogram1D h = null;
        try{
            Arrays.sort(randomSc);
            int min = (int) Math.floor(randomSc[0].doubleValue());
            int max = (int) Math.ceil(randomSc[randomSc.length-1].doubleValue());
            h =   new Histogram1D("hist", 1, min, max);
      for(int k=0; k<randomSc.length; k++){
    	 h.fill(randomSc[k]);
      }
           
           //h = EVD.makeH(randomSc);
        
      //  log.println(h);
        }catch(Exception exc){}
        if(true || randomSc.length>=120){
        try{
       //    if(h!=null) evd= new EVD(h);
           
            
        }catch(Exception exc){
           // exc.printStackTrace();
        }
        int j=0;
        int len1 = randomSc.length;
        double min=  1 / (double)randomSc.length;
        for(; j<len1; j++){
            double rsc = randomSc[len1-j-1];
            if(rsc < sc){
               
                break;
            }
        }
       
           pval = Math.max(min, (double) (j) / (double)len1);
       //    eval = evd==null ? 1.0 : evd.extremeValueP(sc);
          // Logger.global.info("h");
        }
        
    }
    public double[] run(){
        double sc = this.calcLD();
        this.getRandomScores(sc);
        return new double[] {sc, pval, eval};
    }
    List<Double> left ;
    List<Double> right;
    public double calcLD(){
        left = sdt.calcLD(lower, data,  true, start,(int)lenthresh);
        right = sdt.calcLD(lower, data,  false, this.end,(int)lenthresh);
        int maxL = max(left);
        int maxR = max(right);
        if(left.get(maxL) > right.get(maxR)){
            Logger.global.info("best is "+sdt.snpid.get(start - maxL));
            return left.get(maxL);
        }
        else{
            Logger.global.info("best is "+sdt.snpid.get(end + maxR));
            return right.get(maxR);
        }
    }
    public void setCore(int st1, int end1) {
       // this.data.set(st1, immutable.getElement(i));
        for(int i=st1; i<=end1; i++){
            this.data.set(i, immutable.getElement(i));
        }
    //   super.setCore(st1, end1, "A", "_", new Double[2][2], new boolean[] {true, true}, 
    //           new List[2][2]);
       Comparable comp = ComparableArray.make(Emiss.A, Emiss.N);
       for(int i=st1; i<=end1; i++){
           this.data.set(i, comp);
       }
       EmissionStateSpace emStSp = maf1.getEmissionStateSpace();
       maf = maf1.score(emStSp.get(Emiss.N()), start);
       exp = (int) Math.round(maf *sdt.getKeys().size());
     System.err.println("data is "+data);
        
    }
}