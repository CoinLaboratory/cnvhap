package lc1.possel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.util.Constants;

public class EHH extends AbstractEHH{
  
    
    public static void main(String[] args){
        try{
            File user = new File(System.getProperty("user.dir"));
        File[] f = user.listFiles(new FilenameFilter(){
            public boolean accept(File arg0, String arg1) {
               return arg1.endsWith(".zip");
            
            }
           });
        for(int i=0; i<f.length; i++){
            ZipFile zf = new ZipFile(f[i]);
            BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("phased2.txt_1"))));
            SimpleDataCollection sdt = 
                SimpleDataCollection.readFastPhaseOutput(br,Emiss.class, Emiss.getEmissionStateSpace(1));
           int st = 10;
           Map<String, List<PIGData>> m =sdt.getAllHaplotypes(st, st);
           String[] cores =  m.keySet().toArray(new String[0]);
           double[] maf = new double[cores.length];
           double sum =0;
           for(int ik=0; ik<maf.length; ik++){
               maf[ik] = m.get(cores[ik]).size();
               sum+=maf[ik];
           }
           for(int ik=0; ik<maf.length; ik++){
               maf[ik] = maf[ik]/sum;
           }
           
            EHH ehh = new EHH(sdt,f[i].getName());
            Double[][]res = new Double[2][cores.length];
            List<Double>[][] in1 = new List[2][cores.length];
            for(int k=0; k<in1.length; k++){
                for(int j=0; j<in1.length; j++){
                    in1[k][j] = new ArrayList<Double>();
                }
            }
            
            ehh.setCore(st, st, cores, res, new boolean[] {true, true}, in1);
            GraphIHH.plot(in1, 0, 0, st, st, sdt.loc, f[i], maf[0],  true);
            Logger.global.info("h");
        
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    static double thresh = 0.05;  //below which we do not extend any further
    
   
    Map<String, EHHCalculator> map = new HashMap<String, EHHCalculator>();

    public EHH(DataCollection coll,String f) throws Exception{
      super(coll,f);
    }
  
    public int noCharInCommon(String st, String st1, boolean left){
       int len = st.length();
       if(len!=st1.length()) throw new RuntimeException("!!");
    
       for(int i=0; i<len; i++){
           if(!left){
               if(st.charAt(i)!=st1.charAt(i)) return i;
           }
           else if(left){
               if(st.charAt(len -(i+1) ) != st1.charAt(len-(i+1))){
                   return i;
               }
           }
       }
       return len;
   }
   
    /** this class calculates the EHH probability */
    class EHHCalculator{
        boolean left; // 
        List<Double> scores = new ArrayList<Double> ();
        List<Integer> loc1 = new ArrayList<Integer>();
        SortedMap<String, List<PIGData>> hapl_1;
        int originalNumber;
        String core;
        int pos =0;
        StringBuffer ch = new StringBuffer();
        String mostCommon(){
            Iterator<Entry<String, List<PIGData>>> it =hapl_1.entrySet().iterator();
            String st = "";int count=0;
            while(it.hasNext()){
                Entry<String, List<PIGData>> ent = it.next();
                int cnt1 = ent.getValue().size();
                if(cnt1>=count){
                    st = ent.getKey();
                    count = cnt1;
                }
            }
            if(left) return st.substring(0, st.length()-(end-start+1));
            else return  st.substring((end-start+1), st.length());
        }
   
      
        EHHCalculator(String core, boolean left) {
            this.core = core;
            this.left = left;
            List<PIGData> dat = sdt.getAllHaplotypes(start, end).get(core);
            originalNumber = dat.size();
          //  if(originalNumber>=2){
                hapl_1 = new TreeMap<String, List<PIGData>>();
                hapl_1.put(core, dat);
               // haplotypes.add(hapl_1);
                pos++;
                scores.add(numEquiv(hapl_1));
              //for(int i=start; i<=end; i++){
                this.loc1.add(sdt.loc.get(left ? start : end));
                this.ch.append(left ? core.charAt(0) : core.charAt(core.length()-1));
              //
            //}
        }
        
        public void printSum(Logger  log){
                log.info(name+" "+start+" "+end+" "+core);
                log.info(scores+" ");
              //  log.info(haplotypes+" ");
        }
        public double sum(){
            //List<Integer> loc = sdt.loc;
            double sum=0;
            for(int i=0; i<scores.size()-1; i++){
            //    int pos = left ? -1 :1;
            //    int start_ = start+pos*(i+1);
           //     int end_ =  start+pos*(i);
                double midpoint = (scores(i)+scores(i+1))/2.0;
                double width = Math.abs(loc1.get(i+1) - loc1.get(i));
                sum+=width*midpoint;
            }
            return sum;
        }
        public double scores(int i){
            return scores.get(i)==0 ? 0 : scores.get(i)/scores.get(0);
        }
        public int getNextStart(){
            if(left){
               // int pos = haplotypes.size();
                return start - pos; 
            }
            else return start;   
        }
        public int getNextEnd(){
            if(left) return end;
            else{
              //  int pos = haplotypes.size();
                return end+pos;
            }
        }
        
        public boolean skip(int i){
        	return false;
/*           // if(true) return false;
            if(sdt.maf.getFixedInteger(i)!=null) return true;
            else{
                double[] d = sdt.maf.getEmiss(i);
                for(int k=0; k<d.length; k++){
                    if(d[k]>0 && sdt.maf.getEmissionStateSpace().getCN(k)!=1) {
                        return true;
                    }
                }
            }
            return false;*/
        }
        public boolean extend(){
          //  int pos = haplotypes.size();
            int start_ = getNextStart();
            int end_ =  getNextEnd();
            pos++;
            if(left){
                if(skip(start_)) {
                    return true;
                }
            }
            else{
                if(skip(end_)) return true;
            }
            SortedMap<String, List<PIGData>> map_1 =this.hapl_1;// haplotypes.get(pos-1);
            SortedMap<String, List<PIGData>> m_pos = new TreeMap<String, List<PIGData>>();
        
          
            
            this.loc1.add(left ? sdt.loc.get(start_) : sdt.loc.get(end_));
           int max_size =0;
           String commonHapl = "";
            for(Iterator<Entry<String, List<PIGData>>> it = map_1.entrySet().iterator();it.hasNext();){
                Entry<String, List<PIGData>> nxt = it.next();
                for(int i=0; i<nxt.getValue().size(); i++){
                    PIGData dat_i = nxt.getValue().get(i);
                    String hapl = dat_i.getStringRep(start_, end_);
                    List<PIGData> l_i = m_pos.get(hapl);
                    if(l_i==null){
                        m_pos.put(hapl, l_i = new ArrayList<PIGData>());
                    }
                    l_i.add(dat_i);
                    if(l_i.size()>max_size){
                        max_size = l_i.size();
                        commonHapl = hapl;
                    }
                }
            }
            scores.add(numEquiv(m_pos));
            this.ch.append(left ? commonHapl.charAt(0) : commonHapl.charAt(commonHapl.length()-1));
            hapl_1 = m_pos;
            return true;
        }
        /*  public double ehh(int i){
           return numEquiv(i);
           double numEquiv_0 = numEquiv(0);
            double res =  numEquiv_i / numEquiv_0;
            if(numEquiv_0==0){
                throw new RuntimeException("!!");
            }
            return res;
        }*/
        public double combinations(double numHapl){
            return  
          //  withReplacement ? 
          //          numHapl*(numHapl+1.0) / 2.0:    
                        numHapl*(numHapl-1.0) / 2.0;
        }
        public double numEquiv(SortedMap<String, List<PIGData>> map){
            double sum=0.0;
           
            for(Iterator<List<PIGData>> it = map.values().iterator(); it.hasNext();){
                //if(numHapl==1) throw new RuntimeException("!!");
                sum+=  combinations((double) it.next().size());
            }
            return sum;
        }

        public double ehhMin() {
           return scores(this.scores.size()-1);
           
        }
        /** is a next possible? */
        public boolean canExtend(){
            if(getNextStart()<0 || getNextEnd()>=sdt.length()) return false;
            else return true;
        }
        public boolean hasNext(){
            if(canExtend() && ehhMin()>=thresh){
                return true;
            }
            else return false;
        }

        public char ch(int k) {
            // TODO Auto-generated method stub
             return ch.charAt(k);
        }
       
        
       
    
    }
public double minEHH=0;
static boolean[] left = new boolean[] {true, false};



public Integer findAncestral(int start, int end, String derived, String[] ancestral, boolean[] doLR){
    this.start = start;
    this.end = end;
    Integer[] anc_index = new Integer[2];
    for(int i=0; i<left.length; i++){
        if(!doLR[i]) continue; 
        EHHCalculator coredDerived = new EHHCalculator(derived, left[i]);//map.remove(derived);
        EHHCalculator[] coreAncestral = new EHHCalculator[ancestral.length];
        String st_der;
        int[] st_anc = new int[coreAncestral.length];
        for(int j=0; j<coreAncestral.length; j++){
            coreAncestral[j] = new EHHCalculator(ancestral[j], left[i]);
        }
        while(coredDerived.canExtend()){
            coredDerived.extend();
            st_der = coredDerived.mostCommon();
            String[] anc_com = new String[coreAncestral.length]; 
            for(int j=0; j<coreAncestral.length; j++){
                coreAncestral[j].extend();
               anc_com[j] = coreAncestral[j].mostCommon();
                st_anc[j] = this.noCharInCommon(st_der,anc_com[j] , left[i]);
            }
            int max_index = Constants.getMax(st_anc);
            boolean greaterThanAll = st_anc[max_index]>0;
            for(int j=0; j<coreAncestral.length && greaterThanAll; j++){
                if(j==max_index) continue;
                if(st_anc[max_index] <= st_anc[j]) greaterThanAll = false;
            }
            if(greaterThanAll){
                anc_index[i] = max_index; 
            }
        }
    }
    if(anc_index[0]==null) return anc_index[1];
    else if(anc_index[1]==null) return anc_index[0];
    else if(anc_index[0]==anc_index[1]) return anc_index[0];
    else return null;
}
//EHHCalculator[][]in = new EHHCalculator[Constants.modelCNP()][2]; //first index is allele, second is left or right

/** note - this is all based on nC2, which is without replacement
 * returns [[left_d, left_a], [right_d, right_a]
 * 
 * ancestral is a list of potential ancestral alleles
 *  */
    public void setCore(int start, int end, String[] cores, Double[][]res, boolean[] doLR, List<Double>[][] in1){
     //  Double[][] res = new Double[2][2];
        super.setCore(start, end, cores, res, doLR, in1);
        //derived, ancestral, res, doLR, in1);
        EHHCalculator[] ehh = new EHHCalculator[cores.length];
       // boolean[] docalculation = new boolean[cores.length];
         
        for(int i=0; i<left.length; i++){
            if(!doLR[i]) continue; 
            int numNonFalse = 0;
            for(int j=0; j<cores.length; j++){
               ehh[j] = new EHHCalculator(cores[j], left[i]);//map.remove(derived);
               // EHHCalculator coreAncestral =new EHHCalculator(cores[j], left[i]);
                if(ehh[j].originalNumber<2)ehh[j] = null;
                else numNonFalse++;
            }
            if(numNonFalse<=1) continue;  //no point if only one core can 
            while(hasNextCore(ehh)){
                extendCore(ehh);
            }
            for(int j=0; j<ehh.length; j++){
                res[i][j] =  ehh[j].sum();
                in1[i][0] = ehh[j].scores;
            }
        
            
        /*
            if(coredDerived.canExtend()){
                coredDerived.extend();
                coreAncestral.extend();
            }*/
        }
       /*  if(res[0][0]+res[1][0]) ==0 ||(res[0][1]+res[1][1]) ==0){
             throw new RuntimeException("!!");
         }*/
        
    }
   
   private void extendCore(EHHCalculator[] ehh) {
       for(int i=0; i<ehh.length; i++){
          ehh[i].extend();
        }
    
}

private boolean hasNextCore(EHHCalculator[] ehh) {
   for(int i=0; i<ehh.length; i++){
      if(ehh[i].hasNext()) return true;
   }
   return false;
}

static boolean plot = true;

   public double score(Double[][] scores) {
       return Math.log((scores[0][0]+scores[1][0])/(scores[0][1]+scores[1][1]));
   } 
    
   
    
   
}
