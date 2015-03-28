package lc1.dp;

import java.io.File;
import java.util.Arrays;

import lc1.dp.genotype.ExponentialHaplotypeHMM;
import lc1.dp.genotype.GenotypeHMM;
import lc1.dp.genotype.io.FastPhaseFormat;
import lc1.dp.genotype.io.scorable.PhasedIntegerGenotypeData;
import lc1.dp.genotype.io.scorable.ScorableObject;

public class GenotypeTest extends DPTest{

GenotypeTest(int no, int no_copies) throws Exception{
    super(no, no_copies);
    File dir = new File(System.getProperties().getProperty("user.dir"));
    FastPhaseFormat.write(Arrays.asList(data), new File(dir, "simulated.dick.txt"));
}
static int noSites =50;
MarkovModel getHMM(){
    File dir = new File(System.getProperties().getProperty("user.dir"));
    try{
        ExponentialHaplotypeHMM hmm = FastPhaseFormat.readHMM(dir);
        return new GenotypeHMM(hmm,no_copies);
        }catch(Exception exc){
            exc.printStackTrace();
            return null;
        }
  //  ExponentialHaplotypeHMM hmm = new ExponentialHaplotypeHMM("geno", 2, noSites,pseudocount, 2);
   
}

protected void addPoint(int k, int i) {
    PairEmissionState t = (PairEmissionState)states[k][i];
    ((PhasedIntegerGenotypeData)data[k]).addPoint(t.sample(i));//(Boolean) t.dist1.sample(i), (Boolean) t.dist2.sample(i));
   }
ScorableObject makeScoreableObject(int k) {
    return new PhasedIntegerGenotypeData(k+"", 5,2);
  /*  return new ScorableObject(){
        public String toString(){
            return l.toString();
        }
        List<Integer> l = new ArrayList<Integer>();

        public int length() {
           return l.size();
        }

        public Object getElement(int i) {
         return l.get(i);
        }

        public void addPoint(Object obj) {
          l.add((Integer)obj);
        }
        
    };*/
    
}


}
