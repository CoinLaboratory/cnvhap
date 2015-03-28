package lc1.dp;

import java.io.File;

import lc1.dp.genotype.ExponentialHaplotypeHMM;
import lc1.dp.genotype.io.FastPhaseFormat;
import lc1.dp.genotype.io.scorable.HaplotypeData;
import lc1.dp.genotype.io.scorable.ScorableObject;

public class HaplotypeTest extends DPTest{

HaplotypeTest(int no) throws Exception{
    super(no,1);
}
MarkovModel getHMM(){
    File dir = new File(System.getProperties().getProperty("user.dir"));
    try{
    ExponentialHaplotypeHMM hmm = FastPhaseFormat.readHMM(dir);
        
        //new ExponentialHaplotypeHMM("geno",2, 100,pseudocount,2.0);
    //hmm.increaseStates();
   // hmm.doubleStates(0.5, 20);
    return hmm;
    }catch(Exception exc){
        exc.printStackTrace();
        return null;
    }
}
ScorableObject makeScoreableObject(int k) {
    
    return new HaplotypeData(k+"",5);
    
}


}
