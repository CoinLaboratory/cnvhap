package lc1.dp;

import java.util.Iterator;
import java.util.Map.Entry;

import lc1.dp.genotype.io.scorable.ScorableObject;
import lc1.dp.genotype.io.scorable.SimpleScorableObject;

public class DPDiceTest extends DPTest{

    static double d1 = 0.16;
    static double d2 = 0.7;
    static double sixth = (1.0 - d1)/5.0;
    static double tenth = (1 - d2)/5.0;
    
    static double[] dist1 = new double[] {sixth, sixth,sixth,sixth, sixth,d1};
    static double[] dist2 = new double[] {tenth, tenth, tenth, tenth, tenth, d2};
DPDiceTest(int no, int no_copies1) throws Exception{
    super(no, no_copies1);
   
}

class DiceMarkovModel extends FastMarkovModel{
    final int noSnps;
    public DiceMarkovModel(String name, DotState magic, double pseudocount, int noSnps) {
        super(name, magic, pseudocount);
        this.noSnps = noSnps;
    }
    public DiceMarkovModel(DiceMarkovModel dmm ){
        super(dmm);
        this.noSnps = dmm.noSnps;
    }
   
    public Object clone(){
        return new DiceMarkovModel(this);
    }
  
    
    public double getTransitionScore(State from, State to, int indexOfToEmission) {
        if(indexOfToEmission>=noSnps) return 0;
        if(from==MAGIC){
            return this.transProb.getTransition(from, to);//statesOut(from).dist.get(to);
        }
        else if(to==MAGIC){
            if(indexOfToEmission==noSnps-1) return 1.0;
            else return 0;
        }
        else if(indexOfToEmission==0) return 0;
        else return this.transProb.getTransition(from, to);//statesOut(from).dist.get(to);
    }
    
}

MarkovModel getHMM(){
    State st1 = new DiceEmissionState("F", 1,dist1, this.pseudocount);
    State st2 = new DiceEmissionState("L", 1,dist2, pseudocount);
    final int noSnps = 10;
    final SimpleDistribution toMagic = new SimpleDistribution();
    toMagic.dist.put(magic, 1.0);
    final SimpleDistribution magicIn = new SimpleDistribution();
    magicIn.dist.put(st1, 0.5);
    magicIn.dist.put(st2, 0.5);
    FastMarkovModel hmm = new FastMarkovModel("dice", magic,pseudocount);
      double mag = 1e-3;
    hmm.addState(st1);
    hmm.addState(st2);
    hmm.setTransition(st1, st1, 0, 0.95-mag);
    hmm.setTransition(st1, st2, 0, 0.05);
    hmm.setTransition(st1,magic, 0,mag);
    
    hmm.setTransition(st2, st2,  0,0.9-mag);
    hmm.setTransition(st2, st1,  0,0.1);
    hmm.setTransition(st2,magic, 0, mag);
    
    hmm.setTransition(hmm.MAGIC, st1, 0,1.0);
    if(no_copies==1) return hmm;
    else return new DicePairModel(hmm, no_copies);
}
ScorableObject makeScoreableObject() {
    return new SimpleScorableObject("", 10, Integer.class);
}

}
