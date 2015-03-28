/**
 * 
 */
package lc1.dp;

import java.util.List;

class DicePairModel extends PairMarkovModel{
DicePairModel(MarkovModel m1, int no_copies){
    super(m1,no_copies);
}
    public Object clone() {
      //  if(m1!=m2) throw new RuntimeException("!!");
        MarkovModel m = (MarkovModel)this.m1.clone();
      return new DicePairModel(m, this.no_copies);
    }
    @Override
    public PairEmissionState constructPair(List st1) {
        // TODO Auto-generated method stub
        return new DicePairState(st1);
    }
}