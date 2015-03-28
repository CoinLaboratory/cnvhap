/**
 * 
 */
package lc1.dp;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

class DicePairState extends PairEmissionState{
    DicePairState(List<EmissionState> st1){
        super(st1);
      
    }
   

     public Object combine(Object[] o1) {
         int res = 0;
         for(int j=0; j<o1.length; j++){
             res+=((Integer)o1[j]).intValue();
         }
         return res;
     }
    @Override
    public Object clone() {
       return new DicePairState(copy(this.dist));
    }
    

     
}