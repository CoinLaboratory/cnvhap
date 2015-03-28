package lc1.dp.genotype.io.trio;

import lc1.dp.data.representation.ComparableArray;

public class TrioComparableArray extends TrComparableArray {
  
   public TrioComparableArray(ComparableArray comp){
        super();
        if(comp.size()!=3) throw new RuntimeException("!!");
        for(int i=0; i<2; i++){
            this.add(comp.get(i));
        }
        third = (ComparableArray) comp.get(2);
        ComparableArray third_pos = new ComparableArray(false);
        ComparableArray th0= (ComparableArray) comp.get(0);
        ComparableArray th1 = (ComparableArray) comp.get(1);
       
        int index1 = (Integer)third.get(0);
        int index2 = (Integer)third.get(1);
        third_pos.add(th0.get(index1-1)); third_pos.add(th1.get(index2-1));
        this.add(third_pos);
    }
 
}
