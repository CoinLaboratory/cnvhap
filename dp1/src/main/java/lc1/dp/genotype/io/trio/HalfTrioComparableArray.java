package lc1.dp.genotype.io.trio;

import lc1.dp.data.representation.ComparableArray;

public class HalfTrioComparableArray extends TrComparableArray {
    
    public HalfTrioComparableArray(ComparableArray comp){
        super();
        if(comp.size()!=2) throw new RuntimeException("!!");
        this.add(comp.get(0));
        this.third = (ComparableArray) comp.get(1);
        ComparableArray third_pos = new ComparableArray(false);
        ComparableArray th0= (ComparableArray) comp.get(0);
        int index1 = (Integer)third.get(0);
        third_pos.add(th0.get(index1-1)); 
        this.add(third_pos);
    }
 
}
