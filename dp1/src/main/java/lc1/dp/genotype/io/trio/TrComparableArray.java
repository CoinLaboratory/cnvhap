package lc1.dp.genotype.io.trio;

import lc1.dp.data.representation.ComparableArray;

public abstract class TrComparableArray extends ComparableArray{
    public ComparableArray third;  // this is the indices;
    
    TrComparableArray(){
        super(true);
    }
    @Override
    public Comparable getReal(int j){
        if(j==this.size()-1) return third;
        else return this.get(j);
    }
}
