package lc1.possel;

import java.util.List;

import lc1.dp.data.collection.DataCollection;

public abstract class AbstractEHH {
    int start, end;
    final String name;
    DataCollection sdt;
    
    public void setCore(int i, int i2,String[] string,  Double[][] scores, boolean[] doLR, List<Double>[][] in) {
        this.start =i;
        this.end = i2;
    }
    public AbstractEHH(DataCollection coll,String f) throws Exception{
        this.sdt = (DataCollection) coll.clone();
        this.name = f;
    }
    
    public abstract double score(Double[][] scores) ;
}
