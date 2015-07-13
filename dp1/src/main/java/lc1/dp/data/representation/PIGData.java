package lc1.dp.data.representation;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import lc1.CGH.Aberation;
import lc1.dp.states.EmissionState;

public interface PIGData extends CSOData, SSOData, ScorableObject{


    /*public void collapse(){
        for(int i=0; i<this.length(); i++){
            ComparableArray comp =(ComparableArray) this.getElement(i);
            ComparableArray comp_n = new ComparableArray(comp.orderKnown());
            comp_n.add(comp.get(0));
            this.l1.set(i, comp_n);
        }
    }*/

 //   public abstract ScorableObject clone();

    /*
    Double[] certainty;
    
    public void setCertainty(Integer pos, double cert) {
       if(certainty ==null) certainty = new Double[this.length()];
       certainty[pos] = cert;
        
    }*/

    public abstract String toString();

   // public abstract PIGData[] split(int[] is);

    public abstract PIGData recombine(
            Map<Integer, Integer> recSites, int starting_index);

    public abstract Set<Integer>[] getSwitches();

    public abstract Collection<Aberation> getDeletedPositions(EmissionState st, int noCop, Boolean deletion);

    public abstract String getStringRep(int start, int end);

    public abstract void removeAll(List<Integer> toDrop);

    public abstract void reverse();

    public abstract int compareTo(Object arg0);

    public abstract void switchAlleles(int i);

    public abstract void setAsMissing(List<Integer> toDrop, double cn_ratio);

    public abstract void applyAlias(int[] alias);



  //  public abstract int countNull();

    

}