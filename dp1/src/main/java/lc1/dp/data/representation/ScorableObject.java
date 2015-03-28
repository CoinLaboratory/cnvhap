package lc1.dp.data.representation;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Collection;
import java.util.List;

import lc1.stats.PseudoDistribution;

public interface   ScorableObject extends Serializable, Comparable{
    public  int length(); //returns length of object
    public   Comparable getElement(int i);
    public  void print(PrintWriter pw, boolean expand, boolean mark, Collection<Integer> toDrop,
    		List<Character> alleleA, 
            List<Character> alleleB, PseudoDistribution[] ems);
    public  String getName();
    public abstract Object clone();
}
