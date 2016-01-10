package lc1.dp.emissionspace;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import lc1.dp.data.representation.AbstractEmiss;
import lc1.dp.data.representation.ComparableArray;
import lc1.util.Constants;



//Note default list behaviour is over haplo pairs - ie unorder lists of haplotypes!!!!
public  class SimpleEmissionStateSpace extends EmissionStateSpace{


    public SimpleEmissionStateSpace(List<Comparable> list){
        super();
     //   if(list.size()>2) throw newy RuntimeException("!!"); 
        init(list, false);
      
    }
    public List<String>nme = null;
    public SimpleEmissionStateSpace(List<Comparable> list, List<String>nme){
        super();
     //   if(list.size()>2) throw newy RuntimeException("!!"); 
        init(list, false);
        this.nme = nme;
    }
    
    public SimpleEmissionStateSpace(EmissionStateSpace[] stsp){
    	super();
    	Set<Comparable> l = new HashSet<Comparable>();
    	for(int i=0; i<stsp.length; i++){
    		l.addAll(stsp[i].getGenotypeList());
    	}
    	init(new ArrayList<Comparable>(l), false);
    } 
    public String getHaploPairString(Comparable comp){
        return getHaploString(comp);
    }
    public String getGenotypeString(int i){
        return this.getGenotypeString(this.getGenotype(i));
    }
    public String getGenotypeString(Comparable comp){
        if(comp instanceof Integer){
            return Integer.toString((Integer)comp, Constants.radix());
        }
        else if (comp instanceof ComparableArray){
        	return ((ComparableArray)comp).getGenotypeString();
        }
        
        return ((AbstractEmiss)comp).toStringShort().replaceAll("_", "");
    }
    public String getHaploString(Comparable comp){
      return getGenotypeString(comp);
    }
    
  
    
   public SimpleEmissionStateSpace(Comparable[] stateSpace) {
       this(getList(stateSpace));
      
       
    }

private static List<Comparable> getList(Comparable[] stateSpace) {
    List<Comparable> comp = new ArrayList<Comparable> (stateSpace.length);
    for(int i=0; i<stateSpace.length; i++){
        comp.add(stateSpace[i]);
    }
    return comp;
}

public void setNme(int i, String name) {
	this.nme.set(i, name);
	
}

    
    
}
