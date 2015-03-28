/*
 * Created on 17-Aug-2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package lc1.dp.states;

import java.io.Serializable;

import lc1.stats.SimpleDistribution;


/**
 * @author lc
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public abstract class State implements Serializable, Comparable{
	public  String name;
	final public int adv;
    public static int endStateLength = 5;
    
   private   int index=-1;
   
   public void setIndex(int i){
      // if(index!=-1 && index!=i) throw new RuntimeException("changing index");
       index = i;
   }
    public int getIndex(){
       // if(index==-1) throw new NullPointerException("index is null");
        return index;
    }
    
    public int compareTo(Object o){
        State st = (State) o;
        return name.compareTo(st.getName());
    }
    public abstract Object clone();
    
   
    //  public final static Distribution oneOff = new Distribution(new short[] {1}, new double[] {1});//3, 0.01, true);
	public String getName(){
		return name;
	}
    /** this is used for collapsing the map from state to emission score */
   public String getEmissionName(){
       return name;
   }
   // Distribution dist;  //distribution of lengths around adv
    public abstract SimpleDistribution adv(State s);
    
    abstract public void validate() throws Exception;

   
    public String toString(){
		return name+"";
	}
	public State(String name1, int adv){
      //  if(name1=="") throw new RuntimeException("need name");
		this.name = name1;
		this.adv = adv;
	}
    public State(State st1){
        this(st1.name, st1.adv);
        this.index = st1.index;
    }
    
	/*public int hashCode(){
     return name;   
    }*/
    public boolean equals(Object o){
        return ((State)o).getName().equals(this.getName());
    }
   public int hashCode(){
       return this.getName().hashCode();
   }
    
   // public abstract double score(Object element, boolean logspace);
        
}
