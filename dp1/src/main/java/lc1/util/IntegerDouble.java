/**
 * 
 */
package lc1.util;


public class IntegerDouble{
    public Object stateId;
    public double d;
    public IntegerDouble(Object i, double d){
        this.stateId = i;
        this.d = d;
    }
    
    public void offset(int n){
        stateId = ((Integer)stateId) - n;
    }
    IntegerDouble(IntegerDouble i1){
        stateId = i1.stateId;
        d = i1.d;
    }
    
    public void set(Object stateId, double sc){
        this.d =sc;
        this.stateId = stateId;
    }
    public String toString(){
        return stateId+":"+String.format("%5.3g", new Object[] {d});
    }
    public void add(double sc){
        
        d+=sc;
    }
}