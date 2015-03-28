package lc1.util;

import java.io.Serializable;
import java.util.Iterator;
/** enumerates all possible copies of an underlying iterator */
public abstract class CopyEnumerator implements Serializable {
    final  int no_copies;
   
    /** underlying iterator */
    public abstract Iterator<Comparable> getPossibilities(int depth);
    /*what to do with object generated */
    public abstract void doInner(Comparable[] list);
    Comparable[] list;
    public CopyEnumerator(int length){
        list = new Comparable[length];
        this.no_copies = length;
    }
    public void run(){
        inner(0);
    }
    public abstract boolean exclude(Object obj, Object previous);
    
    public abstract boolean exclude(Comparable[] list);
    
    private void inner(int depth){
        for(Iterator<Comparable> it = getPossibilities(depth); it.hasNext();){
            Comparable nxt =  it.next();
            if(exclude(nxt, depth==0 ? null : list[depth-1])) continue;
            list[depth] =nxt;
            if(depth+1 == list.length){
                doInner(list);
            }
            else{
                inner(depth+1);
            }
        }
    }
}
