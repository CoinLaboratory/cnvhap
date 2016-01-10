/**
 * 
 */
package lc1.dp.core;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.StackObjectPool;

public class MyObjectPool extends StackObjectPool{
	public static final boolean print = false;  

    Map<Object, Object> inuse = Collections.synchronizedMap(new HashMap<Object, Object>());
    int max;
    public MyObjectPool(PoolableObjectFactory fact, int idle, int cap){
        super(fact, idle, cap);
        max = cap;
    }
 /*   MyObjectPool(){
        super(;
    }*/
    public  synchronized Object getObj(Object j) throws Exception{
    	//System.err.println("getting1 "+j+ " "+inuse.size());
        Object res = inuse.get(j);
        if(res==null){
        	if(this.getNumActive()>= this.max || this.inuse.size()>=max){
        		throw new RuntimeException("need to return first "+j+" "+inuse.size()+" "+getNumActive()+" "+max);
        	}
            inuse.put(j,res = this.borrowObject());
        }
       if(print) Logger.global.info("getting1 "+j+ " "+inuse.size()+ " "+res.getClass());
        return res;
    }
    public synchronized void returnObj(Object j) throws Exception{
    
        Object obj = this.inuse.remove(j);
        if(print) Logger.global.info("returning1 "+j+ " "+inuse.size()+ "  "+obj.getClass());
      //  System.err.println("removed "+obj);
        this.returnObject(obj);
    }
    
}