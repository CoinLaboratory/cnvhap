/**
 * 
 */
package lc1.dp.core;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.StackObjectPool;

public class MyObjectPool extends StackObjectPool{
    
    Map<Object, Object> inuse = Collections.synchronizedMap(new HashMap<Object, Object>());
    
    public MyObjectPool(PoolableObjectFactory fact, int idle, int cap){
        super(fact, idle, cap);
    }
 /*   MyObjectPool(){
        super(;
    }*/
    public  Object getObj(Object j) throws Exception{
        Object res = inuse.get(j);
        if(res==null){
            inuse.put(j,res = this.borrowObject());
        }
        return res;
    }
    public synchronized void returnObj(Object j) throws Exception{
        Object obj = this.inuse.remove(j);
        this.returnObject(obj);
    }
    
}