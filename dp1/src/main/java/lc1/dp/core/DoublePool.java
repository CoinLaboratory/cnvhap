/**
 * 
 */
package lc1.dp.core;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.StackObjectPool;

public class DoublePool{
    List<StackObjectPool> double_pool = new ArrayList<StackObjectPool>();
    int size =0;
    
  
  //  int idle;
  //  int cap;
    
    public void addPool(final int len){
        for(int i=size; i<=len-2; i++){
           final  int i1 = i+2;
        StackObjectPool pool = new StackObjectPool(new PoolableObjectFactory(){

            public void activateObject(Object arg0) throws Exception {}
            public void destroyObject(Object arg0) throws Exception {}
            public Object makeObject() throws Exception {
               return new double[i1];
            }
            public void passivateObject(Object arg0) throws Exception {}
            public boolean validateObject(Object arg0) {return true;}}   ,4, 4);
        double_pool.add(pool);
        }
        size = double_pool.size();
    }
    
 /*   MyObjectPool(){
        super(;
    }*/
    public  double[] getObj(int len) {
        if(len-2>=size) addPool(len);
        try{
        return (double[]) double_pool.get(len-2).borrowObject();
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
        }
        return null;
        
    }
    public synchronized void returnObj(double[] j) {
     try{
        this.double_pool.get(j.length-2).returnObject(j);
     }catch(Exception exc){
         exc.printStackTrace();
     }
    }
    
}