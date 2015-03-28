/**
 * 
 */
package lc1.dp.core;

import lc1.util.Constants;

import org.apache.commons.pool.PoolableObjectFactory;
import org.apache.commons.pool.impl.StackObjectPool;

public class DoubleDoublePool {
    StackObjectPool[][] double_pool ;
  
  //  int idle;
  //  int cap;
    
    public DoubleDoublePool(int max){
        this.double_pool = new StackObjectPool[max+1][max+1];
        for(int i=0; i<=max; i++){
            final  int i1 = i;
            for(int j=0; j<=max; j++){
           final int j1 = j;
            double_pool[i1][j1] = new StackObjectPool(new PoolableObjectFactory(){

             public void activateObject(Object arg0) throws Exception {}
             public void destroyObject(Object arg0) throws Exception {}
             public Object makeObject() throws Exception {
                return new double[i1][j1];
             }
             public void passivateObject(Object arg0) throws Exception {}
             public boolean validateObject(Object arg0) {return true;}}   ,Math.max(2, Constants.numThreads()), Constants.numThreads());
     }  
        }
    }
    
    
 /*   MyObjectPool(){
        super(;
    }*/
    public  double[][] getObj(int len, int len1) {
        try{
        return (double[][])double_pool[len][len1].borrowObject();
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
        }
        return null;
        
    }
    public synchronized void returnObj(double[][] j) {
     try{
        this.double_pool[j.length][j[0].length].returnObject(j);
     }catch(Exception exc){
         exc.printStackTrace();
     }
    }
    
}