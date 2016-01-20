package lc1.dp.core;

import java.util.logging.Logger;

import lc1.dp.model.MarkovModel;
import lc1.dp.model.WrappedModel;
import lc1.util.Constants;

import org.apache.commons.pool.PoolableObjectFactory;

public class DPPool extends MyObjectPool{

	public DPPool(final MarkovModel hmm, int max, int numThreads) {
		super(new PoolableObjectFactory(){

	         public void activateObject(Object arg0) throws Exception {}
	         public void destroyObject(Object arg0) throws Exception {}
	         public Object makeObject() throws Exception {
	             System.err.println("MAKING NEW DP OBJECT ");
	             return new DP(hmm,"",Constants.isLogProbs(), hmm.noSnps , false);
	         }
	         public void passivateObject(Object arg0) throws Exception {}
	         public boolean validateObject(Object arg0) {return true;}}
	         , max, numThreads);
	}
	@Override
	 public synchronized void returnObj(Object j) throws Exception{
      

	        Object obj = this.inuse.remove(j);
	        ((DP)obj).inuse = false;
	        this.returnObject(obj);
	       
	        if(print) Logger.global.info("returning1 "+j+ " "+inuse.size()+ "  "+obj.getClass());
	    }

	public DPPool(final MarkovModel hmm1, final int noSamples, final boolean unwrapToSample, int max, int numThreads) {
		  super(new PoolableObjectFactory(){

              public void activateObject(Object arg0) throws Exception {}
              public void destroyObject(Object arg0) throws Exception {}
              public Object makeObject() throws Exception {
                  System.err.println("MAKING NEW DP OBJECT -SAMPLER");
                  MarkovModel hmm = hmm1;
                  if(unwrapToSample){
                      while(hmm instanceof WrappedModel){
                       hmm = ((WrappedModel)hmm).getHMM();
                      } 
                  }
                  return new DP(hmm,"", noSamples<=1, hmm.noSnps, noSamples>1 );
              }
              public void passivateObject(Object arg0) throws Exception {}
              public boolean validateObject(Object arg0) {return true;}}   ,max, numThreads);
    
	}
	
}
