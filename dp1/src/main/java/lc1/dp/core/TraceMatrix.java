/**
 * 
 */
package lc1.dp.core;

import java.util.Arrays;

import lc1.util.Constants;


public class TraceMatrix{
    final public AbstractTerm[][] trace;
    final private double[] logscale;
    double overall; 
    final AbstractTerm nullAbstractTerm ;
    final int seqLength;
    final boolean forward;
    
    
    public AbstractTerm getTrace(int j, int i){
        if(i>=0 && i<seqLength){
            AbstractTerm[] term = trace[j];
            if(term==null) return null;
            else return term[i];
        }
        else if(j==0 && i==-1 && forward) return nullAbstractTerm;
        else if(j==0 && i==seqLength && !forward) return nullAbstractTerm;
        else return null;
    }
    
 public double getScore(int j, int i, boolean logspace){
     AbstractTerm t = getTrace(j,i);
     if(t==null) return logspace ? Double.NEGATIVE_INFINITY : 0;
     else return t.score();
 }
    
    public double getLogScale(int i){
        if((i==-1 && forward)||( i==seqLength && !forward)) return 0;
        else return logscale[i];
    }
    
    public void setLogscale(int i, double d){
        this.logscale[i] = d;
    }
    
    
    public void clear(){
        Arrays.fill(logscale, 0);
        /*for(int j=0; j<trace.length; j++){
           for(int i=0; i<seqLength; i++){
               this.trace[j][i].clear();
           }
        }*/
        overall = 0;
    }
    
    
    
TraceMatrix(int modelLength, int seqLength,  boolean forward, boolean logspace, Class clazz) throws Exception
{
    this.forward = forward;
    nullAbstractTerm = new Term(-1,-1, logspace ? 0.0 : 1.0);
    trace = new AbstractTerm[modelLength][];
    for(int i=0; i<trace.length; i++){
        trace[i] = new AbstractTerm[seqLength];
        for(int j=0; j<seqLength; j++){
            trace[i][j] =(AbstractTerm) clazz.getConstructor(new Class[] {int.class, double.class, int.class}).newInstance(new Object[] {0,0, modelLength});
        }
    }
    this.seqLength = seqLength;
    logscale = new double[seqLength];
    Arrays.fill(logscale, 0.0);
}
    double[] minScore(int i){
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for(int j=0; j<trace.length; j++){
            AbstractTerm[] trace_j = trace[j];
            if(trace_j==null || trace_j[i]==null ) continue;
        //    if(trace_j[i].sc)
            if(trace_j[i].scaleScore() < min && trace_j[i].scaleScore()>0){
                min = trace_j[i].scaleScore();
            }
            if(trace_j[i].scaleScore()>max){
                max = trace_j[i].scaleScore();
            }
        }
       
        return new double[] {Math.log(min),Math.log(max)};
    }
    
    boolean  allZero(int i){
       
        for(int j=0; j<trace.length; j++){
            AbstractTerm[] trace_j = trace[j];
            if(trace_j==null || trace_j[i]==null ) continue;
            if( trace_j[i].scaleScore()>0){
               return false;
            }
          
        }
       return true;
    }
    
    Double[] getTr(int i){
    	Double[] res = new Double[trace.length];
    	 for(int j=0; j<trace.length; j++){
             AbstractTerm[] trace_j = trace[j];
             if(trace_j==null || trace_j[i]==null ) continue;
          
                 res[j] = trace_j[i].scaleScore();
            
         }
    	 return res;
    }
    
    void  scale(double d, int i){
    	//if(Double.isNaN(d)) throw new RuntimeException("!!");
        for(int j=0; j<trace.length; j++){
            AbstractTerm[] trace_j = trace[j];
            if(trace_j!=null && trace_j[i]!=null) trace_j[i].scale(d); 
        }
       if(Constants.CHECK && this.allZero(i)){
    	   throw new RuntimeException("!!");
       }
    }

    public void setTrace(int j, int i, int max_j, int max_i, double sc) {
    	//if(Constants.CHECK && Double.isNaN(sc)) {
    	//	throw new RuntimeException("is nan "+sc+" "+j+" "+i);
    	//}
         AbstractTerm term = trace[j][i];
         term.i = max_i;
       
         term.score = sc;
         term.setj(max_j);
      /* if(sc==0){
        	 throw new RuntimeException("!!");
         }*/
    }
    
    public void setTrace(int j, int i,int max_i, double sc) {
    	//if(Constants.CHECK && Double.isNaN(sc)) throw new RuntimeException("is nan "+sc+" "+j+" "+i);
        
    	AbstractTerm term = trace[j][i];
        term.i = max_i;
        term.score = sc;
   }

    public double[] getDoubleArray(int j, int i) {
        // TODO Auto-generated method stub
        ComplexTerm term = (ComplexTerm) trace[j][i];
        return term.prob;
    }

   /* public void setTrace(int j, int i, AbstractTerm term) {
        AbstractTerm[] t = trace[j];
        if(t==null){
            trace[j]= t = new AbstractTerm[seqLength];
        }
        t[i]= term;
    }*/

 }