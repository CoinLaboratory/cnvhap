package lc1.dp.states;

import java.io.PrintWriter;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;

public class SiteEmissionState extends EmissionState {

    
    public double[] emissions; //prob of emitting true
   final int noSnps;
  boolean logspace;
  @Override
  public boolean transferCountsToProbs(double pseudo) {
   throw new RuntimeException("!!");
      
  }
   /** alias is site to site pattern */
   public SiteEmissionState(String name, double[] logprob, int[] alias, boolean logspace){
       super(name, 1);
       this.noSnps =alias.length;
       this.emissions =new double[alias.length];
       for(int i=0; i<noSnps; i++){
           emissions[i] = logspace ? logprob[alias[i]] : Math.exp(logprob[alias[i]]);
           if(emissions[i]<1e-300) throw new RuntimeException("danger of underflow");
       }
   }
  
   
   public void print(PrintWriter pw, String prefix){
       StringBuffer sb1= new StringBuffer(prefix);
       for(int i=0; i<this.noSnps; i++){
          sb1.append("%8.2g ");
       }
      // pw.print(Format.sprintf(sb1.toString(), emissions(this.emissions)));   
   }
 
/* public double[] score(StateIndices obj, boolean logspace){
       if(logspace!=this.logspace) throw new RuntimeException("mismatch");
       return emissions;
   }*/
  
   public Object clone(){
      throw new RuntimeException("!!");
   }
   
   /*[true, false, null] */
    public void addCount(int  obj_index, double value, int i) {
        throw new RuntimeException("!!");
    }
    public void addCountDT( double data_index, int phen_index, double value, int i) {
        throw new RuntimeException("!!");
    }
   
    public void initialiseCounts(){
        throw new RuntimeException("!!");
        }

   
    
    public int sample(int i) {
       
        throw new RuntimeException("!!");
    }
    public void reverse(){
        
    }
   //[true, false, null] is order of state space
public double score(int obj_index, int i) {
    throw new RuntimeException("!!");
  
}  
    public void setRandom(double emiss, boolean restart){
        throw new RuntimeException("!!");
    }

   



public String toString(int i){
    return this.getName();
}

  
   
    public void validate() throws Exception {
        throw new RuntimeException("!!");
        
    }
 


    


    @Override
    public int noSnps() {
        // TODO Auto-generated method stub
        return 0;
    }

    public  void append(EmissionState emissionState){
        throw new RuntimeException("not implemented");
        }
    @Override
    public EmissionStateSpace getEmissionStateSpace() {
        // TODO Auto-generated method stub
        return null;
    }


    @Override
    public void print(PrintWriter pw, String st, List<Integer> columns) {
        // TODO Auto-generated method stub
        
    }


    @Override
    public double[] getEmiss(int i) {
        // TODO Auto-generated method stub
        return null;
    }


    @Override
    public Integer getFixedInteger(int i) {
        // TODO Auto-generated method stub
        return null;
    }


    


    @Override
    public int mostLikely(int pos) {
        // TODO Auto-generated method stub
        return -1;
    }
    @Override
    public int getParamIndex() {
        return 1;
    }

  
   
   
   
    
}
