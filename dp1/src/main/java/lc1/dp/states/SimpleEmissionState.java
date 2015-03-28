package lc1.dp.states;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;


public  class SimpleEmissionState extends EmissionState{
    
    protected PseudoDistribution emissions;
    
   
   
    public  SimpleEmissionState(String name, int adv, PseudoDistribution em){
        super(name, adv);
        this.emissions = new SimpleExtendedDistribution(em);
    }
    protected Collection getStateSpace() {
      throw new RuntimeException("!!");
    }

    public void addCount(int element,  double value, int i) {
        this.emissions.addCount(element,value);
    }
    public void addCountDT(double element, int phen_index,  double value, int i) {
       
    }
   
    
    public SimpleEmissionState(SimpleEmissionState st1){
        this(st1.name, st1.adv, st1.emissions);
    }
   public void reverse(){
       
   }

    public SimpleEmissionState(String string, int adv) {
      super(string, adv);
    }
    public void initialiseCounts(){
      this.emissions.initialise();
    }
    
   
   
     public int sample(int i){
       throw new RuntimeException("!!");
    }
     public double score(int obj, int i){
        return emissions.probs(obj);
    }
  
  
 

    
    
    public void setRandom(double u, boolean restart){
        throw new RuntimeException("!!");
    }
    public String getEmissionString(){
        return ""+ this.emissions;
    }
    public void transfer(){
        throw new RuntimeException("!!");
       
    }

    public void validate(){
        if(Math.abs(this.emissions.sum()-1.0) > 0.001) throw new RuntimeException("!!");
        this.lengthDistrib.validate();
    }
  
    @Override
    public Object clone() {
       return new SimpleEmissionState(this);
    }
    @Override
    public int noSnps() {
        return 0;
    }
    @Override
    public void print(PrintWriter pw, String st, List<Integer> columns) {
       throw new RuntimeException("");
        
    }
    @Override
    public EmissionStateSpace getEmissionStateSpace() {
        // TODO Auto-generated method stub
        return null;
    }
    @Override
    public double[] getEmiss(int i) {
        // TODO Auto-generated method stub
        return null;
    }
    @Override
    public int mostLikely(int pos) {
        // TODO Auto-generated method stub
        return -1;
    }
    @Override
    public Integer getFixedInteger(int i) {
        // TODO Auto-generated method stub
      return this.emissions.fixedInteger();
    }
    
    @Override
    public boolean  transferCountsToProbs(double pseudo) {
     throw new RuntimeException("!!");
        
    }
    @Override
    public int getParamIndex() {
       return 1;
    }
    
}
