package lc1.dp.states;

import java.io.PrintWriter;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.PseudoDistribution;

/** A class to use if there is an inconsistency in emission state spaces */
public class WrappedEmissionState1 extends EmissionState{
    
    //TODO - genotype consistent issue - 
    //what about [AA, _] vs [A,A] - will go through to inner state to deal with - should check this
 public EmissionState inner;

    public WrappedEmissionState1(EmissionState inner){
        super(inner.name, inner.adv);
       
        this.inner = inner;
      
    }
    @Override
    public PseudoDistribution emissions(int it){
    	return this.inner.emissions(it);
    }
    @Override
    public int noSnps() {
        return inner.noSnps();
    }
    public  void append(EmissionState emissionState){
        this.inner.append(emissionState);
         }
   
   /* @Override
    public void fillLikelihood(Locreader mid, List<Integer> loc) {
      //  System.err.println("class inner is "+inner.getClass());
            if(inner instanceof WrappedEmissionState){
                inner.fillLikelihood(mid, loc);
            }
            else if(!(inner instanceof HaplotypeEmissionState)){
              //  System.err.println(i+" "+this.getName()+" "+inner.getClass());
                   inner = new HaplotypeEmissionState(inner);
                   inner.fillLikelihood(mid, loc);
            }
            else inner.fillLikelihood(mid, loc);// =  new HaplotypeEmissionState(em[i], mid, loc[i].loc);
        
        
    }*/
    
    @Override
    public EmissionStateSpace getEmissionStateSpace() {
       return inner.getEmissionStateSpace();
    }
    @Override
    public void addCount(int obj_index,  double value, int i) {
    	/*if(inner instanceof CachedEmissionState){
    		CachedEmissionState inner1  = (CachedEmissionState)inner;
    		((CachedEmissionState) inner).innerState.addCount(obj_index, value, i);
    	}
    	else*/
    	 inner.addCount(obj_index,value,i);
    }
    
    public void refreshSiteEmissions() {
		inner.refreshSiteEmissions();
	}
    @Override
    public void addCountDT(double data_index,  int phen_index, double value, int i) {
       inner.addCountDT(data_index, phen_index, value, i);
        
    }
  
    @Override
    public String getUnderlyingData(int i){
       return this.inner.getUnderlyingData(i);
    }
    @Override
    public void print(PrintWriter pw, String st, List<Integer> columns) {
       inner.print(pw, st, columns);
        
    }
    @Override
    public void initialiseCounts() {
      inner.initialiseCounts();
        
    }
  
    @Override
    public int sample(int i) {
       throw new RuntimeException("");
    }
    @Override
    public double score(int obj_index, int i) {
      
       return this.inner.score(obj_index, i);
    }
    @Override
    public double[] getEmiss(int i) {
    	return inner.getEmiss(i);
      //  throw new RuntimeException("!!");
    }
    @Override
    public int mostLikely(int pos) {
        return inner.mostLikely(pos);
        // TODO Auto-generated method stub
    }
    @Override
    public boolean transferCountsToProbs(double pseudo) {
       return inner.transferCountsToProbs(pseudo);
        
    }
    public void modifyDirectCounts(double[] d){
    	if(inner instanceof CachedEmissionState){
    		if(Math.abs(d[0]-1.0)>0.1) ((CachedEmissionState)inner).modifyDirectCounts(d[0]);
    		if(Math.abs(d[1]-1.0)>0.1)((CachedEmissionState)inner).modifyInDirectCounts(d[1]);
    	}
    	else if(inner instanceof HaplotypeEmissionState){
    		
    	}
    	//else{
    }
    @Override
    public Integer getFixedInteger(int i) {
       return inner.getFixedInteger(i);
    }
    
    @Override
    public Object clone() {
      return new WrappedEmissionState1((EmissionState) this.inner.clone()
    		  );
    }
    
  
    public Integer noCop(){
    	return inner.noCop();
    }
    @Override
    public void validate() throws Exception {
     inner.validate();
        
    }
    @Override
    public void reverse() {
     this.inner.reverse(); 
    }
    @Override
    public int getParamIndex() {
        return this.inner.getParamIndex();
    }

    
    /*double[] prob1 = st.getEmiss(i);
    double[] prob = new double[emStSp.size()];
    Arrays.fill(prob, 0.0);
    for(int k=0; k<prob1.length; k++){
        String compa = ((ComparableArray)subSp.getHaploPair(k)).toStringShort().replaceAll("_", "");
        int genoIndex = emStSp.getFromString(compa);
        int[] indices = emStSp.getGenotypeConsistent( genoIndex);
        double[] weights = emStSp.getWeights(genoIndex);
        for(int k1=0; k1<indices.length; k1++){
             prob[indices[k1]] = weights[k1]* prob1[k];
        }
    }*/
}
