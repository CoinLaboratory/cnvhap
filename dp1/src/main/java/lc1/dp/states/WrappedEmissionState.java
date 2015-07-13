package lc1.dp.states;

import java.io.PrintWriter;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;

/** A class to use if there is an inconsistency in emission state spaces */
public class WrappedEmissionState extends EmissionState{
    
    //TODO - genotype consistent issue - 
    //what about [AA, _] vs [A,A] - will go through to inner state to deal with - should check this
 EmissionState inner;
final EmissionStateSpace emStSp;
final EmissionStateSpaceTranslation trans;
    public WrappedEmissionState(EmissionState inner, EmissionStateSpace target, EmissionStateSpaceTranslation trans){
        super(inner.name, inner.adv);
        this.trans =trans;
        this.inner = inner;
        this.emStSp = target;
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
       return emStSp;
    }
    @Override
    public void addCount(int obj_index,  double value, int i) {
       inner.addCount(trans.bigToSmall(obj_index),value, i);
        
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
       Integer sm = trans.bigToSmall(obj_index);
       if(sm==null) return 0;
       return this.inner.score(sm, i);
    }
    @Override
    public double[] getEmiss(int i) {
//    if(true)    throw new RuntimeException("");
      //  Logger.global.warning("is inefficient to repeatedly call this!");
        double[] d =new  double[this.emStSp.size()];
        for(int j=0; j<d.length; j++){
            d[j] = score(j, i);
        }
        return d;
      //  throw new RuntimeException("!!");
    }
    @Override
    public int mostLikely(int pos) {
        return trans.smallToBig(pos)[0];
        // TODO Auto-generated method stub
    }
    @Override
    public boolean transferCountsToProbs(double pseudo) {
       return inner.transferCountsToProbs(pseudo);
        
    }
    @Override
    public Integer getFixedInteger(int i) {
        Integer i1 = inner.getFixedInteger(i);
        if(i1==null) return null;
        return trans.smallToBig(i1)[0];
    }
    
    @Override
    public Object clone() {
      return new WrappedEmissionState((EmissionState) this.inner.clone(), this.emStSp, this.trans);
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
