package lc1.dp.states;

import java.util.List;

public  class AlleleCopyPairEmissionState extends PairEmissionState {

  
    
   
    @Override
 public Integer calculateIndex(int i){
        EmissionState[] states = getMemberStates(false);
     
        Integer[] internal_indices = new Integer[states.length];
        for(int j=0; j<states.length; j++){
            if(   states[j] instanceof CompoundState){
                internal_indices[j] = ((CompoundState)states[j]).calculateIndex(i);
                if(internal_indices[j]==null){
                  throw new RuntimeException("!!");
                }
            }
            else{
              //  for(int i=0; i<internal_indices[j].length; i++){
                    internal_indices[j] = states[j].getFixedInteger(i);
                    if(internal_indices[j]==null){
                        return null;
                    }
               // }
            }
        }
       // for(int i=0; i<this.noSnps(); i++){
            int[] indices = new int[states.length];
        for(int j=0; j<states.length; j++){
          
           Integer index_i = internal_indices[j];
           if(index_i==null){
               return null;
           }
           indices[j] =  index_i;
        }
        Integer result = this.emStSp.getIndex(indices);
     if(result ==null){
         throw new RuntimeException("is null");
     }
      //  }
        return result;
    }
    
   
    
    
  
    //** assume they have the same emission space */
   public  AlleleCopyPairEmissionState(List<EmissionState> dist,  boolean decompose){
        super(dist,  decompose);
       //  if(this.emStSp.getMembers()[0].size()==1) throw new RuntimeException("!!");
        }
    
    
    
   public double score(int key, int i,  boolean recursive, boolean decompose){
       
            double sc = 1;
           int[] indices =  this.emStSp.getMemberIndices(key);
           if(indices[0]==0) sc = getInnerState(0, false).score(0, i);
           else{
            for(int j=0; j<indices.length; j++){
                EmissionState innerSt = getInnerState(j, false);
                    sc*=innerSt.score(indices[j], i);
            }
           }
            return sc;
    }
   public void addCount(int key,  Double value, int i, boolean decompose) {
      
          int[] indices =  this.emStSp.getMemberIndices(key);
          if(indices[0]==0) getInnerState(0, false).addCount(0, value, i);
          else{
              for(int j=0; j<indices.length; j++){
                    getInnerState(j, false).addCount(indices[j], value, i);
              }
          }
  }
    
    @Override
    public void addCount(int key1, double value, int i) {
        if(value==0) return;
        this.addCount(key1, value, i, decomp);
    }
    
    public double score(int key1, int i){
        Comparable key =  this.getEmissionStateSpace().get(key1);
        double sc = score(key1, i,  false, decomp);
        return sc;
    }
    
    public boolean transferCountsToProbs( double pseudo) {
       EmissionState[] stats = this.getMemberStates(true);
       boolean chnged = false;
       for(int i=0; i<stats.length; i++){
           chnged = chnged || ((HaplotypeEmissionState)stats[i]).transferCountsToProbs(pseudo);
           
       }
       return chnged;
    }
  
}
