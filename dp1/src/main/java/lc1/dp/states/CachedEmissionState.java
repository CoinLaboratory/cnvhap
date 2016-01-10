package lc1.dp.states;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.IntegerDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;


public class CachedEmissionState extends AbstractCachedEmissionState {

    public CompoundState innerState;
   
   public  PseudoDistribution[] emissions;
 //  final double[] pseudo;
   @Override
   public void reverse(){
       List<PseudoDistribution> l =Arrays.asList(emissions);
       Collections.reverse(l);
       emissions = l.toArray(new PseudoDistribution[0]);
   }
   @Override
   public double[] getEmiss(int i) {
       return this.emissions[i].probs();
   }
   public Integer noCop(){
	   return innerState.noCop();
   }
  
   public  void append(EmissionState emissionState){
   	int len = emissions.length;
   	PseudoDistribution[] ems = new PseudoDistribution[this.emissions.length+emissionState.length()];
   	System.arraycopy(this.emissions, 0, ems, 0, len);
   	
   	for(int k=0 ; k<emissionState.length(); k++){
   		ems[len+k] = 
   		((CachedEmissionState)emissionState).emissions(k);
   	}
   	this.emissions =ems;
   }
   
  
   public CachedEmissionState(CompoundState state,  int stateSpaceSize) {
        super(state);
        int noSites = state.noSnps();
        ProbabilityDistribution[][] dists = null;//((HaplotypeEmissionState)state.getMemberStates(true)[0]).emissionsDatatype;
       // int[] numLevels = Constants.numLevels();
        this.emissions = new PseudoDistribution[noSites];
       
        this.innerState = state;
        for(int i=0; i<noSites; i++){
           Integer ind = this.calculateIndex(i);
           if(ind!=null){
               emissions[i] = new IntegerDistribution(ind, this.getEmissionStateSpace());   
           }
           else
               emissions[i] = new SimpleExtendedDistribution(stateSpaceSize);
       }
       refreshSiteEmissions();
    }
	public PseudoDistribution emissions(int i) {
		return this.emissions[i];
	}

 
   public EmissionStateSpace getEmissionStateSpace(){
       return innerState.getEmissionStateSpace();
   }
   
   public boolean transferCountsToProbs( double pseudo) {
	   this.transferCountsToMemberStates();
	   return true;
//       throw new RuntimeException("!!! "+this.getClass());
   }
   
   public void modifyDirectCounts(double d) {
		innerState.modifyDirectCounts(d);
		
	}
   
   public void modifyInDirectCounts(double d) {
	   PseudoDistribution[] dists = this.emissions;
		for(int i1=0; i1<dists.length; i1++){
			dists[i1].multiplyCounts(d);
		}
		
	}
  
   
   /* transfers the counts to the member states which make up this compound state */
   public void transferCountsToMemberStates(){
      for(int i=0; i<emissions.length; i++){
          if(emissions[i].fixedInteger()==null){
              double[] counts = emissions[i].counts();
             
              for(int k=0; k<counts.length; k++){
                  if(counts[k]==0) continue;
                      innerState.addCount(k, counts[k], i);
              }
             
          }
          else{
        	  innerState.addCount(emissions[i].fixedInteger(), ((IntegerDistribution)emissions[i]).cnt, i);
          }
      }
      if(innerState instanceof CachedEmissionState){
          ((CachedEmissionState)innerState).transferCountsToMemberStates();
      }
   }
   

   @Override
   public double score(int object_index, int i1){
       //int i= emissions.length==1 ? 0 : i1;
      return emissions[i1].probs(object_index);
       //return probs[object_index];
   }
   
   
   
    
   public Object clone() {
       return new CachedEmissionState((CompoundState) innerState.clone(),emissions[0].probs().length);
   }
   
   @Override
   public void initialiseCounts(){
       for(int i=0; i<this.emissions.length; i++){
           this.emissions[i].initialise();
       }
   }
   
  /*@Override 
  public int mod(int j, int di){
	  return innerState.mod(j, di);
  }*/
   
 /*  @Override
   public double score(ComparableArray comp_a, int i,  boolean recursive, boolean decompose) {
      return this.innerState.score(comp_a, i, recursive, decompose);
   }*/
   @Override
   public double score(int j, int i,  boolean recursive, boolean decompose) {
      return this.innerState.score(j, i, recursive, decompose);
   }
   
   @Override
   public void addCount(int obj_index,  double value, int i1) {
       int i= emissions.length==1 ? 0 : i1;
    //   if(this.getName().equals("1_1c")){
    //	   System.err.println("p")
;  //     }
       this.emissions[i].addCount(obj_index, value);
   }
   
    
   
int paramIndex = 1;
 
    
    public void refreshSiteEmissions() {
        if(innerState instanceof CachedEmissionState){
            ((CachedEmissionState)innerState).refreshSiteEmissions();
        }
        this.paramIndex = innerState.getParamIndex();
        int noSites = this.emissions.length;
        for(int i=0; i<noSites; i++){
            PseudoDistribution em = emissions[i];
           if( em instanceof IntegerDistribution) continue;
           {
               
                   double[] probs = emissions[i].probs();
                   boolean fixed = false;
                   for(int j=0; j<probs.length; j++){
                       double sc = innerState.score(j,i);
                       if(sc > 0.999){
                           fixed= true;
                         //  break;
                       }
                       probs[j]=sc;
                   }
                   if(fixed){
                       Integer fixedIndex = this.calculateIndex(i);
                       if(fixedIndex!=null){
                           emissions[i] = new IntegerDistribution(fixedIndex, this.getEmissionStateSpace());
                       }
                   }
                      
           }
          if(Constants.CHECK && Math.abs(1-emissions[i].sum())>0.001) {
              try{
              innerState.validate();
              EmissionStateSpace emStSp = innerState.getEmissionStateSpace();
               throw new RuntimeException("!! "+emissions[i].sum()+" "+this.name);
              }catch(Exception exc){
                  exc.printStackTrace();
                  System.exit(0);
              }
         }
         
        }
    }

   

   
   

    @Override
    public int mostLikely(int pos) {
        return this.innerState.mostLikely(pos);
    }

    @Override
    public void print(PrintWriter pw, String st, List<Integer>columns) {
        this.innerState.print(pw, st, columns);
    }

    @Override
    public int sample(int i) {
        throw new RuntimeException("!!");
    }

 

    @Override
    public void validate() throws Exception {
        this.innerState.validate();
    }

    @Override
    public EmissionState[] getMemberStates(boolean real) {
      return innerState.getMemberStates(real);
    }

public EmissionState innerState(){
	return innerState;
}
    /*  @Override
   public void initialise(StateIndices dat) {
       this.innerState.initialise(dat);
        
    }*/

    @Override
    public Integer calculateIndex(int i) {
      return this.innerState.calculateIndex(i);
    }

    @Override
    public Integer getFixedInteger(int i) {
      return this.emissions[i].fixedInteger();
    }
    @Override
    public int getParamIndex() {
       return this.paramIndex;
    }
	

    
   

    
   
  
}
