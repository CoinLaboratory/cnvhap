package lc1.dp.states;

import java.util.ArrayList;
import java.util.List;

import lc1.dp.data.representation.ComparableArray;
import lc1.stats.ProbabilityDistribution;



public abstract class CompoundState extends EmissionState{
    public CompoundState(String name, int adv) {
        super(name, adv);
        // TODO Auto-generated constructor stub
    }
    /* (non-Javadoc)
     * @see lc1.dp.states.CompState#reverse()
     */
    public void reverse(){
        EmissionState[] emStates = getMemberStates(true);
        for(int i=0; i<emStates.length; i++){
            emStates[i].reverse();
        }
    }
    /* (non-Javadoc)
     * @see lc1.dp.states.CompState#isFixed()
     */
   
   protected boolean isFixed(int k) {
        EmissionState[] emStates = getMemberStates(true);
        for(int i=0; i<emStates.length; i++){
            if(emStates[i].getFixedInteger(k)==null) return false;
        }
        return true;
    }
    /*
     * @parameter real - do we want the real underlying member states or the 'equivalent' one.  They are the same, except for
     * HalfTrioEmissionState and TrioEmissionState
     */
    /* (non-Javadoc)
     * @see lc1.dp.states.CompState#getMemberStates(boolean)
     */
    public abstract EmissionState[] getMemberStates(boolean real);
   /* (non-Javadoc)
 * @see lc1.dp.states.CompState#getMemberStatesRecursive(boolean)
 */
public  ComparableArray getMemberStatesRecursive(boolean real){
        EmissionState[] membs = getMemberStates(real);
        List<Comparable> res = new ArrayList<Comparable>();//(membs[0] instanceof CompoundState);
        for(int i=0; i<membs.length; i++){
          /*  if(membs[i] instanceof CompoundState){
                res.add(((CompoundState)membs[i]).getMemberStatesRecursive(real));
            }
            else{*/
                res.add(membs[i]);
            //}
        }
        return new ComparableArray(res);
    }
    public CompoundState(CompoundState st1){
        this(st1.name, st1.adv);
    }
   // public abstract Integer[] calculateIndex();
    /* (non-Javadoc)
     * @see lc1.dp.states.CompState#calculateIndex(int)
     */
    public abstract Integer calculateIndex(int i);

 //   public abstract double score(ComparableArray comp_a, int i,  boolean recursive, boolean decompose);
    public abstract double score(int j, int i, boolean recursive, boolean decompose) ;
       
     /* (non-Javadoc)
     * @see lc1.dp.states.CompState#noSnps()
     */
    public int noSnps(){
         return getMemberStates(false)[0].noSnps();
     }
    /* (non-Javadoc)
     * @see lc1.dp.states.CompState#memberStatesAreCompound()
     */
    public int memberStatesAreCompound(){
         EmissionState[] memberStates = getMemberStates(false);
         int cntMoreThan1 =0;
         for(int i=0; i<memberStates.length; i++){
             if(memberStates[i] instanceof CompoundState 
                     && ((CompoundState)memberStates[i]).getMemberStates(false).length>1){
                 cntMoreThan1++;
             }
         }
         return cntMoreThan1;
     }
  //public abstract void initialise(StateIndices dat) ;
  
    ProbabilityDistribution[] tmp;
    public ProbabilityDistribution calcAverageDistributions(int k, int i) {
      
       EmissionState[] st = this.getMemberStates(true);
       if(tmp==null) tmp = new ProbabilityDistribution[st.length];
       for(int i1 =0; i1<tmp.length; i1++){
           tmp[i1] = ((HaplotypeEmissionState)st[i1]).emissionsDatatype[k][i];
       }
      return tmp[0].clone();
     
    }
    
    public void setAverageDistributions(ProbabilityDistribution res, int k, int i) {
        
        EmissionState[] st = this.getMemberStates(true);
        if(tmp==null) tmp = new ProbabilityDistribution[st.length];
        for(int i1 =0; i1<tmp.length; i1++){
            tmp[i1] = ((HaplotypeEmissionState)st[i1]).emissionsDatatype[k][i];
        }
        res.setParamsAsAverageOf(tmp);
      
     }
    
	public abstract void modifyDirectCounts(double d) ;
	
  
    
   
}
