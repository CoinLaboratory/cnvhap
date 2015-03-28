package lc1.dp.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;

import lc1.dp.data.representation.Emiss;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.dp.transition.BetweenWithinTransitionProbs1;
import lc1.dp.transition.BetweenWithinTransitionProbs3;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.stats.Dirichlet;
import lc1.stats.PermutationSampler;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;

public class FreeExpSiteTrans extends SiteTransitions {
     final List<State> states;
public FreeExpSiteTrans(List<Integer> loc, List<State> states, Double[] exp_p1, Double[] r, int length, Object[] transProbtype) {
        super(loc, states.size()-1, exp_p1, r, length,0, 
        		SiteTransitions.getCN(states), true);
        this.states = states;
        
        this.transProbType = transProbtype;
        // TODO Auto-generated constructor stub
    }

final Object[] transProbType;
   
/** if pseudocount=0 we do not use this as pseudo */
public FreeExpSiteTrans(FreeExpSiteTrans trans_init, boolean swtch){
    super(trans_init, swtch);
    this.states = trans_init.states;
    this.transProbType = trans_init.transProbType;
//    this(trans_init.loc,trans_init.states, trans_init.exp_p1, trans_init.r,  trans_init.transProbs.length,trans_init.transProbType);
}
double[] u = null;
@Override
public SiteTransitions clone(boolean swtch) {
    return new FreeExpSiteTrans(this, swtch);
}

    private static final int[] group(List<State> states, int[] stateToGroup, int[] stateToIndexWithinGroup){
    	int[] group  = Constants.group();
    	if(group==null) return group_1(states, stateToGroup, stateToIndexWithinGroup);
    	SortedMap<Integer, Integer> cnt = new TreeMap<Integer, Integer>();
    	stateToGroup[0] = 0;
    	stateToIndexWithinGroup[0] =0;
    	for(int i=0; i<group.length; i++){
    		Integer count = cnt.get(group[i]);
    		if(count==null){
    			count=0;
    		}
    		stateToIndexWithinGroup[i+1] = count;
    		stateToGroup[i+1] = group[i];
    		cnt.put(group[i], count+1);
    		
    	}
    	int[] res = new int[cnt.lastKey()+1];
    	Arrays.fill(res,-1);
    	res[0] = 1;
    	for(Iterator<Entry<Integer, Integer>> it = cnt.entrySet().iterator(); it.hasNext();){
    		Entry<Integer, Integer> nxt = it.next();
    		res[nxt.getKey()] = nxt.getValue();
    	}
    	return res;
    }
    private static int[] group_1(List<State> states, int[] stateToGroup, int[] stateToIndexWithinGroup) {
        Map<String, List<Integer>> map = new TreeMap<String, List<Integer>>();
        map.put(-1+"", Arrays.asList(new Integer[] {0}));
        Set<Integer> commonClass = new HashSet<Integer>();
        for(int i=0; i<Constants.cnStatesInCommonClass().length; i++){
            commonClass.add(Constants.cnStatesInCommonClass()[i]);
        }
        for(int i=1; i<states.size(); i++){
            EmissionState st = (EmissionState) states.get(i);
            int best_index = st.getBestIndex(0);
            Emiss em = (Emiss) st.getEmissionStateSpace().get(best_index);
          //  String str = em.toStringShort();
            int num = (em).noCopies();
            //    Constants.getMax(st.emissions[0].probs);
           for(int j=1; j<st.noSnps(); j++){
               int num1 = ((Emiss) st.getEmissionStateSpace().get(st.getBestIndex(j))).noCopies();
               if(num1!=num) throw new RuntimeException("!!");
           }
        
           String key = num+"";
           if(num!=1 && !commonClass.contains(num)){
               key= key+"_"+i;
           }
            List<Integer> l = map.get(key);
            if(l==null){
                map.put(key, l = new ArrayList<Integer>());
            }
            l.add(i);
        }
        List<Integer>[] membs=map.values().toArray(new List[0]);
        for(int j=0; j<membs.length; j++){
            for(int i=0; i<membs[j].size(); i++){
                int k = membs[j].get(i);
                stateToGroup[k] = j;
                stateToIndexWithinGroup[k] = i;
            }
        }
        int[] res = new int[membs.length];
        for(int i=0; i<res.length; i++){
            res[i] = membs[i].size();
        }
        return res;
    }
    
    
  
 /*   private String getOutString(FreeTransitionProbs1 probs1) {
        if(Constants.CHECK && ( probs1.transitionsOut[1]!=null)) throw new RuntimeException("!!");
        double[] counts =  probs1.transitionsOut[0].counts;
      StringBuffer sb = new StringBuffer();
        for(int i=0; i<counts.length; i++){
            sb.append(counts[i]+"\t");
        }
        return sb.toString();
    }*/
   // boolean permute = true;
    
 
    /** dist is initial distribution over states, used to build between group transitions */
    @Override
    public void initialise(double[] dist, double permute, double u_g) throws Exception{
            int[] stateToGroup = new int[states.size()];
            int[] stateToIndexWithinGroup = new int[states.size()];
            int[] groupSizes = group(states, stateToGroup, stateToIndexWithinGroup);
            int[][] groupToState = new int[groupSizes.length][];
            for(int i=0; i<groupSizes.length; i++){
            	groupToState[i] = new int[groupSizes[i]];
            	
            }
            for(int i=0; i<stateToGroup.length; i++){
            	groupToState[stateToGroup[i]][stateToIndexWithinGroup[i]] = i;
            }
       //  if(alpha_overall.length!=groupSizes.length) throw new RuntimeException("!!");
          
           int[][] statesToGroupTrans = new int[stateToGroup.length][];
           int[][]  statesToWithinGroupTrans = new int[groupToState.length][];
           FreeSiteTrans1.makeStateToGroupTrans(stateToGroup,stateToIndexWithinGroup, 
            		groupToState, Constants.transitions(0), 
            		statesToGroupTrans, statesToWithinGroupTrans);  
           
                Double[] d = new Double[groupSizes.length];
                Arrays.fill(d, 0.0);
                for(int i=0; i<states.size(); i++){
                    int g = stateToGroup[i];
                    d[g]+=dist[i];
                }
                
                //for(int k=1; k<d.length; k++){
                 //   d[k] = (double) groupSizes[k] / (double) numFounders;
               // }
                double[] u0 = Constants.u_global(0);
            Dirichlet   dir= new Dirichlet(d, Constants.u_global(0)[1]);
            Sampler[] samplers = new Sampler[groupSizes.length];  // within group alpha for transitions in
            for(int k=0; k<groupSizes.length; k++){
                Double[] d1 = new Double[groupSizes[k]];
                double incr = 1.0/(double)groupSizes[k];
                if(Constants.samplePermute()>0 ){
                    int j=0;
                    for(j=0; j<d1.length; j++){
                        d1[j] = //(double) j+1; //
                        Math.pow((double)Constants.samplePermute(), j);
                    }
                //   d1[j-1] = d1[j-1]*Constants.samplePermute();
                    SimpleExtendedDistribution.normalise(d1);
                    samplers[k] = new PermutationSampler(d1, Constants.u_global(1)[1]);
                }
                else{
                    Arrays.fill(d1, incr);
                    samplers[k] = new Dirichlet(d1, Constants.u_global(1)[1]);
                }
            }
            double[] start = new double[states.size()];
            start[0] = 1.0;
         double[] probs = new double[states.size()];
            //Double[] r_0 = new Double[] {1e-60,r[1]};
            transProbs[0] = new BetweenWithinTransitionProbs1(dir, samplers, stateToGroup, stateToIndexWithinGroup, new Dirichlet[] {null, null}, new Object[] {
            	FreeTransitionProbs1.class, FreeTransitionProbs1.class	, (Class[])transProbType[2]
            },samplers.length, groupToState);
            FreeSiteTrans1.fillProbs(transProbs[0], probs, start);
            for(int i=1; i<transProbs.length; i++){
              //  Double[] r_ = special==null || this.special.contains(i) ?  r : r_0;
                double[] expp = 
                    
                   
                    new double[] {
                        loc==null || loc.size()==0 ? exp_p1[0] : Math.exp(-1*r[0][0]*(loc.get(i)-loc.get(i-1))),
                        loc==null || loc.size()==0 ? exp_p1[1] : Math.exp(-1*r[1][0]*(loc.get(i)-loc.get(i-1)))
                    };
              /* if(special!=null && !special.contains(i)){
                   expp[0] = 1.0;  // 0 is between groups
               }*/
             
                Dirichlet[] exp_p = new Dirichlet[] {
                        new Dirichlet(new double[] {expp[0], 1-expp[0]}, Constants.u_global(0)[2]),
                        new Dirichlet(new double[] {expp[1], 1-expp[1]}, Constants.u_global(1)[2]),
                }; //exp_p[0] is between groups
                
              if(transProbs[i]==null)  {transProbs[i] = 
                    Constants.trans1() ? 
                    new BetweenWithinTransitionProbs3(dir, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState):
                        new BetweenWithinTransitionProbs1(dir, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType,samplers.length, groupToState);
              }
                        else{
                        	transProbs[i] = 
                                Constants.trans1() ? 
                                new BetweenWithinTransitionProbs3(probs,transProbs[i], null, samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState,
                                		//statesToGroupTrans, statesToWithinGroupTrans
                                		i):
                                    new BetweenWithinTransitionProbs1(transProbs[i],null,  samplers, stateToGroup, stateToIndexWithinGroup, exp_p, transProbType, groupToState
                                    		//statesToGroupTrans, statesToWithinGroupTrans
                                    		);
                        }
              System.arraycopy(probs,0, start,0,start.length);
              FreeSiteTrans1.fillProbs(transProbs[i], probs, start);
                if(Constants.CHECK){
                    validate(transProbs[i], states.size(), i);
                }
            }
       
    }
   
  
}
