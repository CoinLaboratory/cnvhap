package lc1.dp.states;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import lc1.util.Constants;

public class MergedEmissionState extends HaplotypeEmissionState {

  public   EmissionState[] em;
    
   private  List<int[]> transformation ; //first is index of state, second is index within list
    public void reverse(){
      throw new RuntimeException("!!");
    }
    static String getName(EmissionState[] em){
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<em.length; i++){
            sb.append(em[i].getName());
        }
        return sb.toString();
    }
    /* main is the index to choose in case of conflicts */
    public static void getTransformation(List<Integer>[]locs, 
            List<String> []snp_id,
            List<int[]> transformation, 
            List<Integer> overallLocs, int main){
       transformation.clear();
        overallLocs.clear();
        SortedSet<PosTrans> l = new TreeSet<PosTrans>();
        for(int i=0; i<locs.length; i++){
          for(int j=0; j<locs[i].size(); j++){
             PosTrans pt = new PosTrans(locs[i].get(j), i, j, snp_id[i].get(j));
             if(l.contains(pt) ){
                 if(i==main){
                     PosTrans first = l.tailSet(pt).first();
                     first.replace(pt);
                 }
             }
             else{
                 l.add(pt);
             }
          }
        }
       // Collections.sort(l);
        for(Iterator<PosTrans> it = l.iterator(); it.hasNext();){
            PosTrans nxt = it.next();
            transformation.add(nxt.ind_loc);
            overallLocs.add(nxt.pos);
        }
    }
    
    public static void adjustCoords(List<Integer>[]locs, 
            List<String> []snp_id
          ){
        SortedSet<PosTrans1> l = new TreeSet<PosTrans1>();
        for(int i=0; i<locs.length; i++){
          for(int j=0; j<locs[i].size(); j++){
             PosTrans1 pt = new PosTrans1(locs[i].get(j), i, j, snp_id[i].get(j));
             if(l.contains(pt) ){
                 PosTrans1 first = l.tailSet(pt).first();
                 if(!first.snp_id.equals(pt.snp_id)){
                     int incr = locs[i].get(j)+1;
                     if(locs[i].size() > j+1 && locs[i].get(j+1)==incr) throw new RuntimeException("!!");
                     locs[i].set(j, incr);
                      pt =  new PosTrans1(locs[i].get(j), i, j, snp_id[i].get(j));
                     if(l.contains(pt)) throw new RuntimeException("!!");
                     l.add(pt);
                 }
             }
             else{
                 l.add(pt);
             }
          }
        }
       
    }
   
    public MergedEmissionState(HaplotypeEmissionState[] em, List<int [] > transformation){
        super(em[0].getName(), transformation.size(), em[0].getEmissionStateSpace(), (short)-1);
        if(Constants.CHECK){
            for(int i=1; i<em.length; i++){
                if(em[0].getEmissionStateSpace()!=em[i].getEmissionStateSpace()) throw new RuntimeException("!!");
              //  if(em[0].isFixed()!=em[i].isFixed()) throw new RuntimeException("!!");
            }
        }
        this.em = em;
        this.transformation = transformation;
        for(int i=0; i<this.noSnps; i++){
            int[] inde = transformation.get(i);
            this.emissions[i] = em[inde[0]].emissions[inde[1]];
        }
    }
    
    public MergedEmissionState(MergedEmissionState state) {
        super(state);
        this.em = new EmissionState[state.em.length];
        for(int i=0; i<em.length; i++){
            em[i] =(EmissionState) state.em[i].clone();
        }
        this.transformation = new ArrayList<int[]>(state.transformation);
    }

    

    
   /* @Override
    public boolean transferCountsToProbs(double pseudo) {
    // if(true) throw new RuntimeException("!!");
        boolean chnged = false;
       for(int i=0; i<em.length; i++){
           if(em[i].transferCountsToProbs(pseudo)){
               chnged = true;
           }
       }
       if(chnged) paramIndex++;
       return chnged;
        
    }*/

    @Override
    public Object clone() {
       return new MergedEmissionState(this);
    }

    int paramIndex = 1;
    @Override
    public int getParamIndex() {
        return paramIndex;
    }
  

}
