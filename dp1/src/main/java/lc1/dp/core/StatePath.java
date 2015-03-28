package lc1.dp.core;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import lc1.dp.states.DotState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;



public class StatePath implements Comparable {

    public class StatePosEmission{
        public StatePosEmission(State st, Object em, int pos2, int j) {
            this.state = st;
            this.emission = em;
            this.pos = pos2;
            this.j = j;
        }
        public int j;
       public State state;
       public Integer pos;
       public  Object emission;
        public String toString(){
           String str = state.getName();
          String[] split  = str.split("_");
          if(split.length==3) return split[1]+split[2];
          else return str;
        }
    }
    public final List<StatePosEmission> spe;
    public List<State> getStates(){
        List<State>st = new ArrayList(this.size());
        for(int i=0; i<size(); i++){
            st.add(this.getState(i));
        }
        return st;
    }
    public StatePosEmission[] getSubList(Set<State> include){
        List<StatePosEmission> subl = new ArrayList<StatePosEmission>();
        for(int i=0; i<len(); i++){
            StatePosEmission spe_i = this.getSPE(i);
            if(include.contains(spe_i.state)) subl.add(spe_i);
        }
        return subl.toArray(new StatePosEmission[0]);
    }
   // public final List<Integer>emissionPos;
   // public final List emissions;
    final String name;
   // public final List<State> states;
    public double sc;
private final int[] posToIndex;
   public StatePath(StatePath sp){
        this.name = sp.name;
        this.sc = sp.sc;
        this.spe = new ArrayList<StatePosEmission>();
        this.posToIndex = new int[sp.posToIndex.length];
   }
   public StatePath(String name, double sc, int length){
    this.name = name;
    this.sc = sc;
    this.spe = new ArrayList<StatePosEmission>();
    this.posToIndex = new int[length];
    //emissionPos = new ArrayList<Integer>();
    //states = new ArrayList<State>();
}
   public StatePosEmission getStatePosForSeqIndex(int i){
	   return spe.get(posToIndex[i]);
   }
  
   
  // int len() =0;
   public void add(State st, Object em, int pos,int j){
	   this.posToIndex[pos] = spe.size();
       this.spe.add(new StatePosEmission(st, em, pos,j));
      
    //   emissions.add(em);
     //  emissionPos.add(pos);
      // states.add(st);
    //   len()++;
   }
   public int len(){
       return spe.size();
   }
   public StatePosEmission  getSPE(int i){
       return this.spe.get(len()-i-1);
   }
   public Object getEmission(int i){
       return getSPE(i).emission;
   }
   public Integer getIndex(int i){
       return getSPE(i).pos;
   }
   public State getState(int i){
       return getSPE(i).state;
   }
   public void setState(int k, State st) {
       getSPE(k).state = st;
        
    }
   
   public List<State> getSublist(Set<State> states){
       List<State> list = new ArrayList<State>();
       for(int i=0; i<this.len(); i++){
           State st = this.getState(i);
           if(states.contains(st)){
               list.add(st);
           }
       }
       return list;
   }
   
    public int[] find(EmissionState match1,EmissionState match2) {
        int[] res = new int[] {-1,-1};
        State match = match2;
       for(int i=0; i<spe.size(); i++){
           if(spe.get(i).state==match){
               if(match==match2){
                   res[1] = i;
                   match = match1;
               }
           else{
               res[0] = i;
               return res;
           }
           }
       }
       return null;
    }
  
    public int getEnd(){
        return this.getIndex(len()-1);
    }
        public Object getFirstEmission() {
            if(this.size()==0) return null;
           return this.getEmission(0);
        }
        public Integer getFirstEmissionPos(){
               if(this.spe.size()==0) return null;
               else return this.getIndex(0);
           }
        
        public State getFirstState(){
               if(this.spe.size()==0) return null;
               else return getState(0);
           }
 
    public Object getLastState() {
        if(this.spe.size()==0) return null;
        else return this.getState(len()-1);
    }
    public String getName(){
        return name;
    }
    
    public int getStart(){
        return getIndex(0);
    }
    
  
    
    public String[]getWord(){
        String[] sb = new String[spe.size()];
        for(int i=0; i<spe.size();i++){
            int i1 = spe.size()-i-1;
           sb[i1] = spe.get(i).state.getName(); 
        }
        return sb;
    }
   /* public String getWord() {
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<states.size();i++){
           sb.insert(0,states.get(i).getName()); 
        }
        return sb.toString();
    }*/
    
    public double percCorrect(State[] st){
       // if(st.len()gth!=states.size()) throw new RuntimeException("not right");
       double no_correct=0;
       double no_wrong=0;
       int k=0;
        for(int i=spe.size()-1; i>=0; i--){
            State st0 = spe.get(i).state;
            State st1 = st[k];
           if(st0 instanceof EmissionState){
               if(st0.equals(st1)) no_correct++;
               else no_wrong++;
               k++;
           }
        }
        return no_correct/ (no_wrong+no_correct); 
    }
    
    public int size(){
        return spe.size();
    }
    public String toString(){
        StringBuffer sb1 = new StringBuffer();
        //StringBuffer sb2 = new StringBuffer();
        //StringBuffer sb3 = new StringBuffer();
        for(int i=0; i<len(); i++){
            sb1.append(spe.get(len()-i-1));
          //  sb2.append(getIndex(i));
           // sb3.append(this.getEmission(i));
        }
        return //"state path\n"+emissions.reverse().toString()+"\n"+
          sb1.toString()+"\n";
        //+sb2.toString()+"\n"+sb3.toString()+"\n";
    }
    public int count(State state) {
        int count=0;
       for(int i=0; i<this.size(); i++){
          if(this.getState(i).equals(state)) count++;
       }
       return count;
    }
    public int compareTo(Object o) {
       StatePath sp1 = (StatePath) o;
       if(sc ==sp1.sc) return 0;
       else if(sc <sp1.sc) return 1;
       else return -1;
       
    }
    
    public int getIdenticalStretch(int start){
        State obj =this.getState(start);
        for(int i1 =start+1; i1 < this.size(); i1++){
            if(!this.getState(i1).equals(obj)) return i1;
        }
        return size();
    }
    public int getRecombineIndex(){
        int max_start =0;
        int max_end = 0;
        for(int start=0; start<this.size()-1; start++){
            int end = getIdenticalStretch(start);
            if(end-start > (max_end - max_start)){
                max_end = end;
                max_start = start;
            }
        }
        int max_index = max_start + (int) Math.floor((max_end - max_start)/2.0);
        return max_index;
    }
    public void merge(StatePath spi, DotState begin, DotState end) {
     //for(int i=0; i <spi.)
        throw new RuntimeException("need to reimplement this");
    }
   
   
  
}
