package lc1.stats;

import java.io.Serializable;
import java.util.Arrays;

import lc1.dp.states.State;

public class StateDistribution implements Serializable {
   public  double[] dist;
  //  final State[] states;
    //final int numPairs;
    public StateDistribution(int modelLength){
        this.dist = new double[modelLength];
        Arrays.fill(dist, 0);
      //  numPairs  = (int) Math.round(((double)(dist.length-1))/2.0);
    }
    public StateDistribution(StateDistribution dist1){
        this.dist = new double[dist1.dist.length];
        System.arraycopy(dist1.dist, 0,dist, 0, dist.length );
       // numPairs  = (int) Math.round(((double)(dist.length-1))/2.0);
    }
    public String toString(){
        StringBuffer sb = new StringBuffer();
        Double[] d = new Double[dist.length];
        for(int i=0; i<d.length; i++){
            d[i] = dist[i];
            sb.append("%5.3g  ");
        }
        return String.format( sb.toString(),d);
    }
    
    public void multiplyValues(double sum){
        for(int j=0; j<this.dist.length; j++){
            dist[j]*=sum;
        }
    }
    
    /** returns original odd sum
    public StateDistribution[] split(){
        StateDistribution res = new StateDistribution(this.numPairs+1);
        StateDistribution res1 = new StateDistribution(3);
        res.dist[0] = this.dist[0];
        res1.dist[0] = this.dist[0];
        for(int j=0; j<this.numPairs; j++){
            res.dist[j+1] = dist[j*2+1] + dist[j*2+2];
            res1.dist[1] += dist[j*2+1];
            res1.dist[2] +=dist[j*2+2];
        }
        return new StateDistribution[]  {res, res1};
   } */
    public void addCounts(StateDistribution distribution) {
        if(distribution==null) return;
        if(distribution.dist.length!=this.dist.length) throw new RuntimeException("!!");
        for(int j=0; j<this.dist.length; j++){
        //outer: for(Iterator<Entry<Object, Double>> it = distribution.dist.entrySet().iterator(); it.hasNext();){
           // Entry<Object, Double> entry   = it.next();
   //         Double val = dist.get(entry.getKey());
            dist[j]+=distribution.dist[j];
     ///       dist.put(entry.getKey(),entry.getValue() + (val==null ? 0 : val));
        }
    }
    public void setLength(int length){
        double[] newdist = new double[length];
        System.arraycopy(dist, 0, newdist, 0, dist.length);
        this.dist = newdist;
    }
    public void put(State st, double d){
        dist[st.getIndex()] = d;
    }
    public void put(int st, double d){
        dist[st] = d;
    }
    public Double get(State state_j) {
       Double res =  dist[state_j.getIndex()];
       if(res==null) return zero;
       else return res;
    }
    public Double get(int j) {
        Double res =  dist[j];
        if(res==null) return zero;
        else return res;
     }
    
    public void validate(){
        if(Math.abs(1.0 - sum()) > tolerance){
            if(Math.abs(1.0 - sum()) > 0.001)
            throw new RuntimeException("sum is wrong "+Arrays.asList(dist)+" "+sum());
            else{
                this.multiplyValues(1.0/sum());
            }
        }
   }
    public int nonZero(){
        int nonZero =0;
        for(int i=0; i<dist.length; i++){
            if(dist[i]>0) nonZero++;
        }
        return nonZero;
    }
    public static double tolerance = 0.001;
    
    public double sum(){
        double sum=0;
        for(int i=0; i<dist.length; i++){
            Double res = dist[i];
            if(res!=null) sum+=res;
        }
        return sum;
    }
    public void normalise(){
        double sum = this.sum();
        for(int i=0; i<dist.length; i++){
            //if(dist[i]!=null)
                dist[i] = dist[i]/sum;
        }
    }
        public void setRandom(double u){
            if(Math.abs(1.0-this.sum()) > 0.001) throw new RuntimeException("not valid");
            Dirichlet dir = new Dirichlet(dist, u);
            Double[] res = dir.sample();
            for(int i=0; i<dist.length; i++){
                dist[i] = res[i];
            }
            this.normalise();
        }
        Double zero = 0.0;
        public double KLDistance(StateDistribution d2){
            double sum =0;
            for(int i=0; i<this.dist.length; i++){
                Double num1 = dist[i];
                Double num2 = d2.dist[i];
                if(num1==null) num1=zero;
                if(num2==null) num2=zero;
                if(num1!=0){
                    sum+=num1* Math.log(num1 / num2); 
                }
            }
            return sum;
        }
        public void addCount(Object key, Double value) {
          int ind  = ((Integer)key).intValue();
           //if(dist[ind]==null) dist[ind] = value;
           //else 
               dist[ind]+=value;
            
        }
          public int getMax(){
            int max =0;
       for  (int j=1; j<dist.length; j++){
           if(dist[j]>max){
                max = j;
            }
        }
        return max;
    }
        public synchronized void add(int j, double sc) {
            dist[j] += sc;
        }
        public void reset() {
            Arrays.fill(dist, 0);
            
        }
}
