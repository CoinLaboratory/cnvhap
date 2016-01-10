/**
 * 
 */
package lc1.stats;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import lc1.dp.states.State;
import lc1.util.Constants;

public class SimpleDistribution implements Serializable {
   
    
    public final static transient SimpleDistribution noOffset = new SimpleDistribution(new int[] { 0 }, new double[] { 1 });
    public final static transient SimpleDistribution oneOffset = new SimpleDistribution(new int[] { -1, 0, 1 }, new double[] { 0.01, 0.98, 0.01 });
    public static transient double tolerance = 1e-3;
    
    static transient Double zero = new Double(0);
    
    public static Object sample(Iterator<Entry<Object, Double>> it ) {
        double pr = Constants.rand.nextDouble();
        double sum = 0;
       while( it.hasNext()){
            Entry<Object, Double> entry = it.next();
            sum+=entry.getValue();
            if(sum>=pr){
                return entry.getKey(); 
            }
        }
        throw new RuntimeException("did not find sample "+sum+" "+pr);
    }
    public HashMap<Object, Double> dist;
    public SimpleDistribution(){
        dist = new HashMap<Object, Double>();
    }
   
    public SimpleDistribution(int[] dist1, double[] d){
        this();
        for(int i=0; i<dist1.length; i++){
            put(dist1[i], d[i]);
        }
      //  validate();
    }
    
    public SimpleDistribution(Object[] o, double[] d){
        this();
        for(int i=0; i<o.length; i++){
            dist.put(o[i], d[i]);
        }
    }
    public SimpleDistribution(SimpleDistribution dist1){
        this(dist1, false, 0);
    }
    /*Applies offset first and finally inversion */
    SimpleDistribution(SimpleDistribution dist1, boolean invert, int offset){
        this();
         for(Iterator<Entry<Object,Double>> it = dist1.iterator(); it.hasNext();){
             Entry<Object, Double> entry = it.next();
             Object key =entry.getKey(); 
             if(offset!=0) key = ((Integer)entry.getKey()) - offset;
             if(invert) key = -1*((Integer)entry.getKey());
            dist.put(key, entry.getValue());
        }
    }
    public SimpleDistribution (SimpleDistribution left, SimpleDistribution right){
       this(left, true, 0);
        addAll(right);
    }
   
    
   public void addAll(SimpleDistribution dist){
    this.dist.putAll(dist.dist);
}
    public void addCount(Object key, double d){
            Double val = dist.get(key);
            dist.put(key,d + (val==null ? 0 : val));
        }
    //incr counts of this distribution for all equivlant things in the other distribution
    public void addCounts(SimpleDistribution distribution) {
        outer: for(Iterator<Entry<Object, Double>> it = distribution.dist.entrySet().iterator(); it.hasNext();){
            Entry<Object, Double> entry   = it.next();
            Double val = dist.get(entry.getKey());
            dist.put(entry.getKey(),entry.getValue() + (val==null ? 0 : val));
        }
    }
    
    public void addCounts(StateDistribution distribution, List<State> states) {
        outer: for(int j=0; j<states.size(); j++){
           double val1 = distribution.get(j);
           if(val1==0) continue;
           Double val = this.get(states.get(j));
            dist.put(states.get(j),val1 + (val==null ? 0 : val));
        }
    }
    public boolean different(SimpleDistribution d2) {
        for(Iterator it = this.dist.keySet().iterator(); it.hasNext();){
            Object key = it.next();
            double num1 = dist.get(key);
            double num2 = d2.dist.get(key);
            double diff = num1-num2;
            if(num1==0 || num2==0){
                return true;
            }
            if(Math.abs(diff)>0.001){
                return true;
            }
        }
        return false;
    }
    public double get(Object obj){
        Double d = dist.get(obj);
        return d==null ? 0 : d;
    }
    public Entry<Object, Double> getMax(){
        Iterator<Entry<Object, Double>> it = iterator();
        Entry<Object, Double> max = it.next();
        while( it.hasNext() ){
            Entry<Object, Double> entry = it.next();
            if(entry.getValue()>max.getValue()){
                max = entry;
            }
        }
        return max;
    }
    
    
    public Iterator<Entry<Object, Double>> iterator(){
           return  dist.entrySet().iterator();
    }
   
    
    public double KLDistance(SimpleDistribution d2){
        double sum =0;
        for(Iterator it = this.dist.keySet().iterator(); it.hasNext();){
            Object key = it.next();
            Double num1 = dist.get(key);
            Double num2 = d2.dist.get(key);
            if(num1==null) num1=zero;
            if(num2==null) num2=zero;
            if(num1!=0){
                sum+=num1* Math.log(num1 / num2); 
            }
        }
        return sum;
    }
    
   
    public int length(){
        return dist.size();
    }
    
    public Object mostLikely() {
       Iterator<Entry<Object, Double>> it = this.dist.entrySet().iterator();
       Entry<Object, Double> best = it.next();
       while(it.hasNext()){
           Entry<Object, Double> nxt = it.next();
           if(nxt.getValue()>best.getValue()){
               best = nxt;
           }
       }
       return best.getKey();
    }
    
    public void multiplyValues(double sum){
        for(Iterator<Entry<Object, Double>> it = iterator(); it.hasNext();){
            Entry entry = it.next();
            entry.setValue(((Number)entry.getValue()).doubleValue()*sum);
        }
    }
    
    
    
    
    public void normalise(){
        double sum= this.sum();
        if(sum==0) throw new RuntimeException("!!");
        this.normalise(sum);
    }
 
     public void normalise(double sum1){
        for(Iterator<Entry<Object,Double>> it = dist.entrySet().iterator(); it.hasNext();){
            Entry<Object, Double> entry = it.next();
            entry.setValue(entry.getValue()/sum1);
        }
    }
   
    public void put(Object obj, double d){
        dist.put(obj,d);
    }
   
    public Object sample(){
       return sample(dist.entrySet().iterator());
    }
    public void setRandom(double u, boolean restart){
        if(Math.abs(1.0-this.sum()) > 0.001) throw new RuntimeException("not valid");
        Double[] d = new Double[dist.values().size()] ;
        if(restart) Arrays.fill(d, 1.0/((double)d.length));
        else dist.values().toArray(d);
        Dirichlet dir = new Dirichlet(d,u);
        Double[] res = dir.sample();
        int i=0;
        for(Iterator<Entry<Object, Double>> it = iterator(); it.hasNext();i++){
            Entry <Object, Double>entry = it.next();
         //   if(entry.getValue()!=d[i]) throw new RuntimeException("!!");
            entry.setValue(res[i]);
        }
        this.normalise();
    }
    public int size(){
        return dist.size();
    }
    public double sum(){
        return sum(Object.class);
       
    }
    public double sum(Class clazz){
        double sum1=0;
        for(Iterator<Entry<Object, Double>> it = dist.entrySet().iterator(); it.hasNext();){
            Entry<Object, Double> entry = it.next();
            if(clazz.isInstance(entry.getKey())){
            //    if(Double.isNaN(entry.getValue())) throw new ArithmeticException("!!");
                sum1+=entry.getValue();
            }
         }
        return sum1;
    }
    
 public void removeSmallCounts(double d) {
     for(Iterator<Entry<Object, Double>> it = dist.entrySet().iterator(); it.hasNext();){
         Entry<Object, Double> entry = it.next();
         if(entry.getValue()<d) it.remove();
      }
        
    }
    public String toString(){
        StringBuffer sb = new StringBuffer("{");
        boolean first = true;
        for(Iterator<Entry<Object, Double>> it = iterator(); it.hasNext();){
            Entry<Object, Double> id = it.next();
            if(id.getValue()<1e-5) continue;
            if(!first) sb.append(",");
            else first = false;
            sb.append(id.getKey());
            sb.append("->");
            sb.append(Double.isNaN(id.getValue()) ? "NaN" : String.format("%5.3f", new Object[] {id.getValue()}));
        }
        sb.append("}");
        return sb.toString();
    }
    public void validate(){
        
         if(Math.abs(1.0 - sum()) > tolerance)throw new RuntimeException("sum is wrong "+Arrays.asList(dist)+" "+sum());
    }

   
  
}