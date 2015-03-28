/**
 * 
 */
package lc1.dp;

import java.io.PrintWriter;
import java.util.Iterator;



class DiceEmissionState extends SimpleEmissionState{
    public  DiceEmissionState(String name, int adv, double[] dist, double pseudo){
        super(name, adv, getEmission(dist),pseudo);
    }
   
    static SimpleDistribution getEmission(double[] dist){
        SimpleDistribution d = new SimpleDistribution();
        for(int i=0; i<dist.length; i++){
            d.put(i, dist[i]);
        }
        return d;
    }

    public DiceEmissionState(DiceEmissionState st){
        super(st);
    }
    @Override
    public Object clone() {
      return new DiceEmissionState(this);
    }

    @Override
    public void print(PrintWriter pw, String prefix) {
        pw.print(prefix+" "+emissions.toString());
    }

  
  
    
   
    
    /*public double score(Object element, boolean logspace) {
       Integer val = (Integer) element;
       return logspace ? Math.log(dist[val]) : dist[val];
    }*/
}