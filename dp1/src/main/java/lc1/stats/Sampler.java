package lc1.stats;

import java.io.Serializable;


public abstract class Sampler implements Serializable {
    
    
    
    public Double[] dist;
 
     
     public Sampler(double[] dist1, double u){
         this.dist = new Double[dist1.length];
         for(int i=0; i<dist.length; i++){
             dist[i] = dist1[i];//Math.max(1e-5, Math.min(1.0-1e-5, dist1[i]));
         }
     }
     public Sampler(float[] dist1, double u){
         this.dist = new Double[dist1.length];
         for(int i=0; i<dist.length; i++){
             dist[i] = (double) dist1[i];//Math.max(1e-5, Math.min(1.0-1e-5, dist1[i]));
         }
     }
     
     public Sampler(Double[] dist1, double u){
         this.dist = new Double[dist1.length];
         for(int i=0; i<dist.length; i++){
             dist[i] = dist1[i];
         }
      
     }
     
   
     
     
   
    
     public abstract Double[] sample();

   
}
