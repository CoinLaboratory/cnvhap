package lc1.stats;

import java.util.ArrayList;
import java.util.List;

import lc1.util.Constants;


public class PermutationSampler extends Sampler {
    
   Dirichlet dir;
     public PermutationSampler(double[] dist1, double u){
         super(dist1, u);
         for(int i=0; i<dist.length; i++){
             l.add(i);
         }
         this.dir = new Dirichlet(dist1,u);
     
     }
     
     public PermutationSampler(Double[] dist1, double u){
         super(dist1, u);
         for(int i=0; i<dist.length; i++){
             l.add(i);
         }
         this.dir = new Dirichlet(dist1,u);
  
     }
     
     List<Integer> s = new ArrayList<Integer>();
     List<Integer> l = new ArrayList<Integer>();
    /* public Double[] sample(){
        // if(u==Double.POSITIVE_INFINITY) return dist;
         Double[] res = new Double[dist.length];
         Double[] dist = dir.sample();
         s.addAll(l);
         for(int i=0; i<res.length; i++){
             int j = s.remove(Constants.rand.nextInt(s.size()));
             res[i] = dist[j];
         }
         return res;
     }*/

     
    public  Double[] sample(){
       
        Double[] dist = dir.sample();
        Double[] res = (Double[]) dist.clone();
     //    int[] b = (int[])a.clone();
         for (int k = res.length - 1; k > 0; k--) {
             int w = Constants.rand.nextInt(k+1);//(int)Math.floor(Math.random() * (k+1));
             Double temp = res[w];
             if(temp>0 && res[k]>0){
                 res[w] = res[k];
                 res[k] = temp;
             }
         }
         return res;
         //printArray(b);
         }
   
}
