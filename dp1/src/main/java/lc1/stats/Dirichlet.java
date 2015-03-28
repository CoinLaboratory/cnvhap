package lc1.stats;

import java.util.Arrays;

import lc1.util.Constants;
import cern.jet.random.Gamma;
import edu.cornell.lassp.houle.RngPack.RandomElement;


public class Dirichlet extends Sampler {
    
    public static void main(String[] args){
        Double[] dist = new Double[] {1.0/3.0, 1.0/3.0, 1.0/3.0};
        Dirichlet d = new Dirichlet(dist, 1);
        for(int i=0; i<10; i++){
            Double[] res = d.sample();
            System.err.println(Arrays.asList(res));
        }
       
    }
    
   final  double u;
     Gamma[] g;
     
     public Dirichlet(double[] dist1, double u){
         super(dist1, u);
         this.u = u;
         this.g = new Gamma[dist.length];
         set(u);
     }
     
     public Dirichlet(Double[] dist1, double u){
         super(dist1, u);
         this.u = u;
         this.g = new Gamma[dist.length];
         set(u);
     }
     
     public Dirichlet(int length, double u) {
		this(getUniform(length),u);
	}

	private static double[] getUniform(int length) {
		double[] res = new double[length];
		Arrays.fill(res,1.0/(double)length);
		return res;
	}

	public void set(double u){
         if(u==Double.POSITIVE_INFINITY){
       //      throw new RuntimeException("should not set u to pos inf");
            // return;
         }
         for(int i=0; i<g.length; i++){
             if(dist[i]>0){
                 if(g[i]!=null) g[i].setState(dist[i]*u, 1);
                 g[i] = new Gamma(dist[i]*u, 1, re) ;
             }
         }
     }
     
     
    RandomElement re = new RandomElement(){
         public double raw(){
             return Constants.rand.nextDouble();
         }
     };
     public Double[] sample(){
         if(u==Double.POSITIVE_INFINITY) return dist;
         Double[] res = new Double[dist.length];
         double sum=0;
         for(int i=0; i<res.length; i++){
             res[i] = g[i]==null ? 0 : Math.max(1e-5,g[i].nextDouble());
             sum+=res[i];
         }
         for(int i=0; i<res.length; i++){
             res[i] = res[i]/sum;
         }
         return res;
     }
   

    public double u() {
       return u;
    }
}
