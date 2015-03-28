package data;

import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;
public class Dirichlet {
    
   
   final  double u;
   GammaDistributionImpl[] g;
     final double[] dist1;
    
     public Dirichlet(double[] dist1, double u){
    	 double sum=0;
    	 for(int k=0; k<dist1.length; k++){
    		 sum+=dist1[k];
    	 }
       this.dist1 = dist1;
       for(int k=0; k<dist1.length; k++){
  		 this.dist1[k] = dist1[k]/sum;
  	 }
         this.u = u;
         this.g = new GammaDistributionImpl[dist1.length];
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
             if(dist1[i]>0){
                 if(g[i]!=null){
                	 g[i].setAlpha(dist1[i]*u);
                	 g[i].setBeta(1);
                 }
                 g[i] = new GammaDistributionImpl(dist1[i]*u, 1) ;
             }
         }
     }
     
     
   // RandomElement re = new RandomElement(){
        /* public double raw(){
             return Constants.rand.nextDouble();
         }
     };*/
     public double[] sample() throws MathException{
         if(u==Double.POSITIVE_INFINITY) return dist1;
         double[] res = new double[dist1.length];
         double sum=0;
         for(int i=0; i<res.length; i++){
        	 double d = Math.random();
        	 double d1 = g[i]==null ? 0 : g[i].inverseCumulativeProbability(d);
             res[i] = g[i]==null ? 0 : Math.max(1e-5,d1);
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
