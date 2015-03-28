package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.Dirichlet;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

public class MultiExpProbs  implements ExpTransProb{
    
   private SimpleExtendedDistribution1[] exp_rd;
   public  MultiExpProbs(Sampler dir, int len){
      
        this.exp_rd = new SimpleExtendedDistribution1[len];
        for(int i=0; i<exp_rd.length; i++){
            exp_rd[i] = new SimpleExtendedDistribution1(dir);
        //    if(Constants.CHECK && dir.dist.length==2 && exp_rd[i].probs[0]<0.5){
        //        throw new RuntimeException("!!");
        //    }
        }
    }
   public double[] getPrior() {
		double[] d = this.exp_rd[0].probs();
		double[] res = new double[d.length];
		System.arraycopy(d, 0, res, 0, d.length);
	return d;
	}
  
   public MultiExpProbs(MultiExpProbs probs) {
       int len = probs.exp_rd.length;
       this.exp_rd = new SimpleExtendedDistribution1[len];
       for(int i=0; i<exp_rd.length; i++){
           Dirichlet dir = new Dirichlet(probs.getExp(i).probs, Double.POSITIVE_INFINITY);
           exp_rd[i] = new SimpleExtendedDistribution1(dir);
       }
       
  }
   public  MultiExpProbs(Sampler dir, int len, boolean exp, double[] hs){
		 // double[] hs = Constants.expModelIntHotSpot1();
		  
		   if(exp ){
			   Double[] probs = dir.dist;
			   double log = Math.log(probs[0]);
			   this.exp_rd = new SimpleExtendedDistribution1[len];
		       for(int i=0; i<exp_rd.length; i++){
		    	 
		    		   double exp1 = Math.exp(log*hs[i]);
		    		   Dirichlet dir1 = new Dirichlet(new double[] {exp1,1-exp1}, Double.POSITIVE_INFINITY);
		    		   
		    	  
		           exp_rd[i] = new SimpleExtendedDistribution1(dir1);
		       //    if(Constants.CHECK && dir.dist.length==2 && exp_rd[i].probs[0]<0.5){
		       //        throw new RuntimeException("!!");
		       //    }
		       }
		   }
		   else{
			//   double[] mod = Constants.modifyFrac0;
		       this.exp_rd = new SimpleExtendedDistribution1[len];
		       for(int i=0; i<exp_rd.length; i++){
		    	//   double[] dist = new double[mod.length+1];
			      //  		 System.arraycopy(mod, 0, dist, 1,mod.length);
		           exp_rd[i] = new SimpleExtendedDistribution1(dir);
		       //    if(Constants.CHECK && dir.dist.length==2 && exp_rd[i].probs[0]<0.5){
		       //        throw new RuntimeException("!!");
		       //    }
		       }
		   }
	   }
   public MultiExpProbs(ExpTransProb probs, int[] statesToGroup, double[] u) {
       exp_rd = new SimpleExtendedDistribution1[statesToGroup.length];
       for(int i=0; i<exp_rd.length; i++){
           exp_rd[i] = (SimpleExtendedDistribution1)probs.getExp(statesToGroup[i]).clone(u[i]);
       }
}
public ExpTransProb clone(int[] statesToGroup, double[] u) {
       return new MultiExpProbs(this, statesToGroup, u);
    
   }
   
public ExpTransProb clone(boolean swtch, int noStates) {
      return new MultiExpProbs(this);
  }
/* each row is for different state, within each row we have between then within */
    public double evaluateExpRd(double pseudoC, double[] d) {
        double sum=0;
        double[] tmp = new double[2];
       
        for(int i=1; i<this.exp_rd.length; i++){
        	double v = Math.exp(d[i]);
        	 tmp[0] = pseudoC*v;
        	 tmp[1] = pseudoC*(1 - v);
             sum+=this.exp_rd[i].evaluate(tmp);
            
        }
      //  if(exp_rd[2].probs[0]<0.5){
      // 	 Logger.global.info("h");
      //  }
        return sum;
     }
    
    public double evaluateExpRd(double[][]pseudoC) {
        double sum=0;
       
       
        for(int i=0; i<this.exp_rd.length; i++){
        	
             sum+=this.exp_rd[i].evaluate(pseudoC[i]);
        }
        return sum;
     }
    
    public double evaluateExpRd(double[][]pseudoC, double mod) {
        double sum=0;
       
       
        for(int i=0; i<this.exp_rd.length; i++){
        	
             sum+=this.exp_rd[i].evaluate(pseudoC[i], mod);
        }
        return sum;
     }
    public SimpleExtendedDistribution1 getExp(int groupFrom){
       return exp_rd[groupFrom];
    }
    public Collection getExpRdColl(){
        return Arrays.asList(exp_rd);
    }
public void initialiseExpRd(){
   for(int i=0; i<exp_rd.length; i++){
       this.exp_rd[i].initialise();
   }
   }
  /* public void transferExp(double[] pseudoExp){
       for(int i=0; i<exp_rd.length; i++){
       exp_rd[i].transfer(pseudoExp);
       }
   }*/
public void transferExp(double pseudoExp){
    for(int i=0; i<exp_rd.length; i++){
        exp_rd[i].transfer(pseudoExp);
    }
}
   public void validateExp(){
       for(int i=0; i<exp_rd.length; i++){
       exp_rd[i].validate();
       }
   }

    public void printExp(PrintWriter pw, double dist, String pref){
        pw.print(pref+" ");
        for(int i=0; i<exp_rd.length; i++){
          pw.print(i+" : "); 
        // if(exp_rd[i].probs.length==2) pw.print(AbstractTransitionProbs.transform(exp_rd[i].probs[0], dist));
        //else{
              exp_rd[i].printSimple(pw, "", "", 0.3);
       //  }
          pw.println("; ");
        }
    }
    public int noStates() {
       return this.exp_rd[0].probs.length;
    }
    public String info() {
       StringBuffer sb = new StringBuffer(this.getClass().toString());
       return sb.toString();
    }
    public double[] getProbs(int i){
    	return this.exp_rd[i].probs();
    }
	public void modify(int[] is) {
		double skew = Constants.skewTransitions(1);
		double skew2 = Math.pow(skew, 2);
		for(int i=1; i<is.length; i++){
			if(is[i]>=0){
				double[] probs = this.exp_rd[i].probs;
				if(probs.length>1){
					Arrays.fill(probs,skew2);
					probs[is[i]] = skew;
					Constants.normalise(probs);
				}
			}
		}
		// TODO Auto-generated method stub
		
	}
   
   
}
