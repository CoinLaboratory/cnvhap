package lc1.dp.transition;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;

import lc1.stats.Dirichlet;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

public class SimpleExpTransProb implements ExpTransProb {
    SimpleExtendedDistribution exp_rd1;
    public String info() {
        StringBuffer sb = new StringBuffer(this.getClass().toString());
        return sb.toString();
     }
    public ExpTransProb clone(boolean swtch, int noStates){
        if(swtch){
          return new MultiExpProbs(new Dirichlet(exp_rd1.probs, Constants.switchU()), noStates);
        }
        else return new SimpleExpTransProb(new Dirichlet(exp_rd1.probs, Double.POSITIVE_INFINITY), noStates);
    }
    public SimpleExpTransProb(Sampler samplerFirst, int len) {
        exp_rd1 = new SimpleExtendedDistribution1(samplerFirst);
      //  if(Constants.CHECK) this.validateExp();
    }
    public SimpleExpTransProb(Sampler samplerFirst, int len, boolean exp, double[] hs) {
        exp_rd1 = new SimpleExtendedDistribution1(samplerFirst);
      //  if(Constants.CHECK) this.validateExp();
    }
    
    
    public ExpTransProb clone(int[] statesToGroup, double[] u) {
        return new MultiExpProbs(this, statesToGroup, u);
     
    }
    
    
    public void printExp(PrintWriter pw, double dist, String pref){
        pw.print(pref+" ");
        //  if(exp_rd1.probs.length==2) pw.print(AbstractTransitionProbs.transform(exp_rd1.probs[0], dist));
        //  else{
              exp_rd1.printSimple(pw, "", "", 0.3);
       //   }
          pw.print("; ");
    }
    public SimpleExpTransProb(ExpTransProb exp) {
        this.exp_rd1=  new SimpleExtendedDistribution(exp_rd1);
    }
    public void initialiseExpRd(){
        this.exp_rd1.initialise();
    }
       /*public void transferExp(double[] pseudoExp){
           exp_rd1.transfer(pseudoExp);
       }*/
       public void transferExp(double pseudoExp){
           exp_rd1.transfer(pseudoExp);
       }
       public void validateExp(){
           exp_rd1.validate();
       }

       public Collection getExpRdColl(){
           return Arrays.asList(new SimpleExtendedDistribution[] {exp_rd1});
       }


   
    public double evaluateExpRd(double pseudoC, double[] d) {
    	 double[] tmp = new double[2];
    	 double v = Math.exp(d[0]);
    	 tmp[0] = pseudoC*v;
    	 tmp[1] = (1 - v)*pseudoC;
        return this.exp_rd1.evaluate(tmp);
     }
    public double evaluateExpRd(double[][]pseudoC) {
      
        return this.exp_rd1.evaluate(pseudoC[0]);
     }
    
    public double evaluateExpRd(double[][]pseudoC,
    	double mod) {
        
        return this.exp_rd1.evaluate(pseudoC[0], mod);
     }
    public SimpleExtendedDistribution getExp(int groupFrom) {
       return exp_rd1;
    }
   
    public int noStates() {
        return this.exp_rd1.probs.length;
     }
	public double[] getPrior() {
		double[] d = this.exp_rd1.probs();
		double[] res = new double[d.length];
		System.arraycopy(d, 0, res, 0, d.length);
	return d;
	}
}
