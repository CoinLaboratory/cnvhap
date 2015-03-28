package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Collection;

import lc1.stats.SimpleExtendedDistribution;



public interface ExpTransProb extends Serializable {
    
    public  SimpleExtendedDistribution getExp(int groupFrom);
    public void initialiseExpRd();
   // public void transferExp(double[] pseudoExp);
   public void transferExp(double pseudoExp);
    public void validateExp();
    public double evaluateExpRd(double pseudoC, double[] d) ;
    public double evaluateExpRd(double[][] pseudoC) ;
    
    public double evaluateExpRd(double[][] pseudoC, double mod) ;
    public Collection getExpRdColl();
    public void printExp(PrintWriter pw, double dist, String pref);
    public ExpTransProb clone(boolean swtch, int noStates);
    public int noStates();
    public ExpTransProb clone(int[] statesToGroup, double[] u);
    public String info();
   // public double evaluateExpRd(double[]pseudoC);
	public double[] getPrior();
}
