package lc1.stats;

import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;


public class UniformDist extends PseudoDistribution {
   int len;
public UniformDist(int len){
   this.len = len;
}

public boolean isMeasured() {
	return false;
	}
public  double[] calcDistribution(
        double[] distribution, EmissionStateSpace emStSp, int pos) {
   Arrays.fill(distribution, 1.0/len);
   return distribution;
 }

   

    public void addCount(int obj_index, double value) {
      
    }
    
  
    @Override
    public void setFixedIndex(int k) {
       throw new RuntimeException("!!");
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
        throw new RuntimeException("!!");
    }
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
       throw new RuntimeException("!!");
    }
   public PseudoDistribution clone(double swtch) {
       return new UniformDist(this.len);
    }
   public PseudoDistribution clone(){
       return new UniformDist(this.len);
   }

    public double[] counts() {
        throw new RuntimeException("!!");
    }
    public String getUnderlyingData(EmissionStateSpace emStSp){
        return "U";
    }
    public double logProb() {
        // TODO Auto-generated method stub
        return 0;
    }

    public int getMax() {
       return 0;
    }

    public String getPrintString() {
        return "U";
    }

    public void initialise() {
      
    }

    public void print(PrintWriter pw, boolean b, String printString,
            String string) {
        pw.print(printString+" U");

    }

    public void printSimple(PrintWriter pw, String name, String newLine, double thresh){
        pw.print(name+"->{");
      //  for(int i=0; i<this.probs.length; i++){
        //    if(probs[i] > thresh){
                pw.print("U,");
         //   }
      //  }
        pw.print("}"+newLine);
    }

    public double[] probs() {
        throw new RuntimeException("!!");
    }

    public double probs(int obj_i) {
       return 1.0/(double)len;
    }

    public double sample() {
       return Constants.nextInt(len);
    }

    public void setCounts(int i1, double cnt) {
        throw new RuntimeException("!!");

    }

    public void setProb(double[] prob) {
        throw new RuntimeException("!!");

    }

    public void setProbs(int to, double d) {
        throw new RuntimeException("!!");

    }

    public double sum() {
        return 1.0;
    }

    public void transfer(double pseudoC1) {
//      throw new RuntimeException("!!");

    }

    public void validate() {
        // TODO Auto-generated method stub

    }
    public Integer fixedInteger() {
       return null;
    }
    public PseudoDistribution swtchAlleles() {
      return this;
        
    }
  
   
    @Override
    public double totalCount() {
     return 0;
    }
    public void mix(EmissionStateSpace emStSp) {
        throw new RuntimeException ("!!");
    }
    public  String getCompressedDataString(EmissionStateSpace emStSp){
        
        return  "";
         //throw new RuntimeException("!!");
     }
     
    public String toString(){
    	return "U";
    }
	@Override
	public double scoreB(int j, int i) {
		return 1.0/len;
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emstsp, int j, int i) {
		return 1.0/(double)len;
	}


	

}
