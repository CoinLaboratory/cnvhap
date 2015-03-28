package lc1.stats;

import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;


public class SoftenedIntegerDistribution extends PseudoDistribution {
    Integer i;
    Integer noCop;
     Integer noB;
     int len;
     double prob;
     
public SoftenedIntegerDistribution(int i, Integer noCop, Integer noB, int len, double prob){
    this.i = i;
    this.noB = noB;
    this.noCop = noCop;
    this.len = len;
    this.prob = prob;
}
public void getB(double[] b, double mult){
	throw new RuntimeException("!!");
}
public SoftenedIntegerDistribution(int res, EmissionStateSpace emStSp, double prob) {
	this.i= res;
	this.noB = emStSp.getBCount(res);
	this.noCop = emStSp.getCN(res);
	this.prob = prob;
	this.len = emStSp.size();
}
public boolean isMeasured() {
	return true;
	}
public  double[] calcDistribution(
        double[] distribution, EmissionStateSpace emStSp, int pos) {
   Arrays.fill(distribution, (1-prob)*1.0/(double)len);
   distribution[this.i] += prob;
   return distribution;
 }

    public double KLDistance(PseudoDistribution d2) {
       throw new RuntimeException("!!");
    }

    public void addCount(int obj_index, double value) {
       cnt+=value;
    }
    
  
    @Override
    public void setFixedIndex(int k) {
       
        this.i = k;
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
        throw new RuntimeException("!!");
    }
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
       throw new RuntimeException("!!");
    }
   public PseudoDistribution clone(double swtch) {
       return new SoftenedIntegerDistribution(i, this.noCop, this.noB, this.len, this.prob);
    }
   public PseudoDistribution clone(){
       return new SoftenedIntegerDistribution(i, this.noCop, this.noB, this.len, this.prob);
   }

    public double[] counts() {
        throw new RuntimeException("!!");
    }
    public String getUnderlyingData(EmissionStateSpace emStSp){
        return emStSp.get(i.intValue())+"_"+1.0;
    }
    public double logProb() {
        // TODO Auto-generated method stub
        return 0;
    }

    public int getMax() {
       return i;
    }

    public String getPrintString() {
        return i+"";
    }
int cnt =0;
    public void initialise() {
      cnt=0;
    }

    public void print(PrintWriter pw, boolean b, String printString,
            String string) {
        pw.print(printString+" "+i);

    }

    public void printSimple(PrintWriter pw, String name, String newLine, double thresh){
        pw.print(name+"->{");
      //  for(int i=0; i<this.probs.length; i++){
        //    if(probs[i] > thresh){
                pw.print(i+":"+100+",");
         //   }
      //  }
        pw.print("}"+newLine);
    }

    public double[] probs() {
        throw new RuntimeException("!!");
    }

    public double probs(int obj_i) {
        if(obj_i==i) return this.prob +(1-prob)*(1.0/(double)(len)) ;
        else return (1-prob)*(1.0/(double)(len));
    }

    public double sample() {
    	if(Constants.rand.nextDouble()<prob)
       return i;
    	else{
    		return Constants.nextInt(len);
    	}
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
       return i;
    }
    public PseudoDistribution swtchAlleles() {
       this.i = null;//switchTranslation[i];
        this.noB = noCop-noB;
        if(true) throw new RuntimeException("!!");
        return this;
    }
  
   
    @Override
    public double totalCount() {
     return this.cnt;
    }
    public void mix(EmissionStateSpace emStSp) {
        if(Constants.rand.nextBoolean()){
            this.i = emStSp.flip( emStSp.getHaploPairFromHaplo(i));
        }
    }
    public  String getCompressedDataString(EmissionStateSpace emStSp){
        
        return  emStSp.getHapl(i.intValue()).toString().replaceAll("[,\\s+\\[\\]]","");
         //throw new RuntimeException("!!");
     }
     
    public String toString(){
    	return this.i+"";
    }
	@Override
	public double scoreB(int j, int i) {
	    if(j==this.i) return this.prob +(1-prob)*(1.0/(double)(len)) ;
        else return (1-prob)*(1.0/(double)(len));
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emstsp, int j, int i) {
		int noCop = emstsp.getCN(j);
		int noB = emstsp.getBCount(j);
		if(true) throw new RuntimeException("!!");
		if(noCop == this.noCop && noB == this.noB) return this.prob +(1-prob)*(1.0/(double)(len)) ;
		else return (1-prob)*(1.0/(double)(len));
	}


	

}
