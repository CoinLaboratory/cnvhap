package lc1.stats;

import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;


public class IntegerDistribution extends PseudoDistribution {
    Integer i;
   public Integer noCop;
   public  Integer noB;
     public EmissionStateSpace emstsp;
public IntegerDistribution(int i, Integer noCop, Integer noB){
    this.i = i;
    this.noB = noB;
    this.noCop = noCop;
}

public void setEmStSp(EmissionStateSpace stsp) {
	this.emstsp = stsp;
	}
public int[] getOrder(){
   	int[] res = new int[emstsp.defaultList.size()];
   
   	for(int i=0; i<res.length; i++){
   		res[i] = i;
   	}
	res[0] = i;
   	res[i] = 0;
	return res;
   }

public void getB(double[] b, double mult) {
	b[0]+=this.noB.doubleValue()*mult;
	b[1]+=this.noCop.doubleValue()*mult;
};
public IntegerDistribution(int res, EmissionStateSpace emStSp) {
	this.i= res;
	this.noB = emStSp.getBCount(res);
	this.noCop = emStSp.getCN(res);
	this.emstsp = emStSp;
}
public EmissionStateSpace getEmissionStateSpace(){
	return this.emstsp;
}
public boolean isMeasured() {
	return true;
	}
public  double[] calcDistribution(
        double[] distribution, EmissionStateSpace emStSp, int pos) {
   Arrays.fill(distribution, 0);
   distribution[this.i] = 1.0;
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
       return new IntegerDistribution(i, this.noCop, this.noB);
    }
   public PseudoDistribution clone(){
       IntegerDistribution res =  new IntegerDistribution(i, this.noCop, this.noB);
       res.emstsp = emstsp;
       return res;
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
public double cnt =0;
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
        if(obj_i==i) return 1.0;
        else return 0.0;
    }

    public double sample() {
       return i;
    }

    public void setCounts(int i1, double cnt) {
        throw new RuntimeException("!!");

    }

    public void setProb(double[] prob) {
        throw new RuntimeException("!!");

    }

    public void setProbs(int to, double d) {
    	if(true)  throw new RuntimeException("!!");
    	/*if(d==1){
    		this.i = to;
    		this.noCop = emstsp.getCN(i);
    		this.noB = emstsp.getBCount(i);
    	}else if(d>0){
        throw new RuntimeException("!!");
    	}*/
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
    	PseudoDistribution after = ((CompoundEmissionStateSpace)this.emstsp).getIntDist(
    			emstsp.getSwitchTranslation()[i]);
    	return after;
    
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
		if(j==this.i) return 1.0;
		else return 0.0;
	}
	
	@Override
	public final double scoreBR(EmissionStateSpace emstsp, int j, int i) {
	//	Comparable compa = emstsp.get(j);
		int noCop = emstsp.getCN(j);
		//double exp = ;
		//if(Constants.expon)
		if(noCop==this.noCop){
	//	if(Constants.CHECK && (this.noB==null || this.noCop==null)) throw new RuntimeException("");
		if(emstsp.get(j).equals(this.emstsp.get(this.i.intValue()))){
			return 1.0;
			
			
		}
		else{
			return 0;
//			return Constants.exponentB(this.data_index) <=0  ? 1.0 : 0.0;
		}
		}
		else return 0;//Constants.exponentB(this.data_index) <=0  ? 1.0 : 0.0;
	}


	

}
