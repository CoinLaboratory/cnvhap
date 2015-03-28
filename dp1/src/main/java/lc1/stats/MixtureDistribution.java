package lc1.stats;

import java.io.PrintWriter;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.PairEmissionState;
import lc1.util.Constants;


public class MixtureDistribution extends PseudoDistribution {
   public PseudoDistribution[] dist;
   double[ ] mix;
   public  short getDataIndex(){
       return dist[0].getDataIndex();
   }
   
   
public MixtureDistribution(double[] mix, PseudoDistribution[] dist){
	if(dist[0]==null) throw new RuntimeException("!!");
   this.mix = mix;
   this.dist = dist;
   //this.distribution = new double[mix.length];
}

public EmissionStateSpace getEmissionStateSpace(){
	return dist[0].getEmissionStateSpace();
}
public void getB(double[] b, double mult) {
	
	for(int i=0; i<mix.length; i++){
		dist[i].getB(b,mix[i]*mult);
	}
	
};
public MixtureDistribution(double[] mix2, PseudoDistribution[] dist2,
		double swtch) {
	this.mix = mix2;
	this.dist = new PseudoDistribution[dist2.length];
	for(int  i=0; i<dist.length; i++){
		dist[i] = dist2[i].clone(swtch);
	}
	// this.distribution = new double[mix.length];
	// TODO Auto-generated constructor stub
}

public boolean isMeasured() {
	return true;
	}
public  double[] calcDistribution(
        double[] distribution, EmissionStateSpace emStSp, int pos) {
  
  throw new RuntimeException("!!");
 }

    public double KLDistance(PseudoDistribution d2) {
       throw new RuntimeException("!!");
    }

    public void addCount(int obj_index, double value) {
    	if(true) throw new RuntimeException("!!");
     //  cnt+=value;
    }
    
  
    @Override
    public void setFixedIndex(int k) {
       
        for(int i=0; i<dist.length; i++){
        	dist[i].setFixedIndex(k);
        }
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
        throw new RuntimeException("!!");
    }
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
       throw new RuntimeException("!!");
    }
   public PseudoDistribution clone(double swtch) {
       return new MixtureDistribution(this.mix, this.dist, swtch);
    }
   public PseudoDistribution clone(){
       return new MixtureDistribution(this.mix, this.dist);
   }

    public double[] counts() {
        throw new RuntimeException("!!");
    }
    public String getUnderlyingData(EmissionStateSpace emStSp){
       StringBuffer sb = new StringBuffer(this.dist[0].getUnderlyingData(emStSp));
       for(int i=1; i<dist.length; i++){
    	   sb.append(dist[i].getUnderlyingData(emStSp));
       }
       return sb.toString();
    }
    public double logProb() {
        // TODO Auto-generated method stub
        return 0;
    }

   

    public String getPrintString() {
    	StringBuffer sb = new StringBuffer(this.dist[0].getPrintString());
        for(int i=1; i<dist.length; i++){
     	   sb.append(dist[i].getPrintString());
        }
        return sb.toString();
    }
//int cnt =0;
    public void initialise() {
    //  cnt=0;
    }

    public void print(PrintWriter pw, boolean b, String printString,
            String string) {
        pw.print(printString+" "+this.getPrintString());

    }

    public void printSimple(PrintWriter pw, String name, String newLine, double thresh){
        pw.print(name+"->{");
        for(int i=0; i<dist.length; i++){
        dist[i].printSimple(pw, name, newLine, thresh);
        }
      //  for(int i=0; i<this.probs.length; i++){
        //    if(probs[i] > thresh){
               
         //   }
      //  }
        pw.print("}"+newLine);
    }

    public double[] probs() {
        throw new RuntimeException("!!");
    }
    
protected Double getBEst(EmissionStateSpace emstsp) {
		return this.dist[0].getBEst(emstsp);
	}

    public double probs(int obj_i) {
    	double p = 0;
    	for(int i=0; i<dist.length; i++){
    		p+=dist[i].probs(obj_i)*mix[i];
    	}
       return p;
    }

    public double sample() {
    	throw new RuntimeException("!!");
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
    
     for(int i=0; i<dist.length; i++){
    	 dist[i] = dist[i].swtchAlleles();
     }
        return this;
    }
  
   
    @Override
    public double totalCount() {
    	throw new RuntimeException("!!");
    }
    public void mix(EmissionStateSpace emStSp) {
        throw new RuntimeException("!!");
    }
    public  String getCompressedDataString(EmissionStateSpace emStSp){
        
        return  "";//emStSp.getHapl(i.intValue()).toString().replaceAll("[,\\s+\\[\\]]","");
         //throw new RuntimeException("!!");
     }
     
    public String toString(){
    	return "mix";
    }
	@Override
	public double scoreB(int j, int i) {
	  throw new RuntimeException("!!");
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emstsp,int j, int i) {
		double p = 0;
		//double[] d = new double[emstsp.genoListSize()];
    	for(int ik=0; ik<dist.length; ik++){
    		//d[ik] = dist[ik].scoreBR(emstsp,j, i)*mix[ik];
    		p+=dist[ik].scoreBR(emstsp,j, i)*mix[ik];
    	}
       return p;
	}
	//final double[] distribution;
	
	public double[]  getDist(EmissionStateSpace emstsp, int j, int i, double[] distribution){
	
		for(int ik=0; ik<dist.length; ik++){
    		distribution[ik]=dist[ik].scoreBR(emstsp, j, i)*mix[ik];
    	}
		return distribution;
	}
	public void addRBCount(EmissionStateSpace emstsp, int j, double val, int i) {
		//if(Constants.CHECK && Double.isNaN(d)){
			//throw new RuntimeException("!!");
	//	}
		{
			double[] distribution = PairEmissionState.pool.getObj(dist.length);
		 this.getDist(emstsp, j, i,distribution);
			double sum = Constants.sum(distribution);
			for(int ik=0; ik<distribution.length; ik++){
				dist[ik].addRBCount(emstsp,j, val*(distribution[ik]/sum), i);
			}
			PairEmissionState.pool.returnObj(distribution);
		}
		
	}

	@Override
	public int getMax() {
		// TODO Auto-generated method stub
		return 0;
	}


	

}
