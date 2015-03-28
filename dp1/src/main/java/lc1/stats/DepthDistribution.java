package lc1.stats;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.util.Constants;

public class DepthDistribution extends IlluminaRDistribution{
// public Number depth=null;
 
//PoissonDistribution dist;
//TrainableNegBinomial tnb;
//SkewNormal sn;
final int index_indiv;
public  double average_coverage_per_chromosome;
 

    //boolean isnan =false;
    public DepthDistribution(short data_index, Number depth,  double average_coverage, int index_indiv){
    	super(data_index);
    	this.index_indiv = index_indiv;
//       this.data_index = data_index;
       this.r = depth.doubleValue();
       this.average_coverage_per_chromosome = average_coverage;
    if(this.average_coverage_per_chromosome < 0.0001){
    	throw new RuntimeException("!!");
    }
      // this.sn = new SkewNormal();
      // if(data_index==0) System.exit(0);
    }
    public boolean isMeasured() {
	return this.r!=null;
	}
    
    public String getUnderlyingData(EmissionStateSpace emStSp){
        return this.getUnderlyingData();
    }
    public boolean equals(PseudoDistribution dist){
    	if(dist.getClass().equals(this.getClass())){
    		if(super.equals(dist) && this.r.equals(((DepthDistribution)dist).r())) return true;
    	}
    	return false;
    }
  //}
	public Boolean probeOnly() {
		if(r==null) return null;
		// TODO Auto-generated method stub
		return true;
	}
    public String getPrintString(){
       return "";
    }
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
        throw new RuntimeException("!!");
     }
    public DepthDistribution(DepthDistribution ill) {
    	super(ill);
    	this.index_indiv = ill.index_indiv;
      // this.r= ill.depth;
       this.average_coverage_per_chromosome = ill.average_coverage_per_chromosome;
   //  this.data_index = data_index;
      // this.isnan = ill.isnan;
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
        throw new RuntimeException("!!");
    }

    /*public void applyCorrection(double median) {
       if(r!=null) this.setR(r+median);
        
    }*/
    /*public void applyStdErrCorrection(double stderrmult){
    	if(r!=null) this.r = r*stderrmult;
    }*/
    	
    
    public String getUnderlyingData(){
       String r1 = r==null ? "n":""+Math.round(r.doubleValue()*100.0)/100.0;
        return r1;
    }
   
    
    @Override
    public void addCount(int obj_index, double value) {
        throw new RuntimeException("!!");
        
    }

   
 
    

    public Number r() {
        double d =  r.doubleValue()/ (this.average_coverage_per_chromosome);
        return d;
//       
        //return d;
    }
   @Override
    public Number r1(){
    	 return r();//Math.max(-3,Math.log(r().doubleValue()/2.0)/log2);
    }
   double log2 = Math.log(2);
    
   
    
    
   /* @Override
    public double[] calcDistribution(IlluminaProbB probB, double[] distribution, EmissionStateSpace emStSp, int nocop,int pos){
        double sum=0;
        Arrays.fill(distribution, 0.0);
        for(int j=0; j<distribution.length; j++){
                   if(emStSp.getCN(j)!=nocop){
                       distribution[j] = 0;
                   }
                   else{
                       distribution[j] = this.score(j, probB,pos)*emStSp.getWeight(j);  //do we include the weight term???;
                   }
                   sum+=distribution[j];
                }
        for(int i=0; i<distribution.length; i++){
            distribution[i] = distribution[i]/sum;
        }
        return distribution;
    }*/
    
   /* @Override
    public double[] calcDistribution(IlluminaProbR probR, IlluminaProbB probB,
            double[] distribution, EmissionStateSpace emStSp, int pos) {
        double sum=0;
        Arrays.fill(distribution, 0.0);
        for(int j=0; j<distribution.length; j++){
               
                  
                   distribution[j] = this.score(j, probB,probR, emStSp,pos)*emStSp.getWeight(j);  //do we include the weight term???;
                   sum+=distribution[j];
                }
        for(int i=0; i<distribution.length; i++){
            distribution[i] = distribution[i]/sum;
        }
        return distribution;
    }*/
   /* private double score(int j, IlluminaProbB probB, IlluminaProbR probR, EmissionStateSpace emStsp, int pos) {
    
    
          double sc1 = r==null ? emStsp.invCNVSize():
              probR.calcR(Constants.backgroundCount(),emStsp.getCN(j), r, pos);
          double sc2;
          if(probB.hasB(j)) {
      		return sc2 = 0.0;
      	}
      	else sc2 = 1.0;
         
       return sc1*sc2;
     
    }*/

    

    @Override
    public PseudoDistribution clone() {
       return new DepthDistribution(this);
    }

    @Override
    public PseudoDistribution clone(double swtch) {
     return clone();
    }

    @Override
    public double[] counts() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public double logProb() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Integer fixedInteger() {
        // TODO Auto-generated method stub
        return null;
    }

   

    @Override
    public int getMax() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public void initialise() {
        // TODO Auto-generated method stub
        
    }

    @Override
    public double[] probs() {
    //	if(true) throw new RuntimeException("!!");
     return null;//
     //this.calcDistribution(distribution, emStSp, pos);
    }
    public  double[] calcDistribution(
            double[] distribution, EmissionStateSpace emStSp, int pos) {
    	if(true) throw new RuntimeException("!!");
       return null;
     }

    @Override
    public double probs(int obj_i) {
        if(true) throw new RuntimeException("!!");// TODO Auto-generated method stub
        return 0;
    }

    public double sample() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public void setCounts(int i1, double cnt) {
        // TODO Auto-generated method stub
        
    }

    @Override
    public void setProb(double[] prob) {
        // TODO Auto-generated method stub
        
    }

    @Override
    public void setProbs(int to, double d) {
        // TODO Auto-generated method stub
        
    }

    @Override
    public double sum() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public PseudoDistribution  swtchAlleles() {
     return this;
        
    }

    @Override
    public void transfer(double pseudoC1) {
        // TODO Auto-generated method stub
        
    }

    @Override
    public void validate() {
        // TODO Auto-generated method stub
        
    }

   
   /* public void setR(Double r2){
        if(r2==null || Double.isNaN(r2) ){
            r = null;
        }
        else{
            this.r = r2;
        }
    }*/
    @Override
    public void setFixedIndex(int k) {
       throw new RuntimeException("!!");
        
    }
    @Override
    public double totalCount() {
        throw new RuntimeException("!!");
    }
    @Override
    public String getCompressedDataString(EmissionStateSpace emStSp) {
       return this.r+"";
    }
    public String compressedStringHeader(EmissionStateSpace emStSp) {
        return "Log R";
       }
  
    
public void compareTo(PseudoDistribution pseudoDistribution) {
    if(pseudoDistribution instanceof DepthDistribution){
    
     if(Math.abs(this.r.doubleValue() - ((DepthDistribution)pseudoDistribution).r.doubleValue())>0.75){
         throw new RuntimeException("r values are very different ");
     }  
    }
}
public void dropROutlier(double[] ds) {
if(r==null || this.r.doubleValue() < ds[0] || this.r.doubleValue()>ds[1]) r = null;
	
}
@Override
public Double getIntensity(boolean r) {
	if(r) return this.r.doubleValue();
	else return Double.NaN;
}
@Override
public void addBCount(int j, double val, int i) {
	
	
}



@Override
public void addRBCount(EmissionStateSpace emstsp, int j, double val, int i) {
	if(r!=null){
		int noB = emstsp.getBCount(j);
		int noCop = emstsp.getCN(j);
	if(noB==0){
		int i1 = this.index_indiv;
		DistributionCollection.dc.addRBCount(data_index, noCop,noB, val, r/this.average_coverage_per_chromosome, 0.0, i1);
	   //   getDistributionRB(this.data_index, noCop, noB,i).addCount(r,0, val);
	}
	}
}


@Override
/** j is index, i is position */
public double scoreB(int j, int i) {
	if(true) throw new RuntimeException("!!");
	//if(Emiss.getSpaceForNoCopies(Constants.backgroundCount).getBCount(j)>0) return 0;
	return 1.0;

}

@Override
/** j is index, i is position
 * we use the individual index rather than position here
 *  */
public double scoreBR(EmissionStateSpace emstsp, int j,   int i) {
	if(r==null) return 1.0;
	int noB = emstsp.getBCount(j);
	int noCop = emstsp.getCN(j);
	if(noB>0) return 0;
	else {
		int i1 = this.index_indiv;
	ProbabilityDistribution dist = 	((lc1.stats.OrthogonalProbabilityDistribution)DistributionCollection.dc.getDistributionRB(this.data_index, noCop,noB, i1)).distx;
	//dist.setCoverage(noCop);//*);	
	double sc =  Math.pow(
    		dist.probability(this.r.doubleValue()/this.average_coverage_per_chromosome)
    			,Constants.exponentB(data_index))
    			;
	if(Double.isNaN(sc)){
		throw new RuntimeException("!!");
	}
	return sc;
		//dist.setCoverage(noCop*this.average_coverage_per_chromosome);
		
	}

}



/** j is index, i is position */
@Override
public double scoreBR(EmissionStateSpace emstsp, int j,  int i, int mixComponent) {
	return scoreBR(emstsp,j,i);

}
@Override
public void addRCount(int bg, int no_cop, double val, int i) {
	//super.ad
	if(r!=null && val >Constants.countThresh()){
		int i1 = this.index_indiv;
	DistributionCollection.dc.getDistribution(this.data_index, bg,no_cop,i1).addCount(r/this.average_coverage_per_chromosome, val);
	}
	
}

@Override
public double scoreR(int bg, int no_cop, int i) {
	if(true) throw new RuntimeException("!!");
	return 0;
	
}


}
