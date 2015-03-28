package lc1.stats;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.util.Constants;


public class WeightedIlluminaRDistribution extends PseudoDistribution{
 protected Double r=null;
 protected Double numReads=null;
 protected Double averageNumReads=null;
 protected PoissonDistribution currPoisson=null;
    //boolean isnan =false;
    public WeightedIlluminaRDistribution(short data_index){
       this.data_index = data_index;
      // if(data_index==0) System.exit(0);
    }
    
    
    public WeightedIlluminaRDistribution(short data_index, double averageNumReads){
        this.data_index = data_index;
        this.averageNumReads = averageNumReads;
        currPoisson = new PoissonDistribution(averageNumReads);
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
    		if(super.equals(dist) && this.r.equals(((WeightedIlluminaRDistribution)dist).r())) return true;
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
    public WeightedIlluminaRDistribution(WeightedIlluminaRDistribution ill) {
       this.r = ill.r;
     this.data_index = ill.data_index;
      // this.isnan = ill.isnan;
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
        throw new RuntimeException("!!");
    }

    public void applyCorrection(double median) {
       if(r!=null) this.setR(r+(-1)*median);
        
    }
    public void applyStdErrCorrection(double stderrmult){
    	if(r!=null) {
    		this.r = r*stderrmult;
    	}
    }
    	
    
    public String getUnderlyingData(){
       String r1 = r==null ? "n":""+Math.round(r*100.0)/100.0;
        return r1;
    }
   /* @Override
    public double score(int j, IlluminaProbB[] probB, int pos){
    	if(probB[data_index].hasB(j)) {
    		return 0.0;
    	}
    	else return 1.0;
    }
    
    @Override
    public  double score(int j, IlluminaProbB probBState, int pos){
    	if(probBState.hasB(j)) {
    		return 0.0;
    	}
    	else return 1.0;
//        double sc =  probBState.calcB(j, 0);
 //       return sc;
    }*/
    
    @Override
    public void addCount(int obj_index, double value) {
        throw new RuntimeException("!!");
        
    }

   
 /*   public void addCount(IlluminaProbB[] probB, Integer index, double val,int pos){
      
    }*/
    
    

    public Number r() {
        return r;
    }
    
    public Number numReads() {
        return this.numReads;
    }
   
    
    
    
   /* @Overrider
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
       return new WeightedIlluminaRDistribution(this);
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
	protected Double getBEst(EmissionStateSpace emstsp) {
		return null;
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
    public PseudoDistribution swtchAlleles() {
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

   
    public void setR(Double r2){
        if(r2==null || Double.isNaN(r2) ){
            r = null;
        }
        else{
            this.r = r2;
        }
    }
    
    public void setNumReads(Double numReads){
        if(numReads==null || Double.isNaN(numReads) ){
            this.numReads = null;
        }
        else{
            this.numReads = numReads;
        }
    }
    
    public void setAverageNumReads(Double averageNumReads){
        if(averageNumReads==null || Double.isNaN(averageNumReads) ){
            this.averageNumReads = null;
        }
        else{
            this.averageNumReads = averageNumReads;
        }
    }
    
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
       return String.format("%5.3g",this.r);
    }
    public String compressedStringHeader(EmissionStateSpace emStSp) {
        return "Log R";
       }
  
    
public void compareTo(PseudoDistribution pseudoDistribution) {
    if(pseudoDistribution instanceof WeightedIlluminaRDistribution){
    
     if(Math.abs(this.r - ((WeightedIlluminaRDistribution)pseudoDistribution).r)>0.75){
         throw new RuntimeException("r values are very different ");
     }  
    }
}
public void dropROutlier(double[] ds) {
if(r==null || this.r < ds[0] || this.r>ds[1]){
	r = null;
}
	
}
@Override
public Double getIntensity(boolean r) {
	if(r) return this.r;
	else return Double.NaN;
}
@Override
public void addBCount(int j, double val, int i) {
	
	
}



@Override
public void addRBCount(EmissionStateSpace emstsp, int j, double val, int i) {
	if(r!=null){
	//if(noB==0){
		int noB = emstsp.getBCount(j);
		int noCop = emstsp.getCN(j);
		DistributionCollection.dc.addRBCount(data_index, noCop, noB, val, r, 0.0, i);
	   //   getDistributionRB(this.data_index, noCop, noB,i).addCount(r,0, val);
	//}
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
/** j is index, i is position */
public double scoreBR(EmissionStateSpace emstsp, int j,  int i) {
	if(r==null) return 1.0;
	int noB = emstsp.getBCount(j);
	int noCop = emstsp.getCN(j);
	if(noB>0) return 0;
	else {

	try {
		if(true){
			throw new RuntimeException("need to fix this");
		}
		return 1.0;
//		Math.pow(
//					DistributionCollection.dc.getDistributionRB(this.data_index, noCop,noB, i).probability(this.r, 0)
//					,Constants.exponentB(data_index)*this.currPoisson.cdf((Double)this.numReads()))
//					;
	} catch (Exception e) {
		// TODO Auto-generated catch block
		System.err.println(e.getMessage());
		return Math.pow(
				DistributionCollection.dc.getDistributionRB(this.data_index, noCop,noB, i).probability(this.r, 0)
				,Constants.exponentB(data_index));
	}
		/*return scoreR(Constants.backgroundCount(), EmissionStateSpace.emStSp.noCopies(),i);*/
	}

}



/** j is index, i is position */
@Override
public double scoreBR(EmissionStateSpace emstsp, int j,  int i, int mixComponent) {
	if(r==null) return 1.0;
	int noB = emstsp.getBCount(j);
	if(noB>0) return 0;
	
	else {
	int noCop = emstsp.getCN(j);
		//return 
	try {
		if(true){
			throw new RuntimeException("need to fix this");
		}
		return 1.0
		//Math.pow(
		//			DistributionCollection.dc.getDistributionRB(this.data_index, noCop,noB, i).probability(this.r, 0)
		//			,Constants.exponentB(data_index)*this.currPoisson.cdf((Double)this.numReads()))
					;
	} catch (Exception e) {
		// TODO Auto-generated catch block
		System.err.println(e.getMessage());
		return Math.pow(
				DistributionCollection.dc.getDistributionRB(this.data_index, noCop,noB, i).probability(this.r, 0)
				,Constants.exponentB(data_index));
	}
		/*return scoreR(Constants.backgroundCount(), EmissionStateSpace.emStSp.noCopies(),i);*/
	}

}
@Override
public void addRCount(int bg, int no_cop, double val, int i) {
	if(r!=null && val >Constants.countThresh()){
	DistributionCollection.dc.getDistribution(this.data_index, bg,no_cop,i).addCount(r, val);
	}
	
}

@Override
public double scoreR(int bg, int no_cop, int i) {
	if(true) throw new RuntimeException("!!");
	if(r==null) return 1.0;
	return 
	Math.pow(
	DistributionCollection.dc.getDistribution(this.data_index, bg,no_cop,i).probability(this.r),
	Constants.exponentR(data_index));
}
public Number depth() {
	// TODO Auto-generated method stub
	return this.r;
}
public void setDepth(double r){
	this.r = r;
}
public Number r1() {
	return r;
}

}
