package lc1.stats;

import java.util.logging.Logger;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.util.Constants;

public class IlluminaDistribution extends IlluminaRDistribution{

  Double b=null;
    //boolean isnan =false;
    public IlluminaDistribution(short data_index){
      super(data_index);
      this.data_index = data_index;
      // if(data_index==0) System.exit(0);
    }
    public void getB(double[] b, double mult) {
		if(b!=null)b[0]+=mult*this.b;
		b[1]+=mult;
	};
    public Boolean probeOnly() {
		// TODO Auto-generated method stub
		return false;
	}
    
   
    @Override
	protected Double getBEst(EmissionStateSpace emstsp) {
	//	if(this.r!=null && (r<-0.5 || r> 0.5)) return null;
    	return this.b;
	}
    public IlluminaDistribution(IlluminaDistribution ill) {
    	super(ill);
     
       this.b = ill.b;
      // this.isnan = ill.isnan;
    }
    public Double getIntensity(boolean r) {
    	if(r) return this.r;
    	else return this.b;
    }
    
    public String getUnderlyingData(){
       String r1 = r==null ? "n":""+Math.round(r*100.0)/100.0;
        String b1 = b==null ? "b": ""+Math.round(b*100.0)/100.0;
        return r1+"_"+b1;
    }
  /*  @Override
    public double score(int j, IlluminaProbB[] probB, int pos){
        return  probB[this.data_index].calcB(j, b,pos);
    
    }
    
    @Override
    public  double score(int j, IlluminaProbB probBState, int pos){
      return  probBState.calcB(j, b, pos);
    }
  

   
    public void addCount(IlluminaProbB[] probB, Integer index, double val, int pos){
        if(b!=null)probB[this.data_index].addBCount(index, val, b, pos);
    }*/
    
    @Override
    public void addBCount(int j, double val, int i) {
    	  if(val>Constants.countThresh()){
    	DistributionCollection.dc.getDistribution(this.data_index, j,i).addCount(b, val);
    	  }
    	
    }
    
    @Override
    public void addRBCount(EmissionStateSpace emstsp, int j, double val, int i) {
    	int noCop = emstsp.getCN(j);
    	int noB = emstsp.getBCount(j);
    	/*if(b>0.9 && noB == 2){
    		System.err.println("h!!!!h "+val);
    	}*/
    	
    	DistributionCollection.dc.addRBCount(data_index, noCop,noB, val, r, b, i);
    	//.getDistributionRB(this.data_index, noCop, noB,i).addCount(r,b, val);
    	
    }
    @Override
    /** j is index, i is position */
    public double scoreB(int j, int i) {
    	return 
    	
 //  	Math.pow(
    			DistributionCollection.dc.getDistribution(this.data_index, j,i).probability(this.b)
    	//		,Constants.exponentB(data_index))
    			;

    }
    
    @Override
    /** j is index, i is position */
    public double scoreBR(EmissionStateSpace emstsp, int j, int i) {
    	try{
    		int no = emstsp.getCN(j);
        	int noB = emstsp.getBCount(j);
        	/*if(no>0){
        		System.err.println('h');
        	}*/
    	return 
    	
    	Math.pow(
    	
    			DistributionCollection.dc.getDistributionRB(this.data_index, no,noB,i).probability(this.r, this.b)
    			,Constants.exponentB(data_index))
    			;
    	}catch(Exception exc){
    		Logger.global.info("problem at "+i+" "+this.data_index+" exiting");
    		exc.printStackTrace();
    		System.exit(0);
    	}
    	return 0;
    	//		,Constants.exponentB(data_index));

    }
    
    
    
    public double scoreBR(int no,int noB, int i, int mixComponent) {
    	return 
    	
 	Math.pow(
    			DistributionCollection.dc.getDistributionRB(this.data_index, no,noB,i).probability(this.r, this.b,mixComponent)
    		,Constants.exponentB(data_index))
    			;

    }
    
    public boolean equals(PseudoDistribution dist){
    	if(dist.getClass().equals(this.getClass())){
    		if(super.equals(dist) && this.b.equals(((IlluminaDistribution)dist).b())) return true;;
    	}
    	return false;
    }

   
    public Double b() {
        return b;
    }
    

    
public Double b(int i) {
	if(DataCollection.datC.dc instanceof MatchedDistributionCollection){
		double refc = ((MatchedDistributionCollection)DataCollection.datC.dc).refCount(i);
		return b/(b+refc);
	}else{
		return b;
	}
	
}
   
    /*@Override
    public double[] calcDistribution(
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
    }
    private double score(int j, IlluminaProbB probB, IlluminaProbR probR, EmissionStateSpace emStsp, int pos) {
    
    	if(r==null && b==null) return emStsp.invSize();
      else{
          double sc1 = r==null ? emStsp.invCNVSize():
              probR.calcR(Constants.backgroundCount(),emStsp.getCN(j), r,pos);
          double sc2 = b==null ? emStsp.bSpaceSize(j) :
                  probB.calcB(j, b,pos) ;
       return sc1*sc2;
      }
    }*/

    

    @Override
    public PseudoDistribution clone() {
       return new IlluminaDistribution(this);
    }

    @Override
    public PseudoDistribution clone(double swtch) {
     return clone();
    }

    @Override
    public PseudoDistribution swtchAlleles() {
    if(Constants.format[this.data_index].startsWith("alleleCount")){
    	this.b = r-b;
    }else{
       this.b = 1-b;
    }
        return this;
    }

    public EmissionStateSpace getEmissionStateSpace(){
		throw new RuntimeException("!!");
	}
    
    public void setBR(Double b2, Double r2){
        if(r2==null || Double.isNaN(r2) || b2==null || Double.isNaN(b2)  ){
            r = null;b = null;
        }
        else if(r!=null){// && Constants.cumulativeR(this.data_index)>1){
        	  
        	   this.r = r2+ this.r;	
        	   this.b = b2 + this.b;
        	  
        	   
        }
       /*else {
        		this.b = b2;
        		this.r = r2;
        		
        	}*/
    }
    
    public void setBR1(Double b2, Double r2){
        
        		this.b = b2;
        		this.r = r2;
        		
        	
    }
    
   
    public void setB(Double b2) {
      // System.err.println("set B "+b2);
      if(b2==null || Double.isNaN(b2) ){
          b = null;
      }
      else if(b!=null && Constants.cumulativeR(this.data_index)>1){
   	   this.b = b2 + this.b;
   }
      else {
          this.b = b2;
      }
        
    }
   
    
    @Override
    public String getCompressedDataString(EmissionStateSpace emStSp) {
       return String.format("%5.3g",this.r)+"\t"+String.format("%5.3g",this.b);
    }
    public String compressedStringHeader(EmissionStateSpace emStSp) {
        return "Log R\tB allele";
       }
  
    
public void compareTo(PseudoDistribution pseudoDistribution) {
    if(pseudoDistribution instanceof IlluminaDistribution){
     if(Math.abs(this.b - ((IlluminaDistribution)pseudoDistribution).b)>0.75){
         throw new RuntimeException("b values are very different ");
     }
     if(Math.abs(this.r - ((IlluminaDistribution)pseudoDistribution).r)>0.75){
         throw new RuntimeException("r values are very different ");
     }  
    }
}
public void transform() {
	// TODO Auto-generated method stub
	
}



}
