package lc1.dp.states;

import java.util.List;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.WeightedIlluminaRDistribution;
import lc1.util.Constants;
/** warning sum[i] is only set by getEmiss()! */
public class IlluminaNoBg extends  HaplotypeEmissionState{


  @Override
  public int getParamIndex() {
     return 0;
  }
//  IlluminaProbR probR;
 // IlluminaProbB probB;
 // int[] cn;
   public IlluminaNoBg(final String name,
           EmissionStateSpace emStSp, 
       //    IlluminaProbR r, IlluminaProbB b ,
         List<String> snpid, int cumul, short index
           ){
       super(name, (int) Math.floor((float)snpid.size()/(float) cumul), emStSp, index);
   //    this.probB = b;
     //  this.probR = r;
       Boolean probeOnly = Constants.probeOnly(index);
       boolean po = probeOnly!=null && probeOnly;
   //    int cumul = Constants.cumulativeR(index);
       String format = Constants.format[index];
      for(int i=0; i<this.noSnps; i++){
    	  String id = snpid.get(i*cumul);
    	  if(format.equals("matcheddepth")){
    		  emissions[i] = new IlluminaDistribution(index){
    			  
    			  public Double b(int i) {
    				  MatchedDistributionCollection mdc = ((MatchedDistributionCollection)DataCollection.datC.dc);
    				  double p = b()/r().doubleValue();
    				  double p1 = mdc.refCount(i)/mdc.totnormal();
    				  double mult = p/p1;
    				  
    				 return mdc.inverseMult(mult, i,name);//(mult - (1-cellularity))* (ratio/cellularity);
    						 //Constants.ratioAsLevels()!=null && Constants.ratioAsLevels() ? p/p1 : (p/p1)*ratio; //CHANGED 01/12/2014
    			  }
    			  @Override
    			  public Double r1(){
    				  return this.b();
    			  }
    		  };
    	  }
    	  else if(po || id.startsWith("A_") || id.startsWith("cnv")){
        	 emissions[i] = new IlluminaRDistribution(index);
         }
         else if (id.startsWith("ins")){
        	 emissions[i] = 
    			 new WeightedIlluminaRDistribution(index);
         }
         else {
        	 emissions[i] = 
        			 new IlluminaDistribution(index);
         }
          
      }
     
   }
   
   public IlluminaNoBg(String name, int noSnps, EmissionStateSpace emStSp, short data_index){
       super(name, noSnps, emStSp, data_index);
       for(int i=0; i<this.noSnps; i++){
     	
         	 emissions[i] = 
         		
         		 new IlluminaDistribution(data_index);
         
           
       }
   }
public IlluminaNoBg(IlluminaNoBg st) {
    super(st);
    for(int i=0; i<this.noSnps; i++){
        emissions[i] = st.emissions[i].clone();
        
    }
  //  this.outOfRangePos = st.outOfRangePos;
}






	

public void setB(int i, Double b2) {
 // if(Double.isNaN(b2)) ((IlluminaDistribution)this.emissions[i]).setNaN(true);
   // System.err.println("setting B "+this.getName()+i+" "+b2);
	if(emissions[i] instanceof IlluminaDistribution)
    ((IlluminaDistribution) this.emissions[i]).setB(b2);
   
 }



public void setBR(int i,  Double b2, Double r2) {
	 // if(Double.isNaN(b2)) ((IlluminaDistribution)this.emissions[i]).setNaN(true);
	   // System.err.println("setting B "+this.getName()+i+" "+b2);
		if(emissions[i] instanceof IlluminaDistribution)
	    ((IlluminaDistribution) this.emissions[i]).setBR(b2,r2);
	   
	 }

public void setNumReads(int i, Double b2) {
	 // if(Double.isNaN(b2)) ((IlluminaDistribution)this.emissions[i]).setNaN(true);
	   // System.err.println("setting B "+this.getName()+i+" "+b2);
		if(emissions[i] instanceof WeightedIlluminaRDistribution)
	    ((WeightedIlluminaRDistribution) this.emissions[i]).setNumReads(b2);
	   
	 }


 public void setR(int i, double r2) {
	PseudoDistribution dist = null;
	 if (this.emissions[i] instanceof IlluminaRDistribution)
		((IlluminaRDistribution) this.emissions[i]).setR(r2);
	else if (this.emissions[i] instanceof WeightedIlluminaRDistribution)
		((WeightedIlluminaRDistribution) this.emissions[i]).setR(r2);
	//if(dist.r()!=null)
   //  if(Double.isNaN(r2)) ((IlluminaDistribution)this.emissions[i]).setNaN(true);
//    dist.setR(r2);
   
  }

/*public void remove(int j) {
  r.remove(j);
  b.remove(j);
    
}*/
@Override
public EmissionStateSpace getEmissionStateSpace() {
   return emStSp;
}
public void set(int i, Double sampleR, Double sampleB) {
   // if(true) throw new RuntimeException("!!!!");
   this.setR(i, sampleR);
   this.setB(i, sampleB);
    
}


public double[] getEmiss(int i){
	 throw new RuntimeException("!!");
//    return this.emissions[i].calcDistribution( distribution, emStSp, i);
}
}


