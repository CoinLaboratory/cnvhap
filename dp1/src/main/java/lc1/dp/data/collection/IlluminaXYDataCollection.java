package lc1.dp.data.collection;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;

import lc1.dp.data.representation.Emiss;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.PhasedDataState;
import lc1.stats.IlluminaRDistribution;
import lc1.util.Constants;

public class IlluminaXYDataCollection extends IlluminaRDataCollection {
	 public  IlluminaXYDataCollection (File f, short index,
			 int no_copies, int[][] mid,  File bf, Collection<String> snpidrest) throws Exception{
	        super(f, index, no_copies, mid, bf,snpidrest);
	        
	    }
	 public IlluminaXYDataCollection(IlluminaRDataCollection collection) {
	      super(collection);
	      this.dc = collection.dc;
	   //   this.distR = collection.distR;
	    //  this.distB = collection.distB;
	    }
	 
	 
	 /*NOTE THIS IS A BIT OF AN ABUSE OF FLIP STRAND TO DEAL WITH SEQUENOM ISSUES.  FLIPPING STRAND SHOULD NOT 
	  * NECESSARILY CHANGE BAF(non-Javadoc)
	  * @see lc1.dp.data.collection.DataCollection#flipStrand(int)
	  */
	 public void flipStrand(int i){
    	if(alleleA.size()>0){ 
    		alleleA.set(i, compl(this.alleleA.get(i)));
    	    alleleB.set(i, compl(this.alleleB.get(i)));
    	}
    	 this.baf.set(i, 1-this.baf.get(i));
    	 Boolean strand_ = strand.get(i);
    	 if(strand_!=null && strand_) {
    		 throw new RuntimeException("should only flip to plus strand");
    	 }
    	 for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
    		 HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
    		 hes.emissions[i].swtchAlleles();
    	 }
    	 strand.set(i, true);
    }
	 
	 @Override
	    public  Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy, double[] missing){
	        boolean doneB = false;
	        boolean nullB = false;
	        boolean nullR = false;
	        IlluminaNoBg state = (IlluminaNoBg)this.dataL.get(indiv);
	 //       boolean forward = this.strand.get(i);
	      if(header[4].startsWith("area")){
	    	  double x = Double.parseDouble(geno[geno.length-2]);
	    	  double y = Double.parseDouble(geno[geno.length-1]);
	    	  if(x+y < 1){
	    	   ((IlluminaNoBg)state).emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);
	    	   return false;
	    	  }
	    	  ((IlluminaNoBg)state).setR(i,Math.log((x+y)/2.0));
	    	    ((IlluminaNoBg)state).setB(i, y/(y+x)  );
	    	    
	    	    missing[2]+=x;
	    	    missing[3]+=y;
	    	    return false;
	      }
	        
	        
	       for(int k=0; k<header.length; k++){
	           if(geno[k].equals("null")) geno[k] = "NaN";
	             String hk = header[k].toLowerCase();
	             if(hk.indexOf("y raw")>=0   || hk.indexOf("y_raw")>=0){
	                 Double b = Double.parseDouble(geno[k])/1000.0;
	                 if(b==null || Double.isNaN(b)) nullB  = true;
	                ((IlluminaNoBg)state).setB(i,b);
	                doneB =true;
	                if(b<Constants.minB[index]) b = Constants.minB[index];
	                 if(b>Constants.maxB[index]) b = Constants.maxB[index];
	                 
	            }
	             else if((hk.indexOf("x raw")>=0 || hk.indexOf("x_raw")>=0)){
	                 double r = Double.parseDouble(geno[k])/1000.0;
	                 if( Double.isNaN(r)) nullR  = true;
	               
	                // if(this.snpid.get(i).startsWith("cnv1950")) r+=0.6;
	                 state.setR(i,r);
	                 if(r<Constants.minR[index]) r = Constants.minR[index];
	                 if(r>Constants.maxR[index]) r = Constants.maxR[index];
	                
	                 
	             }
	             else if(false && header[k].indexOf("AgilentPred")>=0){
	              /*   EmissionStateSpace emStSp = state.getEmissionStateSpace();
	                 PhasedDataState pig = (PhasedDataState)this.data.get(indiv);
	                 if(geno[k].equals("0")){
	                     pig.emissions[i] = new IntegerDistribution(super.trans("AA"));
	                 }
	                 else{
	                     IlluminaRDistribution dist = 
	                    	 Constants.suppressB() ?
	                    			 new IlluminaRDistribution(this.index):
	                    	 new IlluminaDistribution(this.index);
	                     double r = Double.parseDouble(geno[k]);
	                     dist.setR(r);
	               //      IlluminaProbR probR = this.getR();
	                //     IlluminaProbB probB = this.getB();
	                    // double[] d = dist.calcDistribution( state.distribution, emStSp);
	                  //   int max = Constants.getMax(d);
	                   
	                   
	                  //  pig.emissions[i] = new IntegerDistribution(max);
	                 }*/
	             }
	                    }
	       if(nullR ){
	    	   ((IlluminaNoBg)state).emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);
	    	   
	       }
	       else if(!doneB
	    		   || snpid.get(i).toLowerCase().indexOf("cnv")>=0
	    		   || snpid.get(i).toLowerCase().indexOf("A_")>=0
	    		 || nullB
	    		   ){
	         Number r =  ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).r();
	        ((IlluminaNoBg)state).emissions[i] = new IlluminaRDistribution(this.index);
	        ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).setR(r.doubleValue());
	         //  throw new RuntimeException("!!");
	       }
	       try{
	    	      //  boolean doneGeno = false;
	    	        PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
	    	        for(int k=0; k<header.length; k++){
	    	        /*if(header[k].toLowerCase().indexOf("geno")>=0){
	    	            if(data!=null && data.emissions[i]==null){
	    	            int ind = trans(geno[k]);
	    	        
	    	                 data.emissions[i] = new IntegerDistribution(ind);
	    	            }
	    	        }*/
	    	        
	    	    }
//	    	        if(data.emissions[i]==null){
	  //  	            EmissionState sta = this.dataL(indiv);
	    //	             data.emissions[i] = new IntegerDistribution(sta.getBestIndex(i));
	    	//         }
	    	  //      data.emissions[i].setDataIndex(this.index);
	    	       
	    	    }catch(Exception exc){
	    	        exc.printStackTrace();
	    	    }
	    	    Boolean res =  state.emissions[i].probeOnly();
	    	 /*  if(res==null){
	    		   int pos = this.loc.get(i);
	    	   	Logger.global.info("all NC at "+i+" "+this.index+" "+pos);
	    	   }*/
	    	    if(nullR) missing[0]++;
	    	    if(nullB) missing[1]++;
	    	    return res;
	    }
}
