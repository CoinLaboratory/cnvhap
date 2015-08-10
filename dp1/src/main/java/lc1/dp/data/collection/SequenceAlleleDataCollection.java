package lc1.dp.data.collection;

import java.io.File;
import java.util.Collection;

import lc1.dp.data.representation.Emiss;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.IlluminaDistribution;
import lc1.util.Constants;
import cern.jet.stat.Probability;

public class SequenceAlleleDataCollection extends SequenceDataCollection {

	
	

    public SequenceAlleleDataCollection(File file, short i, int j, int[][] mid,
			File buildF, Collection<String>str) throws Exception{
		super(file,i,j,mid,buildF, str);
	}

	public int[] countAB(String ab){
    	int[] count=new int[]{0,0};
    	for(int i=0; i<ab.length(); i++){
    		if(ab.charAt(i)=='A') count[0]++;
    		else count[1]++;
    	}
    	if(count[0]==0 && count[1]==0) throw new RuntimeException("no alleles");
    	return count;
    }

    public double homoProb(char errorAllele, String quality, String base){
    	double prob=1.0;
    	for(int i=0; i<base.length(); i++){
    		if(base.charAt(i)==errorAllele) prob*=Math.pow(10, (-((int)quality.charAt(i)-33)/10));
    	}
    	return prob;
    }
    
    public double aveProb(String quality){
    	double prob=0.0;
    	for(int i=0; i<quality.length(); i++){
    		prob+=(1-Math.pow(10, (-((int)quality.charAt(i)-33)/10)));
    	}
    	return prob/(double)quality.length();
    }
    public double[] normalize(double[] val){
    	double[] nor=new double[val.length];
    	double sum=0.0;
    	for(int i=0; i<val.length; i++){
    		sum+=val[i];
    	}
    	for(int i=0; i<nor.length; i++){
    		nor[i]=val[i]/sum;
    	}
    	return nor;
    }
	public double[] getGenoProb(String[] geno, int stspSize){
    	double[] probs=new double[stspSize];
        	double[] temp=new double[stspSize];
    		int[] countAB=this.countAB(geno[0]);
    		if(countAB[0]==0 && countAB[1]==0) throw new RuntimeException("no allele data");
    		if(countAB[0]!=0 && countAB[1]!=0) { //AB
    			temp[0]=this.homoProb('B', geno[1], geno[0]);
    			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5)-Probability.binomial(countAB[0]-1, countAB[0]+countAB[1], 0.5);
    			temp[2]=this.homoProb('A', geno[1], geno[0]);
    		}
    		else if(countAB[0]==0){ //BB
    			temp[0]=this.homoProb('B', geno[1], geno[0]);
    			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5);
    			temp[2]=this.aveProb(geno[1]);    			
    		}
    		else if(countAB[1]==0){//AA
    			temp[0]=this.aveProb(geno[1]); 
    			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5)-Probability.binomial(countAB[0]-1, countAB[0]+countAB[1], 0.5);
    			temp[2]=this.homoProb('A', geno[1], geno[0]);    			
    		}
    		else throw new RuntimeException("!!");
    		probs=this.normalize(temp);

    	return probs;
    }
	@Override
	protected void setMinMaxRValues() {
		for(int i=0; i<this.minR.size(); i++){
		 if(dc!=null) ((DistributionCollection)dc).probRB.setMinMax(minR.get(i), maxR.get(i),
				 Constants.minB(this.index), Constants.maxB(this.index),i);
		}
		
	}
	 //@Override
	 public Boolean process(String indiv, String[] header, String[] geno, int i,
			 int ploidy,double averageDepth, int ind_indiv) {

		/* boolean donea = false;
		 boolean nulla = false;
		 boolean doneb = false;
		 
		 int numPoints=0;*/
		 HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);

		 Integer aCount=null;
		 Integer bCount = null;
		 if(geno.length==1 && geno[0].equals("./.")){
			aCount=0;bCount=0; 
		 }
		 else{
		 for(int k=0; k<header.length; k++){	   		
			 String hk = header[k].toLowerCase();
			 if(geno[k].startsWith("nul")) geno[k] = "NaN";
			 if(hk.indexOf("counta")>=0){
				
				aCount = Integer.parseInt(geno[k]);
				 // if(depth>0){
				 //	 System.err.println(depth);
				 // }
			 }
			 else if(hk.indexOf("countb")>=0){
					bCount = Integer.parseInt(geno[k]);
			 }else if(hk.equals("ad")){
				 if(geno[k].equals(".")) geno[k] = "0,0";
				 String[] str = geno[k].split(",");
				 aCount = Integer.parseInt(str[0]);
				 bCount = Integer.parseInt(str[1]);
			 }
		 }
		 }
		 int totDepth = aCount+bCount;
		 double lrr = Double.NaN; double baf = Double.NaN;
		 if(totDepth>=Constants.depthThresh()){
			 //lrr = Constants.adjustR( Math.log((double)totDepth / ((double)averageDepth))/Math.log(2),index);
			 //baf = totDepth==0 ? Math.random() : (double)bCount/(double)totDepth;
			
			 lrr = totDepth;
			 baf = bCount;
			// System.err.println(lrr+" "+baf);
			 if(lrr <minR.get(i)) minR.set(i, lrr);
			 if(lrr> maxR.get(i)) maxR.set(i,lrr);
		 } 
		 
		 if(Double.isNaN(totDepth) || Double.isNaN(lrr)){
			
			 state.emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(null);

		 }
		 else{
			 IlluminaDistribution ill = new IlluminaDistribution(index);
			 state.emissions[i] = ill;//
			 ill.setR(lrr);
			 ill.setB(baf);
				 //new lc1.stats.AlleleCountDistribution(this.index,aCount,bCount,averageDepth/(double)ploidy, ind_indiv);

		 }
		 return false;		




		 //		return super.process(indiv, header, geno, i, ploidy);
	 }
	
	/* @Override
	    public Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy, double[] missing){
	        boolean doneB = false;
	        boolean doneR = false;
	        boolean nullB = true;
	        boolean nullR = true;
	        IlluminaNoBg state = (IlluminaNoBg)this.dataL.get(indiv);
	      
	       for(int k=0; k<geno.length; k++){
	           if(geno[k].startsWith("nul")) geno[k] = "NaN";
	             String hk = header[k].toLowerCase();
	             
	            if(!doneB && (hk.indexOf("countb")>=0)){
	                 Double b = Double.parseDouble(geno[k]);
	                 nullB  = (b==null || Double.isNaN(b));
	              
	                doneB =true;
	                if(!nullB){
		                 if(b<Constants.minB[index]){
		                	b = Constants.minB[index];
		                 }
		                 if(b>Constants.maxB[index]) {
		                	b =  Constants.maxB[index];
		                 }
		                 ((IlluminaNoBg)state).setB(i,b);
		                }
	                 
	            }
	             else if(
	            		 !doneR && 
	(            		 (header[k].indexOf("counta")>=0 ))){
	                 double r = Double.parseDouble(geno[k]);
	                 nullR  = Double.isNaN(r);
	                if(!nullR){
	                 if(r<Constants.minR[index]){
	                	r = Constants.minR[index];
	                 }
	                 if(r>Constants.maxR[index]) {
	                	r =  Constants.maxR[index];
	                 }
	                 state.setR(i,r);
	                }
	                doneR = true;
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
	       if(r!=null) ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).setR(r.doubleValue());
	         //  throw new RuntimeException("!!");
	       }
	       }
	      
	    	    Boolean res =  state.emissions[i].probeOnly();
	    	
	    	 if(nullB) missing[1]++;
	    	 if(nullR) missing[0]++;
	    	    return res;
	    }*/
	
	    
	
}
