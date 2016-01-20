package lc1.dp.model;

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.AlleleCopyPairEmissionState;
import lc1.dp.transition.FreeRateTransitionProbs1;
import lc1.dp.transition.FreeTransitionProbs1;
import lc1.dp.transition.MultiExpProbs;
import lc1.dp.transition.SimpleExpTransProb;
import lc1.stats.ProbabilityDistribution;
import lc1.util.Constants;

public class HaplotypeHMMIterator implements Iterator<MarkovModel> {
final int noSites;
final int count;
int i=0;
final double[] init;
//final EmissionStateSpace emStSp;
 int numF;
final String[][] modify; 
//final double[][] meanvarskew;
final double[] modifyFrac, modifyFracStart;
//1.0 means that all emissions are 1.0 for a particular category
//-1.0 means that it all emissions are within a copy class, but distributed over that copy class
//0.0 means that emissions are random
//Math.abs(x) > 0 is intermediate

final List<Integer> locs;
final DataCollection datac;
public HaplotypeHMMIterator(String name, int noSites, int count,  DataCollection datac, int numF, int[] mode,
        String[][] modify, double[] modifyFrac, 
        double[] modifyFracStart,
        ProbabilityDistribution[] numLevels, List<String> model){
    this.noSites = noSites;
    this.modelnames = model;
    this.datac = datac;
  //  this.meanvarskew = meanvarskew;
    this.name = name;
    this.count = count;
  //  this.emStSp =emStSp;
    this.locs =datac.getLocations();
    this.modifyFracStart = modifyFracStart;
    EmissionStateSpace emStSp = Emiss.mergedSpace;
    init = new double[emStSp.size()]; 
    Arrays.fill(init, 1.0 / (double)emStSp.size());
    setMode(mode);
    this.numF = numF;
    this.modify = modify;
    this.modifyFrac = modifyFrac;
    if(modify[0].length!=modifyFrac.length) {
  	  throw new RuntimeException("!!" +modify.length+" "+modify.length);
    }
    this.numLevels = numLevels;
}
public static     Iterator<MarkovModel> mergeModel(final Iterator<MarkovModel> counts,
		final Iterator<MarkovModel> alleles) {
	  return  new Iterator<MarkovModel>(){

         public boolean hasNext() {
             return alleles.hasNext() && counts.hasNext();
         }

         public MarkovModel next() {
             MarkovModel cnt=counts.next();
             MarkovModel allele = alleles.next();
            
            return 
          //  new CachedHMM(
            new PairMarkovModel(new MarkovModel[] {cnt,allele }, new int[] {0,1}, AlleleCopyPairEmissionState.class,
         false ,false      
            //)
            );
         }

         public void remove() {}
            
        };
     }


/*
 0 == ExponentialTransitionsProbs
 1 = FreeTransitionProbs1
  2 = ExponentialAndFreeTransitionsProbs
 3 == ExponentialMultiExpTransitionProbs
 4 = ExponentialAndFreeTransitionProbsMulti
 5 = FreeTransitionProbs
*/
Object[] clazz;
//int index;
double[] u_g = Constants.u_global(0);//new double[] {Constants.u_global(),Constants.u_global(), Double.POSITIVE_INFINITY, 100}; //emiss, trans

public static Class[][] getMode(int [] transMode){
	Class[][] clazz = new Class[transMode.length][2];
    for(int i=0; i<transMode.length; i++){
        if(transMode[i]==0){
            clazz[i] =new Class[] { SimpleExpTransProb.class, SimpleExpTransProb.class };  //first is exp, second is alpha
        }
    else if
       ( transMode[i]==1){
        clazz[i] =  new Class[] {FreeRateTransitionProbs1.class};
    }
    else if
       ( transMode[i]==2){
        clazz[i] =  new Class[] { SimpleExpTransProb.class, MultiExpProbs.class };
        
    }
    else  if(transMode[i]==3  ){
        clazz[i] =  new Class[] {  MultiExpProbs.class , SimpleExpTransProb.class};
    }
    else  if(transMode[i]==4  ){
        clazz[i] =   new Class[] {  MultiExpProbs.class , MultiExpProbs.class};
    }
    else  if(transMode[i]==5  ){
    	 clazz[i] =  new Class[] {FreeTransitionProbs1.class};
       // throw new RuntimeException("!!");
 //       clazz[i] =  FreeTransitionProbs.class;
    }
        Logger.global.info("clazz "+Arrays.asList(clazz));
    }	
    return clazz;
}

public void setMode(int[] transMode){
    clazz = getMode(transMode);
 /* if(transMode==5){
      u[1] = Double.POSITIVE_INFINITY;
      u[2] = Double.POSITIVE_INFINITY;
  }
    else{
        u[1] = Constants.u_global;
        u[2] = Double.POSITIVE_INFINITY;
    }*/
}
final String name;
final ProbabilityDistribution[]  numLevels;

List<String> modelnames=null; 
public  MarkovModel next(){
     i++;
   // int index;
   //  Class clazz;
      Double[] r = 
              new Double[] {   Constants.expModelIntHotSpot(0)*Constants.probCrossOverBetweenBP,
              Constants.expModelIntHotSpot(1)*Constants.probCrossOverBetweenBP};
    FreeHaplotypeHMM res = null;
    if(modelnames!=null){
    	/*res = new FreeHaplotypeHMM(datac, modelnames, clazz, r, 
                null,
                modifyFrac, modifyFracStart) ;*/
    	throw new RuntimeException("!!");
    }
    res = new FreeHaplotypeHMM(name+"_free_haplotype_"+i, numF, noSites, init, 
    		Emiss.mergedSpace, clazz, locs, r, 
            null,
            modify, modifyFrac, modifyFracStart, false, ((DataCollection)datac).probeOnly, numLevels) ;
 if(Constants.useHMMFile(0)!=null){
	 
	
	 try{
		 if(Constants.useHMMFile(0).endsWith(".zip")){
	 ZipFile hmmFile = new ZipFile(new File(Constants.useHMMFile(0)));
	 res.modify(hmmFile, Constants.transferHMM(), this.locs);
	 hmmFile.close();
		 }
		 else{
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(new File(Constants.useHMMFile(0)))); 
			res = (FreeHaplotypeHMM) ois.readObject();
		 }
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
 }
    if(Constants.CHECK){
      try{
      res.validate(res.noSnps);
      }catch(Exception exc){
          exc.printStackTrace();
          System.exit(0);
      }
  }
 
        
               return  res;
    }

public boolean hasNext() {
   return i<count;
}

public void remove() {
    // TODO Auto-generated method stub
    
}
}
