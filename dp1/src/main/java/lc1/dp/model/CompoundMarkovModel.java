package lc1.dp.model;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.util.Constants;

public abstract class CompoundMarkovModel extends MarkovModel{
    public abstract MarkovModel getMarkovModel(int i);
    public abstract int noCopies();
  public CompoundMarkovModel(MarkovModel m) {
        super(m);
    }
  @Override
  public String info() {
      StringBuffer sb = new StringBuffer();
      MarkovModel[] m = this.getMemberModels();
      for(int i=0; i<m.length; i++){
         sb.append(m[i].info()+"\n");
      }
      return sb.toString();
   }
  @Override
  public EmissionStateSpace getStateSpace(){
 	
 		  CompoundMarkovModel hmm = this;
 		 while(hmm instanceof WrappedModel){
             hmm = ((WrappedModel)hmm).getHMM();
            } 
 		
 		return  ((PairMarkovModel)hmm).emstsp;
 	//
 		
 	  
   }
  
  public static EmissionStateSpace getStateSpace(MarkovModel[] mm, int[] nocops, boolean onlyRepeats){
	 
	//	this.
		int sum =0;
		for(int i=0; i<nocops.length; i++){
			sum+=nocops[i];
		}
		EmissionStateSpace[] em = new EmissionStateSpace[sum];
	    sum =0;
		for(int i=0; i<nocops.length; i++){
			for(int j=0; j<nocops[i]; j++){
				em[sum+j] = mm[i].getStateSpace();
			}
			sum+=nocops[i];
		}
		return  new CompoundEmissionStateSpace(em, onlyRepeats, Constants.parentobj!=null);
  }
  @Override
  public boolean converged(){
      MarkovModel[] m = this.getMemberModels();
      for(int i=0; i<m.length; i++){
          if(!m[i].converged()) return false;
      }
      return true;
  }
  @Override
  public void validate(int length) throws Exception{
      try{
          super.validate(length);
      }
      catch(Exception exc){
         exc.printStackTrace();
        for(int i=0; i< this.noCopies(); i++){
            getMarkovModel(i).validate(length);
        }
         System.exit(0);
      }
      
  }
  public boolean allOneLength(){
     return this.emstsp.allOneLength();
  }
  
    public CompoundMarkovModel(String name, int noSnps){
        super(name,   noSnps);
    }
   /* public double[] getPseudoCountWeights(){
        return this.getMarkovModel(0).getPseudoCountWeights();
    }*/
    
    public abstract CompoundState getCompoundState(State[] res);
  //  public abstract Set<Integer> getHemizygous(int emissionIndex);
    public abstract void refresh();
  //  public abstract CompoundState disambiguate(CompoundState state, CompoundState prev) ;
    public abstract State[] disambiguate(State[] memberStates, State[] prev, int index, boolean sample) ;
    public abstract EmissionState disambiguate(EmissionState state, EmissionState previous, int indexOfToEmiss, boolean sample, int j);
    public abstract MarkovModel[] getMemberModels();
    
    @Override
    public void print(PrintWriter pw, List<Integer>columns, int popsize){
        try{
       //    for(int i=0; i<Constants.format.length; i++){
         //      this.probB(i).print(pw);
         //  }
        Map<String, MarkovModel> m = new HashMap<String, MarkovModel>();
            for(int i=0; i<this.noCopies(); i++){
                MarkovModel mm =  this.getMarkovModel(i);
                m.put(mm.getName(), mm);
            }
            for(Iterator<MarkovModel> it = m.values().iterator(); it.hasNext();){
                MarkovModel mm = it.next();
                if(mm instanceof HaplotypeHMM && Constants.writeHMM()>=1  && Constants.writeFree()){
                     mm = new FreeHaplotypeHMM((FreeHaplotypeHMM)mm, true);
                }
                 mm.print(pw, columns, popsize);
                pw.println("####################################################");
            }
        }catch(Exception exc){
            exc.printStackTrace();
        }
     }
    public CompoundMarkovModel unwrapModel() {
       return this;
    }
    
	/* public void printProbR(PrintWriter pw) throws Exception{
    String fStr = "%5.3g";
    	int len = probR.get(0).length();
    	for(int k=0; k<len; k++){
    	
    		{
		    		pw.println("index = "+Constants.inputDir[k]+" probR");
		    	Double[] mean = new Double[probR.size()];
		    	Double[] var = new Double[probR.size()];
		    	Double[] skew = new Double[probR.size()];
		    	 for(int i=0;i<mean.length; i++){
		    		 IlluminaRArray pr1 = probR.get(i);
		    		
		    			 SkewNormal sn =(SkewNormal) ((Mixture)pr1.get(k)).dist[0];
		    			 if(sn!=null){
		    			 mean[i] = sn.location;
		    			 var[i] = sn.scale;
		    			 skew[i] = sn.shape;
		    			 }
		    			 else{
		    				 mean[i] = Double.NaN;
		    				 var[i] = Double.NaN;
		    				 skew[i] = Double.NaN;
		    			 }
		    		 
		    	 }
		    		 pw.println(sprintf(fStr, mean,";",""));
		    		 pw.println(sprintf(fStr, var,";",""));
		    		 pw.println(sprintf(fStr, skew,";",""));
		    		 pw.println();
    		}
    		 
    		IlluminaProbB probB = probB(k);
    		 pw.println("index = "+Constants.inputDir[k]+" single R");
    		List<State> states = this.getMarkovModel(0).states;
    		 Double[] mean = new Double[states.size()-1];
		    	Double[] var = new Double[states.size()-1];
		    	Double[] skew = new Double[states.size()-1];
    		    	for(int i=1; i<states.size(); i++){
    		    		HaplotypeEmissionState hes = (HaplotypeEmissionState)states.get(i);
    		    		
    		    			SkewNormal sn = (SkewNormal) ((Mixture)hes.probR().get(k)).dist[0];
    		    			 if(sn!=null){
    		    			mean[i-1] = sn.location;
    		    			var[i-1] = sn.scale;
    		    			skew[i-1]  = sn.shape;
    		    			 }
    		    			 else{
    		    				 mean[i-1] = Double.NaN;
    		    				 var[i-1] = Double.NaN;
    		    				 skew[i-1] = Double.NaN;
    		    			 }
    		    		
    		    	}
    		    	 pw.println(sprintf(fStr, mean,";",""));
		    		 pw.println(sprintf(fStr, var,";",""));
		    		 pw.println(sprintf(fStr, skew,";",""));
		    		 pw.println();
    	    
    		 
    		 pw.println("index = "+Constants.inputDir[k]+" probB");
    		 
    		 
    		 
    		 probB.print1(pw);s
    		 
    		 
    	
    	}
    }	*/
 public static String sprintf(String format, Double[] mean, String spl, String preface) {
		StringBuffer sb = new StringBuffer(preface);
	
		for(int i=0; i<mean.length; i++){
			if(i>0)sb.append(spl);
			if(mean[i]<0 && spl.equals(":")) sb.append("^");
			sb.append(String.format(format,mean[i]));
		}
		return sb.toString().replaceAll(" ", "");
	}
	//List<IlluminaRArray> probR ;

   
   // List<IlluminaRArray> probB = new ArrayList<IlluminaRArray>();
  //  public abstract SimpleExtendedDistribution[] probHomoIsHemi();
   // public abstract IlluminaProbB probB(int i) ;
       
       /* public void  probR(List<List<Integer>> stateIndices, List<Integer> cnvIndex, List<ProbabilityDistribution> s, List<ProbabilityDistribution[]> s1, int ind){
        	probR = new ArrayList<IlluminaRArray>();
        	for(int i=1; i<this.states.size(); i++){
                    EmissionState st = (EmissionState) states.get(i);
                    IlluminaRArray pR = ((EmissionState)st).probR();
                    this.probR.add(pR);
                    ProbabilityDistribution probR = pR.get(ind);
                 //   if(probR==null) continue;
                    int index = s.indexOf(probR);
                    if(index<0){
                        s.add(probR);
                        List<Integer> l = new ArrayList<Integer>();
                        l.add(i);
                        if(stateIndices!=null){
                            stateIndices.add(l);
                            int nocop =  ((EmissionState)states.get(i)).noCop();
                            cnvIndex.add(nocop);
                        }
                        if(s1!=null && st instanceof CompoundState){
                            s1.add( ((CompoundState)st).getProbabilityDist(ind));
                          
                        }
                    }
                    else if(stateIndices!=null){
                        stateIndices.get(index).add(i);
                        if(((EmissionState)states.get(i)).noCop().intValue()!=cnvIndex.get(index) && Constants.trainEnsemble()>0){
                          throw new RuntimeException("!!");
                        }
                    }
            }
        	Collections.sort(probR);
        	//System.err.println(probR);
        }*/
       
}
