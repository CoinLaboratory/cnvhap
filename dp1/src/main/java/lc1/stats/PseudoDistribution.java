package lc1.stats;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public abstract class  PseudoDistribution implements ProbabilityDistribution{
public void setCoverage(double d){
    	throw new RuntimeException("!!");
    }
    short data_index = -1;
    
    public abstract double[] probs();
    public  double[] pseudo(){
    	throw new RuntimeException("!!");
    }
    public void setPriors(lc1.stats.ProbabilityDistribution distx, int type){
		//this.dist[0].setPriors(((Mixture)distx).dist[0]);
	}
    public int fillVariance(DoubleMatrix2D y,  int numObs, double pseudo){
    	throw new RuntimeException("!!");
    }
    public abstract void setProb(double[] prob);
    public void getInterval(double[] in, double[] res_inR){}
    public abstract void transfer(double pseudoC1);
   
    public abstract void addCount(int obj_index, double value);

    public abstract void initialise();

   public double scale(){
	   return 0.0;
   }

   public int[] getOrder(){
   	throw new RuntimeException("!!");
   }
	@Override
	public void setParam(int type, double rho) {
		// TODO Auto-generated method stub
		
	}
    public abstract double probs(int obj_i);
public void print(PrintWriter pw){
	pw.print(this.toString()+"\t");
}
    public abstract double sum();

    public abstract int getMax();

    public abstract double[] counts();
    public void variance( double[] sum){}
    public abstract PseudoDistribution clone();

    public abstract void validate();
    public  void validate(boolean normalise){
    	this.validate();
    }
    @Override
	public int numObs() {
		return 0;
	}

    public abstract void setProbs(int to, double d);

    public abstract String getPrintString();

    public  void printSimple(PrintWriter pw, String string, String string2, double thresh){
        pw.print(string+" "+string2+" "+getPrintString());
    }

    public void print(PrintWriter pw, boolean b, String printString, String string){
        pw.print(printString+" "+string+" "+getPrintString());
    }
    public abstract PseudoDistribution clone(double swtch) ;

    public abstract double logProb();

    public abstract void setCounts(int i1, double cnt);

    public abstract Integer fixedInteger();

    public  abstract PseudoDistribution swtchAlleles();

    
   

    public  short getDataIndex(){
        return data_index;
    }

    public abstract double[] calcDistribution(
           double[] distribution, EmissionStateSpace emStSp, int pos) ;//{
     // return this.probs();
    //}
    
    public double[] calcDistribution( double[] distribution, EmissionStateSpace emStSp, int nocop, int pos){
        double sum=0;
        double[] probs = probs();
        Arrays.fill(distribution, 0.0);
        for(int j=0; j<distribution.length; j++){
                   if(emStSp.getCN(j)!=nocop){
                       distribution[j] = 0;
                   }
                   else{
                       distribution[j] = probs[j];
                   }
                   sum+=distribution[j];
                }
        for(int i=0; i<distribution.length; i++){
            distribution[i] = distribution[i]/sum;
        }
        return distribution;
    }

    public abstract void  setFixedIndex( int k);

    
    public String getUnderlyingData(EmissionStateSpace emStSp){
        int ind = Constants.getMax(this.probs());
        Comparable compa = emStSp.get(ind);
        if(compa instanceof ComparableArray){
            ComparableArray arr = (ComparableArray)compa;
            StringBuffer sb = new StringBuffer();
            for(int i1=0; i1<arr.size(); i1++){
                sb.append(arr.get(i1));
            }
           return sb.toString();
        }
        else return ((Emiss)compa).toStringPrint();
        //+"_"+this.probs()[ind];
    }
   /*public double calcDistribution(PseudoDistribution hmm_state, 
                            SkewNormal[] probRState, 
                            IlluminaProbB[] probBState, 
                            double[] distribution, EmissionStateSpace emStSp){
        double sum=0;
        Arrays.fill(distribution, 0.0);
        for(int j=0; j<distribution.length; j++){
                double prob_j =hmm_state.probs(j);
                if(prob_j>0){
                  
                   distribution[j] =  prob_j 
                       * this.score(j, probBState,probRState)*emStSp.getWeight(j);  //do we include the weight term???;
                   sum+=distribution[j];
                }
        }
        return sum;
    }*/

    public void setDataIndex(short data_index2) {
    	if(data_index==-2) return;
    
       this.data_index = data_index2;
        
    }
    
    public void transfercounts(EmissionState innerState, int phen_index, int i){
        double[] countsData = counts();
        for(int k=0; k<countsData.length; k++){
            if(countsData[k]==0) continue;
                innerState.addCountDT((double)k,phen_index,  countsData[k], i);
        }
    }

    public abstract double totalCount() ;

    public void applyCorrection(double parseDouble) {
        // TODO Auto-generated method stub
        
    }
    
    public void applyStdErrCorrection(double stderrmult){
    
    }
       
    public void addCount(double obj_index, double value) {
        this.addCount((int) Math.round(obj_index), value);
        // TODO Auto-generated method stub
        
    }
    public double probability(double x) {
       return this.probs((int)Math.round(x));
    }

    public double[] getCount(double[] angle){
        throw new RuntimeException("!!");
    }

   
    public  String getCompressedDataString(EmissionStateSpace emStSp){
      //  return emStSp.get(this.getMax()).toString().replaceAll("[,\\s+\\[\\]]","");
         throw new RuntimeException("!! "+this.getClass());
     }

    public String compressedStringHeader(EmissionStateSpace emStSp) {
        return "Genotype";
       }

    public void compareTo(PseudoDistribution pseudoDistribution) {
        // TODO Auto-generated method stub
        
    }

	public void getDataIndices(Set<Short> res) {
		if(data_index>=0) res.add(this.data_index);
		
	}

	public void change(Sampler dir) {
		// TODO Auto-generated method stub
		
	}

	public double getUncertainty() {
		throw new RuntimeException("!!");
	}

	public final double[] fillCN(EmissionStateSpace emissionStateSpace,  double[] d) {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean containsIndex(int index) {
		return this.getDataIndex()==index;
	}

	
   
            
    
	public double getMean() {
		// TODO Auto-generated method stub
		return 0;
	}


	public int getParamIndex() {
		// TODO Auto-generated method stub
		return 0;
	}


	public double getParamValue(int n1) {
		// TODO Auto-generated method stub
		return 0;
	}


	public String id() {
		// TODO Auto-generated method stub
		return null;
	}


	public void maximise(double d, double e, double f) {
		// TODO Auto-generated method stub
		
	}


	public String name() {
		// TODO Auto-generated method stub
		return null;
	}


	public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch) {
		// TODO Auto-generated method stub
		return 0;
	}


	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		// TODO Auto-generated method stub
		
	}


	public double prior() {
		// TODO Auto-generated method stub
		return 0;
	}


	public void recalcName() {
		// TODO Auto-generated method stub
		
	}


	public void setParamValue(int n1, double val) {
		// TODO Auto-generated method stub
		
	}


	public void updateParamIndex() {
		// TODO Auto-generated method stub
		
	}


	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return 0;
	}


	public double getLowerBound(int n) {
		// TODO Auto-generated method stub
		return 0;
	}


	public int getNumArguments() {
		// TODO Auto-generated method stub
		return 0;
	}


	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}


	public double getUpperBound(int n) {
		// TODO Auto-generated method stub
		return 0;
	}
    
	 public double evaluate(double[] argument) {
			// TODO Auto-generated method stub
			return 0;
		}

	public  Double getIntensity(boolean r){
		return Double.NaN;
	}
	public abstract boolean isMeasured();
	
	public abstract double  scoreB(int j, int i) ;
	
	public void addBCount(int j, double val, int i) {
		// TODO Auto-generated method stub
		
	}
	public void addRBCount(EmissionStateSpace emstsp, int j, double val, int i) {
		//throw new RuntimeException("!!");
	}


	public void addRCount(int bg, int no_cop, double val, int i) {
		// TODO Auto-generated method stub
		
	}

	
	public double scoreR(int bg, int no_cop, int i) {
		// TODO Auto-generated method stub
		return 1.0;
	}
	public abstract double scoreBR(EmissionStateSpace emStSp, int j,   int i);
	//	// TODO Auto-generated method stub
	//	return 1.0;
	//}
	public Boolean probeOnly() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public double  probability(double x, int mixComponent){
		throw new RuntimeException("!!");
	}
	
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[] noCop, double pseudo){
	throw new RuntimeException("!!");
	}
	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1){
		throw new RuntimeException("!!");
	}
	public double scoreBR(EmissionStateSpace emstsp, int j, int i, int mixComponent) {
		return scoreBR(emstsp,j, i);
	}
	public void getB(double[] b, double mult) {
		// TODO Auto-generated method stub
		//return 0;
	}
	public void multiplyCounts(double d) {
		// TODO Auto-generated method stub
		
	};

	public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, ProbabilityDistribution disty){
		throw new RuntimeException("!!");
	}
	public void setMinMax(double min, double max){
		throw new RuntimeException("!!");
	}
	public double evaluate(double d) {
		// TODO Auto-generated method stub
		return 0;
	}
	public double evaluate(DoubleMatrix1D viewRow, double d) {
		// TODO Auto-generated method stub
		return 0;
	}
	public final  void check(PseudoDistribution dist2, EmissionStateSpace emstsp){
		
		Double b2 = this.getBEst(emstsp);
		    	if(b2!=null){
		    	
		    		
		    		Double b1 = dist2.getBEst(emstsp);
		    		
		    		if(b1!=null && Math.abs(b2 - b1)>Constants.alleleDiffThresh()){
						throw new RuntimeException("!! problem alleles seem to be inconsistent");
					}
		    		
		    	}
		    	
		
	}
	protected Double getBEst(EmissionStateSpace emstsp) {
		
		Integer fixed = fixedInteger();
		if(fixed==null){
			fixed = Constants.getMax(probs());
		}
		return (double)emstsp.getBCount(fixed)/(double)emstsp.getCN(fixed);
	}
	public double weight() {
		// TODO Auto-generated method stub
		return 1.0;
	}
	public void setEmStSp(EmissionStateSpace stsp) {
	throw new RuntimeException("!!");
		
	}
	public EmissionStateSpace getEmissionStateSpace(){
		throw new RuntimeException("!!");
	}
	public void modify(PseudoDistribution dist1) {
		//SimpleExtendedDistribution dist1 = pseudoDistribution;
		//EmissionStateSpace emstsp = emstsp;
		EmissionStateSpace emstsp1 = dist1.getEmissionStateSpace();
		EmissionStateSpace emstsp = this.getEmissionStateSpace();
		for(int i=0; i<emstsp.genoListSize(); i++){
			this.setProbs(i, 0);
		}
		double sum = 0;
		for(int i1=0; i1<emstsp1.genoListSize(); i1++){
			if(dist1.probs(i1)>0){
				sum+=dist1.probs(i1);
				int cn = emstsp1.getCN(i1);
				int b = emstsp1.getBCount(i1);
				int i = emstsp.getByAlias(cn, b);
				this.setProbs(i, dist1.probs(i1));
			}
		}
		if(sum<1) {
			throw new RuntimeException("sum less than 1 "+sum);
		}
		
	}
	public void standardise(double mean, double sd) {}
	public String getProbString(List<String> nme) {
		  int[] order = this.getOrder();
		  StringBuffer sb = new StringBuffer();
		  for(int i=0; i<order.length; i++){
			  int k = order[i];
			  double p = this.probs(k);
			  if(p>1e-3){
			    sb.append(nme.get(k)+"="+String.format("%5.3g",p)+";");	  
			  }
		  }
		return sb.toString();
	}
	
		
    
}
