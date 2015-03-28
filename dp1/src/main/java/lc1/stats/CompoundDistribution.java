package lc1.stats;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import lc1.dp.emissionspace.EmissionStateSpace;

public class CompoundDistribution extends PseudoDistribution {
	@Override
public boolean isMeasured(){
	for(int i=0; i<l.size(); i++){
		if(l.get(i).isMeasured()) return true;
	}
	return false;
}
	
	public Boolean probeOnly() {
		boolean po1 = true;
		for(int i=0; i<l.size(); i++){
			Boolean po = l.get(i).probeOnly();
			if(po==null) continue;
			po1 = po1 && po;
		}
		return po1;
	}
	
	 public  double[] calcDistribution(
	            double[] distribution, EmissionStateSpace emStSp, int pos) {
	      if(true) throw new RuntimeException("!!");// for(int i=0; i<distribution.length; i++){
	    	//   calcDistribution(distribution,emStSp,pos);
	       //}
	       return distribution;
	     }
	
   public  List<PseudoDistribution> l = new ArrayList<PseudoDistribution>(2);
    public void getDataIndices(Set<Short> res) {
		for(int i=0; i<l.size(); i++){
			l.get(i).getDataIndices(res);
		}
		
	}
    public void getB(double[] b, double mult) {
    	double s = 0;
		for(int j=0; j<l.size(); j++){
			l.get(j).getB(b, mult/(double)l.size());
		}
		
	};
   
    public short[] getDataIndices() {
		short[] res = new short[l.size()];
		for(int i=0; i<l.size();i++){
			res[i] = l.get(i).getDataIndex();
		}
		return res;
	}
   
    public  String getCompressedDataString(EmissionStateSpace emStSp){
    	StringBuffer sb = new StringBuffer();
    	for(int i=0; i<l.size(); i++){
      	sb.append(l.get(i).getCompressedDataString(emStSp));
    	}
    	return sb.toString();
       }

    public CompoundDistribution(PseudoDistribution dist1, PseudoDistribution dist2, EmissionStateSpace emstsp){
        this.data_index = -1;
        dist1.check(dist2, emstsp);
     //   if(true) throw new RuntimeException("!!");
     
      
        this.l.add(dist1);
        this.l.add(dist2);
     }
     public CompoundDistribution(CompoundDistribution cd) {
       for(int i=0; i<cd.l.size(); i++){
           this.l.add(cd.l.get(i).clone());
       }
    }
     public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
         throw new RuntimeException("!!");
      }
     public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
         throw new RuntimeException("!!");
     }
    public String getUnderlyingData(EmissionStateSpace emStSp){
         StringBuffer sb = new StringBuffer();
         for(int i=0; i<l.size(); i++){
             sb.append(l.get(i).getUnderlyingData(emStSp)+";");
         }
         return sb.toString();
     }
     public void applyCorrection(double parseDouble) {
         for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
             it.next().applyCorrection(parseDouble);
         }
     }
     public String getPrintString(){
        return "";
     }
   
    
     @Override
     public void addCount(int obj_index, double value) {
         for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
             it.next().addCount(obj_index, value);
         }
         
     }

    
    
     
     
    
    

     

     @Override
     public PseudoDistribution clone() {
        return new CompoundDistribution(this);
     }

     @Override
     public PseudoDistribution clone(double swtch) {
      return clone();
     }

     @Override
     public double[] counts() {
         throw new RuntimeException("!!");
     }

     @Override
     public double logProb() {
         throw new RuntimeException("!!");
     }

     @Override
     public Integer fixedInteger() {
    	 for(int i=0; i<l.size(); i++){
    		 Integer fix = l.get(i).fixedInteger();
    		 if(fix!=null) return fix;
    	 }
    	 return null;
     }

    

     @Override
     public int getMax() {
         throw new RuntimeException("!!");
     }

     @Override
     public void initialise() {
         for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
             it.next().initialise();
         }
         
     }

     @Override
     public double[] probs() {
         throw new RuntimeException("!!");
     }

     @Override
     public double probs(int obj_i) {
         throw new RuntimeException("!!");
     }

     public double sample() {
         throw new RuntimeException("!!");
     }

     @Override
     public void setCounts(int i1, double cnt) {
        throw new RuntimeException("!!");
         
     }

     @Override
     public void setProb(double[] prob) {
         throw new RuntimeException("!!");
         
     }

     @Override
     public void setProbs(int to, double d) {
         throw new RuntimeException("!!");
         
     }

     @Override
     public double sum() {
         throw new RuntimeException("!!");
     }

     @Override
     public PseudoDistribution swtchAlleles() {
    	 CompoundDistribution dist = new CompoundDistribution(l.get(0).swtchAlleles(), l.get(1).swtchAlleles(), null);
         for(int i=2; i<l.size(); i++){
        	 dist.addDist(this.l.get(i).swtchAlleles(), null);
         }
         
    	 return dist;
         
     }

     @Override
     public void transfer(double pseudoC1) {
          for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
              it.next().transfer(pseudoC1);
          }
         
     }

     @Override
     public void validate() {
         for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
             it.next().validate();
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
    public void addDist(PseudoDistribution pseudoDistribution, EmissionStateSpace emstsp) {
            for(int i=0; i<l.size(); i++){
        pseudoDistribution.compareTo(this.l.get(i));
            }
            l.get(0).check(pseudoDistribution, emstsp);
       this.l.add(pseudoDistribution);
        
    }
    public short getDataIndex(){
        return this.l.get(0).getDataIndex();
    }
    public boolean containsIndex(int index) {
    	for(int i=0; i<l.size(); i++){
		 if(l.get(i).data_index==index) return true;
    	}
    	return false;
	}
    public PseudoDistribution getForIndex(short index) {
        for(Iterator<PseudoDistribution> it = l.iterator(); it.hasNext();){
            PseudoDistribution dist = it.next();
            if(dist.getDataIndex()==index) return dist;
        }
        return null;
    }

	public boolean hasIlluminDist() {
		for(int i=0; i<l.size(); i++){
			if(l.get(i) instanceof IlluminaRDistribution ) return true;
		}
		return false;
	}

	@Override
	public Double getIntensity(boolean r) {
		double res = 0;
		for(int i=0; i<l.size(); i++){
		res+= this.l.get(i).getIntensity(r);
		}
		return res/(double)l.size();
	}

	
	public void addBCount(int j, double val, int i) {
		for(int k=0; k<this.l.size(); k++){
			l.get(k).addBCount(j, val, i);
		}
		
	}


	@Override
	public void addRCount(int bg, int no_cop, double val, int i) {
		for(int k=0; k<this.l.size(); k++){
			l.get(k).addRCount(bg, no_cop, val, i);
		}
	}

	@Override
	public double scoreR(int bg, int no_cop, int i) {
		double prob = l.get(0).scoreR(bg,no_cop, i);
		for(int k=1; k<this.l.size(); k++){
			prob = prob *l.get(k).scoreR(bg,no_cop, i);
		}
		return prob;
	}
	@Override
	public double scoreB(int j, int i) {
		double prob = l.get(0).scoreB(j, i);
		for(int k=1; k<this.l.size(); k++){
			prob = prob *l.get(k).scoreB(j, i);
		}
		return prob;
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emstsp,int j,  int i) {
		double prob = l.get(0).scoreBR(emstsp,j, i);
		for(int k=1; k<this.l.size(); k++){
			prob = prob *l.get(k).scoreBR(emstsp,j, i);
		}
		return prob;
	}

	public int length() {
		return l.size();
	}

	public short getDataIndex(int k) {
		return l.get(k).getDataIndex();
	}


	

	public void setMinMax(double min, double max){
		throw new RuntimeException("!!");
	}

	

	
    
   
}