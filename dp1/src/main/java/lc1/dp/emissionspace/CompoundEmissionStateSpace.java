package lc1.dp.emissionspace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.stats.IntegerDistribution;
import lc1.stats.MixtureDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import lc1.util.CopyEnumerator;

import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;


public class CompoundEmissionStateSpace extends EmissionStateSpace{
     EmissionStateSpace[] members;
     static class ComparaNew implements Comparator<int[]>, Serializable{
    	 public int compare(int[] o1, int[] o2) {
    	    	if(o1.length!=o2.length) return o1.length<o2.length ? -1 :1;
    	      for(int i=0; i<o1.length; i++){
    	          if(o1[i]!=o2[i]) return o1[i]<o2[i] ? -1 :1;
    	      }
    	      return 0;
    	    }
     };
  public static final  Comparator<int[]> comp1 = new ComparaNew();

   
       

   protected int[][] haploToMember;  //maps haplotype list to its member haplotypes
   public SortedMap<int[], Integer> membersIndexToIndex = new TreeMap<int[] , Integer>(comp1);
   public EmissionStateSpace[] getMembers(){
       return members;
   }

  // double[] hwe_dist = null;
   
   Map<Integer,PseudoDistribution> zeroDistProb = null;
   
   public PseudoDistribution getZeroDistProb(int coverage, int ploidy) {
	   if(zeroDistProb==null){
		   zeroDistProb = new HashMap<Integer, PseudoDistribution>();
		   
		double[] dist = new double[defaultList.size()];
		BinomialDistribution bin = new BinomialDistributionImpl(ploidy,0);
		PoissonDistribution poiss = new PoissonDistributionImpl(coverage);
		double[] d = new double[2*coverage];
		//double[] p1 = new double[2*coverage];
		double sum =0;
		Arrays.fill(d,0);
		for(int i=0; i<d.length; i++){
			d[i] = poiss.probability(i);
			sum+=d[i];
		}
		for(int i=0; i<d.length; i++){
			d[i] = d[i]/sum;
		}
		sum=0;
		for(int i=0; i<dist.length; i++){
			if(this.getCN(i)==ploidy){
				//Arrays.fill(p1,0);
				int b = this.getBCount(i);
				bin.setProbabilityOfSuccess((double)b/(double)ploidy);
				double p=0;
				
				for(int k=0; k<d.length; k++){
					double p0=1;
					if(k>0){
						bin.setNumberOfTrials(k);
						p0 = bin.probability(0);
					}
				//	p1[k] = p0;
					p+=d[k] *p0;
				}
				dist[i] = p;
				sum+=p;
			}else dist[i] =0;
		}
		for(int i=0; i<dist.length; i++){
			dist[i] = dist[i]/sum;
		}
		zeroDistProb.put(coverage,  new SimpleExtendedDistribution(dist, Double.POSITIVE_INFINITY, this));
	   }
	   return zeroDistProb.get(coverage);
	}
   
   private  double[] getHWEDist(Double bfrac, int stage) {
	 
	  //double[] d = //new double[] {0,1};
	//	 (double[]) Constants.modifyFrac(2+stage).clone();
	//  double[] d1 = new double[d.length];
	  EmissionStateSpace[] mems =  getMembers();
	//  EmissionStateSpace mem = mems[0];
	
	   
       double[] hwe_dist = new double[defaultList.size()];
    /*   double[] dist1 = new double[mem.defaultList.size()];
       for(int k=0; k<mem.cnLength(); k++){
    	  int [] geno =  mem.getGenoForCopyNo(k);
    	  double d1 = d[k]/(double)geno.length;
    	  for(int k1=0; k1<geno.length; k1++){
    		 dist1[geno[k1]] =d1;  
    	  }
       }*/
       
//    	   Arrays.fill(dist1, 1.0/(double)dist1.length);
       
    	   
       if(true || Constants.ignoreHWEInData()|| bfrac==null){
    	 Arrays.fill(hwe_dist,1.0/(double)hwe_dist.length);
       }
       else{
    	  // calcExpected(hwe_dist, dist1, 1.0);
       }
      
    	   
       if(bfrac!=null){
    	   for(int i=0; i<hwe_dist.length; i++){
    		   double bfrac_ =(double) getBCount(i)/(double)getCN(i); 
   	   		if(!Double.isNaN(bfrac_) && Math.abs(bfrac_-bfrac )>0.001 ){
   	   			double expb = bfrac*getCN(i);
   	   		   if( false && Math.IEEEremainder(expb, 1.0)<1e-9){
	   				int al = getByAlias(getCN(i), (int) expb );
	   				hwe_dist[al]+=hwe_dist[i];
	   			}
   	   			hwe_dist[i] = 0.0;
   	   			
   	   		}
   	   //	if(getCN(i)==mems.length) hwe_dist[i] = 0;
   		   }
    	   Constants.normalise(hwe_dist);
       }
       else{
    	   Constants.normalise(hwe_dist);
       }
  if(bfrac==null && false) {
    	   for(int i=0; i<hwe_dist.length; i++){
    		   if(getCN(i)!=2) hwe_dist[i] = 0.0;
    		 
    	   }
       }
//  Constants.normalise(hwe_dist);
//       else 
	  return hwe_dist;
     
}
   
	
   
   CompoundEmissionStateSpace(){
	   
   }
   public CompoundEmissionStateSpace(final EmissionStateSpace[] stateSpaces, boolean onlyRepeats){
     //   super( noCop);
        this.members = stateSpaces;
        init(initStateSpace(stateSpaces, onlyRepeats));
      this.typeIndices = new int[stateSpaces.length];
      for(int i=0; i<typeIndices.length; i++){
    	  typeIndices[i] = i;
      }
        haploToMember = new int[haploList.size()][];
    //    haploToTypeIndex = new int[haploList.size()][];
        for(int i=0; i<this.haploList.size(); i++){
            ComparableArray compa = (ComparableArray) backTranslate(this.haploList.get(i));
            haploToMember[i] = new int[compa.size()];
      //      haploToTypeIndex[i] = new int[compa.size()];
            for(int ij=0; ij<compa.size(); ij++){
            	String va = members[ij].getHaploString(compa.get(ij));
            	Integer val = members[ij].getHapl(va);
                haploToMember[i][ij] =  val;
              //  haploToTypeIndex[i][ij] = ij;
            }
        }
        for(int i=0; i<this.haploToMember.length; i++){
            int[] members = haploToMember[i];
            membersIndexToIndex.put(members, i);
        }
        for(int i=0; i<this.defaultList.size(); i++){
            int[] hapL = this.getHaps(i);
            for(int i1=1; i1<hapL.length; i1++){
                if(this.haploList.get(i1).toString().equals(this.haploList.get(0).toString())) throw new RuntimeException("!!");
            }
        }
        SortedSet<Integer> s = new TreeSet<Integer>();
        for(int i=0; i< this.haploPairToHaplo.length; i++){
            s.add(haploPairToHaplo[i].length);
        }
  
    }
 
   protected ComparableArray backTranslate(Comparable comparable) {
       return (ComparableArray) comparable;
    }



    public int[] getMemberIndices(int haploIndex){
        return haploToMember[haploIndex];
    }
    
    public void calcExpected(double[] overall, double[] d_allele, double sum){
    	
    	Arrays.fill(overall,0.0);
		for(int j=0; j<haploSize(); j++){
    		int geno = getGenoForHaplopair(getHaploPairFromHaplo(j));
    		int[] ind = this.getMemberIndices(j);
    		double prob = 1.0;
    		for(int k=0; k<ind.length; k++){
    			prob*=d_allele[ind[k]];
    		}
    		overall[geno]+=prob*sum;
    	}
    	
	}
    
    protected Comparable translate(ComparableArray array) {
        return array;
      }
     List<Comparable> initStateSpace(final List<Comparable>[] stateSpaces, final boolean onlyRepeats){
        final List<Comparable> set = new ArrayList<Comparable>();
        
        
              CopyEnumerator posm = new CopyEnumerator(stateSpaces.length){
                  public Iterator<Comparable> getPossibilities(int depth) {
                    return stateSpaces[depth].iterator();
                  }
                  @Override
                  public void doInner(Comparable[] list) {
                     if( !exclude(list)){ 
                      set.add(
                              translate(new ComparableArray(Arrays.asList(list))));//PairEmissionState.perm.get(res).get(0));
                     }
                  }
                
                  @Override
                  public boolean exclude(Object obj1, Object prev) {
                	return false;
                  }
                  
                @Override
                  public boolean exclude(Comparable[] list) {
                	if(!onlyRepeats)
                      return false;
                	else{
                		for(int i=1; i<list.length; i++){
                			if(!list[i].equals(list[0])) return true;
                		}
                		return false;
                	}
                  }
             };
             posm.run();
          //  List<ComparableArray> ss =  new ArrayList<Comparable>(set);
             return set;
    }
    
    
 



	public int getIndex(int[] indices) {
      if(indices.length!=members.length){
    	  throw new RuntimeException("mismatch");
      }
		return this.haploToHaploPair[membersIndexToIndex.get(indices)];
      /*  Comparable[] index = new Comparable[indices.length];
        for(int i=0; i<index.length; i++){
            index[i] = this.members[i].get(indices[i]);
           
        }
        ComparableArray compA = new ComparableArray(Arrays.asList(index));
        
        Integer res = this.get(translate(compA));
        if(res==null){
            throw new RuntimeException("");
        }
       // if(indices.length==2 && indices[0]==2 && indices[1] ==3){
       // }
        return res;*/
    }



    @Override
    public String getGenotypeString(Comparable comp) {
    	String res = null;
     if(comp instanceof ComparableArray) res =  ((ComparableArray)comp).getGenotypeString();
     else if(comp instanceof Integer) res =  Integer.toString((Integer)comp, Constants.radix());
     else res =  ((Emiss)comp).toStringShort();
     return res;
    }
    @Override
    public String getHaploString(Comparable comp) {
    	String res = null;
     if(comp instanceof ComparableArray) res =  ((ComparableArray)comp).getHaplotypeString();
     else if(comp instanceof Integer) res =  Integer.toString((Integer)comp, Constants.radix());
     else res =  ((Emiss)comp).toStringPrint();
     return res;
    }


    @Override
    public String getHaploPairString(Comparable comp) {
    	String res = "";
        if(comp instanceof ComparableArray) res =  ((ComparableArray)comp).getHaploPairString();
        else if(comp instanceof Integer) res =  Integer.toString((Integer)comp, Constants.radix());
        else res =  ((Emiss)comp).toStringPrint();
        return res;
    }


private  int[] typeIndices;
	public int[] getMemberTypeIndices(int key) {
		return typeIndices;
	}

	//IntegerDistribution baseDist;
	Map<Integer, IntegerDistribution> baseDists;
	public IntegerDistribution getBaseDist(int st) {
		if(baseDists==null) baseDists = new HashMap<Integer, IntegerDistribution>();
		IntegerDistribution res= baseDists.get(st);
		if(res==null) baseDists.put(st,res = new IntegerDistribution(this.getByAlias(st, 0), this) );
		return res;
		
		
	}
	
	List<Double>bfrac_inds = new ArrayList<Double>();
	public PseudoDistribution getHWEDist1(final Double bfrac) {
		
			if(hweDist==null){
				Set<Double> bfrac1 = new TreeSet<Double>();
				for(int i=0;i<this.genoListSize(); i++){
					double frac =(double) this.getBCount(i)/(double)this.getCN(i);
					if(!Double.isNaN(frac)){
						bfrac1.add(frac);
					}
				}
				bfrac_inds = new ArrayList<Double>(bfrac1);
				bfrac_inds.add(null);
				hweDist = new PseudoDistribution[bfrac_inds.size()];
			}
			int ind = bfrac_inds.indexOf(bfrac);
			if(hweDist[ind]==null){
				hweDist[ind] = new SimpleExtendedDistr(this.getHWEDist(bfrac,0), Double.POSITIVE_INFINITY, this, bfrac);
				
				hweDist[ind].setDataIndex((short)-2);
			}
			return hweDist[ind];
		}
	
	public void recalcHWEDist(){
		for(int i=0; i<bfrac_inds.size(); i++){
			if(hweDist[i]!=null){
				double[] probs = this.getHWEDist(bfrac_inds.get(i),1);
				((SimpleExtendedDistr)hweDist[i]).update(probs);
			}
		
		}
	}
	public class SimpleExtendedDistr extends SimpleExtendedDistribution{
		Double bfrac;
		public SimpleExtendedDistr(double[] dist, double u,
				CompoundEmissionStateSpace stsp, Double bfrac) {
			super(dist, u, stsp);
			this.bfrac = bfrac;
			// TODO Auto-generated constructor stub
		}
		public void update(double[] probs1) {
		 System.arraycopy(probs1, 0, this.probs, 0, probs1.length);
			
		}
		@Override
		public double weight(){
			return 0;
		}
		@Override
		public PseudoDistribution swtchAlleles() {
			if(bfrac==null) return this;
	    	return CompoundEmissionStateSpace.this.getHWEDist1(1-bfrac);
	    }
	
		
	}
	
		
	private PseudoDistribution[] hweDist;
	private PseudoDistribution[] softenedDist;
	
	
	public PseudoDistribution getSoftenedDist(int ind, double prob) {
		double bfrac = (double) this.getBCount(ind)/(double) this.getCN(ind);
		//int index = probeOnly ? 0 :1;
		if(prob==0) return this.getHWEDist1(bfrac);
			//PseudoDistribution dist = 
		 PseudoDistribution softenedDist = new MixtureDistribution(
				  new double[] {prob, 1-prob},
				 new PseudoDistribution [] {getIntDist(ind), 
				this.getHWEDist1(null)});
		return softenedDist;
	}

	

	
	
	
	



	


	



  /* private Integer translateToTarget(Integer comp, CompoundEmissionStateSpace target, int memInd) {
     //   EmissionStateSpace[] emStSp_0 = emStSp.getMembers();
        //int memInd = 0;
        if(target.size()==size()) return comp;
        int[] memberIndices =getMemberIndices(comp);
        int[] targetMemberIndices = new  int[memberIndices.length];
//        EmissionStateSpace[] targetSp = new EmissionStateSpace[memberIndices.length];
        for(int ij=0; ij<memberIndices.length; ij++){
            targetMemberIndices[ij] = ((CompoundEmissionStateSpace)members[ij]).getMemberIndices(memberIndices[ij])[memInd];
            EmissionStateSpace targ1 = ((CompoundEmissionStateSpace)members[ij]).getMembers()[memInd];
            EmissionStateSpace targ2 =  target.getMembers()[memInd];
            if( targ1!=targ2) throw new RuntimeException("!! "+targ1+" "+targ2);
        }   
        return target.getGenotype(target.get(target.getIndex(targetMemberIndices)));
    }*/



	//int[][] haploToTypeIndex;
 


/*double [][]calcArray ;
public double[] getCalcArray(int length) {
    return calcArray[length-2];
    // TODO Auto-generated method stub
  //  return null;
}*/




    
    
    
    
}
