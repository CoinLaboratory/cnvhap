package lc1.dp.states;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import lc1.CGH.Aberation;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.Info;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.data.representation.CSOData;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.CompoundDistribution;
import lc1.stats.DepthDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.IntegerDistribution;
import lc1.stats.MixtureDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import cern.jet.math.Arithmetic;

public class HaplotypeEmissionState extends EmissionState {
   // IlluminaRArray probR;

	public void addAllCounts(HaplotypeEmissionState st){
		PseudoDistribution[] dist1 = st.emissions;
		for(int i=0; i<this.emissions.length; i++){
			EmissionStateSpace emstsp = emissions[i].getEmissionStateSpace();
			for(int j=0; j<emstsp.size(); j++){
				emissions[i].addCount(j, dist1[1].probs(j));
			}
		}
	}
	
	
	
	public void adjustDepth(){
		if(emissions[0] instanceof DepthDistribution){
		double mediand =this.getMedianDepth();
		int k_start = 0;
		int k_end = emissions.length;
		/*	double sumL =0;
		double sumR = 0;
		for(int k=0; k<k_start; k++){
			sumL+= ((DepthDistribution)emissions[k]).depth.doubleValue();
		}
		for(int k=k_start; k<k_start+len; k++){
			sumR+= ((DepthDistribution)emissions[k]).depth.doubleValue();
		}*/
		for(int k=k_start; k<k_end; k++){
			double bground = mediand;//(sumR+sumL)/(2.0*(double)len);
			if(mediand<0.001) throw new RuntimeException("median depth is wrong");
			((DepthDistribution)emissions[k]).average_coverage_per_chromosome = bground/2.0;
		/*	double curr = ((DepthDistribution)emissions[k]).depth.doubleValue();
			sumL+= curr- ((DepthDistribution)emissions[k-len]).depth.doubleValue();
			sumR+=((DepthDistribution)emissions[k+len]).depth.doubleValue() - curr;*/
		}
		}
	}
	
	public List<PseudoDistribution> subSample(int freq, int start, int end){
		int newlen = (int)Math.ceil((double) (end-start)/ (double) freq);
		
		PseudoDistribution[] dist = new PseudoDistribution[newlen];
		for(int k=0; k<newlen; k++){
			int i=start+k*freq;
			dist[k] = this.emissions[i].clone();
			double depth=0;
			double cnt=0;
			for(int j=i; j<start+(k+1)*freq && j<this.noSnps; j++){
				Number dep = ((IlluminaRDistribution)emissions[j]).depth();
				if(dep!=null){
					depth +=dep.doubleValue();
					cnt++;
				}
			}
			if(cnt>0){
				(((IlluminaRDistribution)dist[k])).setDepth(  depth/cnt);
			}
			//else{
				//System.err.println("h");
			//}
		}
		return Arrays.asList(dist);
		
	}
	public void subSample(int[] freq, int[] start, int[] end){
		List<PseudoDistribution> l = new ArrayList<PseudoDistribution>();
		for(int k=0; k<freq.length; k++){
			l.addAll(subSample(freq[k],start[k],end[k]));
		}
		
		this.noSnps = l.size();
		this.emissions = l.toArray(new PseudoDistribution[0]);
	}
	
	public double getAvgDepth(){
		double d =0;
		double cnt=0;
		//List<Integer> d = new ArrayList<Integer>();
		for(int k=0; k<this.emissions.length; k++){
			Number dept = ((IlluminaRDistribution)emissions[k]).depth();
			if(dept!=null){
			d+=dept.doubleValue();
			cnt++;
			}
		}
		return d / cnt;
	}
	public double getMedianDepth(){
		List<Double> d = new ArrayList<Double>();
		for(int k=0; k<this.emissions.length; k++){
			Number dept = ((IlluminaRDistribution)emissions[k]).depth();
			if(dept!=null) d.add(dept.doubleValue());
		}
		Collections.sort(d);
		double len = d.size();
		double rem = Math.IEEEremainder(len, 2);
		if(rem==0) return d.get((int) (len/2.0)).doubleValue();
		else{
			double d1 = d.get((int) ((len+1.0)/2.0)).doubleValue();
			double d2 = d.get((int) ((len-1.0)/2.0)).doubleValue();
			return (d1+d2)/2.0;
		}
	}
	
   Integer noCop= null;
   @Override
   public void setIndex(int i){
	    super.setIndex(i);
	    this.data_index = (short)i;
	    for(int ik=0; ik<emissions.length; ik++){
	    	if(emissions[ik]!=null) emissions[ik].setDataIndex((short) i);
	    }
	   }
   public void standardiseVariance(double stderrmult){
	   for(int i=0; i<this.emissions.length; i++){
		   PseudoDistribution dist = emissions[i];
		   if(dist instanceof IlluminaRDistribution){
			  dist.applyStdErrCorrection(stderrmult);
		   }
	   }
   }
   public PseudoDistribution emissions(int i) {
		return this.emissions[i];
	}
   public Comparable getElement(int i) {
       return this.emStSp.get(this.getBestIndex(i));
   }
   /* (non-Javadoc)
    * @see lc1.dp.data.representation.PIGData#getDeletedPositions(lc1.dp.states.EmissionState)
    */
   public Collection<Aberation> getDeletedPositions(List<Integer> geno, List<Integer> loc,  PrintWriter pw, int cn, double thresh, double thresh1) {
   	  List<Aberation> l = new ArrayList<Aberation>();
   	  Aberation ab=null;
   	  double prod=1;
   	  double count = 0;
   	double[] scaleLoc = Constants.scaleLoc();
       for(int i=0; i<length(); i++){
       	       PseudoDistribution dist = this.emissions[i];
			
				double prob = 0;
				double prob1 = 0;
				if(i<length()-1){
					 PseudoDistribution dist1 = this.emissions[i+1];
					 for(int j=0; j<geno.size(); j++){
							prob1+=dist1.probs(geno.get(j));//score(geno[j], i);
					 }
				}
				
					for(int j=0; j<geno.size(); j++){
						prob+=dist.probs(geno.get(j));//score(geno[j], i);
					}
					if(prob>thresh){
						//int pos = loc.get(i);
						String loci = loc.get(i)+"";
						if(scaleLoc!=null){
							double a = ((double)loc.get(i)/scaleLoc[1]);
							double chr = Math.floor(a);
							double pos = (a -chr)*scaleLoc[0];
							loci = chr+":"+pos;
						}
						pw.println(this.getName()+"\t"+loci+"\t"+cn+"\t"+prob);
						if(ab==null & prob1>thresh1){
							ab = new Aberation(this.getName(), i, cn);
							count=1;
							prod = prob;
						}
						else{
							count++;
							prod *=prob; 
						}
					}
					else{
						if(ab!=null && prob1<=(1-thresh1)){
							ab.end = i-1;
							ab.certainty = Math.pow(prod, 1.0/count);
							l.add(ab);
							ab = null;
						}
					}
		}
       if(ab!=null){ab.end = length()-1;
		ab.certainty = Math.pow(prod, 1.0/count);
		
    	   l.add(ab);
       }
       pw.flush();
       return l;
   }
  /*  @Override 
    public Integer noCop(int di){
        return noCop;
        
    }*/
    
    public Integer noCop(){
    	return noCop;
    }
    
    public void setName(String name) {
        this.name = name;
        
    }
    public void mix() {
        for(int i=0; i<this.emissions.length; i++){
           ((IntegerDistribution)emissions[i]).mix(this.getEmissionStateSpace());
        }
          
      }
  
   // final public boolean train_j = true; //whether or not to train
    final EmissionStateSpace emStSp;
    public void reverse(){
        List<PseudoDistribution> l =Arrays.asList(emissions);
        Collections.reverse(l);
        emissions = l.toArray(new PseudoDistribution[0]);
    }
    
    
    public  PseudoDistribution[] emissions; //prob of emitting true
    
    
    public  void append(EmissionState emissionState){
    	int len = emissions.length;
    	PseudoDistribution[] ems = new PseudoDistribution[this.emissions.length+emissionState.length()];
    	System.arraycopy(this.emissions, 0, ems, 0, len);
    	
    	for(int k=0 ; k<emissionState.length(); k++){
    		ems[len+k] = 
    		((HaplotypeEmissionState)emissionState).emissions(k);
    	}
    	this.emissions =ems;
    }
    
    public final ProbabilityDistribution[][] emissionsDatatype;
  //  final boolean[] train;
 //   final double[] pseudo;
    int noSnps;
   
  /* public boolean same(HaplotypeEmissionState st1, double thresh){
       for(int i=0; i<emissions.length; i++){
           double[] d1 = emissions[i].probs();
           double[] d2 = st1.emissions[i].probs();
           for(int j=0;j<d1.length; j++){
               if(Math.abs(d2[j] - d1[j])> thresh) return false; 
           }
       }
       return true;
   }*/
   public void restrictSites(int i) {
       PseudoDistribution[] em1 = new PseudoDistribution[emissions.length];
       System.arraycopy(emissions, 0, em1, 0, (int) Math.min(i, emissions.length));
       this.emissions = em1;
   }
   /* (non-Javadoc)
    * @see lc1.dp.data.representation.SSOData#restrictSites(int, int)
    */
   public void restrictSites(int min, int max) {
       PseudoDistribution[] em1 = new PseudoDistribution[max-min+1];
       System.arraycopy(emissions, min, em1, 0, em1.length);
       this.emissions = em1;
        
    }
   
   @Override
   public String getUnderlyingData(int i) {
       
       return this.emissions[i].getUnderlyingData(this.emStSp);
   }
  
   public void setTheta(double[] prob,  int i){
       emissions[i].setProb(prob);
     
      // System.arraycopy(prob, 0, emissions[i].pseudo, 0, prob.length);
   }
   
   public void setTheta(int val, int i2) {
       emissions[i2] = new IntegerDistribution(val, this.emStSp);
   }
   
   @Override
   public boolean transferCountsToProbs( double pseudoC) {
       if(pseudoC > 1e3 ) return super.transferCountsToProbs(pseudoC);
      // double pseudoC1 = pseudoC/(double) this.emStSp.size();
       paramIndex++;
           for(int i=0; i<this.emissions.length; i++){
               PseudoDistribution em = emissions[i];
               if(em==null) continue;
               if(em instanceof SimpleExtendedDistribution){
            	//  Logger.global.info(i+" before "+em);
            	   if(prior!=null){
            		   ((SimpleExtendedDistribution)em).transfer(this.prior.emissions[i].probs(), pseudoC);
            	   }
            	   else{
            		   em.transfer(pseudoC);
            	   }
                 //  Logger.global.info("state probs "+this.name+" "+i+" after "+em);
                   PseudoDistribution fixed = ((SimpleExtendedDistribution)em).makeFixed();
                   if(fixed!=null){
                       emissions[i] = fixed;
                   }
               }
           }
           if(this.emissionsDatatype!=null){
               for(int k=0; k<this.emissionsDatatype.length; k++){
                   for(int i=0; i<emissions.length; i++){
                       ProbabilityDistribution em = emissionsDatatype[k][i];
                       if(em==null) continue;
                           em.transfer(pseudoC);
                           
                   }
               }
           }
        return true;
    }
   @Override 
   public void switchAlleles(int i){
      if(emissions[i]!=null) emissions[i] = this.emissions[i].swtchAlleles();
   }
  short data_index;
  
  //used for data
   public HaplotypeEmissionState(String name, int noSnps, EmissionStateSpace emStSp, short data_index){
       super(name, 1);
       this.data_index = data_index;
     //  this.distribution = new double[emStSp.defaultList.size()];
      // cn = new int[distribution.length];
      // for(int j=0; j<distribution.length; j++){
       //	cn[j] = emStSp.getCN(j);
      // }
       this.noCop = null;
       this.emStSp = emStSp;
       this.noSnps = noSnps;
       this.emissions = new PseudoDistribution[noSnps];
       this.emissionsDatatype =  null;
     
   }
   
   
   public HaplotypeEmissionState(String name, int noSnps, int len, EmissionStateSpace emStSp, Integer noCop, double[][] meanvarskew){
       this(name, noSnps, len, emStSp, noCop, meanvarskew, (short)-1);
   }
  public HaplotypeEmissionState(String name, int noSnps, int len, EmissionStateSpace emStSp, Integer noCop, 
          double[][] meanvarskew, short data_index){
    this(name, noSnps, emStSp, data_index);
    this.noCop = noCop;
   // this.probR = noCop==null ? null : getProbRGroup(name, noCop,1.0, meanvarskew);
    //int[] numLevels = Constants.numLevels();
      for(int i=0; i<noSnps; i++){
          emissions[i] = new SimpleExtendedDistribution(len,emStSp);//{0.475, 0.475, 0.05}
          emissions[i].setDataIndex(data_index);
         /* for(int k=0; k<this.emissionsDatatype.length; k++){
              emissionsDatatype[k][i] = 
                  numLevels[k] >0 ? 
                          (ProbabilityDistribution)   new SimpleExtendedDistribution1(numLevels[k]) :
                        (ProbabilityDistribution) new TrainableNormal(0,1.0, 1000, 100);
                 
          }*/
      }
  }
   public void setPhenotype(Double[] phen){
       this.phenValue = phen;
   }
  
  // int[] offset;
   
   public HaplotypeEmissionState(String name, int noSnps, double u, 
            double[][] init, EmissionStateSpace emStSp, Integer[] noCop,  ProbabilityDistribution[]  numLevels
       ){
       super(name, 1);
   //    this.offset = offsets;
       this.data_index = -1;
      double[][] init1 = new double[init.length][0];
      int[] index = new int[init.length];
      boolean[] fixed = new boolean[init.length];
      this.noCop = noCop[0];
     // int[] noCops = new int[init.length];
      for(int k=0; k<init1.length; k++){
    	  if(noCop[k]!=noCop[0]) noCop = null;
    	  init1[k] = new double[init[k].length];
      for(int i=0; i<init1[k].length; i++){
    	  if(emStSp.getBCount(i)>0){
    		  init1[k][i] = 0;
    	  }
    	  else{
    		  init1[k][i] = init[k][i];
    	  }
      }
       Constants.normalise(init1[k]);
      
       index[k] = Constants.getMax(init[k]);
       fixed[k] =  init[k][index[k]]>0.9999;
      }
     
   
     
     
     
     //  this.probR = getProbRGroup(name, noCop,1.0, meanvarskew);
   //    this.train = new boolean[noSnps];
   //    Arrays.fill(train, false);
       this.emStSp = emStSp;
       this.noSnps = noSnps;
       BinomialDistribution binom = new BinomialDistributionImpl(init[0].length-1,0.5);
       this.emissions = new PseudoDistribution[noSnps];
       this.emissionsDatatype = new ProbabilityDistribution[Constants.countDT ?numLevels.length : 0][noSnps];
       Info[][] map = null;
       if((DataCollection.datC) instanceof MergedDataCollection){
    	 map =   ((MergedDataCollection) DataCollection.datC).map;
       }
       for(int i=0; i<noSnps; i++){
    	   int di =map==null ? 0 :getIndex(map[i],index);
    	  // System.err.println("i " +i);
    	   if(fixed[di]){
    		   emissions[i] =new IntegerDistribution(index[di], emStSp);
    	   }
    	   else{
    		   Boolean prO =  DataCollection.datC.probeOnly(i);
    		 
//    		   double p = ;
  //  		   double q = 1-p;
    	   double[] init2 ;//=  ? init : init1;
    	   if(prO ==null || !prO) {
    		   binom.setProbabilityOfSuccess(DataCollection.datC.baf(i));
    		   init2 = new double[init[0].length];
    		   for(int k=0; k<init2.length; k++){
    			   init2[k] = binom.probability(k);
    			   
    		   }
    	   }
    	   else{
    		   init2 = init1[di];
    	   }
    	   int index1 = Constants.getMax(init2);
           emissions[i] = 
        	   init2[index1]>0.9999 ? new IntegerDistribution(index1,emStSp) : 
               new SimpleExtendedDistribution1(init2, u, emStSp);//{0.475, 0.475, 0.05}
        	  
    	   }
               for(int k=0; k<this.emissionsDatatype.length; k++){
                   emissionsDatatype[k][i] = numLevels[k].clone(100);
                   
                   emissionsDatatype[k][i].initialise();
               }
        //   emissions[i].transferProbToPseudo();
       }
       if(Constants.CHECK && Double.isNaN(emissions[0].probs(0))){
    	   throw new RuntimeException("!!");
       }
      // this.pseudo = new double[init.length];
     // System.arraycopy(init, 0, pseudo, 0, init.length);
   }
  
  
  private int getIndex(Info[] infos,int[] index) {
	  int ind = -1;
	for(int k=0; k<infos.length; k++){
		if(infos[k]!=null){
			if(ind>=0 && index[ind]!=index[k]) {
				throw new RuntimeException("data types overlapping with diferent cn");
			}
			ind = k;
		}
	}
	return ind;
}

/* public HaplotypeEmissionState(String name, int noSnps, double[] d, EmissionStateSpace emStSp, int noCop, double[][] meanvarskew, double r_prior){
       super(name, 1);
       this.noCop = noCop;
       this.data_index = -1;
       int index = Constants.getMax(d);
       boolean fixed = d[index]>0.999;
       this.probR = getProbRGroup(noCop,1.0, meanvarskew, r_prior);
      // this.train = new boolean[noSnps];
      // Arrays.fill(train, false);
     //  this.nullM = null;
       this.emStSp = emStSp;
       this.noSnps = noSnps;
       this.emissions = new PseudoDistribution[noSnps];
       this.emissionsDatatype = new ProbabilityDistribution[Constants.numLevels().length][noSnps];
       for(int i=0; i<noSnps; i++){
           emissions[i] =
               fixed ? new IntegerDistribution(index) : 
               new SimpleExtendedDistribution1(d, Double.POSITIVE_INFINITY);//{0.475, 0.475, 0.05}
       }
     //  this.pseudo = new double[d.length];
     //  System.arraycopy(d, 0, pseudo, 0, d.length);
   }*/
   
   
   
 private double[] modify(double[] init, EmissionStateSpace emStSp2, double a_only, int[] noA, int[] noB) {
	double[] res = new double[init.length];
	Arrays.fill(res, 0.0);
	if(a_only>0 && a_only <1) return init;
	for(int i=0; i<res.length; i++){
		if(init[i]>0){
			if(noA[i] +noB[i]==0) {
				res[i] = 1.0;
			}
			else{
			double bin = Arithmetic.binomial(noA[i]+noB[i], noA[i]);
			res[i] =bin * Math.pow(a_only, noA[i]) * Math.pow(1-a_only, noB[i]);
			}
		}
	}
	Constants.normalise(res);
	return res;
}
//  public void setFixed(boolean f){
 //      for(int i=0; i<emissions.length; i++){
 //          emissions[i].setFixed(f);
 //      }
 //  }
   public void print(PrintWriter pw, String prefix){
       StringBuffer sb1= new StringBuffer(prefix);
     //  StringBuffer sb2= new StringBuffer(prefix);
       for(int i=0; i<this.noSnps; i++){
           sb1.append("%8.2g ");
        }
       Double[] em = new Double[this.noSnps];
       for(int i=0; i<this.getEmissionStateSpace().size(); i++){
           for(int j=0; j<this.noSnps; j++){
               em[j] = emissions[j]==null ? Double.NaN : this.emissions[j].probs(i);
           }
           pw.println(String.format(sb1.toString(), em));  
       }
      
      
   }
   public void print(PrintWriter pw, String prefix, List<Integer>columns){
      // for(int i=0; i<probR.length(); i++){
      //     pw.println("prob R_"+i+this.probR.get(i).toString());
     //  }
       StringBuffer sb1= new StringBuffer(prefix);
       StringBuffer sb2= new StringBuffer(prefix);
       Object[] em = emissions(this, columns);
       for(int i=0; i<((Double[])em[0]).length; i++){
          sb1.append("%8.2g ");
          sb2.append("%8s ");
       }
       pw.println(String.format(sb2.toString(),(String[]) em[1]));  
       pw.println(String.format(sb1.toString(),(Double[]) em[0]));  
   }
   public HaplotypeEmissionState(HaplotypeEmissionState st_to_init,  String name) {
	   this(st_to_init,name, Double.POSITIVE_INFINITY);
   }
   
   /** u controls the sampling.  A high value results in a closer copy
    * use st_to_init to initialise and st_to_pseudo to set pseudo counts
    *  */
   public HaplotypeEmissionState(HaplotypeEmissionState st_to_init,  String name, double d) {
      super(name, 1);
      this.prior = st_to_init;
      this.data_index = st_to_init.data_index;
   //   this.probR = st_to_init.probR;
      this.noCop = st_to_init.noCop;
   //   this.train = st_to_init.train;
      this.emStSp = st_to_init.emStSp;
    if(st_to_init.phenValue!=null){
         this.phenValue = new Double[st_to_init.phenValue.length];
         System.arraycopy(st_to_init.phenValue, 0, phenValue, 0, phenValue.length);
     }
    //  this.pseudo = st_to_init.pseudo;
    //  this.nullM = st_to_init.nullM;
     // if(st_to_pseudo!=null && !st_to_init.getName().equals(st_to_pseudo.getName())) throw new RuntimeException("!!");
      this.noSnps = st_to_init.noSnps;
      this.emissions =new PseudoDistribution[noSnps];
          this.emissionsDatatype =
        	  st_to_init.emissionsDatatype==null ? null : 
        	  new ProbabilityDistribution[st_to_init.emissionsDatatype.length][noSnps];
       for(int i=0; i<noSnps; i++){
           emissions[i] =  st_to_init.emissions[i]==null ? null : st_to_init.emissions[i].clone(d);
         
          if(emissionsDatatype!=null){
               for(int k=0; k<emissionsDatatype.length; k++){
                   emissionsDatatype[k][i] = st_to_init.emissionsDatatype[k][i].clone(d);
               }
          }
       }
   }
   
   HaplotypeEmissionState prior = null;
  
public HaplotypeEmissionState(EmissionState state_j) {
    super(state_j.getName(), 1);
    this.data_index = ((HaplotypeEmissionState)state_j).data_index;
   this.noCop = state_j.noCop();
 //  this.probR = state_j.probR();
   this.emStSp = state_j.getEmissionStateSpace();
   this.noSnps = state_j.noSnps();
   this.emissions =new PseudoDistribution[noSnps];
   if(((HaplotypeEmissionState)state_j).emissionsDatatype==null){
       emissionsDatatype=null;
   }
   else{
   this.emissionsDatatype = 
       new ProbabilityDistribution[((HaplotypeEmissionState)state_j).emissionsDatatype.length][noSnps];
       for(int k=0; k<emissionsDatatype.length; k++){
           for(int i=0; i<noSnps; i++){
               emissionsDatatype[k][i] =((HaplotypeEmissionState)state_j).emissionsDatatype[k][i].clone();
           }
       }
   }
   if(((HaplotypeEmissionState)state_j).phenValue!=null){
       this.phenValue = new Double[((HaplotypeEmissionState)state_j).phenValue.length];
       System.arraycopy(((HaplotypeEmissionState)state_j).phenValue, 0, phenValue, 0, phenValue.length);
   }
   for(int i=0; i<noSnps; i++){
       Integer fixed_i = state_j.getFixedInteger(i);
       if(fixed_i!=null){
           emissions[i] = new IntegerDistribution(fixed_i,emStSp);
       }
       else{
           double[] probs = new double[emStSp.size()];
           for(int j=0; j<probs.length; j++){
               double sc = state_j.score(j,i);// ****VERY DANGEROUS *****   
               probs[j]=sc;
           }
           emissions[i] =  new SimpleExtendedDistribution(probs, Constants.switchU(), emStSp);
          
       }
       
   }
   
}


public Object clone(State pseudo){
       return new HaplotypeEmissionState(this, this.getName());
   }
   
   public Object clone(){
       return new HaplotypeEmissionState(this,this.getName());
   }
   
   /*[true, false, null] */
    public void addCount(int obj_index, int data_index, double value, int i) {
       // if(train_j){
          emissions[i].addCount(obj_index, value);
          if(prior!=null){
        	  prior.addCount(obj_index, value, i);
          }
         //  emissionsDatatype[i].addCount(data_index, value);
       // }
    }
    
    public void addCount(int obj_index,  double value, int i) {
     //   if(train_j){
    	/*if(i==0){
    		System.err.println("added "+this+" "+obj_index+" "+value);
    	}*/
          emissions[i].addCount(obj_index, value);
          if(prior!=null){
        	  prior.addCount(obj_index, obj_index, value, i);
          }
      //  }
    }
    public void addCountDT(double phen,   int phen_index,  double value, int i) {
        if( emissionsDatatype!=null){
          //  for(int k=0; k<obj_index.length; k++){
                emissionsDatatype[phen_index][i].addCount(phen, value);
          //  }
        }
    }
     Double[] phenValue = new Double[0];
    @Override
    public Double[] phenValue(){
        return phenValue;
    }
  
    public void dataIndices(Set<Short>s){
    	for(int i=0; i<this.emissions.length;i++){
    		emissions[i].getDataIndices(s);
    	
    	}
    }
    
    public int dataIndex(){
        return this.data_index;
    }
    @Override
    public int dataIndex(int i){
    	return this.emissions[i].getDataIndex();
    }
    public void initialiseCounts(){
          for(int i=0; i<this.noSnps; i++){
              if(emissions[i]!=null) this.emissions[i].initialise();
            }
          if(this.emissionsDatatype!=null){
              for(int k=0; k<emissionsDatatype.length; k++){
                  for(int i=0; i<this.emissions.length; i++){
                      this.emissionsDatatype[k][i].initialise();
                  }
              }
          }
        }

    public double KLDistance(EmissionState st ){
        double sum=0;
        HaplotypeEmissionState hes = (HaplotypeEmissionState)st;
        for(int i=0; i<emissions.length; i++){
          if(emissions[i].probs()[0]!=0)
               sum+= emissions[i].probs()[0]*Math.log(emissions[i].probs()[0]/hes.emissions[i].probs()[0]);
        }
        return sum/(double)emissions.length;
    }
    
    public int sample(int i) {
        return (int) emissions[i].sample();
    }
   
    public double score(int obj_i, int i1) {
        //int i= emissions.length==1 ? 0 : i1;
           return emissions[i1].probs(obj_i);
      
    }  
    
  
   /* public void setRandom(double emiss, boolean restart){
        for(int i=0; i<noSnps;i++){
           emissions[i].setRandom(emiss, restart);
        }
    }*/

   
public static Object[] emissions(EmissionState d,  List<Integer> columns){
    int len = d.noSnps();
   EmissionStateSpace emStSp = d.getEmissionStateSpace();
    Double[] res = new Double[columns==null ? len : columns.size()];
    String[] res1 = new String[columns==null ?len : columns.size()];
    for(int i=0; i<res.length; i++){
        int pos = columns==null ? i : columns.get(i);
        if(len==1) pos = 0;
        int k = d.getBestIndex(pos);
            res[i] =  d.score(k, pos);
            if(res[i] < 0.01) res[i] = 0.0;
            else if(res[i]>0.99) res[i] = 1.0;
          //  Emiss[] emiss = Emiss.stateSpace;
            Comparable compa =  emStSp.get(k);
            res1[i] =compa instanceof Emiss ?((Emiss)compa).toStringShort() : ((ComparableArray)compa).toStringShort();
      
    }
    return new Object[] {res,res1};
}
public static Object[] emissions(SimpleExtendedDistribution[] d,  List<Integer> columns){
    int len = d.length;
   //EmissionStateSpace emStSp = d.getEmissionStateSpace();
    Double[] res = new Double[columns==null ? len : columns.size()];
    String[] res1 = new String[columns==null ?len : columns.size()];
    for(int i=0; i<res.length; i++){
        int pos = columns==null ? i : columns.get(i);
        if(len==1) pos = 0;
        double[] probs = d[i].probs;
        int k = Constants.getMax(probs);
            res[i] =  probs[k];
            if(res[i] < 0.01) res[i] = 0.0;
            else if(res[i]>0.99) res[i] = 1.0;
          //  Emiss[] emiss = Emiss.stateSpace;
            res1[i] =  k+"";
      
    }
    return new Object[] {res,res1};
}


public String toString(int i){
    return this.getName();
}

  
   
    public void validate() throws Exception {
        this.lengthDistrib.validate();
        for(int i=0; i<noSnps; i++){
            double sum = emissions[i].sum();
          if(Math.abs(1.0-sum)>0.001){
           throw new RuntimeException("invalid! "+emissions[i]);
          }
        }
        
    }

    public int mostLikely(int pos) {
        return this.emissions[pos].getMax();
     
      //  int k1 = nullM.mostLikely(pos);
      // if(emissions[pos].probs[0]>0.5) return Emiss.B;
      //  else return Emiss.A;
    }
   
   /* public void fix(){
        this.setChangedParams(true);
        for(int i=0; i<this.emissions.length; i++){
           emissions[i].fix();
           this.train_j = false;
        }
    }
   
    public void reinitialise(Double[] maf_cases,double u) {
       for(int i=0; i<this.emissions.length; i++){
           if(u==Double.POSITIVE_INFINITY){
               emissions[i].probs[0] = maf_cases[i];
               emissions[i].probs[1] = 1-maf_cases[i];
           }
           else{
               Dirichlet dir = new Dirichlet(new Double[] {maf_cases[i], 1-maf_cases[i]}, u);
               Double[] d = dir.sample();
               emissions[i].probs[0]  = d[0];
               emissions[i].probs[1] = d[1];
           }
       }
        
    }*/
   
 /*   public int length() {
       return this.emissions.length;
    }*/
    public void print(PrintWriter pw, List<Integer>cols) {
       this.print(pw, "", cols);
        
    }
 /*   public void recombine(HaplotypeEmissionState to, int i) {
        PseudoDistribution[] fromEm = emissions;
        PseudoDistribution[] toEm = to.emissions;
        PseudoDistribution[] fromEm1 = new  PseudoDistribution[fromEm.length];
        int len = fromEm.length;
        System.arraycopy(fromEm, 0, fromEm1, 0, i);
        System.arraycopy(toEm, i, fromEm1, i, len-i);
        System.arraycopy(fromEm, i, toEm, i,len-i );
       emissions = fromEm1;
    }*/

    @Override
    public int noSnps() {
       return this.emissions.length;
    }


    

    @Override
    public EmissionStateSpace getEmissionStateSpace() {
       return emStSp;
    }


    @Override
    public Integer getFixedInteger(int i) {
        return emissions[i].fixedInteger();
    }

    public void setFixedIndex(int i, int k) {
        emissions[i].setFixedIndex( k);
       
    }
    
    public boolean hasIlluminaDist(){
        for(int i=0; i<this.emissions.length; i++){
            if(emissions[i] instanceof IlluminaRDistribution  || 
            		(emissions[i] instanceof CompoundDistribution && ((CompoundDistribution)emissions[i]).hasIlluminDist())) return true;
        }
        return false;
    }
    
    int paramIndex =1;
    @Override
    public int getParamIndex() {
       
      int res =  paramIndex;
     // if(probR!=null){
     // for(int i=0; i<probR.length(); i++){
    	  
    	//  ProbabilityDistribution pr1 = this.probR.get(i);
    	//  if(pr1!=null)
    	//  res+=pr1.getParamIndex();
     // }
     // }
      return res;
    }
   @Override
    public void removeAll(List<Integer> toDrop) {
	   //if(true) throw new RuntimeException("!!");
        List<PseudoDistribution> newDist = new ArrayList<PseudoDistribution>();
        int k=0;
        for(int i=0; i<this.length(); i++){
            if(toDrop.contains(i)) continue;
            else{
                newDist.add( this.emissions[i]);
                k++;
            }
        }
        this.emissions = newDist.toArray(new PseudoDistribution[0]);
        this.noSnps = emissions.length;
   //     if(noSnps<55) throw new RuntimeException("!!");
      }

   /* public void fixLikelihood(boolean X, int i){
        double[] em = ((SimpleExtendedDistribution) this.emissions[i]).probs;
        double max = em[Constants.getMax(em)];
        double[] probs = emissions[i].probs();
        for(int k=0; k<probs.length; k++){
        //    if(!X){
         //       if(this.emStSp.getCN(k)!=2) probs[k] =0;
        //    }
         //   else{
                if(this.emStSp.getCN(k)!=((CompoundEmissionStateSpace)this.getEmissionStateSpace()).noCopies()) probs[k] =0;
          //  }
        }
        SimpleExtendedDistribution.normalise(probs);
    }*/
   

  

   
@Override
public void setAsMissing(List<Integer> toDrop, double cn_ratio){
    EmissionStateSpace subSp = this.getEmissionStateSpace();
    List<String> genoList = new ArrayList<String>();
    for(int k=0; k<subSp.getGenotypeList().size(); k++){
        char[] ch =  ((ComparableArray)subSp.getGenotype(k)).toStringShort().replaceAll("_", "").toCharArray();
        Arrays.sort(ch);
         genoList.add(new String(ch));
    }
    for(int i1=0; i1<toDrop.size(); i1++){
        int i = toDrop.get(i1);
        if(this.emissions[i] instanceof IntegerDistribution){
            emissions[i] = new SimpleExtendedDistribution(subSp.defaultList.size());
        }
        double[] prob =((SimpleExtendedDistribution)this.emissions[i]).probs;
        Arrays.fill(prob, 0.0);
        for(int k=0; k<genoList.size(); k++){
            String compa = genoList.get(k);
          
            int genoIndex = emStSp.getFromString(compa);
            int cn= compa.length();
            int[] indices = emStSp.getGenotypeConsistent( genoIndex);
            double[] weights = emStSp.getWeights(genoIndex);
            double val = cn==2 ?cn_ratio :1;
                    //10 : (cn==1 ? 1 : (cn==0 ? 0.1 : 0.01));
            for(int k1=0; k1<indices.length; k1++){
                 prob[indices[k1]] = weights[k1]*val;
            }
        }
        SimpleExtendedDistribution.normalise(prob);
        setTheta(prob, i);
    }
}
    @Override
public void applyAlias(int[] alias) {
    DataCollection.reorder( alias, this.emissions);
     
 }
   
 
   
    
    
  
   
    
    
    /* returns probability */
    public double calcDistribution(EmissionState hmm_state, int i,  double p, boolean log){
        double sum=0;
        double[] distribution = hmm_state.distribution;
        Arrays.fill(distribution,log ? Double.NEGATIVE_INFINITY: 0.0);
        HaplotypeEmissionState data_state = this;
    //    int di = data_state.dataIndex(i);
        Integer nc = hmm_state.noCop();
      
        EmissionStateSpace emStSp = hmm_state.getEmissionStateSpace();
     
        if(distribution.length!=emStSp.defaultList.size()) {
        	throw new RuntimeException("!!");
        }
        PseudoDistribution[] emissions = data_state.emissions;
        //Integer bg = DataCollection.datC.getBGCount(emissions[i].getDataIndex(), i);
      
	        for(int j=0; j<distribution.length; j++){
	        //	Integer nc = hmm_state.noCop(di);
	        	
	            if(nc==null || emStSp.getCN(j)==nc.intValue()){
	          
	                double prob_j =hmm_state.score(j,i)*p;
	                if(prob_j>Constants.countThresh3()){
	                 // 	int j1 = hmm_state.mod(j,di);
	                	double sc1 = emissions[i].scoreBR(emStSp,j, i);
	                  double v = log ? Math.log(prob_j*emStSp.getWeight(j)) + sc1 : prob_j*sc1*emStSp.getWeight(j);
	               
	                    /*(Constants.joint ?
	                    		sc1 : 
	                        emissions[i].scoreR(DataCollection.datC.getBGCount(emissions[i].getDataIndex(), i), emStSp.getCN(j),i)
	                       *emissions[i].scoreB(j,i));*/
	                       //*;  //do we include the weight term???;
	                  distribution[j] =v;
	                   sum+= v;
	                }
	            }
	        }
	        if(log) {
	        	sum=0;
	        	int maxind =  Constants.getMax(distribution);
	        	double v = distribution[maxind];
	        	 for(int j=0; j<hmm_state.distribution.length; j++){
	        		distribution[j]=Math.exp(distribution[j]-v);
	        		sum+=distribution[j];
	        	 }
	         }
      
        return sum;
    }
    public double calcDistribution(EmissionState hmm_state, int i,  double p, boolean log, int mixComp){
        double sum=0;
        double[] distribution = hmm_state.distribution;
        Arrays.fill(distribution, log ? Double.NEGATIVE_INFINITY: 0.0);
       // int di = this.dataIndex(i);
    	//Integer nc = hmm_state.noCop(di);
    	Integer nc = hmm_state.noCop();
        HaplotypeEmissionState data_state = this;
        EmissionStateSpace emStSp = hmm_state.getEmissionStateSpace();
        if(distribution.length!=emStSp.defaultList.size()) {
        	throw new RuntimeException("!!");
        }
        PseudoDistribution[] emissions = data_state.emissions;
      
	        for(int j=0; j<distribution.length; j++){
	        
	            if(nc==null || emStSp.getCN(j)==nc.intValue()){
	                double prob_j =hmm_state.score(j,i)*p;
	               
	                if(prob_j>Constants.countThresh3()){
	                 // 	int j1 = hmm_state.mod(j,di);
	                	double sc1 =    (Constants.joint ?
	                    		emissions[i].scoreBR(emStSp,j, i,mixComp) : 
	    	                        emissions[i].scoreR(DataCollection.datC.getBGCount(emissions[i].getDataIndex(), i), emStSp.getCN(j),i)
	    	                       *emissions[i].scoreB(j,i));
		                double v = log ? Math.log(prob_j*emStSp.getWeight(j)) + sc1 : prob_j*sc1*emStSp.getWeight(j);
	                   //do we include the weight term???;
	                  distribution[j] +=v;
	                   sum+=v;
	                }
	            }
	        }
	        if(log) {
	        	sum=0;
	        	int maxind =  Constants.getMax(distribution);
	        	double v = distribution[maxind];
	        	 for(int j=0; j<hmm_state.distribution.length; j++){
	        		distribution[j]=Math.exp(distribution[j]-v);
	        		sum+=distribution[j];
	        	 }
	         }
      
        return sum;
    }
    public  double calcDistribution( int i, double[] distribution, EmissionStateSpace emStSp){
    	if(true){
    		throw new RuntimeException("!!");
    	}
        double sum=0;
        HaplotypeEmissionState data_state = this;
       // EmissionStateSpace emStSp = data_state.emStSp;
       // IlluminaProbB[] probBState = hmm_state.probB();
      //  IlluminaRArray probRState =  hmm_state.probR();
        PseudoDistribution[] emissions = data_state.emissions;
       // int bg = Constants.backgroundCount();//DataCollection.datC.getBGCount(emissions[i].getDataIndex(), i);
        if(true){
	        for(int j=0; j<distribution.length; j++){
	        //	int j1 = hmm_state.mod(j,di);
	            //if(this.emStSp.getCN(j)==hmm_state.noCop().intValue()){
	               // double prob_j =hmm_state.score(j,i)*p;
	               // if(prob_j>0){
	                  double v = 
	                	 // (Constants.joint ? 
	                			  emissions[i].scoreBR(emStSp,j,i) 
	                			 
	                   //   emissions[i].scoreR(bg, emStSp.getCN(j),i)  * emissions[i].scoreB(j,i))
	                    
	                       *emStSp.getWeight(j);  //do we include the weight term???;
	                  distribution[j] +=v;
	                   sum+=v;
	                //}
	            //}�
	        }
	        return sum;
        }
        else{
        	PseudoDistribution dist= DataCollection.datC.getBGDist(emissions[i].getDataIndex(), i);
        	  for(int j=0; j<distribution.length; j++){
  	            //if(this.emStSp.getCN(j)==hmm_state.noCop().intValue()){
  	               // double prob_j =hmm_state.score(j,i)*p;
  	               // if(prob_j>0){
        		  for(int k=0; k<emStSp.cnLength(); k++){
  	                  double v = 
  	                  //  prob_j 
  	                       emissions[i].scoreR(k, emStSp.getCN(j),i)
  	                      *emissions[i].scoreB(j,i)
  	                       *dist.probability(k)
  	                       *emStSp.getWeight(j);  //do we include the weight term???;
  	                  distribution[j] +=v;
  	                   sum+=v;
        		  }
  	                //}
  	            //}�
  	        }
  	        return sum;
        }
    	
    }
    /** distribution is distribution over cn states */
    public synchronized  double calcDistribution( int bg, int i, double[] distribution, double[] fg_cn){
    	if(true){
    		throw new RuntimeException("!!");
    	}
        double sum=0;
        HaplotypeEmissionState data_state = this;
        EmissionStateSpace emStSp = data_state.emStSp;
        PseudoDistribution[] emissions = data_state.emissions;
    //	Arrays.fill(distribution,0.0);
       // IlluminaProbB[] probBState = hmm_state.probB();
      //  IlluminaRArray probRState =  hmm_state.probR();
      //  int bg = DataCollection.datC.getBGCount(emissions[i].getDataIndex(), i);
        for(int j=0; j<distribution.length; j++){
            //if(this.emStSp.getCN(j)==hmm_state.noCop().intValue()){
             
               // if(prob_j>0){
        	 distribution[j] = 
                  //  prob_j 
                       emissions[i].scoreR(bg, j,i)
                       *fg_cn[j]
                     //  *this.emissions[i].scoreB(j,i)
                       *emStSp.getWeight(j);  //do we include the weight term???;
                sum+=distribution[j];
                //}
            //}
        }
        return sum;
    }
    
    
    
    

    @Override
    public double scoreEmiss(Double[] object_index, int i1){
        double sc = 1.0;
        for(int k=0; k<object_index.length; k++){
            if(object_index[k]==null) continue;
            sc*=emissionsDatatype[k][i1].probability(object_index[k]);
        }
        //int i= emissions.length==1 ? 0 : i1;
       return sc;
        //return probs[object_index];
    }
    
    /*public static PrintWriter scores;
    {
        try{
        scores= new PrintWriter(new BufferedWriter(new FileWriter("scores")));
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }*/
   
    
   
    
    /** warning - for speed we recycle the same d */
    public double[] getEmiss(int i) {
        return emissions[i].probs();
       
    }
    public Set<Short> getDataIndices(Set<Short>res) {
      //List<Short> res = new ArrayList<Short>();
      for(int i=0; i<emissions.length; i++){
         emissions[i].getDataIndices(res);
      }
      return res;
    }
    public void applyCorrection(double parseDouble) {
        for(int i=0; i<this.emissions.length; i++){
            if(emissions[i]!=null) emissions[i].applyCorrection(parseDouble);
        }
        
    }
    public void applyLoess(Double[] loess) {
      for(int i=0; i<emissions.length; i++){
          if(emissions[i]!=null) emissions[i].applyCorrection(-1*(loess[i]));
      }
        
    }
 
    
    public final boolean updateEmissions(int i, PseudoDistribution pseudoDistribution) {
    
      if(emissions[i]==null){
    	
          emissions[i] = pseudoDistribution;
         
      }
      else if (emissions[i] instanceof CompoundDistribution){
          ((CompoundDistribution)emissions[i]).addDist(pseudoDistribution, emStSp);
      }
      else if(emissions[i].equals(pseudoDistribution)){
    	//  System.err.println("distribution is the same");
      }
     /* else if(emissions[i] instanceof IntegerDistribution){
          System.err.println(emissions[i].fixedInteger()+" vs "+pseudoDistribution.fixedInteger());
      }*/
      else{//if(!(emissions[i] instanceof IntegerDistribution)){

    	 
    		  if(pseudoDistribution.getDataIndex()>=-1){
    			  if(emissions[i].getDataIndex()<-1) emissions[i] = pseudoDistribution;
    			  else{
    			   emissions[i] = new CompoundDistribution(this.emissions[i], pseudoDistribution, this.emStSp);
    			  }
    		  }
       
          
      }
      return true;
        
    }
    public void setDataIndex(int dat_ind) {
    	if(dat_ind <0) throw new RuntimeException("!!");
       if(this.data_index<0) this.data_index =(short) dat_ind;
        
    }
    public String getCompressedDataString(int i) {
        return this.emissions[i].getCompressedDataString(this.emStSp);
    }
    public String getCompressedStringHeader() {
    	int min = Math.min(100,emissions.length);
    	String st = "";
    	for(int k=0; k<emissions.length; k++){
    		String str = emissions[0].compressedStringHeader(this.getEmissionStateSpace());
    		if(str.length()>st.length()) st = str;
    	}
    	return st;
    }
      
    
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.CSOData#sampleGenotype(lc1.dp.states.HaplotypeEmissionState, java.util.List)
     */
        public void sampleGenotype( List<CSOData> spList){
          EmissionStateSpace emStSp = getEmissionStateSpace();
          HaplotypeEmissionState emst = this;
          emst.initialiseCounts();
          //emst.initialiseCounts();
//            int[] m = new int[emStSp.size()];
            for(int i=0; i<length(); i++){
                short dat_index = emst.data_index;
                //emst.emissions[i] = new SimpleExtendedDistribution(emStSp.defaultList.size());
                //emst.emissions[i].setDataIndex(dat_index);
      //        Arrays.fill(m, 0);    
                for(int j=0; j<spList.size(); j++){
                    PhasedDataState obj = (PhasedDataState) spList.get(j);
                    double[] probs = obj.emissions[i].probs();
                    for(int k=0; k<probs.length; k++){
                        this.emissions[i].addCount(k, probs[k]);
                    }
                   /* Integer compa_i = emStSp.getHaploPair(obj.getElement(i));
                    if(compa_i==null){
                        try{
                        throw new RuntimeException("no element "+obj.getElement(i)+" in state space \n"+emStSp.defaultList);
                        }catch(Exception exc){
                            exc.printStackTrace();
                        }
                    }
                    else{
                        emst.addCount(compa_i, obj.dataIndex(i), 1.0, i);
                    }*/
                        //m[compa_i]++;
                }
               
               
            }
            emst.transferCountsToProbs(0.0);
          
        }
    	/*public static double[] fillCN( EmissionStateSpace emstsp, double[] d,  int fixed) {
    		if(true) throw new RuntimeException("!!");
    		double[] res = new double[emstsp.size()];
    		double[] probs = new double[emstsp.size()];
    		Arrays.fill(probs, 0.0);
    		ComparableArray compa = ((ComparableArray)emstsp.get(fixed));
    		double fracB = compa.noB()/ (compa.noCopies(true));
    		//if(!compa.noCopies(expandEmiss))
    		for(int k=0; k<probs.length; k++){
    			ComparableArray compa1 = ((ComparableArray)emstsp.get(k));
    			if(compa1.noCopies(true)!=0){
	    			double fracB1 = compa1.noB()/ (compa1.noCopies(true));
	    			probs[k] = 
	    				
	    				Math.abs(fracB - fracB1) < 0.001 ? 1.0 : 0.0;
    			}
    		}
	    	Arrays.fill(res,0.0);
    		for(int cn=0; cn<d.length; cn++){
				
				int cn1 = emstsp.getCN(fixed);
					if(cn1==cn){
						res[i] = d[cn];
					}
					else{
						int[] inds = emstsp.getGenoForCopyNo(cn);
					
						
						double p = d[cn]/(double) inds.length;
						for(int j=0; j<inds.length; j++){
							res[inds[j]] = p;
						}
					}
			}	
			return res;
    	}*/
        /* for background state only */
    	public static double[] fillCN( EmissionStateSpace emstsp, double[] d) {
    		double[] res = new double[emstsp.size()];
	    	Arrays.fill(res,0.0);
    		for(int cn=0; cn<d.length; cn++){
    			
    			
	    			
	    			int[] inds = emstsp.getGenoForCopyNo(cn);
	    			
	    			
	    				//double p = d[cn]/(double) inds.length;
						inner: for(int j=0; j<inds.length; j++){
							if(emstsp.getBCount(inds[j])==0){
								res[inds[j]]=d[cn];
								break inner;
							}
//							res[inds[j]] = p;
						}
	    			
    			}
    	//	double sum = Constants.sum(res);
    		return res;
    	}
        
    	public static double[] fillCN( EmissionStateSpace emstsp, double[] d, double[] probs) {
    		double[] res = new double[emstsp.size()];
	    	Arrays.fill(res,0.0);
    		for(int cn=0; cn<d.length; cn++){
    			
    			
	    			
	    			int[] inds = emstsp.getGenoForCopyNo(cn);
	    			double sum=0;
	    			for(int j=0; j<inds.length; j++){
	    				sum+=probs[inds[j]];
	    			}
	    			if(sum>0){
		    			for(int j=0; j<inds.length; j++){
		    				res[inds[j]] = d[cn]*(probs[inds[j]]/sum);
		    			}
	    			}
	    			else{
	    				double p = d[cn]/(double) inds.length;
						for(int j=0; j<inds.length; j++){
							res[inds[j]] = p;
						}
	    			}
    			}
    	//	double sum = Constants.sum(res);
    		return res;
    	}
    	public void fill( boolean probeOnly, double[] val){
    		CompoundEmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(this.noCop);
    		for(int i=0; i<emissions.length; i++){
    			if(emissions[i]!=null){
    			if(val[1]>0)
    			if(val[1]>=1) emissions[i] = emstsp.getHWEDist1(probeOnly ? new Double(0) : null);
    			else if(val[1]>0){
    				
    				this.emissions[i] = 
    				new MixtureDistribution(
    					val,
    						new PseudoDistribution[] {this.emissions[i],emstsp.getHWEDist1(probeOnly ? new Double(0) : null)});
    			}
    			}
    		}
    	}
    	
    	public void fill( boolean probeOnly, double[] val, List<Integer> sites){
    		CompoundEmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(this.noCop);
    		for(int i1=0; i1<sites.size(); i1++){
    			int i = sites.get(i1);
    			if(emissions[i]!=null){
    			if(val[1]>0)
    			if(val[1]>=1) emissions[i] = emstsp.getHWEDist1(probeOnly ? new Double(0) : null);
    			else if(val[1]>0){
    				
    				this.emissions[i] = 
    				new MixtureDistribution(
    					val,
    						new PseudoDistribution[] {this.emissions[i],emstsp.getHWEDist1(probeOnly ? new Double(0) : null)});
    			}
    			}
    		}
    	}
    	/** keep balance with cn state */
		public  double[] fillCN(int i, EmissionStateSpace emstsp, double[] d) {
		
		   
		    
		    	Integer fixed = this.getFixedInteger(i);
		    	if(fixed==null){
		    		return fillCN( emstsp, d, this.getEmiss(i));
		    	}
		    			else{
		    				ComparableArray compa = (ComparableArray)emstsp.get(i);
		    	    		int copies = compa.noCopies(true);
		    	    		double b = compa.noB();
		    				return fillCN(emstsp, d,getDistWithFrac(
		    						copies ==0  ? 0 :
		    						b / (double)copies, emstsp));
		    					
		    			}
		    		
		    		
		    	
		    	
		   
		}
		public static double[] getDistWithFrac(double fracB, EmissionStateSpace stsp){
	    	double[] res = new double[stsp.defaultList.size()];
	    	for(int i=0; i<res.length; i++){
	    		ComparableArray compa = (ComparableArray)stsp.get(i);
	    		int copies = compa.noCopies(true);
	    		double b = compa.noB();
	    		if(copies ==0 || Math.abs(fracB - b / (double)copies)<0.001){
	    			res[i] = 1.0;
	    		}
	    	}
	    	Constants.normalise(res);
	    	return res;
	    }
		public void setNoCop(int parseInt) {
			this.noCop = parseInt;
			
		}
		public void append(HaplotypeEmissionState nxt1) {
			this.noSnps+=nxt1.noSnps;
			this.emissions = (PseudoDistribution[]) Constants.append(this.emissions, nxt1.emissions,
					
					PseudoDistribution.class);
		}

		public void modify(HaplotypeEmissionState state) {
			this.name = state.name;
		for(int i=0; i <this.length(); i++){
			PseudoDistribution ems1 = state.emissions[i];
			
			if(state.emissions[i] instanceof MixtureDistribution){
				ems1 = 	((MixtureDistribution)state.emissions[i]).dist[0];
			}
			EmissionStateSpace emstsp_ = (emissions[i].getEmissionStateSpace());
			CompoundEmissionStateSpace emstsp1_ = ((CompoundEmissionStateSpace)ems1.getEmissionStateSpace());
			
			 if(ems1 instanceof IntegerDistribution){
				int ind1 = ((IntegerDistribution)ems1).fixedInteger();
				int cn = emstsp1_.getCN(ind1);
				int b = emstsp1_.getBCount(ind1);
				int ind = emstsp_.getByAlias(cn, b);
				this.emissions[i] = emstsp_.getIntDist(ind);
			}
			else{
				
				  emissions[i] = new SimpleExtendedDistribution(emstsp_.genoListSize(),emStSp);//{0.475, 0.475, 0.05}
		          emissions[i].setDataIndex(data_index);
		          emissions[i].modify(ems1);
			}
		}
	}
		
		public  void modify(List<HaplotypeEmissionState> next) {
			HaplotypeEmissionState first = this;
			PseudoDistribution[] ems = first.emissions;
			for(int i=0; i<ems.length; i++){
				EmissionStateSpace emstsp = ems[i].getEmissionStateSpace();
				double[] vals = PairEmissionState.pool.getObj(emstsp.size());
				Arrays.fill(vals,0);
				for(int k=0; k<next.size(); k++){
					HaplotypeEmissionState em = next.get(k);
					PseudoDistribution dist = em.emissions[i];
					//if(dist instanceof SimpleExtendedDistribution){
					//	throw new RuntimeException("!!");
					//}
		//			if(dist instanceof MixtureDistribution) dist = ((MixtureDistribution)dist).dist[0];
					for(int j=0; j<vals.length; j++){
						vals[j]+=dist.probs(j);
					}
				}
				for(int j=0;j<vals.length; j++){
					vals[j] = vals[j]/(double)next.size();
				}
				int max_ind = Constants.getMax(vals);
				PseudoDistribution newDist = null;
				if(1-vals[max_ind]<1e-10){
					newDist = new IntegerDistribution(max_ind, emstsp);;
				}else{
					newDist = new SimpleExtendedDistribution1(vals.clone(),Double.POSITIVE_INFINITY,emstsp);
				}
				if(ems[i] instanceof MixtureDistribution){
					((MixtureDistribution)ems[i]).dist[0] = newDist;
				}else{
					ems[i] =newDist;
				}
				PairEmissionState.pool.returnObj(vals);
			}
		}


		public void insert(int[] old_index, int[] new_index,
				PseudoDistribution baseDist) {
			int len = old_index.length+new_index.length;
			PseudoDistribution[] ems_new = new PseudoDistribution[len];
			for(int k=0; k<this.emissions.length; k++){
				ems_new[old_index[k]] = emissions[k];
			}
			for(int k=0; k<new_index.length; k++){
				ems_new[new_index[k]] = baseDist;
			}
			this.emissions =  ems_new;
		}
		
		
    
}
