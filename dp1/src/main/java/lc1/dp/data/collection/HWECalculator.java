package lc1.dp.data.collection;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.ChiSq;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;

public class HWECalculator {
	public DataCollection dc1;
	//int case_index, int_control_index;
	
	    
	public ChiSq chisq = new ChiSq();
	   // public ArmitageTrendTest att;
	    public double[][] []numCases; //[pos][pheno_cat][copy _number]
	//    public Double[] res;
	  //  public Double[] odds;
	//    public Double[] sum;
	////	private final Double[] odds1;
	//	private Double[] sum1;
	    ////public Double[] cases;
	   // public Double[] controls;
	   // public String[] formatOdds;
	    
	    
	   
	   // private static final double log2 = Math.log(2);
	    //public static double maxOdds = 20;
	   
	  
	  //public int type_len;
	  
	 
	 
	    public static void printResults(File f, DataCollection dc
	    		
		  ){
	  try{
		f.mkdir();
		
		  if(! Constants.calcHWE) return;
			  HWECalculator[] ac_ = (dc).getHWECalculator();
			 if(ac_==null) return; 
			 for(int kk=0; kk<ac_.length; kk++){
				 if(ac_[kk]==null) continue;
				 HWECalculator ac = ac_[kk];
			 PrintWriter[] pw = new PrintWriter[ac.types.size()];
			
		  int n = ac.noCat();
		 
		
		
		//  int n1 = (int)(Math.pow(n, 2)-((double)(n*(n+1)))/2.0);
		  for(int j=0; j<ac.types.size(); j++){
			  pw[j] = new PrintWriter(new BufferedWriter(new FileWriter(new File(f, ac.types.get(j)+".txt"))));
			  pw[j].print("loc\tsnpid");
			  pw[j].print("\t"+getAlleleName(ac.inner));
		  for(int ind1 =0; ind1<ac.type_len; ind1++){
			
			//  for(int ind2 =ind1+1; ind2<n; ind2++){
				  String type = ind1==1 ? "CN_only": 
					  (ind1==0  ? "overall" : "within_cn="+(ind1-1));
				  pw[j].print("\tpval_"+type);
				 
			//  }
			  } 
		  pw[j].println();
		  }
		 
		  
		  //Double[] res = new Double[];
		  for(int i=0; i<dc.loc.size(); i++){
			//if(Constants.savePhasedConfiguration()) ac.scoreChi1(i, true);
			  int k=0;
			  
			
			 
			  for(int j=0; j<ac.types.size(); j++){
				  pw[j].print(dc.loc.get(i)+"\t"+dc.snpid.get(i));
				  ac.getAlleles(i,j);
				  pw[j].print("\t"+getAlleleName(ac.count_allele1.probs));
				  Double[] sig = ac.getSignificance(i, j);
			  for(int ind1 =0; ind1<sig.length; ind1++){
				 // for(int ind2 =ind1+1; ind2<n; ind2++){
					  pw[j].print("\t"+String.format("%5.3g", sig[ind1]));
					  k++;
				  //}
			  }
			  pw[j].println();
			  }
			
		  }
		  for(int j=0; j<pw.length; j++){
			  pw[j].close();
		  }
			 }
	  }catch(Exception exc){
		  exc.printStackTrace();
	  }
	}
	  
	  
	 private static String getAlleleName(EmissionStateSpace sp) {
		StringBuffer sb = new StringBuffer(sp.getGenotype(0).toString());
		
		for(int i=1; i<sp.size(); i++){
			sb.append(";");
			sb.append(sp.getGenotype(i).toString());
		}
		return sb.toString();
	}
	 private static String getAlleleName(double[] d) {
			StringBuffer sb = new StringBuffer(String.format("%5.3g", round(d[0])).trim());
			
			for(int i=1; i<d.length; i++){
				sb.append(";");
				sb.append(String.format("%5.3g", round(d[i])).trim());
			}
			return sb.toString();
		}
	 public static double round(double d){
		 if(d < 0.001){
			 return 0;
		 }
		 else return d;
	 }
	 


	private int noCat() {
		return this.types.size();
	}


	public int type_len;
	  
	
	
	public HWECalculator(DataCollection mdc, int ploidy){
	    
		this.ploidy = ploidy;
           this.dc1 = mdc;
       this.ces = (CompoundEmissionStateSpace)Emiss.getSpaceForNoCopies(ploidy); 
    
     
       this.inner = ces.getMembers()[0];
    
       this.type_len = 1 +ces.cnLength();
       count_allele1 = new SimpleExtendedDistribution(inner.size());
       expected_count1 = new double[ces.size()];
    //   this.count_allele = new SimpleExtendedDistribution[type_len];
       this.expected_count = new double[2][];//[ces.size()];
       this.observed_count = new double[2][];//[ces.size()];
     //  count_allele[0] = new SimpleExtendedDistribution(inner.cnLength());
       expected_count[0] = new double[ces.cnLength()];
       observed_count[0] = new double[ces.cnLength()];
     
    	 //  count_allele[i] =new SimpleExtendedDistribution(inner.size());
    	   expected_count[1] = new double[ces.size()];
    	   observed_count[1] = new double[ces.size()];
      
    		 this.types.add(Constants.experiment());
           if(mdc instanceof MergedDataCollection){
        	
        	   for(int i=0; i<((MergedDataCollection)mdc).ldl.length; i++){
        		   this.types.add(((MergedDataCollection)mdc).ldl[i].name);
        	   }
           }
           numCases = new double[mdc.loc.size()][types.size()][ces.size()]; //controls
           
	    }
	
	
	
	
	//  public int type_len;
	   
	  int currentIndex;
	  
	  
	  /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	   
	    public void scoreChi1(int i, boolean useEmSt){
	    	if(currentIndex==i){
	    		return ;
	    	}
	    	else{
	    		currentIndex=i;
	    	}
	        //EmissionStateSpace emStSp = this.maf.getEmissionStateSpace();
	        CompoundEmissionStateSpace emStSp1 = (CompoundEmissionStateSpace) dc1.getEmStSpace();
	        
	       
	       
	        for(Iterator<String> it= dc1.getKeys().iterator(); it.hasNext();){
	            String key = it.next();
	            HaplotypeEmissionState emv = (HaplotypeEmissionState) dc1.dataL.get(key);
	       int index = emv.dataIndex();
	            double[][] counts =  numCases[i][index];
	        
	                HaplotypeEmissionState st = (HaplotypeEmissionState)dc1.dataL.get(key);
	                PseudoDistribution dist = st.emissions[i];
	                EmissionStateSpace emstsp = st.getEmissionStateSpace();
	          
	                Integer fixed = dist.fixedInteger();
	                if(fixed!=null){
	                    add(counts, (ComparableArray)st.getEmissionStateSpace().get(fixed.intValue()), 1.0);
	                }
	                else{
	                    double[] emiss = st.emissions[i].probs();
	                    double sum = Constants.sum(emiss);
//	                    if(!(st.emissions[i] instanceof SimpleExtendedDistribution  || st.emissions[i] instanceof IntegerDistribution)) throw new RuntimeException("!!");
	                    for(int jj =0; jj<emiss.length; jj++){
	                        add(counts, (ComparableArray)emstsp.get(jj),emiss[jj]/sum);
	                    }
	                }
	           
	            
	        }
	      
	      
	        
	       // return res;
	    }  */
	//  double[] cn_count;
	//  int[][] cn_alias;
//	public final int len1;
	  public List<String> types = new ArrayList<String>();
	
	  public final int ploidy;
	  /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public synchronized void scoreChi1(double[] dist, int category1,  int i){
	    	double sum = Constants.sum(dist);
	    	int category = category1 + 1;
	    	
	    	{
	    		if(category <numCases[i].length){
	               double[] counts =  numCases[i][category];
	               for(int jj =0; jj<dist.length; jj++){
	                    	counts[jj]+=dist[jj]/sum;
	               }
	    		}
	    	}
	    	{
	    		 double[] counts =  numCases[i][0];
	               for(int jj =0; jj<dist.length; jj++){
	                    	counts[jj]+=dist[jj]/sum;
	               }
	    	}
	        
	    }
	      
	   public void initialise(){
		  for(int i=0; i<this.numCases.length; i++){
			  for(int j=0; j<numCases[i].length; j++){
				  Arrays.fill( numCases[i][j],0);
			  }
		  }
	   }
	        
	       // return res;
	  CompoundEmissionStateSpace ces;
	  EmissionStateSpace inner;
	// final SimpleExtendedDistribution[] count_allele;
	 final SimpleExtendedDistribution count_allele1;
	 final double[] expected_count1;
	  final double[][] expected_count, observed_count;
	  //j is type
	    public Double[] getSignificance(int i,int category){
	    	getAlleles(i, category);
	    	Arrays.fill(expected_count1,0.0);
	    	double sum = Constants.sum(numCases[i][category]);
	    	for(int j=0; j<ces.haploSize(); j++){
	    		int geno = ces.getGenoForHaplopair(ces.getHaploPairFromHaplo(j));
	    		int[] ind = this.ces.getMemberIndices(j);
	    		double prob = 1.0;
	    		for(int k=0; k<ind.length; k++){
	    			prob*=count_allele1.probs[ind[k]];
	    		}
	    		expected_count1[geno]+=prob*sum;
	    	}
	    	
	    	
	    	return getSignificance(expected_count1, numCases[i][category]);
	   
	    	
//	  
	    }
	   
	    
	    public void getAlleles(int i, int category){
	    	count_allele1.initialise();
	    	
	    	for(int j=0; j<this.numCases[i][category].length; j++){
	    		int[] ind = this.ces.getMemberIndices(j);
	    		for(int k=0; k<ind.length; k++){
	    			count_allele1.addCount(ind[k],numCases[i][category][j]);
	    		}
	    	}
	    	count_allele1.transfer(0);
	    }
	   public Double[] getSignificance(double[] expected, double[] observed){
		   {
			   double s1 = Constants.sum(expected);
			   double s2 = Constants.sum(observed);
			   if(Math.abs(s2-s1)>0.01) throw new RuntimeException("!! "+s1+" "+s2);
		   }
		  
		   Double[] res = new Double[this.type_len];
		   {
				  
				  res[0] = sig(expected, observed,1.0);
		   }
		   
		   {
			   Arrays.fill(expected_count[0], 0.0);
			   Arrays.fill(observed_count[0], 0.0);
			 
			   for(int j=0; j<expected.length; j++){
				   int cn = ces.getCN(j);
				   expected_count[0][cn]+=expected[j];
				   observed_count[0][cn]+=observed[j];
				   
			   }
			  res[1] = sig(expected_count[0], observed_count[0], 1.0);
			 
		   }
		  
		   for(int k=1; k<ces.cnLength(); k++){
			 //now within each cn state - we ignore cn state = 0 because it is trivially in hwe
			   Arrays.fill(expected_count[1], 0.0);
			   Arrays.fill(observed_count[1], 0.0);
			   double sum_obs =0;
			   double sum_exp =0;
			   for(int j=0; j<expected.length; j++){
				   int cn = ces.getCN(j);
				   if(cn==k){
					   expected_count[1][j] = expected[j];
					   observed_count[1][j] = observed[j];
					   sum_exp+=expected[j];
					   sum_obs+=observed[j];
				   }
			   }
			   res[k+1] = sig(expected_count[1], observed_count[1], sum_obs/sum_exp);
		   }
		 	return res;
	   }

	/* scaling scales expected to observed  */
	   private Double sig(double[] expected, double[] observed, double scaling) {
		double sum=0;
		/*if(Constants.CHECK && Math.abs(Constants.sum(expected)-Constants.sum(observed))>0.01) {
			throw new RuntimeException("!! "+sum);
		}*/
		int degf = 0;
		for(int i=0; i<expected.length; i++){
			if(expected[i]>0){
				double exp_i = expected[i]*scaling;
				sum+= Math.pow(observed[i]-exp_i, 2)/exp_i;
				degf++;
			}
		}
		double res =  this.chisq.chi2prob(degf, sum); 
	//	if(res==0){
		//	throw new RuntimeException("!!");
		//}
		return res;
	}
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#scoreChi(double[][], double[][])
 
public Double[] scoreChi(double[][] ns, double[][] nc){
    Double[] result = new Double[ns[0].length];
    for(int i=0; i<result.length; i++){
          Double R=ns[0][i]+ns[1][i]+ns[2][i];
          Double N=nc[0][i]+nc[1][i]+nc[2][i]+ns[0][i]+ns[1][i]+ns[2][i];;
          Double r1=ns[1][i], r2=ns[2][i], n1=ns[1][i]+nc[1][i], n2=ns[2][i]+nc[2][i];
          Double chiSqaureValues1=N*Math.pow((N*(r1+2*r2)-R*(n1+2*n2)),2);
          Double chiSqaureValues2=R*(N-R)*(N*(n1+4*n2)-Math.pow((n1+2*n2),2));
          if(chiSqaureValues2==0.0) result[i]=0.0;
          else result[i]=chiSqaureValues1/chiSqaureValues2;
    }
    return result;
}*/

//   Regression lr;
  //  Double[] regression = new Double[4] ;
  //  Double[] regressionSig = new Double[4] ;
//  List<String>[][] caseIndiv;
 //  List<String>[][] controlIndiv;
   /* private void scoreRegression(int pos_index,  int phenIndex) {
       if(lr==null){
    	   int t = pheno.type[phenIndex];
           lr = this.pheno.type[phenIndex]==0 ?
        		   new LinearRegression(this):new LogisticRegression(this);
          regression =  new Double[4];
          regressionSig = new Double[4];
          Arrays.fill(regression, 0.0);
          Arrays.fill(regressionSig, 1.0);
          try{
          this.log = new PrintWriter(new BufferedWriter(new FileWriter(new File(this.chrom+"_log.txt"))));
          }catch(Exception exc){
        	  exc.printStackTrace();
          }
          
       }
       for(int i=0; i<regression.length; i++){
       if(i<3){
          regressionSig[i] = lr.calcSignificance(lr.calcLogLDiff(pos_index, phenIndex,i));
          regression[i] = lr.slope();
          if(regressionSig[i]<1e-5){
        	  log.println(this.pheno.phen.get(phenIndex)+" "+this.loc.get(pos_index)+" "+i+" "+regressionSig[i]);
        	  log.flush();
          }
       }
       else{
           regression[i]  =0.0;
           regressionSig[i] = 1.0;
       }
       }
      // System.err.println("he");
       
      // }
       //System.err.println(Arrays.asList("armitage "+Arrays.asList(getSignificance(true))));
       //System.err.println(Arrays.asList("regress "+Arrays.asList(regression))+" sig "+Arrays.asList(regressionSig));
    }*/
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#scoreChi(int)
    
    public Double scoreChi(int i){
        EmissionStateSpace emStSp = this.maf.getEmissionStateSpace();
        CompoundEmissionStateSpace emStSp1 = (CompoundEmissionStateSpace) this.getEmStSpace();
        double[][] ns = new double[3][emStSp.size()]; //cases [0,1,2] copies of allele
        double[][] nc = new double[3][ns[0].length]; //controls
        for(int k=0; k<ns.length; k++){
            Arrays.fill(ns[i],0);
            Arrays.fill(nc[i],0);
        }
        for(Iterator<String> it = this.getKeys().iterator(); it.hasNext();){
            String key = it.next();
            PhasedDataState dat = (PhasedDataState) data.get(key);
          if(true)   throw new RuntimeException("!!");
                      boolean cse = false;//this.phenotypes[dat.emissions[i].getDataIndex()];
                      double[][] counts = cse? ns : nc;

            ComparableArray comp = (ComparableArray)dat.getElement(i);
            int[] memberIndices = emStSp1.getMemberIndices(emStSp1.get(comp));
            for(int k=0; k<counts[0].length; k++){
                counts[count(memberIndices, k)][k]++;
            }
        }
        Double[] res = scoreChi(ns , nc);
        ChiSq ch = new ChiSq();
        return ch.chi2prob(1, res[Constants.getMax(res)]);
    }*/
/*    int currentPhenScIndex = -1;
int currentPosScIndex = -1;
int currentType = -1;
//public static  boolean scoreRegression = false;
public synchronized String getPhenInfo(String string , int pos_index, int phenIndex, int type){
  int numCl = numClasses(phenIndex);
 if(pos_index!=currentPosScIndex || currentPhenScIndex!=phenIndex ){
 	Arrays.fill(regression, 0.0);
 	Arrays.fill(regressionSig, 1.0);
    // System.err.println(scoreRegression);
     if(OptionBuild.scoreRegression())
     		this.scoreRegression(pos_index, phenIndex);
     if(OptionBuild.scoreChi())   
     	this.scoreChi1(pos_index, true, phenIndex);
     
     currentPhenScIndex = phenIndex;
     currentPosScIndex = pos_index;
     currentType = type;
 }
  if(string.startsWith("chisq")){
      
      if(numCl==2) return Format.sprintf("%5.3g", new Double[] {getSignificance(false, type)});
      else return "";
  }
  else if(string.startsWith("armitage")){
  //    this.scoreChi1(pos_index, true, string.endsWith("state"), phenIndex);
      if(numCl==2)  return Format.sprintf("%5.3g",new Double[] { getSignificance(true, type)});
      else return "";
  }
  else if(string.startsWith("regress")){
      if(regression==null) return "null";
     return Format.sprintf("%5.3g", new Double[] {this.regression[type]});
     
  }
  else if(string.startsWith("regrP")){
      if(regressionSig==null) return "null";
      return Format.sprintf("%5.3g", new Double[] {this.regressionSig[type]});
      
   }
  else{
     
      if(string.startsWith("odds_")){
      ///  System.err.println(Arrays.asList(odds[i]));
          if(numCl==2)  return Format.sprintf(formatOdds[type], this.odds[type]);
          else return "";
      }
      else if(string.startsWith("cases_")){
        
          if(numCl==2)  return Format.sprintf(formatOdds[type], this.cases[type]);
          else return "";
      }
      else if(string.startsWith("controls_")){
          if(numCl==2)  return Format.sprintf(formatOdds[type], this.controls[type]);
          else return "";
      }
      else throw new RuntimeException( "!! " +string);
  }
 
}*/

