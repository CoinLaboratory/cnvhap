package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.ArmitageTrendTest;
import lc1.stats.BetaLikelihood;
import lc1.stats.ChiSq;
import lc1.stats.Contigency;
import lc1.stats.PseudoDistribution;
import lc1.stats.StateDistribution;
import lc1.util.Constants;
import lc1.util.PhenoGroup;

public class ArmitageCalculator extends AssociationCalculator {
	
	int case_index, int_control_index;
	
	    
	    public BetaLikelihood ct;
	    public Contigency ct1 ;
	    public ArmitageTrendTest att;
	    public double[][][][] numCases; //[pos][pheno_cat][type][copy _number]
        
	    public Double[] res;
	    public Double[] odds;
	    public Double[] sum;
		private final Double[] odds1;
		private Double[] sum1;
	    //public Double[] cases;
	   // public Double[] controls;
	   // public String[] formatOdds;
	    public static void add(double[][] counts, ComparableArray comp, double prob){
	        int noCop= comp.noCopies(true);
	       // int noB = (int) comp.noB();
	        //int noA = noCop - noB;
	        counts[0][noCop]+=prob;
	      //  counts[1][noA]+=prob;
	      //  counts[2][noB]+=prob;
	      //  if(noCop==2){
	        //    counts[3][noA]+=prob;
	       // }
	       
	    }
	    
	    public List<String> types = null;


		private List<String> type_assoc_header;
	    public String type(int ind1){
	    	/*	if(this.dc1 instanceof MergedDataCollection){
	    			return ((MergedDataCollection)dc1).ldl[ind1].name;
	    		}
	    		else*/
	    	  return types.get(ind1);
	    	  }
	 
	    public Double[][] oddsRatio(int i1, int overall,int j){
	        
	       int[] ind = this.type_decomp.get(overall);
	      
	         // double sum_cases = Constants.sum(numCases[i1][ind1][j]);
	          //double sum_controls = Constants.sum(numCases[i1][ind2][j]);
	          //double sum= sum_cases+sum_controls;
	    	Double[] odds = j==0 ? this.odds : odds1;
	    	Double[] sum = j==0 ? this.sum : sum1;
	              for(int i=0; i<odds.length; i++){
	                 double  cases  = numCases[i1][ind[0]][j][i];
	                double controls = numCases[i1][ind[1]][j][i];
	                  sum[i] = cases+controls;
	                  double  p = sum[i]==0 ? 0 :  cases/ sum[i];
	                  if(p==0) odds[i] = -maxOdds;
	                  if(p==1) odds[i] = maxOdds;
	                  else{
	                	  double odds_ = Math.max(-maxOdds,Math.log(p / (1-p))/log2);
	                	  odds[i] = Math.min(maxOdds,odds_);
	                  }
	              }
	            
	              return new Double[][]{odds, sum};
	    }
	  
	  
	  
	 
	 
	 
	  
	    public  void printResults1(File dir
		
			  ){
		  try{
			  dir.mkdir();
			//  BufferedReader br = new BufferedReader(new FileReader(pheno));
			//  String[] types = br.readLine().split("\t");
			//  
					 ArmitageCalculator ac = this;// (ArmitageCalculator)ac_[ii];
					//if(ac==null) continue;
				//  phenoFile==null ? 
					//	  new ArmitageCalculator((MergedDataCollection)dc):
						//  new ArmitageCalculator(dc, phenoFile, header);
			  int n = ac.noCat();
			  DataCollection dc = this.dc1;
			  int len1 = (int) ((double)(n*(n-1))/2.0);
			  PrintWriter[] pw = new PrintWriter[len1];
			  PrintWriter[][] pw_counts = new PrintWriter[len1][ac.type_len];
			// File[] dirs = new File[len1];
			//  int n1 = (int)(Math.pow(n, 2)-((double)(n*(n+1)))/2.0);
			 int k=0;
			  for(int ind1 =0; ind1<n; ind1++){
				
				  for(int ind2 =ind1+1; ind2<n; ind2++){
					  String nme = ac.types.get(ind1)+"_"+ac.types.get(ind2);
					  nme.replace('=', '_');
					  pw[k] = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, nme+".txt"))));
					  pw[k].print("chrom\tloc\tsnpid\tquality");
					  File dir_out = new File(dir, nme);
					  dir_out.mkdir();
					  for(int j=0; j<ac.state_in; j++){
						
						
						 
   						  String type = ac.type_assoc.get(j);
   						 pw_counts[k][j] 
						               = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir_out, type+".txt"))));
					  		pw[k].print("\tpval_"+type);
					  		pw_counts[k][j].println("snpid\t"+ac.types.get(ind1)+"_"+ac.type_assoc_header.get(j)+"\t"+
					  				ac.types.get(ind2)+"_"+ac.type_assoc_header.get(j));
					 
					  }
					  for(int j=ac.state_in; j<ac.type_len; j++){
						  String type =  "state_"+(j-ac.state_in);
					  		pw[k].print("\tpval_"+type);
					  		 pw_counts[k][j] 
							               = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir_out, type+".txt"))));
					  		 String head ="0;1;2"; 
					  		pw_counts[k][j].println("snpid\t"+ac.types.get(ind1)+"_"+head+"\t"+
					  				ac.types.get(ind2)+"_"+head);
					  }
					  pw[k].println();
					  k++;
				  }
			  }
			  
			  
			  //Double[] res = new Double[];
			  for(int i=0; i<dc.loc.size(); i++){
				if(Constants.savePhasedConfiguration()) ac.scoreChi1(i, true);
				   k=0;
				  
				  
				
				  for(int ind1 =0; ind1<n; ind1++){
					  for(int ind2 =ind1+1; ind2<n; ind2++){
						  Double minQ = dc.dc==null ?  null : dc.dc.minQuality(i);
						  pw[k].print(dc.chrom+"\t"+dc.loc.get(i)+"\t"+dc.snpid.get(i)+"\t"+
								(minQ==null ? "NA":  String.format("%5.3g", minQ)));
										  
										 
						  for(int j=0; j<ac.type_len; j++){
							  pw[k].print("\t"+String.format("%5.3g", ac.getSignificance(i, k, j)));
							  pw_counts[k][j].println(dc.snpid.get(i)+"\t"+ac.getCountString(i,k,j));
						 
						  }
						  pw[k].println();
						  k++;
				  }
				  }
				 
			  }
			  for(int j=0; j<pw.length; j++){
				  pw[j].close();
				  for(int k1=0; k1<pw_counts[j].length; k1++){
					  pw_counts[j][k1].close();
				  }
			  }
		//		}
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
		}
	  
	  
	 
	    protected  void reset(){
			cn_alias = null;
		}
	    public void init(){
		  int modelLength = hmm.modelLength();
		  
		  
		  
		  EmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(this.ploidy);
		  int modelLength1 =  ((CompoundMarkovModel)hmm).getMemberModels()[0].modelLength();
          int new_type_len =modelLength1+this.state_in-1;
     //     int len1 = Constants.backgroundCount()+1;
          if(Constants.saveStates() ){
          	this.type_len=new_type_len;// + (Constants.saveStates() ? Constants.modify0.length : 0);
          	for(int i=0; i<numCases.length; i++){
          		for(int j=0; j<numCases[i].length; j++){
          			int len = numCases[i][j][0].length;
          			numCases[i][j] = new double[type_len][];
          			
          				//if(this.type_assoc.get(k).equals("overall"))
          				numCases[i][j][0] = new double[ emstsp.genoListSize()];
          				numCases[i][j][1] = new double[ emstsp.cnLength()];
          				
          				for(int k=2; k<state_in; k++){
          					int[][] alias_al = emstsp.aliasNB();
           					int[]res = emstsp.aliasNB()[k-1];
          					numCases[i][j][k] = new double[res.length];
          					
          				//else if(typ)
          				}
          			for(int k=state_in; k<type_len; k++) {
          				numCases[i][j][k]=new double[len1];
          			}
          			}
          			
          		}
          
          	res = new Double[type_len];
        
          this.cn_alias = new int[type_len][hmm.modelLength()];
         
          for(int k=state_in; k<type_len; k++){
        	  cn_alias[k][0] =-1;
          }
        	  //this.cn_count = new double[this.dc1.stSp1[1].copyNumber.size()];
          for(int ik=0; ik<this.state_in; ik++){
              cn_alias[ik] = new int[emstsp.size()];
     
          }
      
          for(int i=0; i<emstsp.size(); i++){
          	 cn_alias[0][i] = emstsp.getGenoForHaplopair(i);
          	 cn_alias[1][i] = emstsp.getCN(i);
          	for(int k=1; k<=this.maxC; k++){
          		if(cn_alias[1][i]!=k){
          			cn_alias[k+1][i] =-1;
          		}
          		else{
          		  cn_alias[k+1][i] = emstsp.getBCount(i);
          		}
          	}
          }
          
          int[] count = new int[modelLength1];
         
         
        	  for(int j=1; j<modelLength; j++){
        		  Arrays.fill(count, 0);
        		EmissionState[] st =  ((CompoundState) hmm.getState(j)).getMemberStates(true);
        		
        		for(int k1=0; k1<st.length; k1++){
        			count[st[k1].getIndex()]++;
        		}
        		for(int k=state_in; k<type_len; k++){
        			cn_alias[k][j] = count[k-state_in+1];
        		}
        	  }
          }
         // }
	  }
	 List<int[]> type_decomp = new ArrayList<int[]>();
	  public ArmitageCalculator(MergedDataCollection mdc, int ploidy){
		  super(Arrays.asList(mdc.getUnderlyingDataSets()).toString());
	    	 ct = new BetaLikelihood();
	    	 double[] d1 =new double[mdc.ldl.length];
	    	 for(int i=0; i<d1.length; i++){
	    		d1[i] =  mdc.ldl[i].dataL.size();
	    	 }
	    	 this.ploidy = ploidy;
	    	 init(ploidy);
             att = new ArmitageTrendTest();
           this.dc1 = mdc;
          len1 = ploidy+1;
          types = new ArrayList<String>();
     //     type_assoc_header = new ArrayList<String>();
          for(int i=0; i<mdc.ldl.length; i++){
        	  types.add(mdc.ldl[i].name);
          }
          for(int i=0; i<types.size(); i++){
        	  for(int j=i+1; j<types.size(); j++){
        		  types_all.add(types.get(i)+"_"+types.get(j));
        		  type_decomp.add(new int[] {i,j});
        	  }
          }
           int len = Emiss.getSpaceForNoCopies(ploidy).copyNumber.size();
         for(Iterator<EmissionState> it = mdc.dataLvalues(); it.hasNext();){
        	 HaplotypeEmissionState nxt = (HaplotypeEmissionState) it.next();
        	 this.keysToIndex.put(nxt.getName(), nxt.dataIndex());
         }
           this.type_len=1;// + (Constants.saveStates() ? Constants.modify0.length : 0);
             numCases = new double[mdc.loc.size()][mdc.ldl.length][type_len][len]; //controls
            res = new Double[type_len];
           
 	          
 	           odds = new Double[len];
 	         sum = new Double[len];
 	         this.odds1 = new Double[len1];
 	         this.sum1 = new Double[len1];
 	      
           
	    }
	

	  public void init(int ploidy){
			maxC = Constants.maxCopies*ploidy;
			type_assoc = new ArrayList<String>();
			type_assoc_header =  new ArrayList<String>();
			EmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(ploidy);
			type_assoc_header.add(emstsp.defaultList.toString());
			type_assoc.add("overall");
			type_assoc.add("CN");
			type_assoc_header.add(emstsp.copyNumber.toString());
			for(int i=1; i<=maxC;i++){
				type_assoc.add("within CN = "+i);
				type_assoc_header.add(getString(emstsp,emstsp.getGenoForCopyNo(i)));
			}
			state_in = type_assoc.size();
			this.rank = new boolean[type_assoc.size()];
			Arrays.fill(rank, true);
			rank[0] = false;
			for(int i=1; i<=maxC; i++){
				rank[i] = false;
			}
			
		}
	  private String getString(EmissionStateSpace emstsp, int[] genoForCopyNo) {
		StringBuffer b = new StringBuffer(emstsp.getBCount(genoForCopyNo[0]));
		for(int i=1; i<genoForCopyNo.length;i++){
			b.append(",");
			b.append(emstsp.get(genoForCopyNo[i]));
		}
		return b.toString();
	}

	int maxC;
	
	  
	  
	  
	 
	  
	  public static class Pheno{
		  public Pheno(File[] pheno, String header){
			  
		  }
	  }
	  
	  
	  public ArmitageCalculator(DataCollection mdc, File[] pheno, 
			  String header, int ploidy,
			  String stratHead, String stratVal) throws Exception{
		  super(header+"_"+stratHead+"_"+stratVal);
	    	 ct = new BetaLikelihood();
	    	 init(ploidy);
	    	 this.ploidy = ploidy;
          att = new ArmitageTrendTest();
        this.dc1 = mdc;
       len1 = ploidy+1;
       types = new ArrayList<String>();
       boolean isnumeric = true;
       System.err.println("pheno files are "+Arrays.asList(pheno));
       Map<String, String> keysToValue = new HashMap<String, String>();
      // Map<String, Integer> type_map = new HashMap<String, Integer>();
       Map<String,String> [] rename = mdc.getRenaming();
     inner1: for(int ii=0; ii<pheno.length; ii++){
    	 if(pheno[ii]==null || !pheno[ii].exists()) continue inner1;
    	 BufferedReader br = new BufferedReader(new FileReader(pheno[ii]));
       List<String> l = Arrays.asList(br.readLine().split("\t"));
       int patient_id = l.indexOf("PATIENT");
       int pheno_id = l.indexOf(header);
     if(pheno_id<0) {
    	  // throw new RuntimeException("no index of "+header+"\n"+l);
    	   continue inner1;
      }
       int strat_id = stratHead==null ? -1 : l.indexOf(stratHead);
       String st = "";
      outer: for(int i=0; (st= br.readLine())!=null; i++){
    	   if(st.length()==0) continue;
    	   String[] str = st.split("\t");
    	   String pat = str[patient_id];
    	   if(stratHead!=null && !str[strat_id].equals(stratVal)){
    		   continue outer;
    	   }
    	   if(rename!=null && rename[ii]!=null && rename[ii].containsKey(pat)){
    		   pat = rename[ii].get(pat);
    	   }
    	   EmissionState st1 = mdc.dataL.get(pat);
    	  if(st1==null){
    		  
    		//  if(mdc.name.startsWith("SLEGEN") && pat.indexOf("Inf")>=0){
    			  inner: for(Iterator<String> it = mdc.dataL.keySet().iterator(); it.hasNext();){
    				 String key =  it.next();
    				 if(key.startsWith(pat) || pat.startsWith(key)){
    					 pat = key;
    					 st1 = mdc.dataL.get(pat);
    					 break inner;
    				 }
    			  }
    		//  }
    		 
    	  }
    	 if(st1==null || st1.noCop()!=ploidy) continue outer;
    	   keysToValue.put(pat, str[pheno_id]);
    	   if(i==0){
    		   try{
    			  Double.parseDouble(str[pheno_id]);
    		   }catch(Exception exc){
    			   isnumeric = false;
    		   }
    	   }
    	  
       }
       
       br.close();
     }
       Set<Object> vals = new HashSet<Object>(keysToValue.values());
       if(isnumeric && vals.size()>4 ) throw new RuntimeException("should use linear assoc her");
     //  Collection vals1 = keysToValue.values();
       this.keysToIndex =
    	  convertToIndex(keysToValue, types);
        int len = Emiss.getSpaceForNoCopies(ploidy).copyNumber.size();
      
       
        
        this.type_len=1;// + (Constants.saveStates() ? Constants.modify0.length : 0);
          numCases = new double[mdc.loc.size()][types.size()][type_len][len]; //controls
         res = new Double[type_len];
        
	          
	           odds = new Double[len];
	         sum = new Double[len];
	         this.odds1 = new Double[len1];
	         this.sum1 = new Double[len1];
	         for(int i=0; i<types.size(); i++){
	        	  for(int j=i+1; j<types.size(); j++){
	        		  types_all.add(types.get(i)+"_"+types.get(j));
	        		  type_decomp.add(new int[] {i,j});
	        	  }
	          }
        
	    }
	  private Map<String, Integer> convertToIndex(
			Map<String, String> keysToValue, List<String> types2) {
		types2.clear();
		types2.addAll(new HashSet(keysToValue.values()));
		 Map<String, Integer> keysToIndex = new HashMap<String, Integer>();
		 for(Iterator<String> it = keysToValue.keySet().iterator(); it.hasNext();){
	      		String key = it.next();
	      		keysToIndex.put(key, types2.indexOf(keysToValue.get(key)));
		 }
		 return keysToIndex;
	}


	public static Map<String, Integer> convertToQuantiles(Map<String, String>keysToValue, List<String> types, PhenoGroup[] quantiles){
		  Map<String, Integer> keysToIndex = new HashMap<String, Integer>();
	
		  List<Double> l1 = new ArrayList<Double>();
		 
      	for(Iterator<String> it = keysToValue.values().iterator(); it.hasNext();){
      		String val = it.next();
      		double v = val.length()==0 || val.equals("NA") ? Double.NaN : Double.parseDouble(val);
      		if(!Double.isNaN(v)) l1.add(v);
      	}
      	Collections.sort(l1);
      	
    //  	type_map.clear();
      	for(Iterator<String> it = keysToValue.keySet().iterator(); it.hasNext();){
      		String key = it.next();
      		String val = keysToValue.get(key);
      		
      	//	double v = Double.parseDouble(va));
      		Integer v = getPos(val, quantiles);
      		if(v!=null)
	        		keysToIndex.put(key, v);
      		
      	}
      	types.clear();
    	for(int pos = 0; pos<quantiles.length; pos++){
    	
    	 types.add(quantiles[pos].toString());
    		
    	}
    
    
    //	if(keysToValue.values().contains(Double.NaN)){
    		types.add("NaN");
    	//}
      	return keysToIndex;
      	
	  }
	private static Integer getPos(String v, PhenoGroup[] mid) {
		
		for(int i=0; i<mid.length; i++){
			if(mid[i].within(v))   return i;
		}
		return null;
	}




	
	  int noCat(){
		  return numCases[0].length;
	  }
	   
	  int currentIndex;
	  public void update(StateDistribution dist){
		  
	  }
	  
	  /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public void scoreChi1(int i, boolean useEmSt){
	    	if(currentIndex==i){
	    		return ;
	    	}
	    	else{
	    		currentIndex=i;
	    	}
	        //EmissionStateSpace emStSp = this.maf.getEmissionStateSpace();
	        //CompoundEmissionStateSpace emStSp1 = 
	        //	Emiss.getSpaceForNoCopies(Constants.backgroundCount());
	        	//(CompoundEmissionStateSpace) dc1.getEmStSpace();
	        
	       
	       
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
	    }
	//  double[] cn_count;
	   
	  
	  /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public synchronized void scoreChi1(StateDistribution dist, int i, boolean useEmSt, String key){
	    	if(cn_alias==null) init();
	    	double sum = Constants.sum(dist.dist);
	        //EmissionStateSpace emStSp = this.maf.getEmissionStateSpace();
	       // CompoundEmissionStateSpace emStSp1 = (CompoundEmissionStateSpace) dc1.getEmStSpace();
	        
	      Integer index = this.keysToIndex.get(key);
	  // if(Constants.CHECK){
	//	   EmissionState st = ( (MergedDataCollection)this.dc1).ldl[index].dataL.get(key);
	  
	 //   System.err.println(st.name);
	//   }
	     if(index==null) return;
	            double[][] counts =  numCases[i][index];
	        for(int k=state_in; k<this.type_len; k++){
	               for(int jj =1; jj<dist.dist.length; jj++){
	                    	counts[k][cn_alias[k][jj]]+=dist.dist[jj]/sum;
//	                        add(counts, (ComparableArray)emstsp.get(jj),);
	                    }
	        }
	        
	        }
	    
	    /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public synchronized void scoreChi1(double[] dist, int i, boolean useEmSt, String key){
	    	double sum = Constants.sum(dist);
	    	if(cn_alias==null) init();
	   // 	if(dist[Constants.getMax(dist)] < Constants.imputedThresh()*sum) return;
	        //EmissionStateSpace emStSp = this.maf.getEmissionStateSpace();
	       // CompoundEmissionStateSpace emStSp1 = (CompoundEmissionStateSpace) dc1.getEmStSpace();
	    //    if(Constants.getMax(dist))
	      Integer index = this.keysToIndex.get(key);
	 //  if(Constants.CHECK){
	//	   EmissionState st = ( (MergedDataCollection)this.dc1).ldl[index].dataL.get(key);
	  
	//    System.err.println(st.name);
	//   }
	     if(index==null) return;
	            double[][] counts =  numCases[i][index];
	           for(int k=0; k<state_in; k++){
	               for(int jj =0; jj<dist.length; jj++){
	            	   if(cn_alias[k][jj]>=0){
	                    	counts[k][cn_alias[k][jj]]+=dist[jj]/sum;
	            	   }
//	                        add(counts, (ComparableArray)emstsp.get(jj),);
	                    }
	           }
	        
	        }
	      
	   public void initialise(){
		  for(int i=0; i<this.numCases.length; i++){
		   for(int j=0; j< this.numCases[i].length; j++){
			   for(int k=0; k< this.numCases[i][j].length; k++){
				   Arrays.fill( numCases[i][j][k],0);
			   }
		   }
		  }
	   }
	        
	       // return res;
	  //final int assocTest;
	  
	

	  //j is type
	    public Double getSignificance(int i,  int overall_ind, int j){
	    	
	    	int assocTest =  Constants.assocTest();
	    	if(assocTest==3) assocTest=beta;
	    	if(j==0 ){//&& assocTest == armitage){
	    		assocTest = beta;
	    	}
	    	int[] ind = this.type_decomp.get(overall_ind);
	    
	    	double[] r1 = numCases[i][ind[0]][j];
	    	double[] r2 = numCases[i][ind[1]][j];
	        if(assocTest==armitage){
	            att.set(r1, r2);
	           res[j] = ChiSq.chi2prob1(1, att.chisq());
	        }
	        else if(assocTest==beta){
	            ct.setMatrix(new double[][] {r1, r2});
	            
	            res[j] =ct.getSig(); 
	            	//ChiSq.chi2prob1(1,ct.chisq());//getSig();
	        }
	        else if(assocTest==chisq){
	        	if(ct1==null){
	        		ct1= new Contigency();
	        	}
	        	ct1.setMatrix(new double[][] {r1, r2});
	        	res[j] = ct1.getSig();
	        }
	    /*    if(res[j]<0.0001 && ind[0]==0 && ind2==2){
	        	System.err.println(Constants.print(r1));
	        	System.err.println(Constants.print(r2));
	        }*/
	    return res[j];
	    }
	    
	    private int[] getInt(double[] r1) {
		int[] res = new int[r1.length];
		for(int i=0; i<res.length; i++){
			res[i] =(int) Math.round(r1[i]);
		}
		return res;
	}

		public String getCountString(int i,  int overall_ind, int j){
	    	int[] ind = this.type_decomp.get(overall_ind);
	    
	    	double[] r1 = numCases[i][ind[0]][j];
	    	double[] r2 = numCases[i][ind[1]][j];
	    	return getString(r1)+"\t"+getString(r2);
	    }

		private String getString(double[] r1) {
			StringBuffer sb = new StringBuffer(String.format("%5.3g", r1[0]).trim());
			for(int j=1; j<r1.length; j++){
				sb.append(",");
				sb.append(String.format("%5.3g", r1[j]).trim());
			}
			return sb.toString();
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

