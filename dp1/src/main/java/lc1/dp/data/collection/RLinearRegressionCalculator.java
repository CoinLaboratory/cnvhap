package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RList;
import org.rosuda.JRI.Rengine;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class RLinearRegressionCalculator extends AssociationCalculator {
	
	private double[] x;
	private double[] pheno;
	int numPhen;
//	public DoubleMatrix2D Y,Y1; //multiple phenotypes
//	public DoubleMatrix2D certainty; 
	
	Rengine re=new Rengine(new String[] {"--vanilla"}, false, null);	
	
	public void init(int ploidy){
	//	maxC = Constants.maxCopies*ploidy;
		type_assoc = new ArrayList<String>();
		type_assoc.add("CN");
		type_assoc.add("B-allele");
		/*for(int i=1; i<maxC;i++){
			type_assoc.add("within CN = "+i);
		}*/
		state_in = type_assoc.size();
		this.rank = new boolean[type_assoc.size()];
		Arrays.fill(rank, true);
		rank[0] = true;
		rank[1] = false;
	//	this.oddssum = new Double[2][type_assoc.size()];
	}
	    
	   
	   
	protected  void reset(){
		re.eval("rm(x)");
		//X_all = null;
	}
	  
	
	  
	 
	 
	  public  void printResults1(File dir1
		
			  ){
		  try{
			  
			  String[] top = new String[] {"pval","beta","se"};
			  dir1.mkdir();
			  File dir = new File(dir1, this.name);
			  dir.mkdir();
			//  BufferedReader br = new BufferedReader(new FileReader(pheno));
			//  String[] types = br.readLine().split("\t");
			 
			//	for(int ii=0; ii<ac_.length; ii++){
					RLinearRegressionCalculator ac =this;
				//	if(ac==null) continue;
				//  phenoFile==null ? 
					//	  new ArmitageCalculator((MergedDataCollection)dc):
						//  new ArmitageCalculator(dc, phenoFile, header);
			
			  int len1 = numPhen;
			  PrintWriter[][] pw = new PrintWriter[len1][top.length];
			 
			//  int n1 = (int)(Math.pow(n, 2)-((double)(n*(n+1)))/2.0);
			 int k=0;
			 DataCollection dc = this.dc1;
			  for(int ind1 =0; ind1<len1; ind1++){
				  	  File dir2 = new File(dir, ac.types_all.get(ind1));
				  	  dir2.mkdir();
				 
					 
					  for(int kk=0; kk<top.length; kk++){
						  pw[k][kk] = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir2, top[kk]+".txt"))));
						  pw[k][kk].print("chrom\tloc\tsnpid\tquality");
						  for(int j=0; j<ac.state_in; j++){
							  String type = ac.type_assoc.get(j);
						  		pw[k][kk].print("\t"+top[kk]+"_"+type);
						  }
						  for(int j=ac.state_in; j<ac.type_len; j++){
							  String type =  "state_"+(j-ac.state_in);
								pw[k][kk].print("\t"+top[kk]+"_"+type);
						  }
						  pw[k][kk].println();
					  }
					 
					  k++;
				 
			  }
			  
			  
			  //Double[] res = new Double[];
			  for(int i=0; i<dc.loc.size(); i++){
				if(Constants.savePhasedConfiguration()) throw new RuntimeException("!!");// ac.scoreChi1(i, true);
				   k=0;
				  
				  
				
				  for(int ind1 =0; ind1<len1; ind1++){
					 
						  
						  for(int kk=0; kk<top.length; kk++){
							  Double minQ = dc.dc==null ? null : dc.dc.minQuality(i);
							  pw[k][kk].print(dc.chrom+"\t"+dc.loc.get(i)+"\t"+dc.snpid.get(i)+"\t"+
										(  minQ==null ? "NA" :	  String.format("%5.3g",minQ)));
						  for(int j=0; j<ac.type_len; j++){
							
							 if(kk==0)pw[k][kk].print("\t"+String.format("%5.3g", ac.getSignificance(i,ind1, j)));
							 else if(kk==1) pw[k][kk].print("\t"+String.format("%5.3g", ac.oddssum[0]));
							 else if(kk==2) pw[k][kk].print("\t"+String.format("%5.3g", ac.oddssum[1]));
							  }
						 
						 // }
						  pw[k][kk].println();
						  } 
						  k++;
				 
				  }
				 
			  }
			  for(int j=0; j<pw.length; j++){
				  for(int kk=0; kk<pw[j].length; kk++){
				  pw[j][kk].close();
				  }
			  }
		  }catch(Exception exc){
			  exc.printStackTrace();
		  }
		}
	  
	  
	  
	 int Yrows; //number of individual
	
	
	  public void init(){
		 
		  
		  int modelLength = hmm.modelLength();
		  
		  EmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(this.ploidy);
		  int modelLength1 =  ((CompoundMarkovModel)hmm).getMemberModels()[0].modelLength();
		
          int new_type_len =modelLength1+this.state_in-1;
       //  todo.clear();
          //this.certainty = new DenseDoubleMatrix2D(Y.rows(), this.state_in);
          if(Constants.saveStates()){
        	  re.eval("X_all = array(0,c("+hmm.noSnps+","+Yrows+","+new_type_len+1+")");
    //    	  this.X_all = new DoubleMatrix2D[hmm.noSnps];
             
          	this.type_len=new_type_len;
          this.cn_alias = new int[type_len][hmm.modelLength()];
         
          for(int k=state_in; k<type_len; k++){
        	  cn_alias[k][0] =-1;
          }
        	  //this.cn_count = new double[this.dc1.stSp1[1].copyNumber.size()];
          for(int ik=0; ik<this.state_in; ik++){
              cn_alias[ik] = new int[emstsp.size()];
     
          }
      
          for(int i=0; i<emstsp.size(); i++){
          //	 cn_alias[0][i] = emstsp.getGenoForHaplopair(i);
          	 cn_alias[0][i] = emstsp.getCN(i);
          	for(int k=1; k<this.state_in; k++){
          		if(this.type_assoc.get(k).indexOf("within")>=0 && cn_alias[0][i]!=k){
          			cn_alias[k][i] =-1;
          		}
          		else{
          		  cn_alias[k][i] = emstsp.getBCount(i);
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
	  
	
	 
	
	
	
	  public RLinearRegressionCalculator(DataCollection mdc, File[] pheno, String[] header, int ploidy,
	String stratHead, String stratVal
	  ) throws Exception{
		  super(Arrays.asList(header)+"_"+stratHead+"_"+stratVal);
	    	 init(ploidy);
	    	 this.ploidy = ploidy;
      
        this.dc1 = mdc;
       len1 = ploidy+1;
       
      // Map<String, Integer> type_map = new HashMap<String, Integer>();
       List<int[]> pheno_id = new ArrayList<int[]>()  ;
       int[]patient_id = new int[pheno.length];
       int[] strat_ind = new int[pheno.length];
       List<List<Integer>> row_l = new ArrayList<List<Integer>>();
       int rows = 0;
       
       Map<String, String>[]trans = mdc.getRenaming();
     for(int ii=0; ii<pheno.length; ii++){
    	 if(pheno[ii]==null) continue;
       BufferedReader br = new BufferedReader(new FileReader(pheno[ii]));
       List<String> l = Arrays.asList(br.readLine().split("\t"));
       
	   if(ii==0){
		   
		    	for(int j=0; j<header.length; j++){
		    		 outer: for(int i=0; i<l.size(); i++){
		    		if(l.get(i).startsWith(header[j])){
		    			int[] id_ = new int[pheno.length];
		    			id_[0] = i;
		    			pheno_id.add(id_);
		    			types_all.add(l.get(i));
		    			row_l.add( new ArrayList<Integer>());
		    			continue outer;
		    		}
		    	}
		    }
	   }
	   else{
		   for(int i=0; i<pheno_id.size(); i ++){
			   pheno_id.get(i)[ii] = l.indexOf(types_all.get(i));
		   }
	   }
	this.Yrows = rows;
       patient_id[ii] = l.indexOf("PATIENT");
     strat_ind[ii] = stratHead==null ? -1 : l.indexOf(stratHead) ;
      
       String st = "";
     
       while((st=br.readLine())!=null){
    	   String[] str = st.split("\t");
    	   
    	   String pat =  str[patient_id[ii]];
    	 EmissionState emst =  mdc.dataL.get(pat);
    	   if(emst!=null){
    		 if( emst.noCop()==ploidy &&(stratHead==null || str[strat_ind[ii]].equals(stratVal))){
    		   rows++;
    		 }
    		   
    	   }
    	   else{
    		   if(stratHead==null || str[strat_ind[ii]].equals(stratVal)){
    	    		  if(pat.indexOf("Inf")>=0 && mdc.name.startsWith("SLEGEN")){
    	    			  inner: for(Iterator<String> it = mdc.dataL.keySet().iterator(); it.hasNext();){
    	    				 String key =  it.next();
    	    				 if(key.startsWith(pat)){
    	    					
    	    			    		//   rows++;
    	    			    		
    	    					trans[ii].put(pat, key);
    	    					 rows++;
    	    					
    	    					 break inner;
    	    				 }
    	    			  }
    	    		 
    	    		     
    	    		  }
    		   }
    	    Logger.global.warning("no data for "+pat);
    	    	 
    	   }
    	  
       }
       br.close();
     }
    this.numPhen = pheno_id.size();
    re.eval("Y = array(0,c("+rows+","+numPhen+")");
    re.eval("Y1 = array(0,c("+rows+","+numPhen+")");
   //    this.Y = new DenseDoubleMatrix2D(rows,pheno_id.size()); //first is CN
	//   this.Y1 =new DenseDoubleMatrix2D(rows, pheno_id.size());
	   {
	   int i=0;
	   for(int ii=0; ii<pheno.length; ii++){
		   if(pheno[ii]==null) continue;
	   BufferedReader br = new BufferedReader(new FileReader(pheno[ii]));
     String st =   br.readLine();
  //     boolean isnumeric = true;
      
      inner: for(; (st= br.readLine())!=null; ){
    	   if(st.length()==0) continue;
    	   String[] str = st.split("\t");
    	//   System.err.println(st);
    	 //  int ind = keysToIndex.g
    	   String pat = str[patient_id[ii]];
    	   if(trans[ii]!=null&& trans[ii].containsKey(pat)) pat = trans[ii].get(pat);
    		 EmissionState emst =  mdc.dataL.get(pat);
      	   if(emst!=null){
      		 if( emst.noCop()==ploidy &&(stratHead==null || str[strat_ind[ii]].equals(stratVal))){
    	   //if(mdc.dataL.containsKey(pat)){
    		//   if(stratHead==null || str[strat_ind[ii]].equals(stratVal)){
    			 //  if(pat.equals("471111")){
    				//   System.err.println(pat);
    			   //}
    			   if(keysToIndex.containsKey(pat)){
    				   continue inner;
    				   //throw new RuntimeException("trying to add "+pat+" twice!!!");
    			   }
	    	   keysToIndex.put(pat, i);
	    	   for(int k=0; k<pheno_id.size(); k++){
	    		   String st_ = str[pheno_id.get(k)[ii]];
	    		  
		    		   double d =  st_.length()==0  || st_.equals("NA") || st_.equals("null")? Double.NaN:
	    				   Double.parseDouble(st_);
		    		   re.eval("Y["+i+","+k+"] = "+d);//.setQuick(i, k, d);
		    		   if(!Double.isNaN(d)){
		    			   row_l.get(k).add(i);
		    		   }
	    		   
	    	   }
	    	   i++;
    		   }
    	   }
    	  
    	  
       }
       br.close();
	   }
	   }
       this.rows = new int[row_l.size()][];
       for(int i=0; i<this.rows.length; i++){
    	   this.rows[i] = new int[row_l.get(i).size()];
    	   for(int j=0; j<this.rows[i].length; j++){
    		   this.rows[i][j] = ((Integer)row_l.get(i).get(j)).intValue();
    	   }
       }
      
    //  int len = Emiss.getSpaceForNoCopies(ploidy).copyNumber.size();
     /*  this.logPNull = new double[pheno_id.size()];
      for(int ind =0; ind<pheno_id.size(); ind++){
    	  this.logPNull[ind] =this.logPNull(ind);
      }*/
       
        
        this.type_len=1;// + (Constants.saveStates() ? Constants.modify0.length : 0);
        
	    }
	/*  private Map<String, Integer> convertToIndex(
			Map<String, String> keysToValue, List<String> types2) {
		types2.clear();
		types2.addAll(keysToValue.values());
		 Map<String, Integer> keysToIndex = new HashMap<String, Integer>();
		 for(Iterator<String> it = keysToValue.keySet().iterator(); it.hasNext();){
	      		String key = it.next();
	      		keysToIndex.put(key, types2.indexOf(keysToValue.get(key)));
		 }
		 return keysToIndex;
	}


	public static Map<String, Integer> convertToQuantiles(Map<String, String>keysToValue, List<String> types, int no){
		  Map<String, Integer> keysToIndex = new HashMap<String, Integer>();
	
		  List<Double> l1 = new ArrayList<Double>();
		 
      	for(Iterator<String> it = keysToValue.values().iterator(); it.hasNext();){
      		String val = it.next();
      		double v = val.length()==0 || val.equals("NA") ? Double.NaN : Double.parseDouble(val);
      		if(!Double.isNaN(v)) l1.add(v);
      	}
      	Collections.sort(l1);
      	double[] interval = new double[no-1];
      	for(int i=0; i<interval.length; i++){
      		interval[i] = (double)(i+1)/(double)no;
      	}
      
      	double[] mid = new double[interval.length];
      	for(int k=0; k<mid.length; k++){
      		 mid[k] = l1.get((int)Math.round((double)l1.size()*interval[k]));
      	}
      	
    //  	type_map.clear();
      	for(Iterator<String> it = keysToValue.keySet().iterator(); it.hasNext();){
      		String key = it.next();
      		String val = keysToValue.get(key);
      		double v = val.length()==0 || val.equals("NA") ? Double.NaN : Double.parseDouble(val);
      	//	double v = Double.parseDouble(va));
      		
	        		keysToIndex.put(key, getPos(v, mid));
      		
      	}
      	types.clear();
    	for(int pos = 0; pos<mid.length; pos++){
    	
    		if(pos==0) types.add("<"+mid[pos]);
    		else types.add(mid[pos-1]+"<= x <"+mid[pos]);
    	}
    
    	types.add(">= "+mid[mid.length-1]);
    //	if(keysToValue.values().contains(Double.NaN)){
    		types.add("NaN");
    	//}
      	return keysToIndex;
      	
	  }
	private static Integer getPos(double v, double[] mid) {
		if(Double.isNaN(v)) return mid.length+1;
		for(int i=0; i<mid.length; i++){
			if(v<mid[i]) return i;
		}
		return mid.length;
	}*/



	
	 
	   
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
	    }   */
	//  double[] cn_count;
	
	
	
	
	
	  
	  /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public synchronized void scoreChi1(StateDistribution dist, int i, boolean useEmSt, String key){
	    	double sum = Constants.sum(dist.dist);
	     
	        
	      Integer index = this.keysToIndex.get(key);
	  //    double[] tmp = new double[type_len];
	     if(index==null) return;
	 //	X_all[i].setQuick(index, 0, 1);
	 	//tmp[0] = 1;
	 	  re.eval("X_all["+i+","+index+","+0+"] = 1");
	        for(int k=state_in; k<this.type_len; k++){
	        	double avg =0;
	        	double sum_ =0;
	               for(int jj =1; jj<dist.dist.length; jj++){
	            	   	avg += cn_alias[k][jj]*(dist.dist[jj]/sum);
	                    sum_+=	dist.dist[jj]/sum;
	                    }
	               re.eval("X_all["+i+","+index+","+(1+k)+"] = "+avg/sum_);
	               //tmp[1+k] = avg/sum_;
	   //            this.X_all[i].setQuick(index, 1+k,avg/sum_);
	        }
	       // re.assign("tmp", tmp);
	        re.eval("X_all["+i+","+index+",] = tmp");
	     //   if(index.intValue()==463){
	   	   // 	System.err.println("val1 "+X_all[i].getQuick(index, 0));
	   	    //}
	        }
	    
	  
	  //  Set<Integer> done = new TreeSet<Integer>();
	    /* (non-Javadoc)
	     * @see lc1.dp.data.collection.DataC#scoreChi1(int, boolean)
	     */
	    public synchronized void scoreChi1(double[] dist, int i, boolean useEmSt, String key){
	    //  if(X_all==null) init();
	    	
	    	double sum = Constants.sum(dist);
	    	//boolean add = (dist[Constants.getMax(dist)] >= Constants.imputedThresh()*sum);
	        
	      Integer index = this.keysToIndex.get(key);
	    
	     if(index==null) return;
	//    System.err.println(index+" "+key);
	   //  if(index.intValue()==0 && name.endsWith("_1")){
	   //	   	System.err.println("val "+X_all[i].getQuick(index, 0));
	   	    	
	  // 	 }
	     double[] tmp = new double[state_in];
	     tmp[0] = 1;
	     re.eval("X_all["+i+","+index+","+0+"] = 1");
//	     X_all[i].setQuick(index, 0, 1);
	           for(int k=0; k<state_in; k++){
	        	   double sum_=0;
	        	   double avg=0;
	               for(int jj =0; jj<dist.length; jj++){
	            	   if(cn_alias[k][jj]>=0){
	                    	sum_+=dist[jj]/sum;
	                    	avg+=(dist[jj]/sum)*cn_alias[k][jj];
	            	   }
//	                        add(counts, (ComparableArray)emstsp.get(jj),);
	                    }
	               re.eval("X_all["+i+","+index+","+(1+k)+"] = "+avg/sum_);
	               //this.X_all[i].setQuick(index, 1+k,avg/sum_);
	           }
	          // if(index.intValue()==10){
	   	    //	System.err.println("val "+X_all[i].getQuick(index, 0));
	   	    	
	   	   // }
	       //    if(i==0) done.add(index);
	           //System.err.println("done");
	        }
	  
	      
	
	        
	       // return res;
	  
	  //j is type
	  final  int[][] rows;
	  static DoubleMatrix2D Q = new DenseDoubleMatrix2D(2,2);
	  static DoubleMatrix2D x0 = new DenseDoubleMatrix2D(2,1);
	  static{
		  for(int i=0; i<2; i++){
				 Q.setQuick(i, i,1e-6);
			
				 x0.setQuick(i, 0, 0);
			 }
			// Q.setQuick(1, 1, 1e10);
	  }
	 // double[] avg;
	  static lc1.stats.NormalDistribution normal;
	// double[] logPNull;
	 
	/* public double logPNull(int ind){
		 double var0 = 0;
	    	double mean =0;
	    	for(int i1=0; i1<rows[ind].length; i1++){
		    	   mean+=Y.getQuick(rows[ind][i1], ind);	
		    	}
	    	mean = mean / (double)rows[ind].length;
	    	for(int i1=0; i1<rows[ind].length; i1++){
	    		//if(i1<10) System.err.println(i1+" "+var0);
		    	  var0+=Math.pow(Y.getQuick(rows[ind][i1], ind) -mean ,2);	
		    }
	    	var0 = Math.sqrt(var0/(double)rows[ind].length);
	    	double logpNull = 0;
	    	for(int i1=0; i1<rows[ind].length; i1++){
		    	  logpNull+=normal.logpdf(Y.getQuick(rows[ind][i1], ind) ,mean ,var0);	
		    }
	    	return logpNull;
	 }*/
	 double beta;
	 double se;
	//public static boolean CHECK = true;
	 
	 public Double getSignificance(int i, int ind, int j){
		 	re.assign("indeX", rows[ind]);
		 	re.assign("indeY", new int[] {0,j+1});
			re.assign("indeY1", new int[] {ind});
		 	re.eval("X_ = X_all[i,indeX,indeY");
		 	re.eval("Y_ = Y[indeX,indeY1]");
		 
//	    	DoubleMatrix2D X_ = X_all[i].viewSelection(rows[ind], new int[] {0,j+1});
//	    	DoubleMatrix2D Y_ = Y.viewSelection(rows[ind], new int[] {ind});
	    
	//    	 DoubleMatrix2D AT = lg.transpose(X_);
	//		   DoubleMatrix2D prod = lg.mult(AT, X_);
			  
	  //  	DoubleMatrix2D B = solve(X_,Y_,Q,x0,  AT,prod);
	    //	DoubleMatrix2D Yha = lg.mult(X_, B);
/*	    	double S = 0;
	    	
	    	for(int i1=0; i1<Y_.rows(); i1++){
	    		if(Constants.CHECK){
	    		double xv = X_.getQuick(i1, 0);
	    		if(xv!=1){
	    			boolean cont = keysToIndex.values().contains(ind);
	    			throw new RuntimeException("!! "+ind);
	    		}
	    		}
	    		//if(i1<10) System.err.println(i1+" "+S);
	    		double y1 = Y_.getQuick(i1, 0);
	    	//	double y11 =Y.getQuick(rows[ind][i1], ind);
	    		double y2 = Yha.getQuick(i1, 0);
	    		double p1 = Math.pow(y1-y2,2);
	    	
		    	  S+=p1;
		    	// var0+=p2;
		    }
	    	 
	    	
	         double p = X_.columns();
	         double n =  X_.rows();
	         double var = Math.sqrt(S/(double)n);
	      double   se = Math.sqrt((S/(n-p-1)) * lg.inverse(prod).getQuick(1,1));
	      double  beta = B.getQuick(1, 0);*/
		 	re.eval("summ= summary(Glm(Y_~X_))");
			REXP summ = re.eval("summ");
			RList names= re.eval("dimnames(summ$coefficients)").asList();
			double[][]coeff = summ.asList().at("coefficients").asMatrix();
			int indrow = Arrays.asList(names.at(0).asStringArray()).indexOf("x");
			String[] names1 = names.at(1).asStringArray();
	      this.oddssum[0][0] = coeff[indrow][0];
	      this.oddssum[1][0]= 1.0/coeff[indrow][1];
	    	/*double logPAlt = 0;
	    	for(int i1=0; i1<Y_.rows(); i1++){
		    	  logPAlt+= normal.logpdf(Y_.getQuick(i1, 0), Yha.getQuick(i1, 0) ,var);	
		    }
	    	double chisq = 2*(logPAlt - this.logPNull[ind]);
	    if(chisq<0){
	    	if(chisq<-0.1) {
	    		throw new RuntimeException("!! "+chisq);
	    	}
	    	 return 1.0;
	    }*/
	        double pv =  coeff[indrow][3];// ChiSq.chi2prob1(1, chisq);
	        
	       // if(pv<1 && j==0){
	       // 	Logger.global.info("h");
	        //}
	        
	        
	          return pv;
	       
	       
	  
	    }

Algebra lg = new Algebra();
	    private DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b,
				   DoubleMatrix2D Q, DoubleMatrix2D x0, 
				   DoubleMatrix2D AT,
				   DoubleMatrix2D prod) {
				
				  
				   
				   for(int i=0; i<prod.rows(); i++){
					   for(int j=0; j<prod.columns(); j++){
						   prod.setQuick(i, j, prod.getQuick(i, j)+Q.getQuick(i, j));
					   }
				   }
				   DoubleMatrix2D b1 = b.copy();
				   DoubleMatrix2D Ax0 = lg.mult(A, x0);
				   for(int i=0; i<b1.rows(); i++){
					   b1.setQuick(i, 0, b1.getQuick(i, 0)-Ax0.getQuick(i, 0));
				   }
				   DoubleMatrix2D res =  lg.solve(prod,lg.mult(AT,b1));
				   for(int i=0; i<res.rows(); i++){
					   double v1 = res.getQuick(i, 0);
					  
					   res.setQuick(i, 0, v1+x0.getQuick(i, 0));
				   }
				   
				   return res;
			 }
		
	    /*private DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b
				  ) {
				
				   DoubleMatrix2D AT = lg.transpose(A);
				   DoubleMatrix2D prod = lg.mult(AT, A);
				
				   DoubleMatrix2D res =  lg.solve(prod,lg.mult(AT,b));
				  
				   return res;
			 }*/

final Double[][] oddssum  = new Double[2][1];
		@Override
		public Double[][] oddsRatio(int i, int overall_ind, int j2) {
	//	if(Double.isNaN(oddssum[0][0]) || Double.isNaN(oddssum[1][0]) 
	//	|| 		Double.isInfinite(oddssum[0][0]) || Double.isInfinite(oddssum[1][0]) 
	//	) {
	//		throw new RuntimeException("!!");
	//	}
			return oddssum;
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

