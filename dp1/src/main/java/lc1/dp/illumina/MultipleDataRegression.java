/**
 * 
 */
package lc1.dp.illumina;

import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import lc1.stats.ProbabilityDistribution2;
import lc1.stats.TrainableBinomialDistr;
import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;
import pal.math.UnivariateFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MultipleDataRegression{
	
	protected final IlluminaProbR illuminaProbRB;
	ProbabilityDistribution2[][] prob2 ;
	   int numObsT;
	   DoubleMatrix2D X;
	   public DoubleMatrix2D Y;
	   public DoubleMatrix2D YB, covar;
	   DoubleMatrix2D resultR;
	  int[][] rows;
	 // int[] cols, cols_var;
	  final double[] pseudo; // loc var
	  
	//  int[][]
	  
	//  int noEl, noElVar;
	 // int first;
	// int[][] all_fix_rows;
	//  int[][] lrr_fix_rows;
	 // int[][]non_skew_rows;
	 // int[][] b_fix_rows;
	  int[] pos;
	   
	 
	  double[] priorMod = null;
	  /** note if pseudoFill is less than 0 we use prior */
	  public MultipleDataRegression(IlluminaProbR[] illuminaProbRB, 
			   int[] pos, 
			   double[] pseudo1, double pseudoFill, Prior prior
			 ){
		  this.pos = pos;
		  this.index = illuminaProbRB[0].index;
		 this.priorMod = new double[illuminaProbRB.length];
		 for(int k=0; k<priorMod.length; k++){
			 priorMod[k] = Constants.pseudoMod1(illuminaProbRB[k].index);
		 }
		  this.prior  =prior;
		  this.pseudoFill = pseudoFill;
		  boolean addPseudoFill = Math.abs(pseudoFill)>1e-5 ? true: false;
		   this.len_extra = pos.length-1;
	//	   this.pseudo1 = new double[len_extra+illuminaProbRB[0].basis_mean.length()];
		    this.illuminaProbRB = illuminaProbRB[0];
			prob2 = new ProbabilityDistribution2[illuminaProbRB.length][];
			for(int i=0; i<prob2.length; i++){
				prob2[i] = 
					pos[i]>=0 ? illuminaProbRB[i].r[pos[i]] : 
						illuminaProbRB[i].rGlobal;
			}
			this.dists = new ProbabilityDistribution2[illuminaProbRB.length];
		    this.pseudo = pseudo1;//*Constants.pseudoMod(illuminaProbRB[0].index);
		//    this.probeOnly =pos[0]>=0 ? illuminaProbRB[0].probeOnly[pos[0]] : (pos[0]==-2 ? true : false);	
		//    all_fix_rows = new int[this.illuminaProbRB.all_indices.length][];
		//    this.lrr_fix_rows = new int[this.illuminaProbRB.reg_indices.length][];
		   // int[] lrr_fix_len = new int[lrr_fix_rows.length];
		   // int[] b_fix_len = null;
		    int[] all_fix_len = new int[this.illuminaProbRB.all_indices.length];
		   // int[] nonSkewRowCount=new int[this.illuminaProbRB.covar_indices.length];
		   /*if(!probeOnly){
			   this.b_fix_rows = new int[this.illuminaProbRB.bfrac_indices.length][];
			  b_fix_len = new int[b_fix_rows.length];
		   }*/
		    numObsT = 0; //number of different observations assigned to different cn/b states
			rows = new int[this.illuminaProbRB.all_indices.length][];
		    for(int i=0; i<prob2[0].length; i++){
			   int noCop = this.illuminaProbRB.emstsp.getCN(i);
			  // int k_b = this.illuminaProbRB.bfrac_alias(i);
			   //int k_r = this.illuminaProbRB.reg_alias[i];
			   int k_all = this.illuminaProbRB.all_alias[i];
			   //int skew = this.illuminaProbRB.covar_alias[i];
			   if(toinclude.contains(noCop)){
				   for(int k=0; k<prob2.length; k++){
					   int num = ((ProbabilityDistribution2) prob2[k][i]).numObs()+
					   (addPseudoFill ? 1 : 0);
					 //  lrr_fix_len[k_r]+=num;
					   //if(!probeOnly) b_fix_len[k_b]+=num;
					   //nonSkewRowCount[skew]+=num;
				       numObsT+=num;
				       all_fix_len[k_all]+=num;
			       }
			   }
		  }
		//    len_X = probeOnly ? 4 : 7;
		    len_X = this.illuminaProbRB.basis_mean.length();
	//	  this.res = new double[len_X-1];
		   	X = new DenseDoubleMatrix2D(numObsT,len_X +len_extra);
			initialiseMats(numObsT);//, b_fix_len, nonSkewRowCount);
		/*	for(int i=0; i<this.lrr_fix_rows.length; i++){
				this.lrr_fix_rows[i] = new int[lrr_fix_len[i]];
			}*/
			for(int i=0; i<this.rows.length; i++){
				rows[i] = new int[all_fix_len[i]];
			}
			
			
	   }
	   
	   
	   protected void initialiseMats(int numObsT) {
			Y = new DenseDoubleMatrix2D(numObsT,1);
		
	 }


//	final boolean probeOnly;
	   
	 final int len_extra,len_X;
	//private final DoubleMatrix2D x0R_probeOnly;
//	private final DoubleMatrix2D x0R;
	//private final DoubleMatrix2D x0B;
	 
	 static final SortedSet<Integer> toinclude = new TreeSet<Integer>();
	 static final SortedSet<Integer> toincludeMLzero = new TreeSet<Integer>();
	 static public final SortedSet<Integer> toincludeML = new TreeSet<Integer>();
	 static {
		 int[] regr = Constants.regress();
		 int[] ml = Constants.mltrain();
		 if(regr!=null){
		 for(int i=0; i<regr.length; i++){
			 toinclude.add(regr[i]);
		 }
		 }
		 if(ml!=null){
		 for(int i=0; i<ml.length; i++){
			 toincludeML.add(ml[i]);
		 }
		 }
		 if(toincludeML.remove(0)){
			 toincludeMLzero.add(0);
		 }
		
		 //toincludeML.remove(0);
	 }
		public static void switchSet() {
			  Logger.global.info("swithc");
			 toincludeML.addAll(toinclude);
			 toinclude.clear();
			
		}
	 
	 public void setAbsValues(int col){
		 if(true) throw new RuntimeException("!!");
		 int len = X.rows();
		 for(int i=0; i<len; i++){
			 double v=  X.getQuick(i, col);
			 X.setQuick(i, col, Math.abs(v));
		 }
	 }
	   public int  fillMatrices(double pseudoForFill){
		   int numObs=0;
		   int[] all_fix_len = new int[this.rows.length];
		   for(int i=0; i<prob2[0].length; i++){
			   
			   double  noCop =this.illuminaProbRB.emstsp.getCN(i);
			   double bCop =  this.illuminaProbRB.emstsp.getBCount(i); 
			  if(toinclude.contains((int)noCop) ){
				   double[] vals = this.illuminaProbRB.basis_mean.getVals(noCop, bCop);
				   int k_all = this.illuminaProbRB.all_alias[i];
				   for(int k=0; k<prob2.length; k++){
		  		      int numObs1=((ProbabilityDistribution2) (prob2[k][i])).fill(X, Y, YB,numObs, 
		  		    		 vals,pseudoForFill);
		  		      if(k>0){
		  		    	for(int j=numObs; j<numObs+numObs1; j++){
		  		    		  X.setQuick(j, this.len_X+(k-1), 1.0); //sets values based on which group
		  		    	  }
		  		      }
		  		    for(int j=0; j<numObs1; j++){
		  		    	   this.rows[k_all][all_fix_len[k_all]+j] = numObs+j;
		  		       }
		  		      all_fix_len[k_all]+=numObs1;
		  		     numObs+=numObs1;
				   }
		  		 
			   }
  	       
		   }
		   return numObs;
	   }
	 //  final double[] res ;
	  
	   
	   
	   public int  fillMatricesVariance(double pseudoFill){
		   int numObs=0;
		 //  int[] lrr_fix_len = new int[lrr_fix_rows.length];
		  //  int[] b_fix_len = null;
		   // int[] nonSkewRowCount=new int[this.illuminaProbRB.covar_indices.length];
		    //if(!probeOnly){
		    //b_fix_len = new int[b_fix_rows.length];
		   // }
		   for(int i=0; i<prob2[0].length; i++){
			   double  noCop =this.illuminaProbRB.emstsp.getCN(i);
			  if(toinclude.contains((int)noCop) ){
				 // int k_b = this.illuminaProbRB.bfrac_alias(i);
				//	 int k_r = this.illuminaProbRB.reg_alias[i];
				//	int skew_ind= this.illuminaProbRB.covar_alias[i];
				   for(int k=0; k<prob2.length; k++){
		  		      int numObs1=((ProbabilityDistribution2) (prob2[k][i])).fillVariance(Y,YB, covar, numObs,pseudoFill);//YB,numObs, nocopt,
		  		
		  		      //if(!probeOnly){
		  		  //     for(int j=0; j<numObs1; j++){
		  		   // 	   this.b_fix_rows[k_b][b_fix_len[k_b]+j] = numObs+j;
		  		    //   }
		  		     //  b_fix_len[k_b]+=numObs1;
		  		    
		  		     
		  		    //	 for(int j=0; j<numObs1; j++){
		  		    		 
		  		    //	   this.non_skew_rows[skew_ind][nonSkewRowCount[skew_ind]+j] = numObs+j;
		  		    //	 }
		  		    //	 nonSkewRowCount[skew_ind]+=numObs1;
		  		    
		  		     // }
		  		      //for(int j=0; j<numObs1; j++){
		  		    //	   this.lrr_fix_rows[k_r][lrr_fix_len[k_r]+j] = numObs+j;
		  		     //  }
		  		      //lrr_fix_len[k_r]+=numObs1;
		  		      numObs+=numObs1;
		  	
				   }
		  		  
			   }
  	       
		   }
		   return numObs;
	   }
	   
	   public DoubleMatrix2D regress(int startPos, 
			   RegressParams regr,
			   DoubleMatrix2D Y, int[] rows
			   //, int i
			  ){
		 
		   int numObs=startPos;
		 
		   if(numObs!=numObsT) {
			   throw new RuntimeException("!!!");
		   }
		   DoubleMatrix2D res = regr.regress(startPos, this, Y,rows);
		 if(Constants.CHECK)  check(regr, this.rows);
		   return res;
	   }
	   
	  
	   
	
	
	   public void setMidPoints(int type,int xOrY, RegressParams params , int[] indices){
		
		 //  System.err.println("regression is "+params.toString());
			  params.setX(0,1.0);
			   for(int i1=0; i1<indices.length; i1++){
				   int i = indices[i1];
				   double noCop = this.illuminaProbRB.emstsp.getCN(i);
				   double  bCop = this.illuminaProbRB.emstsp.getBCount(i); 
				   if(toinclude.contains((int)noCop)){
					   double[] vals =type==0 ? this.illuminaProbRB.basis_mean.getVals(noCop, bCop):
						   this.illuminaProbRB.basis_var.getVals(noCop, bCop);
					   for(int iii=0; iii<vals.length; iii++){
						   params.setX(iii, vals[iii]);
					   }
					 
					  
					  for(int k=0; k<prob2.length; k++){
						  if(i<prob2[k].length){
							  ProbabilityDistribution2 dist =    ((ProbabilityDistribution2) (prob2[k][i]));
							  for(int j=0; j<len_extra; j++){
								  params.setX(this.len_X+j, 0.0);
							  }
							  if(k>0) params.setX( this.len_X+(k-1), 1.0);
							  double val = params.calculate();
							  if(!Double.isNaN(val)){
								  if(type==1){
									 if(xOrY==0) val = Math.max(minbafvar2, val);
									 else if(xOrY==1) val = 	  Math.max(minlrr2, val);
										  
								  }
							   dist.setParam(type, xOrY,val);
							  }
							  else{
								  Logger.global.warning("is nan param");
							  }
						  }
						
					  }
				   }
				   }
			   
	   }
	   
	  /* public void setMidPoints(int type,int xOrY, RegressParams params ,DoubleMatrix2D paramR, DoubleMatrix2D paramB, int[] indices){
			
			 //  System.err.println("regression is "+params.toString());
				  params.setX(0,1.0);
				   for(int i1=0; i1<indices.length; i1++){
					   int i = indices[i1];
					   double noCop = this.illuminaProbRB.emstsp.getCN(i);
					   double  bCop = this.illuminaProbRB.emstsp.getBCount(i); 
					   if(toinclude.contains((int)noCop)){
						   double noCopt = this.illuminaProbRB.trans(noCop);
						   double bfract = this.illuminaProbRB.transform(this.illuminaProbRB.calcFracB(noCop,bCop));
						   params.setX(1,noCopt);
						   params.setX(2, Math.pow(noCopt, 2));
						   params.setX(3, noCop);
						   if(!probeOnly){
							   params.setX(4, bfract);
							   params.setX(5,illuminaProbRB.trans1(bfract, type>=1) );
							   params.setX(6,illuminaProbRB.trans2(bfract, type>=1) );
						   }
						  for(int k=0; k<prob2.length; k++){
							  if(i<prob2[k].length){
								  ProbabilityDistribution2 dist =    ((ProbabilityDistribution2) ((Mixture2)prob2[k][i]).dist[0]);
								  for(int j=0; j<len_extra; j++){
									  params.setX(this.len_X+j, 0.0);
								  }
								  if(k>0) params.setX( this.len_X+(k-1), 1.0);
								  double val = params.calculate();
								  double valR = params.calculate(paramR);
								  double valB = params.calculate(paramB);
								  double val1 = val/Math.sqrt(valR*valB);
								  if(!Double.isNaN(val1)){
							 	   dist.setParam(type, xOrY,val1);
								  }
								  else{
									  Logger.global.warning("is nan param");
								  }
							  }
							
						  }
					   }
					   }
				   
		   }*/
	   
	   
	   double minbafvar2 = Math.pow(Constants.minBafVar(), 2);
	   double minlrr2 = Math.pow(Constants.minLRRVar(), 2);
	   
	/*   public void setMidPointsExp(int type,int xOrY, RegressParams params , int[] indices){
			
			 //  System.err.println("regression is "+params.toString());
				  params.setX(0,1.0);
				   for(int i1=0; i1<indices.length; i1++){
					   int i = indices[i1];
					   double noCop = this.illuminaProbRB.emstsp.getCN(i);
					   double  bCop = this.illuminaProbRB.emstsp.getBCount(i); 
					   if(toinclude.contains((int)noCop)){
						   double noCopt = this.illuminaProbRB.trans(noCop);
						   double bfract = this.illuminaProbRB.transform(this.illuminaProbRB.calcFracB(noCop,bCop));
						   params.setX(1,noCopt);
						   params.setX(2, Math.pow(noCopt, 2));
						   params.setX(3, noCop);
						   if(!probeOnly){
							   params.setX(4, bfract);
							   params.setX(5,illuminaProbRB.trans1(bfract,type>=1) );
							   params.setX(6,illuminaProbRB.trans2(bfract,type>=1) );
						   }
						  for(int k=0; k<prob2.length; k++){
							  if(i<prob2[k].length){
								  ProbabilityDistribution2 dist =    ((ProbabilityDistribution2) ((Mixture2)prob2[k][i]).dist[0]);
								  for(int j=0; j<len_extra; j++){
									  params.setX(this.len_X+j, 0.0);
								  }
								  if(k>0) params.setX( this.len_X+(k-1), 1.0);
								  double val = Math.exp(params.calculate());
								  if(!Double.isNaN(val)){
								   dist.setParam(type, xOrY,val);
								  }
								  else{
									  Logger.global.warning("is nan param");
								  }
							  }
							
						  }
					   }
					   }
				   
		   }*/
	  
	   
	  
	   final ProbabilityDistribution2[] dists ; //new ProbabilityDistribution2[prob2.length];
	   public void trainOnNonRegress(SortedSet<Integer> toincludeML, boolean locx, boolean locy, boolean varx, boolean vary, boolean modPrior,
			   double[][]priors
			){
		   if(toincludeML.size()==0 && this.pseudo[0]>0 && this.pseudo[1]>0)return;
		   double[] pseudoR; 
		 
		  
		  
		   
		   double[] pseudoB = priors[2];// Constants.b_train(1,0);
	
		   double[] pseudoRho; 
		   if((toincludeML.size()==0 && pseudo[0]==0 && pseudo[1]==0) || toincludeML.first()>0){
		        pseudoR= priors[1];//Constants.r_train(1,0);
			   pseudoRho = priors[4];//Constants.rho_train(1);
		   }
		   else{
			
			   pseudoR= priors[0];//Constants.r_train0[0];
				   pseudoRho =priors[4];//Constants.rho_train0;
		   }
//		   if(Constants.getMinValue(pseudoR)>1e4 && Constants.getMinValue(pseudoB)>1e4 
	//			   && Constants.getMinValue(pseudoRho)>1e4
				   
		//   ) return;
		   
		   if(illuminaProbRB instanceof IlluminaProbRBCounts) {
			   maximise(prob2[0],((IlluminaProbRBCounts) this.illuminaProbRB).errors[pos[0]]);
		   }
		   else{
		   for(int i=0; i<prob2[0].length; i++){
			   double noCop = this.illuminaProbRB.emstsp.getCN(i);
			   
			   
			   
			  if(toincludeML.contains((int)noCop) || (pseudo[0]==0 && pseudo[1]==0)){
				
			inner:  for(int k=0; k<this.prob2.length; k++){
				//  if(true) break inner;
				  dists[k] =  prob2[k][i];
				  if(modPrior){
					  
						  if(locx) dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 0, true);
						  if(varx){
							  dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 1, true);
							  dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 2, true);
						  }
						  if(locy) dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 0, false);
						  if(vary){
							  dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 1, false);
							  dists[k].setPriors( ((ProbabilityDistribution2)illuminaProbRB.rGlobal[i]), 2, false);
						  }
				  }
				
				dists[k].maximise(
						  !locx ? 1e10 : pseudoR[0]*pseudo[0]*priorMod[k], 
						  !varx  ? 1e10 : pseudoR[1]*pseudo[1]*priorMod[k],
						  !varx  ? 1e10 : pseudoR[2]*pseudo[1]*priorMod[k],
						  !locy ? 1e10 :  pseudoB[0]*pseudo[0]*priorMod[k],
		    		      !vary  ? 1e10 : pseudoB[1]*pseudo[1]*priorMod[k], 
		    		      !vary  ? 1e10 :pseudoB[2]*pseudo[1]*priorMod[k],
		    		        !vary ? 1e10 :  pseudoRho[0]*pseudo[1]*priorMod[k]);
							  }
			  
			      
			       
			   
			 }
		  
	   }
		   }
	   }
	   
	   
	   private void add(double[] from, double[] to) {
		for(int i=0; i<from.length; i++){
			
				to[i] +=from[i];
			
		}
		
	}
	/*public void setVariance(){
		 
		  double[][] rho = new double[overall_rho_counts.length][prob2.length];
		  double[][] sigma_y1 = new double[baf_var_counts.length][prob2.length];
		  double[][] sigma_x1 = new double[lrr_var_counts.length][prob2.length];
		  if(!probeOnly){
			  for(int k=0; k<baf_var_counts.length; k++){
				  double totx = pseudoB[0][5]*pseudo;
				  double expx = totx *Math.pow(this.illuminaProbRB.stddevPriorB[k],2);
				  for(int j=0; j<baf_var_counts[k].length; j++){
					  totx+=baf_var_counts[k][j][1];
					  expx+=baf_var_counts[k][j][0];
				  }
				  for(int j=0; j<baf_var_counts[k].length; j++){
					  double totx1= pseudoB[0][6]*pseudo;
					  double expx1 = totx1 *(expx/totx);
					  totx1+=baf_var_counts[k][j][1];
					  expx1+=baf_var_counts[k][j][0];
					  sigma_y1[k][j] = Math.sqrt(expx1 / totx1);
				  }
			  }
			  for(int k=0; k<overall_rho_counts.length; k++){
				  double noCop = this.illuminaProbRB.emstsp.getCN(k);
				   if(this.illuminaProbRB.toInclude.contains((int)noCop)){
				   int k_r = this.illuminaProbRB.reg_alias[k];
				   int k_b = illuminaProbRB.bfrac_alias[k];
				   double totxy = pseudoRho[0]*pseudo;
				   double expxy = 0; //fix rho_prior of zero
				   double totx = totxy;
				   double expx = totx * Math.pow(this.illuminaProbRB.stddevPriorR[k_r],2);
				   double toty = totxy;
				   double expy = toty *Math.pow(this.illuminaProbRB.stddevPriorB[k_b],2);
				   for(int j=0; j<overall_rho_counts[k].length; j++){
					   totxy+=overall_rho_counts[k][j][1];
					   expxy+=overall_rho_counts[k][j][0];
					   toty+=baf_var_counts1[k][j][1];
					   expy+=baf_var_counts1[k][j][0];
					   totx+=lrr_var_counts1[k][j][1];
					   expx+=lrr_var_counts1[k][j][0];
				   }
				   double sigma_x_prior = Math.sqrt(expx/totx);
				   double sigma_y_prior = Math.sqrt(expy/toty);
				   double sigma_rho_prior = (expxy/totxy)/(sigma_x_prior*sigma_y_prior);
				   for(int j=0; j<overall_rho_counts[k].length; j++){
					   double totxy1 = pseudoRho[1]*pseudo;
					   double expxy1 = totxy1 * sigma_rho_prior *sigma_x_prior*sigma_y_prior; //fix rho_prior of zero
					   double totx1 = totxy1;
					   double expx1 = totx1 * Math.pow(sigma_x_prior,2);
					   double toty1 = totxy1;
					   double expy1 = toty1 *Math.pow(sigma_y_prior,2);
					   totxy1+=overall_rho_counts[k][j][1];
					   expxy1+=overall_rho_counts[k][j][0];
					   toty1+=baf_var_counts1[k][j][1];
					   expy1+=baf_var_counts1[k][j][0];
					   totx1+=lrr_var_counts1[k][j][1];
					   expx1+=lrr_var_counts1[k][j][0];
					   double sigma_x_prior1 = Math.sqrt(expx1/totx1);
					   double sigma_y_prior1 = Math.sqrt(expy1/toty1);
					   rho[k][j] =  (expxy1/totxy1)/(sigma_x_prior1*sigma_y_prior1);
					   
				   }
				   }
				   
			  }
				   
			  
		  }
		  for(int k=0; k<lrr_var_counts.length; k++){
			  double totx = pseudoR[0][5]*pseudo;
			  double expx = totx * Math.pow(this.illuminaProbRB.stddevPriorR[k],2);
			    for(int j=0; j<lrr_var_counts[k].length; j++){
			    	totx+=lrr_var_counts[k][j][1];
			    	expx+=lrr_var_counts[k][j][0];
			    }
			    for(int j=0; j<lrr_var_counts[k].length; j++){
					  double totx1= pseudoR[0][6]*pseudo;
					  double expx1 = totx1 *(expx/totx);
					  totx1+=lrr_var_counts[k][j][1];
					  expx1+=lrr_var_counts[k][j][0];
					  sigma_x1[k][j] = Math.sqrt(expx1 / totx1);
				  }
			   
		  }
		  
	//	  System.err.println("rho is "+rho);
		   for(int i=0; i<prob2[0].length; i++){
			   double noCop = this.illuminaProbRB.emstsp.getCN(i);
			   if(this.illuminaProbRB.toInclude.contains((int)noCop)){
				
				   for(int k=0; k<prob2.length; k++){
					   double sigmax =sigma_x1[this.illuminaProbRB.reg_alias[i]][k];
					   double sigmay = sigma_y1[this.illuminaProbRB.bfrac_alias[i]][k];
					   double rho_ = rho[i][k];
					 
				  ((ProbabilityDistribution2 ) ((Mixture2)prob2[k][i]).dist[0]).setPriorVar(sigmax, sigmay,rho_);//sigma_xy/(sigma_x*sigma_y));
				 if(pseudoR[1][0]*pseudo < 1e4 || pseudoB[1][1]*pseudo < 1e4 || pseudoRho[2]*pseudo < 1e4){
				
					 ((ProbabilityDistribution2) ((Mixture2)prob2[k][i]).dist[0])
				   .maximise(pseudoR[1][0]*pseudo, pseudoR[1][1]*pseudo, pseudoR[1][2]*pseudo,
            		   pseudoB[1][0]*pseudo, pseudoB[1][1]*pseudo, pseudoB[1][2]*pseudo, pseudoRho[2]*pseudo);
				 }
				   }
				  //System.err.println( "SIGMA"+((TrainableNormal2)  ((ProbabilityDistribution2) ((Mixture2)prob2[i]).dist[0])).sigma_y);
			   }
		   }
		 
	   }*/

	   final static boolean print = true;

	 public RegressParams[] paramR, paramRVar;
	final double pseudoFill;
	public  Prior prior;
	public static  boolean modifyBySign = false;
	
	public  void regress(String type, RegressParams[] paramBaf, 
			DoubleMatrix1D[] prior,
			DoubleMatrix2D yB, int[][] rows, int[][] indices, double[] ds, double pseudoM, int st){
		for(int i=0; i<rows.length; i++){
	 		  if(indices[i].length==0 && indices[i][0]==0){
	 			continue;
	 		  }
	 		 // else if(prior==null) continue;
	 		  else if(rows[i].length>0){
	 		    paramBaf[i] = RegressParams.getRegressParams(type,len_extra, len_X,pseudoM, ds);
	 			if(pseudoM>0 )  illuminaProbRB.transfer(paramBaf[i], prior[i]);
	 			regress(st,  paramBaf[i], yB, rows[i]);
	 		//  if(print)   System.err.println("regress  "+paramBaf[i]);
	 		  }
		}
	}
	
	//final double[] pseudo1;
	
	/*private double[] process(double[] ds) {
		for(int i=0; i<len_extra; i++){
			pseudo1[i] = ds[0];
		}
		System.arraycopy(ds, 1, pseudo1, len_extra, ds.length-1);
		return pseudo1;
	}*/

	public void set(){
		if(Constants.regressMean(0)){
		for(int i=0; i<rows.length; i++){
			if(paramR[i]!=null) setMidPoints(0, 0,paramR[i] , this.illuminaProbRB.all_indices[i]);
		}
		}
	}
	
	public void setVar(){
		if(Constants.regressVariance(0)){
		for(int i=0; i<rows.length; i++){
			 if(paramRVar[i]!=null)  setMidPoints(1, 0,paramRVar[i], this.illuminaProbRB.all_indices[i] );
		}
		}
	}
	
	public void regressMean(int st, double[][] priors, double pseudoM, boolean set){
		  if(Constants.regressMean(0)){
			  paramR = new RegressParams[rows.length];
			  this.regress("R", paramR, 
					  prior.priorR, Y,rows, illuminaProbRB.all_indices, priors[0], pseudoM, st);
			 
		
		  }
	}
	
	public void regressVar(int st, double[][] priors, double pseudoM, boolean set){
		 if(Constants.regressVariance(0)) {
			  paramRVar =new RegressParams[rows.length];
	  	 		regress("RVAR", this.paramRVar, 
	  	 				 prior.priorRVar, Y, rows, illuminaProbRB.all_indices, priors[1], pseudoM, st);
	  	 	   
	  	 	 }
	}

	final int index;
	/**priors are [r/b,cov, type] */
	 public boolean makeRegress(
			 boolean set, 
			 boolean global
			) {
		  double[][]priors =  Constants.priors(this.index);
		   double[] pseudoM = this.pseudo;
		  int st = fillMatrices(pseudoFill);
		this.regressMean(st, priors, pseudoM[0],set);
		  if(set) this.set();	
	      if(Constants.regressVariance(0) || Constants.regressVariance(1)){
	 	  int st1 =  fillMatricesVariance(pseudoFill);
	 	 this.regressVar(st1, priors, pseudoM[1], set);
	 	if(set) this.setVar();	
	 	  if(st!=st1) {
	 		  throw new RuntimeException("!! ");
	 	  }
	      }
	   if(set){
 		   boolean modPrior = true;
					this.trainOnNonRegress(toincludeML, true, true, true, true, modPrior, priors);
					this.trainOnNonRegress(toincludeMLzero, true, true, true, true, modPrior, priors);
	   }
		   return false;
		
	}


	/*private void regressCovVar(double pseudoM, boolean set, int st, double[] ds) {
		  for(int i=0; i<this.non_skew_rows.length-1; i++){
				
			  RegressParams paramCovar_;
			 
			  if(illuminaProbRB.covar_indices[i].length==1 && illuminaProbRB.covar_indices[i][0]==0){
		 		//	 paramCovar_ = paramCovar0;
		 		//	 throw new RuntimeException("!!");
				  continue;
		 			// paramBafVar_.setConstantPrior(Math.pow(illuminaProbRB.zeroBafPrior(),2));
		 		 }
		 		 else{
		 			//if(true) continue;
		 			 paramCov[i] = RegressParams.getRegressParams("Cov",len_extra,probeOnly,len_X,pseudoM,ds);
		 			 paramCovar_= paramCov[i];
		 		
		 			 //paramBafVar_.setConstantPrior(Math.pow(this.illuminaProbRB.stddevPriorB(i),2));
		 		 }
				if(pseudoM>0)   illuminaProbRB.transfer(paramCovar_, prior.priorRho[i]);
			//  DoubleMatrix2D paramBafVar_ = regress(st, paramCovar_, this.YB, this.non_skew_rows[i], pos);
			//  DoubleMatrix2D paramRVar_ = regress(st, paramCovar_, this.Y, this.non_skew_rows[i], pos);
			  regress(st, paramCovar_, this.covar, this.non_skew_rows[i]);
			
		//	  setMidPoints(1, 2,paramCovar_,paramRVar_, paramBafVar_, this.illuminaProbRB.covar_indices[i]);//, r_sigmaX, r_sigmaY);
			  if(set)	  setMidPoints(1, 2,paramCovar_, this.illuminaProbRB.covar_indices[i]);//, r_sigmaX, r_sigmaY);
			 // }
			  }
		
	}*/
	/*private void regressBVar(double pseudoM, boolean set, int st, double[] ds) {
		 for(int i=0; i<this.b_fix_rows.length; i++){
				//	if(b_fix_rows[i].length==0) continue;  
						  
					   RegressParams   paramBafVar_;
					   int[] bfrac_indices = illuminaProbRB.bfrac_indices(i);
					   if(bfrac_indices.length==1 && illuminaProbRB.bfrac_indices(i)[0]==0){
						   continue;
				 		 }
				 		 else{
				 			 paramBafVar[i] = RegressParams.getRegressParams("BAF var", len_extra, probeOnly, len_X,pseudoM,
				 					 ds);
//				 					 Constants.b_train(0, 1));
				 			 paramBafVar_= paramBafVar[i];
				 		//	 illuminaProbRB.
				 			// Variance v1 = this.illuminaProbRB.bafVar_prior(it,i);
			if(pseudoM>0)	 			 illuminaProbRB.transfer(paramBafVar_, prior.priorBafVar[i]);//setBafPriors(paramBafVar_, i);
				 			
				 		 }
					
					   regress(st, paramBafVar_, YB, this.b_fix_rows[i]);
					   if(print) System.err.println("regress B var "+paramBafVar_);
					   if(set)  setMidPoints(1, 1,paramBafVar_, this.illuminaProbRB.bfrac_indices(i) );
				   }
		
	}*/
	/*private void regressRVar(double pseudoM, boolean set, int st, double[] ds) {
		// RegressParams paramRVar0 = null;
	 	 for(int i=0; i<this.lrr_fix_rows.length;i++){
	 		// if(lrr_fix_rows[i].length==0) continue;
	 		//double totVar = (sum+Math.pow(this.illuminaProbRB.stddevPriorR[illuminaProbRB.reg_alias[i]],2)*pseudoVar)/(pseudoVar+sum_denom);
	 		RegressParams paramRVar_;
	 	  		
	 		int[] reg_ind = illuminaProbRB.reg_indices[i];
	 		 if(reg_ind.length==1 && reg_ind[0]==0){
	 			 continue;
	 			 //	 paramRVar_ = paramRVar0;
	 		//	 paramRVar_.setConstantPrior(0,this.illuminaProbRB.lrrVarPrior(it).zeroVar);
	 		//	 throw new RuntimeException("!!");
	 			// paramRVar_.setConstantsPriorWeight(1e-3);
	 		 }
	 		 else{
	 			paramRVar[i] =  RegressParams.getRegressParams("RVar", len_extra, probeOnly, len_X,pseudoM, ds);
	 		 	
	 			 paramRVar_= paramRVar[i];
	 			if(pseudoM>0) illuminaProbRB.transfer(paramRVar_, 
	 					probeOnly ? prior.priorRVarProbeOnly[i] : prior.priorRVar[i]);//setLRRPriors(paramRVar_, i, probeOnly);
	 		
	 			// if(!probeOnly) paramRVar_.setConstantPrior(1,this.illuminaProbRB.lrrVarPrior(it).var.get(1));
	 		 }
	 		
		 	
	 	//	  paramRVar_.setConstantPrior(totVar);
	 		  regress(st, paramRVar_, Y, this.lrr_fix_rows[i]);
	 		  if(print) System.err.println("regress R var "+paramRVar_);
	 		 if(set)  setMidPoints(1, 0,paramRVar_, this.illuminaProbRB.reg_indices[i] );
	 	  }
		
	}*/
	private int sumlength(int[][] rows2) {
		int len =0;
		for(int i=0; i<rows2.length; i++){
			len+=rows2[i].length;
		}
		return len;
	}
	private void trainvariance(Set<Integer> toinclude2) {
		// TODO Auto-generated method stub
		
	}


	private double sum(DoubleMatrix2D y2,int i) {
		// TODO Auto-generated method stub
		double sum=0;
		for(int k=0; k<y2.rows(); k++){
			sum+=y2.getQuick(k, i);
		}
		return sum;
	}


	private void check(RegressParams paramR,  int[][] rows) {
		DoubleMatrix2D regressR = paramR.regressR;
		  if(Double.isNaN(regressR.get(0, 0))){
			//  Logger.global.warning("regression failed "+rows.length+"at" +DataCollection.datC.snpid.get(pos)); 
			  Logger.global.warning("details \t"+paramR.toString());
			  Logger.global.warning(Constants.print(this.sum(this.X, paramR.cols, rows[1])));
			  Logger.global.warning("\n");
			  //throw new RuntimeException("is NaN");
		   }
		
	}


	private double[] sum(DoubleMatrix2D x2, int[] cols, int[] rows2) {
		double[] sum  = new double[cols.length];
		for(int i=0; i<rows2.length; i++){
			for(int k=0; k<sum.length; k++){
				sum[k] +=x2.getQuick(rows2[i], cols[k]);
			}
		}
		return sum;
	}


	public void bst(String[] st) {
	
	}


	public void rst(String[] res) {
		
		if(paramR!=null ){
			StringBuffer rst = new StringBuffer();
	       for(int i=0; i<paramR.length; i++){
	    	   if(paramR[i]!=null && this.rows[i].length>0){
	    	   paramR[i].getFormulae(rst, this.illuminaProbRB.basis_mean);
	    	   rst.append("\n");
	    	   }
	       }
	       res[0] = rst.toString();
		}
		if(paramRVar!=null){
			StringBuffer rst = new StringBuffer();
		       for(int i=0; i<paramRVar.length; i++){
		    	   if(paramRVar[i]!=null && this.rows[i].length>0){
		    	   paramRVar[i].getFormulae(rst, this.illuminaProbRB.basis_var);
		    	   rst.append("\n");
		    	   }
		       }
		       res[1] = rst.toString();
			}
	
//	       return rst.toString();
	}
	

	   
	   private void maximise(final ProbabilityDistribution2[] prob22, double[] initialValue) {
		//   if(true) throw new RuntimeException("!!");
		final double[] noB = new double[prob22.length];
		final double[] noCop = new double[prob22.length];
		/*double[] errorA2B = new double[prob22.length];
		double[] errorB2A = new double[prob22.length];*/
		final boolean[] include = new boolean[prob22.length];
		int sze =0;
		for(int k=0; k<noB.length; k++){
			noB[k] = this.illuminaProbRB.emstsp.getBCount(k);
			noCop[k] = this.illuminaProbRB.emstsp.getCN(k);
			TrainableBinomialDistr dis = (TrainableBinomialDistr) prob22[k];
			if(dis.size()>0){
				include[k] = true;
				sze+=dis.size();
			}
/*			if(noB[k] == 0){
				errorB2A[k] = 1.0;
				errorA2B[k] = 0.0;
			}else if(noB[k] == 0){
				errorB2A[k] = 1.0;
				errorA2B[k] = 0.0;
			}
		*/
		}
		if(sze==0) return;
		   MultivariateFunction funct = new MultivariateFunction(){

			@Override
			/* arg[0] is B2A and arg[1] is A2B */
			public double evaluate(double[] arg0) {
				double val = 0;
				double A2B = arg0[0]; double B2A = arg0[1];
				//double A2A = 1 - A2B;
				double B2B = 1-B2A;
				for(int k=0; k<prob22.length; k++){
					if(include[k]){
					double p = (noB[k] * B2B + (noCop[k] - noB[k])*A2B)/noCop[k];
//							noB[k] ==0 ?  A2B : (noB[k]==2 ? 1-B2A : 0.5+(A2B-B2A)/2.0);
				
					val+= ((UnivariateFunction) prob22[k]).evaluate(p);
					}
				}
		//	System.err.println(A2B+" "+B2A+" "+val);
				return val;
			}

			@Override
			public double getLowerBound(int i) {
				// TODO Auto-generated method stub
				return 0;
			}

			@Override
			public double getUpperBound(int i) {
				// TODO Auto-generated method stub
				return 0.5;
			}

			@Override
			public int getNumArguments() {
				return 2;
			}

			@Override
			public OrthogonalHints getOrthogonalHints() {
				// TODO Auto-generated method stub
				return null;
			}
			
		};
		
		 MultivariateMinimum uv = new OrthogonalSearch();
		 uv.findMinimum(funct, initialValue);
	       for(int k=0; k<prob22.length; k++){
           prob22[k].recalcName();
        //   System.err.println(prob22[k]);
	       }
	       
		
	}

   }