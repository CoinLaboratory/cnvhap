package lc1.dp.illumina;

import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MultipleDataRegressionRB extends MultipleDataRegression {

	 
	
	   DoubleMatrix2D resultB;
	   
		 public  RegressParams[]  paramBaf, paramBafVar, paramCov;
	 
	
	public MultipleDataRegressionRB(IlluminaProbR[] illuminaProbRB, int[] pos,
			double[] pseudo1, double pseudoFill, Prior prior) {
		super(illuminaProbRB, pos, pseudo1, pseudoFill, prior);
		// TODO Auto-generated constructor stub
	}
	@Override
	public void initialiseMats(int numObsT){
		super.initialiseMats(numObsT);
			YB = new DenseDoubleMatrix2D(numObsT,1);
			covar = new DenseDoubleMatrix2D(numObsT,1);
			/*for(int i=0; i<this.b_fix_rows.length; i++){
				this.b_fix_rows[i] = new int[b_fix_len[i]];
			}
			this.non_skew_rows = new int[nonSkewRowCount.length][];
			for(int i=0; i<non_skew_rows.length; i++){
				non_skew_rows[i] = new int[nonSkewRowCount[i]];
			}*/
	}
	
	
	public void set(){
		super.set();
		if(Constants.regressMean(1)){
		for(int i=0; i<rows.length; i++){
		    if(paramBaf[i]!=null) setMidPoints(0, 1,paramBaf[i],this.illuminaProbRB.all_indices[i]);
		}
		}
	}
	
	public void setVar(){
		super.setVar();
		if(Constants.regressVariance(1)){
			for(int i=0; i<this.rows.length; i++){
				  if(paramBafVar[i]!=null) setMidPoints(1, 1,this.paramBafVar[i], this.illuminaProbRB.all_indices[i]);
			}
	     }if(false){
			for(int i=0; i<this.rows.length; i++){
				  if(paramCov[i]!=null)	  setMidPoints(1, 2,this.paramCov[i], this.illuminaProbRB.all_indices[i]);//, r_sigmaX, r_sigmaY);
					 
			}}
	}
	
	@Override
	public void regressMean(int st, double[][] priors, double pseudoM, boolean set){
		super.regressMean(st, priors, pseudoM, set);
	
		if( Constants.regressMean(1)){
			paramBaf = new RegressParams[rows.length];
			
		  regress("BAF", paramBaf, ((PriorRB)prior).priorBaf, YB, rows, illuminaProbRB.all_indices, priors[2], pseudoM, st);
	}
	}

	@Override
	public void regressVar(int st, double[][] priors, double pseudoM, boolean set){
		super.regressVar(st, priors, pseudoM, set);
			   if(Constants.regressVariance(1) ){
				   {
					   this.paramBafVar = new RegressParams[rows.length];
					
					   regress("BVAR", this.paramBafVar, ((PriorRB)prior).priorBafVar, YB, this.rows, 
							   illuminaProbRB.all_indices, priors[3], pseudoM, st);
				   }
				 if(false){//Constants.getMinValue(Constants.rho_train(0))<1e10 && ! Constants.orthogonal()){
					  paramCov = new RegressParams[rows.length];
					 regress("COVAR", this.paramCov, ((PriorRB)prior).priorRho, this.covar, this.rows, this.illuminaProbRB.all_indices, priors[4], pseudoM, st);
				   }
				  
			   }
		   
	}
	
@Override
	public void bst(String[] st) {
		 
	     if(paramBaf!=null ){
	    	 StringBuffer bst = new StringBuffer();
	       for(int i=0; i<paramBaf.length; i++){
	    	   if(paramBaf[i]!=null){
	    	   paramBaf[i].getFormulae(bst, illuminaProbRB.basis_mean);
	    	   bst.append("\n");
	    	   }
	       }
	       st[2] = bst.toString();
	     }
	      if(paramBafVar!=null){
	    	 StringBuffer bst = new StringBuffer();
		      
		       for(int i=0; i<paramBafVar.length; i++){
		    	   if(paramBafVar[i]!=null){
		    	   paramBafVar[i].getFormulae(bst, illuminaProbRB.basis_var);
		    	   bst.append("\n");
		    	   }
		       }
		       st[3] = bst.toString();
		     }

   }
	
}
