/**
 * 
 */
package lc1.dp.illumina;

import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class RegressParams{
	  final  DoubleMatrix2D Q;
	  final  DoubleMatrix2D x0;
	  // double[] pseudoR;
	  final  int noEl;
	  final  int[] cols;
	  final DoubleMatrix2D x;
	   final String type;
	   DoubleMatrix2D regressR;
	   
	 //  final DoubleMatrix2D min;
	 //  final DoubleMatrix2D max;
	   
	   
	   public void getFormulae(StringBuffer bst, BasisFunction bf) {
			// TODO Auto-generated method stub
		   String[] st1 = type.split(":");
		//	String[] str = st1[1].split("/");
			bst.append(st1[0]+"=");
			for(int k=len_extra; k<regressR.rows(); k++){
				double v= regressR.get(k, 0);
				int k1 = k-len_extra;
				if(Math.abs(v)>5e-3){
				bst.append((v<0 ? "-" : "+")+String.format("%5.3g",Math.abs(v))+"*"+bf.getName(k1));
						
				}
			}
		}
	   int len_extra;
	   /**note : the order of pseudoR puts category first, where as for x, Q etc, it is last */
	RegressParams(int noEl, int len_extra, double[] pseudo, double pseudoM, String type

			){
		int total_len = noEl;
		this.len_extra = len_extra;
		 x0 =  new DenseDoubleMatrix2D(noEl+len_extra,1);
		 cols = new int[noEl+len_extra];
		 this.type =  type;
		 double[][] q1 = new double[noEl+len_extra][noEl+len_extra];
		
		 this.noEl = noEl;
		// this.min = new DenseDoubleMatrix2D(1, total_len+len_extra);
		// this.max = new DenseDoubleMatrix2D(1, total_len+len_extra);
		  x = new DenseDoubleMatrix2D(1,total_len+len_extra);
		 
		 {
				 for(int i=0; i<noEl; i++){
					 cols[i] = i;
					 q1[i][i] =  Math.pow(pseudo[i+1]*pseudoM,2.0);
					 x0.setQuick(i, 0, 0);
				//	 this.min.setQuick(i,0,min[i+1]);
				//	 this.max.setQuick(i,0,max[i+1]);
				 }
				 for(int i=0; i<len_extra; i++){
		            	cols[noEl+i] = total_len+i;
		            	q1[i+noEl][ i+noEl] =  Math.pow(pseudo[0]*pseudoM,2.0);
		            	x0.setQuick(noEl+i, 0, 0);
		        //    	 this.min.setQuick(i,0,min[0]);
		    	//		 this.max.setQuick(i,0,max[0]);
		         }
		 
		 }
		 if(pseudoM==0){
			  q1[2][2] = 1;
			  q1[3][3] = 1;
			
				if(q1.length>5){
					q1[5][5] = 1;
					q1[6][6] = 1;
					
				}
			 }
		  Q = new DenseDoubleMatrix2D(q1 );
	//	  System.err.println("he");
	}
	   
	public void setUniformPrior(double prior) {
		for(int i=0; i<Q.size(); i++){
			Q.set(i, i, Math.pow(prior,2));
		}
		
	}
	
	
	
	public static RegressParams getRegressParams(String type, int noCat,  int xlen, double pseudoM, double[] pseudo){
		RegressParams res = new RegressParams(xlen,noCat, pseudo, 
		pseudoM,		
		type
	//(	probeOnly ? ":  const/LRR/LRR^2/LRR^3/cat" : ": const/LRR/LRR^2/LRR^3/BAF/BAF^2/BAF^3, cat")
		);
		
		//res.cols[2] =3;
		//res.cols[3] = 4;
		//res.x0.setQuick(3, 0,1.0);
		return res;
	}
	
	
   final private static int[] row_x = new int[]{0};
	   
	public double calculate() {
			DoubleMatrix2D x1 = x.viewSelection(row_x,cols);
			DoubleMatrix2D x2 = IlluminaProbR.lg.mult(x1,regressR);
		   double res = x2.getQuick(0, 0);
		///   if(Constants.CHECK && Double.isNaN(res)){
		//	   throw new RuntimeException("is na "+res+" "+this.regressR);
		//   }
		   return res;
	}
	
	public double calculate(DoubleMatrix2D mat) {
		DoubleMatrix2D x1 = x.viewSelection(row_x,cols);
		DoubleMatrix2D x2 = IlluminaProbR.lg.mult(x1,mat);
	   return x2.getQuick(0, 0);
}

	public void setX(int i, double d) {
		x.set(0, i, d);
		
	}

	
	public DoubleMatrix2D regress(int startPos, 
			  MultipleDataRegression data, 
			  DoubleMatrix2D Y, int[] rows
			  ){
		   int numObs=startPos;
		  
		 //  if(numObs!=numObsT) throw new RuntimeException("!!!");
		  
		  if(rows.length>0){
		   DoubleMatrix2D X1 = rows==null ? data.X : data.X.viewSelection(rows, cols);
		   DoubleMatrix2D Y1 = rows ==null || rows.length == Y.rows() ? Y : Y.viewSelection(rows, colsY);
		   this.regressR =   solve(X1,Y1,Q,x0); 
		  }
		  else{
			 if(regressR==null) regressR = new DenseDoubleMatrix2D(x0.rows(), x0.columns());
			  for(int i=0; i<this.x0.rows(); i++){
				  for(int j=0; j<this.x0.columns(); j++){
					  regressR.set(i, j, x0.getQuick(i, j));
				  }
			  }
		  }
		
		   return regressR;
	   }
	private static final int[] colsY = new int[] {0};
	public static DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b,
			   DoubleMatrix2D Q, DoubleMatrix2D x0) {
			
			   DoubleMatrix2D AT = IlluminaProbR.lg.transpose(A);
			   DoubleMatrix2D prod = IlluminaProbR.lg.mult(AT, A);
			   
			   for(int i=0; i<prod.rows(); i++){
				   for(int j=0; j<prod.columns(); j++){
					   prod.setQuick(i, j, prod.getQuick(i, j)+Q.getQuick(i, j));
				   }
			   }
			   DoubleMatrix2D b1 = b.copy();
			   DoubleMatrix2D Ax0 = IlluminaProbR.lg.mult(A, x0);
			   for(int i=0; i<b1.rows(); i++){
				   b1.setQuick(i, 0, b1.getQuick(i, 0)-Ax0.getQuick(i, 0));
			   }
			   DoubleMatrix2D res =  IlluminaProbR.lg.solve(prod,IlluminaProbR.lg.mult(AT,b1));
			   for(int i=0; i<res.rows(); i++){
				   double v1 = res.getQuick(i, 0);
				   try{
				   if(Constants.CHECK && Double.isNaN(v1)) {
					   validate(A);
					   validate(b);
					   validate(Q);
					   validate(x0);
					   throw new RuntimeException("regression did not work");
				   }
				   }catch(Exception exc){
					   exc.printStackTrace();
				   }
				   res.setQuick(i, 0, v1+x0.getQuick(i, 0));
			   }
			   
			   return res;
		 }
	
	public static DoubleMatrix2D solve(DoubleMatrix2D A, DoubleMatrix2D b
			  ) {
			
			   DoubleMatrix2D AT = IlluminaProbR.lg.transpose(A);
			   DoubleMatrix2D prod = IlluminaProbR.lg.mult(AT, A);
			   
			   
			
			   DoubleMatrix2D res =  IlluminaProbR.lg.solve(prod,IlluminaProbR.lg.mult(AT,b));
			 /*  for(int i=0; i<res.rows(); i++){
				   double v1 = res.getQuick(i, 0);
				   try{
				   if(Constants.CHECK && Double.isNaN(v1)) {
					   validate(A);
					   validate(b);
					   validate(Q);
					   validate(x0);
					   throw new RuntimeException("regression did not work");
				   }
				   }catch(Exception exc){
					   exc.printStackTrace();
				   }
				   res.setQuick(i, 0, v1+x0.getQuick(i, 0));
			   }*/
			   
			   return res;
		 }
private static void validate(DoubleMatrix2D A) {
	for(int i1=0; i1<A.rows(); i1++){
		   for(int j1=0; j1<A.columns(); j1++){
		   if(Double.isNaN(A.getQuick(i1, j1))) {
			   throw new RuntimeException("!!");
		   }
		   }
	   }
		
	}


public String toString(){
	
	
	return this.type+"\n"+this.regressR+"\n"+this.x0+"\n"+this.Q+"\n"+Constants.print(this.cols);
	
}


public void setConstantPrior(int pos, double d) {
	this.x0.set(pos,0,d);
	
}
public void setConstantsPriorWeight(double d){
	this.Q.set(0, 0, d);
}







}