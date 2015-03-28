/**
 * 
 */
package lc1.dp.transition;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.illumina.RegressParams;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Statistic;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;

class BackCalc implements MultivariateFunction{
	   /**
	 * 
	 */
	private final BetweenWithinTransitionProbs6 between;
	DoubleMatrix2D mat;
	 DoubleMatrix1D n;
	 //DoubleMatrix1D b;
	// SimpleExtendedDistribution[] dist;
	 double[] trans;
	 int[] cnts;
	 int nocols, norows, wtrows, transrows;
	 DoubleMatrix1D v;
	 DoubleMatrix2D nullspace;
	 static Algebra al = new Algebra();
	int numGroups;
	int[][] groupToState;
	 int numMembers;
	//int[] groupToState;
	 int[] grps;
	final int[] grps1;
	
	//final public DoubleMatrix1D bCount;
	final public DoubleMatrix2D vMCount;

	public double update(double[] trans, double[] count, double ent, double pseudo){
		double sumexcl = 0;
		
		this.entropyWeight = ent;
	this.pseudo = pseudo;
		for(int i=0; i<grps1.length; i++){
			sumexcl+=trans[grps1[i]];
		}
		/* keep this intact for checking!!
		 * for(int i1=0; i1<wtrows; i1++){ //wt rows first
			  b.setQuick(i1, this.between.frac(groupFrom).probs(i1)*(1-sumexcl));
			  if(count!=null)  bCount.setQuick(i1, this.between.frac(groupFrom).counts()[i1]*(1-sumexcl));
		  }
		  for(int ik=0; ik<grps.length; ik++){
			  int i2 = grps[ik];
			  int i1 =  ik ;
			  b.setQuick(wtrows+i1, trans[i2]);
			 if(count!=null) bCount.setQuick(wtrows+i1, count[i2]);
		  }*/
		  {
			  int k=0;
				 for(int ik=0; ik<grps.length; ik++){
					 int i= grps[ik];
					 double[] counts = this.between.frac(groupFrom).probs();
					// if(i==groupFrom) continue;
					 double sum = Constants.sum(counts);
					  for(int j=0; j<numMembers; j++){
						  double cnt_j = sum==0 ? 0 : counts[j]/sum;
						 vM.setQuick(k, 0, this.between.frac(groupFrom).probs(j)*trans[i]);
						 if(count!=null) vMCount.setQuick(k, 0, (cnt_j)*count[i]);
						  k++;
					  }
				  }
		  }
		 
		  v = vM.viewColumn(0);
		  
		  if(count!=null){
				DoubleMatrix2D nNew = this.getNCount();
				DoubleMatrix2D xNew  = this.solve(nullspace, this.subtract(nNew, vMCount));
			
				this.x =  xNew.viewColumn(0);
				double sum  = this.eval();
				
				double res;
				/*if(sum>1e-10) {
					res =  1e100*sum;
				}
				else{*/
					this.getFromN();
					 res =  -1.0*this.between.evaluteBack()+sum;
				//}*/
				return res;
				}
		  return 0;
	}
	 
	 DoubleMatrix2D vM;
	double pseudo;
	 public BackCalc(BetweenWithinTransitionProbs6 betweenWithinTransitionProbs5, int groupFrom, double[] trans, int[] grps, int[] grps1
		){
		
		 this.entropyWeight = 0;//entropyWeight;
		 this.probs = new double[trans.length];
		 this.between = betweenWithinTransitionProbs5;
		 this.groupToState = between.groupToState();
		this.trans = trans;
		this.grps = grps;
		this.grps1 = grps1;
		numGroups = groupToState.length;
		numMembers = groupToState[groupFrom].length;
		 this.groupFrom = groupFrom;
			int num = grps.length;
		cnts = new int[num];
		 nocols =grps.length*numMembers;
		wtrows = this.between.frac(groupFrom).probs.length;
		transrows = grps.length;
		 norows = wtrows+transrows;
		DoubleMatrix2D mat1 = new DenseDoubleMatrix2D(nocols,nocols);
		  n = new DenseDoubleMatrix1D(nocols); 
		 
		//  b = new DenseDoubleMatrix1D(norows);
		 // bCount = new DenseDoubleMatrix1D(norows);
//		  double max =1 - trans[groupFrom];
		  for(int i1=0; i1<wtrows; i1++){ //wt rows first
			  int k=0;
			 for(int ik=0; ik<grps.length; ik++){
				  for(int j=0; j<numMembers; j++){
					  if(j==i1) mat1.setQuick(i1, k, 1.0);
					  else mat1.setQuick(i1, k, 0.0);
					  k++;
				  }
			  }
		  }
		  for(int ik=0; ik<grps.length; ik++){
			  int i2 = grps[ik];
			  int i1 =  ik ;
			  int k=0;
				inner: for(int ik1=0; ik1<grps.length; ik1++){
					int i = grps[ik1];
					  for(int j=0; j<numMembers; j++){
						  if(i==i2) mat1.setQuick(wtrows+i1, k, 1.0);
						  else mat1.setQuick(wtrows+i1, k, 0.0);
						  k++;
					  }
				  }
		  }
		  vM = new DenseDoubleMatrix2D(nocols,1);   
		  vMCount = new DenseDoubleMatrix2D(nocols,1);   
		 this.update(trans,null, 0.0,1.0);
		  SingularValueDecomposition svd =   new SingularValueDecomposition(mat1);
		
		  double[] singularvals = svd.getSingularValues();
		  double thresh = singularvals[0]*Math.max(norows, nocols)*2.2e-16;
		  int rank =0;
		  for(int i=0; i<singularvals.length; i++){
			  if(singularvals[i]>thresh) rank++;
		  }
		int[] rows = new int[svd.getV().rows()];
		int[]cols = new int[nocols-rank];
		int[] rows1 = new int[norows];
		int[] cols1 = new int[nocols];
		for(int i=rank; i<nocols; i++){
			cols[i-rank] = i;
		}
		for(int i=0; i<rows.length; i++){
			rows[i] = i;
			if(i<norows)rows1[i] = i;
		}
		for(int i=0; i<nocols; i++){
			cols1[i] = i;
		}
		
		  mat = mat1.viewSelection(rows1,cols1);
		
			nullspace = svd.getV().viewSelection(rows, cols);
		this.x_init = new double[nullspace.columns()];
		List<Integer> colsL = new ArrayList<Integer>();
		for(int kk=0; kk<nullspace.size(); kk++){
			colsL.add(kk);
		}
		if(Constants.rotate()){
			List<DoubleMatrix1D> xv = new ArrayList<DoubleMatrix1D>();
		for(int kk=0; kk< this.grps.length; kk++){
			if(grps[kk]==this.groupFrom) continue;
			int ind = grps[kk]-1; //target index for most of transition to
		
			int[] cols2 = new int[this.numMembers-1];
			int[] rows2 = new int[cols2.length];
			/*{
				double trans_max = 0;
			for(int kj=0; kj<cols2.length; kj++){
				double trans_ = this.v.get(kk*this.numMembers+kj);
				if(trans_>trans_max){
					ind = kj;
					trans_max = trans_;
				}
			}
			}*/
			//	
			
			int kj1 =0;
			for(int kj=0; kj<cols2.length; kj++){
				cols2[kj] = colsL.remove(0);
				if(kj1==ind) kj1++;
				rows2[kj] = kk*this.numMembers + kj1;
				kj1++;
			}
			
			
		
			DoubleMatrix2D nullspace1 = nullspace.viewSelection(rows2, cols2);
		//	EigenvalueDecomposition evd = new EigenvalueDecomposition(nullspace1);
	        x = 
	        	this.mult(this.solve(nullspace1,vM.viewSelection(rows2, new int[]{0})).viewColumn(0),-1);
	        xv.add((DoubleMatrix1D)x.clone());
		}
		  DoubleMatrix2D mat_1 = getGramSchmidt(xv);
	        this.nullspace = al.mult(this.nullspace, mat_1);
		}
		x = new DenseDoubleMatrix1D(nullspace.columns());
		for(int k=0; k<x.size(); k++){
			x.set(k, 1.0);
			int min_ind = addMin(al.mult(nullspace, x),v);
			double min_scale =v.get(min_ind)/ nullspace.get(min_ind,k);
			for(int j=0; j<nullspace.rows(); j++){
				nullspace.set(j,k, nullspace.get(j,k)*min_scale);
			}
			x.set(k,0.0);
		}
	 }
	 private DoubleMatrix2D getGramSchmidt(List<DoubleMatrix1D> x2) {
		 int len = x2.get(0).size();
		DoubleMatrix2D inp = new DenseDoubleMatrix2D(len, len);
		for(int i=0; i<len; i++){
			inp.setQuick(i, i, 1);
		}
		DoubleMatrix2D inp1 = this.getNewNullspace(inp, x2); //first column is target
		for(int i=1; i<len; i++){
			DoubleMatrix1D xi = inp1.viewColumn(i);
			double[] weight = new double[i];
			for(int j=0; j<i; j++){
				DoubleMatrix1D vj = inp1.viewColumn(j);
				weight[j] = al.mult(xi, vj) / al.mult(vj, vj);
				
			}
			for(int j=0; j<i; j++){
				this.subtract(xi, inp1.viewColumn(j), weight[j]);
			}
		}
		SingularValueDecomposition svd = new SingularValueDecomposition(inp1);
	      if(svd.rank()!=len) {
	    	  throw new RuntimeException("!!");
	      }
	 //	 EigenvalueDecomposition eig1= new EigenvalueDecomposition(nullspace2);*/
		return inp1;
	}
	private void subtract(DoubleMatrix1D vi, DoubleMatrix1D vj, double w) {
	for(int i=0; i<vi.size(); i++){
		vi.setQuick(i, vi.getQuick(i) - vj.getQuick(i)*w);
	}
		
	}
	private DoubleMatrix2D getNewNullspace(DoubleMatrix2D nullspace2,
			List<DoubleMatrix1D> xx) {
		// EigenvalueDecomposition eig = new EigenvalueDecomposition(nullspace2);
		// TODO Auto-generated method stub
		 int cols1 = nullspace2.columns();
		 int rows1 = nullspace2.rows();
		 DoubleMatrix2D mat1 = new DenseDoubleMatrix2D(rows1,cols1+xx.size());
		 for(int i=0; i<rows1; i++){
			 
			 for(int j=0; j<cols1; j++){
				 mat1.setQuick(i, j, nullspace2.getQuick(i, j));
			 }
			 for(int j=0; j<xx.size(); j++){
				 mat1.setQuick(i, cols1+j, xx.get(j).get(i));
			 }
		 }
			DoubleMatrix2D matrix1 = Statistic.covariance(mat1);
			Statistic.correlation(matrix1);
	      int[] max_ind =new int[xx.size()];
	      double[] vmax =new double[xx.size()];
	      for(int i=0; i<cols1; i++){	
	    	  boolean allocated = false;
	    	  for(int j=0; j<xx.size(); j++){
		    	  double v = Math.abs(matrix1.get(cols1+j, i));
		    	  if(!allocated && v>vmax[j]){
		    		  allocated = true;
		    		  vmax[j] = v;
		    		  max_ind[j] = i;
		    	  }
	    	  }
	      }
	      for(int i=0; i<rows1; i++){
	    	  for(int j=0; j<xx.size(); j++){
		    	  nullspace2.setQuick(i,max_ind[j], nullspace2.getQuick(i, 0+j));
		    	  nullspace2.setQuick(i, 0+j, xx.get(j).get(i));
	    	  }
	      }
	    /*  SingularValueDecomposition svd = new SingularValueDecomposition(nullspace2);
	      if(svd.rank()!=cols1) {
	    	  throw new RuntimeException("!!");
	      }
	 //	 EigenvalueDecomposition eig1= new EigenvalueDecomposition(nullspace2);*/
		return nullspace2;
	}
	private DoubleMatrix1D mult(DoubleMatrix1D viewColumn, double mult) {
			  for(int i=0; i<viewColumn.size(); i++){
				  viewColumn.set(i, mult*viewColumn.get(i));
			  }
			  return viewColumn;
	}
	private int[] getSmallest(double[] probs, int n) {
		int[] order = Constants.getOrder(probs);
		int[] res = new int[n];
		System.arraycopy(order, 0, res, 0, res.length);
		return res;
	}
	private double[] getNew(double[] prev, int within_ind, double valn) {
		double[] res = new double[prev.length];
		double diff = prev[within_ind]-valn;
		double sum = Constants.sum(prev) - prev[within_ind];
		for(int k=0; k<res.length; k++){
			if(k==within_ind) res[k] = valn;
			else{
				res[k] = prev[k] + (prev[k]/sum)*diff;
			}
		}
		return res;
	}
	DoubleMatrix2D v4 ;
	 private double[] get(DoubleMatrix2D v4, int i) {
		double[] res = new double[v4.rows()];
		for(int j=0; j<res.length; j++){
			res[j] = v4.get(j, i);
		}
		return res;
	}
	private DoubleMatrix2D solve(DoubleMatrix2D nullspace2, DoubleMatrix2D v3) {
		return RegressParams.solve(nullspace2, v3);
	}
	private DoubleMatrix2D subtract_all(DoubleMatrix2D v2, DoubleMatrix1D v3) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(v2.rows(),v2.columns());
		for(int i=0; i<res.rows(); i++){
			for(int j=0; j<res.columns(); j++){
				res.setQuick(i, j, v2.getQuick(i, j) - v3.getQuick(i));
			}
		}
		return res;
	}
	private double[] getBound(DoubleMatrix1D v2, DoubleMatrix1D viewColumn) {
		 double low = -100;
			double high = 100;
			double[] x_ = new double[v.size()];
			for(int i=0; i<v.size();i++){
				double lb = -1.0*v.get(i)/nullspace.getQuick(i, 0);
				x_[i] = lb;
				if(lb>0){
					//if(max>0.01){
					if(lb<high){
						high = lb;
					}
					//}
				}
				else if(lb<0){
					if(lb>low){
						low = lb;
					}
				}
				//double up = b.get(i) - 
			}
			return new double[] {low,high};
	}
	 /*private double getTheta(DoubleMatrix 1D Nrow, double b){
		 
	 }*/
	public static DoubleMatrix2D getColSpace(DoubleMatrix2D mat, int[] rows) {
		 int nocols = mat.columns();
		 int norows = mat.rows();
		
		  SingularValueDecomposition svd = 
			
			  new SingularValueDecomposition(mat);
		  
		  double[] singularvals = svd.getSingularValues();
		  double thresh = singularvals[0]*Math.max(norows, nocols)*2.2e-16;
		  int rank =0;
		  for(int i=0; i<singularvals.length; i++){
			  if(singularvals[i]>thresh) rank++;
		  }
		  int[] cols = new int[rank];
		  for(int i=0; i<rank; i++){
			  cols[i] = i;
		  }
		  return svd.getU().viewSelection(rows, cols);
		  
	}
	public static DoubleMatrix2D getIdentity(int i) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(i,i);
		for(int k=0; k<i;k++){
			res.setQuick(k, k, 1.0);
		}
		return res;
	}
	double[] x_init;
	 
	final double[] probs;

	public void getFromN(){
		 int k=0; 
		// between.validate1();
		
			boolean fix = false;
			// if(groupTo==groupFrom) continue;
			 for(int j=0; j<numMembers; j++){
				 int indexFrom = j;
				 Arrays.fill(probs, 0);
				 for(int i1=0; i1<grps.length; i1++){
					 k = i1*numMembers + j;
					 int groupTo = grps[i1];
					 
					 double prob =  n.getQuick(k)/between.frac(groupFrom).probs[indexFrom];
					 if(prob<0){
						 fix = true;
						 prob =0;
					 }
				//	 Math.max(0,);
					 probs[groupTo] = prob;
				//	 between.setProb(groupFrom, groupTo, indexFrom,prob );
				//	 between.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom].setProbs(groupTo,prob);
					 
				 }
				((SimpleExtendedDistribution) between.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom]).transfer(probs, pseudo);
		 }
			 
		if(fix) between.validate1(false);
	 }
	
	
	 
	public DoubleMatrix2D getNCount(){
		 DoubleMatrix2D res = new DenseDoubleMatrix2D(n.size(),1);
		 int k=0; 
		// between.validate1();
		
			
			// if(groupTo==groupFrom) continue;
			 for(int j=0; j<numMembers; j++){
				 int indexFrom = j;
				 Arrays.fill(probs, 0);
				 double sum=0;
				 for(int i1=0; i1<grps.length; i1++){
					 k = i1*numMembers + j;
					 int groupTo = grps[i1];
					 double prob =   between.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom].counts()[groupTo];
					 probs[groupTo] = prob* between.frac(groupFrom).probs[indexFrom];
					 sum+=prob;
				 }
				 for(int i1=0; i1<grps.length; i1++){
					 k = i1*numMembers + j;
					 int groupTo = grps[i1];
					 res.setQuick(0,k, sum==0 ? 0 : probs[groupTo]/sum);
				 }
				 /*
				  *  for(int i1=0; i1<grps.length; i1++){
			 int groupTo = grps[i1];
			// if(groupTo==groupFrom) continue;
			 for(int j=0; j<numMembers; j++){
				 int indexFrom = j;
				 double prob =between.freeTransBetweenGroups[groupFrom].transitionsOut[indexFrom].counts()[groupTo];
					 
					 //between.getProb(groupFrom, groupTo, indexFrom );
				 res.setQuick( 0, k,prob * between.frac(groupFrom).probs[indexFrom]);
				// double prob =  n.getQuick(k)/between.frac(groupFrom).probs[indexFrom];
//				 between.setProb(groupFrom, groupTo, indexFrom,prob );
				
				 k++;
			 }
		 }
				  */
		 }
			return res;
		 //between.validate1(false);
	 }
	
	
	 public double getEntropy(){
		
		 int k=0; 
		 double sum=0;
		 for(int i1=0; i1<grps.length; i1++){
			 int groupTo = grps[i1];
			// if(groupTo==groupFrom) continue;
			 for(int j=0; j<numMembers; j++){
				 int indexFrom = j;
				 double prob = (between.getProb(groupFrom, groupTo, indexFrom ) 
			  *between.frac(groupFrom).probs[indexFrom])
				 / between.getGroupProb(groupFrom, groupTo);
				 if(prob>0){
				sum+=-1*prob*Math.log(prob);
			/*	if(Double.isNaN(sum)){
					throw new RuntimeException("!!");
				}*/}
				 k++;
			 }
		 }
		 return sum;
		 
	 }
	 
	 public final int groupFrom;
	 
	/* public double[]  getN(){
		 int k=0; 
		 for(int i=1; i<dist.length; i++){
			 if(i==groupFrom) continue;
			 for(int j=0; j<dist[i].probs.length; j++){
				 n.setQuick(k, dist[i].probs(j)*trans[i]);
				 k++;
			 }
		 }
		 double[] res = new double[n.size()];
		 for(int i=0; i<res.length; i++){
			 res[i] = n.getQuick(i);
		 }
		 return res;
	 }*/
DoubleMatrix1D x;

double entropyWeight;	




@Override
	public double evaluate(double[] argument) {
		for(int i=0; i<argument.length; i++){
			x.setQuick(i, argument[i]);
		}
		
		double sum  = this.eval();
		/*keep intact for checking
		 * if(Constants.CHECK){
		double diff = subtract(al.mult(mat, n),b);
			if(diff>0.001) throw new RuntimeException("!! not in nullspace");
		}*/
		double res;
		if(sum>1e-10) {
			res =  1e100*sum;
		}
		else{
			this.getFromN();
			 res =  -1.0*this.between.evaluteBack()+sum;
		}
		double entropy = entropyWeight *this.getEntropy();
		//System.err.println("evals "+Constants.print(argument)+" "+(res+entropy)+"\n" +n);
		return res+entropy;
	}

public double eval(){
	DoubleMatrix1D xx = //argument.length==v4.columns() ? 
		al.mult(nullspace, x)
	//:	al.mult( nullspace,x)
		;
double sum =0;
for(int i=0; i<xx.size(); i++){
	double val = xx.get(i)+v.get(i);
	n.setQuick(i, val);
	if(val<0){
		sum+=Math.pow(val, 2);
	}
}
return sum;
}


	


	private static double subtract(DoubleMatrix1D mult, DoubleMatrix1D b2) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(mult.size(),1);
		double sum=0;
		for(int i=0;i<mult.size(); i++){
			
			sum+=Math.pow(mult.getQuick(i)-b2.getQuick(i),2);
		}
	//	System.err.println("h "+sum);
		return sum;
	}
	
	private static int addMin(DoubleMatrix1D mult, DoubleMatrix1D b2) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(mult.size(),1);
		double min = Double.POSITIVE_INFINITY;
		int min_ind =-1;
		for(int i=0;i<mult.size(); i++){
				
			double v = (mult.getQuick(i)+b2.getQuick(i));
			if(v <min){
				min = v;
				min_ind = i;
			}
		}
	//	System.err.println("h "+sum);
		return min_ind;
	}
	
public  static DoubleMatrix2D subtract(DoubleMatrix2D mult, DoubleMatrix2D b2) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(mult.rows(), mult.columns());
		
		for(int i=0;i<mult.rows(); i++){
			for(int j=0;j<mult.columns(); j++){
			res.setQuick(i,j,mult.getQuick(i,j)-b2.getQuick(i,j));
			}
			
		}
	//	System.err.println("h "+sum);
		return res;
	}
	private DoubleMatrix2D add(DoubleMatrix2D mult, DoubleMatrix1D b2) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(mult.size(),1);
		for(int i=0;i<mult.size(); i++){
			res.setQuick(i,0, mult.getQuick(i,0)+b2.getQuick(i));
		}
		return res;
	}
	private DoubleMatrix1D add(DoubleMatrix1D mult, DoubleMatrix1D b2) {
		DoubleMatrix1D res = new DenseDoubleMatrix1D(mult.size());
		for(int i=0;i<mult.size(); i++){
			res.setQuick(i,mult.getQuick(i)+b2.getQuick(i));
		}
		return res;
	}


	@Override
	public double getLowerBound(int i) {
		this.eval();
		double max = -Double.POSITIVE_INFINITY;
		
		for(int k=0; k<n.size(); k++){
			double Ni = nullspace.get(k, i);
			if(Ni>0){
			double v = -(n.get(k) - x.get(i)*Ni)/Ni;
			if( v>max){
				max = v;
			}
			}
		}
		
		return max*0.999;
	}


	@Override
	public int getNumArguments() {
		return this.nullspace.columns();
	}


	@Override
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public double getUpperBound(int i) {
		// TODO Auto-generated method stub
		this.eval();
		double min = Double.POSITIVE_INFINITY;
		
		for(int k=0; k<n.size(); k++){
			double Ni = nullspace.get(k, i);
			if(Ni<0){
			double v = -(n.get(k) - x.get(i)*Ni)/Ni;
			if( v<min){
				min = v;
			}
			}
		}
		
		return min*0.999;
	}
	//final double[] x_high, x_low;
	public void initialise(double[] xvec) {
		int len = xvec.length;
		for(int i=0; i<len; i++){
			double frac = Constants.rotate() ? 0 :
				1/(1+Math.exp(normal.quantile(Math.random(), 0, 1)));
			double v = this.getLowerBound(i)+ frac*(this.getUpperBound(i) - this.getLowerBound(i));
			this.x.set(i, v);
			xvec[i] = v;
		}
		this.eval();
		this.getFromN();
	}
	
	static lc1.stats.NormalDistribution normal;
   }