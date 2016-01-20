/**
 * 
 */
package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Arrays;
import java.util.logging.Logger;

import lc1.stats.SimpleDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;

public class MatrixExp  implements Serializable{
	
	public static void main(String[] args){
		try{
			double[][] d = new double[][] {new double[] {0.8, 0.5}, new double[] {-0.1, 1.0}};
			DoubleMatrix2D rates = new DenseDoubleMatrix2D(d);
			EigenvalueDecomposition evd = new EigenvalueDecomposition(rates);
			DoubleMatrix2D U = evd.getV();
			DoubleMatrix2D Uinv = alg.inverse(U);
			DoubleMatrix1D eig = evd.getRealEigenvalues();
			DoubleMatrix1D 	eigIm = evd.getImagEigenvalues();
			DoubleMatrix2D D = new SparseDoubleMatrix2D(U.rows(), U.rows());
			double distance = 5;
			for(int i=0; i<D.rows(); i++){
				double mui = eig.get(i);
				double expmui = Math.exp(mui*distance);
				D.setQuick(i, i,Math.exp(mui*distance));
			}
			DoubleMatrix2D M = alg.mult(U,alg.mult(D, Uinv));
			//MatrixExp mat = new MatrixExp(d,1.0);
			//mat.setDistance(5.0);
			System.err.println(M);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
		public DoubleMatrix2D rates;
		public  Double currentRate ;//= Constants.expModelIntHotSpot[0]*Constants.probCrossOverBetweenBP;
		DoubleMatrix2D U;
		DoubleMatrix2D Uinv;
		DoubleMatrix2D D;
		public DoubleMatrix2D M;
		DoubleMatrix1D eig;
		public DoubleMatrix2D pi;
		DoubleMatrix2D J;
	//	DoubleMatrix2D Iab;
		double distance;
		final public int len; 
		double[][]countE;
		boolean changed = false; //if rates, J, UInv have changed and need to be recalculated
		
		
		
		public MatrixExp(double [] pi, double rate){
			this(getBaseMatrix(pi.length), pi, rate);
			//this.getNullSpace();
		}
		
		public double getRate(int i, int j) {
			return rates.getQuick(i, j);
		}
		
		
		public double getM(int i, int j) {
			return M.getQuick(i, j);
		}
	
	
		
		public MatrixExp(MatrixExp exp){
			this.currentRate = new Double(exp.currentRate);
			this.rates = exp.rates.copy();
			this.pi = exp.pi==null ? null: (DoubleMatrix2D) exp.pi.clone();
			this.U = exp.U.copy();
			this.D = exp.D.copy();
			this.Uinv = exp.Uinv.copy();
			this.eig = exp.eig.copy();
			changed = true;
			//this.update(null);
			this.len = exp.len;
			countE = new double[len][len];
			this.J = new DenseDoubleMatrix2D(len,len);
		//	changed = true;
		}
		public MatrixExp clone(){
			return new MatrixExp(this);
		}
		
		/**should take advantage of symetric matrix here?*/
		public MatrixExp(DoubleMatrix2D ratesS, double[] pi1, double rate) {
			double[] pi = new double[pi1.length];
			System.arraycopy(pi1, 0, pi, 0, pi.length);
			boolean zeroPi = false;
			for(int k=0; k<pi.length; k++){
				if(pi[k]==0){
					zeroPi = true;
					pi[k] = 1e-10;
				}
			}
			this.currentRate  = rate;
			double[] pi_bg = new double[pi.length];
			this.len = pi.length;
			countE = new double[len][len];
	    	Arrays.fill(pi_bg, 1.0/(double)pi.length);
			this.pi = new DenseDoubleMatrix2D(1,pi.length);
			for(int i=0; i<pi.length; i++){
				this.pi.setQuick(0, i, pi[i]);
			}
			this.rates = getMod(0.5, ratesS, pi, pi_bg);
			
			//multiplyRate(rates, rate);
			this.J = new DenseDoubleMatrix2D(len,len);
			this.makeValid();
			this.calcPiRate();
			this.resetRate(rate);
			this.update();
			if(zeroPi){
				for(int to=0; to<pi.length; to++){
					if(pi1[to]==0){
						for(int i=0; i<pi.length; i++){
							double diff = 0;
							if(i!=to) {
							    diff+=rates.getQuick(i,to);
								rates.set(i, to, 0);
							}
							rates.set(i, i, rates.get(i,i)+diff);
						}
					}
				}
			}
			System.err.println("here");
		//	changed = true;
		
			
		}
		
		public MatrixExp(double[][]ma, double mod){
			//this.currentRate = rate;
			this.len = ma.length;
			countE = new double[len][len];
			this.rates = new DenseDoubleMatrix2D(ma);
		    this.multiplyRate(rates, mod);
			this.J = new DenseDoubleMatrix2D(len,len);
			this.makeValid();
			this.calcPiRate();
			this.update();
		}
		
		public void updateRates(double[][] mat, double[] pi) {
			for(int i=0; i<pi.length; i++){
				pi[i] = this.pi.getQuick(i, 0);
			}
			double max =0;
			for(int i=0; i<mat.length; i++){
				for(int j=0; j<mat[i].length; j++){
					double res = this.rates.getQuick(j, i);
					mat[i][j]  = res;
					if(Math.abs(res) > max){
						max = res;
					}
				}
			}
			/*for(int i=0; i<mat.length; i++){
				for(int j=0; j<mat[i].length; j++){
					
					mat[i][j]  = mat[i][j]/max;
					
				}
			}*/
			
			
		}
		/** rate is target average rate */
		public void update(){
			EigenvalueDecomposition evd = new EigenvalueDecomposition(rates);
			U = evd.getV();
			Uinv = alg.inverse(U);
			this.eig = evd.getRealEigenvalues();
			DoubleMatrix1D 	eigIm = evd.getImagEigenvalues();
			if(nonZero(evd.getImagEigenvalues())){
				System.err.println("error");
				//throw new RuntimeException("not expecting imag eigv");
			}
			this.D = new SparseDoubleMatrix2D(U.rows(), U.rows());
			changed = true;
		}
		
		private boolean nonZero(DoubleMatrix1D imagEigenvalues) {
			for(int i=0; i<imagEigenvalues.size(); i++){
				if(Math.abs(imagEigenvalues.get(i))>0){
					return true;
				}
			}
			return false;
		}

		
		
		public void calcPiRate(){
			this.pi = this.getNullSpace();
			this.currentRate = calcRate();
			
		}
		
		public void resetRate(double rate){
		//Logger.global.info("old rate "+currentRate);
				multiplyRate(rates, currentRate==0 ? 0 : rate/this.currentRate);
				this.currentRate = rate;
				this.changed = true;
	//			Logger.global.info("new rate "+currentRate);
			
		}
		
		public double calcRate(){
			double sum=0;
			for(int i=0; i<this.rates.rows(); i++){
				sum+= - this.pi.get(i, 0)*this.rates.getQuick(i, i);
			}
			return sum;
		}
		
		static Algebra alg = new Algebra();
		
		public void addCount(int a, int b, double v, double distance){
			//if(distance < 1e-9  || v<=0) return;
			this.setDistance(distance);
			double trans = this.M.getQuick(a, b);
			if(trans ==0 && v>0) {
				Logger.global.warning("trans prob was zero "+distance+" "+v);
				return ;
			}
			double abv  =v==0 ? 0 : v *(1/trans);
			/*if(Double.isNaN(abv)){
				System.err.println("problem\n"+M.toString());
				System.err.println(v);
				System.exit(0);
			}*/
			for(int k=0; k<len; k++){
				for(int l=0; l<len; l++){
					double toA = abv * this.J.getQuick(k, l) * this.U.getQuick(a, k) * this.Uinv.getQuick(l, b);
					countE[k][l]+= toA;
				//	if(Double.isNaN(countE[k][l])){
					//	throw new RuntimeException("!!"+J+"\n"+U+"\n"+Uinv+"\t"+toA+"\t"+abv+"\t"+distance);
				//	}
				}
			}
			
		}
		public void validate(DoubleMatrix2D matr){
			for(int i=0; i<matr.rows(); i++){
				for(int j=0; j<matr.columns(); j++){
					double d = matr.getQuick(i, j);
					if(Double.isNaN(d)){
						System.err.println("prob with\n"+matr);
					}
				}
			}
		}
		public void initialise(){
			for(int i=0; i<this.countE.length; i++){
				Arrays.fill(countE[i],0.0);
			}
		}
	
		public void transferInner(){
			
			
			double[][]counts = new double[len][len];
			double[] w = new double[len];
			for(int i=0; i<len; i++){
				for(int j=0; j<len; j++){
					double v = this.rates.getQuick(i, j);
					for(int k=0; k<len; k++){
						double v1 = v * Uinv.getQuick(k, i);
						for(int l=0; l<len; l++){
							counts[i][j] += v1 * U.getQuick(j, l) * countE[k][l];
						}
					}
				}
				 w[i] = counts[i][i]/ this.rates.getQuick(i, i);
				
			}
			for(int i=0; i<len; i++){
				for(int j=0; j<len; j++){
					if(i!=j){
						double r = counts[i][j]/w[i];
						if(r<0){
							Logger.global.info("problem with rates");
							r=1e-10;//
							//throw new RuntimeException("!!");
						}
						this.rates.setQuick(i, j, r);
					}
			}
			}
		}
		
		public void addPseudoCounts(double pseudo, MatrixExp expected, double distancePseudo){
			if(distancePseudo>1e-7){
				expected.setDistance(distancePseudo);
				for(int i=0; i<len; i++){
					for(int j=0; j<len; j++){
						this.addCount(i, j, expected.M.getQuick(i, j)*pseudo, distancePseudo);
					}
				}
			}
			
		}
		
		
		
		public void transfer(double pseudo,  double pseudoAlpha, double pseudoRate, MatrixExp expected, double distancePseudo){
			if(pseudo>0) addPseudoCounts(pseudo, expected, distancePseudo);
			this.transferInner();
			double rate = this.currentRate;
	//	this.tr
		try{
			this.makeValid();
			this.calcPiRate();
		//	Logger.global.info("rate is "+this.currentRate+"\n"+this.pi);
			if(pseudoRate > 1e-3){ //rate is fixed
			//	if(true) throw new RuntimeException("!!");
				 this.resetRate(rate);
				
				if(pseudoRate < 1e5 && false){
					this.update();
					
					addPseudoCounts(-1*pseudo, expected, distancePseudo);
					addPseudoCounts(pseudoRate, this, distancePseudo);
					this.transferInner();
					this.makeValid();
					this.calcPiRate();
					Logger.global.info("NEW RATE "+this.currentRate);
				}
			}
			 if( pseudoAlpha>1e-3 && false) { //alpha is fixed
				// if(true) throw new RuntimeException("!!");
				 double targetRate = this.currentRate;
				 this.rates = expected.rates.copy();
				 this.calcPiRate();
				 this.resetRate(targetRate);
			//	 if(Math.abs(distancePseudo-2719)<1.0){
				//	 System.err.println("TARGET "+targetRate+" "+this.currentRate);
			//	 }
				// if(targetRate > expected.currentRate){
				//	 System.err.println("TARGET "+targetRate+" "+this.currentRate);
				// }
				 if(pseudoAlpha < 1e5 && false){
					    this.update();
					    addPseudoCounts(-1*pseudo, expected, distancePseudo);
					    addPseudoCounts(pseudoAlpha, this, distancePseudo);
						this.transferInner();
						this.makeValid();
						this.calcPiRate();
						
				 }
			}
		//	Logger.global.info("new global rate is "+this.currentRate+"\n"+this.pi);
			
			this.update();//expected.currentRate);
		}catch(Exception exc){
			exc.printStackTrace();
			System.exit(0);
		}
	}
		public double currentRate(){
			return this.currentRate;
		}
	
	public String toString(){
		DoubleMatrix2D rates1 = new DenseDoubleMatrix2D(len, len);
		for(int i=0; i<len; i++){
			for(int j=0; j<len; j++){
				rates1.setQuick(i, j, rates.getQuick(i, j)/this.currentRate);
			}
		}
	//	DoubleMatrix2D nulls = this.getNullSpace();
		return rates1.toString()+"__"+this.pi;
	}
	
	public DoubleMatrix2D getNullSpace() {
		/*double d = this.distance;
		this.setDistance(1e10, true);
		DoubleMatrix2D ones =  ones(M.rows());
		DoubleMatrix2D res = this.alg.mult(alg.transpose(ones),M);
		this.setDistance(d, true);
	//	return res;
	*/
		//this.multiplyRate(rates, 1e10);
		 int nocols = rates.columns();
		 int norows = rates.rows();
		 int[] rows = new int[rates.rows()];
		 for(int i=0; i<rows.length; i++){
			 rows[i] = i;
		 }
		
		  SingularValueDecomposition svd =   new SingularValueDecomposition(alg.transpose(rates));
			
			
		  
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
		/*  DoubleMatrix2D u1 = svd.getV().viewSelection(rows, cols);
		  DoubleMatrix2D	nullspace=BackCalc.getColSpace(BackCalc.subtract(BackCalc.getIdentity(rows.length),alg.mult(u1, alg.transpose(u1))),rows);
		*/
		 DoubleMatrix2D res =  svd.getV().viewSelection(rows, new int[] {rank}).copy();
		 double sum = res.zSum();
		 for(int i=0; i<res.rows(); i++){
			 res.setQuick(i, 0, res.getQuick(i, 0)/sum);
		 }
		 sum = res.zSum();
		 return res;
		// for(int i=0; i<res1.size(); i++){
		//	 res1.set(i, 0, pi.get(0, i));
		 //}
	//	 DoubleMatrix2D res3 = alg.mult(alg.transpose(res1),rates);
		 //return res1;
		
	}
		
		private DoubleMatrix2D ones(int rows) {
		DoubleMatrix2D res = new DenseDoubleMatrix2D(rows,1);
		for(int i=0; i<rows; i++){
			res.setQuick(i, 0, 1.0/(double)rows);
		}
		return res;
	}
		public void setDistance(double distance){
			if(Math.abs(distance-this.distance)<2 && ! changed) return;
			System.err.println("setting distance "+distance+" vs "+this.distance);
			this.distance = distance;
			for(int i=0; i<D.rows(); i++){
				double mui = eig.get(i);
				double expmui = Math.exp(mui*distance);
				D.setQuick(i, i,Math.exp(mui*distance));
				for(int j=0; j<len; j++){
					double muj = eig.get(j);
					if(Math.abs(mui-muj)<1e-10){
						J.setQuick(i, j, distance * expmui);
					}
					else{
						J.setQuick(i, j, (expmui-Math.exp(muj*distance))/(mui-muj));
					}
				}
			}
	
			this.M = alg.mult(U,alg.mult(D, Uinv));
			
			/*next bit is a bit of a hack due to errors in exponentiating - if matrix has low rate, then prob should be set to exactly zero 
			for(int i=0; i<len; i++){
			double resid = 0;
			for(int j=0; j<len; j++){
				if(this.rates.get(i, j)<1e-10){
					resid+=M.get(i,j);
					M.set(i, j, 0);
				}
			}
			M.set(i, i, M.get(i, i)+resid);
		}
			*/
			
			changed = false;
			fixM();
				if(true || Constants.CHECK ) {
					
					this.validate();
					this.validate(M);
					
				}
		    }
		public void fixM(){
			for(int i=0; i<M.rows(); i++){
				double s = 0;
				
				for(int j=0; j<M.columns();j++){
					double p = M.getQuick(i,j);
					if(p<0){
						s+=-p;
						M.setQuick(i,j,0);
					}
					
				}
				if(s>0){
					
					int max_id = 0;
					for(int j=1; j<M.columns();j++){
						if(M.getQuick(i,j)>M.getQuick(i,max_id)){
							max_id = j;
						}
					}
					M.setQuick(i,max_id,M.getQuick(i,max_id)-s);
					if(s>1e-7){
						Logger.getAnonymousLogger().warning("pathcing M");
					}
				}
			}
		}
		
		private boolean allzero(DoubleMatrix2D j2) {
			for(int i=0; i<j2.rows(); i++){
				for(int j=0; j<j2.rows(); j++){
					if(Math.abs(j2.getQuick(i, j))>1e-9){
						return false;
					}
				}
			}
			return true;
		}

		public static void multiplyRate(DoubleMatrix2D rates, double rate){
			int len = rates.rows();
			for(int i=0; i<len; i++){
	    		for(int j=0; j<len; j++){
	    			rates.set(i,j, rates.getQuick(i, j)*rate);
	    		}
	    		
			}
		}
		 public static DoubleMatrix2D getBaseMatrix(int len){
		    	DoubleMatrix2D rate = new DenseDoubleMatrix2D(len, len);
		    	for(int i=0; i<len; i++){
		    		for(int j=0; j<len; j++){
		    			if(j!=i){
		    				rate.set(i, j,(1.0/ ((double) len-1)) );//*Math.pow(pi[j], 1) );
		    			}
		    			else{
		    				rate.set(i, i, -1.0);
		    			}
		    		}
		    	}
		    	
		    	return rate;
		    }
		 protected static DoubleMatrix2D  getMod(double F, DoubleMatrix2D dbrates, double[] frequency, double[] dbfreq) {
		        //  java.util.logging.Logger.global.info("QR "+this.getClass());
		            //System.exit(0);
		  
		         int len = dbrates.rows();
		     	DoubleMatrix2D rates = new DenseDoubleMatrix2D(len, len);
		            for (int i = 0; i < len; i++)
		                {
		                    for (int j = i + 1; j < len; j++)
		                        {
		                            double mod = StrictMath.pow(frequency[j]/dbfreq[j],1- F)*
		                                                    StrictMath.pow(dbfreq[i]/frequency[i], F);
		                            rates.setQuick(i, j, dbrates.getQuick(i, j)*mod);
		                            double mod_1 = StrictMath.pow(frequency[i]/dbfreq[i],1- F)*
		                                                            StrictMath.pow(dbfreq[j]/frequency[j], F);
		                            rates.setQuick(j, i, dbrates.getQuick(i, j)*mod_1);
		                        }
		                }
		                   
		                    return rates;
		                
		                   
		        }
		public void validate() {
			for(int i=0; i<rates.rows(); i++){
				double sum=0;
				double sum1 = 0;
				
				for(int j=0; j<rates.rows(); j++){
					sum1+=M.getQuick(i, j);
					if(j!=i) sum+=rates.getQuick(i,j);
				}
				if(Math.abs(rates.getQuick(i, i) +sum)>SimpleDistribution.tolerance) {
					throw new RuntimeException("!!");
				}
				if(Math.abs(sum1-1)>SimpleDistribution.tolerance) {
					throw new RuntimeException("!! "+sum1);
				}
			}
			
		}
		  public  void makeValid(){
				for(int i=0; i<rates.rows(); i++){
					double sum=0;
					
					for(int j=0; j<rates.rows(); j++){
						//
						if(j!=i) {
							double r = rates.getQuick(i,j);
							sum+=r;
						}
					}
					rates.setQuick(i, i, -sum);
				}
		    }
//		  double[] sumC;
		public void transfer(double pseudo, double[][] counts, MatrixExp expected, double distance){
	    if(true) throw new RuntimeException("!!");
			expected.setDistance(distance);
	    	for(int i=0; i<rates.rows(); i++){
	    		double sum =0;
	    		double sumC = Constants.sum(counts[i]) + expected.M.viewRow(i).zSum()*pseudo;
	    		
	    		for(int j=0; j<rates.rows(); j++){
	        		if(i!=j){
	        			double v =  (counts[i][j] + pseudo*expected.M.getQuick(i, j))/(sumC *distance);
	        			sum+=v;
	        			rates.setQuick(i, j,v);
	        		}
	        		
	        	}
	    	//	rates.setQuick(i, i, -sum);
	    	}
	    	//this.makeValid();
	    	this.makeValid();
	    	this.calcPiRate();
	    	this.update();
	    		
	    	this.setDistance(distance);
	    	//System.err.println("updated");
	    }

		public double evaluate(double[][] counts) {
			double sum=0;
	    	for(int i=0; i<len; i++){
	    		for(int j=0; j<len; j++){
	        			sum += counts[i][j] * Math.log(M.getQuick(i, j));
	        	}
	    	}
	    	return sum;
		}

		public void print(PrintWriter pw) {
			pw.print("--transitionMatrix  ");
			pw.println(print(rates,this.currentRate).toString());
			pw.print("STATIONARY_DIST\t");
			pw.println(print(this.getNullSpace(),1 ).toString());
		}
		static StringBuffer print(DoubleMatrix2D rates, double mod){
			StringBuffer sb = new StringBuffer(String.format("%5.3g", mod)+":");
			
			for(int i=0; i<rates.rows(); i++){
				for(int j=0; j<rates.columns(); j++){
					sb.append(String.format("%5.3g",rates.getQuick(i, j)/mod).trim());
					if(j<rates.columns()-1) sb.append(";");
				}
				if(i<rates.rows()-1) sb.append(":");
			}
			return sb;
		}


		
	}