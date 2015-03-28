package lc1.stats;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.states.HaplotypeEmissionState;
import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import pal.math.UnivariateFunction;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class LinearRegression extends Regression{
	 MultivariateFunction link1 = new MultivariateFunction() {

	        public double evaluate(double x, double phen) {
	            double prob = 1 / (1 + Math.exp(-x));
	            double res = Math.pow(prob, phen) * Math.pow(1 - prob, 1 - phen);
	            return res;
	        }

	        public double getLowerBound(int i) {
	            // TODO Auto-generated method stub
	            return Double.NEGATIVE_INFINITY;
	        }

	        public double getUpperBound(int i) {
	            return Double.POSITIVE_INFINITY;
	        }

	        public double evaluate(double[] argument) {
	            return evaluate(argument[0], argument[1]);
	        }

	        public int getNumArguments() {
	            // TODO Auto-generated method stub
	            return 2;
	        }

	        public OrthogonalHints getOrthogonalHints() {
	            // TODO Auto-generated method stub
	            return null;
	        }

	    };

	    MultivariateFunction link = new MultivariateFunction() {

	        public double evaluate(double[] argument) {// double x, double phen,
	                                                    // double sd) {
	            double phen = argument[1];
	            double sd = argument[2];
	            double x = argument[0];
	            return pal.statistics.NormalDistribution.pdf(phen, x, sd);
	        }

	        public double getLowerBound(int i) {
	            // TODO Auto-generated method stub
	            return Double.NEGATIVE_INFINITY;
	        }

	        public double getUpperBound(int i) {
	            return Double.POSITIVE_INFINITY;
	        }

	        public int getNumArguments() {
	            // TODO Auto-generated method stub
	            return 3;
	        }

	        public OrthogonalHints getOrthogonalHints() {
	            // TODO Auto-generated method stub
	            return null;
	        }

	    };
	    
	    final int[][] cols; //indexed by type
	    Algebra lg = new Algebra();
	    final DoubleMatrix2D X,XT;
	    TrainableNormal tn ;
	    public LinearRegression(DataCollection dc){
	    	super(dc);
	    	 X = new DenseDoubleMatrix2D(dc.indiv().size(),4);
	         this.cols = new int[][] {new int[] {0,1}, new int[] {0,2}, new int[] {0,3}};
	         XT = new DenseDoubleMatrix2D(4,dc.indiv().size());
	    }
	    public void setGeno(int position, int i, HaplotypeEmissionState hes){
	    	super.setGeno(position, i , hes);
	        	// String key = indiv.get(i);
	 		   // HaplotypeEmissionState hes = (HaplotypeEmissionState) dc.dataL.get(key);
	 		//    PhasedDataState ds1 = (PhasedDataState)dc.data.get(key);
	 		    X.setQuick(i, 0,1);
	 		    XT.setQuick(0, i, 1);
	 		    for(int k=0; k<vals.length; k++){
	 		    
	 		    	//Comparable c = ds1.getElement(position);
	 		    	double val = exp(hes.emissions[position], vals[k]);
	 		    	X.setQuick(i,1+k, val);
	 		    	XT.setQuick(1+k, i,val);
	 		    }
	 		   
	    	//}
	    }
	    
	    @Override
	    public double calcLogLDiff(int pos_index, int phen_index,int  type){
	    	super.calcLogLDiff(pos_index, phen_index, type);
	    	regModel.setType(type, phen_index);
	        slope = this.fitModel(type);
	          double logL3 = -1*regModel.evaluate(slope);
	          double logL2 = logL3;
	          if( logL2-logL1[phen_index] >3 ){
	          	
	         regModel.lb = slope - 0.1;
	         regModel.ub = slope+0.1;
	        double slope1 = slope;
	                 slope = os1.findMinimum(slope, regModel, 3);
	                  logL2 =  -1*regModel.evaluate(slope);
	                 
	             System.err.println("comp "+" "+(logL2-logL1[phen_index])+"<-"+(logL3-logL1[phen_index])+" ; "+slope+"<-"+slope1);    
	          }
	          try {
	         
	              if (logL2 < logL1[phen_index]) {
	                  return  0.0;
	               
	              }
	              return logL2 - logL1[phen_index];
	          } catch (Exception exc) {
	              exc.printStackTrace();
	          }
	          return 0.0;
	    }
	    class RegressionModel implements UnivariateFunction {
	        double beta = 0;
	        private int type =0;
	        int[] rowsR = null;
	        DoubleMatrix2D Z1=null;
	        void setType(int type, int phen_index){
	        	this.type = type;
	        	this.rowsR = rows[phen_index];
	        	this.Z1 = Z[phen_index];
	        }
	        TrainableNormal tn1 = new TrainableNormal(null,0,1,10000, 1.0);
	        // double stddev = 0;

	        public double evaluate(double beta) {
	           tn1.initialise();
	            this.beta = beta;
	            for(int kk=0; kk<rowsR.length; kk++){
	            	int ind = rowsR[kk];
	            	double phenotype = Z1.getQuick(ind, 0);
	            	
	                PseudoDistribution dist = covariate[ind];
	                Integer fixed = dist.fixedInteger();
	              //  double pred = 0;
	                if (fixed != null) {
	                    double pred = fixed *beta;
	                    tn1.addCount(phenotype-pred, 1.0);
	                } else {
	                    double[] probs = dist.probs();
	                    for (int i = 0; i < probs.length; i++) {
	                        double pred = vals[type][i] *beta;
	                        tn1.addCount(phenotype-pred, probs[i]);
	                       // pred += probs[i] * beta * i;
	                    }
	                }
	              
	            }
	            tn1.maximise(0, 0, 0);
	            double lh = tn1.calcLH();
	            //System.err.println("comp "+beta+" "+lh);
	            return -1*lh;
	           
	        }

	       
	double lb=-10;
	double ub = 10;
	        public double getLowerBound() {
	           
	            return lb;
	        }

	       

	        public OrthogonalHints getOrthogonalHints() {
	            // TODO Auto-generated method stub
	            return null;
	        }

	        public double getUpperBound() {
	           return ub;
	        }
	    }
	    final RegressionModel regModel = new RegressionModel();
public  double calcNull(int phen_index){
	if(tn==null)tn= new TrainableNormal(null,0,1,10000, 1.0);
	tn.initialise();
	DoubleMatrix2D Yi = Y[phen_index];
	for(int k=0; k<Yi.rows(); k++){
		tn.addCount(Yi.getQuick(k, 0), 1.0);
	}
	 tn.maximise(0, 0, 0);
      //  tn1.maximise(0, 0, 0);
     return tn.calcLH();
}
	    public double fitModel(int type){
	   	 
	  	  DoubleMatrix2D Xi = X.viewSelection( rows[this.phen_index],cols[type]);
	  		 DoubleMatrix2D XTi = //lg.transpose(Xi);
	  		 XT.viewSelection(cols[type], rows[this.phen_index]);
	  		
	  		 try{
	         DoubleMatrix2D res = lg.solve(lg.mult(XTi, Xi),lg.mult(XTi,Y[phen_index]));
	        return res.getQuick(1, 0);
	  	 }catch(Exception exc){
	  		 
	  		 return 0.0;
	  	 }
	      }
}