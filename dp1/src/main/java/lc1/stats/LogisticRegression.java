package lc1.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.TimeUnit;

import lc1.dp.data.collection.DataCollection;
import pal.math.ConjugateDirectionSearch;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;
import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;
import cern.colt.matrix.DoubleMatrix2D;

public class LogisticRegression  extends Regression{
    DataCollection dc;
//    List<PseudoDistribution> data = new ArrayList<PseudoDistribution>();
 //  final  List<Double> phenotype = new ArrayList<Double>();
  // final List<PseudoDistribution> covariate = new ArrayList<PseudoDistribution>();
    UnivariateFunction link = new UnivariateFunction(){

        public double evaluate(double x) {
          return  1/(1+Math.exp(-x));
        }

        public double getLowerBound() {
            // TODO Auto-generated method stub
            return Double.NEGATIVE_INFINITY;
        }

        public double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
        
    };
    
 
  //  double stddev;
   
    
   
    
    public LogisticRegression(DataCollection dc){
       super(dc);
       
    }
  
    final UnivariateMinimum os = new UnivariateMinimum();
    final MultivariateMinimum os1 = new ConjugateDirectionSearch();
    final RegressionModel regModel = new RegressionModel();
    final double[] xvec =  new double[] {0};
    final double[] xvec1 =  new double[] {0,0};
    @Override
    public  double calcNull(int phen_index){
    	double cnt0=0;
    	double cnt1=0;
    	if(alpha==null) alpha = new double[this.Y.length];
    	DoubleMatrix2D m = this.Y[phen_index];
    	for(int i=0; i<m.rows(); i++){
    		if(m.getQuick(i, 0)>0.5){
    			cnt1++;
    		}
    		else{
    			cnt0++;
    		}
    	}
    	double p = cnt1/(cnt0+cnt1);
    	alpha[phen_index] =p==0 ? -10 : (p==1 ? 10 : Math.log( p/(1-p)));
    	if(p==0 || p==1){
    		return 0.0;
    	}
    	try{
    	//this.binomial.setNandP((int) (cnt0+cnt1), p);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    	//return B
    	
    	double res = Math.log(Math.pow(p,cnt1))+Math.log(Math.pow(1-p,cnt0));
    	return res;
    	
    }
    double[] alpha;
    public double calcLogLDiff(int pos_index, int phen_index, int type){
    	//System.err.println(pos_index+"_"+phen_index+"_"+ type);
    	super.calcLogLDiff(pos_index, phen_index, type);
       // this.position = pos_index;
     
  
      
        Arrays.fill(xvec1, 0);
        List call= new ArrayList();
       regModel.setType(type,phen_index);
      xvec1[0] = alpha[phen_index];
        call.add(
         new Callable(){
           public Object call(){
              // os1.findMinimum(regModel, xvec1, 2, 2);
        os1.optimize(regModel, xvec1, 1, 1);
              return null;
          }
        });
        try{
        es.invokeAll(call, 5, TimeUnit.SECONDS);
       //((Callable)call.get(0)).call();
       
        double logL2 = logL1[phen_index]==0 ? 0 : -1*regModel.evaluate(xvec1);
        if(logL2 <= logL1[phen_index]) {
            System.err.println("warning = poor maximisation");
            return 0.0;
//           double[] xvec2 = new double[] {0,xvec[0]};
 //          logL2 = -1 *regModel.evaluate(xvec2);
        }
        else if(logL2>logL1[phen_index]+1);
        {        	call.clear();
        	call.add(
        	         new Callable(){
        	           public Object call(){
        	              // os1.findMinimum(regModel, xvec1, 2, 2);
        	        os1.optimize(regModel, xvec1, 0.1, 0.1);
        	              return null;
        	          }
        	        });
        	es.invokeAll(call, 5, TimeUnit.SECONDS);
        	logL2 = -1*regModel.evaluate(xvec1);
        }
      //  System.err.println("minimizer1 "+xvec[0]);
        System.err.println("minimizer2 "+xvec1[0]+" "+xvec1[1]);
        return logL2-logL1[phen_index];
        }catch(Exception exc){
            exc.printStackTrace();
        }
       return 0.0;
    }
    
   
    
    static double[] max = new double[] {100,100};
    class RegressionModel implements MultivariateFunction{
        double beta=0;
        double alpha=0;
        public double evaluate(double[] argument) {
            alpha = argument[0];   
            beta = argument[1];
            double logL=0;
            for(int kk=0; kk<rowsR.length; kk++){
            	int ind = rowsR[kk];
            	double phenotype = Z1.getQuick(ind, 0);
            	
                PseudoDistribution dist = covariate[ind];
                Integer fixed = dist.fixedInteger();
              
                   logL+=Math.log(eval(dist, phenotype));  
                 
              
           }
              
          // System.err.println(alpha+" "+beta+" "+logL);
               return -1*logL;
        }
        private int type =0;
        int[] rowsR = null;
        DoubleMatrix2D Z1=null;
        void setType(int type, int phen_index){
        	this.type = type;
        	this.alpha = LogisticRegression.this.alpha[phen_index];
        	this.rowsR = rows[phen_index];
        	this.Z1 = Z[phen_index];
        }
        
		protected double eval(PseudoDistribution dist, double phen) {
            double pred = alpha;
            Integer fixed = dist.fixedInteger();
            if(fixed!=null){
               pred +=beta*vals[type][fixed];
            }
            else{
                double[] probs = dist.probs();
                for(int i=0; i<probs.length; i++){
                    pred +=probs[i]*beta*vals[type][i];
                }
            }
            double prob = link.evaluate(pred);
            //double prob = normal.pdf(phen, linkL, stddev);
            double res = Math.pow(prob, phen)*Math.pow(1-prob, 1-phen);
            return res;
        }
        public double getLowerBound(int n) {
         return -max[n];
        }
        public int getNumArguments() {
           return 2;
        }
     
        public OrthogonalHints getOrthogonalHints() {
            // TODO Auto-generated method stub
            return null;
        }
        public double getUpperBound(int n) {
           return max[n];
         }
    }
}
