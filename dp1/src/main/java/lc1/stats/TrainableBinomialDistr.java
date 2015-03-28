package lc1.stats;

import org.apache.commons.math.MathRuntimeException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import pal.math.UnivariateFunction;
import pal.math.UnivariateMinimum;
import cern.colt.matrix.DoubleMatrix2D;


public class TrainableBinomialDistr extends TrainableNormal2 implements UnivariateFunction{

	public TrainableBinomialDistr(String name, double meanx, double stddevx,
			double stddevxprior, double meany, double stddevy,
			double stddevyprior) {
		super(name, meanx, stddevx, stddevxprior, meany, stddevy, stddevyprior);
		this.name = name;
		binom = new BinomialDistributionImpl(1,meany);
		
		// TODO Auto-generated constructor stub
	}
	@Override
	public double probability(double r, double b, int mixComponent) {
		// TODO Auto-generated method stub
	 return this.probability(r, b);
	}

	public TrainableBinomialDistr(TrainableBinomialDistr dist){
		super(dist);
		binom = new BinomialDistributionImpl(1,meany);
	}
	
	@Override
	public ProbabilityDistribution2 clone(double u,
			SimpleExtendedDistribution1 dist1) {
		return new TrainableBinomialDistr(this);
	}
	
	static lc1.stats.NormalDistribution normal;
	BinomialDistribution binom;
	@Override
	    public double probability(double x1, double y1) {
	    	binom.setNumberOfTrials((int)x1);
	    	double p = binom.probability(y1);
	    	//if(y1 < 0.75 * x1);
	    	return p;
	  
	    }
	
	 public double probabilityLog(double x1, double y1) {
	    	binom.setNumberOfTrials((int)x1);
	    	return Math.log(binom.probability(x1));
	    	
	  
	    }
	 
	 public int fill(DoubleMatrix2D x, DoubleMatrix2D y, DoubleMatrix2D yB, int numObs, double[] noCop,   double pseudo){
		 if(true) throw new RuntimeException("!!");
		 return 0;
	 }
	 
	 @Override
	  public void maximise(double pseudo,  double pseudoSD, double pseudoSkew, 
	    		double pseudo1, double pseudoSD1, double pseudoSkew1, double pseudoCov
	    		) {
		 if(this.x_obs.size()>0){
	    	        UnivariateMinimum uv = new UnivariateMinimum();
	    	       meany = uv.findMinimum(this.meany, this, 5);
	                this.recalcName();
	            // System.err.println("afteer "+location+" "+scale+" "+shape);
	 }
	       
	      this.paramIndex++;
	        // TODO Auto-generated method stub
	        
	    }

	@Override
	public double evaluate(double arg0) {
		if(x_obs.size()==0) return 0;
		double l =0;
		if(arg0 > getUpperBound() || arg0  < getLowerBound()) return Double.POSITIVE_INFINITY; 
		try{
		//	System.err.println(arg0);
        binom.setProbabilityOfSuccess(arg0);
		}catch(MathRuntimeException exc){
			exc.printStackTrace();
			return Double.POSITIVE_INFINITY;
		}
        for(int i=0; i<x_obs.size(); i++){
            
             double valx =x_obs.get(i);
             double valy = y_obs.get(i);
             binom.setNumberOfTrials((int)valx);
             
             double weighti =weight.get(i);
         //    if(weight>1e-5){
                 double prob = Math.log(binom.probability(valy));
              //  if(prob==Double.NEGATIVE_INFINITY){
                  //  System.err.println("doubleNeg");
                //  throw new RuntimeException("!!");
                //}
                 l+=prob*weighti;
            // }
         }
        return -l;
	}

	@Override
	public double getLowerBound() {
		// TODO Auto-generated method stub
		return 0.0;
	}

	@Override
	public double getUpperBound() {
		// TODO Auto-generated method stub
		return 1.0;
	}
	    

	
}
