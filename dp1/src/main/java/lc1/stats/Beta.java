package lc1.stats;

import java.util.Iterator;

import lc1.util.Constants;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import pal.math.ConjugateDirectionSearch;
import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;

public class Beta extends TrainableNormal implements MultivariateFunction{
	
	//double alpha, beta;
 
    public Beta clone(){
        return new Beta(this);
    }
    public Beta clone(double u){
        return new Beta(this, u);
    }

@Override
public void setParam(int type,  double d){
	super.setParam(type, d);
	modifyBeta();
	}





private void modifyBeta() {
	double loc1 = Math.max(1e-10, Math.min((1-1e-10), this.location));
	
	double v = Math.pow(this.scale, 2.0);
	double locminus = (loc1*(1-loc1));
	double n = Math.max(1e-20, locminus/v -1);
	
	//if(n<0.0001){
	//	System.err.println("h");
	//}
	double alpha = n * loc1;
	double beta = n*(1-loc1);
	this.dist = new BetaDistributionImpl(alpha,beta);
	//double[] d = new double[] {this.probability(0.0),this.probability(0.5), this.probability(1.0)};
	//System.err.println(d);
}
public Beta(String name, double mean, double stddev, double round, double priorMod){
    super(name, mean, stddev, round,priorMod); 
     this.modifyBeta();   
   
}
public Beta(String name, double mean, double stddev, double stddevprior, double round, double priorMod){
    super(name, mean, stddev, stddevprior,
            round, priorMod
           );
   this.modifyBeta();
}

public Beta(String string, int i) {
   this(null, 0,0,0, 1.0);
    throw new RuntimeException("!!");
}

public Beta(Beta trainableNormal) {
    super(trainableNormal);
    this.modifyBeta();
    // TODO Auto-generated constructor stub
}

public Beta(Beta skewNormal, double u) {
    this(skewNormal);
    
   // this.setPrior(u, u, u);
    this.location =normal.quantile(Constants.rand.nextDouble(),location, 1.0/u);
    this.modifyBeta();
  }
  
public double dsn (double x, boolean log)
{
 
 return log ? Math.log(this.probability(x)) : this.probability(x);
}




    @Override
    public double cumulative(double arg0) {
    	try{
    	return dist.cumulativeProbability(arg0);
    	}catch(MathException exc){
    		return Double.NaN;
    	}
       // return normal.cdf(arg0, this.location, this.scale);
    }

    @Override
    public double inverse(double arg0) {
    	try{
    	return dist.inverseCumulativeProbability(arg0);
    }catch(MathException exc){
		return Double.NaN;
	}
    //    return normal.quantile(arg0,location, scale);
//        throw new RuntimeException("!!");
//       return normal.inverse(arg0);
       
    }

   
    BetaDistribution dist = new BetaDistributionImpl(1,1);
    static double incr = 1e-10;
    @Override
    public double probability(double arg0) {
    	try{
    		//this.dist.
    		double min = Math.max(0, arg0-incr);
    		double max = Math.min(1, arg0+incr);
    	double res =   (this.dist.cumulativeProbability(arg0+incr) - 
    		this.dist.cumulativeProbability(arg0-incr))/(2*incr);
    //	double res = dist.density(arg0);
//    	double res = normal.pdf(arg0, location, scale);
    	
    return res;
    	}catch(Exception exc){
    		exc.printStackTrace();
    		return Double.NaN;
    	}
   //  return sc;
  /*      if(prob==0){
        	throw new RuntimeException("!!");
        }
       return prob;*/
    }

  
        @Override
    public void maximise(double pseudo, double pseudoSD, double pseudoSkew) {
      super.maximise(pseudo, pseudoSD, pseudoSkew);
      this.modifyBeta();
    //    	this.maximise1(pseudo, pseudoSD, pseudoSkew);
   // scale = Math.max(scale, 0.0001);
    }
 
     
    
    public double probability(double r, double offset) {
    	
      return probability(r-offset);
//       normal.setState(mean-offset, stddev);
    }
    public double rsn() {
        return this.inverse( Constants.rand.nextDouble());
    }
public double max(){
	double m = 0;
	for(Iterator<Double> it = obsv.iterator(); it.hasNext();){
		Double d = it.next();
		if(d>m) m = d;
	}
	return m;
	}
    public void maximise1(double pseudoM, double pseudoSd, double pseudoSk) {
   // double max = this.max();	
        if(this.max() < Constants.trainThresh() ) return;
                final ConjugateDirectionSearch os = new ConjugateDirectionSearch();
                final OrthogonalSearch os1 = new pal.math.OrthogonalSearch();
                final double[] xvec =  new double[] {Math.log(this.dist.getAlpha()), Math.log(this.dist.getBeta())};
                final double[] init = new double[xvec.length];
                System.arraycopy(xvec, 0, init, 0, xvec.length);
                Runnable run = new Runnable(){
                public void run(){
                	try{
                  os1.optimize(Beta.this,xvec, 0.01, 0.01);
                	}catch(Exception exc){
                		exc.printStackTrace();
                		//os1.optimize(Beta.this,xvec, 0.01, 0.01);
                	}
                }
                };
                Thread th = new Thread(run);
                th.run();
                try{
                for(int i=0; i<100; i++){
                    Thread.sleep(100);//this.wait(100);
                    if(!th.isAlive()) break;
                }
                if(th.isAlive()){
                    
                    th.stop();
                    System.arraycopy(init, 0, xvec, 0, xvec.length);
                }
                }catch(Exception exc){
                    exc.printStackTrace();
                }
                //location = xvec[0];
                //scale = xvec[1];
               // shape = xvec[2];
                this.recalcName();
            // System.err.println("afteer "+location+" "+scale+" "+shape);
           
       
      this.paramIndex++;
        // TODO Auto-generated method stub
        
    }

    public double calcLH(){
            double l =0;
            double shape =0;
          
            for(int i=0; i<obsx.size(); i++){
                 double val =obsx.get(i);
                 double weight =obsv.get(i);
             //    if(weight>1e-5){
                     double prob = this.dsn(val, true);
                   
                    if(prob==Double.NEGATIVE_INFINITY){
                    	  if(weight > Constants.trainThresh()){
                     /*  System.err.println("doubleNeg");
                       double res =  -1e6*(
                                  Math.pow(location - this.meanPrior, 2) + 
                                  Math.pow(scale -  this.stddevPrior, 2)
                              //    Math.pow(shape - skewPrior, 2)
                                );
                       prob = res;*/
                    }
                    	  else{
                         	 prob=0;
                          }
                     }
                     l+=prob*weight;
                // }
             }
            return l;
        }
    
    
    public double evaluate(double[] argument) {
          if(Double.isNaN(argument[0])) throw new RuntimeException("!!"+argument[0]+" "+argument[1]+" "+argument[2]);
         this.dist.setAlpha(Math.exp(argument[0]));
         this.dist.setBeta(Math.exp(argument[1]));
         try{
        this.location = this.dist.inverseCumulativeProbability(0.5);
         }catch(MathException exc){
        	 exc.printStackTrace();
         }
         double res = -1*this.calcL();
        
//         System.err.println(argument[0]+" "+argument[1]+" "+res);
         return res;
      }
        public double calcL(){
           return calcLH();
        }

        public double getLowerBound(int n) {
            return -100;

        }
        public int getNumArguments() {
          return 2;
        }
        public OrthogonalHints getOrthogonalHints() {
            return null;
        }
       
        public double getUpperBound(int n) {
            return 100;
       
        }

}
