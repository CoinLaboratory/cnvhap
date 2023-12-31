package lc1.stats;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.ConjugateDirectionSearch;
import pal.math.MultivariateFunction;
import pal.math.OrthogonalHints;
import pal.math.OrthogonalSearch;
import JSci.maths.statistics.ProbabilityDistribution;
import cern.colt.matrix.DoubleMatrix2D;

public class SkewNormal extends ProbabilityDistribution implements lc1.stats.ProbabilityDistribution, MultivariateFunction, Comparable{
   int paramIndex = 1;
  // String id;
    public int getParamIndex(){
        return this.paramIndex;
   }
    public void setCoverage(double d){
    	throw new RuntimeException("!!");
    }
    
    public SkewNormal clone(){
        return new SkewNormal(this);
    }
    public double  probability(double x, int mixComponent){
		return this.probability(x);
	}
    public lc1.stats.ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1){
		return this.clone(u);
	}
    public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, lc1.stats.ProbabilityDistribution disty){
		throw new RuntimeException("!!");
	}
    public void setPriors(lc1.stats.ProbabilityDistribution distx, int type){
    	SkewNormal dist1 = (SkewNormal)distx;
    	if(type==1)this.stddevPrior = dist1.scale;
    	else if (type==2) throw new RuntimeException("!!");//this.skewPrior[0] = dist1.shape;
    	else if(type==0) this.meanPrior = dist1.location;
		//this.dist[0].setPriors(((Mixture)distx).dist[0]);
	}
    public SkewNormal clone(double u){
        return new SkewNormal(this, u);
    }
   /* public void setPriorVar(double e){
		this.stddevPrior = e;
		this.scale = e;
	}*/
    public void transfercounts(EmissionState innerState, int phen_index, int i){
       throw new RuntimeException("!!");
       /* for(Iterator<Map.Entry<Double, Double>> it = this.observations.entrySet().iterator(); it.hasNext();){
           Map.Entry<Double, Double> nxt = it.next();
                innerState.addCountDT(nxt.getKey(),phen_index,  nxt.getValue(), i);
        *}*/
    }
    public void print(PrintWriter pw){
    	this.recalcName();
    	pw.print(this.toString()+"\t");
    }
    @Override
	public int numObs() {
		return this.obsx.size();
	}
    public double scale(){
    	return this.scale;
    }
    public double[] getCount(double[] angle){
        throw new RuntimeException("!!");
    }
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
        SkewNormal n = (SkewNormal) probabilityDistribution;
       
        for(int i=0; i<n.obsx.size(); i++){
        	this.addCount(n.obsx.get(i), n.obsv.get(i));
        }
       /* for(Iterator<Map.Entry<Double,Double>> it = n.observations.entrySet().iterator(); it.hasNext();){
            Map.Entry<Double,Double> nxt = it.next();
            this.addCount(nxt.getKey(),  nxt.getValue());
        }*/
    }
    public static void main(String[] args){
     
       SkewNormal sn = new SkewNormal(0.0,0.05,1e10, -5, 5, 100, 1);
      // SkewNormal sn1 = new TrainableNormal(0.0,0.05,1e10, -5, 5, 100, 1);
       double[] x = new double[] {0.0, 0.5, 1.0};
       for(int i=0; i<x.length; i++){
       System.err.print(sn.dsn(x[i], false)+"\t");
       }
       System.err.println();
       
        sn = new SkewNormal(-1.0,0.2,0, -5, 5, 100,1);
       for(int i=0; i<x.length; i++){
       System.err.print(sn.dsn(x[i], false)+"\t");
       }
       System.err.println();
       
       sn = new SkewNormal(-1.0,0.2,10, -5, 5, 100,1);
       for(int i=0; i<x.length; i++){
       System.err.print(sn.dsn(x[i], false)+"\t");
       }
       System.err.println();
    
    }
    public  double getParamValue(int i){
        if(i==0) return this.location;
        if(i==1) return this.scale;
       // if(i==2) return this.shape;
        throw new RuntimeException("!!");
    }
    public double probability(double r, double offset) {
        return this.probability(r);
          
      }
    private String name;
    public String name(){
        return id;
    }
    public boolean equals(Object obj){
        return this.name.equals( ((SkewNormal)obj).name);
    }
    public int hashCode(){
        return this.name.hashCode();
    }
    public double sum(){
      return sum;
    }
    double sum;
    
    public List<Double> obsx = new ArrayList<Double>();
    List<Double> obsv = new ArrayList<Double>();
    //public SortedMap<Double, Double> observations = new TreeMap<Double, Double>();
    
    
    /*public void log(PrintWriter pw){
        
    }*/
    
    final double round;
public double round(double d){
    return Math.round(d*round)/round;
}
    public void addObservation(double d1, double weight){
        //double d1 = round(d);
    	this.obsx.add(d1);
    	this.obsv.add(weight);
//        Double w = observations.get(d1);
      //  observations.put(d1, w==null ? weight : weight+w);
        sum+=weight;
    }
  static lc1.stats.NormalDistribution normal;
public double location, scale;//, shape;//origLoc, origScale, origShape;
static double pi = Math.PI;
//final double min, max;
double[] lower = new double[] {-5, 0.001, -1e10};
double[] upper = new double[] {5, 1e3, 1e10};

public  double meanPrior, stddevPrior;// skewPrior;
public SkewNormal(double mean, double stddev, double skew, double min, double max, double round, double priorModifier) {
    this(null, mean, stddev, skew, new double[] {min, 0.001, -1e10},  new double[] {max, 1e3, 1e10}, round, priorModifier);
 }
public SkewNormal(String name, double mean, double stddev, double skew, double min, double max, double round, double priorModifier) {
    this(name, mean, stddev, skew, new double[] {min, 0.001, -1e10},  new double[] {max, 1e3, 1e10}, round, priorModifier);
 }
public void setParams(double mean, double stddev, double skew) {
    this.location = mean;
    this.scale = stddev;
    //this.shape = skew;
    this.meanPrior = location;
    this.stddevPrior = scale;
//    this.skewPrior = new double[] {shape,1};
    
}
final public String id; //
    public SkewNormal(String name1, 
    		double mean, double stddev, double skew, double[] lower , double[] upper, double round, double priorModifier){
       this.priorModifier = priorModifier;
        this.round = round;
        this.lower =lower;
        this.upper = upper;
         this.id = name1;
     setParams(mean, stddev, skew);
    this.recalcName();
        if(stddev<0) throw new RuntimeException("!!");
    }
    public SkewNormal(String name1, double mean, double stddev,double stddevprior, double skew, double[] lower , double[] upper, double round, double priorModifier){
        this.priorModifier = priorModifier;
         this.round = round;
         this.lower =lower;
         this.upper = upper;
         this.id = name1;
      setParams(mean, stddev, skew);
      this.stddevPrior = stddevprior;
     this.recalcName();
         if(stddev<0) throw new RuntimeException("!!");
     }
   
    public SkewNormal(String string, double mean_i_r, double var_i_r, double skew,
			double round2, double d) {
    	this(string, mean_i_r, var_i_r, skew, new double[] {mean_i_r-10, 1e-3,0}, new double[] {mean_i_r+10,10,1e10 }, round2, d  );
		// TODO Auto-generated constructor stub
	}
    
    
    public SkewNormal(SkewNormal skewNormal) {
        this.priorModifier = skewNormal.priorModifier;
        this.round = skewNormal.round;
        this.lower =skewNormal.lower;
        this.upper = skewNormal.upper;
        this.id = skewNormal.id;
     setParams(skewNormal.location, skewNormal.scale, 0.0);//skewNormal.shape);
     this.name = skewNormal.name;
   obsx = new ArrayList<Double>();
   obsv = new ArrayList<Double>();
    }
    
    
    public SkewNormal(SkewNormal skewNormal, double u) {
      this(skewNormal);
     this.location =   normal.quantile(Constants.rand.nextDouble(),location, 1.0/u);
    }
    
   
	
	public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution[] tmp) {
        // TODO Auto-generated method stub
        
       this.location =0;       this.scale =0;
      // this.shape = 0;
     
       for(int i=0; i<tmp.length; i++){
           location += ((SkewNormal)tmp[i]).location;
           this.scale += ((SkewNormal)tmp[i]).scale;
         //  this.shape += ((SkewNormal)tmp[i]).shape;
       }
       location = location / (double)tmp.length;
       scale = scale/ (double)tmp.length;
    //   shape = shape / (double)tmp.length;
       
    }
    @Override
    public double cumulative(double arg0) {
       throw new RuntimeException("!!");
    }

    @Override
    public double inverse(double arg0) {
       throw new RuntimeException("!!");
    }

    @Override
    public double probability(double arg0) {
    	
    double res =  this.dsn(arg0, false);
    if(Constants.CHECK && Double.isNaN(res)){
        throw new RuntimeException("!!"+arg0+" "+this.location+" "+this.scale);//+" "+this.shape);
    }
   return res;
    }

    public double dsn (double x, boolean log)
    {
    	double shape=0;
    	//if(true) throw new RuntimeException("!!");
     double z  =  (x-location)/scale;
     
     double y;
      if(log)
        y = (-0.9189385332046727-Math.log(scale)-Math.pow(z, 2)/2+zeta(0,shape*z));
      else
        y = 2*dnorm(z)*pnorm(z*shape, false)/scale;
      if(scale <=0) y = Double.NaN;
//      replace(y, scale<= 0, NaN)
      return y;
    }
   private static double dnorm(double z){
       return normal.pdf(z,0,1);
   }
   private static double pnorm(double z, boolean log){
      if(log) return Math.log(normal.cdf(z,0,1));
      else return normal.cdf(z,0,1);
       //return 0;
   }
   
   public double getMean(){
       return this.location;
   }
   public double getStdDev(){
       return this.scale;
   }
    public static Double zeta(int k, double x){// function(k,x){# k integer in (0,4)
       
        if(k<0 | k>4 | k != Math.round(k)) return null;
        k  = Math.round(k);
        boolean na= Double.isNaN(x);
        if(na) x = 0;
        Double z=null;
        if(k==0){
            z =   pnorm(x, true)+ Math.log(2); //logp.p=true
        }
        else if(k==1){
          z =   (x>-20) ?  dnorm(x)/pnorm(x, false) : 
                 (  (x>-200) ? Math.exp(-Math.pow(x,2)/2-0.5*Math.log(2*pi) - pnorm(x,true)) : 
                            -x*(1+1/Math.pow(x,2)-2/Math.pow(x,4)));
        }
       
        else if(k==2){
           z =  (-zeta(1,x)*(x+zeta(1,x)));
        }
        else if(k==3){
           z =  (-zeta(2,x)*(x+zeta(1,x))-zeta(1,x)*(1+zeta(2,x)));
        }
        else if(k==4){
           z =  (-zeta(3,x)*(x+2*zeta(1,x))-2*zeta(2,x)*(1+zeta(2,x)));
        }
        if(x==Double.NEGATIVE_INFINITY){
            if(k==0){}// z = z;
            else if(k==1) z = z.doubleValue() == Double.NEGATIVE_INFINITY ? Double.POSITIVE_INFINITY : z;
            else if(k==2) z =z.doubleValue() == Double.NEGATIVE_INFINITY ? 1.0 : z;
            else if (k==3 || k==4)z = z.doubleValue() == Double.NEGATIVE_INFINITY ? 0.0 : z;
            else if(k>4) z = null;
            if(k>1) z = x ==Double.POSITIVE_INFINITY ? 0.0:z;
            z = na ? Double.NaN :z;
      }
    return z;
    }
    //static final double cntThresh = Constants.countThresh();
    
    public  synchronized  void addCount(double d, double w) {
        if(Constants.CHECK && (Double.isNaN(d) || Double.isInfinite(d))) throw new RuntimeException("!!");
       // if( w > cntThresh){ //throw new RuntimeException("!! "+d);
          this.addObservation(d, w);
       // }
    }
    public void initialiseCounts() {
        obsx.clear();obsv.clear();
        sum =0;
    }
    
    @Override
    /** if pseudo < 0 we take prior value */
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[]  noCop,
 double pseudo) {
    	int i;
    	int nocols = noCop.length;
    	 for( i=0; i<this.obsx.size(); i++){
    		 int k = numObs+i;
    		 double v = this.obsv.get(i);
    		 for(int kk=0;kk<nocols; kk++){
    			 x.setQuick(k, kk, v * (double)noCop[kk]);
    		 }
    		
    		
    		
    		 y.setQuick(k,0, v*( this.obsx.get(i) ));
    	 }
    	if(Math.abs(pseudo)>0){
    		double v = Math.abs(pseudo);
    		int k = numObs+i;
	   		 for(int kk=0;kk<nocols; kk++){
	   			 x.setQuick(k, kk, v * (double)noCop[kk]);
	   		 }
	   		
	   		
	   		
	   		 y.setQuick(k,0, v*(pseudo < 0  ? this.meanPrior : this.location ));
	   		i++;
    	}
    	 return i;
	}
    @Override
    public int fillVariance( DoubleMatrix2D y, int numObs, double pseudo) {
    	int i;
    	
    	 for( i=0; i<this.obsx.size(); i++){
    		 int k = numObs+i;
    		 double v = this.obsv.get(i);
    		
    		 y.setQuick(k,0, v*( Constants.transformVariance(this.obsx.get(i) - this.location )));
    	 }
    	 if(Math.abs(pseudo)>0){
    		 int k = numObs+i;
    		 double v = Math.abs(pseudo);
    		
    		 y.setQuick(k,0, v*( Constants.transformVariance(pseudo<0 ? this.stddevPrior : this.scale )));
    		 i++;
    	 }
    	 return i;
	}
    
    
    public String getObsString(){
    //  if(obsx.g) return median()+"";
    //  else return "-";
    	return this.name;
    }
    public double median(){
     throw new RuntimeException("!!");
       /* double sum=0;
        double tot = sum();
      //  StringBuffer sb = new StringBuffer(tot+"_sum ");
       // Histogram hist = new Histogram(observations.firstKey(), observations.lastKey(), 10);
        for(Iterator<Double> it = this.observations.keySet().iterator(); it.hasNext();){
           Double d =  it.next();
           sum+=observations.get(d)/tot;
           if(sum>0.5) return d;
        //   sb.append("("+d+":"+sum+"),");
        }
        return observations.lastKey();//  */    //  return sb.toString();
    }
    
    public void transfer(double ps){
    	if(this.obsx.size()==0) return;
        this.maximise(ps, ps, ps);
        
    }
    public void maximise(double pseudoM, double pseudoSd, double pseudoSk) {
    	//if(!Constants.trainSkew()) return ;
    //	if(true) throw new RuntimeException("!!");
        if(this.sum() < Constants.trainThresh() ) return;
      
        setPrior(pseudoM, pseudoSd, pseudoSk);
       
      
         
              //  origLoc = location;
             //   origScale = scale;
              //  origShape = shape;
             //   System.err.println("initial "+location+" "+scale+" "+shape+"\n Priors"+this.meanPrior.getMean()+" "+this.stddevPrior.getMean()+" "+this.skewPrior.getMean());
            //   Logger.global.info(this.getObsString());
                final ConjugateDirectionSearch os = new ConjugateDirectionSearch();
                final OrthogonalSearch os1 = new pal.math.OrthogonalSearch();
                final double[] xvec =  new double[] {location, scale, 0};
                final double[] init = new double[xvec.length];
                System.arraycopy(xvec, 0, init, 0, xvec.length);
                Runnable run = new Runnable(){
                public void run(){
                	try{
                  os.optimize(SkewNormal.this,xvec, 0.01, 0.01);
                	}catch(Exception exc){
                		exc.printStackTrace();
                		os1.optimize(SkewNormal.this,xvec, 0.01, 0.01);
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
                location = xvec[0];
                scale = xvec[1];
               // shape = xvec[2];
                this.recalcName();
            // System.err.println("afteer "+location+" "+scale+" "+shape);
           
       
      this.paramIndex++;
        // TODO Auto-generated method stub
        
    }
    
    
  /*  public void setPrior(SkewNormal bground,double pseudo, double pseudoSD, double pseudoS  ){
        this.meanPrior = new NormalDistribution(bground.location,1.0/pseudo);
        this.stddevPrior = new NormalDistribution(bground.scale, 1.0/pseudoSD);
        this.skewPrior = new NormalDistribution(bground.shape, 1.0/pseudoSk);
        this.pseudoM = pseudo;
        this.pseudoSD = pseudoSD;
        this.pseudoSk = pseudoSk;
        if(pseudoM > 1e4) this.lower[0] = this.upper[0] =  meanPrior.getMean();
        if(pseudoSD > 1e4) this.lower[1] = this.upper[1] = stddevPrior.getMean();
        if(pseudoSk > 1e4) this.lower[2] = this.upper[2] = skewPrior.getMean();
       } */
    
    public final double priorModifier;
    
    public  void setPrior(double pseudo, double pseudoSD, double pseudoSk) {
    //	if(true) throw new RuntimeException("!!");
        this.meanPrior = 1.0/(pseudo*priorModifier);
        this.stddevPrior  =  1.0/(pseudoSD*priorModifier);
      //  this.skewPrior = 1.0/pseudoSk;
        this.pseudoM = pseudo;
        this.pseudoSD = pseudoSD;
        this.pseudoSk = pseudoSk;
        if(pseudoM > 1e4) this.lower[0] = this.upper[0] =  meanPrior;
        if(pseudoSD > 1e4) this.lower[1] = this.upper[1] =stddevPrior;
      //  if(pseudoSk > 1e4) this.lower[2] = this.upper[2] = skewPrior;
        
    }
    double pseudoM,  pseudoSD,  pseudoSk;
    public double calcLH(){
    //	if(true) throw new RuntimeException("!!");
        double l =0;
        double shape =0;
        double skewPrior =0;
        for(int i=0; i<obsx.size(); i++){
            
             double val =obsx.get(i);
             double weight =obsv.get(i);
         //    if(weight>1e-5){
                 double prob = this.dsn(val, true);
                if(prob==Double.NEGATIVE_INFINITY){
                  //  System.err.println("doubleNeg");
                   double res =  -1e6*(
                              Math.pow(location - this.meanPrior, 2) + 
                              Math.pow(scale -  this.stddevPrior, 2)+
                              Math.pow(shape - skewPrior, 2)
                            );
                   return res;
                }
                 l+=prob*weight;
            // }
         }
        return l;
    }
    public double prior(){
        double l = 0 ;
      //  if(true) throw new RuntimeException("!!");
        /*
        if(pseudoM>1e-5) l+= normal.logpdf(location,this.meanPrior, meanPrior[1]);
          
        //if(l==Double.NEGATIVE_INFINITY) 
        if(pseudoSD>1e-5) l+=normal.logpdf(scale,this.stddevPrior, stddevPrior[1]);
       // if(l==Double.NEGATIVE_INFINITY) 
        if(pseudoSk>1e-5) l+=  normal.logpdf(shape,this.skewPrior[0],skewPrior[1]);
        if(l==Double.NEGATIVE_INFINITY){
          //  Logger.global.warning("!! "+stddevPrior.getMean()+" "+stddevPrior.getVariance()+" "+scale+" "+meanPrior.probability(this.location));
           // Logger.global.warning("!! "+meanPrior.getMean()+" "+meanPrior.getVariance()+" "+location+" "+stddevPrior.probability(this.scale));
          //  Logger.global.warning("!! "+skewPrior.getMean()+" "+skewPrior.getVariance()+" "+shape+" "+stddevPrior.probability(this.scale));
                //  System.err.println("doubleNeg");
                 double res =  -1e6*(
                            Math.pow(location - this.meanPrior[0], 2) + 
                            Math.pow(scale -  this.stddevPrior[0], 2)+
                            Math.pow(shape - this.skewPrior[0], 2)
                          );
                 return l;
        }
            */
            
           ;//pseudo==0 ? 0 : Math.log(dist.probability(arg0));
        
       
    //    if(Double.isNaN(l)) throw new RuntimeException("!!");
      // System.err.println(this.location+" "+this.scale+" "+this.shape+" "+l);
        return l;
    }
    public double calcL(){
       return calcLH()+prior();
    }
  public String toString(){
	  return name;
  }
    //public String toString1(){
       
     
       // return "c("+this.location+","+this.scale+","+this.shape+");";
    //}
   /* double psn(double x)
    {
     
     double z = (x-location)/scale;
     double p = pmin(1, pmax(0, pnorm(z, false) - 2*tOwen(z, shape)));
     if(scale<=0) p = Double.NaN;
//      replace(p, scale<= 0, NaN)
    }*/
    public double evaluate(double[] argument) {
      //  
   //    if(true) System.exit(0);
        if(Double.isNaN(argument[0])) throw new RuntimeException("!!"+argument[0]+" "+argument[1]+" "+argument[2]);
       this.location = argument[0];
       this.scale = argument[1];
      // this.shape = argument[2];
     
       double res = -1*this.calcL();
       
  //  System.err.println("eval "+argument[0]+" "+argument[1]+" "+argument[2]+" "+res);
       return res;
    }
   
   public void ci(Double[] d, double[] r){
       for(int i=0; i<d.length; i++){
           d[i] = normal.quantile(r[i], location, scale);
       }
   }
    
    public void  plotTheoretical(String name1, boolean cum, XYSeries newD) {
    	//  double maxV =0;
       //   double maxX =0;
        double sum =0;
        double minq = 1e-7;
        double maxq = (1.0-1e-7) ;
        double min = normal.quantile(minq, location, scale);
        double max = normal.quantile(maxq, location, scale);
        double incr = 0.001;
        double sum1 =0;
        for(double j =min; j<max; j+=incr){
            sum1+= this.probability(j);
            
        }
        for(double j =min; j<max; j+=incr){
            double probj = this.probability(j)/sum1;
            sum+=probj*incr;
            double res = cum ? sum : probj;
            if(res>1e-11){
               
                newD.add(j, res);
            }
            
        }
       
      //  return newD;
    }
    
  //  @Override
    public void  getInterval(double[] input, double[] res) {
    	for(int i=0; i<res.length; i++){
    	res[i] =  normal.quantile(input[i], this.location, this.scale);
    	}
      //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
    }
   /* public void  getInterval(double[] in, double[] res) {
    	//  double maxV =0;
       //   double maxX =0;
        double sum =0;
        double minq = 1e-7;
        double maxq = (1.0-1e-7) ;
        double min = normal.quantile(minq, location, scale);
        double max = normal.quantile(maxq, location, scale);
        double incr = 0.001;
        double sum1 =0;
        for(double j =min; j<max; j+=incr){
            sum1+= this.probability(j)*incr;
            
        }
    //  int i=0;
        int min_i=0;
       loop: for(double j =min; j<max; j+=incr){
            double probj = this.probability(j)/sum1;
            sum+=probj*incr;
           // if(i>=in.length) break loop;
            //if(sum>in[i]){
            	//i++;
            //}
          for(int i=min_i; i<in.length; i++){
        	if  (sum <= in[i]){
          
        	   res[i] = j;
           }
        	else{
        		min_i =i;
        	}
          } 
           if(min_i>=in.length) break loop; 
        }
     //  System.err.println("d");
      //  return newD;
    }*/
  
    public double plotObservations(String name1, boolean cum, XYSeries newD, boolean swtch) {
      // double tot = sum();
     //  if(tot<Constants.trainThresh()) return tot;
        double sum1 =0;
      // Double[] res1 = null;
     /*   for(Iterator<Map.Entry<Double, Double>> it = this.observations.entrySet().iterator(); it.hasNext();){
            Map.Entry<Double, Double> nxt = it.next();
         
            sum+=nxt.getValue();
            if(!swtch){
	            double res = cum ? sum : nxt.getValue().doubleValue();
	            if(res>IndividualPlot.plotThresh){
	              
	               newD.add(nxt.getKey().doubleValue(), res);
	            }
            }
           
        }
        double sum1 = sum;
        if(swtch){
        	for(Iterator<Map.Entry<Double, Double>> it = this.observations.entrySet().iterator(); it.hasNext();){
                Map.Entry<Double, Double> nxt = it.next();
             
                sum-=nxt.getValue();
                double res = cum ? sum : nxt.getValue().doubleValue();
                if(res>IndividualPlot.plotThresh){
                  
                   newD.add(nxt.getKey().doubleValue(), res);
                //   if(res1==null && sum>=0.5*tot){
                //	   res1 = new Double[]  {sum/tot,this.mean()};
                //   }
                }
               
            }
        }*/
        return sum1;
       // return res1;
       
    }
  
    public double getLowerBound(int n) {
        return lower[n];

    }
    public int getNumArguments() {
      return 3;
    }
    public OrthogonalHints getOrthogonalHints() {
        return null;
    }
   
    public double getUpperBound(int n) {
        return upper[n];
   
    }
   
    public double rsn()
    {
    	double shape =0;
      
      double u1 =normal.quantile(Constants.rand.nextDouble(),0,1);
      double u2  = normal.quantile(Constants.rand.nextDouble(),0,1);
      boolean id  = u2 > shape*u1;
      if(id) u1 = -u1;
//      u1[id] <- (-u1[id])
      return location+scale*u1;
     
      }
   /* public Map<Double, Double> getObservations() {
       return this.observations;
    }*/
  
    public void transferParams(ProbabilityDistribution pdist0){
        this.location = ((SkewNormal)pdist0).location;
        this.scale = ((SkewNormal)pdist0).scale;
      //  this.shape = ((SkewNormal)pdist0).shape;
    }
    public void setParamValue(int n1, double val) {
      if(n1==0) location = val;
      else if(n1==1){
    	  scale = val;
    	  if(Constants.CHECK && scale==0){
          	   throw new RuntimeException("!!");
             }
      }
     // else if(n1==2) shape = val;
      else throw new RuntimeException("!!");
    }
    public void recalcName() {
    	name =   id+"_"+String.format("%5.3g,", location)+String.format("%5.3g,", scale);//+String.format("%5.3g", shape);
       
    }
    public void appendToName(String string) {
        this.name = this.name+";"+string;
        
    }
    public void updateParamIndex(){
        this.paramIndex++;
    }

   /* public double shape() {
      return shape;
    }*/
    public double sample() {
       throw new RuntimeException("!!");
    }
    public void initialise() {
        this.initialiseCounts();
        
    }
	public int compareTo(Object o) {
		return this.name.compareTo(((SkewNormal)o).name);
	}
/*	public void setName(String string) {
		this.name =string;
		
	}*/
	public String id() {
		return id;
	}
public void setParam(int type,  double d){
	//if(true) throw new RuntimeException("!!");
	if(type==0){
		this.meanPrior = d;
		this.location = d;
	}
	else if(type==1){
		double d1 =Math.sqrt(Math.max(1e-10,d));
		this.stddevPrior = d1;
		this.scale = d1;
	}
	else{
		throw new RuntimeException("!!");
//		this.skewPrior[0] = d;
//		this.shape = d;
	}
//	this.recalcName();
	}

public void variance( double[] sum){
  	if(true)throw new RuntimeException("!!");
        
  }

public int shape() {
	throw new RuntimeException("!!");
}
public void setMinMax(double min, double max){
	//throw new RuntimeException("!!");
}
}

/*
 * # R package for the Skew-Normal (SN) and the skew-t (ST) distributions
# 
# Author: A.Azzalini 
# Home-page: http://azzalini.stat.unipd.it/SN
# major updates: 29/8/1997, 10/12/1997, 1/10/1998, 12/10/1998, 01/04/1999, 
# 15/06/2002, 01/04/2006
# It requires R 2.2.0, plus package mnormt for a few functions
#
#------- 

dsn <- function(x, location=0, scale=1, shape=0, dp=NULL, log=FALSE)
 {
   if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
   z <- (x-location)/scale
   if(log)
     y <- (-0.9189385332046727-logb(scale)-z^2/2+zeta(0,shape*z))
   else
     y <- 2*dnorm(z)*pnorm(z*shape)/scale
   replace(y, scale<= 0, NaN)
 }

psn <- function(x, location=0, scale=1, shape=0,  dp=NULL, ...)
 {
   if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
   z <- (x-location)/scale
   p <- pmin(1, pmax(0, pnorm(z) - 2*T.Owen(z, shape,...)))
   replace(p, scale<= 0, NaN)
 }
 
rsn <- function(n=1, location=0, scale=1, shape=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    }
  u1 <- rnorm(n)
  u2 <- rnorm(n)
  id <- (u2 > shape*u1)
  u1[id] <- (-u1[id])
  y <- location+scale*u1
  attr(y,"parameters") <- c(location,scale,shape)
  return(y)
  }

qsn <- function (p, location = 0, scale = 1, shape = 0, dp=NULL, 
           tol = 1e-08, ...) 
{   if(!is.null(dp)) {
      if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      }
    max.q <- sqrt(qchisq(p,1))
    min.q <- -sqrt(qchisq(1-p,1))
    if(shape == Inf) return(location + scale * max.q)
    if(shape == -Inf) return(location + scale * min.q)
    na <- is.na(p) | (p < 0) | (p > 1)
    zero <- (p == 0)
    one <- (p == 1)
    p <- replace(p, (na | zero | one), 0.5)
    cum <- sn.cumulants(0,1,shape, n=4)
    g1 <- cum[3]/cum[2]^(3/2)
    g2 <- cum[4]/cum[2]^2
    x <- qnorm(p)
    x <- (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
          x * (2 * x^2 - 5) * g1^2/36)
    x <- cum[1] + sqrt(cum[2]) * x
    max.err <- 1
    while (max.err > tol) {
        x1 <- x - (psn(x, 0, 1, shape, ...) - p)/dsn(x, 0, 1, shape)
        x1 <- pmin(x1,max.q)
        x1 <- pmax(x1,min.q)
        max.err <- max(abs(x1 - x)/(1 + abs(x)))
        x <- x1
    }
    x <- replace(x, na, NA)
    x <- replace(x, zero, -Inf)
    x <- replace(x, one, Inf)
    return(location + scale * x)
}

sn.cumulants <- function(location = 0, scale = 1, shape = 0, dp=NULL, n=4)
 {
   cumulants.half.norm <- function(n=4){
     n <- max(n,2)
     n <- as.integer(2*ceiling(n/2))
     half.n  <-  as.integer(n/2)
     m <- 0:(half.n-1)
     a <- sqrt(2/pi)/(gamma(m+1)*2^m*(2*m+1))
     signs <- rep(c(1,-1),half.n)[1:half.n]
     a <- as.vector(rbind(signs*a,rep(0,half.n)))
     coeff <- rep(a[1],n)
     for (k in 2:n) {
        ind <- 1:(k-1)
        coeff[k] <- a[k]-sum(ind*coeff[ind]*a[rev(ind)]/k)
        }
     kappa <- coeff*gamma((1:n)+1)
     kappa[2] <- 1+kappa[2]
     return(kappa)
    }
  if(!is.null(dp)) {
      if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      }
  par <- cbind(location,scale,shape)
  delta <- par[,3]/sqrt(1+par[,3]^2)
  kv <- cumulants.half.norm(n)
  if(length(kv)>n) kv<-kv[-(n+1)]
  kv[2] <- kv[2]-1
  kappa <- outer(delta,1:n,"^")*matrix(rep(kv,nrow(par)),ncol=n,byrow=TRUE)
  kappa[,2] <- kappa[,2]+1
  kappa <- kappa * outer(par[,2],(1:n),"^")
  kappa[,1] <- kappa[,1]+par[,1]
  kappa[,,drop=TRUE]
}

# lambda.of <- function(delta) delta/sqrt(1-delta^2)

# delta.of <- function(lambda) {
#   inf <- (abs(lambda)==Inf)
#   delta <-lambda/sqrt(1+lambda^2)
#   delta[inf] <- sign(lambda[inf])
#   delta
#}

T.Owen <- function(h, a, jmax=50, cut.point=6)
{
 T.int <-function(h,a,jmax,cut.point)
   {
     fui<- function(h,i) (h^(2*i))/((2^i)*gamma(i+1)) 
     seriesL <- seriesH <- NULL
     i  <- 0:jmax
     low<- (h<=cut.point)
     hL <- h[low]
     hH <- h[!low]
     L  <- length(hL)
     if (L>0) {
       b    <- outer(hL,i,fui)
       cumb <- apply(b,1,cumsum)
       b1   <- exp(-0.5*hL^2)*t(cumb)
       matr <- matrix(1,jmax+1,L)-t(b1)
       jk   <- rep(c(1,-1),jmax)[1:(jmax+1)]/(2*i+1)
       matr <- t(matr*jk) %*%  a^(2*i+1)
       seriesL  <- (atan(a)-as.vector(matr))/(2*pi)
     }
     if (length(hH) >0) 
       seriesH <- atan(a)*exp(-0.5*(hH^2)*a/atan(a))*
                  (1+0.00868*(hH^4)*a^4)/(2*pi)
     series <- c(seriesL,seriesH)
     id <- c((1:length(h))[low],(1:length(h))[!low]) 
     series[id] <- series  # re-sets in original order
     series
  }
  if(!is.vector(a) | length(a)>1) stop("a must be a vector of length 1")
  if(!is.vector(h)) stop("h must be a vector")
  aa <- abs(a)    
  ah <- abs(h)
  if(is.na(aa)) stop("parameter 'a' is NA") 
  if(aa==Inf) return(0.5*pnorm(-ah))
  if(aa==0)   return(rep(0,length(h)))
  na  <- is.na(h)
  inf <- (ah==Inf)
  ah  <- replace(ah,(na|inf),0)
  if(aa<=1)
    owen <- T.int(ah,aa,jmax,cut.point)
  else
    owen<-0.5*pnorm(ah)+pnorm(aa*ah)*(0.5-pnorm(ah))- 
               T.int(aa*ah,(1/aa),jmax,cut.point)
  owen <- replace(owen,na,NA)
  owen <- replace(owen,inf,0)
  return(owen*sign(a))
}


cp.to.dp <- function(param){
  # converts centred parameters cp=(mu,sigma,gamma1)
  # to direct parameters dp=(xi,omega,lambda)
  # Note:  mu can be m-dimensional, the other must be scalars
  b <- sqrt(2/pi)
  m <- length(param)-2
  gamma1 <- param[m+2]
  if(abs(gamma1)> 0.995271746431) stop("abs(gamma1)> 0.995271746431")
  A <- sign(gamma1)*(abs(2*gamma1/(4-pi)))^(1/3)
  delta <- A/(b*sqrt(1+A^2))
  lambda <- delta/sqrt(1-delta^2)
  E.Z  <- b*delta
  sd.Z <- sqrt(1-E.Z^2)
  location    <- param[1:m]
  location[1] <- param[1]-param[m+1]*E.Z/sd.Z
  scale <- param[m+1]/sd.Z
  dp    <- c(location,scale,lambda)
  names(dp)[(m+1):(m+2)] <- c("scale","shape")
  if(m==1)  names(dp)[1] <- "location"
  dp
  }

dp.to.cp <- function(param){
# converts 'direct' dp=(xi,omega,lambda) to 'centred' cp=(mu,sigma,gamma1)
  m <- length(param)-2
  omega <-param[m+1]
  lambda<-param[m+2]
  mu.Z  <- lambda*sqrt(2/(pi*(1+lambda^2)))
  s.Z   <- sqrt(1-mu.Z^2)
  gamma1<- 0.5*(4-pi)*(mu.Z/s.Z)^3
  sigma <- omega*s.Z
  mu    <- param[1:m]
  mu[1] <- param[1]+sigma*mu.Z/s.Z
  cp    <- c(mu,sigma,gamma1)
  names(cp)[(m+1):(m+2)]<-c("s.d.","skewness")
  if(m==1) names(cp)[1] <- "mean"
  cp
}

zeta <- function(k,x){# k integer in (0,4)
  if(k<0 | k>4 | k != round(k)) return(NULL)
  k <- round(k)
  na<- is.na(x)
  x <- replace(x,na,0)
  z <- switch(k+1,
            pnorm(x, log.p=TRUE)+ log(2),
            ifelse(x>(-20), dnorm(x)/pnorm(x), 
              ifelse(x>(-200), exp(-x^2/2-0.5*log(2*pi) - pnorm(x,log.p=TRUE)),
                               -x*(1+1/x^2-2/x^4))),              
            (-zeta(1,x)*(x+zeta(1,x))),
            (-zeta(2,x)*(x+zeta(1,x))-zeta(1,x)*(1+zeta(2,x))),
            (-zeta(3,x)*(x+2*zeta(1,x))-2*zeta(2,x)*(1+zeta(2,x))),
            NULL)
  neg.inf<- (x == -Inf)
  if(any(neg.inf))
    z <- switch(k+1,
                z,
                replace(z, neg.inf, Inf),
                replace(z, neg.inf, 1),
                replace(z, neg.inf, 0),
                replace(z, neg.inf, 0),
                NULL)
  if(k>1) z<- replace(z, x==Inf, 0)
  replace(z,na,NA)
}


sn.em <-function(X, y, fixed, p.eps=1e-4, l.eps=1.e-2, trace=FALSE, data=FALSE)
{
#
#  1/10/1998 (elaborando dal em.lm.sn del 2-12-97)
#
#  EM per caso con uno/due/tre parametri ignoti, parametrizzando in modo 
#  "diretta" con (xi, omega, lambda); internamente usa peraltro 'delta'.
#  Le componenti ignote sono i termini NA di fixed,  ma per semplicit�
#  assumiamo che un NA implica che le componenti alla sua sx sono NA
#  (e quindi il primo elemento di fixed � sempre NA).
#
  n <- length(y)
  if(missing(X)) X<-matrix(rep(1,n),n,1)
  nc <- ncol(X)
  if(missing(fixed)) fixed <- rep(NA,3)
  if(all(!is.na(fixed))) stop("all parameter are fixed") 
  if(is.na(fixed[3])) iter<-(1-log(l.eps,10)) else iter<-1 
  qrX<-qr(X)
  beta<-qr.coef(qrX,y)
  xi  <- m <-qr.fitted(qrX,y)
  omega  <- fixed[2] 
  lambda <- fixed[3]
  # delta  <- delta.of(lambda)
  delta <- lambda/sqrt(1+lambda^2)
  s<-sqrt(sum((y-xi)^2)/n)
  if(is.na(fixed[3])) {
    gamma1 <- sum((y-m)^3)/(n*s^3)
    a  <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
    delta<-sqrt(pi/2)*a/sqrt(1+a^2)
    if(abs(delta)>=1) delta<-sign(delta)/(1+1/n)
    # lambda<-lambda.of(delta)
    lambda<-delta/sqrt(1-delta^2)
    }
  mean.Z <- sqrt(2/pi)*delta
  sd.Z <- sqrt(1-mean.Z^2)
  if(is.na(fixed[2])) omega  <- s/sd.Z
  if(is.na(fixed[1])) xi     <- m-s*mean.Z/sd.Z
  old.par   <- c(beta,omega,lambda)
  diverge   <- 1
  incr.logL <- Inf
  logL      <- -Inf
  while(diverge>p.eps | incr.logL>l.eps){
    # E-step
    v  <- (y-xi)/omega
    p  <- zeta(1,lambda*v)
    u1 <- omega*(delta*v+p*sqrt(1-delta^2))
    u2<-omega^2*((delta*v)^2+(1-delta^2)+p*v*delta*sqrt(1-delta^2))
    # M-step
    for(i in 1:iter){
      beta<-qr.coef(qrX,y-delta*u1)
      xi <- qr.fitted(qrX,y-delta*u1)
      d  <- y-xi
      Q  <- sum(d^2-2*delta*d*u1+u2)
      if(is.na(fixed[2])) omega <-sqrt(Q/(2*n*(1-delta^2)))
      r  <- 2*sum(d*u1)/Q
      if(is.na(fixed[3])) delta<-(sqrt((2*r)^2+1)-1)/(2*r)
      }
    # convergence?    # lambda<-lambda.of(delta)
    lambda <- delta/sqrt(1-delta^2)
    param <- c(beta,omega,lambda)
    names(param)[(nc+1):(nc+2)] <-c("scale","shape")
    if(nc==1 & all(X==1)) names(param)[1] <- "location"
    else  names(param)[1:nc] <- colnames(X)
    diverge<-sum(abs(param-old.par)/(1+abs(old.par)))/(nc+2)
    old.par<-param
    new.logL <- sum(dsn(y,xi,omega,lambda, log=TRUE))
    incr.logL<- new.logL-logL
    logL <- new.logL
    if(trace) print(c(param,logL),digits=5)
  }
  cp <- dp.to.cp(param)
  result<-list(dp=param, cp=cp, logL=logL)
  if(data) result$data <- list(X=X, y=y, residuals=d/omega)
  result
}

#-------------------

gamma1.to.lambda<- function(gamma1){
  max.gamma1 <- 0.5*(4-pi)*(2/(pi-2))^1.5
  na <- (abs(gamma1)>max.gamma1)
  if(any(na)) warning("NAs generated") 
  gamma1<-replace(gamma1,na,NA)
  a    <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta<- sqrt(pi/2)*a/sqrt(1+a^2)
  lambda<-delta/sqrt(1-delta^2)
  as.vector(lambda)
}

# Examples:
#  a<-sn.2logL.profile(y=otis)
#  a<-sn.2logL.profile(y=otis, use.cp=FALSE)
#  a<-sn.2logL.profile(X=cbind(1,lbm), y=bmi, npts=50)
#  a<-sn.2logL.profile(y=frontier,param.range=c(0.8,1.6,10,30), use.cp=FALSE, npts=11)

sn.2logL.profile<-function(X=matrix(rep(1,n)), y, 
      param.range=c(sqrt(var(y))*c(2/3, 3/2), -0.95, 0.95),
      use.cp=TRUE, npts= 31 %/% d, plot.it=TRUE, ...)
{# plot 1D or 2D profile deviance (=-2logL) using either parameters
   # if(plot.it & !exists(.Device)) stop("Device not active")
   n<-length(y)
   d<- round(length(param.range)/2)
   if((d!=1)&(d!=2)) stop(" length(param.range) must be either 2 or 4")
   if(d==1){
      param1 <- seq(param.range[1],param.range[2],length=npts)
      llik <- param2 <- rep(NA,npts)}
   else{ 
      param1 <- seq(param.range[1],param.range[2],length=npts)
      param2 <- seq(param.range[3],param.range[4],length=npts)
      llik   <- matrix(NA,npts,npts)}
   if(use.cp){
      if(d==1){
        gamma1<-param1
        sigma <-param2 
        xlab <- "gamma1"
        ylab <- ""}
      else {
        sigma <-param1
        gamma1<-param2
        xlab <- "sigma"
        ylab <- "gamma1"
        }
      if(max(abs(gamma1))>0.9952717) stop("abs(gamma1)>0.9952717")
      lambda <- gamma1.to.lambda(gamma1)
      sc<-sqrt(1-(2/pi)*lambda^2/(1+lambda^2))      
      }
   else{ # use dp 
      if(d==1) {
        lambda<-param1
        omega<-param2
        xlab <- "alpha"
        ylab <- ""}
      else {
         omega<-param1 
     sc <- rep(1,npts)
         lambda<-param2
         xlab <- "omega"
         ylab <- "alpha"
         }
      }
   cat(c("Running until",npts,":"))
   for(i in 1:npts){
     cat(" ");cat(i)
     if(d==1) {
       a<-sn.em(X, y, fixed=c(NA,NA,lambda[i]), ...)       
       llik[i]<-a$logL
       # print(c(i,lambda[i],a$logL))
       }
     else{
     for(j in 1:npts){     
       a<-sn.em(X, y, fixed=c(NA,param1[i]/sc[j],lambda[j]), ...)
       llik[i,j] <- a$logL
       # print( c(i,j, param1[i]/sc[j], lambda[j], a$logL))
     }}
   }
  cat("\n")
  #if(plot)
  f<-2*(llik-max(llik)) 
  if(plot.it){
    if(d==1) plot(param1, f, type="l", 
            xlab=xlab, ylab="profile relative 2(logL)")
    else contour(param1, param2, f, labcex=0.75, 
            xlab=xlab, ylab=ylab,
            levels=-c(0.575, 1.386, 2.773, 4.605, 5.991, 9.210))
            # qchisq(c(0.25,0.5,0.75,0.90,0.95,0.99),2)
    title(main=paste("dataset:", deparse(substitute(y)),
        "\nProfile relative 2(logLikelihood)", sep= " "))   
  }
  invisible( list(param1=param1, param2=param2,
      param.names=c(xlab,ylab), two.logL=f, maximum=2*max(llik)))
}

  
sn.mle <- function(X, y, cp, plot.it=TRUE, trace=FALSE, method="L-BFGS-B",
               control=list(iter.max=100, abs.tol=1e-5)) 
{
  xlab<-deparse(substitute(y))
  if(!is.null(dim(y))) {
    if(min(dim(y))==1) y<-as.vector(y)
    else stop("y must be a vector")
    }
  n<-length(y)
  if(missing(X)) {
     X <-as.matrix(rep(1,n))
     cp.names <- "mean"
   }
  else{
    if(is.null(colnames(X)))
       cp.names<-  outer(deparse(substitute(X)),as.character(1:ncol(X)),
                         paste, sep=".")
    else  cp.names<- colnames(X)
     }
  cp.names<- c(cp.names,"s.d.","skewness")
  m<-ncol(X)
  if(missing(cp)) {
    qrX <- qr(X)
    s <- sqrt(sum(qr.resid(qrX, y)^2)/n)
    gamma1 <- sum(qr.resid(qrX, y)^3)/(n*s^3)
    if(abs(gamma1) > 0.99527) gamma1<- sign(gamma1)*0.95
    cp <- c(qr.coef(qrX,y), s, gamma1)
    }
  else{ 
    if(length(cp)!= (m+2)) stop("ncol(X)+2 != length(cp)")}
  opt<- optim(cp, fn=sn.dev, gr=sn.dev.gh, method=method,
          lower=c(-rep(Inf,m), 10*.Machine$double.eps, -0.99527), 
          upper=c(rep(Inf,m), Inf, 0.99527), 
          control=control, X=X, y=y, trace=trace, hessian=FALSE)
  cp <- opt$par
  if(trace) {
     cat(paste("Message from optimization routine:", opt$message,"\n"))
     cat("estimates (cp): ", cp, "\n")
     }
  if(abs(cp[m+2])> 0.9952717){
     if(trace) cat("optim searched outside admissible range - restarted\n")
      cp[m+2]<- sign(cp[m+2])*runif(1)
      mle <- sn.mle(X, y, cp, plot.it, trace, method, control)
      cp  <- mle$cp
      }
  logL <- (-opt$value)/2
  info <- attr(sn.dev.gh(cp, X, y, trace=FALSE, hessian=TRUE),"hessian")/2
  # se <- sqrt(diag(solve(info)))
  if(all(is.finite(info))) 
    {
      qr.info <- qr(info)
      info.ok <- (qr.info$rank == length(cp))
     }
  else info.ok <- FALSE
  if(info.ok) {
    se2 <- diag(solve.qr(qr.info))
    se <- sqrt(ifelse(se2 >= 0, se2, NA))
    }
  else
    se <- rep(NA, length(cp))
  if(plot.it) {
    dp0<-cp.to.dp(cp)
    if(all(X==rep(1,n))) 
      y0<-y        
    else {
      y0<- as.vector(y - X %*% dp0[1:m])
      dp0<-c(0,dp0[m+1],dp0[m+2])
      xlab<-"residuals"
      }
    x<-seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
    pdf.sn <- dsn(x,dp0[1],dp0[2],dp0[3])
    if(exists("sm.density",mode="function"))
      {
      a<-sm.density(x=y0, eval.points=x, h=hnorm(y0)/1.5, display="none")
      a<-sm.density(x=y0, eval.points=x, h=hnorm(y0)/1.5, xlab=xlab, 
                    lty=3, ylim=c(0,max(a$estimate,pdf.sn)))
      }
    else 
      {
      h <- hist(y0, prob=TRUE, breaks="FD", plot=FALSE)
      hist(y0, prob=TRUE, breaks="FD", xlim=c(min(x),max(x)),
           xlab=xlab, main=xlab, ylim=c(0, max(pdf.sn, h$density)))
      }
    if(n<101) points(y0,rep(0,n),pch=1)
    # title(deparse(substitute(y)))
    curve(dsn(x, dp0[1], dp0[2], dp0[3]), add=TRUE, col=2)
  }
  names(cp)<-  names(se)<- cp.names
  list(call=match.call(), cp=cp,  se=se, info=info, logL=logL, optim=opt)
}


sn.dev <- function(cp, X, y, trace=FALSE)
{ # -2*logL for centred parameters  
  m <- ncol(X)
  if(abs(cp[m+2])> 0.9952717){
    warning("optim search in abs(cp[m+2])> 0.9952717, value adjusted")
    cp[m+2] <- 0.9952717*sign(cp[m+2])
  }
  dp <- as.vector(cp.to.dp(cp))
  location <- as.vector(X %*% as.matrix(dp[1:m]))
  logL <- sum(dsn(y, location, dp[m+1], dp[m+2], log=TRUE))
  if(trace) {cat("sn.dev: (cp,dev)="); print(c(cp,-2*logL))}
  return(-2*logL)
}

sn.dev.gh <- function(cp, X, y, trace=FALSE, hessian=FALSE)
{
  # computes gradient and hessian of dev=-2*logL for centred parameters 
  # (and observed information matrix);
  m  <- ncol(X)
  n  <- nrow(X)
  np <- m+2
  if(abs(cp[m+2])> 0.9952717){
    warning("optim search in abs(cp[m+2])> 0.9952717, value adjusted")
    cp[m+2] <- 0.9952717*sign(cp[m+2])
  }
  score <- rep(NA,np)
  info  <- matrix(NA,np,np)
  beta  <- cp[1:m]
  sigma <- cp[m+1]
  gamma1 <- cp[m+2]
  lambda <- gamma1.to.lambda(gamma1)
  # dp<-cp.to.dp(c(beta,sigma,gamma1))
  # info.dp <- sn.info(dp,y)$info.dp
  mu <- as.vector(X %*% as.matrix(beta))
  d  <- y-mu
  r  <- d/sigma
  E.Z<- lambda*sqrt(2/(pi*(1+lambda^2)))
  s.Z<- sqrt(1-E.Z^2)
  z  <- E.Z+s.Z*r
  p1 <- as.vector(zeta(1,lambda*z))
  p2 <- as.vector(zeta(2,lambda*z))
  omega<- sigma/s.Z
  w    <- lambda*p1-E.Z
  DE.Z <- sqrt(2/pi)/(1+lambda^2)^1.5
  Ds.Z <- (-E.Z/s.Z)*DE.Z
  Dz   <- DE.Z + r*Ds.Z
  DDE.Z<- (-3)*E.Z/(1+lambda^2)^2
  DDs.Z<- -((DE.Z*s.Z-E.Z*Ds.Z)*DE.Z/s.Z^2+E.Z*DDE.Z/s.Z)
  DDz  <- DDE.Z + r*DDs.Z
  score[1:m] <- omega^(-2)*t(X) %*% as.matrix(y-mu-omega*w) 
  score[m+1] <- (-n)/sigma+s.Z*sum(d*(z-p1*lambda))/sigma^2
  score[m+2] <- score.l <- n*Ds.Z/s.Z-sum(z*Dz)+sum(p1*(z+lambda*Dz))
  Dg.Dl <-1.5*(4-pi)*E.Z^2*(DE.Z*s.Z-E.Z*Ds.Z)/s.Z^4
  R <- E.Z/s.Z
  T <- sqrt(2/pi-(1-2/pi)*R^2)
  Dl.Dg <- 2*(T/(T*R)^2+(1-2/pi)/T^3)/(3*(4-pi))
  R. <- 2/(3*R^2 * (4-pi))
  T. <- (-R)*R.*(1-2/pi)/T
  DDl.Dg <- (-2/(3*(4-pi))) * (T./(R*T)^2+2*R./(T*R^3)+3*(1-2/pi)*T./T^4)
  score[m+2] <- score[m+2]/Dg.Dl  # convert deriv wrt lamda to gamma1 
  gradient <- (-2)*score
  if(hessian){
     info[1:m,1:m] <- omega^(-2) * t(X) %*% ((1-lambda^2*p2)*X)
     info[1:m,m+1] <- info[m+1,1:m] <- 
            s.Z* t(X) %*% as.matrix((z-lambda*p1)+d*(1-lambda^2*p2)*
            s.Z/sigma)/sigma^2
     info[m+1,m+1] <- (-n)/sigma^2+2*s.Z*sum(d*(z-lambda*p1))/sigma^3 +
            s.Z^2*sum(d*(1-lambda^2*p2)*d)/sigma^4
     info[1:m,m+2] <- info[m+2,1:m] <- 
            t(X)%*%(-2*Ds.Z*d/omega+Ds.Z*w+s.Z*(p1+lambda*p2*(z+lambda*Dz)
            -DE.Z))/sigma 
     info[m+1,m+2] <- info[m+2,m+1] <- 
            -sum(d*(Ds.Z*(z-lambda*p1)+s.Z*(Dz-p1-p2*lambda*(z+lambda*Dz))
             ))/sigma^2
     info[m+2,m+2] <- (n*(-DDs.Z*s.Z+Ds.Z^2)/s.Z^2+sum(Dz^2+z*DDz)-
            sum(p2*(z+lambda*Dz)^2)- sum(p1*(2*Dz+lambda*DDz)))
     info[np,] <- info[np,]/Dg.Dl # convert info wrt lamda to gamma1 
     info[,np] <- info[,np]*Dl.Dg # an equivalent form of the above
     info[np,np] <- info[np,np]-score.l*DDl.Dg
     attr(gradient,"hessian") <- 2*info
     }
  if(trace) {cat("sn.dev.gh: gradient="); print(-2*score)}
  return(gradient)
}

# 

dmsn <- function(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, log=FALSE)
{
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else {if(!is.null(dp$beta)) xi <- as.vector(dp$beta)}
      Omega <- dp$Omega
      alpha <- dp$alpha
    }
    d <- length(alpha)
    Omega <- matrix(Omega,d,d)
    x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x) 
    if(is.vector(xi)) xi <- outer(rep(1,nrow(x)), xi)
    X <- t(x - xi)
    Q <- apply((solvePD(Omega) %*% X) * X, 2, sum)
    L <- as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
    logDet <- sum(logb(abs(diag(qr(Omega)$qr))))
    logPDF <- (logb(2) - 0.5 * Q + pnorm(L, log = TRUE)
               - 0.5 * (d * logb(2 * pi) + logDet))
    if (log) logPDF
    else exp(logPDF)
}


rmsn <- function(n=1, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL)
{# generates SN_d(xi,Omega,alpha) variates using transformation method
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
    }
  d <- length(alpha)
  Z <- msn.quantities(xi,Omega,alpha)
  y <- matrix(rnorm(n*d),n,d) %*% chol(Z$Psi)
  # each row of y is N_d(0,Psi)
  abs.y0 <- abs(rnorm(n))  
  abs.y0 <- matrix(rep(abs.y0,d), ncol=d)
  delta  <- Z$delta
  z <- delta * t(abs.y0) + sqrt(1-delta^2) * t(y)
  y <- t(xi+Z$omega*z)
  attr(y,"parameters") <- list(xi=xi,Omega=Omega,alpha=alpha)
  return(y)
}

pmsn <- function(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, ...)
{
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
  if(!is.null(dp$xi)) xi <- dp$xi
     else
  if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
     Omega <- dp$Omega
     alpha <- dp$alpha
  }   
  pmst(x, xi=xi, Omega=Omega, alpha=alpha, df=Inf, ...)
}


dsn2.plot <- function(x, y, xi, Omega, alpha, dp=NULL, ...)
{# plot bivariate density SN_2(xi,Omega,alpha) computed at (x,y) grid
  if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
      }  
  if(any(dim(Omega)!=c(2,2))) stop("dim(Omega) != c(2,2)")
  nx <- length(x)
  ny <- length(y)
  xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=TRUE)))
  X <- matrix(xoy, nx*ny, 2, byrow=FALSE)
  pdf<-dmsn(X, xi, Omega, alpha)
  pdf<-matrix(pdf, nx, ny)
  contour(x, y, pdf, ...)
  invisible(list(x=x, y=y, density=pdf, xi=xi, Omega=Omega, alpha=alpha))
}

msn.quantities <- function(xi=rep(0,length(alpha)), Omega, alpha, dp=NULL)
{# 21-12/1997; computes various quantities related to SN_d(xi,Omega,alpha)
  if(!is.null(dp)){
    if(any(!missing(xi) | !missing(Omega) | !missing(alpha))) 
      stop("you cat set either dp or its components, but not both")
    xi<- as.vector(dp$xi)
    if(is.null(dp$xi)) xi<- as.vector(dp$beta)
    Omega<- dp$Omega
    alpha<- dp$alpha
   }
  d <- length(alpha)
  Omega<- as.matrix(Omega)
  if(length(xi)!=d | any(dim(Omega)!=c(d,d))) 
       stop("dimensions of arguments do not match")
  omega <- sqrt(diag(Omega))
  O.cor <- cov2cor(Omega)
  tmp <- as.vector(sqrt(1 + t(as.matrix(alpha))%*%O.cor%*%alpha)) 
  delta<- as.vector(O.cor %*%alpha)/tmp
  lambda<- delta/sqrt(1-delta^2)
  D <-diag(sqrt(1+lambda^2))
  Psi <- D %*% (O.cor-outer(delta,delta)) %*% D
  Psi <- (Psi+t(Psi))/2
  O.inv <- solvePD(Omega)
  O.pcor <- -cov2cor(O.inv) 
  O.pcor[cbind(1:d, 1:d)] <- 1
  muZ <- delta*sqrt(2/pi)
  muY <- xi+omega*muZ
  Sigma <- diag(omega,d,d) %*% (O.cor-outer(muZ,muZ)) %*% diag(omega,d,d) 
  Sigma <- (Sigma+t(Sigma))/2
  cv <- muZ/sqrt(1-muZ^2)
  gamma1 <- 0.5*(4-pi)*cv^3
  list(xi=xi, Omega=Omega, alpha=alpha, omega=omega,  mean=muY, variance=Sigma,
       Omega.conc=O.inv, Omega.cor=O.cor, Omega.pcor=O.pcor,
       lambda=lambda, Psi=Psi,  delta=delta, skewness=gamma1)
}

msn.conditional <- function(xi=rep(0,length(alpha)), Omega, alpha,
                            fixed.comp, fixed.values, dp=NULL)
{
# conditional Multivariate SN  (6/11/1997).
# Given a rv Y ~ SN_d(xi,Omega,alpha), this function computes cumulants of
# conditrional distribution, given that the fixed.com take on fixed.values;
# then it finds MSN with matching cumulants.  
  Diag <- function(x) diag(x,nrow=length(x),ncol=length(x))
  msqrt <- function(A) Diag(sqrt(diag(as.matrix(A))))
  imsqrt<- function(A) Diag(1/sqrt(diag(as.matrix(A))))
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi <- dp$xi
    else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
    Omega <- dp$Omega
    alpha <- dp$alpha
  }
  d <- length(alpha) 
  h <- length(fixed.comp)
  if(any(dim(Omega)!=c(d,d)) | length(xi)!=d | h!=length(fixed.values))
        stop("dimensions of parameters do not match")
  fc <- fixed.comp
  O  <- as.matrix(Omega)
  O11<- O[fc,fc, drop=FALSE]
  O12<- O[fc,-fc, drop=FALSE]
  O21<- O[-fc,fc, drop=FALSE]
  O22<- O[-fc,-fc, drop=FALSE]
  o22<- sqrt(diag(O22))
  inv.O11 <- solvePD(O11)
  xi1 <- xi[fc, drop=FALSE]
  xi2 <- xi[-fc, drop=FALSE]
  alpha1 <- matrix(alpha[fc])
  alpha2 <- matrix(alpha[-fc]) 
  O22.1 <- O22 - O21 %*% inv.O11 %*% O12
  O22.b <- imsqrt(O22) %*% O22.1 %*% imsqrt(O22)
  xi.c  <- xi2 + as.vector(O21 %*% inv.O11 %*% matrix(fixed.values-xi1))
  a     <- sqrt(1+as.vector(t(alpha2) %*% O22.b %*% alpha2))
  alpha.b <- (alpha1 + msqrt(O11) %*% inv.O11 %*% O12 %*% (alpha2/o22))/a  
  d2    <- as.vector(O22.b %*% alpha2)/a
  x0    <- sum(alpha.b * (fixed.values-xi1)/sqrt(diag(O11)))
  E.c   <- xi.c + zeta(1,x0)*o22*d2
  var.c <- O22.1+zeta(2,x0)*outer(o22*d2,o22*d2)
  gamma1<- zeta(3,x0)*d2^3/diag(O22.b+zeta(2,x0)*d2^2)^1.5
  cum   <- list(as.vector(E.c),var.c,gamma1)
  # cumulants are computed; now choose SN distn to fit them
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  E.z   <- a/sqrt(1+a^2)
  delta <- E.z*sqrt(pi/2) 
  omega <- sqrt(diag(var.c)/(1-E.z^2))
  O.new <- var.c+outer(omega*E.z,omega*E.z) 
  xi.new<- E.c-omega*E.z
  B   <- diag(1/omega,d-h,d-h)
  m   <- as.vector(solvePD(B %*% O.new %*% B) %*% as.matrix(delta))
  a   <- m/sqrt(1-sum(delta*m))
  # cum2<- msn.cumulants(xi.new,O.new,a)
  list(cumulants=cum, fit=list(xi=xi.new, Omega=O.new, alpha=a, delta=delta))
}


msn.marginal <- function(xi=rep(0,length(alpha)), Omega, alpha,
                         comp=1:d, dp=NULL)
                         
{# calcola parametri della marginale associata a comp di un SN_d 
 # cfr SJS 2003, p.131-2
  if(!is.null(dp)){
    if(!is.null(dp$xi)) xi <- dp$xi
    else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
    Omega <- dp$Omega
    alpha <- dp$alpha
    }
  alpha <- as.vector(alpha)
  d <- length(alpha)
  xi <- as.vector(xi)
  comp <- as.integer(comp)
  if(length(xi) != d) stop("parameter size not compatible")
  if(all(dim(Omega) != c(d,d))) stop("parameter size not compatible")
  if(length(comp)<d){
    if(any(comp>d | comp<1)) stop("comp makes no sense")
    O   <- cov2cor(Omega)
    O11 <- O[comp,comp, drop=FALSE]
    O12 <- O[comp,-comp, drop=FALSE]
    O21 <- O[-comp,comp, drop=FALSE]
    O22 <- O[-comp,-comp, drop=FALSE]
    alpha1 <- matrix(alpha[comp], ncol=1)
    alpha2 <- matrix(alpha[-comp], ncol=1)
    O11_inv <- solvePD(O11)
    O22.1 <- O22 - O21 %*% O11_inv %*% O12
    a.sum <- as.vector(t(alpha2) %*% O22.1 %*% alpha2)
    a.new <- as.vector(alpha1 + O11_inv %*% O12 %*% alpha2)/sqrt(1+a.sum)
    result<- list(xi=xi[comp], Omega=Omega[comp,comp], alpha=a.new)
  }
  else {
    if(any(sort(comp)!=(1:d))) stop("comp makes no sense")
    result <- list(xi=xi[comp], Omega=Omega[comp,comp], alpha=alpha[comp])
  }
  if(!is.null(dp$tau)) result$tau <- dp$tau
  result
}




msn.cond.plot <- function(xi, Omega, alpha, fixed.comp, fixed.values, n=35)
{# fa il grafico di Y_2|Y_1; assumiamo che dim(Y_2)= 2 
  msn.pdf2.aux <- function(x,y,xi,Omega,alpha,fc,fv)
   {
     nx <- length(x)
     ny <- length(y)
     FV <- matrix(rep(fv,nx*ny), nx*ny, length(fv), byrow=TRUE)
     X <- matrix(NA, nx*ny, length(alpha))
     X[,fc] <- FV
     xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=TRUE)))
     X[,-fc] <- matrix(xoy, nx*ny, 2, byrow=FALSE)
     pdf<-dmsn(X,xi,Omega,alpha)
     matrix(pdf,nx,ny)
   }
  dsn2 <- function(x,y,d1,d2,omega)
    {
     u <- (x*(d1-omega*d2)+y*(d2-omega*d1))/
          sqrt((1-omega^2-d1^2-d2^2+2*omega*d1*d2)*(1-omega^2))
     pdfn2 <- exp(-0.5*(x^2-2*omega*x*y+y^2)/(1-omega^2))/
              (2*pi*sqrt(1-omega^2))
     2*pdfn2*pnorm(u)
    }
  fc <- fixed.comp
  fv <- fixed.values
  cond <- msn.conditional(xi,Omega,alpha,fc,fv)
  xi.c <- cond$fit$xi
  O.c  <- cond$fit$Omega
  a.c  <- cond$fit$alpha
  if(any(dim(O.c)!=c(2,2))) stop("length(alpha)-length(fixed.com)!=2")
  scale1<-sqrt(as.vector(O.c[1,1]))
  scale2<-sqrt(as.vector(O.c[2,2]))
  delta <- cond$fit$delta
  omega <-as.vector(O.c[1,2])/(scale1*scale2)
  x<-seq(xi.c[1]-3*scale1, xi.c[1]+3*scale1, length=n)
  y<-seq(xi.c[2]-3*scale2, xi.c[2]+3*scale2, length=n)
  plot(x,y,type="n", main="Conditional multivariate SN pdf")
  z1<-(x-xi.c[1])/scale1
  z2<-(y-xi.c[2])/scale2
  pdf.fit<-outer(z1,z2,dsn2,d1=delta[1],d2=delta[2],omega=omega)/
                (scale1*scale2)
  cond$pdf<-list(x=x,y=y,f.fitted=pdf.fit)
  levels<-pretty(pdf.fit,5)
  contour(x,y,pdf.fit,levels=levels,add=TRUE,col=2)
  # fino a qui per il calcolo della densit� approx;
  # ora otteniamo quella esatta
  numer <- msn.pdf2.aux(x, y, xi, Omega, alpha, fc, fv)
  marg  <- msn.marginal(xi, Omega, alpha, fc)
  denom <- dmsn(fv, marg$xi, marg$Omega, marg$alpha)
  pdf.exact<- numer/as.vector(denom)
  contour(x, y, pdf.exact, add=TRUE, levels=levels, col=1, lty=4, labcex=0)
  legend(x[1],y[length(y)],c("approx","exact"), lty=c(1,4),col=c(2,1))
  cond$pdf$f.exact<-pdf.exact
  cond$rel.error<-summary((pdf.fit-pdf.exact)/pdf.exact)
  cond$abs.error<-summary(abs(pdf.fit-pdf.exact))
  invisible(cond)
}


msn.moment.fit <- function(y)
{# 31-12-1997: simple fit of MSN distribution usign moments
  y     <- as.matrix(y)
  k     <- ncol(y)
  m.y   <- apply(y,2,mean)
  var.y <- var(y)
  y0    <- (t(y)-m.y)/sqrt(diag(var.y))
  gamma1<- apply(y0^3,1,mean)
  out   <- (abs(gamma1)>0.99527)
  gamma1[out] <- sign(gamma1[out])*0.995
  a     <- sign(gamma1)*(2*abs(gamma1)/(4-pi))^0.33333
  delta <- sqrt(pi/2)*a/sqrt(1+a^2)
  m.z   <- delta*sqrt(2/pi) 
  omega <- sqrt(diag(var.y)/(1-m.z^2))
  Omega <- var.y+outer(omega*m.z,omega*m.z) 
  xi    <- m.y-omega*m.z
  O.cor <- cov2cor(Omega)
  O.inv <- solvePD(O.cor)
  tmp   <- as.vector(1-t(delta) %*% O.inv %*% delta)
  if(tmp<=0) {tmp <- 0.0001; admissible <- FALSE} 
        else admissible<-TRUE
  alpha <-as.vector(O.inv%*%delta)/sqrt(tmp)
  list(xi=xi, Omega=Omega, alpha=alpha, Omega.cor=O.cor, omega=omega, 
       delta=delta, skewness=gamma1, admissible=admissible) 
}

msn.fit <- function(X, y, freq, plot.it=TRUE, trace=FALSE, ... )
{
  y.name <- deparse(substitute(y))
  y.names<- dimnames(y)[[2]] 
  y <- as.matrix(y)
  colnames(y)<-y.names
  k <- ncol(y)
  if(missing(freq)) freq<-rep(1,nrow(y))  
  n <- sum(freq)
  if(missing(X)) {
    X <- rep(1,nrow(y))
    missing.X <- TRUE }
  else
    missing.X <- FALSE
  X <- as.matrix(X)
  m <- ncol(X)
  if(length(dimnames(y)[[2]])==0) {
      dimnames(y) <- list(NULL, outer("V",as.character(1:k),paste,sep=""))
      y.names<- as.vector(dimnames(y)[[2]])
      }
  qrX <- qr(X)
  mle<- msn.mle(X=X, y=y, freq=freq, trace=trace, ...)
  mle$call <- match.call()
  # print(mle$nlminb$message)
  beta  <- mle$dp$beta
  Omega <- mle$dp$Omega
  alpha <- mle$dp$alpha
  omega <- sqrt(diag(Omega))
  xi    <- X %*% beta
  if(plot.it & all(freq==rep(1,length(y)))) {
    if(missing.X) { 
      y0  <-y 
      xi0 <- apply(xi,2,mean)} 
    else  {
      y0  <- y-xi 
      xi0 <- rep(0,k)
      }
    if(k>1) {
       opt<-options()
       options(warn=-1)
       pairs(y0, labels=y.names,
        panel=function(x, y, Y, y.names, xi, Omega, alpha)
         {
          for(i in 1:length(alpha)){
            # if(y.names[i]==deparse(substitute(x))) Ix<-i
            # if(y.names[i]==deparse(substitute(y))) Iy<-i
            if(all(Y[,i]==x)) Ix<-i
            if(all(Y[,i]==y)) Iy<-i
            }
          points(x,y)
          marg<-msn.marginal(xi,Omega,alpha,c(Ix,Iy))
          xi.marg<-marg$xi
          Omega.marg<-marg$Omega
          alpha.marg<-marg$alpha     
          x1 <- seq(min(x), max(x), length=30)
          x2 <- seq(min(y), max(y), length=30)
          dsn2.plot(x1, x2, xi.marg, Omega.marg, alpha.marg, add=TRUE, col=2)},  
          # end "panel" function
      Y=y0, y.names=y.names, xi=xi0, Omega=Omega, alpha=alpha)
      options(opt) 
      }
    else{ # plot for case k=1
      y0<-as.vector(y0)
      x<-seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
      if(missing.X){
         dp0<-c(xi0,omega,alpha)
         xlab<-y.name}
      else {
         dp0<-c(0,omega,alpha)
         xlab <- "residuals"}
      hist(y0, prob=TRUE, breaks="FD", xlab=xlab, ylab="density")
      lines(x, dsn(x,dp0[1],dp0[2],dp0[3]))
      if(length(y)<101) points(y0, rep(0,n), pch=1)
      title(y.name)
      }
    cat("Press <Enter> to continue..."); readline()
    par(mfrow=c(1,2))
    pp <- qchisq((1:n)/(n+1),k)
    # Xb <- qr.fitted(qrX,y)
    res<- qr.resid(qrX,y)
    rad.n  <- apply(res    * (res    %*% solvePD(var(res))), 1, sum)
    rad.sn <- apply((y-xi) * ((y-xi) %*% solvePD(Omega)), 1, sum)
    plot(pp, sort(rad.n), pch=1, ylim=c(0,max(rad.n,rad.sn)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for normal distribution", sub=y.name)
    plot(pp, sort(rad.sn), pch=1, ylim=c(0,max(rad.n,rad.sn)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for skew-normal distribution", sub=y.name)
    cat("Press <Enter> to continue..."); readline()
    plot((1:n)/(n+1), sort(pchisq(rad.n,k)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for normal distribution", sub=y.name)
    plot((1:n)/(n+1), sort(pchisq(rad.sn,k)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for skew-normal distribution", sub=y.name)
    par(mfrow=c(1,1))
    } # end ploting
  dev.norm<- msn.dev(c(qr.coef(qrX,y),rep(0,k)), X, y, freq)
  test <- dev.norm + 2*mle$logL
  p.value <-  1-pchisq(test,k)
  if(trace) {
    cat("LRT for normality (test-function, p-value): ")
    print(c(test,p.value))
    }
  mle$test.normality <- list(LRT=test, p.value=p.value)
  invisible(mle)
}

msn.mle <-function(X, y, freq, start, trace=FALSE, 
                algorithm=c("nlminb", "Nelder-Mead", "BFGS", "CG",  "SANN"),
                control=list() )
{
  y <- data.matrix(y)
  if(missing(X)) X <- rep(1,nrow(y))
    else {if(!is.numeric(X)) stop("X must be numeric")}
  if(missing(freq)) freq <- rep(1,nrow(y))
  algorithm <- match.arg(algorithm)
  X <- data.matrix(X) 
  k <- ncol(y)  
  n <- sum(freq)
  m <- ncol(X)
  y.names<-dimnames(y)[[2]] 
  x.names<-dimnames(X)[[2]]
  if(missing(start)) {
      fit0  <- lm.fit(X, y, method="qr")
      beta  <- as.matrix(coef(fit0))
      res   <- resid(fit0)
      a     <- msn.moment.fit(res)
      Omega <- a$Omega
      omega <- a$omega
      alpha <- a$alpha
      if(!a$admissible) alpha<-alpha/(1+max(abs(alpha)))
      beta[1,] <- beta[1,]-omega*a$delta*sqrt(2/pi)  
     }
  else{
    beta  <- start$beta
    Omega <- start$Omega
    alpha <- start$alpha
    omega <- sqrt(diag(Omega)) 
    }
  al.om <-alpha/omega
  if(trace){ 
    cat("Initial parameters:\n")
    print(cbind(t(beta),al.om,Omega))
    }
  param<- c(beta,al.om)
  dev <- msn.dev(param,X,y,freq)    
  if(algorithm == "nlminb"){
    opt <- nlminb(param, msn.dev, msn.dev.grad, 
                control=control, X=X, y=y, freq=freq, trace=trace)      
    opt$value<- opt$objective 
    }
  else{
   opt <- optim(param, fn=msn.dev, gr=msn.dev.grad, method=algorithm,
               control=control, X=X, y=y, freq=freq, trace=trace)    
   }
  if(trace) 
    cat(paste("Message from optimization routine:", opt$message,"\n"))
  logL <- opt$value/(-2) 
  beta <- matrix(opt$par[1:(m*k)],m,k)
  dimnames(beta)[2] <- list(y.names)
  dimnames(beta)[1] <- list(x.names)
  al.om  <- opt$par[(m*k+1):(m*k+k)]
  xi    <- X %*% beta
  Omega <- t(y-xi) %*% (freq*(y-xi))/n
  omega <- sqrt(diag(Omega))
  alpha <- al.om*omega
  param <- cbind(omega,alpha)
  dimnames(Omega) <- list(y.names,y.names)
  dimnames(param)[1] <- list(y.names)
  info <- num.deriv2(opt$par, FUN="msn.dev.grad", X=X, y=y, freq=freq)/2
    if (all(is.finite(info))) {
        qr.info <- qr(info)
        info.ok <- (qr.info$rank == length(param))
    }
    else info.ok <- FALSE
    if (info.ok) {
        se2 <- diag(solve(qr.info))
        if (min(se2) < 0) 
          se <- NA
        else {
          se <- sqrt(se2)
          se      <- sqrt(diag(solve(info)))
          se.beta <- matrix(se[1:(m*k)],m,k)
          se.alpha<- se[(m*k+1):(m*k+k)]*omega
          dimnames(se.beta)[2]<-list(y.names)
          dimnames(se.beta)[1]<-list(x.names)
          se  <- list(beta=se.beta, alpha=se.alpha, info=info)
        }
    }
  else 
    se <- NA       
  dp  <- list(beta=beta, Omega=Omega, alpha=alpha)
  opt$name <- algorithm
  list(call=match.call(), dp=dp, logL=logL, se=se, algorithm=opt)
}

 
msn.dev<-function(param, X, y, freq, trace=FALSE)
{
  d <- ncol(y)
  if(missing(freq)) freq<-rep(1,nrow(y))
  n <- sum(freq)
  m <- ncol(X)
  beta<-matrix(param[1:(m*d)],m,d)
  al.om<-param[(m*d+1):(m*d+d)]
  y0 <- y-X %*% beta
  Omega <- (t(y0) %*% (y0*freq))/n  
  D <- diag(qr(2*pi*Omega)[[1]])
  logDet <- sum(log(abs(D)))
  dev <- n*logDet-2*sum(zeta(0,y0 %*% al.om)*freq)+n*d
  if(trace) { 
    cat("\nmsn.dev:",dev,"\n","parameters:"); 
    print(rbind(beta,al.om))
    }
  dev
}

msn.dev.grad <- function(param, X, y, freq, trace=FALSE){
  d <- ncol(y)
  if(missing(freq)) freq<-rep(1,nrow(y))
  n <- sum(freq)
  m <- ncol(X)
  beta<-matrix(param[1:(m*d)],m,d)
  al.om<-param[(m*d+1):(m*d+d)]
  y0 <- y-X %*% beta
  Omega <- (t(y0) %*% (freq*y0))/n
  p1 <- zeta(1,as.vector(y0 %*% al.om))
  Dbeta <- t(X)%*% (y0*freq) %*%solvePD(Omega) - 
            outer(as.vector(t(X*freq)%*%p1), al.om)
  Dal.om <- as.vector(t(y0*freq) %*% p1)
  if(trace){
    cat("gradient:\n")
    print(rbind(Dbeta,Dal.om))}
  -2*c(Dbeta,Dal.om)
}

num.deriv1 <- function(x, FUN, ...)
{# calcola gradiente in modo numerico, se FUN calcola la funzione
   FUN <- get(FUN, inherit = TRUE)
   value <- FUN(x, ...)
   p <- length(x)
   grad <- numeric(p)
   delta <- cbind((abs(x) + 1e-10) * 1e-5, rep(1e-06, p))
   delta <- apply(delta, 1, max)
   for(i in 1:p) {
        x1 <- x
        x1[i] <- x1[i]+delta[i]
    grad[i] <- (FUN(x1, ...) - value)/delta[i]
   }
   grad
}
    
num.deriv2 <- function(x, FUN, ...)
{# derivate seconde numeriche, se FUN calcola il gradiente
   FUN <- get(FUN, inherit = TRUE)
   values <- FUN(x, ...)
   p <- length(values)
   H <- matrix(0, p, p)
   delta <- cbind((abs(x) + 1e-10) * 1e-5, rep(1e-06, p))
   delta <- apply(delta, 1, max)
   for(i in 1:p) {
    x1 <- x
        x1[i] <- x1[i]+delta[i]
    H[, i] <- (FUN(x1, ...) - values)/delta[i]
   }
   (H+t(H))/2
}
 
#---
# skew-t portion

dst <-  function (x, location=0, scale=1, shape=0, df=Inf, dp=NULL, log=FALSE)
{ 
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     df <- dp[4]
    }
  if(df==Inf) return(dsn(x,location,scale,shape, log=log))
  z   <- (x - location)/scale
  pdf <- dt(z, df=df, log=log)
  cdf <- pt(shape*z*sqrt((df+1)/(z^2+df)), df=df+1, log=log)
  if(log)
    log(2) + pdf + cdf -logb(scale)
  else
    2 * pdf * cdf / scale
}


rst <- function (n=1, location = 0, scale = 1, shape = 0, df=Inf, dp=NULL)
{ 
    if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     df <- dp[4]
    }
  z <- rsn(n,location=0,scale,shape)
  if(df==Inf) return(z+location)
  v <- rchisq(n,df)/df
  y <- z/sqrt(v)+location
  attr(y,"parameters")<- c(location,scale,shape,df)
  return(y)
}



pst <- function (x, location=0, scale=1, shape=0, df=Inf, dp=NULL, ...) 
{        
    if(!is.null(dp)) {
      if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
     }
    fp <- function(v, shape, df, t.value) psn(sqrt(v) * t.value, 0, 1, 
                    shape) * dchisq(v * df, df = df) * df
    if (df == Inf) 
        p <- psn(x, location, scale, shape)
    else 
      {
        if(df <= 0) stop("df must be non-negative")
        z <- (x-location)/scale
        p <- numeric(length(z))
        for (i in 1:length(z)){
          if(round(df)==df) 
             p[i] <- pmst(z[i], 0, matrix(1,1,1), shape, df, ...)
          else{
            if(abs(z[i]) == Inf)
               p[i] <- (1+sign(z[i]))/2
            else{
              if(z[i] < 10+50/df)  
                p[i] <- integrate(dst, -Inf, z[i], shape = shape, df = df,
                             ...)$value
              else           
                p[i] <- integrate(fp, 0, Inf, shape = shape, df = df,
                                t.value = z[i], ...)$value
        }}
      }
    pmax(0,pmin(1,p))
  }
}

qst <- function (p, location = 0, scale = 1, shape = 0, df=Inf,  
                 tol = 1e-08, dp = NULL, ...)
{
    if(!is.null(dp)) {
      if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
      }
    if (df == Inf) 
        return(qsn(p, location, scale, shape))
    max.q <- sqrt(qf(p, 1, df))
    min.q <- -sqrt(qf(1 - p, 1, df))
    if (shape == Inf) 
        return(location + scale * max.q)
    if (shape == -Inf) 
        return(location + scale * min.q)
    na <- is.na(p) | (p < 0) | (p > 1)
    zero <- (p == 0)
    one <- (p == 1)
    p <- replace(p, (na | zero | one), 0.5)
    cum <- st.cumulants(0, 1, shape, max(df,5), n=4)
    g1 <- cum[3]/cum[2]^(3/2)
    g2 <- cum[4]/cum[2]^2
    x <- qnorm(p)
    x <- (x + (x^2 - 1) * g1/6 + x * (x^2 - 3) * g2/24 -
           x * (2 *  x^2 - 5) * g1^2/36)
    x <- cum[1] + sqrt(cum[2]) * x
    max.err <- 1
    while (max.err > tol) {
        x1 <- x - (pst(x, 0, 1, shape, df, ...) - p)/dst(x, 0, 1, shape, df)
        x1 <- pmin(x1, max.q)
        x1 <- pmax(x1, min.q)
        max.err <- max(abs(x1 - x)/(1 + abs(x)))
        x <- x1
    }
    x <- replace(x, na, NA)
    x <- replace(x, zero, -Inf)
    x <- replace(x, one, Inf)
    return(location + scale * x)
}

st.cumulants <- function(location=0, scale=1, shape=0, df=Inf, dp=NULL, n=4)
{
  if(!is.null(dp)) {
      if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
      }
  if(df == Inf) return(sn.cumulants(location, scale, shape, n=n))
  n <- min(as.integer(n),4)
  if(df <= n) stop("need df>n")
  par <- cbind(location,scale,shape)
  delta <- par[,3]/sqrt(1+par[,3]^2)
  mu <- delta*sqrt(df/pi)*exp(lgamma((df-1)/2)-lgamma(df/2))
  cum<- matrix(NA, nrow=nrow(par), ncol=n)
  cum[,1]<- mu
  if(n>1) cum[,2] <- df/(df-2) - mu^2
  if(n>2) cum[,3] <- mu*(df*(3-delta^2)/(df-3) - 3*df/(df-2)+2*mu^2)
  if(n>3) cum[,4] <- (3*df^2/((df-2)*(df-4)) - 4*mu^2*df*(3-delta^2)/(df-3)
                      + 6*mu^2*df/(df-2)-3*mu^4)- 3*cum[,2]^2
  cum <- cum*outer(par[,2],1:n,"^")
  cum[,1] <- cum[,1]+par[,1]
  cum[,,drop=TRUE]
}


dmst <- function(x, xi = rep(0, length(alpha)), Omega, alpha,
                 df = Inf, dp=NULL, log = FALSE) 
{
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
    }  
    if (df == Inf) 
        return(dmsn(x, xi=xi, Omega=Omega, alpha=alpha, log = log))
    d <- length(alpha)
    Omega <- matrix(Omega,d,d)
    x <- if(is.vector(x)) matrix(x, 1, d) else data.matrix(x)
    if(is.vector(xi)) xi <- outer(rep(1,nrow(x)), xi)
    X <- t(x - xi)
    Q <- apply((solvePD(Omega) %*% X) * X, 2, sum)
    L <- as.vector(t(X/ sqrt(diag(Omega))) %*% as.matrix(alpha))
    logDet <- sum(logb(abs(diag(qr(Omega)$qr))))
    if(df < 10000)  {
          const<- lgamma((df + d)/2)- lgamma(df/2)-0.5*d*logb(df)
          log1Q <- logb(1+Q/df) 
          }
    else {
         const <- (-0.5*d*logb(2)+ log1p((d/2)*(d/2-1)/df))
         log1Q <- log1p(Q/df)
         }
    log.dmt <- const - 0.5*(d * logb(pi) + logDet + (df + d)* log1Q)                
    log.pt <- pt(L * sqrt((df + d)/(Q + df)), df = df + d, log = TRUE)
    logPDF <-  logb(2) + log.dmt + log.pt
    if (log) logPDF
    else exp(logPDF)
}

rmst <- function(n=1, xi=rep(0,length(alpha)), Omega, alpha, df=Inf, dp=NULL)
{ 
  if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
  if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
     }  
  d <- length(alpha)
  x <- if(df==Inf) 1 else rchisq(n,df)/df
  z <- rmsn(n, rep(0,d), Omega, alpha)
  y <- t(xi+ t(z/sqrt(x)))
  attr(y,"parameters") <- list(xi=xi, Omega=Omega, alpha=alpha, df=df)
  return(y)
}

pmst <- function(x, xi=rep(0,length(alpha)), Omega, alpha, df=Inf,
         dp= NULL, ...)
{
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
      }  
  d <- length(alpha)
  Omega<- matrix(Omega,d,d) 
  omega<- sqrt(diag(Omega))
  Ocor <- cov2cor(Omega)
  O.alpha <- as.vector(Ocor %*% alpha)
  delta <- O.alpha/sqrt(1+sum(alpha*O.alpha))
  Obig <- matrix(rbind(c(1,-delta),cbind(-delta,Ocor)),d+1,d+1)
  x <- c(0,(x-xi)/omega)
  if(df > .Machine$integer.max)  
    2 * pmnorm(x, mean=rep(0,d+1), varcov=Obig, ...)
  else 
    2 * pmt(x, mean=rep(0,d+1), S=Obig, df=df, ...)
}

  


dst2.plot <- function(x, y, xi, Omega, alpha, df, dp=NULL, ...)
{# plot bivariate density ST_2(xi,Omega,alpha,df) computed at (x,y) grid
    if(!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
       stop("You cannot set both component parameters and dp")
    if(!is.null(dp)){
      if(!is.null(dp$xi)) xi <- dp$xi
        else
      if(!is.null(dp$beta)) xi <- as.vector(dp$beta)
      Omega <- dp$Omega
      alpha <- dp$alpha
      df <- dp$df
      }  
    if(any(dim(Omega) != c(2, 2))) stop("dim(Omega) != c(2,2)")
    nx <- length(x)
    ny <- length(y)
    xoy <- cbind(rep(x, ny), as.vector(matrix(y, nx, ny, byrow = TRUE)))
    X <- matrix(xoy, nx * ny, 2, byrow = FALSE)
    pdf <- dmst(X, xi, Omega, alpha, df)
    pdf <- matrix(pdf, nx, ny)
    contour(x, y, pdf, ...)
    invisible(list(x = x, y = y, density = pdf, xi = xi, Omega = Omega,
        alpha = alpha, df=df))
}

mst.fit <- function(X, y, freq, start, fixed.df=NA, plot.it=TRUE, 
                 trace=FALSE, ...)
{
  y.name <- deparse(substitute(y))
  y.names<- dimnames(y)[[2]] 
  y <- as.matrix(y)
  d <- ncol(y)
  if(is.null(d)) d<- 1
  if(d>1){
    if(length(y.names)==0){
      dimnames(y) <-
         list(dimnames(y)[[1]], outer("V",as.character(1:d),paste,sep=""))
      y.names<- as.vector(dimnames(y)[[2]])
      }}
  else 
    colnames(y)<-y.name
  if(missing(freq)) freq <- rep(1,nrow(y))  
  n <- sum(freq)
  if(missing(X)) {
     X <- rep(1,nrow(y)) 
     missing.X <- TRUE }
  else 
     missing.X <- FALSE
  X   <- as.matrix(X)
  qrX <- qr(X)
  m   <- ncol(X)
  mle <- mst.mle(X=X, y=y, freq=freq,  fixed.df=fixed.df, start=start, 
                 trace=trace, ...)
  mle$call <- match.call()
  beta  <- mle$dp$beta
  Omega <- mle$dp$Omega
  alpha <- mle$dp$alpha
  omega <- sqrt(diag(Omega))
  df    <- mle$dp$df
  xi    <- X %*% beta
  if(plot.it & all(freq==rep(1,length(y)))) {
    if(missing.X) { 
      y0  <-y 
      xi0 <- apply(xi,2,mean)} 
    else  {
      y0  <- y-xi 
      xi0 <- rep(0,d)
      }
    if(d>1) {
       opt<-options()
       options(warn=-1)
       pairs(y0, labels=y.names,
        panel=function(x, y, Y, y.names, xi, Omega, alpha)
          {
          for(i in 1:length(alpha)){
            if(all(Y[,i]==x)) Ix<-i
            if(all(Y[,i]==y)) Iy<-i
            }
          points(x,y)
          marg <- msn.marginal(xi, Omega ,alpha, c(Ix,Iy))
          xi.marg <- marg$xi
          Omega.marg <- marg$Omega
          alpha.marg <- marg$alpha     
          x1 <- seq(min(x), max(x), length=30)
          x2 <- seq(min(y), max(y), length=30)
          dst2.plot(x1, x2, xi.marg, Omega.marg, alpha.marg, df, 
                    add=TRUE, col=2)
         },  # end "panel" function
      Y=y0, y.names=y.names, xi=xi0, Omega=Omega, alpha=alpha)
      options(opt)
      }
    else{ # plot for case d=1
      y0<-as.vector(y0)
      x<-seq(min(pretty(y0,10)),max(pretty(y0,10)),length=100)
      if(missing.X){
         dp0<-c(xi0,omega,alpha,df)
         xlab<-y.name}
      else {
         dp0<-c(0,omega,alpha,df)
         xlab <- "residuals"}
      hist(y0, prob=TRUE,  breaks="FD", xlab=xlab, ylab="density", main="")
      lines(x, dst(x,dp0[1],dp0[2],dp0[3],dp0[4]),  col=2)
      if(length(y)<101) points(y0, rep(0,n), pch=1)
      title(y.name)
      }
    cat("Press <Enter> to continue..."); readline()
    par(mfrow=c(1,2))
    pp  <- d*qf((1:n)/(n+1),d,df)
    pp2 <- qchisq((1:n)/(n+1),d)
    # Xb  <- qr.fitted(qrX,y)
    res <- qr.resid(qrX,y)
    rad.n  <- apply(res    * (res %*% solvePD(var(res))), 1, sum)
    rad.st <- apply((y-xi) * ((y-xi) %*% solvePD(Omega)), 1, sum)
    plot(pp2, sort(rad.n), pch=1, ylim=c(0,max(rad.n,rad.st)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for normal distribution", sub=y.name)
    plot(pp, sort(rad.st), pch=1, ylim=c(0,max(rad.n,rad.st)),
        xlab="Percentiles of chi-square distribution", 
        ylab="Mahalanobis distances")
    abline(0,1,lty=3)
    title(main="QQ-plot for skew-t distribution", sub=y.name)
    cat("Press <Enter> to continue..."); readline()
    plot((1:n)/(n+1), sort(pchisq(rad.n,d)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for normal distribution", sub=y.name)
    plot((1:n)/(n+1), sort(pf(rad.st/d,d,df)), xlab="", ylab="")
    abline(0,1,lty=3)
    title(main="PP-plot for skew-t distribution", sub=y.name)
    par(mfrow=c(1,1))
    } # end ploting
  dev.norm<- msn.dev(c(qr.coef(qrX,y),rep(0,d)), as.matrix(X), y, freq)
  test <- dev.norm + 2*mle$logL
  p.value <-  1-pchisq(test,d+1)
  if(trace) {
    cat("LRT for normality (test-function, p-value): ")
    print(c(test,p.value))
    }
  mle$test.normality <- list(LRT=test, df=d+1, p.value=p.value, 
                             normal.logL=dev.norm/(-2))
  invisible(mle)
}

#

st.mle <- function(X, y, freq,  start, fixed.df=NA, trace=FALSE, 
              algorithm = c("nlminb","Nelder-Mead", "BFGS", "CG", "SANN"),
              control=list())
{ 
  y.name  <- deparse(substitute(y))
  y <- data.matrix(y)
  if(missing(X)) X<- matrix(1, nrow=length(y), ncol=1)
  dimnames(y)[[2]] <- list(y.name) 
  if(missing(start)){
    cp0 <- sn.mle(X=X, y=y, plot.it=FALSE, trace=trace)$cp
    m <- length(cp0)-2
    cp0[m+2] <- cp0[m+2]*0.9
    mle0 <- cp.to.dp(cp0)
    start <- list(beta=mle0[1:m], Omega=matrix(mle0[m+1]^2,1,1),
                  alpha=mle0[m+2], df=10)
    }
  else {
    m <- length(start)-3
    if(m<1) stop("bad start vector")
    start<-  list(beta=start[1:m], Omega=matrix(start[m+1]^2,1,1),
                  alpha=start[m+2], df=start[m+3])
    } 
  fit <- mst.mle(X, y, freq, start=start, fixed.df=fixed.df, trace=trace, 
              algorithm=algorithm, control=control)
  mle <- list()
  mle$call<- match.call()
  dp <- fit$dp
  se <- fit$se
  p  <- length(dp$beta)
  dp.names <- c(if(p==1) "location" else dimnames(dp$beta)[[1]],
                "scale","shape","df")
  mle$dp  <- c(dp$beta, sqrt(as.vector(dp$Omega)), dp$alpha, dp$df)
  names(mle$dp) <- dp.names
  mle$se  <- c(se$beta, mle$dp[p+1] * se$internal[p+1], se$alpha,
               dp$df * se$internal[p+3])
  if(length(mle$se) == length(dp.names)) names(mle$se) <- dp.names
  mle$logL <- fit$logL
  mle
}


mst.mle <- function (X, y, freq, start, fixed.df = NA, trace = FALSE, 
                 algorithm=c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
                 control = list()) 
{
    algorithm <- match.arg(algorithm)
    y.name <- deparse(substitute(y))
    y.names <- dimnames(y)[[2]]
    y <- data.matrix(y)
    X <- if (missing(X)) matrix(rep(1, nrow(y)), ncol = 1)        
              else data.matrix(X)
    if (missing(freq)) freq <- rep(1, nrow(y))
    x.names <- dimnames(X)[[2]]
    d <- ncol(y)
    n <- sum(freq)
    m <- ncol(X)
    if (missing(start)) {
        qrX <- qr(X)
        beta <- as.matrix(qr.coef(qrX, y))
        Omega <- matrix(var(qr.resid(qrX, y)), d, d)
        omega <- sqrt(diag(Omega))
        alpha <- rep(0, d)
        df <- ifelse(is.na(fixed.df), 10, fixed.df)
        if (trace) {
            cat("mst.mle: dp=", "\n")
            print(c(beta, Omega, alpha))
            cat("df:", df, "\n")
        }
    }
    else {
        if (!is.na(fixed.df)) 
            start$df <- fixed.df
        if (all(names(start) == c("beta", "Omega", "alpha", "df"))) {
            beta <- start$beta
            Omega <- start$Omega
            alpha <- start$alpha
            df <- start$df
        }
        else stop("start parameter is not in the form that I expected")
    }
    eta <- alpha/sqrt(diag(Omega))
    Oinv <- solvePD(Omega)
    upper <- chol(Oinv)
    D <- diag(upper)
    A <- upper/D
    D <- D^2
    if (d > 1) 
       param <- c(beta, -log(D)/2, A[!lower.tri(A, diag = TRUE)], eta)
    else 
       param <- c(beta, -log(D)/2, eta)
    if (is.na(fixed.df)) 
        param <- c(param, log(df))
    if(algorithm == "nlminb"){
      opt <- nlminb(param, objective = mst.dev, gradient = mst.dev.grad, 
               control = control,  X = X, y = y, freq = freq, 
               trace = trace, fixed.df = fixed.df)
      info <- num.deriv2(opt$par, FUN="mst.dev.grad", X=X, y=y,
                 freq=freq, fixed.df = fixed.df)/2
      opt$value <-  opt$objective
      }
    else{
        opt <- optim(param, fn = mst.dev, gr = mst.dev.grad, 
                method = algorithm, control = control, hessian = TRUE,
                X = X, y = y, freq = freq, trace = trace, fixed.df = fixed.df)
        info <- opt$hessian/2
        }
    dev   <- opt$value
    param <- opt$par
    opt$name <- algorithm
    if (trace) {
        cat("Message from optimization routine:", opt$message, "\n")
        cat("deviance:", dev, "\n")
    }
    beta <- matrix(param[1:(m * d)], m, d)
    D <- exp(-2 * param[(m * d + 1):(m * d + d)])
    A <- diag(d)
    i0 <- m*d+d*(d+1)/2
    if(d>1)  A[!lower.tri(A,diag=TRUE)] <- param[(m*d+d+1):i0]
    eta <- param[(i0 + 1):(i0 + d)]
    if (is.na(fixed.df)) 
        df <- exp(param[i0 + d + 1])
    else df <- fixed.df
    Oinv <- t(A) %*% diag(D,d,d) %*% A
    Omega <- solvePD(Oinv)
    omega <- sqrt(diag(Omega))
    alpha <- eta * omega
    dimnames(beta) <- list(x.names, y.names)
    dimnames(Omega) <- list(y.names, y.names)
    if (length(y.names) > 0) names(alpha) <- y.names
    if (all(is.finite(info))) {
        qr.info <- qr(info)
        info.ok <- (qr.info$rank == length(param))
    }
    else info.ok <- FALSE
    if (info.ok) {
        se2 <- diag(solve(qr.info))
        if (min(se2) < 0) 
            se <- NA
        else {
            se <- sqrt(se2)
            se.beta <- matrix(se[1:(m * d)], m, d)
            se.alpha <- se[(i0 + 1):(i0 + d)] * omega
            dimnames(se.beta)[2] <- list(y.names)
            dimnames(se.beta)[1] <- list(x.names)
            names(se.alpha) <- y.names
            se.df <- df * se[i0 + d + 1]
            se <- list(beta = se.beta, alpha = se.alpha, df = se.df, 
                internal = se, info = info)
        }
    }
    else se <- NA
    dp <- list(beta = beta, Omega = Omega, alpha = alpha, df = df)
    list(call = match.call(), logL = -dev/2, deviance = dev, 
        dp = dp, se = se, algorithm = opt)
}



mst.dev <- function(param, X, y, freq=rep(1,nrow(X)), fixed.df=NA, trace=FALSE)
{
  # Diag <- function(x) diag(x, nrow=length(x), ncol=length(x))
  d <- ncol(y)
  # if(missing(freq)) freq<-rep(1,nrow(y))
  n <- sum(freq)
  m <- ncol(X)
  beta<-matrix(param[1:(m*d)],m,d)
  D <- exp(-2*param[(m*d+1):(m*d+d)]) 
  i0 <- m*d+d*(d+1)/2
  A <- diag(d)
  if(d>1) A[!lower.tri(A,diag=TRUE)] <- param[(m*d+d+1):i0]
  eta <- param[(i0+1):(i0+d)]
  if(is.na(fixed.df)) df <- exp(param[i0+d+1])
     else df <- fixed.df
  Oinv <- t(A) %*% diag(D,d,d) %*% A
  # Omega <- solvePD(Oinv)
  u <-  y - X %*% beta
  Q <- apply((u %*% Oinv)*u,1,sum)
  L <- as.vector(u %*% eta) 
  logDet<- sum(-log(D))
  if(df < 10000)  {
          const<- lgamma((df + d)/2)- lgamma(df/2)-0.5*d*logb(df)
          DQ <-  (df+d) * sum(freq *logb(1+Q/df))
          L. <- L*sqrt((df+d)/(Q+df))
          }
  else {
         const <- (-0.5*d*logb(2)+ log1p((d/2)*(d/2-1)/df))
         DQ <- if(df<Inf) (df+d) * sum(freq *log1p(Q/df)) else sum(freq*Q)
         L. <- L*sqrt((1+d/df)/(1+Q/df))
         }
  dev <- (n*(logDet - 2*const+ d*logb(pi)) + DQ
          -2*sum(freq * (log(2)+log.pt(L., df+d))))
  if(trace) cat("mst.dev: ",dev, "\n")
  dev
}


mst.dev.grad <- function(param, X, y, freq=rep(1,nrow(X)), fixed.df=NA,
                         trace=FALSE)
{
  # Diag <- function(x) diag(x, nrow=length(x), ncol=length(x))
  d <- ncol(y)
  n   <- sum(freq)
  m   <- ncol(X)
  beta<- matrix(param[1:(m*d)],m,d)
  D  <- exp(-2*param[(m*d+1):(m*d+d)])
  A  <- diag(d)
  i0 <- m*d+d*(d+1)/2
  if(d>1) A[!lower.tri(A,diag=TRUE)] <- param[(m*d+d+1):i0]
  eta   <- param[(i0+1):(i0+d)]
  if(is.na(fixed.df)) df <- exp(param[i0+d+1])
     else df <- fixed.df
  Oinv  <- t(A) %*% diag(D,d,d) %*% A
  u     <- y-X %*% beta
  Q     <- as.vector(apply((u %*% Oinv)*u,1,sum))
  L     <- as.vector(u %*% eta)
  sf    <- if(df<10000) sqrt((df+d)/(Q+df)) else sqrt((1+d/df)/(1+Q/df))
  t.    <- L*sf
  dlogft<- (-0.5)*(1+d/df)/(1+Q/df)
  dt.dL <- sf
  dt.dQ <- (-0.5)*L*sf/(Q+df)
  logT. <- log.pt(t., df+d)
  dlogT.<- exp(dt(t., df+d, log=TRUE)- logT.)
  u.freq<- u*freq
  Dbeta <- (-2* t(X) %*% (u.freq*dlogft) %*% Oinv 
            - outer(as.vector(t(X) %*% (dlogT. * dt.dL* freq)), eta)
            - 2* t(X) %*% (dlogT.* dt.dQ * u.freq) %*% Oinv )
  Deta  <- apply(dlogT.*sf*u.freq, 2, sum)
  if(d>1){
     M  <- 2*( diag(D,d,d) %*% A %*% t(u * dlogft
               + u * dlogT. * dt.dQ) %*% u.freq)
     DA <- M[!lower.tri(M,diag=TRUE)]
     }
  else DA<- NULL
  M     <- ( A %*% t(u*dlogft + u*dlogT.*dt.dQ) %*% u.freq %*% t(A))
  if(d>1) DD <- diag(M) + 0.5*n/D
     else DD <- as.vector(M + 0.5*n/D) 
  grad <- (-2)*c(Dbeta,DD*(-2*D),DA,Deta)
  if(is.na(fixed.df)) {
    df0<- if(df<Inf) df else 1e8
    if(df0<10000){
       diff.digamma<-  digamma((df0+d)/2) - digamma(df0/2)
       log1Q<- log(1+Q/df0)
     }
    else
      {
       diff.digamma<- log1p(d/df0)
       log1Q <- log1p(Q/df0)
      }
    dlogft.ddf <- 0.5 * (diff.digamma - d/df0
                        + (1+d/df0)*Q/((1+Q/df0)*df0) - log1Q)
    eps   <- 1.0e-4
    df1 <- df0 + eps
    sf1 <- if(df0<10000) sqrt((df1+d)/(Q+df1)) else sqrt((1+d/df1)/(1+Q/df1))
    logT.eps <- log.pt(L*sf1, df1+d)
    dlogT.ddf <- (logT.eps-logT.)/eps
    Ddf   <- sum((dlogft.ddf + dlogT.ddf)*freq)
    grad <- c(grad, -2*Ddf*df0)
    }
  if(trace) cat("mst.dev.grad: norm is ",sqrt(sum(grad^2)),"\n")  
 return(grad)
}
#-------------

st.2logL.profile<-function(X=matrix(rep(1,n)), y, freq, trace=FALSE,
          fixed.comp = c(ncol(X)+2, ncol(X)+3), 
          fixed.values = cbind(c(-4,4), log(c(1,25))),
          npts=30/length(fixed.comp), plot.it=TRUE, ...)
{# plot2D profile deviance (=2(max.logL-logL)) using either parameters
 # if(plot.it & !exists(.Device)) stop("Device not active")
 #
   if(missing(freq)) freq <- rep(1,length(y))
   n <- sum(freq)
   m <- ncol(X)
   if(length(fixed.comp) == 1){
      param1 <- seq(fixed.values[1], fixed.values[2], length=npts)
      logL <- param2 <- rep(NA,npts)}
   else{ 
      param1 <- seq(fixed.values[1,1], fixed.values[2,1], length=npts)
      param2 <- seq(fixed.values[1,2], fixed.values[2,2], length=npts)
      logL   <- matrix(NA,npts,npts)}
   ls <- lm.fit(X,y)
   omega <- sqrt(var(resid(ls)))
   param <- c(coef(ls), log(omega), 0, log(20))[-fixed.comp]
   max.logL <- (-Inf)
   if(trace) cat(c("Running up to",npts,":"))
   for(i in 1:npts){
     if(trace) cat(" ",i)
     if(length(fixed.comp) == 1) {
       opt  <- optim(param, fn=st.dev.fixed, method="Nelder-Mead",
                 X=X, y=y, freq=freq, trace=trace, 
                 fixed.comp=fixed.comp, fixed.values=param1[i])
       logL[i] <- opt$value/(-2)
       param <- opt$par
       if(logL[i] > max.logL) {
         max.logL<- logL[i]
         param <- numeric(m+3)
         param[fixed.comp]  <- param1[i]
         param[-fixed.comp] <- opt$par
         dp<- c(param[1:m], exp(param[m+1]), param[m+2], exp(param[m+3]))
         best <- list(fixed.comp1=param1[i], fixed.comp2=NA,
                      dp=dp, logL=max.logL, opt=opt)
       }}
     else{
     for(j in 1:npts){     
       opt  <- optim(param, fn=st.dev.fixed, method="Nelder-Mead",
                 X=X, y=y, freq=freq, trace=trace, 
                 fixed.comp=fixed.comp,
                 fixed.values=c(param1[i], param2[j] ))
       logL[i,j] <- opt$value/(-2)
       if(j==1) param0 <- opt$par
       if(j<npts) 
         param <- opt$par
       else
         param <- param0
       if(logL[i,j] > max.logL) {
         max.logL<- logL[i,j]
         param <- numeric(m+3)
         param[fixed.comp]  <- c(param1[i], param2[j])
         param[-fixed.comp] <- opt$par
         dp<- c(param[1:m], exp(param[m+1]), param[m+2], exp(param[m+3]))
         best <- list(fixed.comp1=param1[i], fixed.comp2=param2[j],
                      dp=dp, logL=max.logL, opt=opt)
       }   
     }}
   }
  if(trace) cat("\n") 
  dev <- 2 * (max(logL) - logL)
  if(plot.it){
    if(length(fixed.comp) == 1){ 
      plot(param1, dev, type="l", ...)
      points(x=best$fixed.comp1, y=0, pch=1)
      }
    else{
      contour(param1, param2, dev,   labcex=0.75, 
          levels=c(0.57, 1.37, 2.77, 4.6, 5.99, 9.2), 
          labels=c(0.25, 0.5,  0.75, 0.90,0.95, 0.99),
          ...)
    points(x=best$fixed.comp1, y=best$fixed.comp2, pch=1)  
    }
  }    
  title(main=paste("dataset:", deparse(substitute(y)),
        "\nProfile deviance", sep= " "))
  invisible(list(call=match.call(), param1=param1, param2=param2,
            deviance=dev, max.logL=max.logL, best=best))
}


st.dev.fixed <- function(free.param, X, y, freq, trace=FALSE, 
                fixed.comp=NA,  fixed.values=NA)
{# param components: beta, log(omega), alpha, log(df)
  n <- sum(freq)
  m <- ncol(X)
  param <- numeric(length(free.param)+length(fixed.comp))
  param[fixed.comp]  <- fixed.values
  param[-fixed.comp] <- free.param
  beta  <- param[1:m]
  omega <- exp(param[m+1])
  eta <- param[m+2]/omega
  df  <- exp(param[m+3])
  u <-  y - X %*% beta
  Q <- freq*(u/omega)^2
  L <- u*eta 
  logDet <- 2*log(omega)
    if(df < 10000)  {
          const<- lgamma((df + 1)/2)- lgamma(df/2)-0.5*logb(df)
          log1Q <- logb(1+Q/df) 
          }
  else {
         const <- (-0.5*logb(2)+ log1p((1/2)*(-1/2)/df))
         log1Q <- log1p(Q/df)
         }
  dev <- (n*(logDet - 2*const+ logb(pi)) + (df+1) * sum(freq * log1Q)        
         -2*sum(log(2)+log.pt(L * sqrt((df+1)/(Q+df)),df+1)))
  if(trace) cat("st.dev.fixed (param, dev): ", param, dev,"\n")
  dev
}

#----

sn.SFscore <- function(shape, z, trace=FALSE)
{# Sartori-Firth's modified score function (2nd version)
  a42<- integrate(function(x) dsn(x,0,1,shape) * x^4 * zeta(1,shape*x)^2,
                   -Inf, Inf)$value
  a22<- integrate(function(x) dsn(x,0,1,shape) * x^2 * zeta(1,shape*x)^2,
                   -Inf, Inf)$value
  score <- sum(zeta(1,shape*z)*z)-0.5*shape*a42/a22
  if(trace) cat("sn.SFscore: (shape,score)=", shape, score,"\n")
  score
}

sn.mmle <- function(X, y, plot.it=TRUE, trace=FALSE,...)
   {
     n <- length(y)
     if (missing(X)){
         X <- as.matrix(rep(1, n))
         colnames(X) <- "constant"
       }
     m  <- ncol(X)
     dp <- cp.to.dp(sn.mle(X=X, y=y, plot.it=plot.it, trace=trace,...)$cp)
     z  <- (y - as.vector(X %*% dp[1:m]))/dp[m+1]
     start <- sign(dp[m+2])*min(5000,abs(dp[m+2]))
     a0 <- start/4
     f0 <- sn.SFscore(a0, z, trace=trace)
     a1 <- start
     f1 <- sn.SFscore(a1, z, trace=trace)
     while(f0*f1 > 0){
       a1 <- a0
       f1 <- f0
       a0 <- a0/4
       f0 <- sn.SFscore(a0, z, trace=trace)
       if(trace) cat("interval: ", a0,a1,"\n")
       }
     if(trace)cat("a0, a1: ",a0, a1,"\n")
     a <- uniroot(sn.SFscore, interval=c(a0, a1), z=z, trace=trace)
     dp <- sn.em(X, y, fixed=c(NA,NA,a$root), trace=trace)$dp
     if (plot.it) {
       dp0 <- dp
       if (all(X == rep(1, n)))
         y0 <- y
       else {
         y0 <- y - as.vector(X %*% dp0[1:m])
         dp0 <- c(0, dp0[m + 1], dp0[m + 2])
         xlab <- "residuals"
         }
       curve(dsn(x, dp=dp0), add=TRUE, lty=2, col=3)
       }
     names(dp)[m+2] <- "shape"
     info <- sn.Einfo(dp=dp,x=X)
     list(call=match.call(), dp=dp, se=info$se.dp, Einfo=info$info.dp)
   }

st.SFscore <- function(shape, df, z, trace=FALSE)
{# Sartori-Firth's modified score function for skew-t case  
   U <- function(x,shape,df){
          u <- x*sqrt((df+1)/(x^2+df))
          u * dt(shape*u, df=df+1)/pt(shape*u, df=df+1)      
        }
   J <- function(x,shape,df){
          u <- x*sqrt((df+1)/(x^2+df))
          t <- dt(shape*u, df=df+1)
          T <- pt(shape*u, df=df+1)
          ((df+1)*shape*u^3*t/((df+2)*(1+(shape*u)^2/(df+1)))+(t*u/T)^2)
        }
   EJ <- integrate(function(x, shape=shape, df=df)
                 J(x,shape=shape, df=df) * dst(x,0,1,shape,df),
                 -Inf,Inf, shape=shape, df=df)$value
   nu111 <- integrate(function(x,shape=shape, df=df)
                 U(x,shape=shape, df=df)^3 * dst(x,0,1,shape,df),
                 -Inf,Inf, shape=shape, df=df)$value
   nu1.2 <- integrate(function(x, shape=shape, df=df)
                 U(x,shape=shape, df=df) * J(x,shape=shape, df=df) *
                 dst(x,0,1,shape,df), -Inf, Inf, shape=shape, df=df)$value
   M <- 0.5*(nu111+nu1.2)/EJ
   u <- z*sqrt((df+1)/(z^2+df))
   score <- sum(u * dt(shape*u, df=df+1)/pt(shape*u, df=df+1)) + M
   if(trace) cat("st.SFscore(shape,score):", shape, score,"\n")
   score
}

st.mmle <- function(X, y, df, trace=FALSE)
   {
     n <- length(y)
     if (missing(X)){
         X <- as.matrix(rep(1, n))
         colnames(X) <- "constant"
       }
     m <- ncol(X)
     dp <- st.mle(X=X, y=y, fixed.df=df, trace=trace)$dp
     z <- (y - as.vector(X %*% dp[1:m]))/dp[m+1]
     start <- sign(dp[m+2])*min(5000,abs(dp[m+2]))
     a0 <- start/4
     f0 <- st.SFscore(a0, df, z, trace=trace)
     a1 <- start
     f1 <- st.SFscore(a1, df, z, trace=trace)
     while(f0*f1 > 0){
       a1 <- a0
       f1 <- f0
       a0 <- a0/4
       f0 <- st.SFscore(a0, df, z, trace=trace)
       }
     if(trace) cat("st.mmle: (a0, a1)= ",a0, a1,"\n")
     a <- uniroot(st.SFscore, interval=c(a0, a1), df=df, z=z, trace=trace)
     dp <- c(dp[1:(m+1)], shape=a$root, df=df)
     list(call=match.call(), dp=dp)
   }


sn.Einfo <- function(dp=NULL, cp=NULL, n=1, x=NULL)
{# computes Expected  Fisher information matrix for SN variates
  if(is.null(dp) & is.null(cp)) stop("either dp or cp must be set")
  if(!is.null(dp) & !is.null(cp)) stop("either dp or cp must be set")
  if(is.null(cp)) cp<- dp.to.cp(dp)
  if(is.null(dp)) dp<- cp.to.dp(cp)
  if(is.null(x))
     {
      x <-matrix(rep(1,n),nrow=n,ncol=1)
      xx <- n
      sum.x <- n
      p <- 1
     }
    else
    { if(n!=1) warning("You have set both n and x, I am setting n<-nrow(x)")
      n <- nrow(x)
      p <- ncol(x)
      xx <- t(x) %*% x
      sum.x <- matrix(apply(x,2,sum))
    }
  if(length(cp) != (p+2)| length(dp) != (p+2))
        stop("length(dp) | length(cp) must be equal to ncol(x)+2")
  omega <- dp[p+1]
  alpha <- dp[p+2]
  E.z   <- sqrt(2/pi)*alpha/sqrt(1+alpha^2)
  s.z   <- sqrt(1-E.z^2)
  I.dp  <- matrix(NA,p+2,p+2)
  if(abs(alpha)< 200){
    a0 <- integrate(function(x) dsn(x,0,1,alpha) * zeta(1,alpha*x)^2,
          -Inf, Inf)$value
    a1 <- integrate(function(x) dsn(x,0,1,alpha) *x * zeta(1,alpha*x)^2,
          -Inf, Inf)$value
    a2 <- integrate(function(x) dsn(x,0,1,alpha) *x^2 * zeta(1,alpha*x)^2,
          -Inf, Inf)$value
    }
  else
    {a0 <- sign(alpha)*0.7206/abs(alpha)
     a1 <- -sign(alpha)*(0.6797/alpha)^2
     a2 <- 0.005897/alpha^2 + 30.611/alpha^4
    }
  I.dp[1:p,1:p] <- xx* (1+alpha^2*a0)/omega^2  
  I.dp[p+1,p+1] <- n * (2+alpha^2*a2)/omega^2
  I.dp[p+2,p+2] <- n * a2
  I.dp[1:p,p+1] <- sum.x * (E.z*(1+E.z^2*pi/2)+alpha^2*a1)/omega^2
  I.dp[p+1,1:p] <- t(I.dp[1:p,p+1])
  I.dp[1:p,p+2] <- sum.x * (sqrt(2/pi)/(1+alpha^2)^1.5-alpha*a1)/omega
  I.dp[p+2,1:p] <- t(I.dp[1:p,p+2])
  I.dp[p+1,p+2] <- I.dp[p+2,p+1] <- n*(-alpha*a2)/omega
  # cp <- dp.to.cp(dp)
  sigma <-cp[p+1]
  gamma1<-cp[p+2]
  D  <- diag(p+2)
  R  <- E.z/s.z
  T  <- sqrt(2/pi-(1-2/pi)*R^2)
  Da.Dg <- 2*(T/(T*R)^2+(1-2/pi)/T^3)/(3*(4-pi))
  DE.z <- sqrt(2/pi)/(1+alpha^2)^1.5
  Ds.z <- (-E.z/s.z)*DE.z
  D[1,p+1] <- (-R)
  D[1,p+2] <- (-sigma*R)/(3*gamma1)
  D[p+1,p+1] <- 1/s.z
  D[p+1,p+2] <- (-sigma)* Ds.z* Da.Dg/s.z^2
  D[p+2,p+2] <- Da.Dg
  I.cp  <- t(D) %*% I.dp %*% D
  I.cp <- (I.cp + t(I.cp))/2
  se.dp <- sqrt(diag(solvePD(I.dp)))
  se.cp <- sqrt(diag(solvePD(I.cp)))
  dimnames(I.cp)<- list(names(cp), names(cp))
  dimnames(I.dp)<- list(names(dp), names(dp))
  aux <- list(Deriv=D, a.int=c(a0,a1,a2))
  list(dp=dp, cp=cp, info.dp=I.dp, info.cp=I.cp, se.dp=se.dp, se.cp=se.cp, aux=aux)
}

#----


sn.logL.grouped <- function(param, breaks, freq, trace=FALSE)
{
  cdf <- pmax(psn(breaks, param[1],exp(param[2]), param[3]), 0)
  p <- diff(cdf)
  logL <- sum(freq*log(p))
  if(trace) print(c(param, logL))
  logL
}

sn.mle.grouped <- function(breaks, freq, trace=FALSE, start=NA)
{
  if(any(is.na(start))){
    b <- breaks
    d <- diff(b)
    if(b[1]== -Inf) b[1]<- b[2]-d[2]
    if(b[length(b)]==Inf) b[length(b)] <- b[length(b)-1]+d[length(d)-1]
    mid<- (b[-1]+b[-length(b)])/2
    dp <- msn.mle(y=mid, freq=freq, trace=trace)$dp
    start <- c(dp[[1]], log(sqrt(dp[[2]])), dp[[3]])
    }
  opt <- optim(start,  sn.logL.grouped, 
              control=list(fnscale=-1),
              breaks=breaks, freq=freq, trace=trace)
  param <- opt$par
  dp <- c(param[1], exp(param[2]), param[3]) 
  invisible(list(call=match.call(), dp=dp, logL=opt$value, end=param, opt=opt))
}


st.logL.grouped <- function(param, breaks, freq, trace=FALSE)
{
  if(param[4] > 5.5214609) # 5.5214609=log(250)
      cdf<- psn(breaks, param[1], exp(param[2]), param[3])
    else
      cdf<- pst(breaks, param[1], exp(param[2]), param[3], exp(param[4]))
  p <- pmax(diff(cdf), 1.0e-10)
  logL <- sum(freq*log(p)) 
  if(trace) print(c(param, logL))
  logL
}

st.mle.grouped <- function(breaks, freq, trace=FALSE, start=NA)
{
  if(any(is.na(start))){
    a <- sn.mle.grouped(breaks, freq)
    start <- c(a$end, log(15))
    if(trace)  cat("Initial parameters set to:", format(start),"\n")
    }
  opt <- optim(start,  st.logL.grouped, 
            control=list(fnscale=-1),
            breaks=breaks, freq=freq, trace=trace)  
  param<-opt$par
  dp <- c(param[1],exp(param[2]),param[3], exp(param[4])) 
  logL <- opt$value 
  invisible(list(call=match.call(), dp=dp, logL=logL, end=param, opt=opt))
}

msn.affine <- function(dp, a=0, A, drop=TRUE)
{
# computes distribution of affine transformation of MSN/MST variate, T=a+AY,
# using formulae in Appendix A.2 of Capitanio et al.(2003)
#
  Diag  <- function(x) diag(x,nrow=length(x),ncol=length(x))
  if(is.null(dp$xi))  xi <- dp$beta  else  xi <- dp$xi
  xi.T  <- as.vector(A %*% matrix(xi,ncol=1)+a)
  Omega <- dp$Omega
  O.T   <- as.matrix(A %*% Omega %*% t(A)) 
  oi    <- Diag(1/sqrt(diag(Omega)))
  B     <- oi %*% Omega %*% t(A)
  tmp   <- (oi %*% Omega %*% oi - B %*% solvePD(O.T) %*% t(B)) %*% dp$alpha
  den   <- sqrt(1+sum(dp$alpha*as.vector(tmp)))
  num   <- Diag(sqrt(diag(O.T))) %*% solvePD(O.T) %*% t(B) %*% dp$alpha
  alpha <- as.vector(num/den)
  if(all(dim(O.T)==c(1,1)) & drop)
     dp.T<- list(location=xi.T, scale=sqrt(as.vector(O.T)), shape=alpha)
  else
     dp.T <- list(xi=xi.T, Omega=O.T, alpha=alpha)
  if(!is.null(tau=dp$tau)) dp.T$tau <- dp$tau
  if(!is.null(tau=dp$df)) dp.T$df <- dp$df
  return(dp.T)
}

mst.affine <- function(dp, a=0, A, drop=TRUE) msn.affine(dp, a, A, drop)
                                     
#---

st.cumulants.inversion <- function(cum, abstol=1e-8)
{
  st.cumulants.matching <- function(par, gamma)
    {
      cum <- st.cumulants(shape=par[1], df=exp(par[2])+4, n=4)
      g1  <- cum[3]/cum[2]^1.5
      g2  <- cum[4]/cum[2]^2
      (abs(g1-gamma[1])^1.5/(1+abs(gamma[1])) +
         abs(g2-gamma[2])^1.5/(1+gamma[2]))
    }
  if(length(cum) != 4) stop("cum must be a vector of length 4")
  g1 <- cum[3]/cum[2]^1.5
  g2 <- cum[4]/cum[2]^2
  # if(g2<0) {
  #     warning("cumulants matching may be inaccurate")
  #     return(c(location=cum[1], scale=sqrt(cum[2]), shape=0, df=Inf))
  #   } 
  opt1 <- optim(c(0,1), st.cumulants.matching,
            control=list(abstol=abstol), gamma=c(g1,g2))
  opt2 <- nlminb(c(0,1), st.cumulants.matching,
            control=list(abs.tol=abstol), gamma=c(g1,g2))
  if(opt1$value < opt2$objective) par<- opt1$par else par<- opt2$par
  if(min(opt1$value, opt2$objective) > abstol) 
    warning("cumulants matching may be inaccurate")
  alpha <- par[1]
  df  <- exp(par[2])+4
  cumZ <- st.cumulants(dp=c(0,1,alpha,df))
  omega <- sqrt(cum[2]/cumZ[2])
  c(location=cum[1]-omega*cumZ[1], scale=omega, shape=alpha, df=df)
}


sample.centralmoments <- function(x, w=rep(1,length(x)), order=4)
{ # central moments, but first term is ordinary mean
  if( order < 1 | order != round(order))
    stop("order must be a positive integer") 
  x <- as.vector(x)
  m <- weighted.mean(x, w=w, na.rm = TRUE)
  mom <- rep(0,4)
  mom[1] <- m
  if(order > 1) 
    for(k in 2:order)
      mom[k] <- weighted.mean((x-m)^k, w=w, na.rm = TRUE)
  mom
}

solvePD <- function(x) 
{ # inverse of a symmetric positive definite matrix
    u <- chol(x, pivot = FALSE)
    if(prod(diag(u)) <= 0) stop("matrix not positive definite")
    # ui <- backsolve(u,diag(ncol(x)))
    # ui %*% t(ui)
    chol2inv(u)
}



#---
log.pt <- function(x, df){
   # fix for log(pt(...)) when it gives -Inf 
   # see Abramowitz & Stegun formulae 26.7.8 & 26.2.13)
   # However, new releases of R (>=2.3) seem to have fixed the problem
   if(df == Inf) return(pnorm(x, log.p=TRUE))
   p <- pt(x, df=df, log.p=TRUE)
   ninf <- (p == -Inf) 
   x0 <- (1-1/(4*df))*(-x[ninf])/sqrt(1+x[ninf]^2/(2*df))
   p[ninf] <- dnorm(x0,log=TRUE)-log(x0)+log1p(-1/(x0^2+2))
   p
}

#---

.onAttach <- function(library, pkg)
{
  Rv <- R.Version()
  if(Rv$major < 2 |(Rv$major == 2 && Rv$minor < 2.0))
    stop("This package requires R 2.2.0 or later")
  assign(".sn.home", file.path(library, pkg),
         pos=match("package:sn", search()))
  sn.version <- "0.4-2 (2006-10-26)"
  assign(".sn.version", sn.version, pos=match("package:sn", search()))
  if(interactive())
  {
    cat(paste("Package 'sn', ", sn.version, ". ",sep=""))
    cat("Type 'help(SN)' for summary information\n")
  }
  invisible()

}

 */