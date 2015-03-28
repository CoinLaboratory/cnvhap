package lc1.stats;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.states.EmissionState;
import lc1.util.Constants;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import JSci.maths.statistics.ProbabilityDistribution;
import cern.colt.matrix.DoubleMatrix2D;

public class TrainableNormal extends ProbabilityDistribution implements lc1.stats.ProbabilityDistribution,  Comparable{
	 static lc1.stats.NormalDistribution normal;
	 public double location, scale;//, shape;//origLoc, origScale, origShape;
    public TrainableNormal clone(){
        return new TrainableNormal(this);
    }
    public void setCoverage(double d){
    	throw new RuntimeException("!!");
    }
    public TrainableNormal clone(double u){
        return new TrainableNormal(this, u);
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
    	disty.addCount(b, val);
    	mixe1.addCount(0, val);
		//throw new RuntimeException("!!");
	}
@Override
public void plotTheoretical(String name1, boolean cum, XYSeries newD) {
   // XYSeries newD = new XYSeries(name1);
    double min = normal.quantile(1e-7, this.location, this.scale);
    double max = normal.quantile(1.0-1e-7, this.location, this.scale);
    double incr = 0.001;
    for(double j =min; j<max; j+=incr){
        double res = cum ? normal.cdf(j, this.location, this.scale) : normal.pdf(j, this.location, this.scale);
        if(res > 1e-11){
        newD.add(j, res);
        }
        else{
            System.err.println("res");
        }
    }
   
  //  return newD;
}
public String toString(){
	  return name;
}
public void transfer(double ps){
	if(this.obsx.size()==0) return;
    this.maximise(ps, ps, ps);
    
}

public void setParamValue(int n1, double val) {
	if(true ) throw new RuntimeException("!!");
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
public  double getParamValue(int i){
    if(i==0) return this.location;
    if(i==1) return this.scale;
   // if(i==2) return this.shape;
    throw new RuntimeException("!!");
}
public int getParamIndex(){
    return this.paramIndex;
}
public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution[] tmp) {
   if(true) throw new RuntimeException("!!");
    
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
public int numObs() {
	return this.obsx.size();
}
public String id() {
	return id;
}
public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
    TrainableNormal n = (TrainableNormal) probabilityDistribution;
    for(int i=0; i<n.obsx.size(); i++){
    	this.addCount(n.obsx.get(i), n.obsv.get(i));
    }
   /* for(Iterator<Map.Entry<Double,Double>> it = n.observations.entrySet().iterator(); it.hasNext();){
        Map.Entry<Double,Double> nxt = it.next();
        this.addCount(nxt.getKey(),  nxt.getValue());
    }*/
}
int paramIndex=0;
public void updateParamIndex(){
    this.paramIndex++;
}
public double prior(){
     throw new RuntimeException("!!");
} public double sample() {
    throw new RuntimeException("!!");
}
public double scale(){
	return this.scale;
}
public double getMean(){
    return this.location;
}
public void initialise() {
    obsx.clear();obsv.clear();
    sum =0;
}
public void setParam(int type,  double d){
	//if(true) throw new RuntimeException("!!");
	if(type==0){
		this.meanPrior = d;
		this.location = d;
	}
	else if(type==1){
	/*	if(this.id.indexOf('A')>0 && id.indexOf('B')>0){
			System.err.println("h");
		}*/
		double d1 =Math.sqrt(Math.max(1e-10,d));
		/*if(Constants.CHECK && (Double.isNaN(d1) || d1==0)) {
			throw new RuntimeException("!!");
		}*/
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
public int compareTo(Object o) {
	return this.name.compareTo(((TrainableNormal)o).name);
}

public void setPriors(lc1.stats.ProbabilityDistribution distx, int type){
	
	TrainableNormal dist1 = (TrainableNormal) (distx instanceof Mixture ? ((Mixture)distx).dist[0] : distx);
	if(type==1)this.stddevPrior = dist1.scale;
	//else if (type==2) throw new RuntimeException("!!");//this.skewPrior[0] = dist1.shape;
	else if(type==0) this.meanPrior = dist1.location;
	//else throw new RuntimeException("!!");
	//this.dist[0].setPriors(((Mixture)distx).dist[0]);
}
public void print(PrintWriter pw){
	this.recalcName();
	pw.print(this.toString()+"\t");
}
public void recalcName() {
	name =   id+"_"+String.format("%5.3g,", location)+String.format("%5.3g,", scale);//+String.format("%5.3g", shape);
   
}
@Override
public void  getInterval(double[] input, double[] res) {
	for(int i=0; i<res.length; i++){
	res[i] =  normal.quantile(input[i], this.location, this.scale);
	}
	
  //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
}
public List<Double> obsx = new ArrayList<Double>();
public List<Double>obsv = new ArrayList<Double>();
public int size(){
    return obsx.size();
}
double round;
final public String id; //
public  double meanPrior, stddevPrior;// skewPrior;
private double priorModifier=1.0;
public TrainableNormal(String name, double mean, double stddev, double round, double priorMod){
	this.id = name;
	this.location = mean;
	this.scale = stddev;
	if(Constants.CHECK && scale==0){
		throw new RuntimeException("!!");
	}
	this.round = round;
this.meanPrior = mean;
this.stddevPrior =stddev;
//    super(name, mean, stddev, 0.0, 
 //           new double[] {normal.quantile(1e-3, mean, stddev), 0.001, 0},  new double[] { normal.quantile(0.999, mean, stddev), 1e3, 0},
  //          round, priorMod
   //        );
   
}
public TrainableNormal(String name, double mean, double stddev, double stddevprior, double round, double priorMod){
	this.id = name;
	this.name = name;
	this.location = mean;
	this.scale = stddev;
	this.round = round;
this.meanPrior = mean;
this.stddevPrior =stddevprior;
this.priorModifier = priorMod;
if(Constants.CHECK && scale==0){
	throw new RuntimeException("!!");
}
  /*  super(name, mean, stddev, stddevprior, 0.0, 
            new double[] {normal.quantile(1e-3, mean, stddev), 0.001, 0},  new double[] { normal.quantile(0.999, mean, stddev), 1e3, 0},
            round, priorMod
           );*/
   
}
private String name;
public String name(){
    return id;
}
public TrainableNormal(String string, int i) {
   this(null, 0,0,0, 1.0);
    throw new RuntimeException("!!");
}

public TrainableNormal(TrainableNormal skewNormal) {
	this.id = skewNormal.id;
	  this.priorModifier = skewNormal.priorModifier;
      this.round = skewNormal.round;
   this.name = skewNormal.name;
   this.meanPrior = skewNormal.meanPrior;
   this.location = skewNormal.location;
   this.scale = skewNormal.scale;
   this.stddevPrior = skewNormal.stddevPrior;
    // TODO Auto-generated constructor stub
}

public TrainableNormal(TrainableNormal skewNormal, double u) {
    this(skewNormal);
   // this.setPrior(u, u, u);
    this.location =normal.quantile(Constants.rand.nextDouble(),location, 1.0/u);
    if(Constants.CHECK && scale==0){
		throw new RuntimeException("!!");
	}
  }
  
public double dsn (double x, boolean log)
{
 
 return log ? Math.log(this.probability(x)) : this.probability(x);
}




    @Override
    public double cumulative(double arg0) {
        return normal.cdf(arg0, this.location, this.scale);
    }

    @Override
    public double inverse(double arg0) {
        return normal.quantile(arg0,location, scale);
//        throw new RuntimeException("!!");
//       return normal.inverse(arg0);
       
    }

    @Override
    public double probability(double arg0) {
    	
//    	double res =  ( normal.cdf(arg0+0.001, location, scale) - normal.cdf(arg0-0.001, location, scale))/0.002;
    		
    	//	(this.dist.cumulativeProbability(arg0+0.001) - 
        	//	this.dist.cumulativeProbability(arg0-0.001))/0.002;
   	double res = normal.pdf(arg0, location, scale);
    	
    return res;
   //  return sc;
  /*      if(prob==0){
        	throw new RuntimeException("!!");
        }
       return prob;*/
    }

  
    
   //double tot=0;
   //double exp=0;
    public double average(double pseudo){
        double tot =pseudo;
        double exp =pseudo * this.meanPrior;
        for(int i=0; i<obsx.size(); i++){
           
            double sc = obsx.get(i);
            double w = obsv.get(i);
           
            exp+=sc*w;
            tot+=w;
        }
       // if(Math.abs(mean- Constants.r_mean()[3])<0.01){
        //   System.err.println(this.meanPrior.getMean()+" rvalues "+observations);
        //   System.err.println(this.meanPrior.getMean()+" rvalues "+weights);
       // }
        return exp/tot;
    }
  
    public double variance(double mu, double pseudo){
    	
        double tot =pseudo;
      double exp =Math.pow(this.stddevPrior,2)*pseudo;
      for(int i=0; i<obsx.size(); i++){
      
          double sc =Math.pow(obsx.get(i)-mu, 2);
          double w =obsv.get(i);
            exp+=sc*w;
            tot+=w;
          //  if(exp < 0 || tot < 0) throw new RuntimeException("!! "+w+" "+sc);
        }
        if(Constants.CHECK && Double.isInfinite(exp/tot) || tot==0){
            throw new RuntimeException(" ");//+this.observations+"\n"+this.observations);
        }
       
        return exp/tot;
    }
    public String getObsString(){
    	return this.name();
   /*      if(observations.size()>0) return median()+"";
        else return "-";*/
    }
    public double median(){
    throw new RuntimeException("!!");
/*        double sum=0;
        double tot = sum();
      //  StringBuffer sb = new StringBuffer(tot+"_sum ");
       // Histogram hist = new Histogram(observations.firstKey(), observations.lastKey(), 10);
        for(Iterator<Double> it = this.observations.keySet().iterator(); it.hasNext();){
           Double d =  it.next();
           sum+=observations.get(d)/tot;
           if(sum>0.5) return d;
        //   sb.append("("+d+":"+sum+"),");
        }
        return observations.lastKey();//        return sb.toString();*/
    }
    public double sum(){
        return sum;
      }
    
//static final double cntThresh = Constants.countThresh();
    public synchronized void addCount(double d, double w) {
      //  if(!Double.isNaN(d) &&  !Double.isInfinite(d) 
          //      && d!=Double.NEGATIVE_INFINITY
      //          && d!=Double.POSITIVE_INFINITY && w > Constants.countThresh()){
       //   if( d < 0) System.err.println("less than zero "+d+" "+w+this.name);
       // double cntThresh = 
    //    if(w > Constants.countThresh2()){
    	
        	 if(Constants.CHECK && (Double.isNaN(d) || Double.isInfinite(d))) throw new RuntimeException("!! "+d+" "+w);
        this.addObservation(d, w);
      //  }
   // if(Math.abs(location- -5.0)<0.01){
  //      System.err.println("hj");
   //   }
        }
    double sum=0;
    public double round(double d){
        return Math.round(d*round)/round;
    }
        public void addObservation(double d1, double weight){
            //double d1 = round(d);
        	this.obsx.add(d1);
        	this.obsv.add(weight);
//            Double w = observations.get(d1);
          //  observations.put(d1, w==null ? weight : weight+w);
            sum+=weight;
        }
        @Override
    public void maximise(double pseudo, double pseudoSD, double pseudoSkew) {
     //   if(this.id.startsWith("[A, [A, B]]") ){
        	
        //		Logger.global.info("h");
        //	}
       // 	if(true) throw new RuntimeException("!!");
      this.maximise(pseudo, pseudoSD, pseudoSkew, Constants.trainThresh());
   // scale = Math.max(scale, 0.0001);
    }
 
        public double[] getCount(double[] angle){
            double[] res = new double[angle.length];
            for(int i=0; i<res.length; i++){
              
                res[i] =  this.cumulative(angle[i]) - 
                  ( i==0 ? 0 : this.cumulative(angle[i-1]));
            }
            return res;
        }
        
    public void maximise(double pseudo,  double pseudoSD, double pseudoSkew, double trainThresh) {
    	double sum = sum();
    	if(sum<trainThresh){
    		this.location = this.meanPrior;
    		this.scale = this.stddevPrior;
    		return;
    	}
          double mean =pseudo <10e4 ? average(pseudo):this.meanPrior;
        //  double mean1 = average(0);
          if(pseudoSD<10e4){
          double variance = variance(mean, pseudoSD);
          if(Constants.CHECK && (Double.isNaN(variance) || Double.isNaN(mean))){
              throw new RuntimeException("!!");
          }
         this.updateVariance(variance);
          if(Constants.CHECK && (scale==0 || Double.isNaN(this.scale) || Double.isNaN(this.scale))){
              throw new RuntimeException("!! " +variance);
          }
          }
          else{
        	  this.scale = this.stddevPrior;
          }
          this.location = mean;
         
          this.recalcName();
       /* else{
            Logger.global.info("sum is too small!!! "+sum+" "+getMean()+" "+getStdDev());
        }*/
    }
    public void variance( double[] sum){
   	 double totx =0;
   	 double expx =0;
   	 for(int i=0; i<obsv.size(); i++){
            double w = obsv.get(i);
            expx+= Math.pow(this.obsx.get(i) - this.location,2)*w ;
              totx+=w;
          
          }
        sum[0] +=expx;
        sum[1]+=totx;
         
   }
    
public void maximiseVariance(){
    double mean =getMean();//average(pseudo);
    double variance = variance(mean, 0);
    if(Double.isNaN(variance) || Double.isNaN(mean)){
        throw new RuntimeException("!!");
    }
    Logger.global.info("after "+mean+" "+variance);
    this.updateVariance(variance);
}
  
    /*public String toString(){
        Double[] d = new Double[] {getMean(), getStdDev(), 0.0};
        return 
        Format.sprintf("c(%5.2g, %5.2g, %5.2g)", d);
//        "c("+getMean()+","+getStdDev()+",0);";
    }*/
    public void updateVariance(double var) {
        this.scale = Math.sqrt(var);
        
    }
    public double probability(double r, double offset) {
      return normal.pdf(r-offset, location, scale);
//       normal.setState(mean-offset, stddev);
    }
    public double rsn() {
        return this.inverse( Constants.rand.nextDouble());
    }
   @Override 
   public int getNumArguments(){
       return 2;
   }

public void getCi(Double[] ci, double[] mm) {
    // TODO Auto-generated method stub
    
}
public double[] getAngle(int num) {
   double[] res = new double[num];
   double incr = 1.0 / (double)num;
   for(int i=0; i<res.length; i++){
       res[i] = i==res.length-1 ? this.inverse(0.99999999):
           this.inverse((i+1)*incr);
       
   }
   return res;
}


public double evaluate(double[] argument) {
	throw new RuntimeException("!!");
}
public double getLowerBound(int n) {
	throw new RuntimeException("!!");
}

public OrthogonalHints getOrthogonalHints() {
	throw new RuntimeException("!!");
}
public double getUpperBound(int n) {
	throw new RuntimeException("!!");
}
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
		    	if(Math.abs(pseudo)>1e-5){
		    		double v = Math.abs(pseudo);
		    		int k = numObs+i;
		    		//x.setQuick(k, 0, 1.0*v);
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
	double avg = this.location;// this.average(0.0);
	//double sum =0;
	//double cnt=0;
	 for( i=0; i<this.obsx.size(); i++){
		 int k = numObs+i;
		 double v = this.obsv.get(i);
		double diff = this.obsx.get(i) -avg;
		//sum+=v*Math.pow(diff, 2);
		//cnt+=v;
		 y.setQuick(k,0, v*( Constants.transformVariance(diff)));
	 }
	// double av = Math.sqrt(sum/cnt);
	 if(Math.abs(pseudo)>1e-5){
		 int k = numObs+i;
		 double v = Math.abs(pseudo);
		double vv = Constants.transformVariance(pseudo<0 ? this.stddevPrior : this.scale);
			if(Constants.CHECK && (Double.isNaN(vv) || Double.isNaN(vv))){
	   			throw new RuntimeException("!!");
	   		}
		 y.setQuick(k,0, v*( vv));
		 i++;
	 }
	 return i;
}
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
public void transfercounts(EmissionState innerState, int phen_index, int i){
    throw new RuntimeException("!!");
    /* for(Iterator<Map.Entry<Double, Double>> it = this.observations.entrySet().iterator(); it.hasNext();){
        Map.Entry<Double, Double> nxt = it.next();
             innerState.addCountDT(nxt.getKey(),phen_index,  nxt.getValue(), i);
     *}*/
 }
public double calcLH() {
	throw new RuntimeException("!!");
}
public double getStdDev() {
	throw new RuntimeException("!!");
}
public void setMinMax(double min, double max){
	//throw new RuntimeException("!!");
}
}

