package lc1.stats;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import lc1.dp.illumina.RegressParams;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
public class SkewedUniformDistribution implements ProbabilityDistribution2{
	
	public static void main(String[] args){
		
	}
	
	Comparator<Double[]> compd = new Comparator<Double[]>(){

		@Override
		public int compare(Double[] o1, Double[] o2) {
			int res = o1[0].compareTo(o2[0]);
			if(res!=0) return res;
			else return o1[1].compareTo(o1[1]);
		}
		
	};
	public void addCount(Double r2,  Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG){
		throw new RuntimeException("!!");
	}
	
	
	//DoubleMatrix2D covariance =  new DenseDoubleMatrix2D(2,2);
	List<Double> x_obs = new ArrayList<Double>();
	List<Double> y_obs = new ArrayList<Double>();
	List<Double> weight = new ArrayList<Double>();
	// public SortedMap<Double[], Double> observations = new TreeMap<Double[], Double>(compd);
public 	 double sigma_x;
	
final double miny,maxy;
final double y_prob;

	private double rho;
	public double rho(){
		return rho;
	}
	 double meanx;
    public SkewedUniformDistribution clone(){
        return new SkewedUniformDistribution(this);
    }
    public SkewedUniformDistribution clone(double u){
        return new SkewedUniformDistribution(this, u);
    }
    @Override
	public void print(PrintWriter pw) {
		pw.println(this.toString());
		
	}
    int paramIndex =0;
    public int getParamIndex(){
        return this.paramIndex;
   }
    @Override
    public void setParam(int type, int i, double d){
    	throw new RuntimeException("!!");
    	
	}
    
    @Override
	public int numObs() {
		return this.x_obs.size();
	}
    public void updateParamIndex(){
        this.paramIndex++;
    }
    public int compareTo(Object o) {
		return this.name.compareTo(((SkewedUniformDistribution)o).name);
	}
public int size(){
    return this.weight.size();
}
public String name;
String id;
public void recalcName() {
	name =   id+"_"+String.format("%5.3g,", this.meanx)+
	String.format("%5.3g", sigma_x);
   
}

public SkewedUniformDistribution(String name, double meanx, double stddevx, double stddevxprior,
		double miny, double maxy
		){
	this.id = name;
	this.name = name;
this.sigma_x = stddevx;
this.miny = miny;
this.maxy = maxy;
this.y_prob = 1.0/ (maxy - miny);
this.meanx = meanx;

this.rho = 0;
this.rhoprior=0;
this.stddevpriorx = stddevxprior;

this.meanPriorx = meanx;
x0.set(0, 0, this.meanPriorx);
x0.set(1,0,0);
   
}



public SkewedUniformDistribution(SkewedUniformDistribution trainableNormal) {
	this.id = trainableNormal.id;
	this.sigma_x = trainableNormal.sigma_x;

	this.meanx = trainableNormal.meanx;

	this.rho = trainableNormal.rho;
	this.rhoprior=trainableNormal.rhoprior;
	this.stddevpriorx = trainableNormal.stddevpriorx;
	this.meanPriorx =trainableNormal.meanPriorx;
	this.x0 = trainableNormal.x0;
	this.y_prob = trainableNormal.y_prob;
	this.maxy = trainableNormal.maxy;
	this.miny = trainableNormal.miny;
}

public SkewedUniformDistribution(SkewedUniformDistribution skewNormal, double u) {
    this(skewNormal);
   // this.setPrior(u, u, u);
   
  }
  
public double dsn (double x, double y, boolean log)
{
 
 return log ? this.probabilityLog(x,y) : this.probability(x,y);
}



private double meanx(double b){
	return this.meanx + this.rho*b;
}
static lc1.stats.NormalDistribution normal;

    @Override
    public double probability(double x1, double y1) {
    //	double mu = meanx(y1);
    	double res = normal.pdf(x1, meanx(y1), this.sigma_x) *this.y_prob;
    	return res;
    	
  
    }
    
  
    public double probabilityLog(double x1, double y1) {
    	return normal.logpdf(x1, meanx(y1), this.sigma_x) *this.y_prob;
    }

  double meanPriorx, stddevpriorx, rhoprior;
    
  
  
   
    
    public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yB, DoubleMatrix2D covar, int numObs, double pseudo){
    	int i;
    	throw new RuntimeException("exc");
    }
    
   
    double sum;
    public String getObsString(){
         if(size()>0) return median()+"";
        else return "-";
    }
    public Double[] median(){
     
        double sum=0;
        double tot = sum;
        int mid =(int) Math.floor( (double)this.weight.size()/2.0);
        return new Double[] {this.x_obs.get(mid), this.y_obs.get(mid)};
    }
//static final double cntThresh = Constants.countThresh();
    public synchronized void addCount(double x, double y, double w) {
  
      //  if(w > cntThresh){
        this.addObservation(x,y, w);
     //   }
  
        }
    
    public void addObservation(double d1,double d2,  double weight){
        //double d1 = round(d);
    	this.x_obs.add(d1);
    	this.y_obs.add(d2);
    	this.weight.add(weight);
    
        sum+=weight;
    }
    
     DoubleMatrix2D x0 =  new DenseDoubleMatrix2D(2,1);
	static  DoubleMatrix2D Q = new DenseDoubleMatrix2D(2,2);
   
	

    public void initialiseCounts() {
    	this.x_obs.clear();
    	this.y_obs.clear();
    	this.weight.clear();
//        observations.clear();
        sum =0;
    }
        
  
public double variance(double mu, double pseudo){
    	
        double tot =pseudo;
      double exp =Math.pow(this.stddevpriorx,2)*pseudo;
      for(int i=0; i<this.x_obs.size(); i++){
      
          double sc =Math.pow(x_obs.get(i)-this.meanx(y_obs.get(i)), 2);
          double w =this.weight.get(i);
            exp+=sc*w;
            tot+=w;
          //  if(exp < 0 || tot < 0) throw new RuntimeException("!! "+w+" "+sc);
        }
        if(Constants.CHECK && Double.isInfinite(exp/tot) || tot==0){
            throw new RuntimeException(" ");//+this.observations+"\n"+this.observations);
        }
       
        return exp/tot;
    }
public void maximise(double pseudo,  double pseudoSD, double pseudoSkew, 
		double pseudo1, double pseudoSD1, double pseudoSkew1, double pseudoCov
   
    		) {
    	if(sum<Constants.trainThresh()) {
    		return;
    	}
    	DoubleMatrix2D lrr = new DenseDoubleMatrix2D(this.y_obs.size(),1);
    	DoubleMatrix2D baf = new DenseDoubleMatrix2D(this.y_obs.size(),2);
    	
    	for(int i=0; i<this.x_obs.size(); i++){
    		 int k = i;
    		 double v = this.weight.get(i);
    		
    		 baf.setQuick(k, 0, 1.0*v);
    		
    		 baf.setQuick(k,1, v*( this.y_obs.get(i) ));
    		 lrr.setQuick(k,0, v*( this.x_obs.get(i)));
    	 }
        	this.Q.set(0, 0, Math.pow(pseudo,2));
        	this.Q.set(1, 1, Math.pow(pseudoCov,2));
        DoubleMatrix2D res = RegressParams.solve(baf, lrr, Q, x0);
        this.meanx = res.get(0, 0);
        this.rho = res.get(1,0);
	          double variance = variance(meanx, pseudoSD);
	          this.sigma_x = Math.sqrt(variance);
          this.recalcName();
    }

  
   
    
    
    

  
   
   public String toString(){
	   return this.meanx+"_U:"+this.miny+"-"+this.maxy;
   }
  


public double[] getAngle(int num) {
	 throw new RuntimeException("!!");
  /* double[] res = new double[num];
   double incr = 1.0 / (double)num;
   for(int i=0; i<res.length; i++){
       res[i] = i==res.length-1 ? this.inverse(0.99999999):
           this.inverse((i+1)*incr);
       
   }
   return res;*/
}
@Override
public String id() {
	return id;
}
@Override
public void initialise() {
this.initialiseCounts();
	
}
@Override
public String name() {
	return name;
}

public void transfer(double ps){
	if(true) throw new RuntimeException("!!");
	if(this.size()==0) return;
    this.maximise(ps, ps, ps,ps,ps,ps,ps);
    
}

@Override
public void getInterval(double[] input, DoubleMatrix2D cov, double[] mean) {
	double sigma_y = this.maxy - this.miny;
	cov.setQuick(0, 0, Math.pow(this.sigma_x, 2));
	cov.setQuick(1, 1, Math.pow(sigma_y, 2));
	double offd = sigma_x*sigma_y*rho;
	cov.setQuick(0, 1, offd);
	cov.setQuick(1, 0, offd);
	mean[0] = this.meanx;
	mean[1] = (this.maxy+this.miny)/2.0;
}
  @Override
public int fill(DoubleMatrix2D x, DoubleMatrix2D y, DoubleMatrix2D yB, int numObs, double[] noCop,  double pseudo){
	int i;
	if(true) throw new RuntimeException("exc");
 for( i=0; i<this.x_obs.size(); i++){
	 int k = numObs+i;
	 double v = this.weight.get(i);
	
	 for(int kk=0; kk<noCop.length; kk++){
		 x.setQuick(k, kk, v * (double)noCop[kk]);
	 }
	 
	 y.setQuick(k,0, v*( this.x_obs.get(i) ));
	 yB.setQuick(k,0, v*( this.y_obs.get(i)));
 }
 return i;
}
public void setRho(double v) {
	this.rho = v;
	
}

public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
){
	throw new RuntimeException("!!");
	}
@Override
public void variance(int type, double[] sum) {
	throw new RuntimeException("!!");
	
}

@Override
public ProbabilityDistribution2 clone(double u,
		SimpleExtendedDistribution1 dist1) {
	throw new RuntimeException("!!");
}
@Override
public double probability(double r, double b, int mixComponent) {
	// TODO Auto-generated method stub
	throw new RuntimeException("!!");
}
@Override
public void setMinMax(double minR, double maxR, double min, double max) {
	throw new RuntimeException("!!");
	
}
@Override
public void setToExclude() {
	throw new RuntimeException("!!");
	
}

}
