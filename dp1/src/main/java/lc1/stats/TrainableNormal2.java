package lc1.stats;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix2D;
public class TrainableNormal2 implements ProbabilityDistribution2{
	
	public static void main(String[] args){
		
	}
	@Override
	public void print(PrintWriter pw) {
		pw.println(this.toString());
		
	}
	public void addCount(Double r2,  Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG){
		mixe1.addCount(0, val);
		probDistG.addCount(r2, b, val);
	}
	
	public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
	){
   	throw new RuntimeException("!!");
   	}
	Comparator<Double[]> compd = new Comparator<Double[]>(){

		@Override
		public int compare(Double[] o1, Double[] o2) {
			int res = o1[0].compareTo(o2[0]);
			if(res!=0) return res;
			else return o1[1].compareTo(o1[1]);
		}
		
	};
	
	//DoubleMatrix2D covariance =  new DenseDoubleMatrix2D(2,2);
	List<Double> x_obs = new ArrayList<Double>();
	List<Double> y_obs = new ArrayList<Double>();
	List<Double> weight = new ArrayList<Double>();
	// public SortedMap<Double[], Double> observations = new TreeMap<Double[], Double>(compd);
public 	 double sigma_x;
	 public double sigma_y;
	private double rho;
	public double rho(){
		return rho;
	}
	 double meanx, meany;
    public TrainableNormal2 clone(){
        return new TrainableNormal2(this);
    }
    public TrainableNormal2 clone(double u){
        return new TrainableNormal2(this, u);
    }
    int paramIndex =0;
    public int getParamIndex(){
        return this.paramIndex;
   }
    
    public void setParam(int type, int i, double d){
    	if(Double.isNaN(d)){
    		throw new RuntimeException("!!");
    	}
    	if(type==0){
    		if(i==0 ) {
    			this.meanPriorx = d;
    			this.meanx = d;
    		}
    		else{
    			this.meanPriory = d;
    			this.meany = d;
    		}
    	}
    	else if(type==1){
    		if(i==0){
    			double d1 = Math.sqrt(Math.max(1e-10,d));
    			this.stddevpriorx = d1;
    			this.sigma_x = d1;
    		}
    		else if(i==1){
    			double d1 =Math.sqrt(Math.max(1e-10,d));
    		//	if(//Math.abs(x+0.04397)<0.001 && Math.abs(y-0.4296303)<0.001 && 
    				//	(id.indexOf("_0.33")>=0|| id.indexOf("_0.66")>=0)){
    	    	//	System.err.println(id  +this.meanx+" "+this.meany+" "+this.sigma_x+" "+this.sigma_y);
    	    	//}
    			this.stddevpriory = d1;
    			this.sigma_y =d1;
    		}
    		else{
    		 double d_ = d/(this.sigma_x*sigma_y);
    		 d_ = Math.signum(d_)*Math.min(Math.abs(d_), Constants.rhoMax());
    			this.rho = d_;
    			this.rhoprior = d_;
    		
    	//	if(Math.abs(rho)>0.1 && this.sum<5){
    	//			throw new RuntimeException("!!");
    	//		}
    		}
    	}
    	
	}
    
    @Override
	public int numObs() {
		return this.x_obs.size();
	}
    public void updateParamIndex(){
        this.paramIndex++;
    }
/*@Override
public void  getInterval(double[] input, double[] res) {
	for(int i=0; i<res.length; i++){
	res[i] =  normal.quantile(input[i], this.location, this.scale);
	}
  //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
}*/
    public int compareTo(Object o) {
		return this.name.compareTo(((TrainableNormal2)o).name);
	}
public int size(){
    return this.weight.size();
}
public String name;
String id;
public void recalcName() {
	name =   id+"_"+String.format("%5.3g,", this.meanx)+
	String.format("%5.3g,", meany)+String.format("%5.3g", sigma_x)+String.format("%5.3g", sigma_y);
   
}

public TrainableNormal2(String name, double meanx, double stddevx, double stddevxprior,
		double meany, double stddevy, double stddevyprior
		){
	this.id = name;
	this.name = name;
this.sigma_x = stddevx;
this.sigma_y = stddevy;
if(Double.isNaN(sigma_y) || Double.isNaN(sigma_x)){
	  throw new RuntimeException("!!");
}
this.meanx = meanx;
this.meany = meany;
this.rho = 0;
this.rhoprior=0;
this.stddevpriorx = stddevxprior;
this.stddevpriory = stddevyprior;
this.meanPriorx = meanx;
this.meanPriory = meany;

    //super(name, mean, stddev, 0.0, 
     //       new double[] {normal.quantile(1e-3, mean, stddev), 0.001, 0},  new double[] { normal.quantile(0.999, mean, stddev), 1e3, 0},
      //      round, priorMod
       //    );
   
}



public TrainableNormal2(TrainableNormal2 trainableNormal) {
	this.id = trainableNormal.id;
	this.sigma_x = trainableNormal.sigma_x;
	this.sigma_y = trainableNormal.sigma_y;
	this.meanx = trainableNormal.meanx;
	this.meany = trainableNormal.meany;
	this.rho = trainableNormal.rho;
	this.rhoprior=trainableNormal.rhoprior;
	this.stddevpriorx = trainableNormal.stddevpriorx;
	this.stddevpriory = trainableNormal.stddevpriory;
	this.meanPriorx =trainableNormal.meanPriorx;
	this.meanPriory = trainableNormal.meanPriory;

}

public TrainableNormal2(TrainableNormal2 skewNormal, double u) {
    this(skewNormal);
   // this.setPrior(u, u, u);
   
  }
  
public double dsn (double x, double y, boolean log)
{
 
 return log ? this.probabilityLog(x,y) : this.probability(x,y);
}




 

   
static lc1.stats.NormalDistribution normal;
@Override
    public double probability(double x1, double y1) {
    	if(Constants.suppressB()){
    		return normal.pdf(x1, this.meanx, this.sigma_x);
    	}
    	else if(Constants.suppressR()){
    		return normal.pdf(y1, this.meany, this.sigma_y);
    	}
    	double x = x1-this.meanx;
    	double y = y1 - this.meany;
    	double res = (1.0 / (2*Math.PI*sigma_x*sigma_y*Math.sqrt(1-Math.pow(rho,2))))*
    	Math.exp((-1/(2*(1-Math.pow(rho,2))))*(
    			Math.pow(x,2)/Math.pow(sigma_x, 2) + 
    			Math.pow(y,2)/Math.pow(sigma_y, 2) -
    			(2*rho*x*y)/(sigma_y*sigma_x)
    			)
    			);
    	if(Double.isNaN(res)){
    		double logp = this.probabilityLog(x1, y1);
    		throw new RuntimeException("!!");
    	}
    	return res;
  
    }
    
  
    public double probabilityLog(double x1, double y1) {
    	double x = x1-this.meanx;
    	double y = y1 - this.meany;
    	double resa  = (2*Math.PI*sigma_x*sigma_y*Math.sqrt(1-Math.pow(rho,2)));
    	double resb = -1/(2*(1-Math.pow(rho,2)));
    	double resc = Math.pow(x,2)/Math.pow(sigma_x, 2);
    	double resd =	Math.pow(y,2)/Math.pow(sigma_y, 2);
    	double rese = (2*rho*x*y)/(sigma_y*sigma_x);
    	double res = 
    	
    	Math.log(1.0 / resa)+
    	((resb)*(
    			resc + resd
    		 -rese
    			
    			)
    			);
    	if(Double.isInfinite(res) || res==Double.NEGATIVE_INFINITY){
    		throw new RuntimeException("!!");
    	}
  return res;
    }

  double meanPriorx, meanPriory, stddevpriorx, stddevpriory, rhoprior;
    
   //double tot=0;
   //double exp=0;
    public double[] average(double pseudox, double pseudoy){
        double totx =pseudox;
        double toty = pseudoy;
        double expx =pseudox * this.meanPriorx;
        double expy = pseudoy *this.meanPriory;
        for(int i=0; i<weight.size(); i++){
         //   Map.Entry<Double[], Double> nxt = it.next();
          //  Double[] sc = nxt.getKey();
            double w = weight.get(i);//nxt.getValue();
           
            expx+=this.x_obs.get(i)*w;
            expy+=this.y_obs.get(i)*w;
            totx+=w;
            toty+=w;
        }
       // if(Math.abs(mean- Constants.r_mean()[3])<0.01){
        //   System.err.println(this.meanPrior.getMean()+" rvalues "+observations);
        //   System.err.println(this.meanPrior.getMean()+" rvalues "+weights);
       // }
        return new double[] {expx/totx, expy/toty};
    }
  
    public final double[] variance(double mux, double muy, double pseudox, double pseudoy,double pseudoxy){
    	
        double totx =0;//pseudox;
        double toty = 0;//pseudoy;
        double totxy =0;// pseudoxy; 
      
      double expx =0;//;
      double expy =0;//;
   
      double expxy =0;
      
     // check(expx==0 ? 0 : expx/totx, expy==0 ? 0 :expy/toty, expxy==0 ? 0 : expxy/totxy);
      if(true){
      for(int i=0; i<weight.size(); i++){
         
          double w = weight.get(i);
          double x = this.x_obs.get(i);
          double y = this.y_obs.get(i);
            expx+=Math.pow(x - mux,2)*w;
            expy+=Math.pow(y - muy,2)*w;
            expxy+=(x - mux)*(y - muy)*w;
            totx+=w;
            toty+=w;
            totxy+=w;
        //    check(expx/totx, expy/toty, expxy/totxy);
          //  if(exp < 0 || tot < 0) throw new RuntimeException("!! "+w+" "+sc);
        }
      }
       // if(Constants.CHECK && Double.isInfinite(exp/tot) || tot==0){
       //     throw new RuntimeException(" "+this.observations+"\n"+this.observations);
      //  }
       double sigma_x1 = Math.sqrt((expx+Math.pow(this.stddevpriorx,2)*pseudoxy)/(totx+pseudoxy));
       double sigma_y1 = Math.sqrt((expy+Math.pow(this.stddevpriory,2)*pseudoxy)/(toty+pseudoxy));
       double sigma_xy1 = (expxy+0)/(totxy+pseudoxy);
        return new double[] {(expx+Math.pow(this.stddevpriorx,2)*pseudox)/(totx+pseudox), 
        		(expy+Math.pow(this.stddevpriory,2)*pseudoy)/(toty+pseudoy), 
        		sigma_xy1 /( sigma_x1*sigma_y1)};
    }
    
    
    public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yB, DoubleMatrix2D covar, int numObs, double pseudo){
    	int i;
    	double[] avg = this.average(0, 0);
     for( i=0; i<this.x_obs.size(); i++){
    	 int k = numObs+i;
    	 double v = this.weight.get(i);
    	/* x.setQuick(k, 0, 1.0*v);
    	 x.setQuick(k, 1, v * (double)noCop);
    	 x.setQuick(k, 2, v*Math.pow((double) noCop,2));
    	 x.setQuick(k, 3, v * (double)bCop);
    	 x.setQuick(k, 4, v*Math.pow((double) bCop,2));*/
    	 y.setQuick(k,0, v*( Constants.transformVariance(this.x_obs.get(i) - avg[0]) ));
    	 yB.setQuick(k,0, v*(Constants.transformVariance(this.y_obs.get(i) - avg[1])));
    	 covar.setQuick(k,0, v*(x_obs.get(i) - meanx)*(y_obs.get(i) - meany));
     }
     if(pseudo>0){
    	 int k = numObs+i;
    	 double v = Math.abs(pseudo);
    	 y.setQuick(k,0, v*( Constants.transformVariance(pseudo < 0  ? this.stddevpriorx : this.sigma_x) ));
    	 yB.setQuick(k,0, v*(Constants.transformVariance(pseudo < 0 ? this.stddevpriory : this.sigma_y)));
    	 i++;
     }
     return i;
    }
    
    public void variance(int type, double[] sum){
    	 double totx =0;
    	 double expx =0;//Math.pow(type==0 ? this.stddevpriorx : (type==1 ? this.stddevpriory : this.rhoprior),2)*pseudox;
    	 for(int i=0; i<weight.size(); i++){
    		
             double w = weight.get(i);
             
             expx+=type ==0 ? Math.pow(this.x_obs.get(i) - this.meanx,2)*w :
            	 (type==1 ?Math.pow(this.y_obs.get(i) - this.meany,2)*w :
            		 (x_obs.get(i) - meanx)*(y_obs.get(i) - meany)*w
            		 );
               totx+=w;
           
           }
    	 
         sum[0] +=expx;
         sum[1]+=totx;
         //if(Math.sqrt(sum[0]/sum[1]) > 5){
    		// double sigma = Math.sqrt(expx/totx);
    		 //Logger.global.info("h");
    	 }
    //}
    double sum;
    public String getObsString(){
         if(size()>0) return median()+"";
        else return "-";
    }
    public Double[] median(){
     
        double sum=0;
        double tot = sum;
      //  StringBuffer sb = new StringBuffer(tot+"_sum ");
       // Histogram hist = new Histogram(observations.firstKey(), observations.lastKey(), 10);
       /* for(Iterator<Double[]> it = this.observations.keySet().iterator(); it.hasNext();){
           Double[] d =  it.next();
           sum+=observations.get(d)/tot;
           if(sum>0.5) return d;
        //   sb.append("("+d+":"+sum+"),");
        }
        return observations.lastKey();//        return sb.toString();*/
        int mid =(int) Math.floor( (double)this.weight.size()/2.0);
        return new Double[] {this.x_obs.get(mid), this.y_obs.get(mid)};
    }
static final double cntThresh = Constants.countThresh();
    public synchronized void addCount(double x, double y, double w) {
      //  if(!Double.isNaN(d) &&  !Double.isInfinite(d) 
          //      && d!=Double.NEGATIVE_INFINITY
      //          && d!=Double.POSITIVE_INFINITY && w > Constants.countThresh()){
       // if(Constants.CHECK && (Double.isNaN(d) || Double.isInfinite(d))) throw new RuntimeException("!!");
     //   if( d < 0) System.err.println("less than zero "+d+" "+w+this.name);
       // double cntThresh = 
    	if(Constants.CHECK && Double.isNaN(w)){
        	throw new RuntimeException("!!");
        }
    //	if(true) throw new RuntimeException("need to modify this to be more like cut normal with addmincount and addmaxcount");
        if(w > cntThresh){
        this.addObservation(x,y, w);
        }
  //    if(Math.abs(mean- Constants.r_mean()[3])<0.01 && Math.abs(d +0.6068727)<0.00001){
   //      System.err.println("hj");
  //     }
        }
    
    public void addObservation(double d1,double d2,  double weight){
        //double d1 = round(d);
    	this.x_obs.add(d1);
    	this.y_obs.add(d2);
    	this.weight.add(weight);
    	 
    	//Double[] val = new Double[]{d1,d2};
      //  Double w = observations.get(val);
      //  observations.put(val, w==null ? weight : weight+w);
        sum+=weight;
    }
    
    
     
  /*  public void maximise(double pseudo, double pseudoSD, double pseudoSkew,
    		double pseudo1, double pseudoSD1, double pseudoSkew1) {
      this.maximise(pseudo, pseudoSD, pseudoSkew, Constants.trainThresh());
    }
 */
    public void initialiseCounts() {
    	this.x_obs.clear();
    	this.y_obs.clear();
    	this.weight.clear();
//        observations.clear();
        sum =0;
    }
        
    public void maximise(double pseudo,  double pseudoSD, double pseudoSkew, 
    		double pseudo1, double pseudoSD1, double pseudoSkew1, double pseudoCov
    		) {
    	//if(this.meany==1.0 || meany==0.0 || true) return;
    	//Logger.global.info("maximising "+this.toString()+" "+sum);
    //    Logger.global.info(this.weights+"");
    	if(sum<Constants.trainThresh()) {
    		this.meanx = this.meanPriorx;
    		this.meany = this.meanPriory;
    		this.sigma_x = this.stddevpriorx;
    		this.sigma_y = this.stddevpriory;
    		this.rho = this.rhoprior;
    		return;
    	}
    	double mean_x, mean_y;
    	if(pseudo>=1e4 && pseudo1>=1e4){
    		mean_x = this.meanPriorx;
    		mean_y = this.meanPriory;
    		
    	}
    	else{
    		double[] new_mean = average(pseudo, pseudo1);
    		mean_x = pseudo<1e4 ? new_mean[0] : this.meanPriorx;
    		mean_y = pseudo1<1e4 ? new_mean[1] : this.meanPriory;
    	}
    	  
        
          if(pseudoSD>=1e4 && pseudoSD1>=1e4 && pseudoCov>=1e4){
        	 this.sigma_x = this.stddevpriorx;
        	  this.sigma_y = this.stddevpriory;
        	
        	  this.rho = this.rhoprior;
          }
          else{
        	//  if(true) throw new RuntimeException("!!");
	          double[] variance = variance(mean_x,mean_y, pseudoSD, pseudoSD1, pseudoCov);
	        //  if(Constants.CHECK && (Double.isNaN(variance) || Double.isNaN(mean))){
	         //     throw new RuntimeException("!!");
	          //}
	          this.sigma_x = Math.sqrt(variance[0]);
	        //if(false && this.meany!=0.0 && this.meany!=1.0){
	        	
	        
	          this.sigma_y = Math.sqrt(variance[1]);
	      // }
	   
	    this.meanx = mean_x;
    // if(this.meany!=0.0 && this.meany!=1.0) 
    	 this.meany = mean_y;
        
       if(true){ this.rho = variance[2];
	          if(Math.abs(rho)>Constants.rhoMax()) {
	        	  rho = Math.min(Constants.rhoMax(), Math.max(rho, -Constants.rhoMax()));
	          }
          }
          }
        //  this.updateVariance(variance);
       
          this.recalcName();
       /* else{
            Logger.global.info("sum is too small!!! "+sum+" "+getMean()+" "+getStdDev());
        }*/
    }
    public void check(double var0, double var1, double var2){
    	double sigma_x = Math.sqrt(var0);
    	double sigma_y = Math.sqrt(var1);
    	double rho = var2 / (sigma_x*sigma_y);
    	System.err.println("rho is :       "+rho);
    //	if(Math.abs(rho)>1) {
        	//  rho =0;
	//  throw new RuntimeException("!! "+rho+" "+this.meanx+" "+this.meany+" "+this.sigma_x+" "+this.sigma_y);
       //   }
    }
    
    
    

  
    /*public String toString(){
        Double[] d = new Double[] {getMean(), getStdDev(), 0.0};
        return 
        Format.sprintf("c(%5.2g, %5.2g, %5.2g)", d);
//        "c("+getMean()+","+getStdDev()+",0);";
    }*/
   
   public String toString(){
	   return this.meanx+"_"+meany;
   }
  

public void getCi(Double[] ci, double[] mm) {
    // TODO Auto-generated method stub
    
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
	
	cov.setQuick(0, 0, Math.pow(this.sigma_x, 2));
	cov.setQuick(1, 1, Math.pow(this.sigma_y, 2));
	double offd = sigma_x*sigma_y*rho;
	cov.setQuick(0, 1, offd);
	cov.setQuick(1, 0, offd);
	mean[0] = this.meanx;
	mean[1] = this.meany;
}
  @Override
public int fill(DoubleMatrix2D x, DoubleMatrix2D y, DoubleMatrix2D yB, int numObs, double[] noCop,   double pseudo){
	int i;
	int len = noCop.length;
 for( i=0; i<this.x_obs.size(); i++){
	 int k = numObs+i;
	 double v = this.weight.get(i);
	 for(int kk=0; kk<len; kk++){
		 x.setQuick(k, kk, v * (double)noCop[kk]);
	 }
	 y.setQuick(k,0, v*( this.x_obs.get(i) ));
	 yB.setQuick(k,0, v*( this.y_obs.get(i)));
 }
 if(Math.abs(pseudo)>0){
	 int k = numObs+i;
	 double v = Math.abs(pseudo);
	// x.setQuick(k, 0, 1.0*v);
	 for(int kk=0; kk<len; kk++){
		 x.setQuick(k, kk, v * (double)noCop[kk]);
	 }
	 y.setQuick(k,0, v*( pseudo < 0 ? this.meanPriorx : this.meanx ));
	 yB.setQuick(k,0, v*(pseudo < 0 ? this.meanPriory :  this.meany));
	 i++;
 }
 return i;
}
public void setRho(double v) {
	this.rho = v;
	
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
//	throw new RuntimeException("!!");
	
}
@Override
public void setToExclude() {
	throw new RuntimeException("!!");
	
}
}
