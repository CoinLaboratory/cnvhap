package lc1.stats;

import java.io.PrintWriter;

import lc1.dp.states.EmissionState;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import JSci.maths.statistics.ProbabilityDistribution;
import cern.colt.matrix.DoubleMatrix2D;

public class UniformDistribution extends ProbabilityDistribution implements lc1.stats.ProbabilityDistribution{
     double prob;
    public UniformDistribution(double min, double max){
        this.min = min;
        this.max = max;
        prob = 1.0/ (max- min);
       // if(max==min) throw new RuntimeException("!!");
    }
    @Override
	public int numObs() {
		return 0;
	}
    public void setCoverage(double d){
    	throw new RuntimeException("!!");
    }
    public void setPriors(lc1.stats.ProbabilityDistribution distx, int type){
	}
    public double  probability(double x, int mixComponent){
		return this.probability(x);
	}
    public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, lc1.stats.ProbabilityDistribution disty){
	//	throw new RuntimeException("!!");
	}
    public lc1.stats.ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1){
		return this.clone(u);
	}
    public void variance( double[] sum){
	        
	  }
    public void print(PrintWriter pw){
    	pw.print(this.toString()+"\t");
    }
    public int fillVariance(DoubleMatrix2D y,  int numObs, double pseudo){
    	return -1;
    }
    public void setMinMax(double min, double max){
  //  if(true) throw new RuntimeException("should be immutable");
    	this.min = min;
    	this.max = max;
    	this.prob = 1.0/(max-min);//*/
    }
    double min, max;
    @Override
    public double cumulative(double arg0) {
        if(arg0<min) return 0;
        else if(arg0>max) return 1;
        else return (arg0 - min )/max;
    }

    @Override
    public double inverse(double arg0) {
        return arg0*(max - min) + min;
    }

    @Override
    public double probability(double arg0) {
        if(arg0<min) return 0;
        else if(arg0>max) return 0;
        else return prob;
    }

	public void addCount(double obj_index, double value) {
		// TODO Auto-generated method stub
		
	}

	public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution) {
		// TODO Auto-generated method stub
		
	}

	public lc1.stats.ProbabilityDistribution clone(double u) {
		//return this;//
		return new UniformDistribution(min, max);
	}
	public lc1.stats.ProbabilityDistribution clone() {
		return this;
//		return new UniformDistribution(min, max);
	}

	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}

	public double getMean() {
		// TODO Auto-generated method stub
		return (max - min)/2.0;
	}
	
	public String toString(){
    	return "U:"+min+"-"+max+"_"+prob;
    }

	public int getParamIndex() {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getParamValue(int n1) {
		// TODO Auto-generated method stub
		return 0;
	}

	public String id() {
		// TODO Auto-generated method stub
		return null;
	}

	public void initialise() {
		// TODO Auto-generated method stub
		
	}

	public void maximise(double d, double e, double f) {
		// TODO Auto-generated method stub
		
	}

	public String name() {
		// TODO Auto-generated method stub
		return null;
	}

	public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch) {
		// TODO Auto-generated method stub
		return 0;
	}

	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		// TODO Auto-generated method stub
		
	}

	public double prior() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void recalcName() {
		// TODO Auto-generated method stub
		
	}

	public double sample() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void setParamValue(int n1, double val) {
		// TODO Auto-generated method stub
		
	}

	public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution[] tmp) {
		// TODO Auto-generated method stub
		
	}

	public double sum() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void transfer(double pseudoC) {
		// TODO Auto-generated method stub
		
	}

	public void transfercounts(EmissionState innerState, int phen_index, int i) {
		// TODO Auto-generated method stub
		
	}

	public void updateParamIndex() {
		// TODO Auto-generated method stub
		
	}

	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double evaluate(double[] argument) {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getLowerBound(int n) {
		// TODO Auto-generated method stub
		return 0;
	}

	public int getNumArguments() {
		// TODO Auto-generated method stub
		return 0;
	}

	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}
	 
	public double getUpperBound(int n) {
		// TODO Auto-generated method stub
		return 0;
	}
	public void  getInterval(double[] input, double[] res) {
		for(int i=0; i<res.length; i++){
		res[i] =  input[i]*(max-min)+min;
		}
	  //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
	}
	
	public double scale(){
    	return (max-min);
    }
	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[]  noCop,
		double pseudo) {
		
		return -1;
	}

	public void setPriorMean(double d){
		
	}
	@Override
	public void setParam(int type, double rho) {
		// TODO Auto-generated method stub
		
	}
	
}
