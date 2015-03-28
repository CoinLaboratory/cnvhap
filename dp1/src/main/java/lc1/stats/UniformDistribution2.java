package lc1.stats;

import java.io.PrintWriter;

import lc1.dp.states.EmissionState;

import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix2D;

public class UniformDistribution2  implements lc1.stats.ProbabilityDistribution2{
    double probx, proby;
    public String toString(){
    	return "U:"+minx+"-"+maxx+"_"+miny+"-"+maxy+"_"+probx+"_"+proby;
    }
    @Override
	public void print(PrintWriter pw) {
		pw.println(this.toString());
		
	}
    
    public void setPriors(ProbabilityDistribution2 pr2, int type, boolean x
	){
   	throw new RuntimeException("!!");
   	}
    public void addCount(Double r2,  Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG){
		throw new RuntimeException("!!");
	}
	
    public UniformDistribution2(double minx, double maxx, double miny, double maxy){
        this.minx = minx;
        this.maxx = maxx;
        probx = 1.0/ (maxx- minx);
        this.miny = miny;
        this.maxy = maxy;
        proby = 1.0/ (maxy- miny);
       // if(max==min) throw new RuntimeException("!!");
    }
    
    public UniformDistribution2(UniformDistribution2 un) {
    	this(un.minx, un.maxx, un.miny, un.maxy);
	}
	@Override
	public int numObs() {
		return 0;
	}
    @Override
    public int fill(DoubleMatrix2D x, DoubleMatrix2D y,DoubleMatrix2D yb, int numObs,double[] noCop, double pseudo){
		return 0;
	}
    
    public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
			DoubleMatrix2D covar, int numObs, double pseudo){
	return 0;
 }
    public void setMinMax(double minx, double maxx, double miny, double maxy){
    	this.minx = minx;
    	this.maxx = maxx;
    	this.miny = miny;
    	this.maxy = maxy;
    	this.probx = 1.0/(maxx-minx);
    	this.proby = 1.0/(maxy-miny);
    }
    double minx, maxx, miny, maxy;
    /*@Override
    public double cumulative(double arg0) {
        if(arg0<min) return 0;
        else if(arg0>max) return 1;
        else return (arg0 - min )/max;
    }

    @Override
    public double inverse(double arg0) {
        return arg0*(max - min) + min;
    }*/

    @Override
    public double probability(double argx, double argy) {
        if(argx<minx || argx>maxx) return 0;
        else if(argy < miny || argy>maxy) return 0;
        else return probx*proby;
    }

	public void addCount(double obj_index, double value, double valuey) {
		// TODO Auto-generated method stub
		
	}

	public void addCounts(lc1.stats.ProbabilityDistribution2 probabilityDistribution) {
		// TODO Auto-generated method stub
		
	}

	public lc1.stats.ProbabilityDistribution2 clone(double u) {
		return new UniformDistribution2(minx, maxx, miny, maxy);
	}
	public lc1.stats.ProbabilityDistribution2 clone() {
		return new UniformDistribution2(minx, maxx, miny, maxy);
	}

	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}

	public double[] getMean() {
		// TODO Auto-generated method stub
		return new double[] {0,0};
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

	public void maximise(double d, double e, double f, double d1, double e1, double f1, double g) {
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

	public void setParamsAsAverageOf(lc1.stats.ProbabilityDistribution2[] tmp) {
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
	public void  getIntervalX(double[] input, double[] res) {
		for(int i=0; i<res.length; i++){
		res[i] =  input[i]*(maxx-minx)+minx;
		}
	  //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
	}
	
	public void  getInterval(double[] input, DoubleMatrix2D matr, double[] mean) {
	mean[0] = (maxx - minx)/2.0;
	mean[1] = (maxy - miny)/2.0;
	matr.set(0, 0, maxx - minx);
	matr.set(1, 1, maxy - miny);
		
	}
	public void  getIntervalY(double[] input, double[] res) {
		for(int i=0; i<res.length; i++){
		res[i] =  input[i]*(maxy-miny)+miny;
		}
	  //  res[2] = normal.quantile(greaterThan, this.location, this.scale);
	}

	public void setParam(int type, int i,double d){
		
	}
	@Override
	public void variance( int type, double[] sum) {
		// TODO Auto-generated method stub
		
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
	public void setToExclude() {
		throw new RuntimeException("!!");
		
	}
	

}
