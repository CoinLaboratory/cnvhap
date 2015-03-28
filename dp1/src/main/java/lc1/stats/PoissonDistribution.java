package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;

import lc1.dp.states.EmissionState;

import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import cern.colt.matrix.DoubleMatrix2D;

public class PoissonDistribution implements ProbabilityDistribution, Serializable {

	private double rate;
	static PoissonDistributionImpl poisson;
	
	
	public PoissonDistribution(double rate)
	{
		this.rate=rate;
		poisson = new PoissonDistributionImpl(rate);
	}

	public PoissonDistribution(PoissonDistribution poissonDistribution, double rate) {
		this(rate);
	}

	@Override
	public double probability(double x) {
		return poisson.probability(x);
	}

	@Override
	public double getMean() {
		return poisson.getMean();
	}
	
	public void setCoverage(double depth) {
		poisson.setMean(depth+0.0001);
	}
	
	@Override
	public void addCount(double objIndex, double value) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, ProbabilityDistribution disty) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void addCounts(ProbabilityDistribution probabilityDistribution) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public PoissonDistribution clone(){
		return new PoissonDistribution(rate);
	}
	
	@Override
	public PoissonDistribution clone(double rate){
		return null;
	}

	@Override
	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs,
			double[] noCop,  double pseudo) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int fillVariance(DoubleMatrix2D y, int numObs, double pseudo) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getCount(double[] angle) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void getInterval(double[] in, double[] resInR) {
		// TODO Auto-generated method stub
		
	}



	@Override
	public int getParamIndex() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getParamValue(int n1) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String id() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void initialise() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void maximise(double d, double e, double f) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String name() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int numObs() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double plotObservations(String string, boolean b, XYSeries obs,
			boolean swtch) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void plotTheoretical(String string, boolean b, XYSeries theor) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void print(PrintWriter pw) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double prior() {
		// TODO Auto-generated method stub
		return 0;
	}


	@Override
	public double probability(double x, int mixComponent) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void recalcName() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double sample() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double scale() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setMinMax(double min, double max) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParam(int type, double rho) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParamValue(int n1, double val) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setPriors(ProbabilityDistribution distx, int type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double sum() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void transfer(double pseudoC) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void transfercounts(EmissionState innerState, int phenIndex, int i) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void updateParamIndex() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void variance(double[] sum) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double evaluate(double[] arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getLowerBound(int arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getNumArguments() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public OrthogonalHints getOrthogonalHints() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getUpperBound(int arg0) {
		// TODO Auto-generated method stub
		return 0;
	}


	

}
