package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;

import lc1.dp.states.EmissionState;

import org.jfree.data.xy.XYSeries;

import pal.math.MultivariateFunction;
import cern.colt.matrix.DoubleMatrix2D;


/** interface for a distribution over doubles */
public interface ProbabilityDistribution  extends Comparable, MultivariateFunction,  Serializable {

    public double sample();

    public double probability(double x);
  
    public void addCount(double obj_index, double value);

    public ProbabilityDistribution clone();
    public ProbabilityDistribution clone(double u);
    
    public void transfercounts(EmissionState innerState, int phen_index, int i);

    public void initialise();

    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp);

    public void transfer(double pseudoC);

    public String id();
    
    public void addCounts(ProbabilityDistribution probabilityDistribution);

    public double[] getCount(double[] angle);

	public int getParamIndex();

	public String name();

	public double getMean();

	public double sum();

	public void maximise(double d, double e, double f);

	public void updateParamIndex();

	public void recalcName();

	public double prior();

	public double getParamValue(int n1);

	public void setParamValue(int n1, double val);

	public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch);

	public void plotTheoretical(String string, boolean b, XYSeries theor);

	public void getInterval(double[] in, double[] res_inR);

	public double scale();

	public int numObs();

	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs, double[] noCop,   double pseudo);

	//public void setPriorMean(double d);

	public void variance(double[] sum);

	//public void setPriorVar(double e);

	public int fillVariance(DoubleMatrix2D y,  int numObs, double pseudo);

	public void setParam(int type, double rho);

	public void print(PrintWriter pw);

	public void setPriors(ProbabilityDistribution distx, int type);

	public void setMinMax(double min, double max);

	public double  probability(double x, int mixComponent);

	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1);

	public void addCount(Double b, double val,
			SimpleExtendedDistribution1 mixe1, ProbabilityDistribution disty);

	public void setCoverage(double d);


}
