package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix2D;


/** interface for a distribution over doubles */
public interface ProbabilityDistribution2  extends Comparable, Serializable{
   // public double sample();

    public double probability(double x1, double x2);
  
    public void addCount(double obj_index,double obj_indexy, double value);

    public ProbabilityDistribution2 clone();
    public ProbabilityDistribution2 clone(double u);
    
  //  public void transfercounts(EmissionState innerState, int phen_index, int i);

    public void initialise();

  //  public void setParamsAsAverageOf(ProbabilityDistribution2[] tmp);

    public void transfer(double pseudoC);

    public String id();
    
   // public void addCounts(ProbabilityDistribution2 probabilityDistribution);

   // public double[] getCount(double[] angle);

	public int getParamIndex();

	public String name();

	//public double[] getMean();

	//public double sum();

	public void maximise(double d, double e, double f, double d1, double e1, double f1, double g);

	public void updateParamIndex();

	public void recalcName();

	//public double prior();

	//public double getParamValue(int n1);

//	public void setParamValue(int n1, double val);
	
	public void  getInterval(double[] input,DoubleMatrix2D res, double[] mean) ;

	public int numObs();

	public int fill(DoubleMatrix2D x, DoubleMatrix2D y,
			DoubleMatrix2D yB, int numObs,double[] noCop,  double pseudo);
	 public void variance( int type, double[] sum);
	//public double plotObservations(String string, boolean b, XYSeries obs, boolean swtch);

	//public void plotTheoretical(String string, boolean b, XYSeries theor);

	//public void getIntervalX(double[] in, double[] res_inR);
	//public void getIntervalY(double[] in, double[] res_inR);
	/*order [mean_x, mean_y][ sigma_x, sigma_y, covar][skew_x, skew_y] */
	public void setParam(int type, int i, double d);

	

	public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
			DoubleMatrix2D covar, int numObs, double pseudo);

	public void print(PrintWriter pw);

	public void setPriors(ProbabilityDistribution2 probabilityDistribution2, int type, boolean x
			);

	
	
	public ProbabilityDistribution2 clone(double u, SimpleExtendedDistribution1 simpleExtendedDistribution1);

	public void setMinMax(double minR, double maxR, double min, double max);

	public void setToExclude();

	public double probability(double r, double b, int mixComponent);

	public void addCount(Double r2,  Double b, double val,
			SimpleExtendedDistribution1 mixe1,
			ProbabilityDistribution2 probDistG);

	

	
	
	

}
