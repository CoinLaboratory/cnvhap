package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.states.EmissionState;

import org.apache.commons.math.util.MathUtils;
import org.jfree.data.xy.XYSeries;

import pal.math.OrthogonalHints;
import Jama.Matrix;
import cern.colt.matrix.DoubleMatrix2D;



public class TrainableWeightedPoissonDistribution implements ProbabilityDistribution, Serializable{

	public static final double ERROR=1E-5;
	double alpha, lambda, r;
	public double mean1;
	public double mode;
	
	List<Double> obsv = new ArrayList<Double>();
	List<Double> obsx = new ArrayList<Double>();
	double sum=0;
	
	public TrainableWeightedPoissonDistribution(double alpha, double lambda, double r, double mean, double mode)
	{
		this.alpha=alpha;
		this.lambda=lambda;
		this.r=r;
		this.mean1=mean;
		this.mode=mode;
	}
	
	public TrainableWeightedPoissonDistribution(double mean, double mode)
	{
		this.mean1 = mean;
		this.mode = mode;		
	}
	
	public TrainableWeightedPoissonDistribution(String string, double meanIR,
			double varIR, double round, double d) {
		this.lambda = meanIR;
		this.alpha = 1e-5;
		this.r = 0;
		this.mean1 = meanIR;
	
		this.mode = mean1;
		
	}

	public double limitOfC(double lamda, double dispersion, double displacement)//, double error)
	{
		double currError = Double.POSITIVE_INFINITY;
		double sum=0;
		int k=0;
			while (currError>ERROR)
			{
				currError=sum;
				sum+=Math.pow(lamda, k)*Math.pow((k+displacement),dispersion)/MathUtils.factorialDouble(k);

				currError=Math.abs(currError-sum)/sum;

				k++;
			}
		return sum;
	
	}
	
	public double limitOfD(double lamda, double dispersion, double displacement, int s)
	{
		double currError = Double.POSITIVE_INFINITY;
		double sum=0;
		int k=0;
		while (currError>ERROR)//&& k<=12800)// || k<=3*Math.round(lamda))
		{
			currError=sum;
			sum+=Math.pow(lamda, k)*Math.pow((k+displacement),dispersion)*Math.pow(Math.log(k+displacement),s)/MathUtils.factorialDouble(k);
			currError=Math.abs(currError-sum);///sum;
			k++;
		}
//		System.out.println(k+"\t"+sum);
		return sum;
	
	}
	
	
	public double getTau1(double lamda, double r, double a, double error)
	{
		return lamda*limitOfC(lamda, r, a+1)/limitOfC(lamda, r, a);
	}
	

	
	public double getTau2(double lamda, double r, double a, double error)
	{
		return limitOfD(lamda, r, a, 1)/limitOfC(lamda, r, a);
	}
	
	public double[][] getInverseFisherInformationMatrix(double lamda, double r, double a, double error)
	{
		double[][] fisher = new double[2][2];
//		double[] tempCLimit = new double[]{limitOfC(lamda,r, a+2,error),limitOfC(lamda,r,a,error),limitOfC(lamda,r,a+1,error)};
//		double[] tempDLimit = new double[]{limitOfD(lamda,r,a+1,1,error)};


		fisher[0][0]=((Math.pow(lamda,2)*(limitOfC(lamda,r, a+2)*limitOfC(lamda,r,a)-Math.pow(limitOfC(lamda,r,a+1),2)))+lamda*limitOfC(lamda,r,a+1)*limitOfC(lamda,r,a))/Math.pow(limitOfC(lamda,r,a),2);
		fisher[0][1]=lamda*(limitOfC(lamda,r,a)*limitOfD(lamda,r,a+1,1)-limitOfC(lamda,r,a+1)*limitOfD(lamda,r,a,1))/Math.pow(limitOfC(lamda,r,a),2);
		fisher[1][0]=fisher[0][1];
		fisher[1][1]=(limitOfC(lamda,r,a)*limitOfD(lamda,r,a,2)-Math.pow(limitOfD(lamda,r,a,1),2))/Math.pow(limitOfC(lamda,r,a),2);			

		

		

		
//		System.out.println("fisher matrix");
//		System.out.println(fisher[0][0]+"\t"+fisher[0][1]);
//		System.out.println(fisher[1][0]+"\t"+fisher[1][1]);
		
		Matrix matrixFisher = new Matrix(fisher);
		if (matrixFisher.det()==0)
			return new double[][]{{Double.NaN,Double.NaN},{Double.NaN,Double.NaN}};
		else
		{
			double[][] inverse= matrixFisher.inverse().getArray();
			return inverse;
		}
//		return inverseFisher;
	}
	
	public double[][] getFisherInformationMatrix(double lamda, double r, double a, double error)
	{
		double[][] fisher = new double[2][2];
		
		fisher[0][0]=((Math.pow(lamda,2)*(limitOfC(lamda,r, a+2)*limitOfC(lamda,r,a)-Math.pow(limitOfC(lamda,r,a+1),2)))+lamda*limitOfC(lamda,r,a+1)*limitOfC(lamda,r,a))/Math.pow(limitOfC(lamda,r,a),2);
		fisher[0][1]=lamda*(limitOfC(lamda,r,a)*limitOfD(lamda,r,a+1,1)-limitOfC(lamda,r,a+1)*limitOfD(lamda,r,a,1))/Math.pow(limitOfC(lamda,r,a),2);
		fisher[1][0]=fisher[0][1];
		fisher[1][1]=(limitOfC(lamda,r,a)*limitOfD(lamda,r,a,2)-Math.pow(limitOfD(lamda,r,a,1),2))/Math.pow(limitOfC(lamda,r,a),2);			

		return fisher;
	}
	
	public double getSampleMean(Double[] sample, Double[] probs, double pseudocount)
	{
		double mean=this.mean1*pseudocount;
		double cnt=pseudocount;
		for (int i=0;i<sample.length;i++){
			mean+=probs[i]*sample[i];
			cnt+=probs[i];
		}
		mean/=cnt;
		return mean;
	}
	
	public double getSampleVariance(Double[] sample,double mean)
	{
		double variance=0;
		for (int i=0;i<sample.length;i++)
			variance+=Math.pow((sample[i]-mean),2);
		
		variance/=sample.length;
		return variance;
	}
	public double getLogGeometricShiftedSampleMean(Double[] sample, double a, Double[] weights, double pseudocount)
	{
		double mean=pseudocount*Math.log(this.mean1+a);
		double cnt=pseudocount;
		for (int i=0;i<sample.length;i++){
			mean+=weights[i]*Math.log(sample[i]+a);
			cnt+=weights[i];
		}
		mean/=cnt;
		return mean;
	}
	
	public double empiricalProbability(int k, Double[] sample)
	{
		int count=0;
		for (int i=0;i<sample.length;i++)
		{
			if (sample[i]==k)
				count++;
//			System.out.println(k);
		}
				
//		System.out.println((double) count/sample.length);
		return (double) count/sample.length;
	}
	boolean flagFirst = true;
	public double[] getInitialGuesses(double a, Double[] sample)
	{
		if(flagFirst){
			flagFirst = false;
		
		int k=(int) Math.round(this.mode);
//		double r = Math.log10(probability(k,sample)*probability(k+2,sample)*(k+2)*Math.pow((k+a+1),2))-Math.log10((Math.pow(probability(k+1,sample),2)*(k+1)*(k+a)*(k+a+2)));
		double r = Math.log10(empiricalProbability(k,sample)*empiricalProbability(k+2,sample)*(k+2)/(Math.pow(empiricalProbability(k+1,sample), 2)*(k+1)))/Math.log10((k+a)*(k+a+2)/Math.pow((k+a+1),2));
//		System.out.println(probability(k,sample)*Math.pow((1+1/(double) (k+a)),r));
		double lamda = empiricalProbability(k+1,sample)*(k+1)/(double) (empiricalProbability(k,sample)*Math.pow((1+1/(double) (k+a)),r));
//		System.out.println(lamda+"\t"+r);
		return new double[]{lamda,r};	
		}
		else return new double[] {this.lambda, this.r};
	}
	
	public double[] estimateParameters(double a,Double[] sample,Double[] weights, double pseudocount)
	{
		 double parameterError = Double.POSITIVE_INFINITY;
		
		double t1 = getSampleMean(sample,weights,pseudocount);
		double t2 = getLogGeometricShiftedSampleMean(sample,a,weights,pseudocount);
		double[] parameters = getInitialGuesses(a,sample);
		int iterations=0;
		while (parameterError>=ERROR && iterations<1000)
		{
			parameterError = parameters[0];
			double[][] inverseFisherMatrix = getInverseFisherInformationMatrix(parameters[0], parameters[1], a, ERROR);
			double tau1 = getTau1(parameters[0], parameters[1], a, ERROR);
			double tau2 = getTau2(parameters[0], parameters[1], a, ERROR);
			parameters[0]=parameters[0]+inverseFisherMatrix[0][0]*(t1-tau1)+inverseFisherMatrix[0][1]*(t2-tau2);
			parameters[1]=parameters[1]+inverseFisherMatrix[1][0]*(t1-tau1)+inverseFisherMatrix[1][1]*(t2-tau2);
			parameterError = Math.abs(parameterError-parameters[0]);
			iterations++;
		}
		return parameters;
	}
	
	public double[][] parameterProfile(Double[] sample, Double[] weight, double maxAlpha, double step,double pseudo) {
		double t1 = getSampleMean(sample,weight,pseudo);
		
		double[][] profile = new double[((int) Math.round((maxAlpha-0.01)/step))+1][4];
		
		for (double[] row : profile){
			Arrays.fill(row, Double.NaN);
		}
		int i=0;
//		System.out.println("mean: "+t1);
		for (double alpha=0.01d;alpha<maxAlpha;alpha+=step)
		{
//			System.out.println(getInitialGuesses(i, sample)[1]);

			double t2 = getLogGeometricShiftedSampleMean(sample,alpha,weight,pseudo);

			try
			{
				double[] param = estimateParameters(alpha, sample,weight,pseudo);
				profile[i][0] = logL(param,alpha,sample.length,t1,t2);
				profile[i][1] = alpha;
				profile[i][2] = param[0];
				profile[i++][3] = param[1];
//				System.out.println(logL(param,i,0.0001,sample.length,t1,t2)+"\t"+alpha+"\t"+param[0]+"\t"+param[1]);
			}catch(Exception e)
			{
				System.err.println(e.getMessage());
				i++;
			}
		}
		return profile;
	}
	
	public double[][] parameterProfile(double minAlpha, double maxAlpha, double step,double pseudo) {
		Double[] sample = new Double[this.obsx.size()];
		sample = this.obsx.toArray(sample);
		Double[] weight = new Double[this.obsx.size()];
		weight = this.obsv.toArray(weight);
		double t1 = getSampleMean(sample,weight,pseudo);
		
		double[][] profile = new double[((int) Math.round((maxAlpha-minAlpha)/step))+1][4];
		
		for (double[] row : profile){
			Arrays.fill(row, Double.NaN);
		}
		int i=0;
//		System.out.println("mean: "+t1);
		for (double alpha=minAlpha;alpha<maxAlpha;alpha+=step)
		{
//			System.out.println(getInitialGuesses(i, sample)[1]);

			double t2 = getLogGeometricShiftedSampleMean(sample,alpha,weight,pseudo);

			try
			{
				double[] param = estimateParameters(alpha, sample,weight,pseudo);
				profile[i][0] = logL(param,alpha,sample.length,t1,t2);
				profile[i][1] = alpha;
				profile[i][2] = param[0];
				profile[i++][3] = param[1];
//				System.out.println(logL(param,i,0.0001,sample.length,t1,t2)+"\t"+alpha+"\t"+param[0]+"\t"+param[1]);
			}catch(Exception e)
			{
				System.err.println(e.getMessage());
				i++;
			}
		}
		return profile;
	}
	
	public double logL(double[] params, double a, int sampleSize, double mean, double logMean)
	{
		double theta = Math.log(params[0]);
		double rt2 = params[1]*logMean;
		double constant = limitOfC(params[0],params[1],a);
//		System.out.println("params0: "+params[0]+"\tparams1: "+params[1]+"\ttheta: "+theta+"\tlogmean: "+logMean+"\trt2: "+rt2+"\tconstant: "+constant);
		return sampleSize*(theta*mean+rt2-Math.log(constant));
	}
	
	
	
	@Override
	public double probability(double k)
	{
		double pdf = Math.exp(k*Math.log(lambda)+r*Math.log(k+alpha))/(Math.exp(GammaFunction.lnGamma(k+1))*limitOfC(lambda,r,alpha));
		//System.out.println(pdf);
		if(Double.isNaN(pdf)){
			System.err.println(pdf);
		}
		return pdf;
	}
	
	public double probability(int k){
		double pdf = Math.exp(k*Math.log(lambda)+r*Math.log(k+alpha))/(MathUtils.factorialDouble(k)*limitOfC(lambda,r, alpha));
//		System.out.println(pdf);
		return pdf;
	}
	
	
	
	
	@Override
	public void addCount(double objIndex, double value) {
		this.obsx.add(objIndex);
		this.obsv.add(value);
		
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

	public ProbabilityDistribution clone() {
		return new TrainableWeightedPoissonDistribution(this.alpha,this.lambda,this.r, this.mean1, this.mode);
	}
	@Override
	public ProbabilityDistribution clone(double u) {
		return new TrainableWeightedPoissonDistribution(this.alpha,this.lambda,this.r,this.mean1, this.mode);
	}

	@Override
	public ProbabilityDistribution clone(double u,
			SimpleExtendedDistribution1 dist1) {
		return new TrainableWeightedPoissonDistribution(this.alpha,this.lambda,this.r,this.mean1, this.mode);
	}

	@Override
	public int fill(DoubleMatrix2D x, DoubleMatrix2D y, int numObs,
			double[] noCop, double pseudo) {
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
	public double getMean() {
		// TODO Auto-generated method stub
		return 0;
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
		obsx.clear();obsv.clear();
	    sum =0;
		
	}

	@Override
	public void maximise(double pseudo, double pseudo_var, double pseudo_skew) {
		this.maximise(0.0);
		
	}

	public void maximise(double pseudo){
		double minAlpha = 1E-5;
//		double maxAlpha = mean1;
		double maxAlpha = 1.0;
		double step = 0.05;
		double[][] profile = parameterProfile(minAlpha, maxAlpha, step, pseudo);
		double maxLogL = Double.NEGATIVE_INFINITY;
		int maxLogLindex = -1;
		double bestLambdaDiff = Double.POSITIVE_INFINITY;
		int bestLambdaIndex = -1;
		
			
		for (int i=0; i<profile.length;i++)
		{
			if (profile[i][0]>maxLogL)// && profile[i][0]<0)
			{
				maxLogL = profile[i][0];
				maxLogLindex = i;
			}
			if (Math.abs(this.mean1-profile[i][2])<bestLambdaDiff)
			{
				bestLambdaDiff = Math.abs(this.mean1-profile[i][2]);
				bestLambdaIndex = i;
			}
			
		}
		if (maxLogLindex==-1 && bestLambdaIndex==-1)
		{
			System.err.println("Not well behaved Weighted Poisson. Reverting to simple Poisson");
			this.alpha = 1e-5;
			this.lambda = mean1;
			this.r = 0;
			return;
		}
		if (maxLogLindex!=bestLambdaIndex)
		{
			if (bestLambdaIndex!=-1)
				if (!Double.isNaN(profile[bestLambdaIndex][0]))// && (profile[bestLambdaIndex][3]>=-1))
//					if (profile[bestLambdaIndex][0]<0)
					{
						System.err.println("Chosen best lambda for Weighted Poisson");
						this.alpha = profile[bestLambdaIndex][1];
						this.lambda = profile[bestLambdaIndex][2];
						this.r = profile[bestLambdaIndex][3];
						return;
					}
			
		}
		if (maxLogLindex!=-1)
			if (!Double.isNaN(profile[maxLogLindex][0]))// && (profile[maxLogLindex][3]>=-1))
			{
				System.err.println("Chosen max logL for Weighted Poisson");
				this.alpha = profile[maxLogLindex][1];
				this.lambda = profile[maxLogLindex][2];
				this.r = profile[maxLogLindex][3];
				return;
			}
		System.err.println("Not well behaved... Reverting to simple Poisson");
		this.alpha = 1e-5;
		this.lambda = mean1;
		this.r = 0;
		
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
	public void setCoverage(double d) {
		// TODO Auto-generated method stub
		
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
	public int compareTo(Object arg0) {
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
