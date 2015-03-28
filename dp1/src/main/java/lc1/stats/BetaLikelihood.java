package lc1.stats;

import lc1.util.Constants;
import cern.jet.stat.Gamma;

public class BetaLikelihood {
	
	public static void main(String[] args){
		try{
			
		double[] a = new double[] {52.3,  289,  303};
	double[] b = new double[] { 1.00,0.000318,0.000155};
		
		//double[] a = new double[] {645,	1.95e-05,	1.90e-16};
	//	double[] b = new double[] { 1.00,	1.13e-07,	4.71e-19};

			//double[] a= new double[] {1.3,3.5,6.4};
			//double[] b = new double[] {1.8,5.4, 8.7};
			double[][] m = new double[][] {a,b};
			BetaLikelihood c = new BetaLikelihood();
			c.setMatrix(m);
			double bf = c.bayesFactor();
			double p = c.getSig(bf,0.5);
			double sig = c.getSig();
			System.err.println(bf+" "+p+" "+sig);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
	ChiSq chisq = new ChiSq();
	//int deg ;
	/**
	   * constructor for Contigency table
	   *
	   * @param maxSize is the maximum sum that will be encountered by contigency table
	   */
	  public BetaLikelihood() {
	
	  }
	  
	  public double nullLogLikelihood(){
		  double a  =rsum[0]+1;
		  double b = rsum[1]+1;
		  double logg = logBeta(a,b);
		  return logg;
		  
	  }
	  public double alternativeLogLikelihood(){
		  double prob = 0;
		  for(int i=0; i<cols; i++){
			  double a  =this.contig[0][i]+1;
			  double b = this.contig[1][i]+1;
			double logg =logBeta(a,b);
			//  System.err.println(a+" "+b+" "+logg);
			  prob+=logg;
		  }
		  return prob;
	  }
	  public double logBeta(double a, double b){
		  return Gamma.logGamma(a)+Gamma.logGamma(b) - Gamma.logGamma(a+b);
	  }
	 
	  public double bayesFactor(){
		  double altLH  = this.alternativeLogLikelihood();
		  double nullLH = this.nullLogLikelihood();
		  double bf = Math.exp(altLH - nullLH);
		 return bf;
	  }
	  final public double getSig(double bf, double priorOdds){
		  double po = (priorOdds / (1-priorOdds))*bf;
		  double ppa = po/(1+po);
		 return 1.0-ppa; 
	  }
	  
	  public double getSig(){
		  double bf = this.bayesFactor();
		  return this.getSig(bf,Constants.priorOdds());
	  }
double[][] contig;
double[] rsum;
int rows, cols;
//double[] crow,ccol;// rowDist, colDist;
//double[][] expectation;
	  public void setMatrix(double[][] tcontig)
	  {  //this permutes of rowDist to rapidly do the permutations
	 // int i,j,k,count=0;
	
	     rows=tcontig.length;
	     if(rows!=2) throw new RuntimeException("!!");
	     cols=tcontig[0].length;
	     contig=tcontig;
	     this.rsum = new double[rows];
	     for(int i=0; i<rows; i++){
	    	 for(int j=0; j<cols; j++){
	    		 this.rsum[i] +=contig[i][j];
	    	 }
	     }
	
	  }
//	  rowDist=new int[csum];     //This sets up the row distribution so that the random number requires only one call
	//  colDist=new int[csum];
	  

}
