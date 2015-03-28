package lc1.stats;

public class Contigency {
	
	public static void main(String[] args){
		try{
			double[] a= new double[] {1.3,3.5,6.4};
			double[] b = new double[] {1.8,5.4, 8.7};
			double[][] m = new double[][] {a,b};
			Contigency c = new Contigency();
			c.setMatrix(m);
			double chisq = c.calcchiSquare();
			double p = c.getSig();
			System.err.println(chisq+" "+p);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
	ChiSq chisq = new ChiSq();
	int deg ;
	/**
	   * constructor for Contigency table
	   *
	   * @param maxSize is the maximum sum that will be encountered by contigency table
	   */
	  public Contigency() {
	
	  }
	  final double calcchiSquare()
	  {
	  	int i,j;
	  double chi=0, E;

	  	for(i=0; i<rows; i++)
	            {
	            for(j=0; j<cols; j++)
	                {if (contig[i][j]>0)
	                    {E=expectation[i][j];
	                    chi+=(((float)contig[i][j]-E)*((float)contig[i][j]-E)/E);
	                    }
	                }
	            }
//	  	printf("The chi=%8.4f\n",chi);
	  	return chi;
	  }
	  public double chisq(){
		  return calcchiSquare();
	  }
	  final public double getSig(){
		  double chi2 = calcchiSquare();
		  return this.chisq.chi2prob(deg, chi2);
	  }
double[][] contig;
double csum;
int rows, cols;
double[] crow,ccol;// rowDist, colDist;
double[][] expectation;
	  public void setMatrix(double[][] tcontig)
	  {  //this permutes of rowDist to rapidly do the permutations
	 // int i,j,k,count=0;
	
	  rows=tcontig.length;
	  cols=tcontig[0].length;
	 
	  contig=tcontig;
	  csum=0;
	  crow=new double[rows];
	  ccol=new double[cols];
	//  int no_non_zero_cols=0;
	 
		//  boolean nonZero = false;
	  for(int i=0; i<rows; i++)
		{
		  for(int j=0; j<cols; j++)
	      {
			//  nonZero=true;
			  csum+=contig[i][j];
	          crow[i]+=contig[i][j];
	          ccol[j]+=contig[i][j];
	     
		
		}
	 // if(nonZero) no_non_zero_cols++;
	   }
	  this.deg = (rows-1)*(cols-1);
	  this.expectation = new double[rows][cols];
	  for(int i=0; i<rows; i++)
		{for(int j=0; j<cols; j++)
			
	          {
			expectation[i][j] = ((double) ccol[j]/ (double) csum)*(double) crow[i];
			
	          }
		}
	  }
//	  rowDist=new int[csum];     //This sets up the row distribution so that the random number requires only one call
	//  colDist=new int[csum];
	  

}
