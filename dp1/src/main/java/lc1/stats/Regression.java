package lc1.stats;

import java.util.List;
import java.util.concurrent.ExecutorService;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.Constants;
import lc1.util.Executor;
import pal.math.UnivariateMinimum;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.random.Binomial;

public abstract class  Regression {
   
    static pal.statistics.NormalDistribution normal;
    DataCollection dc;
    final Binomial binomial = new Binomial(10,0.5, null);
    // List<PseudoDistribution> data = new ArrayList<PseudoDistribution>();
    final PseudoDistribution []covariate;// = new PseudoDistribution[3];
   

    // double stddev;
    int position = -1;
    int phen_index = -1;
    int[][] vals;
    int max_cn;

    public Regression(DataCollection dc) {
        this.dc = dc;
        indiv = dc.indiv();
        EmissionStateSpace emStSp =Emiss.getSpaceForNoCopies(Constants.backgroundCount(dc.index));// dc.getEmStSpace();
        List<Integer> l = emStSp.copyNumber;
        max_cn = l.get(l.size() - 1) + 1;
        List<Integer> cn = emStSp.copyNumber;
        vals = new int[3][emStSp.size()];
        for (int i = 0; i < emStSp.size(); i++) {
            ComparableArray comp = (ComparableArray) emStSp.get(i);
            int noCop = comp.noCopies(true);
            int noB = (int) comp.noB();
            int noA = noCop - noB;

            vals[0][i] = noCop;

            vals[1][i] = noA;
            vals[2][i] = noB;
        }
        
        this.covariate = new PseudoDistribution[indiv.size()];
       
        int[] ycols = new int[] {0};
       
       Y = new DoubleMatrix2D[dc.pheno.phen.size()];
      Z = new DoubleMatrix2D[dc.pheno.phen.size()];
       this.rows = new int[Y.length][];
         for(int i=0; i<Y.length; i++){
         	Z[i] = new DenseDoubleMatrix2D(dc.indiv().size(),1);
         }
      
         for (int i=0; i<indiv.size(); i++) {
        	 String key = indiv.get(i);
 		    HaplotypeEmissionState hes = (HaplotypeEmissionState) dc.dataL.get(key);
 		    Double[] phenV = hes.phenValue();
 		   for (int phen_index=0; phen_index<Y.length; phen_index++) {
 		        Double ph = phenV[phen_index];
 		        if (ph != null && !Double.isNaN(ph)) {
 		            Z[phen_index].setQuick(i, 0, ph);
 		        }
 		        else Z[phen_index].setQuick(i, 0, Double.NaN);
 		    }
 		}
         for(int i=0; i<Y.length; i++){
        	 rows[i] = nonNA(Z[i]);
          	Y[i] = Z[i].viewSelection(rows[i], ycols);
          }
         	logL1 = new double[Y.length];
         	this.fillNull();
        
       // DoubleMatrix2D XT =  lg.transpose(X);
    }
     int[] nonNA(DoubleMatrix2D matrix2D) {
	int num= 0;
	for(int i=0; i<matrix2D.rows(); i++){
		if(!Double.isNaN(matrix2D.getQuick(i, 0))) num++;
	}
	int[] res = new int[num];
	 num=0;
	for(int i=0; i<matrix2D.rows(); i++){
		if(!Double.isNaN(matrix2D.getQuick(i, 0))){
			res[num] =i;
			num++;
		}
		
	}
	return res;
	}
	final List<String> indiv;
    public void setGeno(int position){
    	this.position = position;
    	 
    	 for (int i=0; i<indiv.size(); i++) {
        	 String key = indiv.get(i);
 		    HaplotypeEmissionState hes = (HaplotypeEmissionState) dc.dataL.get(key);
 		    this.setGeno(position, i, hes);
 		
 		}
    	//}
    }
    	

    void setGeno(int position2, int i, HaplotypeEmissionState hes) {
    	covariate[i]= hes.emissions[position2];
	}
	public static double exp(PseudoDistribution dist, int[] vals){
 	   Integer i = dist.fixedInteger();
 	   if(i!=null) return vals[i];
 	   else{
 		   double d = 0;
 		   double[] d1 = dist.probs();
 		   for(int k=0; k<d1.length; k++){
 			   d+=vals[k]*d1[k];
 		   }
 		   return d;
 	   }
    }
final DoubleMatrix2D[] Y,Z; //


int[][] rows; //indexed for each Y
  public abstract double calcNull(int phen_index);
 
    final double[] logL1;
    public void fillNull(){
    	 
    	for(int phen_index=0; phen_index < this.Y.length; phen_index++){
    		logL1[phen_index] = calcNull(phen_index);
    	
    	}
    }
    final UnivariateMinimum os1 = new UnivariateMinimum();
   
  
   
	public double calcLogLDiff(int pos_index, int phen_index, int type){
	
    	if(phen_index==0)this.setGeno(pos_index);
        this.phen_index = phen_index;
        return 0.0;
    }

   // double mean;
 //double variance;
   


    


	public double calcSignificance(double logLDiff) {
	    double res = ChiSq.chi2prob(1, 2 * logLDiff);
	    return res;
	}

	  public static ExecutorService es =   Executor.getEs(Regression.class, Constants.numThreads);
public double slope=0;
	public double slope() {
		// TODO Auto-generated method stub
		return slope;
	};
	
	
	
	
}