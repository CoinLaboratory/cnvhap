/**
 * 
 */
package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Collection;

import lc1.stats.PseudoDistribution;


public class FreeExpTransitionProbs  extends AbstractTransitionProbs implements Serializable{
    
	//MatrixExp exp;
	

	double invlen;
	//final double distance;
	
	//double[][]counts;
    @Override
    public int noStates() {
       return this.no_states;
    }
   public MatrixExp mat;
   
    public FreeExpTransitionProbs(double distance, double[] pi, double rate){
    	
    	mat = new MatrixExp( pi,rate);
    	mat.setDistance(distance);
		this.no_states = pi.length+1;
		//this.rates = rates;
		this.invlen = 1.0 / (double)pi.length;
	//	this.distance = distance;
	//	this.updateM(this.distance);
    //  counts = new double[pi.length][pi.length];
    }
    public void updateRates(double[][] mat, double[] pi) {
		this.mat.updateRates(mat, pi);
		
	}
    
  public FreeExpTransitionProbs(double dist, double[][]matr, double rate){
	  double[][]matr1  = new double[matr.length-1][];
	  System.arraycopy(matr, 1, matr1, 0, matr1.length);
    mat = new MatrixExp(matr1, matr[0][0]);
    mat.setDistance(dist);
	this.no_states = matr.length;
	//this.rates = rates;
	this.invlen = 1.0 / (double)matr.length;
	//this.distance = dist;
    }
   
   
    
 /*   protected void normalize(DoubleMatrix2D rates, double[] frequency)
    {
      double subst = 0.0;
int dimension = rates.size();
      for (int i = 0; i < dimension; i++)
      {
        subst += -rate[i][i]*frequency[i];
      }
      for (int i = 0; i < dimension; i++)
      {
        for (int j = 0; j < dimension; j++)
        {
          rate[i][j] = rate[i][j]/subst;
        }
      }
    }*/

    public FreeExpTransitionProbs(double[] modifyFrac0,
			double[] expModelIntHotSpot1, Double double1) {
		this(0,modifyFrac0, double1);
	}
  
  
    public FreeExpTransitionProbs(double[] rel1, double dist, int len,
			double[] expModelIntHotSpot1, Double double1) {
		this(dist,  modify(rel1), double1);//, modify(expModelIntHotSpot1),0.5), );
	}
    private static double[] modify(double[] rel1) {
		double[]res = new double[rel1.length-1];
		System.arraycopy(rel1, 1, res, 0, res.length);
		return res;
	}


	
    
  
    
   
	public FreeExpTransitionProbs(FreeExpTransitionProbs probs){
        no_states = probs.noStates();
     this.mat = probs.mat;
      this.invlen = probs.invlen;
    //  this.distance = probs.distance;
    //  this.counts = probs.counts;
        
    }
	
	public FreeExpTransitionProbs(FreeExpTransitionProbs probs, double dist){
        no_states = probs.noStates();
     this.mat = new MatrixExp(probs.mat);
      this.invlen = probs.invlen;
      mat.setDistance(dist);
    //  this.counts = probs.counts;
        
    }









final int no_states;
   
  public FreeExpTransitionProbs clone(boolean swtch){
     return new FreeExpTransitionProbs(this);
  }
   
   
    
    
   
    @Override
    public void addCount(int from, int to, double val){
    	throw new RuntimeException("!!");
    }
    @Override
    void addCount(int from, int to, double val, double dist) {
    	if(from>0 && to>0  && dist>1e-7){
    		//counts[from-1][to-1] +=val;
    		this.mat.addCount(from-1, to-1, val, dist);
    		
    	}
	}
    @Override
    public MatrixExp mat() {
		return mat;
	}
    
    public double evaluate(){
    	return 0;//this.mat.evaluate();
    }
    
    @Override
    public void addCount(int indexFrom,int from, int to, double val, double dist){
       throw new RuntimeException("!!");
     
    }
    
    @Override
    public double getTransitionToPaint(int from, int to) {
    	if(to==0) return 0;
    	else if(from==0) return invlen;
    	else{
    		double res =Math.log( Math.abs(mat.rates.getQuick(from-1, to-1)));
    		//if(Double.isNaN(res)){
    		//	throw new RuntimeException("!!");
    		//}
    		return res;
    	}
	}
    
    public double getTransition(int from, int to){
    	if(to==0) return 0;
    	else if(from==0) return invlen;
    	else{
    		double res = mat.M.getQuick(from-1, to-1);
    		//if(Double.isNaN(res)){
    		//	throw new RuntimeException("!!");
    		//}
    		/*if(res<0){
    			Logger.global.info("h");
    		}*/
    		return res;
    	}
        
    }
   
    
  
    
   
     public void initialiseCounts(boolean start, boolean end){
    	 mat.initialise();
        /*  for(int j=0; j<this.counts.length; j++){
                 Arrays.fill(counts[j],0);
                 
          }*/
     }
    
    
 
  
    
   
   
  
     
  
    public String toString(){
         return mat.M.toString();
     }
  //  final double[] pseudoCountTrans;
     /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#transfer(double)
   
    
    public void transfer(double[] pseudoC, double[] pseudo_exp){
        transfer(pseudoC[0], pseudo_exp[0]);
    }  */
    
     
    
  
    
   

    
  

    




    public int length() {
        return this.no_states;
    }



  

    public void print(PrintWriter pw, Double[] hittingProb, double dist) {
       
          
      }
	@Override
	public double transfer(double[] pseudoCExp, double[][] d, int pos_index) {
		return 0;
		//	transfer(pseudoCExp[pos_index]);
	//return this.evaluate();
	}

	
	
	
	
	public double  transferAlpha(double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
		 // double[][]v = new double[][] {alpha_overall};
		 return 0;
	  }


	@Override
	public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
		return (new FreeTransitionProbs1(this)).clone(statesToGroup, u);
//		throw new RuntimeException("!!");
	}


	@Override
	public double[] getAlphaPrior() {
		throw new RuntimeException("!!");
	}


	@Override
	public Collection getDistributions() {
		throw new RuntimeException("!!");
	}


	@Override
	public void validate() {
		mat.validate();
		
	}
	
@Override
	public double transferQ(double[] ds, double pseudoAlpha, double pseudoRate, MatrixExp initial, int i, double distance, int it) {
		this.mat.transfer(ds[i], pseudoAlpha, pseudoRate, initial, distance);
		return 0;//this.mat.evaluate(counts);
		
	}
public double getRate(){
	return this.mat.currentRate;
}

	public void addCounts(PseudoDistribution[] observed, double dist,
			int[] stateToGroup) {
		 int no_states = observed.length;
	        for(int j=0; j<no_states; j++){
	            PseudoDistribution dist1 =observed[j];
	            int st = j;
	            double[] counts = dist1.counts();
	            if(dist1==null) continue;
	            for(int j1 = 0; j1<no_states; j1++){
	                int state = j1;
	                Double val =counts[state];
	                if(val==0) continue;
	                addCount(stateToGroup[st], stateToGroup[state], val, dist);
	            }
	        }
		
	}
	
	
	



	/*public void transferQ(double[] ds, MatrixExp initial, int i,
			AbstractTransitionProbs[] transProbs, List<Double> loc) {
		for(int k=1; k<transProbs.length; k++){
			double dist = loc.get(k)-loc.get(k-1);
			this.mat.transfer(pseudo, counts, expected)
		}
		
	}*/
	
	  
	
    
    
  
   
   
       
      
}