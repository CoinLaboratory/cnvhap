package lc1.stats;

import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.State;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class SimpleExtendedDistribution extends PseudoDistribution implements Serializable {
     public double[] probs;
   // public double[] pseudo;
    public  double[]  counts;
  //  public boolean fixed=false;
  //  public void setFixed(boolean f){
  //      this.fixed = f;
  //  }
   
    public void multiplyCounts(double d) {
    	if(true) throw new RuntimeException("!!");
		for(int i=0; i<counts.length; i++){
			counts[i] = d*counts[i];
		}
		
	};
    public void getB(double[] b, double mult) {
    	//double b0=0;
    	
		for(int i=1; i<probs.length; i++){
		  if(probs[i]!=0){
			  double cnt = (double) this.emstsp.getBCount(i);
			  b[0] +=mult*probs[i]* emstsp.getBCount(i);
			  b[1] +=mult*probs[i]*emstsp.getCN(i);
			
		  }
		}
		
	};
   public  EmissionStateSpace emstsp;
    public  double[] calcDistribution(
            double[] distribution, EmissionStateSpace emStSp, int pos) {
       return this.probs;
     }
    public boolean isMeasured(){
    	if (probs[Constants.getMax(this.probs)]>Constants.imputedThresh(0)) return true;
    	return false;
    }
    
    
    public double getUncertainty() {
		return probs[Constants.getMax(probs)];
	}
   
    
    public void setProb(double[] prob){
        System.arraycopy(prob, 0, probs, 0, prob.length);
      
    }
    
    @Override
    public double totalCount() {
     return Constants.sum(counts);
    }
    
    public final void transfer(double[] pseudo){
        double sum=0;
        for(int i=0; i<probs.length; i++){
            double total =  counts[i]+ pseudo[i];
            probs[i] = total;
           // if(total<0){
            //	throw new RuntimeException("!!");
           // }
            sum+=total;
        }
        if(sum>1e-9){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
            }
        }
        if(Constants.CHECK){
            try{
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
        // for(int i=0; i<probs.length; i++){   
           if( Double.isNaN(probs[0]) || probs[0] <0) {
        	//   throw new RuntimeException("!!" );
           }
         //}	
        }
       // else{
        //    Arrays.fill(probs, val)
      //  }
    }
    public final void transfer(double[] pseudo, double mod){
        double sum=0;
        for(int i=0; i<probs.length; i++){
            double total =  counts[i]+ pseudo[i]*mod;
            probs[i] = total;
           // if(total<0){
            //	throw new RuntimeException("!!");
           // }
            sum+=total;
        }
        if(sum>1e-9){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
            }
        }
        else{
        	Arrays.fill(probs, 1.0/(probs.length));
        }
        if(Constants.CHECK){
            try{
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
        // for(int i=0; i<probs.length; i++){   
           if( Double.isNaN(probs[0]) || probs[0] <0) {
        	//   throw new RuntimeException("!!" );
           }
         //}	
        }
       // else{
        //    Arrays.fill(probs, val)
      //  }
    }
    
    public final void transfer(DoubleMatrix1D viewRow, double mod){
        double sum=0;
        for(int i=0; i<probs.length; i++){
            double total =  
            	i> 0 ? 
            	(counts[i]+ viewRow.get(i-1)*mod ): counts[i];
            probs[i] = Math.max(0,total);
           
            sum+=total;
        }
//        if(sum>1e-9){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
               
            }
  //      }
        if(Constants.CHECK){
            try{
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
        // for(int i=0; i<probs.length; i++){   
           if( Double.isNaN(probs[0]) || probs[0] <0) {
        	//   throw new RuntimeException("!!" );
           }
         //}	
        }
       // else{
        //    Arrays.fill(probs, val)
      //  }
    }
    public final void transfer(DoubleMatrix1D viewRow){
       // double sum=0;
        for(int i=0; i<probs.length; i++){
            double total =  
            	i> 0 ? 
            	 viewRow.get(i-1) :0;
            probs[i] = total;
          //  sum+=total;
        }
       Constants.normalise(probs);
        if(Constants.CHECK){
            try{
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
        // for(int i=0; i<probs.length; i++){   
           if( Double.isNaN(probs[0]) || probs[0] <0) {
        	//   throw new RuntimeException("!!" );
           }
         //}	
        }
       // else{
        //    Arrays.fill(probs, val)
      //  }
    }
    
    public void transferWithoutCounts(double[] pseudo){
   	 double sum=0;
        for(int i=0; i<probs.length; i++){
            double total =  pseudo[i];
            probs[i] = total;
            sum+=total;
        }
        if(sum>1e-9){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
            }
        }
        if(Constants.CHECK){
            try{
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
            }
    
 
       // else{
        //    Arrays.fill(probs, val)
      //  }
   }
    
    public final double evaluate(double[] pseudocountWeight){
        //  if(true) throw new RuntimeException("!!");
           //if(Constants.nested() ){
        	   this.transfer(pseudocountWeight);
           //}
        	 //  else{
        		//   this.transferWithoutCounts(pseudocountWeight);
        	  // }
            return logProb();
        }
    
    public double evaluate(DoubleMatrix1D viewRow, double mod) {
		this.transfer(viewRow, mod);
		return logProb();
	}
    
    public double evaluate(DoubleMatrix1D viewRow) {
	
		return logProb(viewRow);
	}
    
   
	public final double evaluate(double[] pseudocountWeight, double mod){
        //  if(true) throw new RuntimeException("!!");
           //if(Constants.nested() ){
        	   this.transfer(pseudocountWeight, mod);
           //}
        	 //  else{
        		//   this.transferWithoutCounts(pseudocountWeight);
        	  // }
            return logProb();
        }
    
    public final double evaluate(double pseudocountWeight){
        //  if(true) throw new RuntimeException("!!");
           //if(Constants.nested() ){
        	   this.transfer(pseudocountWeight);
           //}
        	 //  else{
        		//   this.transferWithoutCounts(pseudocountWeight);
        	  // }
            return logProb();
        }

    public double sample(){
        double r = Constants.rand.nextDouble();
        double cum=0;
        for(int i=0; i<probs.length; i++){
            cum+=probs[i];
            if(r < cum) return i;
        }
        throw new RuntimeException("!!");
    }
    
    /** convention is that order is [true, false] */
    /* prob is left undefined */
    public SimpleExtendedDistribution(int len){
        this.probs = new double[len];
        this.counts = new double[len];
      //  this.pseudo = new double[len];
        Arrays.fill(counts, 0);
       
    }
    
   public SimpleExtendedDistribution(int len, double u){
       this(len);
       double[] pseudo = new double[len];
       Arrays.fill(pseudo, 1.0/(double)len);
      Dirichlet dir = new Dirichlet(pseudo, u);
      arraycopy(dir.sample(), 0, probs, 0, len);
    }
   
   
    
   public static void normalise(double[] probs){
       double sum=0;
       for(int i=0; i<probs.length; i++){
         //  sum_ps+=pseudo[i];
           sum+=probs[i];
       }
       if(sum==0 && probs.length>0 ) {
    	   Arrays.fill(probs, 1.0/(double)probs.length);
//    	   throw new RuntimeException("!!");
       }
       else{
       for(int i=0; i<probs.length; i++){
         //  pseudo[i] = pseudo[i]/sum_ps;
           probs[i] = probs[i]/sum;
           if(Constants.CHECK && Double.isInfinite(probs[i])){
        	   throw new RuntimeException("");
           }
       }
       }
       //System.err.println("done");
   }
   public static void normalise(Double[] probs){
       double sum=0;
       for(int i=0; i<probs.length; i++){
         //  sum_ps+=pseudo[i];
           sum+=probs[i];
       }
       if(sum==0 && probs.length>0 ) throw new RuntimeException("!!");
       for(int i=0; i<probs.length; i++){
         //  pseudo[i] = pseudo[i]/sum_ps;
           probs[i] = probs[i]/sum;
       }
   }
   public SimpleExtendedDistribution(double[] init, double u){
        this( init.length);
       // if(init[Constants.getMax(init)] >0.999)throw new RuntimeException("!!");
      int len = init.length;
    //    System.arraycopy(init, 0, pseudo, 0, len);
        if(u==Double.POSITIVE_INFINITY){
          //  throw new RuntimeException("!!");
            System.arraycopy(init, 0, probs, 0, len);
        }
        else{
            Dirichlet dir = new Dirichlet(init, u);
            arraycopy(dir.sample(), 0, probs, 0, len);
        }
      }
   
   public SimpleExtendedDistribution(Double[] init){
       this( init.length);
      // if(init[Constants.getMax(init)] >0.999)throw new RuntimeException("!!");
     int len = init.length;
   //    System.arraycopy(init, 0, pseudo, 0, len);
      
         
           arraycopy(init, 0, probs, 0, len);
       
     }

 
   public double logProb(){
       double logL = 0;
       for(int j=0; j<counts.length; j++){
           if(counts[j]!=0)
              logL+=counts[j]*Math.log(probs[j]);
       }
      if(Double.isNaN(logL)){
    	  return Double.NEGATIVE_INFINITY;
      }
       return logL;
   }
   public double logProb(DoubleMatrix1D viewRow){
       double logL = 0;
       for(int j=1; j<counts.length; j++){
           if(counts[j]!=0)
              logL+=counts[j]*Math.log(viewRow.getQuick(j-1));
       }
      if(Double.isNaN(logL) ){
    	 return Double.NEGATIVE_INFINITY;
      }
       return logL;
   }
    
    
    
    public double prob(){
    	double res = 1.0;
    	for(int i=0; i<probs.length; i++){
    		res*=Math.pow(probs[i], counts[i]);
    	}
    	return res;
    }
   /* @Override
    public double evaluate(double[] pseudocountWeight){
        System.arraycopy(pseudocountWeight, 0,probs,0, probs.length);
      //  this.transfer(pseudocountWeight);
        return logProb();
    }*/
    
    public SimpleExtendedDistribution(PseudoDistribution dist_init) {
        this(dist_init.probs().length);
        int len = dist_init.probs().length;
    //    System.arraycopy(dist_init.pseudo, 0, pseudo, 0, len);
        this.setProb(dist_init.probs());
    }

  /*  public SimpleExtendedDistribution(SimpleExtendedDistribution dist_init) {
        this(dist_init.probs.length);
        int len = dist_init.probs.length;
      //  System.arraycopy(dist_init.pseudo, 0, pseudo, 0, len);
        System.arraycopy(dist_init.probs, 0, probs, 0, len);
    }*/
    
    @Override
    public void change(Sampler dir) {
    	Double[] sample=  dir.sample();
        arraycopy(sample, 0, probs, 0, probs.length);
		
	}
    
    public SimpleExtendedDistribution(Sampler dir) {
        this(dir.dist.length);
      //  arraycopy(dir.dist, 0, pseudo, 0, pseudo.length);
        Double[] sample=  dir.sample();
        arraycopy(sample, 0, probs, 0, probs.length);
    }

    public SimpleExtendedDistribution(double[] dist, double u,
			CompoundEmissionStateSpace stsp) {
		this(dist, u);
		this.emstsp=  stsp;
	}
	public SimpleExtendedDistribution(int len, EmissionStateSpace emStSp2) {
		this(len);
		this.emstsp = emStSp2;
	}
	public SimpleExtendedDistribution(double[] probs2, double switchU,
			EmissionStateSpace emStSp2) {
		this(probs2, switchU);
		this.emstsp = emStSp2;
	}
	public void initialise() {
       Arrays.fill(counts,  0);
    }

   /* public void setRandom(double emiss, boolean restart) {
        Dirichlet d = new Dirichlet(
                restart ? 
                      pseudo : 
                          probs, emiss);
        arraycopy(d.sample(), 0, probs, 0, probs.length);
        
    }*/
   protected void arraycopy(Double[] doubles, int i, double[] probs2, int j, int length) {
      for(int k=0; k<length; k++){
          probs2[k+j] = doubles[k+i];
      }
        
    }

  
    
    public void transfer(double pseudo){
        double sum=0;
       /* if(this.pseudo==null){*/
            for(int i=0; i<probs.length; i++){
                double total =  counts[i]+ (probs[i]==0 ? 0 : pseudo);
                probs[i] = total;
                sum+=total;
            }
      /*  }
        else{
            for(int i=0; i<probs.length; i++){
                double total =  counts[i]+ pseudo*(this.pseudo[i]);
                probs[i] = total;
                sum+=total;
            }
        }*/
        if(sum>0){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
            }
        }
        else{
            Arrays.fill(probs, 1.0/(double)probs.length);
        }
        if(Constants.CHECK){
        try{
        this.validate();
        }catch(Exception exc){
            exc.printStackTrace();
        }
        }
        if(Constants.CHECK && probs.length>2 && Double.isNaN(probs[2])) {
        	throw new RuntimeException("!!");
        }
    }
    
    public double KLDistance(PseudoDistribution d2){
        double sum =0;
        for(int j=0; j<probs.length; j++){
            double num1 = probs[j];
            double num2 = d2.probs(j);
            if(num1!=0){
                sum+=num1* Math.log(num1 / num2); 
            }
        }
        return sum;
    }

    public double sum() {
        double sum=0;
        for(int i=0; i<probs.length; i++){
            sum+=probs[i];
        }
        return sum;
    }

  /*  public void revertToPseudo() {
        System.arraycopy(this.pseudo,0, this.probs,0, probs.length);
        
    }

   public void resample(double u, boolean restart) {
        Dirichlet dir = new Dirichlet(restart ? pseudo : probs, u);
        arraycopy(dir.sample(), 0, probs, 0, probs.length); 
     }
    */
    public void addCounts(StateDistribution distribution) {
       double[] dist = distribution.dist;
        for(int j=0; j<dist.length; j++){
            this.counts[j]+=dist[j];
        }
    }
    public String toString(){
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        print(pw, false, getPrintString(), "\n");
        pw.close();
        return sw.getBuffer().toString();
    }
    
    public String getPrintString(){
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<probs.length; i++){
            sb.append("%5.3g ");
        }
        return sb.toString();
    }
    
    
    public void printSimple(PrintWriter pw, String name, String newLine, double thresh){
        pw.print(name+"->{");
        int[] l = Constants.getMax2(probs);
      //  for(int i=0; i<this.probs.length; i++){
        //    if(probs[i] > thresh){
        for(int k=0; k<l.length; k++){
                pw.print(l[k]+":"+Math.round(100*(probs[l[k]]))+",");
            }
      //  }
        pw.print("}"+newLine);
    }
    public void print(PrintWriter pw, boolean probsOnly, String st, String newLine){
       
        Double[] d= new Double[probs.length];
        for(int i=0; i<probs.length; i++){
            d[i] = probs[i];
        }
        pw.print(String.format(st, d)+";   ");
          if(!probsOnly){  
                  for(int i=0; i<probs.length; i++){
                   d[i] = counts[i];
                }
                pw.print(String.format(st, d)+";   ");
              /*  for(int i=0; i<pseudo.length; i++){
                   d[i] = pseudo[i];
                }
                pw.print(String.format(st, d)+";   ");*/
          }
        pw.print(newLine);
    }



    public void multiplyValues(double d) {
        for(int j=0; j<probs.length; j++){
            probs[j] = probs[j]*d;
        }
        
    }



    public void put(State newState, double small) {
        probs[newState.getIndex()]= small;
        
    }



   /* public void doubleSt() {
       double[] newp = new double[(probs.length-1)*2+1];
       this.counts = new double[newp.length];
       double[] newps = new double[newp.length];
       newp[0] = probs[0];
       newps[0] = pseudo[0];
       for(int i=1; i<probs.length; i++){
           newp[i] = probs[i]/2.0;
           newp[i+probs.length-1] = probs[i]/2.0;
           newps[i] = pseudo[i]/2.0;
           newps[i+probs.length-1] = pseudo[i]/2.0;
       }
        this.probs = newp;
        this.pseudo = newps;
       // if(Math.abs(1.0-this.sum())>0.001) throw new RuntimeException("!!");
    }


    public void extend(double small, Double[] ps, Dirichlet newStates) {
        Double [] d = newStates.dist;
        this.pseudo = new double[ps.length];
        if(pseudo.length!=d.length+probs.length) throw new RuntimeException("!!");
       arraycopy(ps, 0, pseudo, 0, ps.length);
        for(int i=0; i<probs.length; i++){
            probs[i] = probs[i]*(1-small);
        }
        double[] newprobs = new double[pseudo.length];
        System.arraycopy(probs, 0, newprobs, 0, probs.length);
        Double[] samp = newStates.sample();
        for(int i=0; i<samp.length; i++){
            newprobs[i+probs.length] = samp[i]*small;
        }
        this.probs = newprobs;
        this.counts = new double[probs.length];
        
    }




    public void fix() {
        Arrays.fill(pseudo, 0);
        if(probs[0]<0.5){
            probs[0] = 0;
            probs[1] = 1;
        }
        else {
            probs[0] =1;
            probs[1] =0;
        }
        
    }*/



   /* public void correct() {
        if(fixed) return;
        else if(probs[0] < 0.01) {
            probs[0] =0.01;
            probs[1] = 0.99;
            fixed = true;
        }
        else if(probs[1] < 0.01){
            probs[1] = 0.01;
            probs[0] = 0.99;
            fixed=true;
        }
        
    }*/

    public double getProbs(int[] nullIndices, double[] probs2) {
        double sum=0;
        for(int i=0; i<nullIndices.length; i++){
            probs2[i] = this.probs[nullIndices[i]];
            sum+=probs2[i];
        }
        return sum;
        
    }

    public int mostLikely() {
      int m_i =0;
      double max = this.probs[0];
      for(int j=1; j<probs.length; j++){
          if(probs[j]>max){
              max = probs[j];
              m_i = j;
          }
      }
      return m_i;
    }
  /* public void transferProbToPseudo() {
       pseudo = new double[probs.length];
       System.arraycopy(probs,0, pseudo, 0, probs.length);
        
    }*/
    public  static double logToProb(double[] probs) {
       double max = Double.NEGATIVE_INFINITY;
       for(int i=0; i<probs.length; i++){
           if(probs[i]>max){
               max = probs[i];
           }
       }
      // if(max >= -30) max =0;
      
       for(int i=0; i<probs.length; i++){
           probs[i] = Math.exp(probs[i] - max);
       }
       
        return -max;
    }
    public double getProbs(int to) {
       return probs[to];
    }
 /*   public double getPseudo(int to) {
        return pseudo[to];
    }*/
  
    public double[] probs() {
       return probs;
    }
 /*   public double[] pseudo() {
      return pseudo;
    }
    public void fillPseudo(double d) {
        Arrays.fill(pseudo, 0);
        
    }*/
    public void setProbs(int to, double d) {
      this.probs[to] = d;
        
    }
   /* public void setPseudo(int to, double d_ps) {
      this.pseudo[to] = d_ps;
        
    }*/
    public PseudoDistribution clone() {
        SimpleExtendedDistribution res =  new SimpleExtendedDistribution(this );
        res.emstsp = this.emstsp;
        return res;
    }
    public PseudoDistribution clone(double swtch) {
        SimpleExtendedDistribution dist1 =  new SimpleExtendedDistribution(this.probs, swtch);
        dist1.emstsp = this.emstsp;
        return dist1;
    }
    public void validate() {
       
        double sum = sum();
        if(Math.abs(sum-1.0) > SimpleDistribution.tolerance) {
        	if(Constants.CHECK) {
        		throw new RuntimeException("!!! "+sum);
        	}
        	for(int i=0; i<this.probs.length; i++){
        		probs[i] = probs[i] /sum;
        	}
        	
        }
        
    }
    public void validate(boolean check) {
        if(!check) this.normalise();
    	
        else this.validate();
        
    }

    public void apply(int[] transformation) {
        double[] cp = new double[transformation.length];
        System.arraycopy(probs, 0,cp, 0, probs.length );
        for(int j=0; j<transformation.length; j++){
            probs[transformation[j] ] = cp[j];
        }
        System.arraycopy(counts, 0,cp, 0, probs.length );
        for(int j=0; j<transformation.length; j++){
            counts[transformation[j] ] = cp[j];
        }
        
    }
    public void addCount(int obj_index, double value) {
       // try{
    /*	if(Constants.CHECK && Double.isInfinite(value)) {
    		throw new RuntimeException("!!");
    	}
      */
    	counts[obj_index]+=value;
       // }catch(ArrayIndexOutOfBoundsException exc){
       //     Logger.global.info("h");
       // }
        
    }
    public double probs(int obj_i) {
       return probs[obj_i];
    }
    public int getMax() {
        return  Constants.getMax(probs);
    }
    public int[] getOrder(){
    	return Constants.getOrder(probs);
    }
    public double[] counts() {
       return counts;
    }
    public void setCounts(int i1, double cnt) {
       counts[i1] = cnt;
        
    }
    public void setEmStSp(EmissionStateSpace stsp) {
    	this.emstsp = stsp;
    	}
   
    public Integer fixedInteger() {
        return null;
       
    }
    public PseudoDistribution makeFixed() {
        int ind = Constants.getMax(probs);
        if(probs[ind]>Constants.fixedThresh()){
            return new IntegerDistribution(ind, this.emstsp);
        }
        else return null;
    }
    public void swtch(double[] d, int[] alias){
        double[] tmp = new double[d.length];
        System.arraycopy(d, 0, tmp, 0, d.length);
        for(int i=0; i<d.length; i++){
            d[alias[i]] = tmp[i];
        }
    }
    public PseudoDistribution swtchAlleles() {
    	double[] probs1 = new double[probs.length];
    	System.arraycopy(probs, 0, probs1, 0, probs1.length);
    	int[] switchTranslation = this.emstsp.getSwitchTranslation();
        
    	  swtch(probs1, switchTranslation);
        SimpleExtendedDistribution dist = new SimpleExtendedDistribution(probs1, Double.POSITIVE_INFINITY);
    	dist.setEmStSp(this.emstsp);
        return dist;
    }
   
  
   

    @Override
    public void setFixedIndex( int k) {
        double[] probs = probs();
        int max_index = Constants.getMax(probs);
        double prob_k = probs[k];
        double prob_max = probs[max_index];
        probs[k] = prob_max;
        probs[max_index] = prob_k;
        
    }
    public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
       Arrays.fill(probs, 0.0);
       for(int i=0; i<tmp.length; i++){
           for(int j=0; j<probs.length; j++){
               probs[j]+=((SimpleExtendedDistribution)tmp[i]).probs(j);
           }
          
       }
       this.normalise(probs);
    }
   
    public void addCounts(lc1.stats.ProbabilityDistribution probabilityDistribution){
        SimpleExtendedDistribution dist = (SimpleExtendedDistribution) probabilityDistribution;
        for(int i=0; i<dist.counts.length; i++){
            this.counts[i]+=dist.counts[i] ;
        }
    }
    public void addCounts( SimpleExtendedDistribution dist,
            Map<Integer, Integer> map) {
         for(int i=0; i<dist.counts.length; i++){
            this.counts[map.get(i)]+=dist.counts[i];
        }
        
    }
    
    public double[] getCount(double[] angle){
        if(angle.length==probs.length) return probs;
       int ratio = (int) Math.round((double)probs.length / (double) angle.length);
       if(Math.abs((double)counts.length - ((double)ratio)*angle.length) > 0.0001) throw new RuntimeException("!!");
       double[] res = new double[angle.length];
       Arrays.fill(res, 0.0);
       for(int i=0; i<angle.length; i++){
           for(int j=0; j<ratio; j++){
               res[i]+=probs[i*ratio+j];
           }
       }
       return res;
    }
   
    public  String getCompressedDataString(EmissionStateSpace emStSp){
     //   StringBuffer sb = new StringBuffer();
    	  if(emStSp==null) return "NaN	NaN";
        double[] probs = this.probs;
       // boolean first = true;
        double thresh = Math.min(Constants.printThresh(), probs[Constants.getMax(probs)]);
        List<DoubleInt> l = new ArrayList<DoubleInt>();
      
        for(int i=0; i<probs.length; i++){
            if(probs[i]>=thresh){
            	l.add(new DoubleInt(i,probs[i], emStSp.getHaploPairString(emStSp.get(i))+"="+String.format("%5.3f", probs[i]).trim()));
            }
        }
        if(l.size()==0){
        	System.err.println("was zero for "+Constants.print(probs));
        	int i = Constants.getMax(probs);
        	l.add(new DoubleInt(i,probs[i], emStSp.getHaploPairString(emStSp.get(i))+"="+String.format("%5.3f", probs[i]).trim()));
        }
        Collections.sort(l);
       String res =  getString(l);
       return res;
         //throw new RuntimeException("!!");
     }
    
    private String getString(List<DoubleInt> l) {
    	if(l.size()==0) return "";
		StringBuffer sb = new StringBuffer(l.get(0).toString());
		for(int i=1; i<l.size(); i++){
			sb.append(";");
			sb.append(l.get(i));
		}
		return sb.toString();
	}

	static class DoubleInt implements Comparable{
    	int i1;
    	public DoubleInt(int i, double e, String string) {
			this.d = e;
			this.i = string;
			this.i1 = i;
		}
		double d; String i;
		public String toString(){
			return i;
		}
		@Override
		public int compareTo(Object o) {
		   double d1 = ((DoubleInt)o).d;
		   if(d1<d) return -1;
		   else if(d1>d) return 1;
		   else return 0;
		}
		
    }

    public String compressedStringHeader(EmissionStateSpace emStSp) {
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<emStSp.size(); i++){
            if(i>0) sb.append("\t");
            sb.append( emStSp.getHaploPairString(emStSp.get(i)));
           
        }
        return sb.toString();
       }


	public void normalise() {
		this.normalise(probs);
		
	}


	@Override
	public double scoreB(int j, int i) {
	return this.probs[j];
	}
	
	@Override
	public double scoreBR(EmissionStateSpace emstsp, int j , int i) {
	//	int cn = emstsp.getCN(j);
		//if(!this.emstsp)
		
		Integer v = this.emstsp.get(emstsp.get(j));
		if(v==null) return 0;
			return this.probs[v];
	
	}
	public void setCounts(double d, Double[] sample) {
	for(int i=0; i<counts.length; i++){
		counts[i] = d *sample[i].doubleValue();
	}
		
	}
	//d is for relationship
	public double evaluate( double e, DoubleMatrix2D nullspace, int[][] groupToState) {
		double[] tmp = new double[groupToState.length];
		if(nullspace==null) this.fill(groupToState, tmp);
		else{
			for(int i=1; i<tmp.length; i++){
				tmp[i] = nullspace.get(i-1, 0);
			}
		}
		this.evaluate(e);
		this.normalise1(groupToState, tmp);
		return this.logProb();
	}
	
	public void fill(int[][] groupToState, double[] tmp){
		
		for(int i=0; i<tmp.length; i++){
			for(int j=0; j<groupToState[i].length; j++){
				tmp[i]+=this.probs[groupToState[i][j]];
			}
		}
	}
	
	public void normalise1( int[][] groupToState, double[] tmp){
		
		for(int i=0; i<tmp.length; i++){
			double sum =0;
			for(int j=0; j<groupToState[i].length; j++){
				sum+=this.probs[groupToState[i][j]];
			}
			if(sum>1e-10){
			for(int j=0; j<groupToState[i].length; j++){
				this.probs[groupToState[i][j]] = probs[groupToState[i][j]] *tmp[i]/sum;
			}
			}
		}
	}
	
	public EmissionStateSpace getEmissionStateSpace(){
		return this.emstsp;
	}


	
   
}
