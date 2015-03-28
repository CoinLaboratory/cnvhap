package lc1.dp.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import lc1.stats.ProbabilityDistribution;
import lc1.stats.SkewNormal;
import lc1.util.Constants;
import lc1.util.Executor;
import pal.math.ConjugateDirectionSearch;
import pal.math.MultivariateFunction;
import pal.math.MultivariateMinimum;
import pal.math.OrthogonalHints;

public class ProbMultivariate implements MultivariateFunction {

    final static double base = Constants.base();
    final static double logbase = Math.log(Constants.base());
    public static ExecutorService es = Executor.getCachedEs(ProbMultivariate.class);// numT);//
   public static int timeout = 100;
    public void minimise(double pseudoM, double pseudoSd, double pseudoSk){
    	/*double sum=0;
        for(int i=0; i<this.single.size(); i++){
        	sum+= this.noSamples.get(i);
            ((ProbabilityDistribution)single.get(i)).setPrior(pseudoM, pseudoSd, pseudoSk);
            single.get(i).updateParamIndex();
        }
        if(sum<Constants.trainThresh()) return;*/
    	;
        final MultivariateMinimum os = new ConjugateDirectionSearch();
        final MultivariateMinimum os1 = new ConjugateDirectionSearch();
        final double[] init = new double[this.totParam];
        final double[] xvec =  new double[this.totParam];
        Logger.global.info("minimizing "+totParam);
        for(int i=0; i<xvec.length; i++){
            xvec[i]= this.getParamValue(i);
        }
        System.arraycopy(xvec, 0, init, 0, xvec.length);
       List run = new ArrayList();
       run.add(new Callable(){
	        public Object call(){
	        	try{
	        		Logger.global.info("optimizing ");
	            os.optimize(ProbMultivariate.this,xvec, 0.01, 0.01);
	        	}catch(Exception exc){
	        		exc.printStackTrace();
	        		   os1.optimize(ProbMultivariate.this,xvec, 0.01, 0.01);
	        	}
	            return null;
	        }
        });
        try{
        	if(true){
        		for(int i=0; i<run.size(); i++){
        			try{
        			((Callable)run.get(i)).call();
        			}catch(Exception exc){
        				exc.printStackTrace();
        			}
        		}
        	}
        	else es.invokeAll(run, timeout, TimeUnit.SECONDS);
        }catch(InterruptedException exc){
        	exc.printStackTrace();
        }
/*        try{
        for(int i=0; i<100; i++){
            Thread.sleep(100);//this.wait(100);
            if(!th.isAlive()) break;
        }
        if(th.isAlive()){
            
            th.stop();
           
        }
       //System.arraycopy(init, 0, xvec, 0, xvec.length);
        }catch(Exception exc){
            exc.printStackTrace();
        }*/
        for(int i=0; i<xvec.length; i++){
           this.setParamValue(i, xvec[i]);
        }
        for(int i=0; i<pairs.size(); i++){
            SkewNormal pair = (SkewNormal) pairs.get(i);
            ProbabilityDistribution[] mem = members.get(i);
            transferExp(pair, mem);
            pair.updateParamIndex();
            
        }
        for(int i=0; i<this.single.size(); i++){
        	single.get(i).recalcName();
        }
    	Logger.global.info("done minimizing ");
    }
    
  public  final List<ProbabilityDistribution> pairs;
  public   final List<ProbabilityDistribution[]> members;
  public   final List<ProbabilityDistribution> single = new ArrayList<ProbabilityDistribution>();
  //  final List<Double> prior_multipliers = new ArrayList<Double>()
    final int[] numParam; //cumulative
    int totParam =0;
    final int len;
    
  //
   
   public void updateSampleSize(double[] noSamples){
       Arrays.fill(noSamples, 0.0);
       for(int i=0; i<members.size(); i++){
           ProbabilityDistribution[] dist = members.get(i);
           for(int k=0; k<dist.length; k++){
               int index = single.indexOf(dist[k]);
               
                   double sum = (((ProbabilityDistribution)pairs.get(i)).sum());
                  noSamples[index]+=sum;
               
           }
       }
    //   System.err.println("sample sizes : "+this.noSamples);
   }
   private void makeSingleList(){
      
       for(int i=0; i<members.size(); i++){
           ProbabilityDistribution[] dist = members.get(i);
           for(int k=0; k<dist.length; k++){
               int index = single.indexOf(dist[k]);
               if(index<0){
                   single.add(dist[k]);
               }
               
           }
       }
    //   System.err.println("sample sizes : "+this.noSamples);
   }
   
    public ProbMultivariate(List<ProbabilityDistribution> s1,
            List<ProbabilityDistribution[]> s1_memb) {
       this.pairs = s1;
      // boolean che = s1.get(0) ==s1_memb.get(0)[0];
     
       this.members = s1_memb;
     this.makeSingleList();
       numParam = new int[single.size()];
       
       for(int i=0; i<numParam.length; i++){
           MultivariateFunction mvf = (MultivariateFunction) single.get(i);
           numParam[i] = mvf.getNumArguments()+totParam;
           totParam+=mvf.getNumArguments();
       }
       len = numParam.length;
    }

   private int findIndex(int n){
       int i=0;
       for(; i<len-1; i++){
           if(n < numParam[i]  ) break;
       }
       return i;
   }

  static  double[] params = new double[3];
   
 /* public static ProbabilityDistribution getExp(Mixture[] mem, String name, int ik){
      // System.err.println("b "+pair.toString());
      double[] lower = new double[] {-5, 1e-3, -1e10};
      double[] upper= new double[] {5, 1e3, 1e10};
      double priorMod = 0.0;
       Arrays.fill(params, 0);
       for(int k=0; k<mem.length; k++){
           SkewNormal m_k =(SkewNormal) mem[k].dist[0];
           priorMod+=m_k.priorModifier;
           for(int j=0; j<params.length; j++){
               if(j==0 && Constants.trainEnsemble()==2){
                   params[j] +=Math.pow(base,m_k.getParamValue(j));
               }
               else{
                   params[j] += m_k.getParamValue(j);
               }
           }
       }
       double fir =  params[0]/(double)mem.length;
      ProbabilityDistribution res = 
           params[2]==0 ? new TrainableNormal(name, Constants.trainEnsemble()==2 ?  (Math.log(fir) / logbase):fir, 
                   
                   params[1],  Constants.round(), priorMod/(double)mem.length): 
       
       new SkewNormal(name,Constants.trainEnsemble()==2 ?  (Math.log(fir)/logbase):fir, 
               
               params[1], params[2], lower, upper, Constants.round(), priorMod/(double)mem.length);
           new SimpleExtendedDistribution(mixt, Double.POSITIVE_INFINITY);
      return new Mixture(res, Constants.minR(),Constants.maxR(), Constants.mixCoeff(ik));
   }*/
  
  private static String getName(ProbabilityDistribution[] mem) {
	StringBuffer sb = new StringBuffer();
	for(int i=0; i<mem.length; i++){
		sb.append(mem[i].name()+"-");
	}
	return sb.toString();
}

/* public static void transfer(ProbabilityDistribution pair, ProbabilityDistribution[] mem){
      // System.err.println("b "+pair.toString());
       Arrays.fill(params, 0);
       for(int k=0; k<mem.length; k++){
           ProbabilityDistribution m_k = (ProbabilityDistribution)mem[k];
           for(int j=0; j<params.length; j++){
               params[j] +=m_k.getParamValue(j);
           }
       }
       for(int i=0; i<params.length; i++){
           pair.setParamValue(i, i==0 ? params[i]/(double)mem.length : params[i]);
       }
     //  System.err.println("aft "+pair.toString());
   }*/
  // static boolean exp = true;
   public static void transferExp(SkewNormal pair, ProbabilityDistribution[] mem){
    // System.err.println("b "+pair.toString());
        Arrays.fill(params, 0);
        for(int k=0; k<mem.length; k++){
           SkewNormal m_k = (SkewNormal)mem[k];
            for(int j=0; j<params.length; j++){
                if(j==0 && Constants.trainEnsemble()==2){
                    params[j] +=Math.pow(base, m_k.getParamValue(j));
                }
                else{
                    params[j] += m_k.getParamValue(j);
                }
               
            }
        }
        for(int i=0; i<params.length; i++){
            if(i==0){
                if(Constants.trainEnsemble()==2){
                    pair.setParamValue(i,Math.log(params[i]/(double)mem.length)/logbase);
                }
                else{
                    pair.setParamValue(i,params[i]/(double)mem.length);
                }
            }
            else{
                pair.setParamValue(i, params[i]);
            }
           
        }
        pair.recalcName();
   //  System.err.println("aft "+pair.toString());
    }
   
    public double evaluate(double[] argument) {
        
        for(int i=0; i<argument.length; i++){
            this.setParamValue(i, argument[i]);
        }
        double logL =0;
        for(int i=0; i<pairs.size(); i++){
           SkewNormal pair = (SkewNormal) pairs.get(i);
            ProbabilityDistribution[] mem = members.get(i);
            transferExp(pair, mem);
            logL+= pair.calcLH();
            
        }
        double prior = 0;
        for(int i=0; i<this.single.size(); i++){
            prior+=((ProbabilityDistribution)single.get(i)).prior();
        }
     //   System.err.println("argument "+lc1.util.Constants.print(argument)+" "+-1*(logL+prior)+" "+logL+" "+prior);
      //  if(logL>-1e-5) throw new RuntimeException(" log is zero "+logL);
        return -1*(logL+prior);
    }

    public double getLowerBound(int n) {
        int i = findIndex(n);
       int n1 = i==0 ? n : n - numParam[i-1];
       return ((MultivariateFunction)single.get(i)).getLowerBound(n1);
    }
    public double getParamValue(int n){
        int i = findIndex(n);
        int n1 = i==0 ? n : n - numParam[i-1];
        return ((ProbabilityDistribution)single.get(i)).getParamValue(n1);
    }
  
    public void setParamValue(int n, double val){
        int i = findIndex(n);
        int n1 = i==0 ? n : n - numParam[i-1];
         ((ProbabilityDistribution)single.get(i)).setParamValue(n1, val);
        
    }

    public int getNumArguments() {
       return totParam;
    }

    public OrthogonalHints getOrthogonalHints() {
        // TODO Auto-generated method stub
        return null;
    }

    public double getUpperBound(int n) {
        int i = findIndex(n);
        int n1 = i==0 ? n : n - numParam[i-1];
        return ((MultivariateFunction)single.get(i)).getUpperBound(n1);
    }

    public String getCompoundName(int i) {
        ProbabilityDistribution[] sn = this.members.get(i);
       StringBuffer sb = new StringBuffer(sn[0].name());
       for(int j=1; j<sn.length; j++){
           sb.append(":"+sn[j].name());
       }
       return sb.toString();
    }

	
    
}
