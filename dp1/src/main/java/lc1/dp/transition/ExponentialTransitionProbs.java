package lc1.dp.transition;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import lc1.stats.PseudoDistribution;
import lc1.stats.Sampler;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;

public class ExponentialTransitionProbs extends AbstractTransitionProbs {

   public ExpTransProb exp;
   public ExpTransProb alpha;
   public double[] getAlphaPrior(){
	   return alpha.getPrior();
   }
   public String info() {
       StringBuffer info = new StringBuffer(super.info());
       info.append("exp "+exp.info()+"  alpha "+alpha.info());
       return info.toString();
   }
    /*st is from, state is to 
    public void addCount(int from, int to, double val){
        SimpleExtendedDistribution exp_rd = this.getExp(from);
        SimpleExtendedDistribution alpha = (SimpleExtendedDistribution)this.getAlpha(from);
        if(to!=from){//!state.equals(st)){
            alpha.counts[to]+=val;
            exp_rd.counts[1]+=val;
        }
        else{
           double exp =exp_rd.probs[0];
           double non_jump_prob = exp;
           double jump_prob = (1 - exp)*alpha.probs[to];
           double alloc = val*(jump_prob/(jump_prob+non_jump_prob));
           alpha.counts[to]+= alloc;
           exp_rd.counts[1]+=alloc;
           exp_rd.counts[0]+=(val-alloc);
        }
    }*/
   
  
   
    
    /*st is from, state is to */
    public void addCount(int from, int to, double val){
        SimpleExtendedDistribution exp_rd =exp.getExp(from);
        PseudoDistribution alpha = this.alpha.getExp(from);
        if(to!=from){//!state.equals(st)){
            alpha.addCount(to,val);
            exp_rd.counts[1]+=val;
        }
        else{
           double exp =exp_rd.probs[0];
           double non_jump_prob = exp;
           double jump_prob = (1 - exp)*alpha.probs(to);
           double alloc = val*(jump_prob/(jump_prob+non_jump_prob));
           alpha.addCount(to,alloc);
           exp_rd.counts[1]+=alloc;
           exp_rd.counts[0]+=(val-alloc);
        }
    }
    
    public  void addCount(int index, int from, int to, double val){
        SimpleExtendedDistribution exp_rd =exp.getExp(index);
        PseudoDistribution alpha = this.alpha.getExp(index);
        if(to!=from){//!state.equals(st)){
            alpha.addCount(to,val);
            exp_rd.counts[1]+=val;
        }
        else{
           double exp =exp_rd.probs[0];
           double non_jump_prob = exp;
           double jump_prob = (1 - exp)*alpha.probs(to);
           double alloc = val*(jump_prob/(jump_prob+non_jump_prob));
           alpha.addCount(to,alloc);
           exp_rd.counts[1]+=alloc;
           exp_rd.counts[0]+=(val-alloc);
        }
    }
    
    public ExponentialTransitionProbs(ExpTransProb exp, ExpTransProb alpha){
        this.exp = exp;
        this.alpha = alpha;
    }
    public ExponentialTransitionProbs( Sampler samplerFirst, Sampler exp_p, int len) throws Exception{
        this.alpha = new SimpleExpTransProb(samplerFirst, len);
        this.exp = new SimpleExpTransProb(exp_p, len);
    }

    public ExponentialTransitionProbs(ExponentialTransitionProbs probs, boolean swtch){
       if(probs.exp instanceof SimpleExpTransProb){
           exp = probs.exp.clone(swtch, probs.alpha.noStates());
           alpha = probs.alpha.clone(false, probs.alpha.noStates());
       }
       else if(probs.alpha instanceof SimpleExpTransProb){
           exp = probs.exp.clone(false, probs.alpha.noStates());
           alpha = probs.alpha.clone(swtch, probs.alpha.noStates());
       }
       else{
       
         exp = probs.exp.clone(false, probs.alpha.noStates());
          alpha = probs.alpha.clone(false, probs.alpha.noStates());
       }
    }
    public ExponentialTransitionProbs(ExponentialTransitionProbs probs,  int[] statesToGroup, double[] u) {
        this.alpha = probs.alpha.clone(statesToGroup, u);
        this.exp = probs.exp.clone(statesToGroup,u);
    }
  
    @Override
    public AbstractTransitionProbs clone(boolean swtch){
        if(Constants.useFree() && swtch){
            if((exp instanceof MultiExpProbs) && (alpha instanceof MultiExpProbs)){
                return new FreeTransitionProbs1(this);
            }
           
        }
        return new ExponentialTransitionProbs(this, swtch);
    }
    @Override
    public AbstractTransitionProbs clone(int[] statesToGroup, double[] u) {
        
        return new ExponentialTransitionProbs(this, statesToGroup, u);
    }
   
    /* (non-Javadoc)
     * @see lc1.dp.AbstractTransitionProbs#getDistributions()
     */
    public Collection getDistributions(){
        List l = new ArrayList();
        l.add(alpha.getExpRdColl());
        l.addAll(exp.getExpRdColl());
        return l;
//           return Arrays.asList(new Object[] {alpha, exp_rd});
     }
     
    public double getTransition(int from, int to) {
        SimpleExtendedDistribution exp_rd = exp.getExp(from);
        PseudoDistribution alpha = this.alpha.getExp(from);
        double exp = exp_rd.probs[0];
        double toProb = alpha.probs(to);
        /*if(Constants.limitTransByParent()){
        	int noc = Constants.noCopies[0];
        	toProb = Math.abs(to-from)<noc && (from<=noc && to<=noc || from >noc && to>noc) ? 0.5:0;
        }*/
        /*if(Constants.CHECK && (Double.isNaN(toProb) || Double.isInfinite(toProb))){
          	throw new RuntimeException("!!");
          }
        */
        if(to==from){
            return exp+(1-exp)*toProb;
        }
        else{
            return (1-exp)*toProb;
        }
    }
    /** index allows flexibility of having multiple models */
    public  double getTransition(int index, int from, int to){
        SimpleExtendedDistribution exp_rd = exp.getExp(index);
        PseudoDistribution alpha = this.alpha.getExp(index);
        double exp = exp_rd.probs[0];
        double toProb = alpha.probs(to);
       
        if(to==from){
            return exp+(1-exp)*toProb;
        }
        else{
            return (1-exp)*toProb;
        }
    }
  

    public double getRate(int i) {
		// TODO Auto-generated method stub
    	SimpleExtendedDistribution dist = exp.getExp(0);
    	double v= dist.probs[1];
    	//System.err.println(i+" "+dist.probs[0]+" "+dist.probs[1]);
    	return Math.pow(10,v);
		//return Double.parseDouble(AbstractTransitionProbs.transform(,1000));
	}
  

    public void initialiseCounts(boolean start, boolean end) {
        exp.initialiseExpRd();
        alpha.initialiseExpRd();
       //this.exp_rd.initialise();
//       this.alpha.initialise();

    }
 
    
    
   

   /* public void transfer(double pseudoTrans, double[] pseudoExp) {
        exp.transferExp(pseudoExp);
        alpha.transferExp(pseudoTrans);
        //exp_rd.transfer(pseudoExp);
       
    }*/

    /*@Override
    public void transfer(double[] pseudoTrans, double[] pseudoExp) {
        exp.transferExp(pseudoExp[0]);
        alpha.transferExp(pseudoTrans[0]);
        
    }*/

    @Override
    public double transfer(double[] pseudoC, double[][] d, int pos_index) {
    //	Constants.expModelIntHotSpot1();
    	//double[][] pseudoC = new double[d.length][2] ;
    	/*for(int i =0; i<pseudoCExp.length; i++){
    		pseudoC[i][0] = pseudoCExp[pos_index]*Math.exp(d[i]);
    		pseudoC[i][1] = pseudoCExp[pos_index] -pseudoC[i][0];
    	}*/
    	if(d==null) {
    		exp.transferExp(pseudoC[pos_index]);
    		/*if(((SimpleExpTransProb)exp).exp_rd1.probs[1]>0.3){
    			Logger.global.info("h");
    		}*/
    		return 1.0;
    	}
    	return this.exp.evaluateExpRd(pseudoC[pos_index], d[pos_index]);// + this.alpha.evaluateExpRd(pseudoTrans);
    	 /* for(int kk=0; kk<pseudoCExp.length; kk++){
         pseudoC[kk][0] = pseudoCExp[kk]*Math.exp(d*arg0);
         pseudoC[kk][1] =pseudoCExp[kk]-pseudoC[kk][0];
    }*/
    		
    }
    public double  transferAlpha( double[] pseudoTrans, double[][] alpha_overall, int pos_index) {
	
    	if(alpha_overall!=null)return this.alpha.evaluateExpRd(alpha_overall, pseudoTrans[pos_index]);
    	else {
    		this.alpha.transferExp(pseudoTrans[pos_index]);
    		return 1;
    	}
		
	}

    

    public double transitionDistance(AbstractTransitionProbs probs) {
        throw new RuntimeException("!!");
//        ExponentialTransitionProbs probs1 = (ExponentialTransitionProbs)probs;
 //       return alpha.KLDistance(probs1.alpha) + exp_rd.KLDistance(probs1.exp_rd);
    }

    public void print(PrintWriter pw, Double[] hittingProb, double dist) {
        exp.printExp(pw, dist, "exp");
    // pw.print(transform(((Siexp., dist));
        pw.print("; ");
        alpha.printExp(pw, dist, "alpha");
      
        pw.println();
          
      }
   

    
   

    public String toString(){
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        print(pw, null, 1);
        return sw.getBuffer().toString();
    }

    public void validate() {
      alpha.validateExp();
       exp.validateExp();
     super.validate(false);
//       exp_rd.validate();
    }

  /*  @Override
    public void transfer(double pseudoTrans, double pseudoExp) {
        exp.transferExp(pseudoExp);
        alpha.transferExp(pseudoTrans);
    }*/

  
  
  
 





@Override
public int noStates() {
    return ((SimpleExtendedDistribution)alpha.getExp(0)).probs.length;
}

/*public double evaluateExpRd(double[][] pseudoC) {
  return exp.evaluateExpRd(pseudoC);
}*/
public static AbstractTransitionProbs get(Object clazz, Sampler samplerFirst, Sampler expD, int noStates, double[] hs) throws Exception{
    if(clazz instanceof Class){
        return
   
        (AbstractTransitionProbs) ((Class)clazz).getConstructor(new Class[] {Sampler.class,  int.class}).newInstance(
                new Object[] {samplerFirst,  noStates});
    }
    else{
        Class[] cl = (Class[]) clazz;
        if(cl[0] ==FreeTransitionProbs1.class){
            return new FreeTransitionProbs1(samplerFirst);
        }
    	
        ExpTransProb exp = (ExpTransProb) cl[0].getConstructor(new Class[] {Sampler.class, int.class, boolean.class, double[].class}).newInstance(
                new Object[] {expD, noStates, true, hs});
        ExpTransProb alpha = (ExpTransProb) cl[1].getConstructor(new Class[] {Sampler.class, int.class, boolean.class, double[].class}).newInstance(
                new Object[] {samplerFirst, noStates, false, null});
        return new ExponentialTransitionProbs(exp, alpha);
    }
}







	













  
}
