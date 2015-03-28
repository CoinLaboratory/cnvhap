package lc1.stats;

import java.util.Arrays;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;

public class SimpleExtendedDistribution1 extends SimpleExtendedDistribution {
public  final double[] pseudo;
@Override
public double[] pseudo(){
	return pseudo;
}
public SimpleExtendedDistribution1(Sampler dir) {
    super(dir);
  this.pseudo = new double[dir.dist.length];
  arraycopy(dir.dist,0, pseudo,0,dir.dist.length  );
  //if(pseudo.length==2 && pseudo[1] > 0.5){
   //    throw new RuntimeException("!!");
 // }
  //if(Math.abs(1.0-Constants.sum(pseudo))>0.001) throw new RuntimeException("!! "+Arrays.asList(dir.dist));
        
}
@Override
public void setProb(double[] prob){
	super.setProb(prob);
    System.arraycopy(prob, 0, pseudo, 0, prob.length);
  
}
public void setPriors(lc1.stats.ProbabilityDistribution dist){
	SimpleExtendedDistribution1 dist1 = (SimpleExtendedDistribution1) dist;
	System.arraycopy(dist1.pseudo, 0, pseudo, 0, pseudo.length);
}

public void change(Sampler dir) {
super.change(dir);
arraycopy(dir.dist,0, pseudo,0,dir.dist.length  );
	
}
public void swtchAlleles(int[] switchTranslation) {
    swtch(pseudo, switchTranslation);
}
public PseudoDistribution clone(double swtch) {
    return new SimpleExtendedDistribution1(this.probs, swtch);
}
public PseudoDistribution clone() {
    return new SimpleExtendedDistribution1(this );
}
@Override
public void setProbs(int to, double d) {
    this.probs[to] = d;
      this.pseudo[to] = d;
  }
public SimpleExtendedDistribution1(int len){
    super(len);
    this.pseudo = new double[len];
    Arrays.fill(pseudo, 1.0 / (double)pseudo.length);
}
    public SimpleExtendedDistribution1(double[] init, double u) {
        super(init, u);
      this.pseudo = new double[init.length];
      System.arraycopy(init,0, pseudo,0,init.length  );
    /*  if(pseudo.length==2 && pseudo[1] > 0.1){
          throw new RuntimeException("!!");
     }*/ 
    }
    public SimpleExtendedDistribution1(double[] init, double u, EmissionStateSpace emstsp) {
        super(init, u);
        this.emstsp = emstsp;
      this.pseudo = new double[init.length];
      System.arraycopy(init,0, pseudo,0,init.length  );
    /*  if(pseudo.length==2 && pseudo[1] > 0.1){
          throw new RuntimeException("!!");
     }*/ 
    }
    public SimpleExtendedDistribution1(SimpleExtendedDistribution1 exp_rd) {
       super(exp_rd);
       this.pseudo = new double[exp_rd.pseudo.length];
       System.arraycopy(exp_rd.pseudo, 0, pseudo, 0, pseudo.length);
     /*  if(pseudo.length==2 && pseudo[1] > 0.1){
           throw new RuntimeException("!!");
      }*/
    }
    
   
    public void transfer(double pseudo){
    	/*if(this.probs.length==2 && probs[1] ==0){
    		Logger.global.info("h");
    	}*/
        if(Constants.CHECK && 
                Double.isNaN(probs[0])) {
        	throw new RuntimeException("!! ");//+print(probs)+"\n"+print(counts)+"\n"+print(this.pseudo));
        }
        double sum=0;
    
            for(int i=0; i<probs.length; i++){
            	probs[i]  = counts[i]+pseudo*(this.pseudo[i]);
                sum+=probs[i];
            }
      
        if(Constants.CHECK && sum==0){
        	throw new RuntimeException("sum is zero");
        }
       if(sum>0){
            for(int i=0; i<probs.length; i++){
                probs[i] = probs[i]/sum;
            }
       }
       else{
    	   for(int i=0; i<probs.length; i++){
    		   probs[i] = this.pseudo[i];
    	   }
       }
       
        

        if(Constants.CHECK){
            try{
            	if(Double.isInfinite(sum)) throw new RuntimeException("!!");
            this.validate();
            }catch(Exception exc){
                exc.printStackTrace();
            }
            }
        if(Constants.CHECK  && 
                Double.isNaN(probs[0])){
        	throw new RuntimeException("!! "+pseudo+"\n"+print(probs)+"\n"+print(counts)+"\n"+print(this.pseudo));
        }
    }
    public String print(double[] d){
        StringBuffer sb = new StringBuffer();
        Double[] d1 = new Double[d.length];
        for(int i=0; i<d.length; i++){
            sb.append("%5.3g ");
            d1[i] = d[i];
        }
        return String.format(sb.toString(), d1);
    }
}
