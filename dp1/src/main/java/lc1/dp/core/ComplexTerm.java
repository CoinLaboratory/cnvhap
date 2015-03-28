package lc1.dp.core;

import lc1.util.Constants;

public class ComplexTerm extends Term {

   
    public double[] prob; //should add to 1
    public ComplexTerm(int i, double sc, int modelLength){
        super(i, sc, modelLength);
        this.prob = new double[modelLength];
    }
    /*protected ComplexTerm( double[] prob, int i, double sc){
        super(i,sc);
        this.prob = prob;
        double sum=0;
        for(int k=0; k<prob.length; k++){
            sum+=prob[k];
        }
        if(Math.abs(sum-score)>0) throw new RuntimeException("!!");
    }*/
/* samples a path
    @Override*/
    public int getBestPath() {
        return Constants.sample(prob, score);
    }
   @Override
    public void scale(double d) {
       super.scale(d);
       for(int j=0; j<this.prob.length; j++){
            prob[j] = Math.exp(Math.log(prob[j])+d);
        }
    }
 /*   @Override
    public double scaleScore(){
       
        return prob[this.j];
    }*/
    
}
