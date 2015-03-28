package lc1.dp.core;

import lc1.util.Constants;


public abstract class AbstractTerm {
  protected double score;
    public int i;
    protected AbstractTerm( int i, double sc){
    	//if(Double.isNaN(sc)) throw new RuntimeException("!!");
    	if(sc<0)throw new RuntimeException("!!");
        this.score =sc;
        this.i = i;
    }
    public double score(){
        return score;
    }
    /** used for scaling */
    public double scaleScore(){
        return this.score;
    }
    public abstract int getBestPath();
    public int compare(Term t){
        if(score > t.score()) return -1;  // i.e. this is more likely
        else if(score==t.score()) return 0;
        else return 1;
    }
   
    public String toString(){
        return this.score+" "+i;
//        (j==0 ? '0' :CharacterModel.getChar(this.j))+" "+this.i;
    }
   
   public void scale(double d) {
	//   double scoreP = score;
	   if(score>0){
	       double score1 = Math.exp(Math.log(score)+d);
	       if(Constants.CHECK &&( Double.isInfinite(score) || Double.isNaN(score))) {
	    	   throw new RuntimeException("!!");
	       }
	       score = score1;
	      // if(score<=0){
	    	//   score = 0;
	      // }
	   }
      //  if(Double.isNaN(score)) throw new RuntimeException("!! "+d+" "+scoreP);
   }
   
   public abstract void setj(int j);
public void clear() {
  i = -1;
  score =0;
    
}
}
