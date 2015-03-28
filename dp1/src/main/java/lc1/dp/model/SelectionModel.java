/*
 * Created on 21-Mar-2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package lc1.dp.model;

import lc1.dp.core.DP;
import lc1.dp.states.SiteEmissionState;
import lc1.dp.states.State;
import lc1.stats.PseudoDistribution;


public class SelectionModel extends FastMarkovModel{
    static final long serialVersionUID = 1;
	public static String[] fromStates = new String[] {"start", "purifying", "positive", "pseudogene"};
	final int[] states_in;
	public static String[][] paramNames = new String[][]{
			new String[] {"positive", "pseudogene"}, 
			new String[] {"positive", "pseudogene"},
			new String[] {"purifying", "pseudogene"},
			new String[] {"purifying", "positive"}
	};
	
	private  static double[][] default_vals = new double[][] {
        /*	new double[] {0.01, 0.001},
			new double[] {0.01, 0.0},
			new double[] {0.2, 0.0},
			new double[] {0.0, 0.0}
          */  
          
                    new double[] {0.05, 0.001},
                    new double[] {0.05, 0.0},
                    new double[] {0.2, 0.0},
                    new double[] {0.0, 0.0}
                    
	};
    
   
    final int[] alias;
	private static void check(){
		for(int i=0; i<default_vals.length; i++){
			if(default_vals[i][0]+default_vals[i][1]>=0.99) throw new RuntimeException("params not valid");
		}
	}
	public static String getHMMString(){
	//	if(usingDefault) return "";
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<default_vals.length; i++){
			for(int j=0; j<default_vals[i].length; j++){
				sb.append("_");
				sb.append(default_vals[i][j]); 
			}
		}
		return sb.toString();
	}
	public static double[] getVals(int k){
		return default_vals[k];
	}
	public static void setVals(double[][] vals1){
		default_vals = vals1;
	}
	public static void setVals(int k, int i, double v){
		default_vals[k][i] = v;
	}
	public final State domdom;
	public final State nucdom;
	public final State protdom;
	
	
   
    public SelectionModel(double[] avgProtDomEmiss, double[] avgNucDomEmiss, double[] avgDomDomEmiss, int[] alias){
        this(default_vals, avgProtDomEmiss, avgNucDomEmiss, avgDomDomEmiss, alias);
        
    }
    /** alias maps position in values list to position in alignment */
	public  SelectionModel(double[][] vals, double[] avgProtDomEmiss, double[] avgNucDomEmiss, double[] avgDomDomEmiss, int[] alias ){
		super("");
        this.alias = alias;
		check();
		domdom = addState(new SiteEmissionState("d", avgDomDomEmiss, alias, false));
		protdom =addState(new SiteEmissionState("p",avgProtDomEmiss, alias, false));
		nucdom = addState(new SiteEmissionState("n",avgNucDomEmiss, alias, false));
        super.initialiseTransitions();
		this.states_in = new int[this.modelLength()];
        for(int i=0; i<states_in.length; i++){
            states_in[i] = i;
        }
	//	java.util.logging.Logger.global.info("recaclulating selection posterior "+Print.toString(vals));
		
		setTransition(MAGIC, domdom, new Double(1.0 -vals[0][0] - vals[0][1]));
		setTransition(MAGIC, protdom, new Double(vals[0][0]));
		setTransition(MAGIC, nucdom, new Double(vals[0][1]));
		
		setTransition(domdom, domdom, new Double(0.99 - vals[1][0] - vals[1][1]));
		setTransition(domdom, protdom, new Double(vals[1][0]));
		setTransition(domdom, nucdom, new Double(vals[1][1]));
		setTransition(domdom, MAGIC, new Double(0.01));
		
		setTransition(protdom, domdom, new Double(vals[2][0]));
		setTransition(protdom, protdom, new Double(0.99 - vals[2][0] - vals[2][1]));
		setTransition(protdom, nucdom, new Double(vals[2][1]));
		setTransition(protdom, MAGIC, new Double(0.01));
		
		setTransition(nucdom, domdom, new Double(vals[3][0]));
		setTransition(nucdom, nucdom, new Double(0.99-vals[3][0] - vals[3][1]));
		setTransition(nucdom, protdom, new Double(vals[3][1]));
		setTransition(nucdom, MAGIC, new Double(0.01));
      // super.fix();
		
	}
	public double[][] getVals(){
		return default_vals;
	}
	
	
	
   
    
	/** alias maps the sequence (which we are interested in) to position in the alignment, can be null */
	public  double[][] getSitePosterior(){
      
        // for(int k=0; k<node.length; k++){
	         DP dp = new DP(this, "selection",  false, alias.length, false);
             dp.reset(false);
	         dp.search(true, false);
            return  dp.getPosteriorMatch();
	}
    @Override
    public Object clone(boolean swtch) {
        throw new RuntimeException("!!");
    }
    @Override
    public boolean converged() {
        // TODO Auto-generated method stub
        return false;
    }
    @Override
    public void setPseudoCountWeights(double[][] d){
        
    }
    @Override
    public void addCounts(PseudoDistribution[] transProbs, int i, int numIndiv) {
        // TODO Auto-generated method stub
        
    }
    @Override
    public void initialiseTransitionCounts() {
        // TODO Auto-generated method stub
        
    }
    @Override
    public void transferCountsToProbs(int index) {
        // TODO Auto-generated method stub
        
    }
    @Override
    public int[] statesIn(int j, int i) {
       return states_in;
    }
    @Override
    public int[] statesOut(int j, int i) {
        // TODO Auto-generated method stub
        return states_in;
    }
    @Override
    public boolean trainEmissions() {
        // TODO Auto-generated method stub
        return false;
    }
   
}