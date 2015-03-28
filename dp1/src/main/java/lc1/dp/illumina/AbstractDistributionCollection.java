package lc1.dp.illumina;

import java.awt.Color;
import java.awt.Shape;
import java.io.File;
import java.io.Serializable;
import java.util.List;

import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.PseudoDistribution;


public abstract class AbstractDistributionCollection implements Serializable {

	public abstract ProbabilityDistribution getDistribution(short data_index,
			int cn_bottom, int cn_top, int pos);

	public abstract double scoreB(short data_index, int j,double b,  int i);

	public abstract double scoreR(short data_index, int backgroundCount,
			int no_cop,double r,  int i);

	public abstract void addBCount(short data_index, 
			int j, double weight, double val, int i);

	public abstract void addRCount(short data_index, int backgroundCount,
			int no_cop,double weight,  double val, int i);

	/** gets b distribution */
	public abstract ProbabilityDistribution getDistribution(short data_index,
			int j, int i);

	/** maximisation over r dists
	 * @param pseudo1 */
	public abstract void maximisationStep(final int i, final double []pseudo,
			final double []pseudoGlobal,
			 List tasks);

	public abstract void print(File pw);

	
	public abstract String getBName(int ij);

	public abstract Double sampleR(int data_index, int cn_bg, int cn_fg, int pos);

	public abstract Double sampleB(int data_index, int obj_index, int pos);

	//public abstract Color[] getColB();
	//public abstract Color[] getColR();
	
	/*public  Color[] getColR(int bg){
		return getColR();
	}*/

	public abstract ProbabilityDistribution getDistributionBfrac(short i2,
			int i, int j);

	public abstract ProbabilityDistribution getDistribution1(short data_index,
			int i, int j, int i2);

	public abstract void initialise();

	//public abstract Shape[] getShapeB();

	public abstract void addRBCount(short data_index,  int j,
			double weight, double valR, double valB, int i) ;
	public abstract ProbabilityDistribution2 getDistributionRB(short data_index, int n,int noB,  int i) ;
	public abstract ProbabilityDistribution2 getDistributionRBGlob(short data_index, int n,int noB) ;
	public abstract double  scoreRB(short data_index, int j,double r, double b, int i);

	public abstract Double getFrac(int i, int ii, boolean r);
	public abstract Double getFracGlob(int ind, boolean po, boolean R);

	public abstract Double minQuality(int relative_position);

	public abstract void addRBCount(short data_index, int noCop, int noB, double val,
			Double r, Double b, int i);

	
	public abstract void print(int i);
	
	String[][] r_formulae;
	public String[] getForm(int i2, int i) {
		// TODO Auto-generated method stub
		return null;
	}

	public void addMixture(HaplotypeEmissionState st, Boolean[] booleans) {
		//throw new RuntimeException("!! not implemented");
		
	}

	public abstract String[] getFormGlob(int ind) ;
}