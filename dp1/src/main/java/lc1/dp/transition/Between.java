package lc1.dp.transition;

import lc1.stats.SimpleExtendedDistribution;

public interface Between {

	double[] getTrans(int groupFrom);

	int[][] groupToState();

	SimpleExtendedDistribution frac(int groupFrom);

	void setProb(int groupFrom, int groupTo, int indexFrom, double prob);
	double  getProb(int groupFrom, int groupTo, int indexFrom);

	double evaluteBack();

	public double getGroupProb(int groupFrom, int groupTo);

	void validate1();

}
