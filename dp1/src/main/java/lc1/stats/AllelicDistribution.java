package lc1.stats;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.util.Constants;

public class AllelicDistribution extends IlluminaDistribution {

	public AllelicDistribution(short data_index) {
		super(data_index);
		// TODO Auto-generated constructor stub
	}
	public AllelicDistribution(AllelicDistribution r){
		super(r);
	}
	
	public Double b(int i) {
		if(this.r>=Constants.depthPlotThresh())
		return this.b/this.r().doubleValue();
		else return Double.NaN;
	}
	
	
	
	
}
