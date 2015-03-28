package lc1.dp.appl;

import java.util.List;

import lc1.dp.data.collection.DataCollection;
import lc1.stats.IlluminaDistribution;
import lc1.stats.PseudoDistribution;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class IndelFilter {
DataCollection obj;
	public IndelFilter(DataCollection obj) {
		this.obj =obj;
	}

	public void filter() {
		int lent = obj.length();
		List<String> samples = obj.indiv();
		double[] bar = new double[samples.size()];
		double[] lrr = new double[samples.size()];
		NormalDistribution normal =new NormalDistributionImpl(0.5,1.0);
		double logLikelihood = 0;
		
		for(int i=0; i<lent;i++){
			for(int k=0; k<samples.size(); k++){
				PseudoDistribution dist = obj.dataL.get(samples.get(k)).emissions(i);
				bar[k] = ((IlluminaDistribution)dist).b();
				lrr[k] = ((IlluminaDistribution)dist).r().doubleValue();
				logLikelihood+= Math.log(normal.density(bar[k]));
			}
		}
	}

}
