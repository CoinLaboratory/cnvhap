/**
 * 
 */
package lc1.stats;

import JSci.maths.statistics.ProbabilityDistribution;

public class TruncatedDistribution extends ProbabilityDistribution{
      ProbabilityDistribution dist1;
      final double min, max;
      final double probLessThanMin, probMoreThanMax;
      final double totalWeight;
      TruncatedDistribution(ProbabilityDistribution  dist, double min, double max){
          this.dist1 = dist;
          this.max = max;
          this.min = min;
          this.probLessThanMin = dist.cumulative(min);
          this.probMoreThanMax = 1.0 - dist.cumulative(max);
          totalWeight = dist.cumulative(max) - dist.cumulative(min);
      }
      

    @Override
    public double cumulative(double arg0) {
        return (dist1.cumulative(arg0) - probLessThanMin) / totalWeight;
    }

    @Override
    public double inverse(double arg0) {
        throw new RuntimeException("!!");
    }

    @Override
    public double probability(double arg0) {
       return dist1.probability(arg0) / totalWeight;
    }
      
  }