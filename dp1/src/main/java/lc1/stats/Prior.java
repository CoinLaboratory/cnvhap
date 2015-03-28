package lc1.stats;

import java.io.Serializable;

public interface Prior extends Serializable{
    public double probability(double x);
}
