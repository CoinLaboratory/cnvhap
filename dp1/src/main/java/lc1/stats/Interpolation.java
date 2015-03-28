package lc1.stats;
import java.util.HashMap;
import java.util.Map;


public class Interpolation {
    
    public Interpolation(double[] x, double[] mean){
        this.mean = mean;
        this.x = x;
    }
    double[]  mean; //y_vals;
    double[] x; //x_vals;
    private double getSlope(int x_1, int x_2){
        double y_1 = this.mean[x_1];
        double y_2 = this.mean[x_2];
        return  ((y_2-y_1)/(x_2 - x_1));
    }
    
    Map<Double , Double> sc = new HashMap<Double, Double>();
    
    public double getScore(double x1){
        Double res = sc.get(x1);
        if(res==null){
            sc.put(x1, res = getScore1(x1));
        }
        return res;
    }
    public double getScore2(double x1){
        Double res = sc.get(x1);
        if(res==null){
            return getScore1(x1);
        }
        return res;
    }
    
    private double getScore1(double x1){
        if(x1 < x[0]) throw new RuntimeException(x1+" !! "+x[0]);
        if(x1 > x[x.length-1]) throw new RuntimeException(x1+" !! "+x[x.length-1]);
        if(x1 ==x[0]) return mean[0];
        for(int i=0; i<x.length-1; i++){
            if(x1==x[i+1]) return mean[i+1];
            else if(x1>x[i] && x1<x[i+1]){
                double slope = getSlope(i, i+1);
                return mean[i]+slope*(x1-x[i]);
            }
        }
        throw new RuntimeException("!!");
    }
    
    private double extrapolate(double x, int x1, double slope){
        return mean[x1]+slope * (x - (double) x1);
    }
  
}
