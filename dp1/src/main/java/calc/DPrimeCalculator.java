package calc;

import java.util.Arrays;

public class DPrimeCalculator implements LDCalculator {

        public void normalise(double[][]d){
            double sum=0;
            for(int i=0; i<d.length; i++){
                for(int j=0; j<d.length; j++){
                    sum+=d[i][j];
                }
            }
            for(int i=0; i<d.length; i++){
                for(int j=0; j<d.length; j++){
                     d[i][j] = d[i][j]/sum;
                }
            }
            
        }
        
        
        
        public double calculate(double[][] p_ij) {
            normalise(p_ij);
          
            double[] p_i = new double[p_ij.length]; //row sum
            double[] p_j = new double[p_ij.length]; //col sum
            Arrays.fill(p_i, 0);
            Arrays.fill(p_j, 0);
            for(int i=0; i<p_ij.length; i++){
                for(int j=0; j<p_ij.length; j++){
                    p_i[i]+=p_ij[i][j];
                    p_j[j]+=p_ij[i][j];
                }
            }
            
            double[][] r = new double[p_ij.length][p_ij.length];
            for(int i=0; i<p_ij.length; i++){
                for(int j=0; j<p_ij.length; j++){
                    r[i][j] = p_ij[i][j] - p_i[i]*p_j[j];
                }
            }
            double Dminus = -Math.min(p_i[0]*p_j[0], p_i[1]*p_j[1]);
            double Dplis = Math.min(p_i[0]*p_j[1], p_i[1]*p_j[0]);
            double Dprime = r[0][0] > 0 ? r[0][0] / Dplis :
            r[0][0]/Dminus;
            return Dprime;
//            return r[0][0];
            
        }

}
