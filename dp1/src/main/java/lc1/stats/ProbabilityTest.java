package lc1.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.util.Constants;

public class ProbabilityTest {

    public static void main(String[] args){
        ProbabilityTest pt = new ProbabilityTest();
      double bef = Double.NEGATIVE_INFINITY;
        for(int i=0; i<1000; i++){
            System.err.println(i);
            double n = pt.infer();
            if( n  - bef <0.01) break;
            bef = n;
        }
    }
    
    List<JSci.maths.statistics.ProbabilityDistribution> dist = new ArrayList<JSci.maths.statistics.ProbabilityDistribution>();
    
    List<Double> samples = new ArrayList<Double>();
    
    static int lim = 10000;
    static int grid = 100000;
    double[] prob;
   SimpleExtendedDistribution mix;
    ProbabilityTest(){
//        setUp();
        this.setUp1();
    }
    
    public void setUp(){
        dist.add((new SkewNormal(0,0.02,1e10, -10, 10, grid,1.0)));
        dist.add((new SkewNormal(1.0,0.02,-1e10, -10, 10, grid,1.0)));
        dist.add(new TrainableNormal(null, 0.28, 0.02, grid,1.0));
        dist.add(new TrainableNormal(null,0.33, 0.02, grid,1.0));
        dist.add(new TrainableNormal(null,0.48, 0.12, grid,1.0));
        dist.add(new TrainableNormal(null,0.65, 0.02, grid,1.0));
        dist.add(new TrainableNormal(null,0.78, 0.02, grid,1.0));
        mix = new SimpleExtendedDistribution1(dist.size());
        Arrays.fill(mix.counts, 1.0);
        mix.counts[0] =  2.0;
        mix.counts[1]=2.0;
        mix.counts[2]=2.0;
      mix.transfer(0);
      mix.initialise();
      
       System.err.println("mix "+mix.toString());
      //  for(int i=0; i<dist.size(); i++){
           
            for(int j=0; j<lim; j++){
                SkewNormal dist_i = (SkewNormal) dist.get((int) mix.sample());
                double x = dist_i.rsn();
             //   System.err.println(x);
                samples.add(x);
            }
          mix = new SimpleExtendedDistribution1(dist.size());
       // }
       dist.set(0,(new SkewNormal(0.0,0.2,1e10, -10, 10, grid,1.0)));
      dist.set(1,(new SkewNormal(1.0,0.2,-1e10, -10, 10, grid,1.0)));
      dist.set(2, new TrainableNormal(null,0.25, 0.02, grid,1.0));
      dist.set(3, new TrainableNormal(null,0.33, 0.02, grid,1.0));
      dist.set(4, new TrainableNormal(null,0.5, 0.02, grid,1.0));
      dist.set(5, new TrainableNormal(null,0.66, 0.02, grid,1.0));
      dist.set(6, new TrainableNormal(null,0.75, 0.02, grid,1.0));
        prob = new double[dist.size()];
        System.err.println("old params ");
        for(int j=0; j<prob.length; j++){
            System.err.println(dist.get(j).toString());
        }
        
    }
    //new double[]{ -1.0, -0.53,   0, 0.35, 0.54, 1.09, 1.38};
    
    public void setUp1(){
      //  dist.add((new SkewNormal(-1.0,0.1,-10, -10, 10, grid)));
        dist.add((new TrainableNormal(null,-1,0.2, grid,1.0)));
        dist.add(new TrainableNormal(null,1, 0.2, grid,1.0));
        //dist.add(new TrainableNormal(0.35, 0.1, grid));
        //dist.add(new TrainableNormal(0.54, 0.1,  grid));
      //  dist.add(new SkewNormal(1.09, 0.1, 10, -10, 10 , grid));
        mix = new SimpleExtendedDistribution1(dist.size());
        Arrays.fill(mix.counts, 1.0);
    //    mix.counts[2]=5.0;
      mix.transfer(0);
      mix.initialise();
      
       System.err.println("mix "+mix.toString());
      //  for(int i=0; i<dist.size(); i++){
           
            for(int j=0; j<lim; j++){
                SkewNormal dist_i = (SkewNormal) dist.get((int)mix.sample());
                double x = dist_i.rsn();
             //   System.err.println(x);
                samples.add(x);
            }
          mix = new SimpleExtendedDistribution1(dist.size());
       // }
      //    dist.set(0,(new SkewNormal(-0.80,.1,0, -10, 10, grid)));
          dist.set(0,(new TrainableNormal(null,-1,0.1,grid,1.0)));
          dist.set(1,new TrainableNormal(null,-2, 0.1, grid,1.0));
          
          
//          dist.set(0,(new SkewNormal(-1,0.1,1,-10,10,grid)));
//          dist.set(1,new SkewNormal(1, 0.1, 1, -10, 10, grid));
          dist.set(0,(new TrainableNormal(null,0,0.01,grid,1.0)));
          dist.set(1,new TrainableNormal(null,1, 0.01,   grid,1.0));
          
          //dist.set(3,new TrainableNormal(0.40, 0.1, grid));
          //dist.set(4,new TrainableNormal(0.54, 0.1,  grid));
       //   dist.set(5,new SkewNormal(1.09, 0.1, 0, -10, 10 , grid));
        prob = new double[dist.size()];
        System.err.println("old params ");
        for(int j=0; j<prob.length; j++){
            System.err.println(dist.get(j).toString());
        }
        
    }
    public double infer(){
        double likelihood = 0;
       mix.initialise();
        for(int i=0; i<samples.size(); i++){
            Arrays.fill(prob, 0.0);
            double x = samples.get(i);
            for(int j=0; j<prob.length; j++){
                prob[j] = dist.get(j).probability(x);//*mix.probs[j];
                
            }
            likelihood+=Math.log(Constants.sum(prob));
            try{
            SimpleExtendedDistribution.normalise(prob);
            }
            catch(Exception exc){
                Arrays.fill(prob, 1.0/(double)prob.length);
            }
            for(int j=0; j<prob.length; j++){
                JSci.maths.statistics.ProbabilityDistribution distj = dist.get(j);
                mix.addCount(j, prob[j]);
                ((SkewNormal)distj).addCount(x, prob[j]);
            }
        }
        System.err.println("log likelihood is "+likelihood);
        System.err.println("new params ");
        mix.transfer(0);
        System.err.println(mix);
        for(int j=0; j<prob.length; j++){
            JSci.maths.statistics.ProbabilityDistribution distj = dist.get(j);
            if(j>1){
                if(getDist(j).median() < getDist(j-1).median()) throw new RuntimeException("!!");
            }
                
                        ((SkewNormal)distj).maximise(0, 0,0);
            System.err.println(dist.get(j).toString());
        }
      return likelihood;
        
    }
    public SkewNormal getDist(int j){
        return (SkewNormal ) dist.get(j);
    }
    
}
