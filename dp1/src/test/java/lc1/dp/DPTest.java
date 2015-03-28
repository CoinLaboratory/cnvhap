package lc1.dp;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.DP;
import lc1.dp.core.StatePath;
import lc1.dp.data.representation.ScorableObject;
import lc1.dp.model.MarkovModel;
import lc1.dp.model.PairMarkovModel;
import lc1.dp.states.DotState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.State;
import lc1.stats.SimpleDistribution;
import lc1.stats.StateDistribution;

import com.braju.format.Format;



public abstract class DPTest {
    
    static Logger logger = Logger.getLogger("lc1.dp.DPTest");
    {
            logger.setLevel(Level.ALL);
            logger.getParent().getHandlers()[0].setFormatter( new Formatter(){
                public String format(LogRecord record){
                    return 
                            record.getMessage()+"\n";
                }
            });
    }
    
    double pseudocount=0;
    MarkovModel hmm;
    ScorableObject[] data;
    State[][] states;
   
    final int no_copies;
   static Random rand = new Random();
    
   static int count;
    static double[] sampleDistribution(int n){
        double d = 0.5;
        double[] res = new double[n];
        if(true){
            Arrays.fill(res, d/(double)n);
            res[count]+=1-d;
            
        }
        else{
            double total = 0;
            for(int i=0; i<res.length; i++){
                res[i] = rand.nextDouble();
                total+=res[i];
            }
            for(int i=0; i<res.length; i++){
                res[i] = res[i]/total;
            }
        }
        count++;
        return res;
        
    }
    
   
    DotState magic =  new DotState("!", SimpleDistribution.noOffset, SimpleDistribution.noOffset);
     
    public void run() throws Exception{
        MarkovModel hmm1 = (MarkovModel) hmm.clone();
        hmm1.setRandom(randTrans ? u : Double.POSITIVE_INFINITY,
                       randEmiss ? u : Double.POSITIVE_INFINITY, false, false);
        train(hmm1);
    }
   
   public static void main(String[] args){
         try{
             DPTest pt =//new HaplotypeTest(1);
             //new DPDiceTest(10, 1);
                 new GenotypeTest(100,2);
          //  pt.run();
           pt.search();
         }catch(Exception exc){
             exc.printStackTrace();
         }
     }
   abstract MarkovModel getHMM();
   abstract ScorableObject makeScoreableObject(int k);
     
     double[] bin =new double[] {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
     double[] noright = new double[bin.length];
     double[] nowrong = new double[bin.length];
    
   public  DPTest(int no, int no_copies1) throws Exception{
        this.no_copies = no_copies1;
        data = new ScorableObject[no];
        Arrays.fill(noright, 0);
        Arrays.fill(nowrong, 0);
        states = new State[data.length][];
        hmm = getHMM();
     
        //  hmm.setUniform();
        for(int k=0; k<data.length; k++){
          //  logger.info("adding data point ..."+k+" of "+data.length);
            states[k] = hmm.emitStatePath(hmm.MAGIC, EmissionState.class);
           // logger.info("state path "+Arrays.asList(states[k]));
            data[k] = makeScoreableObject(k);
            for(int i=0; i<states[k].length; i++){
                addPoint( k, i);
            }
        }
        logger.info("done adding data");
    }
   protected void addPoint( int k, int i) {
       EmissionState t = (EmissionState)states[k][i];
       data[k].addPoint( t.sample(i));
    }
    


    boolean randTrans = true;
    boolean randEmiss = true;
    double u = 10 ;
    
    public void print(MarkovModel hmm){
        MarkovModel hmm1 = hmm instanceof PairMarkovModel ? ((PairMarkovModel)hmm).m1 : hmm;
        System.err.println(hmm1);
    }
    
    public void train(MarkovModel hmm1) throws Exception{
        
      this.logger.info("TRAINING ");
        if(printHMM){
            System.err.println("HMM:*************************");
            print(hmm);
            System.err.println("CF HMM1*************************");
            print(hmm1);
        }
        BaumWelchTrainer bwt = new BaumWelchTrainer(hmm1, data);
        int num=0;
        double bestScore = getBestScore();
        double logProb = Double.NEGATIVE_INFINITY;
        this.logger.info("initial kl dist is "+hmm.totalTransitionEmissionDist(hmm1)+" "+hmm1.totalTransitionEmissionDist(hmm));
        
        while(true){
          double lp =   bwt.trainIteration(randTrans,randEmiss);
          if(lp+1 < logProb) throw new Exception("warning: log prob should always be increasing "+lp +" should be > "+logProb);
          this.logger.info("kl dist is "+hmm.totalTransitionEmissionDist(hmm1)+" "+hmm1.totalTransitionEmissionDist(hmm));
          logger.info("\n log prob is "+lp+" cf "+bestScore);
          if(printHMM){
              System.err.println("HMM:*************************");
              print(hmm);
              System.err.println("cf HMM1*************************");
              print(hmm1);
          }
          if(Math.abs(lp-logProb) < 0.01)num++;
          if(num>1) break;
          logProb = lp;
         // break;
        }
    }
    boolean printHMM = true;
public double getBestScore(){
    BaumWelchTrainer bwt = new BaumWelchTrainer(hmm, data);
    return bwt.trainIteration(false, false);
}
    private String getPrint(int length, String str){
        StringBuffer sb = new StringBuffer(length*str.length());
        for(int i=0; i<length; i++){
            sb.append(str);
        }
        return sb.toString();
    }
   public void search() throws Exception{
       MarkovModel hmm1 = hmm;
       for(int l=0; l<data.length; l++){
            String[] sb = new String[states[l].length];
            String[] sb1 = new String[states[l].length];
            int k=0;
            int maxLength =0;
            for(int i=0; i<states[l].length;i++){
                    sb[i] = states[l][i].getName();
                    if(sb[i].length()>maxLength) maxLength = sb[i].length();
                if(states[l][i] instanceof EmissionState) {
                        sb1[i] = (data[l].getStringRep(k));
                        if(sb1[i].length()>maxLength) maxLength = sb1[i].length();
                        k++;
                }
                else{
                    sb1[i] = "_";
                }
            }
            String printStr = getPrint(sb.length, "%"+maxLength+"s ");
                logger.info("Actual states "+Format.sprintf(printStr, sb));
            logger.info("Emissions     "+Format.sprintf(printStr, sb1)); 
            boolean logspace = true;
            DP dp  =new DP(hmm1,"", data[l], logspace) ; 
            double sc1 = logspace ? dp.searchViterbi() : dp.search(true);
           StatePath sp =  dp.getStatePath();
               logger.info("state path is "+Format.sprintf(printStr,sp.getWord()));
         
           if(!logspace){
       for(int i=0; i<data[l].length(); i++){
           try{
               StateDistribution dist = dp.getPosterior(i);
              int max = dist.getMax();
              State maxSt = this.hmm.getState(max);
              double maxVal = dist.get(maxSt);
               for(int j=0; j<bin.length; j++){
                   if(maxVal<bin[j]){
                       if(maxSt==this.states[l][i]) noright[j]+=1.0;
                       else nowrong[j]+=1.0;
                   }
               }
             //  double sum = dist.sum(SimpleEmissionState.class);
               //System.err.println("sum "+sum);
               
           }catch(Exception exc){
               Exception exc1 = new Exception("problem with position "+i+" of "+data[l].length());
               exc1.initCause(exc);
               exc1.printStackTrace();
               // throw exc1;
               
           }
       }
     
       for(int i=-1; i<data[l].length(); i++){
           StateDistribution[] dist = dp.getTransitionPosterior(i);
           double sum = 0;
           for(Iterator<StateDistribution> it = Arrays.asList(dist).iterator(); it.hasNext();){
               StateDistribution dst = it.next();
               sum+=
                   dst.sum();
              
           }
          // logger.info("sum "+sum);
           if(Math.abs(1-sum)>0.001) throw new Exception("posterior sum not right!");
       }
           }
       // logger.info( st);
       }
       for(int i=0; i<bin.length; i++){
           System.err.println(bin[i]+" "+noright[i]/(noright[i]+nowrong[i]));
       }
   }
   
   
      //  logger.info(toString(sum(match)));
    //    lc1.dp.StatePath[] path = dp.getStatePath().splitStatePaths(hmmAlign.hmmA.beginIndex, hmmAlign.hmmA.endIndex);
    static int[] sum(double[][] match){
        int[] res = new int[match[0].length];
            for(int i=0; i<match[0].length; i++){
                double sum =0;
                for(int j=0; j<match.length; j++){
                    sum+=match[j][i];
                }
                res[i] = (int) Math.round(sum);
            }
        return res;
    }
    static String toString(int[] m){
        StringBuffer sb = new StringBuffer(m.length);
        for(int i=0; i<m.length; i++){
            sb.append(m[i]);
        }
        return sb.toString();
    }
   
}
