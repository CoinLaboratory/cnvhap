package lc1.dp.data.classification;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.states.HaplotypeEmissionState;

public abstract class ErrorClassificationAbstract {

    static boolean exclCopyMisMatch =true;
    protected final List<Comparable> emSt;
    protected StringBuffer formatStr = new StringBuffer();
    protected StringBuffer formatStr1 = new StringBuffer();
    private  int[][] errorcase; //first index is different sources of data
    
   String original_name, inferred_name;
    static int no_bins = 10;
    private static double incr  = 1.0/ (no_bins);
   
    private int[][] error ;
  
    protected final double thresh;
    protected final ErrorClassificationAbstract parent;
    public final PrintWriter log;
    public void setNames(String orig, String inf){
    	this.original_name = orig;
    	this.inferred_name = inf;
    }
    public ErrorClassificationAbstract(List<Comparable> ss, 
            double threshold, ErrorClassificationAbstract parent, Integer len1,
         
            PrintWriter log) {
        this.log = log;
        int len = len1==null ? ss.size() : len1;
        this.emSt = ss;
        this.parent = parent;
        this.thresh = threshold;
      
        for(int i=0; i<len; i++){
            //     stateSpaceIndex.put(ss.get(i), i);
                 formatStr.append("%10i ");
                 formatStr1.append("%10s ");
              
             }
             formatStr.append("%7i ");
             formatStr1.append("%7s ");
             this.errorcase = init(len);
             this.error = new int[no_bins+1][2];
           
             
    }


    protected int[][] init(int i) {
      int[][] res = new int[i][i];
           for(int j=0; j<i; j++){
               Arrays.fill(res[j],0);
           }
        return res;
    }
   
 
    public void compare(HaplotypeEmissionState orig, HaplotypeEmissionState pred
    ) {
        int[] fromTo = new int[] {0,0};
        
       // for(int k=0; k<res.length; k++){
            for(int i=0; i<orig.length(); i++){
                double cert = orig.emissions[i].getUncertainty() * pred.emissions[i].getUncertainty();
                if(cert < thresh) continue;
                try{
                    ComparableArray compA = (ComparableArray)orig.getElement(i);
                    ComparableArray compB = (ComparableArray)pred.getElement(i);
                    if(matchAtHigherLevel(compA   , compB)){
                       if( wobbleRoom(compA)>1){
	                        compare(compA,  compB,  fromTo);
	                        errorcase[fromTo[0]][fromTo[1]]++;
	                        int bin  = (int) Math.floor(cert / incr);
	                        error[bin][fromTo[0]==fromTo[1] ? 0 :1] ++;
	                        if(log!=null &&  fromTo[0]!=fromTo[1]){
	                           log.println(" mismatch  "+orig.getName()+" "+i+" "+compA+" "+compB);
	                        }
                       }
                    }
                }
                catch(Exception exc){
                    exc.printStackTrace();
                    System.err.println(getClass());
                    System.err.println(fromTo[0]+" "+fromTo[1]);
                    System.err.println(errorcase.length+" "+errorcase[0].length);
                    System.err.println("problem at "+i);
                    System.err.println(orig);
                    System.err.println(pred);
                    System.exit(0);
                }
            }
        //}
    }
    /** number of possibilities, given that we have fixed at higher level */
   protected abstract int wobbleRoom(ComparableArray compA) ;


    protected final void  conditionalCompare(ComparableArray orig, ComparableArray pred,  int[] fromTo) {
      if(parent!=null) parent.compare(orig, pred, fromTo);
      else{
          fromTo[0] = 0;
          fromTo[1] = 0;
      }
    }
    // Higher level could be copy number
    private boolean matchAtHigherLevel(ComparableArray orig,ComparableArray pred) {
        int[] fromTo = new int[] {0,0};
        conditionalCompare(orig, pred, fromTo);
        return fromTo[0]==fromTo[1];
        // return orig.noCopies(true)== pred.noCopies(true);
    }
protected abstract void compare(ComparableArray orig, ComparableArray pred,  int[] fromTo);
  public abstract String getString(int j);
  
    public void print(PrintWriter pw) {
      //  for(int kk=0; kk<this.errorcase.length; kk++){
        pw.println("Classifiction for "+this.getClass().getName()+" ");
       pw.println("rows are from "+original_name+"; columns are from "+this.inferred_name);
            double[]  origEqMod = new double[] {0,0};//  [correct, false]
            //pw.println(this.emSt);
            int len = errorcase.length+1;
            String[] header = new String[len];
            for(int j=0; j<errorcase.length; j++){
                header[j] = getString(j);
            }
            header[len-1] = "sum";
            pw.println(String.format("%-10s ", new String[] {""})+String.format(formatStr1.toString(), header));
            Integer[] sum = new Integer[len];
            Arrays.fill(sum, 0);
            for(int j=0; j<this.errorcase.length; j++){
                
                Double[] d = new Double[len];
               Arrays.fill(d, 0.0);
                for(int k=0; k<this.errorcase[j].length; k++){
                    d[k] = (double)this.errorcase[j][k];
                    sum[k]+=d[k].intValue();
                    d[len-1]+=d[k];
                    origEqMod[k==j ? 0 :1]+=d[k];
                }
                sum[len-1]+=d[len-1].intValue();
              //  pw.print("\t");
                //System.err.println(Arrays.asList(d));
                if(d[len-1]>0) pw.println(
                        String.format("%-10s ", new String[] {getString(j)})+String.format(formatStr.toString(), d));
               
            }
         //   pw.print("\t");
            pw.println(String.format("%-10s ", new String[]{"sum"})+String.format(formatStr.toString(), sum));
            pw.println(origEqMod[1]+" "+(origEqMod[1]+origEqMod[0])+" "+(origEqMod[1]/(origEqMod[1]+origEqMod[0])));
            for(int i=0; i<error.length; i++){
                int[] origEqMo = error[i];
                double bot =i*incr;
                double top = (double)(i+1)*incr;
                Object[]res = new Object[] {bot, top, origEqMo[1], origEqMo[1]+origEqMo[0],((double)origEqMo[1]/(double)(origEqMo[1]+origEqMo[0]))};;
                if(((Number)res[3]).doubleValue()>0 ){
                    pw.println(String.format("%5.1f - %5.1f  %5.3f %5.3f %5.3f", res));
                }
               // pw.println(bot+"-"+top+" "+origEqMo[1]+" "+(origEqMo[1]+origEqMo[0])+" "+((double)origEqMo[1]/(double)(origEqMo[1]+origEqMo[0])));
            }
       // }
        
    }

}