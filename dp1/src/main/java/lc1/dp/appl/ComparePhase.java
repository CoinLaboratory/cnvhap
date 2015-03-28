package lc1.dp.appl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.util.Constants;

public class ComparePhase {
	
	public static double thresh = 0.9;
    static Logger logger = Logger.getAnonymousLogger();
    DataCollection inferred;
    DataCollection original;
    String name_orig, name_inferred;
    public ComparePhase(DataCollection inferred, DataCollection original, String name_inferred, String name_origin){
        this.inferred=  inferred;
        this.original=  original;
        this.name_inferred = name_inferred;
        this.name_orig = name_origin;
       
    }
 
    private static void add(int[] d, int[]toAdd){
        for(int i=0; i<d.length; i++){
            d[i]+=toAdd[i];
        }
    }
    public static Set getSet(Integer[] i){
        return new HashSet<Integer>(Arrays.asList(i));
    }
    
    
    public void comparePhase(PrintWriter pw, File logDir, List<Integer> loc){
    	if(inferred==original) throw new RuntimeException("!!");
        if(inferred.length()!=original.length()) throw new RuntimeException(inferred.length+" !! "+original.length());
        logger.info("comparing ...");
        Set<Integer>[][] noCopies = new Set[][] {
            new Set[]{ getSet(new Integer[] {0,1,2,3,4,5,6}), getSet(new Integer[] {0,1,2,3,4,5,6})},
           new Set[]{ getSet(new Integer[] {1}), getSet(new Integer[] {1})},
           new Set[]{ getSet(new Integer[] {1}), getSet(new Integer[] {2})},
           new Set[]{ getSet(new Integer[] {1}), getSet(new Integer[] {3})},
           new Set[]{ getSet(new Integer[] {2}), getSet(new Integer[] {1})},
          new Set[]{ getSet(new Integer[] {2}),getSet( new Integer[] {2})},
            new Set[]{ getSet(new Integer[] {2}), getSet(new Integer[] {3})},
            new Set[]{ getSet(new Integer[] {3}), getSet(new Integer[] {1})},
            new Set[]{ getSet(new Integer[] {3}), getSet(new Integer[] {2})},
            new Set[]{ getSet(new Integer[] {3}), getSet(new Integer[] {3})}
        };
        try{
        PrintWriter[] out = new PrintWriter[noCopies.length];
        for(int i=0; i<noCopies.length; i++){
        	out[i] = new PrintWriter(new BufferedWriter(new FileWriter(new File(logDir, Arrays.asList(noCopies[i][0])+"->"+Arrays.asList(noCopies[i][1])))));
        }
        int[][] switch_ = new int[noCopies.length][];
        for(int i=0; i<switch_.length; i++){
            switch_[i] = new int[] {0,0};
            }
        for(Iterator<String>it = original.getKeys().iterator(); it.hasNext();){
            String key = it.next();
            PhasedDataState inference =( PhasedDataState) inferred.get(key);
            PhasedDataState original_no_missing =( PhasedDataState) original.get(inference.getName());
           double[]  cert = inferred.uncertaintyPhase.get(key);
            if(inference.length()!=original_no_missing.length()) {
            throw new RuntimeException("!!");
            }
          
          
            for(int i=0; i<switch_.length; i++){
                add(switch_[i], inference.phaseCorrect( original_no_missing,noCopies[i][0], noCopies[i][1], out[i], loc, cert, thresh));
                }
        }
        {
            pw.println("phasing accuracy");
            for(int i=0; i<switch_.length; i++){
            pw.println(Arrays.asList(noCopies[i][0])+"->"+Arrays.asList(noCopies[i][1]));
            int[] sw = switch_[i];
            Integer[] res = new Integer[] {sw[0], sw[1]};
            pw.println("perc "+String.format("%5i : %5i", res));
            pw.println("perc "+String.format("%5.3g ", new Double[] {(double)res[0]/(double)res[1]}));
            }
            }
        
        for(int i=0; i<noCopies.length; i++){
        	out[i].close();
        }
        }catch(Exception exc){
        	exc.printStackTrace();
        }
    }
 /*   public void compare(
           // String origin,
            double thresh,
            PrintWriter pw,
         // int[][] sources, 
           PrintWriter[] logSumm
           ) {
        EmissionStateSpace[] numberofcopies = original.getNoCopies();
        for(int ii1 =0; ii1<numberofcopies.length; ii1++){
        pw.println("NO COPIES == "+(ii1+1));
        
        
        int[] missing = new int[]{0,0,0, 0,0,0};
        EmissionStateSpace emStSp = numberofcopies[ii1];
        if(emStSp==null) continue;
        ErrorClassificationAbstract[] classif = new ErrorClassificationAbstract[]{
          new ErrorClassificationCNV(emStSp, thresh,  logSumm[0]), new ErrorClassificationGenotype(emStSp, thresh,   logSumm[1]),
          new ErrorClassificationHaplopair(emStSp, thresh,  logSumm[2])};
        for(int k=0; k<classif.length; k++){
        	classif[k].setNames(name_orig, name_inferred);
        }
        for(Iterator<String>it = original.getKeys().iterator(); it.hasNext();){
        String key = it.next();
        HaplotypeEmissionState inference = (HaplotypeEmissionState) inferred.dataL.get(key);
        if(inference==null) continue;
        HaplotypeEmissionState original_no_missing = (HaplotypeEmissionState)original.dataL.get(inference.getName());
        if(original==null) continue;
        if(!(inference.noCop()==ii1+1)) continue;
        if(!original.containsKey(inference.getName())) continue;
        
        // double[] ld1 = data.dataL.get(key).getEmiss(3);
        //   double[] ld2 = orig_no_missing.dataL.get(key).getEmiss(3);
        //     double[] nxt = unc==null || unc.size()==0  ? null : unc.get(key);
        //  PIGData original_with_missing = orig_with_missing.get(inference.getName());
        
        if(Constants.fillGaps()){
        
        for(int ik=0; ik<classif.length; ik++){
          classif[ik].compare(original_no_missing,   inference);// inferredPositions.get(original_no_missing.getName()));
        }
        // add(missing, inference.genotypeCorrect( original_no_missing, original_no_missing, nxt, thresh));
        }
        
        
        
        }
       
        if(Constants.fillGaps()){
        pw.println("missing data accuracy");
        pw.println("'-' called as '?' : total '-' //" +
         "'?' called as '-' : total '?' //" +
         " inferred ? with correct genotype  : total real '?'");
        Integer[] res = new Integer[] { missing[0], missing[1], missing[2], missing[3], missing[4], missing[5]};
        pw.println("perc "+String.format("%5i:%5i // %5i:%5i // %5i:%5i", res));
        for(int ik=0; ik<classif.length; ik++){
        classif[ik].print(pw);
        }
        }
}
pw.flush();
logger.info("... done comparing");
}*/
}
