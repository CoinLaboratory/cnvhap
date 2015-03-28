package lc1.dp.data.representation;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.EmissionState;
import lc1.stats.PseudoDistribution;

public interface CSOData extends SSOData, ScorableObject {


    public abstract int noCopies();

   // public abstract Comparable copyElement(Comparable element);

    public abstract void mkTrCompArray();

   // public abstract void mix();

    public abstract void set(int i, Comparable obj);

    // if elements are TrComparableArrays, it returns the index rather than the state!
    public abstract PIGData[] split();

    /**phase position i relative to i1.  Does not actually change data, but removes all elemenst 
     *  from set poss all those things not compatible with the most likely count of num differences
     * only conisders those elements which have an equal genotype
     *  i1 > i2
     *  phasing i1 relative to i2
     * */
 //   public abstract double samplePhase(EmissionStateSpace emStSp,
  //          List<CSOData> spList, int i1, int i2,
  //          Set<Comparable> poss);

    /** phase positions are relative to these genotypes */
    public abstract int[] phaseCorrect(CSOData original,
            Collection<Integer> noCopiesL, Collection<Integer> noCopiesR, PrintWriter pw, List<Integer>loc, double[] cert, double thresh);

    public abstract void samplePhase(List<CSOData> spList,
            double[] uncertainty, EmissionStateSpace emStSp);

    /** gets positions which need to be phased */
  //  public abstract Integer[] updatePhasePositions(EmissionStateSpace emStSp);


 //  public abstract double[] getUncertainty(EmissionState st);

   public abstract double getUncertainty(EmissionState st, int i);

    /** also updates emst 
    public abstract void sampleGenotype(HaplotypeEmissionState emst,
            List<CSOData> spList);*/

    public abstract void print(PrintWriter pw, boolean expand, boolean mark,
            Collection<Integer> toDrop,List<Character> alleleA, 
            List<Character> alleleB, PseudoDistribution[] unc);

    public abstract void print(PrintWriter pw, boolean idline, boolean expand,
            boolean mark, Collection<Integer> toDrop, List<Character> alleleA, 
            List<Character> alleleB, PseudoDistribution[] unc);

   // public abstract boolean equals(CSOData dat, int i);

  //  public abstract boolean isNested();


}