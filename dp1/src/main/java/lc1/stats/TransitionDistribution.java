package lc1.stats;


public interface TransitionDistribution {

    public abstract void transferProbToPseudo() ;

    public abstract void addCounts(StateDistribution distribution);

    public abstract double getProbs(int to);

    public abstract double getPseudo(int to);

    public abstract void initialise();

  //  public abstract void resample(double u, boolean restart);

    public abstract void transfer(double pseudoCountWeight);

    public abstract double KLDistance(TransitionDistribution d2);

    public abstract double sum();

   // public abstract double[] probs();

  //  public abstract double[] pseudo();

    public abstract void fillPseudo(double d);

    public abstract void setProbs(int to, double d);

    public abstract void setPseudo(int to, double d_ps);

    public abstract TransitionDistribution clone(TransitionDistribution distribution);

    public abstract void validate();

}
