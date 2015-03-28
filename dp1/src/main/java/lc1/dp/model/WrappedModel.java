package lc1.dp.model;


public interface WrappedModel {
    public CompoundMarkovModel getHMM();
    public PairMarkovModel unwrapModel();
}
