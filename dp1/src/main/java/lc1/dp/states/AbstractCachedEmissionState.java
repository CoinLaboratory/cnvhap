package lc1.dp.states;


public abstract class AbstractCachedEmissionState extends CompoundState{
    public AbstractCachedEmissionState(CompoundState state) {
        super(state.getName()+"c", state.adv);
    }
    public abstract void transferCountsToMemberStates();
    public abstract void refreshSiteEmissions() ;
}
