package lc1.dp.states;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import lc1.CGH.Aberation;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;

public class BackgroundEmissionState extends HaplotypeEmissionState {

	
	
	public double[] cn_distribution;
	 public BackgroundEmissionState(String string, Integer length, int size,
			EmissionStateSpace emStSp, short s) {
		super(string, length, size, emStSp, null, null);
		this.data_index = s;
		this.cn_distribution = new double[emStSp.copyNumber.size()];
	}
	 
	 
	public BackgroundEmissionState(
			BackgroundEmissionState bes) {
		super(bes);
		this.cn_distribution = new double[emStSp.copyNumber.size()];
	}


	
	public Object clone(){
	       return new BackgroundEmissionState(this);
	   }
	@Override
	/** just updates background emission and does not update hmm state */
	public  void addCount(HaplotypeEmissionState data_state, Double double1, int i,double weight) {
	        //EmissionState data_state = this;
		HaplotypeEmissionState hmm_state = this;
if(double1>Constants.countThresh()){
		addCount(hmm_state.noCop(),  double1, i); //takes advantage of ordering of emission state space
}
	 }
	/*private double calcDistribution(int noCopies, int i) {
		Arrays.fill(this.distribution,0.0);
		double sum=0;
		for(int j=0; j<this.cn_distribution.length; j++){
			if(this.emStSp.getCN(j)==noCopies){
				distribution[j] = this.emissions[j].probs()[j]*this.emStSp.getWeight(j);
				sum+=distribution[j];
			}
		}
		return sum;
	}*/
	 public double calcDistribution(EmissionState hmm_state, int i, double[] distribution, double p){
	      /*  double sum=0;
	        for(Iterator<EmissionState> it = DataCollection.datC.dataL.values().iterator(); it.hasNext();){
	        	HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
	        if(hes==this) continue;
	        	sum+=	Math.log(hes.calcDistribution(hmm_state.noCop(), i, this.cn_distribution)
	        		);
	        }
	       // Arrays.fill(distribution,0.0);*/
	        double res = this.emissions[i].probs(hmm_state.noCop());//Math.exp(sum)*this.emissions[i].probs(hmm_state.noCop())*p;
	        distribution[hmm_state.noCop()] += res;
	        
	        return res;
	    }
	
	 @Override
	    public  double score(HaplotypeEmissionState data_state , int i_hmm) {
	    	return 1.0;//
	    /*    double sc =  0.0;
	        for(Iterator<EmissionState> it = DataCollection.datC.dataL.values().iterator(); it.hasNext();){
	        	HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
	        	
	        }
	    
	        return Math.exp(sc)*((SimpleExtendedDistribution)this.emissions[i_hmm]).probs[hmm_state.noCop()];//*this.emStSp.getWeight(hmm_state.noCop());
	    */
	    }
	 
	/* (non-Javadoc)
	    * @see lc1.dp.data.representation.PIGData#getDeletedPositions(lc1.dp.states.EmissionState)
	    */
	   public final Collection<Aberation> getDeletedPositions(List<Integer> geno, List<Integer> loc,  PrintWriter pw, int cn) {
	   	  List<Aberation> l = new ArrayList<Aberation>();
	   	 
	       return l;
	   }

	   
	public double score(EmissionState hmm_state, boolean logspace,
			HaplotypeEmissionState hes, double[] fg_cn, int i_hmm) {
		// TODO Auto-generated method stub
		if(hes==this) throw new RuntimeException("!!");
        //	double sc1 =0;
        	//
        	double sc=
        		hes.calcDistribution(hmm_state.noCop(), i_hmm, this.cn_distribution, fg_cn)
        		*   hmm_state.score(hmm_state.noCop(),i_hmm);
        	//	*this.calcDistribution(hmm_state, i_hmm)
        		;
        				
        	//	*this.emissions[i].sc		
        		
        	//sc1+=
        	//	this.emissions[i_hmm].scoreR(hmm_state.noCop(),fg,i_hmm);
		return sc;
	}

	public double[] calcCNDist(Integer noCop, int i_hmm, PseudoDistribution dist){
		 for(int j=0; j<emStSp.cnLength(); j++){
			  this.cn_distribution[j] = ( dist.scoreR(j, noCop,i_hmm))*this.emissions[i_hmm].probs(j);
			   //sc1 +=   d[j];
		   }
		 return this.cn_distribution;
	}
	public double scoreR(Integer noCop, int i_hmm, PseudoDistribution dist) {
		 //double[] d = new double[emStSp.cnLength()];
		this.calcCNDist(noCop, i_hmm, dist);
		return Constants.sum(this.cn_distribution);
	}

}
