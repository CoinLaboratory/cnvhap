package lc1.dp.appl;

import java.util.Arrays;
import java.util.Collection;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LightWeightDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;

public class LightWeightDataCollection2 extends LightWeightDataCollection {

	LightWeightDataCollection[] dc;
	private int default_index;
	//List<String> indiv;
	public LightWeightDataCollection2(LightWeightDataCollection[] dcin1) {
		super(getName(dcin1), dcin1[0].midpoint);
		this.dc = dcin1;
		this.indiv = dcin1[0].indiv();
		EmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(2);
		for (int i = 0; i <indiv.size(); i++) {
			String key = indiv.get(i);
			HaplotypeEmissionState s1 =(HaplotypeEmissionState) dcin1[0].dataL.get(key);
			HaplotypeEmissionState value = createEmissionState(key, s1.noCop());
			value.emissions[0] = new SimpleExtendedDistribution(emstsp.defaultList.size());
			value.setNoCop(s1.noCop());
			dataL.put(key, value);
			data.put(key, SimpleScorableObject.make(key, 1, value
					.getEmissionStateSpace(), this.index));

		}
		this.default_index = emstsp.getByAlias(2, 0);
		this.fixed = new Integer[dcin1.length];
		this.dist = new PseudoDistribution[dcin1.length];
		/*for(int i=0; i<dist.length; i++){
			dist[i] = new SimpleExtendedDistribution(emstsp.defaultList.size());
		}*/
	}

	private static String getName(DataCollection[] dcin1) {
		StringBuffer sb = new StringBuffer(dcin1[0].name);
		for(int i=1; i<dcin1.length; i++){
			sb.append("_"+dcin1[i].name);
		}
		return sb.toString();
	}
	@Override
	public int restricToAlias(Collection<String> alias) {
		int sze=0;
		for(int i=0; i<dc.length; i++){
			sze+=dc[i].restricToAlias(alias);
		}
		return sze;
	}

	public boolean canProcess(String snp_id){
		for(int i=0; i<dc.length; i++){
			if(!dc[i].canProcess(snp_id)) return false;
		}
		return true;
	}
	
	PseudoDistribution[] dist;
	Integer[] fixed;
	
	
	public boolean updateIndex(String rsid) {
		boolean done = true;
		for(int i=0; i<this.dc.length; i++){
			done = done && this.dc[i].updateIndex(rsid);
		}
		if(done){
			for(int j=0; j<this.indiv.size(); j++){
				String key = indiv.get(j);
				HaplotypeEmissionState em = ((HaplotypeEmissionState)this.dataL.get(key));
				
				PseudoDistribution combined = em.emissions[last];
				double[] probs = combined.probs();
				
				Arrays.fill(probs,0.0);
				boolean allFixed = true;
				boolean allEqual = true;
				for(int i=0; i<dc.length; i++){
					dist[i] = ((HaplotypeEmissionState)dc[i].dataL.get(key)).emissions(last);
					fixed[i] = dist[i].fixedInteger();
					allFixed = allFixed && fixed[i]!=null;
					allEqual = allFixed && fixed[i].equals(fixed[0]);
				}
				if(allEqual){
					probs[fixed[0]] = 1.0;
				}
				else{
					EmissionStateSpace emstsp = em.getEmissionStateSpace();
					double sum=0;
					for(int k=0; k<probs.length; k++){
						
						probs[k] = getMin(dist,k);
						sum+=probs[k];
					}
					probs[default_index]+=1.0- sum;
				}
				
			}
		}
		
		
		return done;
	}
	
	

	private double getMin(PseudoDistribution[] dist2, int k) {
		double min = 1.0;
		for(int i=0; i<dist2.length; i++){
			double v = dist2[i].probs(k);
			if(v<min) min = v;
		}
		return min;
	}
}
