 package lc1.dp.model;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.data.collection.ParseTree;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.SimpleEmissionStateSpace;
import lc1.dp.states.CachedEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.State;
import lc1.dp.states.WrappedEmissionState1;
import lc1.stats.IntegerDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Compressor;
import lc1.util.Constants;

public abstract class HaplotypeHMM extends MarkovModel implements Comparable {
	// final EmissionStateSpace stateEmissionStateSpace;

	// public static final int nullIndex = Emiss.indexOf(Emiss.N);
	/*
	 * public List<Integer[]> sameStates(){ List<Integer[]> res = new
	 * ArrayList<Integer[]>(); for(int j=1; j<this.modelLength(); j++){
	 * HaplotypeEmissionState hes = (HaplotypeEmissionState)this.getState(j);
	 * for(int j1=j+1; j1<this.modelLength(); j1++){ HaplotypeEmissionState hes1
	 * = (HaplotypeEmissionState)this.getState(j1); if(hes.same(hes1,0.1)){
	 * res.add(new Integer[] {j, j1}); } } } return res; } public
	 * EmissionStateSpace getStateEmissionStateSpace(){ return
	 * this.emissionStateSpace; }
	 */
	public boolean trainEmissions() {
		return this.pseudocountWeights[0][0] < 1e4;
	}

	public void modify(ZipFile hmmFile, int[] nocop1, List<Integer> locs)
			throws Exception {
		List<Integer> noCop = new ArrayList<Integer>();
		for (int i = 0; i < nocop1.length; i++) {
			noCop.add(nocop1[i]);
		}
		List<String> snps = Compressor.getIndiv(hmmFile, "SNPS", 3);
		List<String> locs1 = Compressor.getIndiv(hmmFile, "SNPS", 1);
		List<String> vals = Compressor.getIndiv(hmmFile, "Samples", 0);
		Map<Integer, List<Integer>> avail = new HashMap<Integer, List<Integer>>();
		for (int i = 0; i < vals.size(); i++) {
			Integer v = Integer.parseInt(vals.get(i));
			if (noCop.contains(v)) {
				List<Integer> avail_n = avail.get(v);
				if (avail_n == null) {
					avail.put(v, avail_n = new ArrayList<Integer>());
				}
				avail_n.add(i);
			}
		}
		int[] fromHMMToRepos = new int[states.size()];
		Arrays.fill(fromHMMToRepos, -1);
		for (int j = 1; j < states.size(); j++) {
			EmissionState st = (EmissionState) states.get(j);
			List<Integer> avail_n = avail.get(st.noCop());
			int sze = avail_n.size();
			if (sze > 0) {
				fromHMMToRepos[j] = avail_n.remove(sze - 1);
			}
		}

		inner: for (int i = 0; i < locs.size(); i++) {
			int ind1 = locs1.indexOf(locs.get(i).toString());
			if (ind1 < 0)
				continue inner;
			List<String> vals1 = Compressor
					.readZipFrom(hmmFile, snps.get(ind1));
			for (int j = 1; j < this.states.size(); j++) {
				int index = fromHMMToRepos[j];
				if (index >= 0) {

					HaplotypeEmissionState st = (HaplotypeEmissionState) states
							.get(j);
					if (!(st.emissions[i] instanceof IntegerDistribution)) {
						String[] str = vals1.get(index).split("\t");
						int nocop = st.noCop();
						int[] indices = st.getEmissionStateSpace()
								.getGenoForCopyNo(nocop);
						double[] d = new double[st.getEmissionStateSpace()
								.size()];
						Arrays.fill(d, 0.0);
						for (int k = 0; k < indices.length; k++) {
							d[indices[k]] = Double.parseDouble(str[k]);
						}
						double[] d1 = ((SimpleExtendedDistribution) st.emissions[i]).probs;
						Constants.normalise(d);
						st.emissions[i] = new SimpleExtendedDistribution1(d,
								Double.POSITIVE_INFINITY);
					}
				}

			}
		}
	}

	/**
	 * first is emissions, second is transitions, third is exponential trans,
	 * fourth is data
	 */
	protected double[][] pseudocountWeights = null;

	@Override
	public void setPseudoCountWeights(double[][] d) {
		this.pseudocountWeights = d;

	}

	public double logp = Double.NEGATIVE_INFINITY; // this is the score

	/** higher logp is less */
	public int compareTo(Object obj) {
		HaplotypeHMM obj1 = (HaplotypeHMM) obj;
		if (obj1.logp == logp)
			return 0;
		else if (logp < obj1.logp)
			return 1;
		else
			return -1;
	}

	public HaplotypeHMM(MarkovModel mm) {
		super(mm.getName() + "x" + mm.getName(), mm.noSnps);
		// this.innerModel = innerModel;
		// trans = new FreeSiteTrans1(this.states, noSnps,
		// FreeTransitionProbs1.class);

		// fixed =
		for (int j = 1; j < mm.modelLength(); j++) {
			EmissionState newState;
			EmissionState state_j = (EmissionState) mm.getState(j);

			newState = (EmissionState) state_j.clone();

			this.addState(newState);

		}
		in = mm.statesIn(1, 1);
	}

	public HaplotypeHMM(HaplotypeHMM hmm) {
		super(hmm);
		// this.stateEmissionStateSpace = hmm.stateEmissionStateSpace;
		in = hmm.in;

	}

	/*
	 * public HaplotypeHMM(CachedHMM hmm){ super(hmm); //
	 * this.stateEmissionStateSpace = hmm.getStateEmissionStateSpace(); in =
	 * hmm.statesIn(1, 1); if(Constants.CHECK){ for(int i=1; i<this.noSnps;
	 * i++){ for(int j=0; j<this.states.size(); j++){ int[] inij =
	 * hmm.statesIn(j, i); for(int k=0; k<inij.length; k++){ if(Math.abs(inij[k]
	 * - in[k])>0.01) throw new RuntimeException("!!"); } }
	 * 
	 * } }
	 * 
	 * }
	 */

	public abstract MarkovModel clone(boolean swtch);

	/*
	 * public HaplotypeHMM(List<HaplotypeEmissionState> hesList,
	 * EmissionStateSpace emStSp) { super("", hesList.get(0).noSnps() , emStSp);
	 * for(int j=0; j< hesList.size(); j++){ this.addState(hesList.get(j)); } in
	 * = new int[this.states.size()-1]; for(int jk=1; jk<states.size(); jk++){
	 * in[jk-1] = jk; }
	 * Logger.global.info("state space size is "+this.states.size
	 * ()+" for "+this.getName()); }
	 */

	/**
	 * num is the num in each state returns [group, index in group]
	 */
	public int[][] expand(int[] num) {
		int num1 = Constants.sum(num);
		List<State> st = new ArrayList<State>(this.states);
		this.states.clear();
		int[][] m = new int[st.size()][];
		char[] ch = new char[num1];
		char[] mod0 = Constants.modify(0);
		int k = 1;
		m[0] = new int[] { 0 };
		this.states.add(st.get(0));
		for (int i = 1; i < st.size(); i++) {

			m[i] = new int[num[i - 1]];
			for (int n = 0; n < num[i - 1]; n++) {
				ch[k - 1] = mod0[i - 1];
				HaplotypeEmissionState orig = (HaplotypeEmissionState) st
						.get(i);
				EmissionState sta = new HaplotypeEmissionState((orig), "" + k,
						Constants.expand_init_prior(0));
				this.addState(sta);
				m[i][n] = k;
				k++;

			}
		}
		in = new int[this.states.size() - 1];
		for (int jk = 1; jk < states.size(); jk++) {
			in[jk - 1] = jk;
		}
		Constants.modify0 = new char[Constants.inputDir.length][ch.length];
		Arrays.fill(Constants.modify0,ch);
//		System.arraycopy(src, srcPos, dest, destPos, k)ch;
		// HaplotypePanel1.calcStatEm();
		// ProbDists.ca_b = new ColorAdapter();
		// HMMPanel.ca = new ColorAdapter();
		// // HMMPanel.ca.getColor(0+"");
		// HMMPanel.ca.getColor(1+"");
		return m;
	}

	public HaplotypeHMM(String name, int numFounders, int noSnps,
			double[] init1, EmissionStateSpace emStSp1, char[][] codes,
			List<Integer> locs, Boolean[] probeOnly,
			ProbabilityDistribution[] numLevels) {
		super(name, noSnps);
		int maxcn = Constants.maxCopies();
		
		List<EmissionState>[] statesByCn = new List[maxcn + 1];
		for (int i = 0; i < statesByCn.length; i++) {
			statesByCn[i] = new ArrayList<EmissionState>();
		}
		double thresh = 1e-8;
		char[] codesSum ;
		if (codes == null) {
			codes = new char[1][numFounders];
			Arrays.fill(codes[0], 'a');
			codesSum = codes[0];
		}else{
			codesSum = new char[numFounders];
			boolean hasA = false;  
			for(int k=0; k<numFounders; k++){
				codesSum[k] = codes[0][k];
				for(int j=1; j<codes.length; j++){
					if(codes[j][k]!=codes[0][k]){
						codesSum[k] = 'a';
						hasA = true;
					}
				}
			}
			if(hasA) Arrays.fill(codesSum,'a'); //if it has one a it must have all a otherwise the pairing goes wrong, the emssion state space is messed up - cant find the right paired index
		}
		if (codes[0].length != numFounders)
			throw new RuntimeException("!!");

		for (int j = 0; j < codesSum.length; j++) {
			Integer i = codesSum[j]=='a' ? null : Integer.parseInt(codesSum[j] + "");
			List<EmissionState> statesByCn_i = i==null ? null : statesByCn[i];

			EmissionStateSpace emStSp = i==null ?  Emiss.mergedSpace:
			Emiss.spaceByCN[i.intValue()];
			double[][] probs = new double[codes.length][];
			for(int k=0; k<codes.length; k++){
				probs[k] = emStSp.getArray(codes[k][j] + "");
			}
			double[] probsSum = emStSp.getArray(codesSum[j]+"");
				

			EmissionState sta = null;
			if(i==null){
				sta = (EmissionState) makeState((j + 1) + "",
						noSnps, probs, Math.abs(probsSum[Constants
								.getMax(probsSum)] - 1.0) < 0.00001,
						numLevels, emStSp);
				for (int ik = 0; ik < probeOnly.length; ik++) {
					HaplotypeEmissionState st1 = (HaplotypeEmissionState) sta;
					PseudoDistribution dist1 = st1.emissions[ik];
					if (probeOnly[ik] != null && probeOnly[ik]
							&& dist1 instanceof SimpleExtendedDistribution) {
						SimpleExtendedDistribution dist2 = (SimpleExtendedDistribution) dist1;
						int nocop = emStSp.getCN(Constants.getMax(dist2.probs));
					     boolean sameCop = true;
					     for(int j1=0;sameCop && j1<dist2.probs.length; j1++){
					    	 if(dist2.probs[j1]>0 && emStSp.getCN(j1)!=nocop){
					    		 sameCop = false;
					    	 }
					     }
					     if(sameCop){
						   st1.emissions[ik] = new IntegerDistribution(emStSp
								.getByAlias(nocop, 0), nocop, 0);
					     }
					}
				}
			}
			else if ((i == 0 || i == 1)) {
				statesByCn_i
						.add(sta = (EmissionState) makeState((j + 1) + "",
								noSnps, probs, Math.abs(probsSum[Constants
										.getMax(probsSum)] - 1.0) < 0.00001,
								numLevels, emStSp));
				// if(sta instanceof HaplotypeEmissionState){
				for (int ik = 0; ik < probeOnly.length; ik++) {
					HaplotypeEmissionState st1 = (HaplotypeEmissionState) sta;
					PseudoDistribution dist1 = st1.emissions[ik];
					if (probeOnly[ik] != null && probeOnly[ik]
							&& dist1 instanceof SimpleExtendedDistribution) {
						SimpleExtendedDistribution dist2 = (SimpleExtendedDistribution) dist1;
						int nocop = emStSp.getCN(Constants.getMax(dist2.probs));
					
						st1.emissions[ik] = new IntegerDistribution(emStSp
								.getByAlias(nocop, 0), nocop, 0);
					}
				}
				// }
			} else {
				int sze = statesByCn_i.size();
				EmissionState[] states = new EmissionState[i];
				Arrays.fill(states, statesByCn[1].get(sze));
				CompoundEmissionStateSpace emstsp = (CompoundEmissionStateSpace) Emiss.spaceByCN[i];
				statesByCn_i.add(sta = new CachedEmissionState(
						new PairEmissionState(Arrays.asList(states), emstsp,
								true), emstsp.size()));
			}

			this.addState(sta);
		}
		in = new int[this.states.size() - 1];
		for (int jk = 1; jk < states.size(); jk++) {
			in[jk - 1] = jk;
		}
		Logger.global.info("state space size is " + this.states.size()
				+ " for " + this.getName());
	}

	public HaplotypeHMM(String name, int numFounders, int noSnps,
			HaplotypeEmissionState orig) {
		super(name, noSnps);

		for (int j = 0; j < numFounders; j++) {
			HaplotypeEmissionState sta = new HaplotypeEmissionState((orig), ""
					+ (j + 1), Constants.expand_init_prior(0));
			this.addState(sta);
		}
		in = new int[this.states.size() - 1];
		for (int jk = 1; jk < states.size(); jk++) {
			in[jk - 1] = jk;
		}
		Logger.global.info("state space size is " + this.states.size()
				+ " for " + this.getName());
	}

	/*public HaplotypeHMM(DataCollection datac,List<String> names) {
		this(datac.name,getSub( datac.dataL, names), datac.loc.size());
		
	}*/
	private static Map<String, List<HaplotypeEmissionState>> getSub(DataCollection datC,
			List<String> names) {
		Map<String, EmissionState> dataL = datC.dataL;
		List<String> indiv = datC instanceof MergedDataCollection ? ((MergedDataCollection)datC).ldl[0].indiv() : datC.indiv();
		ParseTree tree = datC instanceof MergedDataCollection ? ((MergedDataCollection)datC).ldl[0].tree : datC.tree;
		Map<String, List<HaplotypeEmissionState>> sub = new HashMap<String, List<HaplotypeEmissionState>>();
		
		for(int i=0; i<names.size(); i++){
			String nme = names.get(i);
			List<HaplotypeEmissionState> l;
		
			l = new ArrayList<HaplotypeEmissionState>();
			
			for(Iterator<String> it = indiv.iterator(); it.hasNext();){
				String ke = it.next();
				if(//ke.startsWith(nme) || tree!=null &&
						tree.isAncestral(nme, ke)){
					System.err.println(nme+"->"+ke);
					l.add((HaplotypeEmissionState)dataL.get(ke));
				}
			}
			sub.put(nme, l);
		}
		return sub;
		
	}

	

	public HaplotypeHMM(String name, List<State> subList, int noSnps) {
		super(name, noSnps);

		for (int j = 0; j < subList.size(); j++) {

			this.addState(subList.get(j));
		}
		in = new int[this.states.size() - 1];
		for (int jk = 1; jk < states.size(); jk++) {
			in[jk - 1] = jk;
		}
		Logger.global.info("state space size is " + this.states.size()
				+ " for " + this.getName());
	}

	public HaplotypeHMM(String name, int noSnps) {
		super(name, noSnps);

	}

	/** distributes prob at each site among same copy number */
	private static void distributeAmongSame(double[] probs,
			EmissionStateSpace emStSp) {
		double[] copyCount = new double[3];
		int[] count = new int[3];
		Arrays.fill(copyCount, 0.0);
		Arrays.fill(count, 0);
		for (int i = 0; i < probs.length; i++) {
			int cp = ((Emiss) emStSp.get(i)).noCopies();
			copyCount[cp] += probs[i];
			count[cp]++;
		}
		Arrays.fill(probs, 0);
		double sum = 0;//
		for (int i = 0; i < probs.length; i++) {
			int cp = ((Emiss) emStSp.get(i)).noCopies();
			double val = copyCount[cp] / (double) count[cp];
			probs[i] = val;
			sum += val;
		}
		if (Math.abs(sum - 1.0) > 0.0001) {
			throw new RuntimeException("!!");
		}
	}

	/*
	 * public HaplotypeHMM(String name, int numFounders, int noSnps, double[]
	 * init,EmissionStateSpace stateSpace, // EmissionState[]
	 * emissionStateSpaceDist, boolean allowedScore, double[][] meanvarskew){
	 * super(name, noSnps, stateSpace); // this.stateEmissionStateSpace =
	 * Emiss.getStateEmissionStateSpace(numFounders); this.allowedScore =
	 * allowedScore; double[] r_prir = Constants.r_prior; //
	 * this.emissionStateSpaceDist = emissionStateSpaceDist; for(int j=0; j<
	 * numFounders; j++){
	 * 
	 * this.addState(makeState((j+1)+"", noSnps, init,
	 * Math.abs(init[Constants.getMax(init)]-1.0 )< 1e-10,
	 * Constants.meanvarskew(j), r_prir[j])); // } }
	 * 
	 * in = new int[this.states.size()-1]; for(int jk=1; jk<states.size();
	 * jk++){ in[jk-1] = jk; }
	 * Logger.global.info("state space size is "+this.states
	 * .size()+" for "+this.getName()); }
	 */

	boolean allowedScore = true;

	public EmissionState makeState(String st, int noSnps, double[][] init,
			boolean fixed, ProbabilityDistribution[] numLevels,

			EmissionStateSpace emissionStateSpace) {
		//
		int[] index = new int[init.length];
		Integer[] cn = new Integer[init.length];
		for(int j=0; j<init.length; j++){
			index[j] = Constants.getMax(init[j]);
		    cn[j] = emissionStateSpace.getCN(index[j]);
		  inner: for(int k=0; k<init[j].length; k++){
			if(init[j][k]>0 && emissionStateSpace.getCN(k)!=cn[j].intValue()){
				cn[j] = null;
				break inner;
			}
		  }
		}
		return new HaplotypeEmissionState(st, noSnps, Constants.u_global(0)[0],
				init, emissionStateSpace, cn,
				numLevels);

	}

	// abstract public void setNewTransitionsRandom( int i, double small);

	// abstract public boolean increaseStates(double u, double small, int
	// numToAdd, double[] init);

	@Override
	public void initialiseCounts() {
		super.initialiseCounts();
		// if(nullEm!=null) nullEm.initialise();

	}

	/*
	 * @Override public void updateEmissionStateSpaceDist(int stage){
	 * 
	 * if(Constants.inputFile().endsWith("lhood") ){//||
	 * !Constants.modelInaccurateData(stage)){ this.emissionStateSpaceDist=null;
	 * return; } if(emissionStateSpaceDist==null )// ||
	 * !Constants.modelInaccurateData(stage)) return; emissionStateSpaceDist=new
	 * EmissionState[emissionStateSpace.size()]; for(int i=0;
	 * i<emissionStateSpace.size(); i++){ double[] res = new
	 * double[emissionStateSpace.size()]; if(emissionStateSpace.get(i)
	 * ==Emiss.N){ double[] probGapIsNotGap = Constants.probGapIsNotGap();
	 * if(stage >=probGapIsNotGap.length) continue; Arrays.fill(res,
	 * Math.abs(probGapIsNotGap[stage])/ ((double)res.length-1.0)); res[i] = 1.0
	 * -Math.abs(probGapIsNotGap[stage]); emissionStateSpaceDist[i] = new
	 * HaplotypeEmissionState("b"+i,Constants.siteSpecificBground() ? noSnps :
	 * 1,res); if(probGapIsNotGap[stage]<0){ ((HaplotypeEmissionState)
	 * emissionStateSpaceDist[i]).setFixed(true); } } else{ double[]
	 * probNonGapIsGap = Constants.probNonGapIsGap(); if(stage
	 * >=probNonGapIsGap.length) continue; Arrays.fill(res,
	 * Math.abs(Constants.modelError()[stage]) / ((double)res.length-2.0));
	 * res[this.stateSpaceToIndex.get(Emiss.N)] =
	 * Math.abs(probNonGapIsGap[stage]);
	 * 
	 * res[i] = 1.0 - Math.abs(probNonGapIsGap[stage])-
	 * Math.abs(Constants.modelError()[stage]); emissionStateSpaceDist[i] = new
	 * HaplotypeEmissionState("b"+i,Constants.siteSpecificBground() ? noSnps :
	 * 1,res); if(probNonGapIsGap[stage]<0 || Constants.modelError()[stage]<0){
	 * ((HaplotypeEmissionState) emissionStateSpaceDist[i]).setFixed(true); }
	 * //probError = Constants.modelError(); }
	 * 
	 * } }
	 */

	@Override
	public List<State> statesOut(State st, int beforeToEmission) {
		if (beforeToEmission == this.noSnps - 1)
			return states.subList(0, 2);
		else
			return this.states.subList(1, states.size());

	}

	@Override
	public List<State> statesIn(State st, int beforeToEmission) {
		if (beforeToEmission == 0)
			return states.subList(0, 1);
		else
			return this.states.subList(1, states.size());

	}

	public abstract void transferTransitionCountsToProbs(int index);

	public void transferCountsToProbs(int index) {
		if(Constants.trainEmissions()){
		this.transferEmissionCountsToProbs(index);
		}
		/*
		 * if(emissionStateSpaceDist!=null && Constants.modelInaccurateData(0)){
		 * for(int i=0; i<this.emissionStateSpace.size(); i++){
		 * HaplotypeEmissionState hes =
		 * (HaplotypeEmissionState)this.getEmissionStateSpaceDistribution(i);
		 * if(hes!=null)
		 * (hes).transferCountsToProbs(this.pseudocountWeights[1]); } }
		 */
		if (Constants.trainTransitions())
			this.transferTransitionCountsToProbs(index);
		// if(nullEm!=null)
		// nullEm.transferCountsToProbs(this.pseudocountWeights);
	}

	public void transferEmissionCountsToProbs(int index) {
		if (pseudocountWeights[0][0] < 1000 || true) {
			for (int j = this.modelLength() - 1; j > 0; j--) { // backwards so
																// we transfer
																// from paired
																// states first
				EmissionState hes = (EmissionState) getState(j);
				if ((Constants.modifyDirectCounts[0] != 1.0 || Constants.modifyDirectCounts[1] != 1.0)
						&& hes instanceof WrappedEmissionState1) {
					((WrappedEmissionState1) hes).modifyDirectCounts(Constants
							.modifyDirectCounts());
				}
				// System.err.println(j+" "+Constants.print(((HaplotypeEmissionState)hes).emissions[0].counts()));
				// if(Constants.stepwise() <= 0 || index <=
				// Constants.stepwise()*(j-1))
				hes.transferCountsToProbs(this.pseudocountWeights[0][0]);
				
			}
			 for (int j = this.modelLength() - 1; j > 0; j--) { // backwards so
					// we transfer
					// from paired
					// states first
				EmissionState hes = (EmissionState) getState(j);
				hes.refreshSiteEmissions();
				}
			  
			
		}
	}

	public static Double[] subList(Double[] l, List<Integer> cols) {
		if (cols == null)
			return l;
		Double[] res = new Double[cols.size()];
		for (int i = 0; i < cols.size(); i++) {
			res[i] = l[cols.get(i)];
		}
		return res;
	}

	public void print(PrintWriter pw, List<Integer> cols, int popsize) {
		/*StringBuffer sb1 = new StringBuffer();
		for (int i = 0; i < this.modelLength(); i++) {
			sb1.append("%8.2g ");
		}*/
	//	super.print(pw, cols, popsize);
		// if(nullEm!=null) //pw.println(+this.nullEm.probs[0]);
		// this.nullEm.print( pw, "bground ");
	//	pw.println("hitting probs");
		List<Integer> pos = DataCollection.datC.loc;
		List<String>snpid = DataCollection.datC.snpid;
		StringBuffer sb = new StringBuffer("snpid\tpos");
		Emiss.EmissionStateSpaceForNoAlleles ems = Emiss.emiss;
		EmissionStateSpace emstsp1 = ems.mergedSpace;
		//int i1 = emstsp1.get(Emiss.b());
		List<Comparable> list = Arrays.asList(new Comparable[] {Emiss.b()});
		//EmissionStateSpace emstsp1 = ((CompoundEmissionStateSpace)Emiss.emiss.mergedSpace).getMembers()[0];
		for(int i1=0; i1<list.size(); i1++){
			Comparable compa = list.get(i1);
			sb.append("\t");
			sb.append(proc(emstsp1.getHaploString(compa)));
		}
		for(int j=1; j<this.modelLength(); j++){
			 sb.append("\t");
			 sb.append(this.getState(j).getName());
		}
		for(int j=1; j<this.modelLength(); j++){
			sb.append("\t");
			 sb.append(this.getState(j).getName());
		}
		pw.println(sb.toString());
		
		double[][] hp = new double[this.noSnps][this.modelLength()];
		this.getHittingProb(this.noSnps,hp);
		double[] pr = new double[this.states.size()-1];
		PseudoDistribution[] dist = new PseudoDistribution[states.size()-1];
		int[] cn = new int[states.size()-1];
		EmissionStateSpace[] stsp1 = new EmissionStateSpace[states.size()-1];
		EmissionState[] states_ = new EmissionState[states.size()-1];
		for(int i=0; i<this.noSnps; i++){
		//	boolean break_ = snpid.get(i).startsWith("640");
			sb = new StringBuffer(snpid.get(i)+"\t"+pos.get(i));
			{
				for(int i1=0; i1<dist.length; i1++){
					EmissionState state = ((EmissionState)states.get(i1+1));
				
					if(state instanceof WrappedEmissionState1){
						state = ((WrappedEmissionState1)state).inner;
					}
					if(state instanceof HaplotypeEmissionState){
						dist[i1] = ((HaplotypeEmissionState)state).emissions[i];
					}
					else{
						dist [i1]= ((CachedEmissionState)state).emissions[i];
					}
					cn[i1] = state.noCop();
					states_[i1] = state;
					stsp1[i1] = state.getEmissionStateSpace();
				}
				for(int i2=0; i2<list.size(); i2++){
					double sum=0;
					sb.append("\t");
					Comparable compa = list.get(i2);
					int i1 = emstsp1.get(compa);
					int cn1 = emstsp1.getCN(i1);
					
					for(int k=0; k<pr.length; k++){
						if(cn1!=cn[k]) pr[k]=0;
						else {
							int i3 = stsp1[k].get(compa);
							pr[k] = states_[k].score(i3, i)*hp[i][k+1];
						}
						sum+=pr[k];
					}
					for(int k=0; k<pr.length; k++){
						pr[k] = pr[k] /sum;
					}
					int max = Constants.getMax(pr);
					double maxd = pr[max];
					sb.append((max+1)+"="+String.format("%5.3g", maxd).trim());
					if(maxd<0.95){
						int max2 = Constants.getMax2(pr, max);
						double maxd2 = pr[max2];
						sb.append(";"+(max2+1)+"="+String.format("%5.3g", maxd2).trim());
					}
				}
			}
			
			int ii= i;
		
			for(int j=1; j<this.modelLength(); j++){
				 State state =this.getState(j);
				  sb.append("\t");
				 //if(state instanceof HaplotypeEmissionState){
				
					EmissionStateSpace emstsp = ((EmissionState)state).getEmissionStateSpace();
					Integer fixed = dist[j-1].fixedInteger();
					if(fixed!=null){
						sb.append(proc(emstsp.getHaploString(emstsp.get(fixed.intValue()))));
					}
					else{
						double[] probs = dist[j-1].probs();
						int max = Constants.getMax(probs);
						double maxd = probs[max];
						sb.append(proc(emstsp.getHaploString(emstsp.get(max)))+"="+String.format("%5.3g", maxd).trim());
						if(maxd<0.95){
							int max2 = Constants.getMax2(probs, max);
							double maxd2 = probs[max2];
							sb.append(";"+emstsp.get(max2).toString()+"="+String.format("%5.3g", maxd2).trim());
						}
					}
				 //}
			}
			for(int j=1; j<this.modelLength(); j++){
				sb.append("\t");
				sb.append(String.format("%5.2g", hp[ii][j]));
			}
			pw.println(sb.toString());
		}
		pw.flush();
		//pw.close()
		/*Double[] hp1 = new Double[hp[0].length];
		if (hp != null) {
			for (int i = 0; i < hp.length; i++) {
				for (int j = 0; j < hp[i].length; j++) {
					hp[i][j] = (double) Math.round(10000 * hp[i][j]) / 100;
				}
			}
		}
		pw.println("hitting probs ");
		for (int i = 0; i < this.noSnps; i++) {
			for (int k = 0; k < hp1.length; k++) {
				hp1[k] = hp[i][k];
			}
			if (cols == null || cols.contains(i)) {
				pw.println(i + " " + String.format(sb1.toString(), hp1));
			}
		}*/

	}

	private Object proc(String haploString) {
		if(haploString.length()==0) return "_";
		else return haploString;
	}

	/*
	 * public void fix() { for(int j=1; j<this.modelLength(); j++){
	 * HaplotypeEmissionState hes = (HaplotypeEmissionState) this.getState(j);
	 * hes.fix(); }
	 * 
	 * }
	 * 
	 * public void setTrain(boolean emiss, boolean gap) { for(int j=1;
	 * j<this.modelLength(); j++){ HaplotypeEmissionState hes =
	 * (HaplotypeEmissionState) this.getState(j); // hes.train_j = emiss; }
	 * 
	 * }
	 */

	public double getLogP() {
		return logp;
	}

	int[] in;// = new int[0];
	final int[] in0 = new int[] { 0 };

	public int[] statesIn(int j, int i) {
		if (i == 0 && j != 0)
			return in0;
		else
			return in;

	}

	public int[] statesOut(int j, int i) {
		if (i == noSnps - 1 && j != 0)
			return in0;
		else
			return in;
	}

	public void modify(DataCollection datC, List<String> datamodel) {
//		DataCollection datC = datC1 instanceof MergedDataCollection ? ((MergedDataCollection)datC1).ldl[0] : datC1;
		Map<String,List<HaplotypeEmissionState>> l = getSub(datC, datamodel);
		SimpleEmissionStateSpace	emstsp = (SimpleEmissionStateSpace)this.emstsp;
		//Iterator<Map.Entry<String, List<HaplotypeEmissionState>>> it = l.entrySet().iterator();
		for(int i=0; i<datamodel.size(); i++){
		
			String key = datamodel.get(i);
			 List<HaplotypeEmissionState> val = l.get(key);
			HaplotypeEmissionState st = (HaplotypeEmissionState)this.states.get(i+1);
			st.modify(val);
			//nxt.getValue()
			emstsp.setNme(i,key);
		}
	}

	/*
	 * public void transferProbToPseudo() { for(int j=1; j<this.modelLength();
	 * j++){ SimpleExtendedDistribution[] em =
	 * ((HaplotypeEmissionState)getState(j)).emissions; for(int i=0;
	 * i<em.length; i++){ em[i].transferProbToPseudo(); } }
	 * 
	 * }
	 */

}
