	package lc1.dp.data.collection;
	
	import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.AbstractMultiDimMax;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.Multidimmax;
import lc1.stats.MultidimmaxSingle;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

import org.apache.commons.math.special.Gamma;

import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.Beta;
import cern.jet.random.engine.DRand;
	
	
	public class MatchedDistributionCollection extends
			AbstractDistributionCollection {
	
		
		final AbstractMultiDimMax mdf;
		
		int sample_index1 =0;
		@Override
		public void setIndiv(String protName) {
			this.sample_index1 = this.indiv.indexOf(protName);
			// TODO Auto-generated method stub
			
		}
		@Override
		public String getInfo(String key) {
			if(Constants.trainCellularity() >=1){
			  int indivk = indiv.indexOf(key);
				return " ##Estimated cellularity:"+cellularity[indivk];
			}else{
				return  "";
			}
		}
		
		static class BetaBinom {
			 double alpha, beta, trials;
			 public void set(double alpha, double beta, double trials){
				 this.alpha = alpha;
				 this.beta  = beta;
				 this.trials = trials;
			 }
			 public double getDensity(double k) {
	             //int k = (int) Math.rint(x);
	             
	 //    if (k < 0 | k > trials) return 0;
	     
	     return (Gamma.logGamma(k+alpha)+Gamma.logGamma(trials-k+beta)+Gamma.logGamma(alpha+beta)+Gamma.logGamma(trials+2)) - 
	             (Math.log(trials+1)+Gamma.logGamma(alpha+beta+trials)+Gamma.logGamma(alpha)+Gamma.logGamma(beta)+Gamma.logGamma(k+1)+Gamma.logGamma(trials-k+1));
	 }
	
			
		}
		//double cellularity;
		//double ratio;
		double average =0;
		double average_count=0;
		
		static final boolean ratioAsLevels = Constants.ratioAsLevels()!=null && Constants.ratioAsLevels;
		static final boolean cellAsLevels = Constants.ratioAsLevels()!=null && !Constants.ratioAsLevels;
		static final boolean trainPool = Constants.ratioAsLevels()!=null && Constants.trainPool();
		static final int trainCellularity = Constants.trainCellularity();
		static final double mix = 0.5; //essentially controls the starting bands for the ratios.  This is effectively the minimunm value
		
//		double[] vals; // cellularity, ratio
		
		final double[] cellularity, ratio, ratiosLe, ratiosR;

		public double ratios(double back) {
			if(Constants.trainPool())
			{
				if(back==Constants.maxPloidy1) return 1;
			
			if(back<Constants.maxPloidy1) return ratiosLe[(int) back];
			else return ratiosR[(int)back-Constants.maxPloidy1-1];
			}else{
				return (double) back/(double)Constants.maxPloidy1;
			}
		}
		
		private double totnormal;
		//private double prop = 1;
	//	private double tottumour;
		private double priorWeight = Constants.priorWeight();
		static int numIt = Constants.noSampleFromBeta();
		static double betaDownWeight = Constants.betaDownWeight();
		
		double noprobes;
		
		
		//double countProbesTumour=0;
		
		 Double[] refCount;  //this lists the tumour count in each window
		final public HaplotypeEmissionState pool;
		
		 Double[] backCNall;
		
		
		
		public double refCount(int i){
			return refCount[i];
		}
		public double backCN(int i){
			return backCNall[i];
		}
		
		public void drop(List<Integer> toDrop) {
			this.pool.removeAll(toDrop);
			refCount = (Double[]) Constants.drop(refCount, toDrop);
			backCNall = (Double[]) Constants.drop(backCNall, toDrop);
			
			}

		
		final short index;
		
	
		List<String> indiv;// = new ArrayList<String>();
		final int numsamples;
		final int maxCN;		
		
	  //  List<BinomialDistr1> dist = new ArrayList<BinomialDistr1>();
		public MatchedDistributionCollection(int index, 
			 double cellularity, File dir, int maxCN,int maxCN1, int noprobes,HaplotypeEmissionState ref, Double[] ratiosL, List<String> indiv) {
			// TODO Auto-generated constructor stub
			this.index = (short)index;
			mdf = indiv.size()>1 && Constants.trainCellularity()==1 ? new MultidimmaxSingle(this): new Multidimmax(this);
			
			this.indiv = new ArrayList<String>(indiv);
			this.indiv.remove("pool");
			this.indiv.remove(Constants.reference());
			this.numsamples = this.indiv.size();
			System.err.println("indiv: "+Arrays.asList(this.indiv));
			this.cellularity = mdf.add("cellularity", Constants.initialCellularity[0], 0.001, 1.0, numsamples, trainCellularity>=1);
			if(Constants.plasma()){
				this.ratio = mdf.add("ratio", Constants.initialCellularity[1], 0.01, 1.0, 1, true);
			}else{
			this.ratio = mdf.add("ratio", Constants.initialCellularity[1], 0.5, 3.0, numsamples, trainCellularity>=2 && Constants.useAvgDepth());
			
			}
			
					  
			
			this.refCount = new Double[ref.length()];
			
			this.maxCN = maxCN;
			this.ratiosR = new double[maxCN1-Constants.maxPloidy1];
			this.ratiosLe = new double[Constants.maxPloidy1];
			this.ratiosLe[0] = 0.001;
			for(int i=0; i<ratiosLe.length; i++){
				if(Constants.plasma) ratiosLe[i] = (double) i/ (double) Constants.maxPloidy1();
				else ratiosLe[i]=(1-mix)*((double)i/(double)Constants.maxPloidy1())+(mix) ;
			}
			for(int i=Constants.maxPloidy1+1; i<=maxCN1; i++){
				if(Constants.plasma) ratiosR[i-(Constants.maxPloidy1+1)] = (double) i/ (double) Constants.maxPloidy1();
				else ratiosR[i-(Constants.maxPloidy1+1)]=(1-mix)*((double)i/(double)Constants.maxPloidy1())+(mix) ;
			}
			if(trainPool  && ! Constants.plasma()){
				mdf.add("ratiosL", ratiosLe, ratiosLe.clone(), ratiosLe.clone());
				mdf.add("ratiosR", ratiosR, ratiosR.clone(), ratiosR.clone());
			}
			this.pool = new HaplotypeEmissionState("pool", ref.length(), Emiss.getSpaceForNoCopies(Constants.maxPloidy1()),(short)ref.dataIndex());
			pool.setNoCop(Constants.maxPloidy1());
			this.backCNall =new Double[ref.length()];
			if(ratiosL!=null){
				throw new RuntimeException("!!");
				/*for(int k=0; k<ratiosL.length; k++){
				backCNall[k] = find(ratios, ratiosL[k]);
			}
		
			{
				List<Double>[] l = new List[ratios.length];
				for(int k=0; k<l.length; k++){
					l[k] = new ArrayList<Double>();
				}
				for(int k=0; k<backCNall.length; k++){
					l[backCNall[k]].add(ratiosL[k]);
				}
				for(int k=0; k<l.length; k++){
					ratios[k] = avg(l[k],ratios[k],1);
					if(k==Constants.maxPloidy1()) ratios[k] = 1;
					else if(k<Constants.maxPloidy1()) ratios[k] = Math.min(ratios[k], 0.99);
					else if(k>Constants.maxPloidy1()) ratios[k] = Math.max(ratios[k], 1.01);
				}
			}*/
			}
			else{
				Arrays.fill(this.backCNall, (double) Constants.maxPloidy1());
			}
			if(trainPool)mdf.updateBounds();
			for(int i =0; i<refCount.length;i++){
				refCount[i] = ((IlluminaDistribution)ref.emissions[i]).b().doubleValue();
				pool.emissions[i] = new BackgroundDistribution(maxCN, refCount[i],backCNall[i], pool.getEmissionStateSpace(), indiv.size());
			}
			this.noprobes = noprobes;
			this.totnormal =((IlluminaRDistribution)ref.emissions[0]).r().doubleValue();
			mdf.print();
			
		}
		public Double b(Double b, Number r, int i, String name) {
		  double p = b/r.doubleValue();
		  double p1 =refCount(i)/totnormal();
		  double mult = p/p1;
		  
		 return inverseMult(mult, i,name);//(mult - (1-cellularity))* (ratio/cellularity);
		}
			
		
	
		private double avg(List<Double> list, double d1, double pseudo) {
		double d = d1*pseudo;
		for(int k=0; k<list.size(); k++){
			d+=list.get(k);
		}
		return d/(pseudo+(double)list.size());
	}
		private int find(double[] ratios2, double ratiosL) {
		int min = -1;
		double minv = Double.POSITIVE_INFINITY;
		int start =0;
		int end = ratios2.length;
		if(ratiosL<1) end = Constants.maxPloidy1();
		else if(ratiosL>1) start = Constants.maxPloidy1()+1;
		    for(int k=start; k<end; k++){
		    	double diff = Math.abs(ratios2[k] - ratiosL);
		    	if(diff<minv){
		    		minv = diff; min = k;
		    	}
		    }
		    return min;
	     }
		@Override
		public ProbabilityDistribution getDistribution(short data_index,
				int cn_bottom, int cn_top, int pos) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public double scoreB(short data_index, int j, double b, int i) {
			// TODO Auto-generated method stub
			return 0;
		}
	
		@Override
		public double scoreR(short data_index, int backgroundCount, int no_cop,
				double r, int i) {
			// TODO Auto-generated method stub
			return 0;
		}
	
		@Override
		public void addBCount(short data_index, int j, double weight, double val,
				int i) {
			// TODO Auto-generated method stub
	
		}
	
		@Override
		public void addRCount(short data_index, int backgroundCount, int no_cop,
				double weight, double val, int i) {
			// TODO Auto-generated method stub
	
		}
	
		
		
		@Override
		public ProbabilityDistribution getDistribution(short data_index, int j,
				int i) {
			return null;
		}
	int max_index=0;
		@Override
		public void maximisationStep(int i, double[] pseudo, double[] pseudoGlobal,
				List tasks) {
		{
			//for(int k=0; k<pool.emissions.length; k++){
			//	System.err.println(((BackgroundDistribution)this.pool.emissions[k]).ratio1_count/((BackgroundDistribution)this.pool.emissions[k]).cnt);
		//	}
				BackgroundDistribution bd = ((BackgroundDistribution)this.pool.emissions[0]);
				if(bd.cnt>0 && Constants.ratioAsLevels() && Constants.makeNewRef()){
					Logger.global.info("updating backgroun");
					for(int k=0; k<this.backCNall.length; k++){
						 bd = ((BackgroundDistribution)this.pool.emissions[k]);
						 bd.transfer(0);
						 this.backCNall[k] = bd.ratio1_count/bd.cnt;
						System.err.print(String.format("%3.2g,", backCNall[k]));
					}
					//System.err.println("\n");
				}
				if(max_index>0){
			   
				if(trainPool){
					mdf.updateBounds();
				}
				  
				   mdf.minimize();	
				  
				}
			
				
			}
		max_index++;
			mdf.print();
		}
	
		
		
		@Override
		public void print(File pw) {
			// TODO Auto-generated method stub
	
		}
	
		@Override
		public String getBName(int ij) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public Double sampleR(int data_index, int cn_bg, int cn_fg, int pos) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public Double sampleB(int data_index, int obj_index, int pos) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public ProbabilityDistribution getDistributionBfrac(short i2, int i, int j) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public ProbabilityDistribution getDistribution1(short data_index, int i,
				int j, int i2) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public void initialise() {
			this.average =0;
			this.average_count =0;
			this.pool.initialiseCounts();
			/*for(int k=0; k<this.pool.length(); k++){
				pool.emissions[k].initialise();
			}*/
		//	this.countProbesTumour = 0;
	
		}
	
		@Override
		public void addRBCount(short data_index, int j, double weight, double valR,
				double valB, int i) {
			
			System.err.println("h");
			// TODO Auto-generated method stub
	
		}
		
		public double inverseMult(double mult, int i, String name){
			//if(true) return mult;
			int sample_index2 = indiv.indexOf(name);
			double cellularity = this.cellularity[sample_index2]; double ratio = Constants.plasma() ? this.ratio[0] : this.ratio[sample_index2]; 
			double backCN = this.backCNall[i];

			 if(Constants.plasma()){
					double rcn1 = ratios(backCN);
				 return (mult*(rcn1*ratio + (1-ratio)) - (1-cellularity))/cellularity;
			 }else{
			if(ratioAsLevels){
				ratio = ratios(backCN)*ratio;
			}
			else if(cellAsLevels){
					cellularity =ratios(backCN)*cellularity;
			}
			return (ratio/cellularity) * ( mult - (1-cellularity) );
			 }
		}
		//static double mod = 2;
		public double getMult(double cellularity, double ratio, double rcn, double backCN){
			  if(Constants.plasma()){
				//	ratio = vals[1];
					double rcn1 = ratios(backCN);
			    	return ((rcn*cellularity) + (1-cellularity))/(rcn1*ratio + (1-ratio));
				}else{
					if(ratioAsLevels){
						
						ratio = ratios(backCN)*ratio;
					}
					else if(cellAsLevels){
							cellularity = Math.min(1.0, ratios(backCN)*cellularity);
							//if(cellularity>1) cellularity = 1.0/cellularity;
					}
					 return 	((rcn/ratio)*cellularity) + (1-cellularity);//
				}
		}
		
		
		class BinomialDistr1 implements ProbabilityDistribution2 {
		    Beta b = new Beta(0,0, new DRand());
		   
		 //   double[]logprobs = new double[numIt];
			final double rcn;
			final int cn;
			BetaBinom bb = new BetaBinom();
			List<Double> countRk = new ArrayList<Double>();
			List<Double> countRki = new ArrayList<Double>();
			//List<Double> val = new ArrayList<Double>();  //maintains the probabilities
	
			
			//List<Double> countref = new ArrayList<Double>();
			//List<Double> countratio = new ArrayList<Double>();
			
			private final double countnormal;
		  private double backCN= Constants.maxPloidy1();
			
			BinomialDistr1(int cn, double refCount, double backCN){
				this.cn = cn;
				this.rcn = (double)cn/Constants.backgroundCount1(0);
				this.countnormal =refCount;
				this.backCN = backCN;
			}
			public void setRefCount(double backCN) {
			//	this.countnormal =d;
				this.backCN = backCN;
				
			}
			
			public double eval(){
				double lp=0;
				for(int k =0; k<countRk.size(); k++){
					//this.setRefCount(countref.get(k), this.countratio.get(k));
					lp+=this.probability(countRk.get(k), countRki.get(k));//*val.get(k);
				}
				return lp;
			}
			
			@Override
			public int compareTo(Object arg0) {
				// TODO Auto-generated method stub
				return 0;
			}
			public ProbabilityDistribution2 clone() {
				return null;//new BinomialDistr(minx, maxx, miny, maxy);
			}
	//b/(a+b) = baf1  lrr = a+b
			//lrr == totTumour
			//baf = countTumour
			// ratio can be interpreted as an effective copynumber in the population
			
			
			
			
			public double probability(double lrr, double baf) {
				//double cellularity = 
				//double ratio = vals[1];
				double mult = getMult(cellularity[sample_index1], Constants.plasma() ? ratio[0] : ratio[sample_index1], rcn, backCN);
			 // double 	mult = ((rcn/ratio)*cellularity) + (1-cellularity);// this wa sold one
				
	//			double baf;
	/*			if(MatchedSequenceDataCollection.rescale){
					double modf = (totnormal+tottumour)/(2*tottumour);
					baf = baf1/modf;
				}else{
					baf = baf1;
				}*/
				//double countnormal = (1-baf)*lrr;
			   
				//double counttumour = baf;
				
				double k = baf;//counttumour;
				double n = lrr;// tottumour;
				
			    if(numIt==1){
			    	  double p =  (countnormal/totnormal)* mult;
			    	  double v = k *Math.log(p)+(n-k)*Math.log(1-p);
			    	  return  v; 
			    	
			    }else{
			    	//if(true){
			    	bb.set((mult*countnormal+1)/betaDownWeight, ((totnormal-mult*countnormal)+1)/betaDownWeight,n);
	//		    	bb.setBeta(
	//		    	bb.setTrials((int) Math.round(n));
			    	double v =  bb.getDensity(k);
			    //	System.err.println(rcn+" "+v);
			    	return v;
			    	//}
			    /*	 double maxv = Double.NEGATIVE_INFINITY;
			    	 b.setState((countnormal+1)/betaDownWeight, ((totnormal-countnormal)+1)/betaDownWeight);
				    for(int j = 0; j<numIt; j++){
						  double p =  b.nextDouble()* mult;
						  double v = k *Math.log(p)+(n-k)*Math.log(1-p);
						
						  if(v>maxv) maxv = v;
						  logprobs[j] = v;
				    }
					double sum=0;
					for(int j =0; j<numIt; j++){
						sum+=Math.exp(logprobs[j]-maxv);
					}
					return (Math.log(sum/numIt)+maxv);*/
			    }
			    }
			    
	
			
			public void addCount(double lrr, double baf, double value) {
				average += value * this.cn;
				average_count +=value;
				this.countRk.add(lrr); this.countRki.add(baf);
				//this.val.add(value);this.countref.add(countnormal);this.countratio.add(this.ratio1);
			}
	
			
			
			@Override
			public void initialise() {
				
				this.countRk.clear();
				this.countRki.clear();
				//this.countref.clear();
				//this.val.clear();
				
				
			}
	
			@Override
			public void transfer(double pseudoC) {
				this.initialise();
				
			}
	
			@Override
			public String id() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public int getParamIndex() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public String name() {
				// TODO Auto-generated method stub
				return null;
			}
	
			
			public void maximise(double d, double e, double f, double d1,
					double e1, double f1, double g) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public void updateParamIndex() {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public void recalcName() {
				// TODO Auto-generated method stub
				
			}
	
			
			public void getInterval(double[] input, DoubleMatrix2D res,
					double[] mean) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public int numObs() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			
			public int fill(DoubleMatrix2D x, DoubleMatrix2D y, DoubleMatrix2D yB,
					int numObs, double[] noCop, double pseudo) {
				// TODO Auto-generated method stub
				return 0;
			}
	
		
			public void variance(int type, double[] sum) {
				// TODO Auto-generated method stub
				
			}
	
			
			public void setParam(int type, int i, double d) {
				// TODO Auto-generated method stub
				
			}
	
			
			public int fillVariance(DoubleMatrix2D y, DoubleMatrix2D yb,
					DoubleMatrix2D covar, int numObs, double pseudo) {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public void print(PrintWriter pw) {
				// TODO Auto-generated method stub
				
			}
	
			
			public void setPriors(
					ProbabilityDistribution2 probabilityDistribution2, int type,
					boolean x) {
				// TODO Auto-generated method stub
				
			}
	
			
			
			public void setMinMax(double minR, double maxR, double min, double max) {
				// TODO Auto-generated method stub
				
			}
	
		
			public void setToExclude() {
				// TODO Auto-generated method stub
				
			}
	
			
			public double probability(double r, double b, int mixComponent) {
				// TODO Auto-generated method stub
				return 0;
			}
	
			
			public void addCount(Double r2, Double b, double val,
					SimpleExtendedDistribution1 mixe1,
					ProbabilityDistribution2 probDistG) {
				// TODO Auto-generated method stub
				
			}
	
			
	
			@Override
			public ProbabilityDistribution2 clone(double u) {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public ProbabilityDistribution2 clone(double u,
					SimpleExtendedDistribution1 simpleExtendedDistribution1) {
				// TODO Auto-generated method stub
				return null;
			}
	
			
			
		}
	
		
		
		//@Override
		public ProbabilityDistribution2 getDistributionRB(short data_index, int n,
				int noB, int i) {
			// TODO Auto-generated method stub
		      BinomialDistr1 dist = ((BackgroundDistribution)this.pool.emissions[i]).bin[sample_index1][n];//.get(n);
		   //   dist.setRefCount(refCount[i],ratio[i]);
		    //  dist.tumourCount=tumour[i];
		      return dist;
		}
	
		@Override
		public ProbabilityDistribution2 getDistributionRBGlob(short data_index,
				int n, int noB) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public double scoreRB(short data_index, int j, double r, double b, int i) {
			// TODO Auto-generated method stub
			return 0;
		}
	
		@Override
		public Double getFrac(int i, int ii, boolean r) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public Double getFracGlob(int ind, boolean po, boolean R) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public Double minQuality(int relative_position) {
			// TODO Auto-generated method stub
			return null;
		}
	
		@Override
		public void addRBCount(short data_index, int noCop, int noB, double val,
				Double r, Double b, int i) {
			//if(val>0.5){
			//	System.err.println("adding "+noCop+" "+r+" "+b+" "+val);\
			//this.countProbesTumour+=((double)noCop/2.0)*val; 
			 BackgroundDistribution dist = ((BackgroundDistribution)this.pool.emissions[i]);
		
		   //   dist.setRefCount(refCount[i], this.ratio[i]);
			 dist.addCount(r,b, val,noCop);
			//}
		}
	
		@Override
		public void print(int i) {
			// TODO Auto-generated method stub
	
		}
	
		@Override
		public String[] getFormGlob(int ind) {
			// TODO Auto-generated method stub
			return null;
		}
	
		
		public double evaluate() {
			//double lp=-Math.abs(arg0)*priorWeight;
			double lp =0;
			for(int k=0; k<this.pool.length(); k++){
				BackgroundDistribution dists = ((BackgroundDistribution)this.pool.emissions[k]);
				for(int j=0; j <numsamples; j++){//numsamples
					this.sample_index1 = j;
					for(int kk=0; kk<=maxCN; kk++){//noCN
						lp+=dists.bin[j][kk].eval();// dist.get(k).eval();
					}
				}
			}
			//System.err.println(arg0+" "+lp);
			return -lp;
		}
		
		public double evaluate(int j) {
			//double lp=-Math.abs(arg0)*priorWeight;
			double lp =0;
			for(int k=0; k<this.pool.length(); k++){
				BackgroundDistribution dists = ((BackgroundDistribution)this.pool.emissions[k]);
				//for(int j=0; j <numsamples; j++){//numsamples
					this.sample_index1 = j;
					for(int kk=0; kk<=maxCN; kk++){//noCN
						lp+=dists.bin[j][kk].eval();// dist.get(k).eval();
					}
				//}
			}
			//System.err.println(arg0+" "+lp);
			return -lp;
		}
	
	
	
		/*
		public double evaluate(double[] arg0) {
			// TODO Auto-generated method stub
			double lp =0;
			for(int k=0; k<arg0.length; k++){
				lp+=-Math.abs(arg0[k] - initial[k])/1e10;
			}
			vals[0] = arg0[0];
			vals[1] = arg0[1];
			if(trainPool)System.arraycopy(arg0, 2, this.ratios, 0, ratios.length);
			for(int k=0; k<this.pool.length(); k++){
				BackgroundDistribution dists = ((BackgroundDistribution)this.pool.emissions[k]);
				for(int j=0; j <numsamples; j++){
					for(int kk=0; kk<=maxCN; kk++){
					lp+=dists.bin[j][kk].eval();// dist.get(k).eval();
					}
				}
			}
		
		//	System.err.println(arg0[0]+" "+arg0[1]+" "+lp);
			return -lp;
		}*/
	
		
		
		public double totnormal() {
			return this.totnormal;
		}
	
		public class BackgroundDistribution extends PseudoDistribution{
			 double ratio1_count=0;
			 double cnt=0;
			
			
			BinomialDistr1[][] bin; //first is by sample, second is by CN.  All should have same backCN
			
			//following just for plotting
			private double refcount;
			double allr = 0;//new ArrayList<Double>(); //just for plotting
			double allb = 0;//new ArrayList<Double>(); //just for plotting
			
			public Number b1() {
				  MatchedDistributionCollection mdc = ((MatchedDistributionCollection)DataCollection.datC.dc);
				  double p = allb/allr;
				  double p1 = refcount/totnormal();
				  double mult = p/p1;
				  
				 return mult;//mdc.inverseMult(mult, i,name);
				
				//return ratios(bin[0][0].backCN);
			}
			public void addCount(Double r, Double b, double val,int noCop) {
			//	System.err.println(sample_index1);
				 bin[sample_index1][noCop].addCount(r, b, val); //this.pool[i].bin[noCop];
				 allr+=r;
				 allb+=b;
				    ratio1_count+=noCop*val;
				    cnt+=val;
				  //  System.err.println(cnt);
//				this.allr.add(r); this.allb.add(r);
			}
			/*public Number r1(int cn){
				bin[cn].
			}*/
			
			//EmissionStateSpace emStp;
			public BackgroundDistribution(int maxCN, double refCount, double backCN, EmissionStateSpace emstsp, int numsamp) {
				bin = new BinomialDistr1[numsamp][maxCN+1];
				this.refcount = refCount;
				//this.emstsp = emstsp;
				// TODO Auto-generated constructor stub
				for(int j=0; j<numsamp; j++){
				for(int i = 0; i<=maxCN; i++){
					bin[j][i] =new BinomialDistr1(i, refCount, backCN);
				}
				}
			}
			@Override
			public void addRBCount(EmissionStateSpace emStSp, int j, double val, int i) {
				if(true) throw new RuntimeException("!!");
				//double noCop = (double)emStSp.getCN(j)/(double) Constants.maxPloidy1();///Constants.backgroundCount(index);;//
				int noCop =emStSp.getCN(j);///(double) Constants.maxPloidy1();///Constants.backgroundCount(ind
			    ratio1_count+=noCop*val;
			    cnt+=val;
			  
			}
			
			@Override
			public double sample() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public void addCounts(ProbabilityDistribution probabilityDistribution) {
				// TODO Auto-generated method stub
				throw new RuntimeException("!!");
			}
	
			@Override
			public double[] probs() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public void setProb(double[] prob) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public void transfer(double pseudoC1) {
				double r1 = (Constants.maxPloidy()/Constants.backgroundCount(0));
				double ratio1 = (this.ratio1_count/this.cnt) * r1;
				//System.err.println("hh "+ratio1_count+" "+cnt);
				for(int ssi = 0; ssi < bin.length; ssi++){
			    for(int k=0; k<bin[ssi].length; k++){
			    	bin[ssi][k].setRefCount(ratio1);
			    }
				}
				
			}
	
			@Override
			public void addCount(int obj_index, double value) {
				System.err.println("h");
				
			}
	
			@Override
			public void initialise() {
				this.ratio1_count =0;
				this.cnt=0;
				allr =0;
				allb =0;
//					allb.clear();
				for(int k=0; k<this.bin.length; k++){
					for(int j=0; j<this.bin[k].length; j++){
				bin[k][j].initialise();
					}
				}
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public double probs(int obj_i) {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public double sum() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public int getMax() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public double[] counts() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public PseudoDistribution clone() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public void validate() {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public void setProbs(int to, double d) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public String getPrintString() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public PseudoDistribution clone(double swtch) {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public double logProb() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public void setCounts(int i1, double cnt) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public Integer fixedInteger() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public PseudoDistribution swtchAlleles() {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public double[] calcDistribution(double[] distribution,
					EmissionStateSpace emStSp, int pos) {
				// TODO Auto-generated method stub
				return null;
			}
	
			@Override
			public void setFixedIndex(int k) {
				// TODO Auto-generated method stub
				
			}
	
			@Override
			public double totalCount() {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public boolean isMeasured() {
				// TODO Auto-generated method stub
				return false;
			}
	
			@Override
			public double scoreB(int j, int i) {
				// TODO Auto-generated method stub
				return 0;
			}
	
			@Override
			public double scoreBR(EmissionStateSpace emStSp, int j, int i) {
				//double noCop = (double) emStSp.getCN(j)/(double) Constants.maxPloidy1();//Constants.backgroundCount(index);
				int noCop = emStSp.getCN(j);
			//	System.err.println(noCop);
				double lp =0;
				for(int k=0; k<numsamples; k++){
					sample_index1 =k ;
					for(int jj=0; jj<=maxCN; jj++){
					bin[k][jj].setRefCount(noCop);
					
					lp+=bin[k][jj].eval();
					}
				}
				return lp;
			}
			
		}



		
		
		
	}
