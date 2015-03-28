	package lc1.dp.data.collection;
	
	import java.io.File;
	import java.io.PrintWriter;
	import java.util.ArrayList;
	import java.util.Arrays;
	import java.util.List;
	
	import lc1.dp.data.representation.Emiss;
	import lc1.dp.emissionspace.EmissionStateSpace;
	import lc1.dp.illumina.AbstractDistributionCollection;
	import lc1.dp.states.HaplotypeEmissionState;
	import lc1.stats.IlluminaDistribution;
	import lc1.stats.IlluminaRDistribution;
	import lc1.stats.ProbabilityDistribution;
	import lc1.stats.ProbabilityDistribution2;
	import lc1.stats.PseudoDistribution;
	import lc1.stats.SimpleExtendedDistribution1;
	import lc1.util.Constants;
	
	import org.apache.commons.math.special.Gamma;
	
	import pal.math.MultivariateFunction;
	import pal.math.MultivariateMinimum;
	import pal.math.OrthogonalHints;
	import pal.math.OrthogonalSearch;
	import pal.math.UnivariateFunction;
	import pal.math.UnivariateMinimum;
	import cern.colt.matrix.DoubleMatrix2D;
	import cern.jet.random.Beta;
	import cern.jet.random.engine.DRand;
	
	
	public class MatchedDistributionCollection extends
			AbstractDistributionCollection implements UnivariateFunction, MultivariateFunction{
	
		
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
		static final boolean trainPool = Constants.ratioAsLevels()!=null;
		static final int trainCellularity = Constants.trainCellularity();
		static final double mix = 0.5; //essentially controls the starting bands for the ratios.  This is effectively the minimunm value
		
		double[] vals; // cellularity, ratio
		
		final double[] ratios;
		  double[] initial =Constants.initialCellularity();
		 double[] upper = new double[] {initial[0],trainCellularity <2 ? initial[1] : 3.0};
		double[] lower = new double[] {trainCellularity <1 ?  initial[0]: 0.0, trainCellularity<2  ? initial[1] : 0.5};
	  
		
		private double totnormal;
		//private double prop = 1;
	//	private double tottumour;
		private double priorWeight = Constants.priorWeight();
		static int numIt = Constants.noSampleFromBeta();
		static double betaDownWeight = Constants.betaDownWeight();
		
		double noprobes;
		
		
		//double countProbesTumour=0;
		
		final double[] refCount;  //this lists the tumour count in each window
		final public HaplotypeEmissionState pool;
		
		final int[] backCNall;
		
		public double refCount(int i){
			return refCount[i];
		}
		public double ratio(int i){
			double ratio = this.vals[1];
			if(ratioAsLevels) ratio = ratios[backCNall[i]]*ratio;
			return ratio;
			//return ratios[backCNall[i]];
		}
		public double cell(int i){
			double cellularity = vals[0]; 
			
			if(cellAsLevels){
					cellularity = Math.min(1.0, ratios[backCNall[i]]*cellularity);
					//if(cellularity>1) cellularity = 1.0/cellularity;
			}
		  
		
		  return cellularity;
		
		}
		final short index;
		
		void updateBounds(){
		   
			
			for(int k=0; k<ratios.length; k++){
				if(k==Constants.maxPloidy1()){
					lower[2+k]=1.0;
					upper[2+k]=1.0;
				}else{
					  lower[2+k] = k==0 ? 0.001: ratios[k-1]+0.0001;
					  upper[2+k] = k<ratios.length-1 ? ratios[k+1]-0.0001 : ratios[k]+0.5;
					  if(k<Constants.maxPloidy1()){
						  lower[2+k] = Math.min(lower[2+k], 0.99999999);
						  upper[2+k] = Math.min(upper[2+k], 0.99999999);
					  }
					  else if(k>Constants.maxPloidy1()){
						  lower[2+k] = Math.max(lower[2+k], 1.0000001);
						  upper[2+k] = Math.max(upper[2+k], 1.0000001);
					  }
				  if(Math.signum(lower[2+k]-1) != Math.signum(upper[2+k]-1)){
					  throw new RuntimeException("!!");
				  }
				}
			}
			
		}
	  //  List<BinomialDistr1> dist = new ArrayList<BinomialDistr1>();
		public MatchedDistributionCollection(int index, 
			 double cellularity, File dir, int maxCN,int maxCN1, int noprobes,HaplotypeEmissionState ref, Double[] ratiosL) {
			// TODO Auto-generated constructor stub
			this.index = (short)index;
			this.refCount = new double[ref.length()];
			this.ratios = new double[maxCN1+1];
			if(trainPool){
			//	lower[1] = 1.0;
			//	upper[1] = 1.0;
				double[] lower1 = new double[lower.length+ratios.length];
				double[] upper1 = new double[lower.length+ratios.length];
				double[] initial1 =new double[lower.length+ratios.length];
				System.arraycopy(lower, 0, lower1, 0, lower.length);
				System.arraycopy(upper, 0, upper1, 0, upper.length);
				System.arraycopy(initial, 0, initial1, 0, initial.length);
				lower = lower1;
				upper = upper1;
				initial = initial1;
			}
			ratios[0] = 0.001;
			for(int i=0; i<ratios.length; i++){
				ratios[i]=(1-mix)*((double)i/(double)Constants.maxPloidy1())+(mix) ;
			}
			if(trainPool) System.arraycopy(ratios,0,initial,2,ratios.length);
			this.pool = new HaplotypeEmissionState("pool", ref.length(), Emiss.getSpaceForNoCopies(Constants.maxPloidy1()),(short)ref.dataIndex());
			pool.setNoCop(Constants.maxPloidy1());
			this.backCNall =new int[ref.length()];
			if(ratiosL!=null){
			for(int k=0; k<ratiosL.length; k++){
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
			}
			}
			else{
				Arrays.fill(this.backCNall, Constants.maxPloidy1());
			}
			if(trainPool)updateBounds();
			for(int i =0; i<refCount.length;i++){
				refCount[i] = ((IlluminaDistribution)ref.emissions[i]).b().doubleValue();
				pool.emissions[i] = new BackgroundDistribution(maxCN, refCount[i],backCNall[i], pool.getEmissionStateSpace());
			}
			this.noprobes = noprobes;
			this.vals = new double[] {cellularity,initial[1]};
			//this.vals = new double[][] {this.
			//this.ratio =1.0;
			//this.tottumour = tottumour;
			this.totnormal =((IlluminaRDistribution)ref.emissions[0]).r().doubleValue();
			//double pr = tottumour/(totnormal+tottumour);
			//this.cellularity = cellularity;
			print();
			
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
	
		@Override
		public void maximisationStep(int i, double[] pseudo, double[] pseudoGlobal,
				List tasks) {
			if(false){
				UnivariateMinimum uvm = new UnivariateMinimum();
				double cell = Constants.updateCellularity()? uvm.findMinimum(vals[0], this, 3) : vals[0];
				//vals[1] = Math.min( Math.max(noprobes/countProbesTumour,0.5),2.0);
				
				vals[0]= cell;
			}else{
				
				
				
				BackgroundDistribution bd = ((BackgroundDistribution)this.pool.emissions[0]);
				if(bd.cnt>0 && trainPool){
				
					for(int k=0; k<this.backCNall.length; k++){
						 bd = ((BackgroundDistribution)this.pool.emissions[k]);
						 bd.transfer(0);
						this.backCNall[k] = bd.ratio1_count/bd.cnt;
					}
				}
				{
				   MultivariateMinimum mvm =new OrthogonalSearch(); // new ConjugateGradientSearch();//n
				double[] vals1 = new double[this.getNumArguments()];
				System.arraycopy(vals,0,vals1,0,vals.length);
				if(trainPool){
					System.arraycopy(ratios,0,vals1,2,ratios.length);
					 this.updateBounds();
				}else{
				if(!Constants.useAvgDepth()){
					double r = (this.average/this.average_count)/Constants.backgroundCount((int)this.index);
			  this.lower[1] = r;
				  this.upper[1] = r;
				  vals1[1] = r;
				}
				}
				   mvm.findMinimum(this, vals1, 1, 3);
				   System.arraycopy(vals1,0,vals,0,vals.length);
					if(trainPool)System.arraycopy(vals1,2,ratios,0,ratios.length);
			//	this.pool.initialiseCounts();
				
				}
				//this.cellularity = vals[0];
				//this.ratio = vals[1];
				
			}
			//System.err.println(Arrays.asList(ratio));
		
			//System.err.println("avg r " +r);
			print();
		}
	
		
		public void print(){
			System.err.println("new val "+vals[0]+" "+vals[1]);
			StringBuffer  sb = new StringBuffer("initial:  ");
			for(int k=0; k<this.initial.length; k++) sb.append(" "+initial[k]);
			System.err.println(sb.toString());
			sb = new StringBuffer("lower: ");
			for(int k=0; k<this.lower.length; k++) sb.append(" "+lower[k]);
			System.err.println(sb);
			sb = new StringBuffer("upper: ");
			for(int k=0; k<this.upper.length; k++) sb.append(" "+upper[k]);
			System.err.println(sb);
			sb = new StringBuffer("vals:  "+vals[0]+" "+vals[1]);
			for(int k=0; k<this.ratios.length; k++) sb.append(" "+ratios[k]);
			System.err.println(sb);
			
			
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
		
		
		
		class BinomialDistr1 implements ProbabilityDistribution2 {
		    Beta b = new Beta(0,0, new DRand());
		   
		    double[]logprobs = new double[numIt];
			final double rcn;
			final int cn;
			 BetaBinom bb = new BetaBinom();
			List<Double> countRk = new ArrayList<Double>();
			List<Double> countRki = new ArrayList<Double>();
			//List<Double> val = new ArrayList<Double>();  //maintains the probabilities
	
			
			//List<Double> countref = new ArrayList<Double>();
			//List<Double> countratio = new ArrayList<Double>();
			
			private final double countnormal;
			private int backCN= Constants.maxPloidy1();
			
			BinomialDistr1(int cn, double refCount, int backCN){
				this.cn = cn;
				this.rcn = (double)cn/Constants.backgroundCount(0);
				this.countnormal =refCount;
				this.backCN = backCN;
			}
			public void setRefCount(int backCN) {
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
				double cellularity = vals[0]; 
				double ratio = vals[1];
				if(ratioAsLevels) ratio = ratios[backCN]*ratio;
				else if(cellAsLevels){
						cellularity = Math.min(1.0, ratios[backCN]*cellularity);
						//if(cellularity>1) cellularity = 1.0/cellularity;
				}
			    double 	mult = ((rcn/ratio)*cellularity) + (1-cellularity);
			
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
		      BinomialDistr1 dist = ((BackgroundDistribution)this.pool.emissions[i]).bin[n];//.get(n);
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
			
			 BinomialDistr1 dist = ((BackgroundDistribution)this.pool.emissions[i]).bin[noCop]; //this.pool[i].bin[noCop];
		   //   dist.setRefCount(refCount[i], this.ratio[i]);
			 dist.addCount(r,b, val);
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
	
		@Override
		public double evaluate(double arg0) {
			//double lp=-Math.abs(arg0)*priorWeight;
			double lp =0;
			vals[0] = arg0;
			for(int k=0; k<this.pool.length(); k++){
				BackgroundDistribution dists = ((BackgroundDistribution)this.pool.emissions[k]);
				for(int j=0; j <dists.bin.length; j++){
					lp+=dists.bin[j].eval();// dist.get(k).eval();
				}
			}
			//System.err.println(arg0+" "+lp);
			return -lp;
		}
	
		@Override
		public double getLowerBound() {
			// TODO Auto-generated method stub
			return 0;
		}
	
		@Override
		public double getUpperBound() {
			// TODO Auto-generated method stub
			return 1;
		}
	
		@Override
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
				for(int j=0; j <dists.bin.length; j++){
					lp+=dists.bin[j].eval();// dist.get(k).eval();
				}
			}
		
		//	System.err.println(arg0[0]+" "+arg0[1]+" "+lp);
			return -lp;
		}
	
		
		@Override
		public double getLowerBound(int arg0) {
			// TODO Auto-generated method stub
		 return lower[arg0];
			
		}
	
		@Override
		public int getNumArguments() {
			// TODO Auto-generated method stub
			return this.vals.length+(trainPool ? ratios.length :0);
		}
	
		@Override
		public OrthogonalHints getOrthogonalHints() {
			// TODO Auto-generated method stub
			return null;
		}
		
		@Override
		public double getUpperBound(int arg0) {
			// TODO Auto-generated method stub
			return upper[arg0];
		}
	
		public double totnormal() {
			return this.totnormal;
		}
	
		public class BackgroundDistribution extends PseudoDistribution{
			private int ratio1_count=0;
			private int cnt=0;
			
			
			BinomialDistr1[] bin;
			
			
			
			public Number b1() {
				return ratios[bin[0].backCN];
			}
			/*public Number r1(int cn){
				bin[cn].
			}*/
			
			//EmissionStateSpace emStp;
			public BackgroundDistribution(int maxCN, double refCount, int backCN, EmissionStateSpace emstsp) {
				bin = new BinomialDistr1[maxCN+1];
				//this.emstsp = emstsp;
				// TODO Auto-generated constructor stub
				for(int i = 0; i<=maxCN; i++){
					bin[i] =new BinomialDistr1(i, refCount, backCN);
					
				}
			}
			@Override
			public void addRBCount(EmissionStateSpace emStSp, int j, double val, int i) {
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
				int ratio1 = this.ratio1_count/this.cnt;
				//System.err.println("hh "+ratio1_count+" "+cnt);
			    for(int k=0; k<bin.length; k++){
			    	bin[k].setRefCount(ratio1);
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
					
				for(int k=0; k<this.bin.length; k++){
				bin[k].initialise();
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
				for(int k=0; k<this.bin.length; k++){
					bin[k].setRefCount(noCop);
					
					lp+=bin[k].eval();
				}
				return lp;
			}
			
		}
		
		
	}
