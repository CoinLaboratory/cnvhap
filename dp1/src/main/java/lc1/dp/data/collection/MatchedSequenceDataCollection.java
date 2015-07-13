package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;
import cern.jet.random.Binomial;
import cern.jet.random.Normal;

public class MatchedSequenceDataCollection extends IlluminaRDataCollection{
static double genomeLength = 1;//3.2e7;
	static String tum_id = "p0";
	static String ref_id = "ref";
	static String refall = Constants.reference();//"HG"REFALL";
	static int topoutlier =0;// 5;
	static int bottomoutlier = 0;
	//static double[] totals1 = new double[]{0,0};
	
	//final static double min_thresh = Constants.minNormalDepth();//5*Constants.cumulativeR(0);
	//public static boolean rescale = false; //rescaling for baf plots to be centred at 0.5
	
	public MatchedSequenceDataCollection(DataCollection dat) {
		super(dat);
		
		// TODO Auto-generated constructor stub
	}

	/*static String[] str_ = "2515532543:2483280919".split(":");
	public static double totnormal=  Double.parseDouble(str_[0]);
	public static	double tottumour =  Double.parseDouble(str_[1]);
	*/
	@Override    
public final void makeDistributions(int index) {
		double cellularity = Constants.initialCellularity()[0];
		int maxCN = Constants.maxCopies()*Constants.maxPloidy();
		int maxCN1 = Constants.maxCopies()*Constants.maxPloidy1;
		//double[]  geno = (double[])this.avgDepth.get(0);
		//int counta_ind = this.indexOf1(headsLowerCase,ref_id.split(":"));
		//int    countb_ind = this.indexOf1(headsLowerCase, tum_id.split(":"));
		
		
		//double a = (geno[counta_ind]); //normal
		//double b = (geno[countb_ind]); //tumour
		//totals[0] = b * 3.2e6;
		//totals[1] = a *3.2e6;
		HaplotypeEmissionState ref = ((HaplotypeEmissionState)this.dataL.get(refall));
	
		//double tottumour = totals[0];//
		/*double totnormal = totals[1];//
		 if(rescale){
			 double modf = (totnormal+tottumour)/(2*tottumour);
			IlluminaNoBg state = (IlluminaNoBg)this.dataL.get(this.indiv.get(0));
			
			 for(int j =0; j<state.length(); j++){
				 double baf = ((IlluminaDistribution)state.emissions[j]).b()*modf;
				  ((IlluminaNoBg)state).setB(j,baf);
					if(baf < Constants.minB(this.index)) Constants.minB[index] = baf;
					if(baf > Constants.maxB(this.index)) Constants.maxB[index] = baf;
				 
			 }
		 }*/
		
		//System.err.println(tottumour/(totnormal+tottumour));
		Double[] ratio = null;
		
		File f = new File("chrY.txt");
		if(f.exists()){
			ratio = new Double[this.loc.size()];
			double[] overlap = new double[this.loc.size()];
			Arrays.fill(ratio,0.0);
			Arrays.fill(overlap,0.0);
			try{
			BufferedReader br = new BufferedReader(new FileReader(f));
			String st = br.readLine();
			int kmin =0; int kmax = loc.size()-1;
			double base = 2.0;
			outer: while((st = br.readLine())!=null){
				if(st.startsWith("#")) continue outer;
				String[] str = st.split("\\s+");
			//	System.err.println(st);
				int start = Integer.parseInt(str[2]);
				int end = Integer.parseInt(str[3]);
				String[] gen = str[6].split(":");
				double cn=0;
				double indiv = 0;
				for(int k=0; k<gen.length; k++){
					double j =Double.parseDouble(gen[k]); 
					cn+=(k*j)/base; indiv +=j;
							
				}
				cn = cn/indiv;;
				inner: for(int k=kmin; k<kmax; k++){
					int bin_st = loc.get(k);
					int bin_end = loc.get(k+1);
					if(bin_st >= end) break inner;
					double overl = Math.min(bin_end, end) - Math.max(start, bin_st);
					if(overl>5){
						ratio[k] += cn*overl;
						overlap[k]+= overl;
					}
				}
				
				
			}
			br.close();
			for(int k=0; k<ratio.length; k++){
				if(overlap[k]>0){
					ratio[k] = ratio[k]/overlap[k];
				}else{
					ratio[k] = 1.0;
				}
			}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}

		this.dc = new MatchedDistributionCollection(index, cellularity,this.dir, maxCN, maxCN1,this.loc.size(), ref, ratio);
		
	/*	int len = this.indiv().size();
		Boolean[] probeOnly1 = new Boolean[len];
		Arrays.fill(probeOnly1,true);
	
	    	this.dc = new DistributionCollection(this.getBackgroundCN(),index,len,true, probeOnly1, this.dir);
	        for(int k=0; k<this.avgDepth.size(); k++){
	        	HaplotypeEmissionState hes = (HaplotypeEmissionState)this.dataL.get(this.indiv().get(k));
	        	EmissionStateSpace ems = Emiss.getSpaceForNoCopies(hes.noCop());
	        	for(int j=0; j<ems.cnLength(); j++){
	        		ProbabilityDistribution2 pds = dc.getDistributionRB(this.index, j, 0, k);
	        		System.err.println("h");
	        	}
	        }
	    }*/
	}
	protected void setMinMaxRValues() {
		// if(dc!=null) ((DistributionCollection)dc).setMinMax(Constants.minR(index), Constants.maxR(index),Constants.minB(index), Constants.maxB(index));
		
	}
	

	protected double defaultBAF() {
		// TODO Auto-generated method stub
		return 0;
	}
	
	static double resample(double tot, double n, double sub){
		if(Math.abs(sub-1)<0.001) {
			return n; 
		}
		double p = (n/tot)/sub;
		double d;
		//double p1 = Beta.staticNextDouble(n/sub+1,(tot-n/sub)+1);
		//double p2 = Normal.staticNextDouble(tot*p, tot*p*(1-p));
		if(tot>Integer.MAX_VALUE) {
			 d =  Normal.staticNextDouble(tot*p, Math.sqrt(tot*p*(1-p)));
			
		}else{
		  int n1 = (int) tot;
		  d =  Binomial.staticNextInt((int) tot, p);
		}
		return d;
		
	}
	
	public void mergeSamples( int len1) throws Exception{
		int len = Math.abs(len1);
        List<String> vals = new ArrayList<String>(this.indiv);
        List<String> todrop = new ArrayList<String>();
	  //  Map<String, EmissionState> states = new HashMap<String, EmissionState>();
		for(int k=0; k<vals.size(); k++){
		    String v = vals.get(k);
		    int index  =  Math.min(v.length(),len);
		    int ind1 = v.indexOf(".");
		    if(ind1>0){
		    	index = Math.min(ind1,index);
		    }
		    if(index<v.length() && ! v.equals(Constants.reference())){
			String v1 = v.substring(0,index);
			todrop.add(v);
			HaplotypeEmissionState hes = (HaplotypeEmissionState) dataL.get(v1);
			HaplotypeEmissionState hes0 = (HaplotypeEmissionState) dataL.get(v);
			if(hes==null){
				 hes = len1<0 ? (HaplotypeEmissionState) hes0.clone() : hes0;
				 hes.name = v1;
				 dataL.put(v1, hes);
				 this.indiv.add(v1);
			}else{
				PseudoDistribution[] dist = hes.emissions;
			    for(int k1=0; k1<dist.length; k1++){
			    	IlluminaDistribution dist1 = (IlluminaDistribution) hes0.emissions[k1];
			    	((IlluminaDistribution)dist[k1]).setBR(dist1.b(), dist1.r().doubleValue());
			    }
			}
		    }
			
		}
		if(len1>0) this.dropIndiv(todrop.toArray(new String[0]));
	}
	
	
	
	
protected Collection<? extends Integer> findLowDepth() {
	boolean mknewRef = Constants.makeNewRef();
	double[] tot = null;
	HaplotypeEmissionState href =  ((HaplotypeEmissionState)this.dataL.get(refall)) ;
	 PseudoDistribution[] emsref= href.emissions ;
	 if(mknewRef){
		// mknewRef = true;
		 tot = new double[emsref.length];
		 Arrays.fill(tot,0);
	 }else{
		 emsref = href.emissions;;
	 }
	for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
		HaplotypeEmissionState st = (HaplotypeEmissionState) it.next();
		String nme = st.name;
		double totD=0;
		PseudoDistribution[] ems = st.emissions;
        if(Constants.makeNewRef()){
        	if(!nme.equals(refall)){
				for(int k=0; k<ems.length; k++){
					double b = ((IlluminaDistribution) ems[k]).b();
			    	if(mknewRef){
			    		tot[k] +=b;
			    	}
				}
        	}
		}
        if(Constants.useAvgDepth()){
			totD = (Double) this.avgDepth.get(indiv.indexOf(nme));
			totD = totD *(genomeLength/Constants.downsample(1));
//		    System.err.println(totD+" "+nme);
		}else{
			for(int k=0; k<ems.length; k++){
				double b = ((IlluminaDistribution) ems[k]).b();
		    	totD+=b;
			}
		}
		for(int k=0; k<ems.length; k++){
	    	((IlluminaRDistribution) ems[k]).setR(totD);
		}
	}
	if(mknewRef){
		double totD = 0;
		if(Constants.useAvgDepth()){
			for(int k=0; k< this.avgDepth.size(); k++){
				if(!this.indiv.get(k).equals(refall)){
					totD +=(Double)avgDepth.get(k);
				}
			}
			totD = totD *(genomeLength/Constants.downsample(1));
		   // System.err.println(totD+" "+refall);
		}
		else{
		  for(int k=0; k<tot.length; k++){
			 totD+=tot[k];
		 }
		}
		for(int k=0; k<tot.length; k++){
			((IlluminaDistribution) emsref[k]).setBR1(tot[k], totD);
		}
	}
	
  
    List<Integer> res = 
			new ArrayList<Integer>();
    Comparator<Double>cd = new Comparator<Double>(){

		@Override
		public int compare(Double arg0, Double arg1) {
			return -1*arg0.compareTo(arg1);
		}
    	
    };
    double maxt=Double.POSITIVE_INFINITY;
    double maxn=Double.POSITIVE_INFINITY;
    double minn = 0;
    if(topoutlier>0 || bottomoutlier>0){
	//    List<Double>allt= new ArrayList<Double>();;
	    List<Double> alln = new ArrayList<Double>();
	  
	    for(int k=0; k<emsref.length; k++){
	    	//double r = ((IlluminaDistribution) emsref[k]).r().doubleValue();
	    	double b = ((IlluminaDistribution) emsref[k]).b();
	    	//double t = b*r;
	    	//double n = r - t;
	    	
	    		//allt.add(t);
	        	alln.add(b);
	    	
	    }
//	    Collections.sort(allt,cd);
	    Collections.sort(alln,cd);
	 //   maxt = allt.get(topoutlier);
	    maxn = alln.get(topoutlier);
	    minn = alln.get(alln.size()-bottomoutlier);
	    System.err.println(minn+" "+maxn);
    }
    for(int k=0; k<emsref.length; k++){
    	if(emsref[k]!=null){
    
    	double n = ((IlluminaDistribution) emsref[k]).b();
    //	double t = n;
    
    	if( (n <= minn && bottomoutlier >0) || (n>=maxn && topoutlier>0) || n <=Constants.minNormalDepth()){
    		res.add(k);
    	//	totals[0] -=t;
    //		totals[1] -=n;
    		System.err.println("dropped "+this.snpid.get(k)+" "+" "+n);
    	}
    	}
    }
  
		return res;
	}
	
	 @Override
	    public  Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy, double[] miss){
	   //    System.err.println(i);
	       HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);
	    	//int counta_ind = this.indexOf1(headsLowerCase,ref_id.split(":"));
			int    countb_ind;
			if(geno.length==1) countb_ind=0;
			else countb_ind = this.indexOf1(headsLowerCase, tum_id.split(":"));
			
		//	double a = Double.parseDouble(geno[counta_ind]); //normal
			double b = Double.parseDouble(geno[countb_ind]); //tumour
			
	
			double downsample = Constants.downsample(0); // tumour then ref
			//double[] tot = (double[])this.avgDepth.get(0); 
			/*double[] dilute = Constants.dilute();
			
				b = (b*dilute[0]+a*dilute[1])/(dilute[1]+dilute[0]); //dilution of tumour with normal;
			
			a = resample(tot[counta_ind], a, downsample[1]); //ref 
			b = resample(tot[countb_ind], b, downsample[0]); //tumour 
			*/
		/*	if(a+b<min_thresh){
				System.err.println("dropped site ");
				miss[0]+=1;
				return null;
			}*/
			if( Double.isNaN(b)){
				/*((IlluminaNoBg)state).emissions[i] = null;
				miss[0]++;
				miss[1]++;
				return null;*/
				//a=0;
				b=0;
			}
		//	else{
		//		System.err.println(indiv +" " +i);
		//	}
		//	totals[0]+=b; //tumour
		//	totals[1]+=a;  //normal
	
			//System.err.println(totals[0]+" "+totals[1]);
			// +totals1[0]+" "+totals1[1]);
			//double baf = b==0 ? 0 : (b/(a+b));
		//if(b> Constants.maxR(i)) Constants.maxR[index] = b;
		//if(!rescale){
		//	if(baf < Constants.minB(this.index)) Constants.minB[index] = baf;
		//	if(baf > Constants.maxB(this.index)) Constants.maxB[index] = baf;
		//}
			  ((IlluminaNoBg)state).setB(i, b/downsample);
	    	    return false;
	    }

	public MatchedSequenceDataCollection(File f, short index, int noCopies,
			int[][] mid, File bf,
			Collection<String> snpidrest) throws Exception {
		super(f, index, noCopies, mid, bf,snpidrest);
		this.dropIndiv(new String[] {refall});
		if(Constants.mergeSamplesLen()!=0){
			this.mergeSamples(Constants.mergeSamplesLen());
		}
		/*for(Iterator it = this.dataLvalues(); it.hasNext();){
			HaplotypeEmissionState hes = (HaplotypeEmissionState)it.next();
			for(int k=0; k<hes.emissions.length; k++){
				IlluminaDistribution ill = (IlluminaDistribution)hes.emissions[k];
				((MatchedDistributionCollection)this.dc).addRBCount(index, Constants.backgroundCount(index), 0, 1.0, ill.r().doubleValue(), ill.b().doubleValue(), k);
			}
		}*/
	//	this.dropIndiv(this.indiv.toArray(new String[0]));
		if(((MatchedDistributionCollection)this.dc).trainPool){
			this.dataL.put("pool",((MatchedDistributionCollection)this.dc).pool);
			this.indiv.add("pool");
		}
		// TODO Auto-generated constructor stub
	}
	
	
	
	


	}


	

//}
