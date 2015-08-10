 package lc1.dp.data.collection;
//TO DO - fix HWE calculation
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;

import lc1.CGH.Aberation;
import lc1.CGH.Location;
import lc1.CGH.Locreader;
import lc1.dp.data.representation.AbstractEmiss;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;
import lc1.dp.emissionspace.SimpleEmissionStateSpace;
import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.BackgroundEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.ensj.PedigreeDataCollection;
import lc1.possel.PlotLD;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.IntegerDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.stats.SkewNormal;
import lc1.stats.TrainableNormal;
import lc1.util.ApacheCompressor;
import lc1.util.CompressDir;
import lc1.util.Constants;
import lc1.util.PhenoGroup;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipFile;

import calc.DPrimeCalculator;
import calc.LDCalculator;
import calc.R2Calculator;


public  class DataCollection extends DataC implements DataCInterface{
	public static DataCollection datC=null;
	public  boolean strand_represents_ref =false;	
  public String chrom;
  public AbstractDistributionCollection dc=null;
  public String chrom(){
	  return chrom;
  }
//  final double[] totals = new double[]{0,0}; //for matched analysis
  public void subsample(String[] det){
	 
	  int[] freq = new int[det.length];
	  int[] dist = new int[det.length];
	  for(int k=0; k<freq.length; k++){
		  String[] spl = det[k].split("_");
		  freq[k] = Integer.parseInt(spl[0]);
		  dist[k] = spl.length>1 ? 
			  
			  Constants.convert(spl[1]) : Integer.MAX_VALUE;
	  }
	  if(freq.length==1 && freq[0]<=1) return;
	  	int[] mid = Constants.mid()[0];
	  	int[] start = new int[freq.length];
	  	int[] end = new int[freq.length];
	  	List<int[]> li = new ArrayList<int[]>();
	  //	boolean[] include = new boolean[freq.length];
	  	for( int k=0; k<dist.length; k++){
	  		
	  		start[k] =dist[k] ==Integer.MAX_VALUE ? 0 : Math.max(0, lastLessThan(loc,mid[0] - dist[k],0));
	  		end[k] = dist[k] ==Integer.MAX_VALUE ? loc.size() : firstGreaterThan(loc, mid[1]+dist[k]);
	  		if(k==0){
	  			li.add(new int[] {start[k], end[k],freq[k]});
	  		}
	  		else{
	  			if(start[k] < start[k-1]){
	  				li.add(new int[] {start[k], start[k-1],freq[k]});
	  			}
	  			if(end[k]>end[k-1]){
	  				li.add(new int[] {end[k-1], end[k],freq[k]});
	  			}
	  		}
	  		//if(start[k] )
	  	//	if(k>0 && start[k])
	  	}
	  	Collections.sort(li, new Comparator<int[]>(){

			

			@Override
			public int compare(int[] arg0, int[] arg1) {
				if(arg0[0] !=arg1[0]) return arg0[0] < arg1[0]  ? -1 :1;
				return 0;
			}
	  		
	  	});
	  	 start = new int[li.size()];
	  	 end = new int[li.size()];
	  	 freq = new int[li.size()];
	  	for(int k=0; k<start.length; k++){
	  		int[] d = li.get(k);
	  		start[k] = d[0]; end[k] = d[1]; freq[k] = d[2];
	  	}
		for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
			((HaplotypeEmissionState)it.next()).adjustDepth();
		}
	//	int extra = Constants.extra;
		for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
			((HaplotypeEmissionState)it.next()).subSample(freq, start,end);
		}
		snpid = subsample(this.snpid,freq,start,end);
		strand = subsample(this.strand,freq,start,end);
		loc = subsample(this.loc,freq,start,end);
		alleleA = subsample(this.alleleA,freq,start,end);
		alleleB = subsample(this.alleleB,freq,start,end);
		Object[] obj = subsample(Arrays.asList(probeOnly), freq,start,end).toArray(new Object[0]);
		this.probeOnly = new Boolean[obj.length];
		System.arraycopy(obj, 0, probeOnly, 0, obj.length);
		this.makeDistributions(this.index);
		this.length = loc.size();
	//	this.
	}
  private int lastLessThan(List<Integer> loc, int i, int start) {
	for(int j=start; j<loc.size(); j++){
		if(i<loc.get(j)) 
			return j-1;
	}
	return loc.size();
}
public static List subsample(List l, int freq) {
	if(freq==1) return l;
      int newlen = (int)Math.ceil((double) (l.size())/ (double) freq);
	
	List newl = new ArrayList();
	for(int k=0; k<newlen; k++){
		newl.add(l.get((int) k*freq));
	}
	return newl;
	
}
  public static List subsample(List l, int[] freq, int[] start, int[] end){
	  if(l.size()==0) return l;
	  List l1 = new ArrayList();
	  for(int k=0; k<freq.length; k++){
		  l1.addAll(subsample(l.subList(start[k],end[k]),freq[k]));
	  }
	  return l1;
  }
  
  public Map<String, String> rename; //used if overlaps with other datasets, maps from original name to new name
// public  BackgroundEmissionState bg;
// public PIGData bg1;
  
 // double[]cn_count;
  protected BackgroundEmissionState bg;
  protected PhasedDataState bg1;
  
public PseudoDistribution getBGDist(int data_index, int pos){
	return bg.emissions[pos];
}
public Integer getBGCount(int data_index, int pos){
	/*if(Constants.modelbg())
	return bg.getFixedInteger(pos);
	else */
		return Constants.backgroundCount(data_index);
	/*Integer val = bg.getFixedInteger(pos);
	if(val!=null) return bg.getEmissionStateSpace().getCN(val);
	else return null;
	*/
}


  public List<Integer> loc() {
  	return loc;
  }

  public List<Character> alleleA() {
  	return alleleA;
  }

  public List<Character> alleleB() {
  	return alleleB;
  }

  public String name() {
  	return name;
  }

  public List<String> snpid() {
  	return snpid;
  }
  Comparator<String> comp1 = new Comparator<String>(){

		@Override
		public int compare(String arg0, String arg1) {
			return arg0.compareTo(arg1)*-1;
		}
  };
    public  Map<String, PIGData> data = new TreeMap<String, PIGData>(comp1);
    public  Map<String, EmissionState> dataL = new TreeMap<String, EmissionState>(comp1);
    //}
    public  Map<String, EmissionState> dataLSt = new TreeMap<String, EmissionState>(comp1);
  //  public Map<String, Boolean> cc = new HashMap<String, Boolean>();
  // int no_copies =2;
  // EmissionStateSpace[] stSp; //index is no of copies -1 
   // EmissionStateSpace[] stSp1;
    //EmissionStateSpaceTranslation[] trans;
    
    public short index = -1;
    List<Double>[] hwe ;
    /*public void calcHWE(boolean state, boolean exclMissing){
       
       hwe =  new List[this.getNumberDataTypes()];//new ArrayList<Double>();
        for(int i=0; i<this.getNumberDataTypes(); i++){
            hwe[i] = new ArrayList<Double>();
        this.calcHWE(state, exclMissing, hwe[i], null, i);
       
        }
    }*/
    public String getCompressedString(String key, int i, boolean likelihood, boolean state){
    	HaplotypeEmissionState st =null;
    	if(state){
    		
    		if(likelihood) st =  ((HaplotypeEmissionState)this.viterbiL.get(key));
        	else{
        		st = ((HaplotypeEmissionState) this.viterbi.get(key));
        	}
    	}
    	else{
	    	if(likelihood)st =  ((HaplotypeEmissionState) dataL(key));
	    	else{
	    		st = ((HaplotypeEmissionState) data.get(key));
	    	}
    	}
    	if(st==null) return "";
    	else return st.getCompressedDataString(i);
    }
  /* public void dropFixed(){
	  
	   this.calculateMaf(true);
	   List<Integer> toD = new ArrayList<Integer>();
	   for(int i=0; i<loc.size();i++){
		   if(this.maf.getEmiss(i)[this.maf.getBestIndex(i)]>0.9999){
			   toD.add(i);
		   }
	   }
	   this.drop(toD, false);
   }*/
   
    public  void applyVarianceThreshold(ZipFile zf, 
    		boolean standardise,List<Integer> dToInc , int sampleid) throws Exception{
        List<String >indiv = ApacheCompressor.getIndiv(zf, "Samples", sampleid);
        int index1 = header_sample.indexOf("variance");
        if(index1 <0){
            Logger.global.warning("no variance column ");
            return ;
        }
        
        List<String >mc = ApacheCompressor.getIndiv(zf, "Samples", index1);
       double[] mcs = new double[mc.size()];
        List<String> toDel = new ArrayList<String>();
   	 for(int k_=0; k_<dToInc.size(); k_++){
		 int i = dToInc.get(k_);
       
        	mcs[i] = Double.parseDouble(mc.get(i));//Math.sqrt();
            if(mcs[i] > Constants.var_thresh(this.index) || (this.dp!=null && dp.todrop.contains(indiv.get(i)))){
                System.err.println("removing "+indiv.get(i)+" "+mc.get(i));
                String indiv_i = indiv.get(i);//.split("#")[0];
               toDel.add(indiv_i);
            }
           
        }
       
     //   int len = mcs.length-toDel.size();
        if(standardise){
        	System.err.println("STANDARDISING VARIANCE");
        double mid = Constants.median(mcs, mcs.length-toDel.size());
        for(int i=0; i<indiv.size(); i++){
        	HaplotypeEmissionState ems = (HaplotypeEmissionState) this.dataL.get(indiv.get(i));
        	//double stderr =  mcs[i]);
        	double mult = Math.sqrt(mid/mcs[i]);
        	ems.standardiseVariance(mult);
        }
        }
        this.dropIndiv(toDel.toArray(new String[0]));
        this.indiv =indiv();
       
        //return mid;
    }
    
   
   
    
    public void randomizePhenotypes(){
        this.currentPhenScIndex = -1;
        this.currentPosScIndex = -1;
       this.currentType=-1;
       
        List<Double[]> phen = new ArrayList<Double[]>(dataL.size());
        for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
          phen.add( ( (HaplotypeEmissionState)it.next()).phenValue());
        }
        Collections.shuffle(phen);
        int i=0;
        for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();i++){
             ( (HaplotypeEmissionState)it.next()).setPhenotype(phen.get(i));
          }
    }
    public void applyMedianCorrection(ZipFile zf, List<Integer> toinc) throws Exception{
        List<String >indiv = ApacheCompressor.getIndiv(zf, "Samples", header_sample.indexOf("id"));
        int index1 = header_sample.indexOf("median_loess");
        int index2 = header_sample.indexOf("median");
        int index3 = header_sample.lastIndexOf("median");
        if(index3!=index2){
        	index1 = index2+1;
        }
        if(index1 <0){
            Logger.global.warning("no median correction column ");
            return;
        }
        List<String >median_loess = ApacheCompressor.getIndiv(zf, "Samples", index1);
    //    List<String >median = Compressor.getIndiv(zf, "Samples", index2);
        for(int i_=0; i_<toinc.size(); i_++){
        	int i = toinc.get(i_);
           String indiv_i = indiv.get(i);//.split("#")[0];
           HaplotypeEmissionState state =  (HaplotypeEmissionState)dataL.get(indiv_i);
           
           if(state==null){
               throw new RuntimeException("no state for "+indiv_i);
           }
           double medianL = Double.parseDouble(median_loess.get(i));
           /*if(medianL==0){
               medianL = Double.parseDouble(median.get(i));
           }*/
           state.applyCorrection(-1*medianL);
        }
    }
    
   
    
    
   
    
    public Set<String> applyLoess(ZipFile zf, Map<String, String>platem, boolean loess_, boolean gc_,String pref) throws Exception{
        BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(pref+"SNPS"))));
        String st = br.readLine();
        Set<String> todelete = new HashSet<String>();
        List<String> plates = new ArrayList<String>(new HashSet<String>(platem.values()));
        Double[][] loess = new Double[plates.size()][snpid.size()];
        Double[][] gc = new Double[plates.size()][snpid.size()];
        int[][] inds= new int[plates.size()][2];
        int offset=0;
        int len = st.split("\\s+").length;
        if(len<this.header_snp.size()){
     	   Logger.global.warning("SNPS is not same length as name - adjusting");
     	   offset  = header_snp.size()-len;
        }
        else if (len > header_snp.size()){
     	   throw new RuntimeException("!!");
        }
       
        for(int i=0; i<plates.size(); i++){
        	inds[i][0] = header_snp.indexOf(plates.get(i)+
        			(plates.get(i).length()==0 ? "" : "_") +
        					"loess") - offset;
        	inds[i][1] =  header_snp.indexOf(plates.get(i)+
        			(plates.get(i).length()==0 ? "" : "_") +
        			"gc")-offset;
        	 if(inds[i][0]<0 ){
        		 loess[i] = null;
        		 gc[i] = null;
              //  throw new RuntimeException("no loess column for "+plates.get(i)+" "+this.name+"\t This may be because there is no plate.txt file");
             //    return;
             }
        }
        int id_index = this.header_snp.indexOf("id");
     
       
        while((st)!=null){
            String[] str = st.split("\\s+");
            int ind = this.snpid.indexOf(str[id_index]);
            if(ind>=0){
            	for(int j=0; j<inds.length; j++){
            		int loess_index = inds[j][0];
            		if(loess_index>=0){
            		 loess[j][ind] = 
                         str[loess_index].equals("NA") ? Double.NaN : 
                         Double.parseDouble(str[loess_index]);
            		int gc_index = inds[j][1];
               
                   gc[j][ind] = gc_index<0 ? 0 : 
                         Double.parseDouble(str[gc_index]); 
            		}
            	}
            }
            st = br.readLine();
        }
        for(int j=0; j<plates.size(); j++){
        	if(loess[j]!=null){
	        for(int i=0; i<loess[j].length; i++){
	          // if(loess[j][i]==0) throw new RuntimeException("no loess for "+this.snpid.get(i));
	           if(Double.isNaN(loess[j][i])){
	                loess[j][i] = average(loess[j], i);
	            }
	        }
        	}
        }
        for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
        	HaplotypeEmissionState nxt = ((HaplotypeEmissionState)it.next());
        	String key = find(platem,nxt.getName());
        	int j = 
        		
        		plates.indexOf(key);
     //   	System.err.println(nxt.getName()+" "+j);
        	if(j<0){
        		todelete.add(nxt.getName());
        		System.err.println("DELETING SAMPLE "+nxt.getName()+" because we cannot find it in plate file");
        	}
        	else{
           if(loess_) nxt.applyLoess(loess[j]);
           if(gc_)  nxt.applyLoess(gc[j]);
        	}
        }
        return todelete;
    }
    
   private String find(Map<String, String> platem, String name) {
		String res = platem.get(name);//JHP656_1673047051_A	
		if(res==null){
			String key = null;
			for(Iterator<String> it =platem.keySet().iterator(); it.hasNext();){
				String k = it.next();
				if(k.startsWith(name)){
					if(key!=null && !key.equals(k) && !platem.get(key).equals(platem.get(k))) {	
						System.err.println("plate conflict "+k+" "+key+" "+name);
					}
					key = k;
				}
			}
			res = platem.get(key);
		}
		return res;
	}
private Double average(Double[] loess, int i) {
       int before_index = i-1;
       int after_index = i+1;
       
           while(before_index>=0 && Double.isNaN(loess[before_index])){
               before_index--;
           }
      while(after_index<loess.length && Double.isNaN(loess[after_index])){
          after_index++;
      }
      if(before_index<0){
    	  if(after_index>=loess.length) return 0.0;
    	  return loess[after_index];
      }
      if(after_index>=loess.length) return loess[before_index];
      else return (loess[after_index] + loess[before_index])/2.0;
    }

protected int getNumberDataTypes() {
       return 1;
    }

   /* public  void calcHWE(boolean state, boolean exclMissing, List<Double> hwe, List<Double> hwe_Exp, int data_index){
        hwe.clear();
       if(hwe_Exp!=null) hwe_Exp.clear();
        for(int i=0; i<this.length(); i++){
            hwe.add(this.getHWE(i, state, exclMissing, data_index));
            if(hwe_Exp!=null) calculatedExpectedHWE(hwe, hwe_Exp,i);
        }
       
    }*/
    public void calculatedExpectedHWE(List<Double> hwe, List<Double> hwe_Exp, int data_index){
        List<Double[]> l = new ArrayList<Double[]>();
        for(int i=0; i<hwe.size(); i++){
            l.add(new Double[] {hwe.get(i), 0.0, (double) i});
        }
        Collections.sort(l,new Comparator<Double[]>(){

            public int compare(Double[] arg0, Double[] arg1) {
               return arg0[0].compareTo(arg1[0]);
            }
            
        });
        for(int i=0; i<l.size(); i++){
            Double[] d = l.get(i);
            d[1] = (((double) i)+0.5)/(double) l.size();
        }
        Collections.sort(l,new Comparator<Double[]>(){

            public int compare(Double[] arg0, Double[] arg1) {
               return arg0[2].compareTo(arg1[2]);
            }
            
        });
      hwe_Exp.clear();
        for(int i=0; i<l.size(); i++){
            Double[] d = l.get(i);
           hwe_Exp.add(d[1]);
        }
    }
    public EmissionState getState(String key, EmissionStateSpace stSp) {
        if(stSp==null) return null;
        EmissionState st =  dataLSt.get(key);
        if(st!=null){
            if(stSp!=st.getEmissionStateSpace()) throw new RuntimeException("!!");
            return st;
        }
        else{
            dataLSt.put(key, st = new HaplotypeEmissionState(key, length, stSp.size(), stSp,  null, null));
        }
        return st;
    }
    /* first key is haplotype, second key is phenotype value */
    public SortedMap<String,Map<String,List<PIGData> >> getAllHaplotypes(int start, int end, String pheno2) {
        SortedMap<String,List<PIGData> >map =  getAllHaplotypes(start, end);
        SortedMap<String,Map<String,List<PIGData> >> res = new TreeMap<String, Map<String, List<PIGData>>>();
        int phenoIndex =-1;
        if(pheno!=null){
        for(int i=0; i<this.pheno.size(); i++){
            if(pheno.phen.get(i).startsWith(pheno2)){
                phenoIndex = i;
            }
        }
        }
        //if(phenoIndex <0) return null;
        Map<Double, String> dec = pheno==null ? null :  pheno.reverse(phenoIndex);
        for(Iterator<Map.Entry<String,List<PIGData> >> it = map.entrySet().iterator(); it.hasNext();){
            Map.Entry<String,List<PIGData>> nxt = it.next();
            Map<String, List<PIGData>> val = new TreeMap<String, List<PIGData>>();
            res.put(nxt.getKey(), val);
            for(int i=0; i<nxt.getValue().size(); i++){
                PhasedDataState hes = (PhasedDataState) nxt.getValue().get(i);
               
                String phenV = "null";
                    if(phenoIndex >=0){
                        Double[] phenVV = hes.phenValue();
                    phenV = dec==null ? phenVV[phenoIndex].toString() : dec.get(phenVV[phenoIndex]);
                    }
                    List<PIGData> val1 = val.get(phenV);
                    if(val1==null){
                        val.put(phenV, val1 = new ArrayList<PIGData>());
                    }
                val1.add(hes);
            }
        }
        return res;
    }
    
    public SortedMap<String,List<PIGData> >  getAllHaplotypes(int start, int end){
        
        return getHaplotypes(start, end, new HashSet<PIGData>(this.data.values()));
    }
    //first is below and hence DPrime, second is above and hence r2
   
    public SortedSet<Entry<String, List<PIGData>> >
        getHaplotypes(Map<String,List<PIGData> > left)
        {
        SortedSet<Entry<String, List<PIGData>>> s = new TreeSet<Entry<String, List<PIGData>>>(
                new Comparator<Entry<String, List<PIGData>>>(){
           // @Override
            public int compare(Entry<String, List<PIGData>> e1, 
                    Entry<String, List<PIGData>> e2){
                int s1 = e1.getValue().size();
                int s2 =  e2.getValue().size();
                if(s1!=s2){
                    return s1 > s2  ? -1 :1;
                }
                else return e1.getKey().compareTo(e2.getKey());
            }
        });
        s.addAll(left.entrySet());
        return s;
    }
    
    public SortedMap<String, List<String>> getHaplotypes(int st, int end, int mid){
    	final int pos = mid - st;
    	
    	  SortedMap<String, List<String>> left = new TreeMap<String,List<String>>(
                  new Comparator<String>(){

                      public int compare(String o1, String o2) {
                    	  int res1 = -1*o1.substring(pos).compareTo(o2.substring(pos));
                    	  if(res1!=0) return res1;
                         return -1*o1.compareTo(o2);
                      }

                     
                      
                  });
    	  StringBuffer[] sB =new StringBuffer[indiv.size()*Constants.noCopies[this.index]];
          for(int k=0; k<sB.length; k++){
          	sB[k] = new StringBuffer();
          }
          for(int k=st; k<=end; k++){
          	append(k, sB);
          }
         for(int i=0; i<indiv.size(); i++){
        	 int no_copies = dataL.get(indiv.get(i)).noCop();
        	 for(int k = 0; k<no_copies;k++){
	        	 String hap = sB[no_copies*i+k].toString();
	        	 List<String >l = left.get(hap);
	        	 if(l==null) left.put(hap, l = new ArrayList<String>());
	        	 l.add(indiv.get(i));
        	 }
         }
         return left;
    }
    
    public SortedMap<String, List<String>[]> getHaplotypes1( int mid, int phenIndex, Set<String>excl){
    	SortedMap<String, List<String>[]> left = new TreeMap<String, List<String>[]>();
    	Comparable[] compar = new Comparable[indiv.size()];
    	this.getCompa(mid, compar);
         for(int i=0; i<indiv.size(); i++){
        	String key =  indiv.get(i);
        	if(excl.contains(key)) continue;
        	ComparableArray compa = (ComparableArray)compar[i];
        	String st = compa.getGenotypeString();
        	List<String>[] res = left.get(st);
        	if(res==null){
        		left.put(st, res = new List[]{new ArrayList(), new ArrayList(), new ArrayList()});
        		
        	}
        		Double phen = dataL.get(key).phenValue()[phenIndex];
        		if(phen==null){
        			res[2].add(key);
        		}
        		else if(phen > 0.5){
        			res[1].add(key);
        		}
        		else res[0].add(key);
         }
         return left;
    }
    
    public SortedMap<String,List<PIGData> >
        getHaplotypes(int st, int end, Collection<PIGData> set){
        SortedMap<String, List<PIGData>> left = new TreeMap<String,List<PIGData>>(
                new Comparator<String>(){

                    public int compare(String o1, String o2) {
                       return -1*o1.compareTo(o2);
                    }

                   
                    
                });
        StringBuffer[] sB =new StringBuffer[indiv.size()];
        for(int k=0; k<sB.length; k++){
        	sB[k] = new StringBuffer();
        }
        for(int k=st; k<=end; k++){
        	append(k, sB);
        }
       
        for(Iterator<PIGData> it = set.iterator(); it.hasNext();){
            PIGData nxt = it.next();
            String stL = sB[indiv.indexOf(nxt.getName())].toString();//nxt.getStringRep(st, end);
            List<PIGData> cntL = left.get(stL);
            if(cntL==null){
                left.put(stL,cntL = new ArrayList<PIGData>());
            }
            cntL.add(nxt);
        }
      //  EHHFinder.check(left.values());
       return left;
    }
    public void getCompa(int pos, Comparable[] genotypes){
    	if(indiv.size()==0) indiv =indiv();
		for(int i=0; i<indiv.size(); i++){
			genotypes[i] = (data.get(indiv.get(i)).getElement(pos));
			
		}
	}
    public void append(int pos, StringBuffer[] sb){
    	if(indiv.size()==0) indiv = indiv();
    	for(int i=0; i<indiv.size(); i++){
    		
    		ComparableArray comp  = (ComparableArray)(data.get(indiv.get(i)).getElement(pos));
    		int no_copies = comp.size();
    		for(int k=0; k<no_copies; k++){
    			sb[i*no_copies+k].append(comp.get(k).toString());
    		}
    	}
    	
    }
    
    public void insert(List<Info> toInsert){
    	int length_new = toInsert.size();
    
    	int length_old = loc.size();
    	int[] old_index = new int[length_old];  // which indices to map the old positions to
    	int[] new_index = new int[length_new]; // which indices to map the new positions to
    	int j_old=0;
    	int j_new=0;
    	int pos_new = toInsert.get(0).loc;
    	int pos_old = loc.size()==0 ? Integer.MAX_VALUE : loc.get(0);
    	
    	List<Integer> loc_new = new ArrayList<Integer>();
    	List<String> snp_new = new ArrayList<String>();
    	List<Character> alleleA_new = new ArrayList<Character>();
    	List<Character> alleleB_new = new ArrayList<Character>();
    	List<Double> baf_new = new ArrayList<Double>();
    	List<Boolean> probeOnly_new = new ArrayList<Boolean>();
    	List<Boolean > strand_new = new ArrayList<Boolean>();
     	for(int k=0; j_new < length_new || j_old < length_old;k++){
    		if(pos_new==pos_old) throw new RuntimeException("!!");
    		else if(pos_new < pos_old){
    			Info inf = toInsert.get(j_new);
    			loc_new.add(pos_new);
    			snp_new.add(inf.id);
    			alleleA_new.add(inf.alleleA);
    			alleleB_new.add(inf.alleleB);
    			baf_new.add(0.0);
    			probeOnly_new.add(inf.probeOnly);
    			strand_new.add(inf.strand);
    			new_index[j_new] = k;
    			j_new++;
    			pos_new = j_new< length_new ? toInsert.get(j_new).loc : Integer.MAX_VALUE;
    		}else{
    			old_index[j_old] = k;
    			loc_new.add(pos_old);
    			snp_new.add(this.snpid.get(j_old));
    			alleleB_new.add(alleleB.get(j_old));
    			alleleA_new.add(alleleA.get(j_old));
    			baf_new.add(baf.get(j_old));
    			probeOnly_new.add(probeOnly[j_old]);
    			strand_new.add(strand.get(j_old));
    			j_old++;
    			pos_old = j_old < length_old ? loc.get(j_old) : Integer.MAX_VALUE;
    		}
    	}
    	loc = loc_new;
    	this.snpid = snp_new;
    	this.alleleA = alleleA_new;
    	this.alleleB = alleleB_new;
    	this.baf = baf_new;
    	this.probeOnly = probeOnly_new.toArray(new Boolean[0]);
    	this.strand = strand_new;
    	this.length = loc.size();
    	if(this.dc!=null){//remake the distributions
    		this.makeDistributions(this.index);
    	}
    	
    	
    	  for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
              HaplotypeEmissionState data =  (HaplotypeEmissionState) it.next();
 //            PIGData emst = this.data.get(data.getName());
              PseudoDistribution base;
              if(Constants.equaliseGroupMode()!=2){
            	  base = ((CompoundEmissionStateSpace)data.emissions[0].getEmissionStateSpace()).getBaseDist(data.noCop());
              }else{
            	  int ploidy = data.noCop();
            	  int coverage =(int) Constants.coverage(index);
            	 
            	  CompoundEmissionStateSpace emstsp =  Emiss.getSpaceForNoCopies(ploidy);;
            		  //(CompoundEmissionStateSpace) data.getEmissionStateSpace();
            	 // base = (emstsp).getZeroDistProb(coverage, ploidy);
            	  base = //((CompoundEmissionStateSpace)data.emissions[0].getEmissionStateSpace())
            	  emstsp.getHWEDist1(null);
              }
              data.insert(old_index,new_index, base);
      		
             
          }
    	//loc=  makeNew(loc,old_index)
    }
    
    public void reverse(){
    	   //     this.maf.reverse();
    	        Collections.reverse(this.loc);
    	        this.dc.reverse();
    	        for(int i=0; i<loc.size(); i++){
    	            loc.set(i,-1 *loc.get(i));
    	        }
    	        Collections.reverse(this.snpid);
    	        Collections.reverse(this.strand);
    	        Collections.reverse(this.alleleA);
    	        Collections.reverse(this.alleleB);
    	        Collections.reverse(this.baf);
    	        Constants.reverse(this.probeOnly);
    	        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
    	           PIGData nxt = it.next();
    	           if(nxt!=null) nxt.reverse();
    	        }
    	        for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
    	            EmissionState nxt = it.next();
    	            nxt.reverse();
    	        }
    	    }
    public void addCollection(DataCollection datC, int offset){
    	this.dc.addCollection(datC.dc);
    	int offset1 = offset - datC.loc.get(0) + loc.get(loc.size()-1);
    
    	for(int k=0; k< datC.loc.size(); k++){
    		this.loc.add(datC.loc.get(k)+offset1);
    		this.snpid.add(datC.snpid.get(k));
    		if(k<datC.alleleA.size()) this.alleleA.add(datC.alleleA.get(k));
    		if(k<datC.alleleB.size())this.alleleB.add(datC.alleleB.get(k));
    		if(k<datC.baf.size())this.baf.add(datC.baf.get(k));
    		if(k<datC.strand.size()) this.strand.add(datC.strand.get(k));
    	}
    	this.probeOnly = (Boolean[]) Constants.join(this.probeOnly, datC.probeOnly).toArray(new Boolean[0]);//po.toArray(new Boolean[0]);
    	this.length = loc.size();
    	   for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
    		   EmissionState data =  it.next();
    	          // HaplotypeEmissionState dat = (HaplotypeEmissionState) dataL.get(data.getName());
    	          PIGData emst = this.data.get(data.getName());
    	       data.append(datC.dataL.get(data.getName()));
    	       if(emst!=null) ((PhasedDataState)emst).append((PhasedDataState)datC.data.get(data.getName()));
    	   }
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#noCopies(java.lang.String)
     */
    public void drop(List<Integer> toDrop1, boolean setAsMissing) {
    	// if(true) throw new RuntimeException("!!");
    	SortedSet<Integer> toDr = new TreeSet<Integer>(toDrop1);
    	List<Integer> toDrop = new ArrayList<Integer>(toDr);
       // Collections.sort(toDrop);
       List<Boolean>po = probeOnly==null ? null : new ArrayList<Boolean>(Arrays.asList(probeOnly));
        for(int i=toDrop.size()-1; i>=0 ; i--){
            this.loc.remove(toDrop.get(i).intValue());
            if(snpid!=null && snpid.size()>0) this.snpid.remove(toDrop.get(i).intValue());
            if(alleleA!=null && alleleA.size()>0) this.alleleA.remove(toDrop.get(i).intValue());
            if(alleleB!=null && alleleB.size()>0) this.alleleB.remove(toDrop.get(i).intValue());
            if(baf!=null && this.baf.size()>0) this.baf.remove(toDrop.get(i).intValue());
            if(strand!=null && this.strand.size()>0) this.strand.remove(toDrop.get(i).intValue());
            if(this.probeOnly!=null && probeOnly.length>0) po.remove(toDrop.get(i).intValue());
           
        }
        if(this.dc!=null){
        	((DistributionCollection)this.dc).drop(toDrop);
        }
        probeOnly = po==null ? null : po.toArray(new Boolean[0]);
       for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
           EmissionState data =  it.next();
          // HaplotypeEmissionState dat = (HaplotypeEmissionState) dataL.get(data.getName());
          PIGData emst = this.data.get(data.getName());
           if(setAsMissing){
               //data.setAsMissing(toDrop);
               if(emst!=null){
                   emst.setAsMissing(toDrop, Constants.cn_ratio());
                  // dat.setAsMissing(toDrop,Constants.cn_ratio());
               }
           }
           else{
        	  // dat.removeAll(toDrop);
               data.removeAll(toDrop);
               //System.err.println(dat.emissions.length+" "+data.emissions.length);
               //if(dat.emissions.length!=loc.size() || data.emissions.length!=loc.size()) throw new RuntimeException("!!");
               if(emst!=null){
                   emst.removeAll(toDrop);
               }
               
           }
       }
       
       if(!setAsMissing){
          // if(maf!=null) this.maf.removeAll(toDrop);
           this.length=loc.size();
       }
    }
    
    
    
    
    public void append(DataCollection idc){
    	
    	
        for(int i=0; i<idc.loc.size(); i++){
            if(!idc.snpid.get(i).equals(this.snpid.get(i))) throw new RuntimeException("!!");
            if(!idc.loc.get(i).equals(this.loc.get(i))) throw new RuntimeException("!!");
            if(this.alleleA!=null && this.alleleA.size()>i){
                if(!idc.alleleA.get(i).equals(this.alleleA.get(i)))  throw new RuntimeException("!!");
                if(!idc.alleleB.get(i).equals(this.alleleB.get(i)))  throw new RuntimeException("!!");
            }
        }
        if(!idc.loc.equals(this.loc)) throw new RuntimeException("!! not equal");
        for(Iterator<PIGData> it = idc.data.values().iterator(); it.hasNext();){
            PIGData nxt = it.next();
            this.data.put(nxt.getName(), nxt);
        }
        for(Iterator<EmissionState> it = idc.dataL.values().iterator(); it.hasNext();){
            EmissionState nxt = it.next();
            this.dataL.put(nxt.getName(), nxt);
        }
    }
    public void reorder() {
        int[] alias = getAlias(loc);
        apply(alias);
    }
    private void apply(int[] alias) {
       reorder(alias, loc);
        reorder(alias, snpid);
        reorder(alias, alleleA);
        reorder(alias, alleleB);
        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
            PIGData data = it.next();
            EmissionState emst = this.getL(data.getName());
            data.applyAlias(alias);
            emst.applyAlias(alias);
        }
     //  if(maf!=null) maf.applyAlias(alias);
    }

   public static void reorder(int[] alias, List loc) {
      List clone = new ArrayList(loc);
      for(int i=0; i<clone.size(); i++){
          loc.set(alias[i], clone.get(i));
      }
        
    }
   
   public static void reorder(int[] alias, Object[] loc) {
       Object[] clone = new Object[loc.length];
       System.arraycopy(loc, 0, clone, 0, loc.length);
       for(int i=0; i<loc.length; i++){
           loc[alias[i]] = clone[i];
       }
         
     }

    private int[] getAlias(List<Integer> loc) {
        List<Integer> loc1 = new ArrayList<Integer>(loc);
        Collections.sort(loc1);
        int[] res = new int[loc.size()];
        for(int i=0; i<res.length; i++){
            res[i] = loc1.indexOf(loc.get(i));
        }
        return res;
    }

   /* public void fillLikelihoodData(Locreader mid){
        System.err.println("filling in likelihood data");
        for(Iterator<Entry<String, EmissionState>> it = this.dataL.entrySet().iterator(); it.hasNext();){
            Entry<String, EmissionState> key = it.next();
            key.getValue().fillLikelihood(mid, loc);
//            this.dataL.put(key.getKey(), new HaplotypeEmissionState(key.getValue(), mid, loc));
           // it.remove();
        }
    }*/
   public void fillLikelihoodData(boolean probeOnly, double[] v){
	  if(probeOnly){
		  for(int i=0; i<this.probeOnly.length; i++){
			  if(this.probeOnly[i]!=null)this.probeOnly[i] = probeOnly;
		  }
		
	  }
	   for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
		   ((HaplotypeEmissionState)it.next()).fill(probeOnly, v);
	   }
   }
   public void fillLikelihoodData(boolean probeOnly, double[] v, List<Integer> sites){
		  if(probeOnly){
			  for(int i1=0; i1<sites.size(); i1++){
				  int i = sites.get(i1);
				  if(this.probeOnly[i]!=null)this.probeOnly[i] = probeOnly;
			  }
			
		  }
		   for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
			   ((HaplotypeEmissionState)it.next()).fill(probeOnly, v, sites);
		   }
	   }
    public int noCopies(String i) {
      EmissionState st  = this.dataL.get(i);
        return st.noCop();
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#fillLikelihoodData()
     * pseudo is the fraction of prob to split amongst all possibilities
     * keeps allele distribution intact
    
    public void fillLikelihoodData1(double pseudo, int[] indices){
        System.err.println("filling in likelihood data".toUpperCase());
        double p = 0.2; //0
        double q = 0.6; //1
        double r = 0.2; //2
        double[] d = new double[] {p*p, 2*p*q,q*q+2*p*r, 2*q*r, r*r}; //0,1,2,3,4 
     //   Arrays.fill(d, 1.0/(double)d.length);
        int[] se =null;
        if(indices!=null && indices.length>1){
            se = new int[] {
                    this.firstGreaterThan(this.loc, indices[0]),
                     this.firstGreaterThan(loc, indices[1])
            };
            if(loc.contains(indices[1])){
                se[1]--;
            }
        }
        else{
            se = new int[] {0,loc.size()};
        }
       EmissionStateSpace emStSp =  Emiss.getSpaceForNoCopies(Constants.backgroundCount());
       
       double[] init = new double[emStSp.size()];
      
       Arrays.fill(init, 1.0/(double) init.length);
       short di = this.index;
        for(Iterator<Entry<String, PIGData>> it = this.data.entrySet().iterator(); it.hasNext();){
            Entry<String, PIGData> key = it.next();
           HaplotypeEmissionState hes = (HaplotypeEmissionState)  this.dataL.get(key.getKey());
            PhasedDataState emSt = new PhasedDataState(name, this.length(),
            		((HaplotypeEmissionState)key.getValue()).getEmissionStateSpace(),
            		(short) ((HaplotypeEmissionState)key.getValue()).dataIndex());//EmissionState.getEmissionState((PhasedDataState) key.getValue(),false, pseudo, se);
          //  SimpleExtendedDistribution dist = new SimpleExtendedDistribution(init, Double.POSITIVE_INFINITY);
          //  dist.setDataIndex(di);
            this.dataL.put(key.getKey(), emSt);
            for(int i=0; i<emSt.emissions.length; i++){
            	double[] probs = hes.fillCN(i, hes.getEmissionStateSpace(), d);
            	emSt.emissions[i] = new SimpleExtendedDistribution(probs, Double.POSITIVE_INFINITY);
            	emSt.emissions[i].setDataIndex(di);
            }
            //            Arrays.fill(emSt.emissions, dist);
           // if(key.getKey().startsWith("18.5")){
            //    System.err.println(key);
           // }
        }
    } */
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#fillLikelihoodData()
     * pseudo is the fraction of prob to split amongst all possibilities
    
    public void fillLikelihoodData(String[] ids, int[] indices){
        System.err.println("filling in likelihood data".toUpperCase());
        Set<String> toFill = new HashSet<String>(Arrays.asList(ids));
        if(ids.length==1 && (ids[0].equals("1")|| ids[0].equals("all"))){
        	toFill = this.dataL.keySet();
        }
        int[] se =null;
        if(indices!=null && indices.length>1){
            se = new int[] {
                    this.firstGreaterThan(this.loc, indices[0]),
                     this.firstGreaterThan(loc, indices[1])
            };
            if(loc.contains(indices[1])){
                se[1]--;
            }
        }
        else{
            se = new int[] {0,loc.size()};
        }
       EmissionStateSpace emStSp =
    	   
    	   Emiss.getSpaceForNoCopies(Constants.backgroundCount());
       double[] init = new double[emStSp.size()];
      
       Arrays.fill(init, 1.0/(double) init.length);
       short di = this.index;
        for(Iterator<Entry<String, PIGData>> it = this.data.entrySet().iterator(); it.hasNext();){
            Entry<String, PIGData> key = it.next();
            if(toFill.contains(key.getKey())){
	            PhasedDataState emSt = new PhasedDataState(name, this.length(),
	            		((HaplotypeEmissionState)key.getValue()).getEmissionStateSpace(),
	            		(short) ((HaplotypeEmissionState)key.getValue()).dataIndex());//EmissionState.getEmissionState((PhasedDataState) key.getValue(),false, pseudo, se);
	            SimpleExtendedDistribution dist = new SimpleExtendedDistribution(init, Double.POSITIVE_INFINITY);
	            dist.setDataIndex(di);
	            this.dataL.put(key.getKey(), emSt);
	            Arrays.fill(emSt.emissions, dist);
            }
           // if(key.getKey().startsWith("18.5")){
            //    System.err.println(key);
           // }
        }
    } */
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#writeSNPFile(java.io.File, java.lang.String, boolean, java.util.Collection)
     */
    public final void writeSNPFile(File file, String chr, boolean header, Collection<Integer> toD)  throws Exception{
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
       if(header) pw.println("snpid\tchr\t"+Constants.build(index)+"\tA\tB\tstrand");
       for(int i=0; i<loc.size(); i++){
           if(toD!=null && toD.contains(i)) continue;
           String snp_id = this.snpid==null ? i+"" : snpid.get(i);
           pw.print("chr"+chr+"\t");
           pw.print(loc.get(i)+"\t");
           pw.print(loc.get(i)+"\t");
           pw.print(snp_id+"\t");
          if(alleleA.size()>0) pw.print(alleleA.get(i)+"\t");
          if(alleleB.size()>0)  pw.print(alleleB.get(i)+"\t");
           pw.print(this.getTypes(i)+"\t");
         if(strand.size()>0)  pw.println(strand.get(i)+"\t");
         else pw.println();
       }
       pw.close();
        
    }
    protected String getTypes(int i) {
		// TODO Auto-generated method stub
		return this.name;
	}
	/**Used to print original data back out */
    public void printWide(File file1) throws Exception{
    	File file = new File(file1, name);
    	file.mkdir();
    	PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(file, "raw_data_table.txt"))));
    	IlluminaRDistribution dist = (IlluminaRDistribution)((HaplotypeEmissionState)dataL.values().iterator().next()).emissions[0];
    	int colsPerIndiv = dist instanceof IlluminaDistribution ? 2 : 1;
    	int nocols = this.dataL.size()*colsPerIndiv+3;
    	List<String> indiv = this.indiv();
    		pw.print("Name\tChr\tPosition");
    	 	PrintWriter pw_plate = new PrintWriter(new BufferedWriter(new FileWriter(new File(file, "plate.txt"))));
        	PrintWriter pw_pheno = new PrintWriter(new BufferedWriter(new FileWriter(new File(file, "pheno.txt"))));
        	pw_plate.println("PATIENT\tPLATE");
        	pw_pheno.println("PATIENT\tCASE");
    	for(int i=0; i<indiv.size(); i++){
    		if(colsPerIndiv==2) pw.print("\t"+indiv.get(i)+".B Allele Freq");
    		pw.print("\t"+indiv.get(i)+".Log R Ratio");
    		pw_plate.println(indiv.get(i)+"\tplate1");
    		pw_pheno.println(indiv.get(i)+"\t"+Math.round(Math.random()));
    	}
    	pw.println();
    	for(int k=0; k<loc.size(); k++){
    		pw.print(this.snpid.get(k)+"\t"+Constants.chrom0()+"\t"+this.loc.get(k));
    		for(int i=0; i<indiv.size(); i++){
    			IlluminaRDistribution dist1 = (IlluminaRDistribution)((HaplotypeEmissionState)dataL.get(indiv.get(i))).emissions[k];
        		if(colsPerIndiv==2) pw.print("\t"+String.format("%5.3g",((IlluminaDistribution)dist1).b()).trim());
        		pw.print("\t"+String.format("%5.3g",(dist1).r()).trim());
        	}
    		pw.println();
    	}
    
    	pw.close();
    	
    pw_plate.close();
    pw_pheno.close();
    	
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#count(int[], int)
     */
    public int count(int[] indices, int index){
    	int cnt =0;
    	for(int i=0; i<indices.length; i++){
    		if(indices[i]==index) cnt++;
    	}
    	return cnt;
    }
   // public EmissionStateSpace getEmStSpace() {
	//	return this.dataL.values().iterator().next().getEmissionStateSpace();
	//}
	
	
   
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#readCaseControl(java.io.File)
     
    public void readCaseControl(BufferedReader br) throws Exception{
    	if(br!=null){
    		String st = "";
    		while((st = br.readLine())!=null){
    			String[] str = st.split("\\t");
    			String id = str[1];
    			String type = str[2];
    			//if(this.dataL.containsKey(id)){
    				cc.put(id,!type.equals("Desir"));
    			//}
    		}
    		br.close();
    		//if(cc.size()!=this.dataL.size()) throw new RuntimeException("!!!");
    		
    	}
    
    }*/
    
    static int getVal(int i1, EmissionStateSpace emstsp, int type){
    	if(type==0){
			return emstsp.getCN(i1);
		}
		else{
			return (emstsp.getBCount(i1, type));
		}
    	
    	
    }
    static double getAvg(double[] probs, EmissionStateSpace emstsp, int type){
    	//int len = emstsp.defaultList.size();
    	//double[] d = PairEmissionState.pool.getObj(len);
    	if(probs[Constants.getMax(probs)]<Constants.imputedThresh(0)){
    		return Double.NaN;
    	}
    	double val = 0;
		for(int i1=0; i1<probs.length; i1++){
				val+=probs[i1]*getVal(i1, emstsp, type);
		}
		if(val < 1e-3) val=0; 
		return val;
    }
    static void getAvg(double[] probs, EmissionStateSpace emstsp, int type, double[] v){
    	//int len = emstsp.defaultList.size();
    	//double[] d = PairEmissionState.pool.getObj(len);
    	
    	Arrays.fill(v,0.0);
		for(int i1=0; i1<probs.length; i1++){
				v[getVal(i1, emstsp, type)]+=probs[i1];
		}
		
    }
 /* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#getInfo(java.lang.String, java.lang.String, int, boolean)
 */
public String getInfo(String tag, String key, int i, boolean style) throws Exception{
	if(DistributionCollection.dc!=null && 
			DistributionCollection.dc.minQuality(i) !=null &&
			DistributionCollection.dc.minQuality(i)>Constants.qualityThresh()){
	   	  return  "NaN";
	     }
		if(tag.equals("countAll")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			return String.format("%7.3f", getAvg(state.emissions[i], state.getEmissionStateSpace(), 0));
		}
		else if(tag.equals("countB")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			return String.format("%7.3f", getAvg(state.emissions[i], state.getEmissionStateSpace(), 2));
		}
		else if(tag.equals("countA")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			return String.format("%7.3f", getAvg(state.emissions[i], state.getEmissionStateSpace(), 1));
		}
		else if(tag.equals("logR")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			PseudoDistribution dist = state.emissions[i];
			if(dist instanceof IlluminaRDistribution){
				return String.format("%7.3f", ((IlluminaRDistribution)dist).r());
			}
			else return "NaN";
		}
		else if(tag.equals("allele")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			EmissionStateSpace emstsp = state.getEmissionStateSpace();
			//Integer val = state.getFixedInteger(i);
			
			double[] probs = state.emissions[i].probs();
			int i1 = Constants.getMax(probs);
			if(probs[i1]>0.9 && emstsp.getCN(i1)==2){
				return emstsp.get(i1)+"_"+String.format("%5.3g", probs[i1]);
			}
			else return "NaN";
		}
		else if(tag.equals("BAF")){
			HaplotypeEmissionState state = (HaplotypeEmissionState) dataL.get(key);
			PseudoDistribution dist = state.emissions[i];
			if(dist instanceof IlluminaDistribution){
				return String.format("%7.3f", ((IlluminaDistribution)dist).b());
			}
			else return "NaN";
		}
		else if(tag.equals("state")){
			//st_ = Integer.parseInt(tag.substring(5));
			HaplotypeEmissionState state = (HaplotypeEmissionState) this.viterbiL.get(key);
			
		if(state==null) return null;
		//	 ComparableArray arr = (ComparableArray)((PIGData)data.get(key)).getElement(i);
			CompoundEmissionStateSpace emstsp = (CompoundEmissionStateSpace) state.getEmissionStateSpace();
			EmissionStateSpace emstsp1 = emstsp.getMembers()[0];
			double[] count = new double[emstsp1.size()];
			double[] probs = state.emissions[i].probs();
			//double val = 0;
			for(int i1=0; i1<probs.length; i1++){
				//if(probs[i1]>0){
				int[] memInd = emstsp.getMemberIndices(i1);
				for(int k=0; k<memInd.length; k++){
				count[memInd[k]]+=probs[i1];
				}
					
				//}
			}
			for(int j=0; j<count.length; j++){
				if(Math.abs(count[j] - Math.round(count[j]))<Constants.printRoundThresh()) count[j] = Math.round(count[j]);
			}
		//	if(Math.abs(val-3)<0.1)System.err.println(arr+" "+val+" "+Constants.sum(probs));
			return Constants.print(count);
		}
        Object o= ((Map) this.getClass().getField(tag).get(this)).get(key);
        if(o ==null){
            return null;
        }
        else if(o instanceof PIGData){
            ComparableArray arr = (ComparableArray)((PIGData)o).getElement(i);
            StringBuffer sb = new StringBuffer();
            for(int i1=0; i1<arr.size(); i1++){
                sb.append(arr.get(i1));
            }
           String st =  sb.toString();
           return st;
//           if(st.indexOf('_')>=0) return "N";
 //          else return count(st, 'B');
          // return st.replace('_', 'N');
        }
        else if(o instanceof EmissionState){
            return ((EmissionState)o).getUnderlyingData(i);
        }
        else if(o instanceof double[]){
            return String.format("%7.3f", ((double[])o)[i]);//String.format("%7.3f ", new Double[] {});
        }
        
        else return "-";
    }
    
   
    private double getAvg(PseudoDistribution pseudoDistribution,
		EmissionStateSpace emstsp, int type) {
	Integer k = pseudoDistribution.fixedInteger();
	if(k!=null){
		return getVal(k, emstsp, type);
		//this.getAvg(probs, emstsp, type);
	}
	else{
		double v =  this.getAvg(pseudoDistribution.probs(), emstsp, type);
		if(Math.abs(v - Math.round(v))<Constants.printRoundThresh()) v = Math.round(v);
		return v;
	}
	
}
	protected int count(String st, char character) {
    	int cnt =0;
		for(int i=0; i<st.length(); i++){
			if( st.charAt(i)==character) cnt++;
		}
		return cnt;
	}

	private String count(Character character) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#reverse()
     */
   
   
    
    
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getSourcePositions()
     */
    public int[][] getSourcePositions() {
       int[]res = new int[loc.size()];
       for(int i=0; i<res.length; i++){
           res[i] =i;
       }
       return new int[][] {res};
    }
   
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#mix()*/
    
    public void mix(){
        for(Iterator<PIGData> it = data.values().iterator(); it.hasNext();){
           ((HaplotypeEmissionState)it.next()).mix();
        }
    } 
    public PedigreeDataCollection ped;
 
    public static String[] readIndiv(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        List<String> l = new ArrayList<String>();
        String st = "";
        while((st = br.readLine())!=null){
            l.add(st);
        }
        return l.toArray(new String[0]);
            
        
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#makeMafState(lc1.dp.emissionspace.EmissionStateSpace)
     */
/*    public EmissionState makeMafState(EmissionStateSpace emStSp1){
        if(emStSp1 instanceof CompoundEmissionStateSpace){
            CompoundEmissionStateSpace emStSp = (CompoundEmissionStateSpace) emStSp1;
        EmissionStateSpace[] ems = emStSp.getMembers();
        EmissionState[] st = new EmissionState[ems.length];
        for(int i=0; i<st.length; i++){
            st[i] = new HaplotypeEmissionState("maf_"+i, length, ems[i].size(), ems[i], null, null);
        }
        
        this.maf =  new AlleleCopyPairEmissionState(Arrays.asList(st), emStSp, false, null);// emStSp.size(), emStSp);
        }
        else{
           maf =   new HaplotypeEmissionState("maf_", length, emStSp1.size(), emStSp1,  null, null);
        }
       
        maf.initialiseCounts();
        return maf;
    }*/
    
    protected void makeBackgroundState() {
    	if(!Constants.modelbg()) return;
    	  EmissionStateSpace emStSp = null;//this.stSp[no_copies-1];//Emiss.getEmissionStateSpace(1, 3);
    	  short di = this.index;
    	 bg =   new BackgroundEmissionState("!bg", length, emStSp.size(), emStSp, this.index);
    	 bg1 = SimpleScorableObject.make("!bg",loc.size(), emStSp, this.index);
    	    //	
    	// BackgroundEmissionState bg1  = SimpleScorableObject.make(bg,Arrays.asList(new String[] {subst, subst1}), emStSp, (short)-1);
    	//this.cn_count = new double[emStSp.copyNumber.size()];
    	if(Constants.fixedBG()){
    		 for(int i=0; i<bg.emissions.length; i++){
    			 bg.emissions[i] = new IntegerDistribution(emStSp.get(ComparableArray.make(Emiss.A, Emiss.A)), emStSp);
              	bg.emissions[i] = new IntegerDistribution(emStSp.get(ComparableArray.make(Emiss.A, Emiss.A)),emStSp);
              	bg.emissions[i].setDataIndex(di);
              }
    	}
    	else{
    	double q = Constants.bg(); //1
    	double[] d;
    	if(q>0){
    	double p = (1-q)/2.0; //0
         
         double r = 1-q-p; //2
          d = new double[]
        //	 new double[5];//
         {p*p, 2*p*q,q*q+2*p*r, 2*q*r, r*r}; //0,1,2,3,4 
    	}
    	else{
    		d = new double[5];
    		 Arrays.fill(d,(1+q)/(double)(d.length));
    	      d[2] -=q;
    	}
      System.err.println("background start is "+d[0]+" "+d[1]+" "+d[2]+" "+d[3]+" "+d[4]);
    
     
       double[] probs = bg.fillCN(emStSp,  d);
        
             for(int i=0; i<bg.emissions.length; i++){
            	 bg.emissions[i] = new IntegerDistribution(emStSp.get(ComparableArray.make(Emiss.A, Emiss.A)),emStSp);
             	bg.emissions[i] = new SimpleExtendedDistribution1(probs, Double.POSITIVE_INFINITY);
             	bg.emissions[i].setDataIndex(di);
             }
    	}
    	bg.initialiseCounts();
    	dataL.put(bg.name, bg);
    	data.put(bg.name, bg1);
    	indiv.add(bg.name);
           
	}
	public static void main(String[] args){
        String input_f = args[0];
        String output_f = args[1];
        String input_file = args[2];
        String output_file = args[3];
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#clone()
     */
    public DataC clone(){
    	return new DataCollection(this);
    }
    public static  List readPosInfo(File f, int index, boolean header, Class clazz) throws Exception{
    	 List res = new ArrayList();
         readPosInfo(f, new int[] {index}, header, new List[] {res}, new Class[] {clazz});
          return res;
    }
    public static  List<Integer> readPosInfo(File f, int index, boolean header) throws Exception{
        
        List<Integer> res = new ArrayList<Integer>();
       readPosInfo(f, new int[] {index}, header, new List[] {res}, new Class[] {Integer.class});
        return res;
    }
    public static BufferedReader getBufferedReader(File f) throws Exception{
        if(f.exists() && f.length()>0){
            return new BufferedReader(new FileReader(f));
        }
        else{
            File f1 = new File(f.getAbsolutePath()+".gz");
            if(f1.exists())
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f1))));
            else return null;
        }
    }
    
    public static  void readPosInfo(File f, int[] index, boolean header, List[] res, Class[] cl) throws Exception{
        if(!f.exists() || f.length()==0) return;
        readPosInfo(getBufferedReader(f), index, header, res, cl);
    }
    
    public static  void readPosInfo(BufferedReader br, int[] index, boolean header, List[] res, Class[] cl) throws Exception{
        readPosInfo(br, index, header, res, cl, "\\s+");
    }
    public static  void readPosInfo(BufferedReader br, int[] index, boolean header, List[] res, Class[] cl, String spl) throws Exception{
        if(header) br.readLine();
        String st = "";
      
        while((st = br.readLine())!=null){
            String st1 = st.trim();
            String[] str = st1.split(spl);
         //   System.err.println(st1);
            for(int i=0; i<index.length; i++){
                if(str.length>index[i]){
                    try{
                        if(cl[i].equals(String.class)){
                            res[i].add(str[index[i]]);
                        }
                        else{
                            res[i].add(
                                    cl[i].getConstructor(new Class[] {String.class}).newInstance(new Object[] {str[index[i]]}));
                        }
                    }catch(Exception exc){
                   //     System.err.println(Arrays.asList(str));;
                        exc.printStackTrace();
                        System.exit(0);
                    }
                }
            }
        }
        br.close();
    }
    public static SimpleDataCollection getMarchiniFile(File parent, int i, boolean trio, boolean answer, int data_index){
        try{
            File[] f;
            if(answer){
                f = new File[] {new File(parent, "ans.haps."+i)};
            }
            else if(trio){
                f = new File[] {new File(parent, "pgenos.haps."+i), new File(parent, "cgenos.haps."+i)};
            }
            else{
                f = new File[] {new File(parent, "pgenos.haps."+i)};
            }
            SimpleDataCollection coll = readMarchiniFormat(f, new File(parent, "posinfo."+i), data_index);
            coll.name = parent.getName();
            coll.chrom = "MAR";
       //     coll.calculateMaf(false);
            for(int i1=0; i1<coll.loc.size(); i1++){
                coll.snpid.add(i1+"");
            }
            if(answer) coll.setPedigree(null);
            return coll;
        }catch(Exception exc){
            exc.printStackTrace();
            return null;
        }
    }
    
    public static void download(String address, String localFileName) {
        OutputStream out = null;
        URLConnection conn = null;
        InputStream  in = null;
        try {
            URL url = new URL(address);
            out = new BufferedOutputStream(
                new FileOutputStream(localFileName));
            conn = url.openConnection();
            in = conn.getInputStream();
            byte[] buffer = new byte[1024];
            int numRead;
            long numWritten = 0;
            while ((numRead = in.read(buffer)) != -1) {
                out.write(buffer, 0, numRead);
                numWritten += numRead;
            }
            System.out.println(localFileName + "\t" + numWritten);
        } catch (Exception exception) {
            exception.printStackTrace();
        } finally {
            try {
                if (in != null) {
                    in.close();
                }
                if (out != null) {
                    out.close();
                }
            } catch (IOException ioe) {
            }
        }
    }

    
   /* public static DataCollection readHapMapPhasedFormat(File dir, String prefix, double soften,
            double probPair, int[] mid, int[] kb) throws Exception{
        boolean trio = true;
        File[] f = new File[] {
           new File(dir, prefix+"legend.txt"),
           new File(dir, prefix+"phased"),  
           new File(dir, prefix+"sample.txt"),  
        };
      
        List<Integer> loc1 = new ArrayList<Integer>();
        List<String> snps = new ArrayList<String>();
        List<String> major = new ArrayList<String>();
        List<String> minor =  new ArrayList<String>();
        SimpleDataCollection.readPosInfo(f[0], new int[] {0,1, 2, 3}, true, new List[] {snps, loc1, major, minor}, new Class[] {String.class, Integer.class, String.class, String.class});
        
        List<String> ids = new ArrayList<String>();
        SimpleDataCollection.readPosInfo(new File(dir, prefix+"sample.txt"), new int[] {0}, false, new List[] {ids}, new Class[] {String.class});
        if(trio){
            int noFam = ids.size();
            int noParents = (int) Math.round((double)noFam*(2.0/3.0));
            if(noParents*(3.0/2.0)!=noFam) throw new RuntimeException("!!");
            ids = ids.subList(0, noParents);
        }
        int read; int end;
       if(Constants.mid()[0][0]>0 || Constants.mid().length>1){
           int mid_index;
           int mid_index1;
           if(kb.length>1){
               int start_index = mid[0] - (int) kb[0]*1000;
                mid_index = firstGreaterThan(loc1, start_index);
                mid_index1 = firstGreaterThan(loc1, (int) mid[1]+(int) kb[1]*1000);
           }
           else{
               mid_index = firstGreaterThan(loc1, mid[0]) - Constants.restrict()[0];
               mid_index1 = firstGreaterThan(loc1, mid[1])+Constants.restrict()[1];
           }
           read = mid_index;
           end  = mid_index1 ;
       }
       else{ 
           read = IlluminaRDataCollection.firstGreaterThan(loc1,  Constants.offset());
           end  = Math.min(getEndIndex(loc1), read+Constants.restrict()[0]);
         
       }
       read = Math.max(0,read);
       end =  Math.min(end,loc1.size());
       DataCollection sdc = new LikelihoodDataCollection();
       sdc.name = "hapmap";
        sdc.loc = new ArrayList<Integer>(loc1.subList(read, end)) ;
        sdc.snpid = new ArrayList<String>(snps.subList(read, end)) ;
        sdc.majorAllele = new ArrayList<Character>();
        sdc.alleleB = new ArrayList<Character>();
        for(int i=read; i<end; i++){
            if(loc1.get(i).intValue()!=sdc.loc.get(i-read).intValue()) throw new RuntimeException("!!");
            char allA = major.get(i).charAt(0);
            char allB = minor.get(i).charAt(0);
            sdc.majorAllele.add(allA);
            sdc.alleleB.add(allB);
        }
       // System.err.println(sdc.snpid);
      //  System.err.println(sdc.loc);
       // System.err.println(sdc.majorAllele);
      //  System.err.println(sdc.alleleB);
        BufferedReader br = DataCollection.getBufferedReader(f[1]);
        String st = ""; String st1 = "";
        sdc.length = sdc.loc.size();
        int i=0;
       
        for( i=0; (st = br.readLine())!=null; i++){
            String subst = st.replaceAll("\\s+", "").substring(read, end).replace('0','A').replace('1','B');
            String subst1 = br.readLine().replaceAll("\\s+", "").substring(read, end).replace('0','A').replace('1','B');
            EmissionStateSpace emStSp = (CompoundEmissionStateSpace) Emiss.getEmissionStateSpace(2-1);
            if(subst.length()!=sdc.length() || subst1.length()!=sdc.length()) throw new RuntimeException("!!");
            if(probPair>0){
               
                PIGData dat = SimpleScorableObject.make(ids.get(i),Arrays.asList(new String[] {subst, subst1}), emStSp, (short)-1);
                if(dat.length()!=sdc.length()) throw new RuntimeException("!!");
                sdc.data.put(dat.getName(), dat);
            }
            else{
                PIGData dat = SimpleScorableObject.make(ids.get(i)+"_"+0,Arrays.asList(new String[] {subst}), emStSp,(short)-1);
                if(dat.length()!=sdc.length()) throw new RuntimeException("!!");
                sdc.data.put(dat.getName(), dat);
                PIGData dat1 = SimpleScorableObject.make(ids.get(i)+"_"+1,Arrays.asList(new String[] {subst1}),  emStSp,(short)-1);
                if(dat1.length()!=sdc.length()) throw new RuntimeException("!!");
                sdc.data.put(dat1.getName(), dat1);
            }
        }
     if(probPair<1.0) sdc.split(probPair, true);
      //  sdc.fillLikelihoodData(soften, Constants.specialTrans());
        sdc.calculateMaf(false);
     if(Constants.topBottom())   sdc.convertToTopBottom();
        return sdc;
    }*/
    public DataCollection[] splitOnSize(){
        EmissionStateSpace[] m = this.getNoCopies();
        DataCollection[] coll = new DataCollection[m.length];
        for(int i=0; i<coll.length; i++){
            coll[i] = new LikelihoodDataCollection();
            coll[i].loc = new ArrayList<Integer>(this.loc);
             coll[i].snpid = new ArrayList<String>(this.snpid);
            coll[i].alleleA = new ArrayList<Character>(this.alleleA);
            coll[i].alleleB = new ArrayList<Character>(this.alleleB);
            coll[i].length = this.length;
        }
        for(Iterator<String> it = this.getKeyIterator(); it.hasNext();){
            String key = it.next();
            PIGData dat = this.get(key);
            EmissionState emst = this.getL(key);
            int no_cop = dat.noCopies();
            coll[no_cop-1].data.put(key, dat);
            coll[no_cop-1].dataL.put(key, emst);
        }
     //   for(int i=0; i<coll.length; i++){
      //      coll[i].calculateMaf(false);
       // }
        return coll;
       // for(Iterator<Integer> it = )
    }
  
    public int restricToAlias(Collection<String> alias) {
     //   boolean b = alias.contains("22086");
        List<String> keys = new ArrayList<String>(this.getKeys());
        
        for(Iterator<String> it = keys.iterator(); it.hasNext();){
            String key = it.next();
            if(!alias.contains(key)){
               dataL.remove(key);
             recSites.remove(key);
               viterbi.remove(key);
             data.remove(key);
           if(indiv!=null)  indiv.remove(key);
            }
        }
        this.indiv=null;
      return  dataL.size();
     //   System.err.println("sze is "+sze);
        
    }
    public DataCollection[] splitRandom(){
        DataCollection[] coll = new DataCollection[2];
        for(int i=0; i<coll.length; i++){
            coll[i] = new LikelihoodDataCollection();
            coll[i].loc = new ArrayList<Integer>(this.loc);
             coll[i].snpid = new ArrayList<String>(this.snpid);
            coll[i].alleleA = new ArrayList<Character>(this.alleleA);
            coll[i].alleleB = new ArrayList<Character>(this.alleleB);
            coll[i].length = this.length;
        }
        for(Iterator<String> it = this.getKeyIterator(); it.hasNext();){
            String key = it.next();
            PIGData dat = this.get(key);
            EmissionState emst = this.getL(key);
            int no_cop = Constants.rand.nextBoolean() ? 1:2;
            coll[no_cop-1].data.put(key, dat);
            coll[no_cop-1].dataL.put(key, emst);
        }
       // for(int i=0; i<coll.length; i++){
         //   coll[i].calculateMaf(false);
       // }
        return coll;
       // for(Iterator<Integer> it = )
    }
    /* reads marchini format */
    public static SimpleDataCollection readMarchiniFormat(File[] f, File posF, int index) throws Exception{
        List<Integer> locs = readPosInfo(posF,2, false);
        int noSites = Math.min(Constants.restrict()[0],locs.size());
        SimpleDataCollection coll = new SimpleDataCollection(noSites);
        coll.index = (short)index;
        coll.loc = locs.subList(0, noSites);
        //coll.snpid = new ArrayList<String>()
        for(int ik=0; ik<f.length; ik++){
        BufferedReader br = new BufferedReader(new FileReader(f[ik]));
        String st = "";
        //br.readLine();
        EmissionStateSpace emStSp = (CompoundEmissionStateSpace) Emiss.getEmissionStateSpace(1);
        while((st = br.readLine())!=null ){
            if(st.length()==0) break;
            String[] str1 = st.split("\\s+");
            String[] str2 = br.readLine().split("\\s+");
            PhasedDataState dat1 =  SimpleScorableObject.make(str1[3].replace(':', '_'),noSites, emStSp,(short)-1);
            dat1.setDataIndex(index);
            for(int j=0; j<noSites; j++){
               Comparable c1 = translate(str1[0].charAt(j), true);
               Comparable  c2 = translate(str2[0].charAt(j), false);
               
                    dat1.addPoint(j,ComparableArray.make(c1, c2));
             
            }
        //    System.err.println(dat1.dataIndex());
            coll.data.put(dat1.getName(), dat1);
        }
        br.close();
        }
        Map<String, String> mother = new TreeMap<String, String>();
        Map<String, String> father= new TreeMap<String, String>();
       for(int i=1; ;i++){
           String child_name = "FAM"+i+"_CHILD";
           if(!coll.data.containsKey(child_name)) {
               break;
           }
           String moth_name = "FAM"+i+"_MOTH";
           String fath_name = "FAM"+i+"_FATH";
           mother.put(child_name, moth_name);
           father.put(child_name, fath_name);
       }
        coll.setPedigree(new PedigreeDataCollection(mother, father));
        return coll;
    }
    
    private static Comparable translate(char c, boolean first) {
        if(c=='0') return Emiss.N();
        else if(c=='1') return Emiss.b();
        else if(c=='2')return Emiss.a();
        else if(c=='9'){
            if(first) return Emiss.a();
            else return Emiss.b();
        }
        else throw new RuntimeException("!!");
    }
    /** reads dick format */
    public static SimpleDataCollection readDickFormat(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        String[] st = br.readLine().split("\\s+");
        int noIndiv = Integer.parseInt(st[1]);
        int noSites = Math.min(Constants.restrict()[0], Integer.parseInt(st[2]));
        SimpleDataCollection coll = new SimpleDataCollection(noSites);
        Character[] majorAllelle = new Character[noSites];  //Emiss.A
        Character[]minorAllelle = new Character[noSites];   //Emiss.B
        st = br.readLine().split("\\s+"); 
        for(int i=1; i<=noSites; i++){
            coll.loc.add(Integer.parseInt(st[i]));
        }
        st = br.readLine().split("\\s+"); 
        for(int i=1; i<=noSites; i++){
            coll.snpid.add(st[i]);
        }
        for(int i=st.length; i<noSites; i++){
            coll.snpid.add("-");
        }
        for(int i=0; i<noIndiv; i++){
            String[] str = br.readLine().split("\\s+");
            EmissionStateSpace emStSp = (CompoundEmissionStateSpace) Emiss.getEmissionStateSpace(1);
            
            PIGData dat =  SimpleScorableObject.make(str[0],noSites, emStSp,(short)-1);
            for(int j=0; j<noSites; j++){
                char c = str[j+1].charAt(0);
                if(c=='?' || c=='-') dat.addPoint(j, ComparableArray.make(Emiss.N()));
                else{
                    if(majorAllelle[j]==null){
                        majorAllelle[j] = new Character(c);
                        dat.addPoint(j,ComparableArray.make(Emiss.a()));
                    }
                    else if(minorAllelle[j]==null && c!=majorAllelle[j]){
                        minorAllelle[j] = c;
                        dat.addPoint(j,ComparableArray.make(Emiss.b()));
                    }
                    else{
                        dat.addPoint(j,ComparableArray.make(c==majorAllelle[j] ? Emiss.a() : Emiss.b()));
                    }
                }
            }
            coll.data.put(dat.getName(), dat);
        }
      
        return coll;
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#noAllelles()
     */
    public int noAllelles(){
        int res =0;
        for(Iterator<String> it = this.data.keySet().iterator(); it.hasNext();){
            res+=this.noCopies(it.next());
        }
        return res;
    }
    public void switchAlleles(int i){
        Character maj = this.alleleA.get(i);
        Character min = this.alleleB.get(i);
        alleleA.set(i, min);
        alleleB.set(i, maj);
 //    maf.switchAlleles(i);
        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
            it.next().switchAlleles(i);
        }
        for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
        	EmissionState st = it.next();
        	st.switchAlleles(i);
          // System.err.println("h");
        }
//        EmissionState st1 = this.dataL.get("NA18561");
 //       System.err.println("h");
    }
    
    public void flipStrand(int i){
    	if(alleleA.size()>0){
    		{
    			alleleA.set(i, compl(this.alleleA.get(i)));
    			alleleB.set(i, compl(this.alleleB.get(i)));
    		}
    	}
    	// this.baf.set(i, 1-this.baf.get(i));
    	/* Boolean strand_ = strand.get(i);
    	 if(strand_!=null && strand_) {
    		 throw new RuntimeException("should only flip to plus strand");
    	 }*/
    	
    	 if(strand.size()>0) strand.set(i, true);
    }
    
    public void swapAlleles(int i){
    	//if(alleleA.size()>0){
    		{
    			if(alleleA.size()>0 && alleleB.size()>0){
    			Character tmp = alleleA.get(i);
    			alleleA.set(i, alleleB.get(i));
    			alleleB.set(i, tmp);
    			}
    			this.baf.set(i, 1-this.baf.get(i));
    			for(Iterator<EmissionState> st = this.dataLvalues(); st.hasNext();){
    				((HaplotypeEmissionState)st.next()).switchAlleles(i);
    			}
    			strand.set(i, !strand.get(i));
    		}
    	//}
    	// this.baf.set(i, 1-this.baf.get(i));
    	/* Boolean strand_ = strand.get(i);
    	 if(strand_!=null && strand_) {
    		 throw new RuntimeException("should only flip to plus strand");
    	 }*/
    }
    
    public static Character compl(char ch){
    	if(ch=='A') return 'T';
    	else if(ch=='T') return 'A';
    	else if(ch=='G') return 'C';
    	else if(ch=='C') return 'G';
    	else if(ch=='N') return 'N';
    	else if(ch=='-') return '-';
    	else if(ch=='I') return 'I';
    	else if(ch=='D') return 'D';
    	else {
    		throw new RuntimeException("!!"+ch);
    	}
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#extractFromTrioData()
     */
    public void extractFromTrioData(){
     //    ped = new PedigreeDataCollection();
        Map<String, PIGData> l = new HashMap<String, PIGData>();
        Map<String, EmissionState> data_l = new HashMap<String, EmissionState>();
        Map<String, double[]> unc = new HashMap<String, double[]>();
        for(Iterator<PIGData> it = data.values().iterator(); it.hasNext();){
            PIGData dat_i =it.next();
            EmissionState dat_i1 = this.dataL.get(dat_i.getName());
            if(((ComparableArray)dat_i.getElement(0)).isNested()){
                PIGData[] data_i = dat_i.split();
                EmissionState[] data_il = dat_i1==null ? null : dat_i1.split();
                for(int j=0; j<data_i.length; j++){
                    l.put(data_i[j].getName(), data_i[j]);
                    if(dat_i1!=null) data_l.put(data_il[j].getName(), data_il[j]);
                    
                }
              //  ped.setTrio(data_i[0].getName(), data_i[1].getName(), data_i[2].getName());
             
            }
            else{
                l.put(dat_i.getName(), dat_i);
                if(dat_i1!=null) data_l.put(dat_i1.getName(), dat_i1);
                
            }
        }
        this.data = l;
        this.dataL = data_l;
       
       // this.makeDataIndex();
    }
    
  /*  public void calculateMaf(boolean state){
    	//if(false){
    	this.makeBackgroundState();
            if(!state) calculateMaf(state, true);
            else  calculateMaf(state, false);
            
    //	}
    }*/
    public  Phenotypes pheno(){
        return this.pheno;
    }
    
   
    public PrintWriter log;
  //This makes each data state a mixture of a 'null' hwe distribution, as well actual distribution
   public void addMixture(){
	for(Iterator<EmissionState> it  = this.dataLvalues(); it.hasNext();){
	    HaplotypeEmissionState  st =(HaplotypeEmissionState) it.next();
	    this.dc.addMixture(st, this.probeOnly);
	}
	  
  }
    
    public String headSNP(){
    	for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
 	  HaplotypeEmissionState st =(HaplotypeEmissionState)it.next();
 	  if(st.getEmissionStateSpace()!=null){
     //   HaplotypeEmissionState st1 = (HaplotypeEmissionState) dataL.values().iterator().next();
    String headSNP = st==null ? "" :  st.getCompressedStringHeader();
        if(st!=null) {
            headSNP= headSNP+"\t"+"distr";
            //+Emiss.getSpaceForNoCopies(st.noCop()).getHeader();
        }
    	
        return headSNP;
 	  }
    	}
    	return "";
    }
    public String head_snp(){
    	String head_snp = "chr\tstart\tend\tid";
        if(alleleA!=null && alleleA.size()>0){
            head_snp = head_snp+"\tA\tB";
        }
        return head_snp;
    }
 
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#calculateMaf()
   
    public void calculateMaf(boolean state, boolean exclMissing){
       CompoundEmissionStateSpace emStSp=Emiss.getSpaceForNoCopies(Constants.backgroundCount()); 
       DistributionCollection.dc = this.dc;
       DataCollection.datC = this;
      EmissionStateSpace emStSp_1 = emStSp.getMembers()[0];
        this.maf =  this.makeMafState(emStSp_1);
        double[] dist = new double[emStSp.defaultList.size()];
        outer: for(int i=0; i<maf.noSnps();  i++){
            if(!state){
                for(Iterator<PIGData> it = data.values().iterator(); it.hasNext();){
                    PhasedDataState nxt = (PhasedDataState) it.next();
                    CompoundEmissionStateSpace emstsp = (CompoundEmissionStateSpace) Emiss.getEmissionStateSpace(nxt.noCopies()-1);
                    Comparable compa = nxt.getElement(i);
                  
                    if(exclMissing){ 
                        if(compa instanceof ComparableArray ){
                            int cnt = ((ComparableArray)compa).size();
                            int cn = ((ComparableArray)compa).copyNumber();
                            if(cn!=cnt) continue outer;
                            int cntNull = ((ComparableArray)compa).countNull();
                            if(cnt==cntNull) continue outer;
                        }
                        else{
                            if(((Emiss)compa)==Emiss.N()) continue outer;
                        }
                    }
                    Integer ind = emstsp.get(compa);
                    int[] indices = emstsp.getMemberIndices(ind);
                        for(int k=0; k<indices.length; k++){
                            if(indices[k]!=0 && indices[k]!=1) {
                               // System.err.println(compa);
                            }
                            maf.addCount(indices[k], 1.0, i);
                        }
                }
            }
            else{
                for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
                    EmissionState nxt = it.next();
                    CompoundEmissionStateSpace emStsp = (CompoundEmissionStateSpace)nxt.getEmissionStateSpace();
                    Integer fixed = nxt.getFixedInteger(i);
                    if(fixed!=null){
                        int ind = nxt.getFixedInteger(i);
                        int[] indices = emStsp.getMemberIndices(ind);
                            for(int k=0; k<indices.length; k++){
                                maf.addCount(indices[k],  1.0, i);
                            }
                    }
                    else{
                    	
                        //double[] prob = nxt.getEmiss(i);
                    	Arrays.fill(dist,0);
                        double sum = ((HaplotypeEmissionState)nxt).calcDistribution(i, dist,emStSp);
               
		                        for(int ind =0; ind<dist.length; ind++){
		                            int[] indices =emStSp.getMemberIndices(ind);
		                            for(int k=0; k<indices.length; k++){
		                              
		                                maf.addCount(indices[k], dist[ind]/sum, i);
		                            }
		                        }
                      
                    }
                    
                }
            }
            
        }
        this.maf.transferCountsToProbs(0.0);
        this.maf.initialiseCounts();
    }  */
    
   
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#calculateMaf1()
     
    public EmissionState calculateMaf1(){
       EmissionStateSpace[] m = this.getNoCopies();
        EmissionStateSpace emStSp_1 = null;
        for(int i=0; i<m.length; i++){
            if(m[i]!=null){
                emStSp_1 = ((CompoundEmissionStateSpace)m[i]).getMembers()[0];
            }
        }
  //     (CompoundEmissionStateSpace) emSt
       EmissionState maf =  this.makeMafState(emStSp_1);
         for(int i=0; i<length();  i++){
             inner: for(Iterator<PIGData> it = data.values().iterator(); it.hasNext();){
                 PhasedDataState nxt = (PhasedDataState) it.next();
                 Comparable compa = nxt.getElement(i);
                 int no_cop = nxt.noCopies()-1;
                 CompoundEmissionStateSpace emStsp = (CompoundEmissionStateSpace) m[no_cop];
                 Integer ind = emStsp.get(compa);
                 if(ind==null) continue inner;
               //  maf.addCount(ind, 1.0, i);
                 int[] indices = emStsp.getMemberIndices(ind);
                     for(int k=0; k<indices.length; k++){
                         maf.addCount(indices[k], 1.0, i);
                     }
             }
             
         }
        
         maf.transferCountsToProbs(0);
         maf.initialiseCounts();
         return maf;
     }*/
    
   static  class FastphaseFormatReader{
        String currentString;
        BufferedReader br;
        Class clazz;
        String st;
        final int ploidy;
        EmissionStateSpace emStSp;
        SimpleDataCollection coll =new SimpleDataCollection();
        int i=0;
        FastphaseFormatReader(BufferedReader br, Class clazz, EmissionStateSpace emStSp, int ploidy) throws Exception{
            this.br = br;
            coll.chrom = Constants.chrom0();
            this.clazz = clazz;
            this.ploidy = ploidy;
            this.emStSp = emStSp;
           // st = br.readLine();
        }
        public PIGData readSingleFastPhaseLine(String name, EmissionStateSpace em, int ploidy) throws Exception{
            int noCopies = 2;
            List<String>[] dat1 =new List[ploidy];
            for(int i=0; i<dat1.length; i++){
            	dat1[i] = new ArrayList<String>();
            }
            EmissionStateSpace emStSp =  Emiss.getEmissionStateSpace(noCopies-1); //source
            EmissionStateSpace emStSp1 = Emiss.getEmissionStateSpace(noCopies-1); //target
            EmissionStateSpaceTranslation trans = new EmissionStateSpaceTranslation(emStSp, emStSp1, true);
            for(int k=0; k<ploidy && st!=null; ){
		            if(st!=null &&!st.startsWith("#") && !st.startsWith("END")  &&!st.startsWith("certainty") && !st.startsWith(">") ){
		               // st =  st.replaceAll("\\s+", "");
		                dat1[k].add(st);//new String(st.replaceAll("\\s+","")));
		               
		                if(st!=null) st = st.replaceAll("\\s+","");
		               
		                //System.err.println(st);
		            }
		            else if(st!=null && st.startsWith("#")){
		            	k++;
		            	//if(k==ploidy) break;
		            }
		           if(k<ploidy){
		        	   st = br.readLine();
		           }
            }
         /*  if(dat1.size()<=2){
                return SimpleScorableObject.make(name, dat1,  em);
            }*/
          //  else if(dat1.size()==3){*/
           
                int len = dat1[0].get(0).length();
                PhasedDataState dat =  new PhasedDataState(name.trim(), dat1[0].get(0).length(), em, (short)-1);
              dat.setNoCop(ploidy);
              Comparable[] compa = new Comparable[ploidy];
               // System.err.println(dat.getName()+" "+len);
                for(int i=0; i<len; i++){
                	Arrays.fill(compa, null);
                    StringBuffer sb = new StringBuffer();
                    for(int j=0; j<ploidy; j++){
                    	StringBuffer sb1 = new StringBuffer();
                    	for(int k=0; k<dat1[j].size(); k++){
	                        char ch = dat1[j].get(k).charAt(i);
	                        if(ch!=' ' && ch!='-' && ch!='_')sb1.append(ch);
	                        if(k<dat1[j].size()-1) sb1.append(",");
                    	}
                    	String st1 =  sb1.toString();
                    	char[] ch1 = st1.toCharArray();
                    	if(Constants.phaseInner()){
                    	String st2 = st1.replaceAll(",", "");
                    	if(st2.length()>1){
                    		
	                    	EmissionStateSpace inn =  ((CompoundEmissionStateSpace)  em).getMembers()[j];
	                    	String st1_ = new String(ch1);
	                    	Integer comp1 = inn.getHapl(st1_);
	                    	if(comp1==null){
	                    		throw new RuntimeException("prob with hap");
	                    	}
	                    	compa[i] =inn.get(comp1.intValue());
                    	}
                    	}
                    	Arrays.sort(ch1);
                    	for(int k=0; k<ch1.length; k++){
                    		if(ch1[k]!=','){
                    			sb.append(ch1[k]);
                    		}
                    	}
                    	if(j<dat1.length-1) {
                    		sb.append(",");
                    	}
                    }
                    char[] ch = sb.toString().toCharArray();
                 //   Arrays.sort(ch);
                    Integer comp =em.getHapl(new String(ch));
                    if(comp==null){
                    	throw new RuntimeException("prob with hap");
                    }
                    ((PhasedDataState)dat).emissions[i] = new IntegerDistribution(comp,em);
                    ComparableArray compa_ = (ComparableArray) em.get(comp.intValue());
                    if(Constants.phaseInner() && compa_.noCopies(true)>ploidy){
                    	ComparableArray compa1=  compa_.copy();
                    	for(int j=0; j<compa.length; j++){
                    		if(compa[j]!=null){
                    			compa1.set(j, compa[j]);
                    		}
                    	}
                    }
                    
                    ((PhasedDataState)dat).phased[i] =  compa_;
                  //  if(dat.getName().startsWith("21939") && (i==224 || i==225)){
                  //      System.err.println(dat.getName()+" "+i+" "+comp);
                  //  }
                //    dat.addPoint( i, emStSp1.get(comp.intValue()));
                }
                return dat;
           // }
        }
        public  void readLocation() throws Exception{
            coll.loc = new ArrayList<Integer>();
            while((st = br.readLine())!=null && !st.startsWith("#")){
                String[] str = st.trim().split("\\s+");
                for(int i=0; i<str.length; i++){
                    if(str[i].length()==0) continue;
                    coll.loc.add(Integer.parseInt(str[i]));
                }
            }
           Collections.sort(coll.loc);
           st = br.readLine();
        }
        public void readMultiLine() throws Exception{
            PIGData dat = null;
            int index = st.indexOf("id");
            //String[] str = st.split("\\s+");
            String name = new String(st);
            if(index>=0)
                name = name.substring(index+2);
           
            if(name.indexOf('|')>=0) name = name.substring(0, name.indexOf('|'));
            name = new String(name.trim().toCharArray());
            st = br.readLine();
            while(st.startsWith("#")){
                st = br.readLine();
            }
            if(st.startsWith(">")){
                boolean trCompArray = st.endsWith("true");
                List<PIGData> l = new ArrayList<PIGData>();
                for(int i=0; st!=null && st.startsWith(">"); i++){
                    st = br.readLine();
                    l.add(readSingleFastPhaseLine(name+i, emStSp,ploidy));
                }
                dat =  SimpleScorableObject.make(l.toArray(new PhasedDataState[0]), false , "_", this.emStSp);
                if(trCompArray){
                    dat.mkTrCompArray();
                }
            }
            else{
                dat = readSingleFastPhaseLine(name,emStSp,ploidy);
            }
            coll.data.put(name,dat);
            if(st!=null && st.startsWith("certainty")){
                String st1 = br.readLine();
                String st2 = br.readLine();
                double[] cer = new double[st1.length()];
                for(int ij=0; ij<cer.length; ij++){
                    String string = new String(new char[] {st1.charAt(ij), st2.charAt(ij)});
                    if(string.equals("**")) cer[ij] =1.0;
                    else cer[ij] = Double.parseDouble(string)/100.0;
                }
                st = br.readLine();
            }
            if(st!=null && st.startsWith("certainty")){
                String st1 = br.readLine();
                String st2 = br.readLine();
                double[] cer = new double[st1.length()];
                for(int ij=0; ij<cer.length; ij++){
                    String string = new String(new char[] {st1.charAt(ij), st2.charAt(ij)});
                    if(string.equals("**")) cer[ij] =1.0;
                    else cer[ij] = Double.parseDouble(string)/100.0;
                }
                coll.uncertaintyPhase.put(name, cer);
                st = br.readLine();
            }
            i++;
        }
        
        public void read() throws Exception{
            while(!(st = br.readLine()).startsWith("#")){}
            if(st!=null && st.startsWith("# locs")){
                this.readLocation();
            }
            while(st!=null &&  !st.startsWith("#")){
                st = br.readLine();
            }
            while(st!=null && st.startsWith("#")){
                this.readMultiLine();
            }
            coll.length = coll.data.values().iterator().next().length();
        }
    }
    public static PhasedDataState readBasic(BufferedReader br, Class clazz, EmissionStateSpace emStSp, int noSnps) throws Exception{
        PhasedDataState dat =  new PhasedDataState("", noSnps, emStSp, (short)-1);
        String st = "";
        int i=0;
        while((st = br.readLine())!=null){
            
            String[] str = st.split("\\s+");
            if(str.length!=emStSp.haplopairListSize()) {
            	throw new RuntimeException("!!");
            }
            double[] probs = new double[str.length];
            for(int k=0; k<probs.length; k++){
                probs[k] = Double.parseDouble(str[k]);
            }
           /* int max = Constants.getMax(probs);
            if(probs[max]>0.9999){
                dat.emissions[i] = new IntegerDistribution(max);
            }
            else*/
            dat.emissions[i] = new SimpleExtendedDistribution(probs, Double.POSITIVE_INFINITY);
            
            i++;
        }
        return dat;
        
    }

    public static SimpleDataCollection readFastPhaseOutput(BufferedReader br, Class clazz, EmissionStateSpace emStSp) throws Exception{
       int ploidy = ((CompoundEmissionStateSpace)emStSp).getMembers().length;
    	FastphaseFormatReader r = new FastphaseFormatReader(br, clazz, emStSp, ploidy);
       r.read();
       br.close();
       return r.coll;
     }
    
   
    
   
    SortedMap<Integer, Integer> getRecSites(double mult, String pos){
        SortedMap<Integer, Integer> res = new TreeMap<Integer, Integer>();
//        PIGData res = SimpleScorableObject.make(pos, this.length());
        for(int i=1; i<this.length(); i++){
            double d =  1 - Math.exp(-1*(this.loc.get(i) - this.loc.get(i-1))* mult* Constants.probCrossOverBetweenBP);
            if(Constants.rand.nextDouble()<d){
                res.put(i,1);
            }
          
        }
      
      return res;
    }
    
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#printTrioData(java.io.PrintWriter)
     */
    public void printTrioData(PrintWriter pw){
        int i=0;
        for(Iterator<PIGData> it = data.values().iterator(); it.hasNext();i++){
            PIGData dat_i = it.next();
            if(((ComparableArray)dat_i.getElement(0)).isNested()){
                PIGData[] data_i = dat_i.split();
                     int childs = data_i[2].noCopies()==1 ? 1 : 2;
                     pw.println(i+"\t 1 \t 0 \t 0 \t 1 \t urn \t urn:lsid:dcc.hapmap.org:Sample:"+data_i[1].getName()+":1");
                     pw.println(i+"\t 2 \t 0 \t 0 \t 2 \t urn \t urn:lsid:dcc.hapmap.org:Sample:"+data_i[0].getName()+":1");
                     pw.println(i+"\t 3 \t 1 \t 2 \t "+childs+"\t urn \t urn:lsid:dcc.hapmap.org:Sample:"+data_i[2].getName()+":1");
               
            }
    }
    }
  

   /* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#get(java.lang.String)
 */
public PIGData get(String st) {
       PIGData dat = data.get(st);
       //    this.get(name_index.get(st));
   //    if(!dat.getName().equals(st)) throw new RuntimeException("!! " +dat.getName()+" "+st);
       return dat;
    }
    
   
    
 
 
   
     
   public Map<String, double[]> uncertaintyVitPhase = new HashMap<String, double[]>();
   public Map<String, double[]> uncertaintyPhase = new HashMap<String, double[]>();
     public Map<String, SortedMap<Integer, Integer>[]> recSites =  new HashMap<String, SortedMap<Integer, Integer>[]>();
     public Map<String, double[]> noSwitches = new HashMap<String, double[]>();
     public Map<String, PIGData> viterbi = new HashMap<String, PIGData>();
     public Map<String, HaplotypeEmissionState> viterbiL = new HashMap<String, HaplotypeEmissionState>();
   //  public EmissionStateSpace emStSp =
  //     Emiss.stateSpace;
     DataCollection(int length){
         this.length = length;
        //maf = new HaplotypeEmissionState("maf", length, emStSp.size(), emStSp);
     }
     DataCollection(){
     }
     
    /* private static boolean equals(ComparableArray comp, ComparableArray prev){
         for(int j=0; j<comp.size(); j++){
             if(comp.get(j) instanceof ComparableArray){
                 if(!equals((ComparableArray)comp.get(j), (ComparableArray)prev.get(j))) return false;
             }
             else{
                 if(!comp.get(j).equals(prev.get(j))){
             }
         }
     }*/
     
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#countSwitches()
     */
    public int countSwitches(){
        int sum=0;
        for(Iterator<PIGData> keys = this.viterbi.values().iterator(); keys.hasNext();){
            PIGData dat = keys.next();
            ComparableArray prev = (ComparableArray)dat.getElement(0);
            for(int i=1; i<dat.length(); i++){
                ComparableArray comp = (ComparableArray) dat.getElement(i);
                for(int j=0; j<comp.size(); j++){
                    if(!comp.get(j).equals(prev.get(j))){
                        sum++;
                    }
                }
                prev = comp;
            }
        }
        return sum;
    }
     
   /*  public List<Clusters> getClusters(int numF){
         List<Integer> bound = getBlockBoundaries();
         List<Clusters> result = new ArrayList<Clusters>();
         for(int ik=0; ik<bound.size(); ik++){
             Clusters res = new Clusters(numF, bound.get(ik), ik==bound.size()-1 ? length() : bound.get(ik+1));
             int i = bound.get(ik);
             for(Iterator<PIGData> it = this.viterbi.values().iterator(); it.hasNext();){
                 PIGData nxt = it.next();
                 ComparableArray dat= (ComparableArray) nxt.getElement(i);
               Integer[] val = new Integer[dat.size()];
                 for(int j=0; j<dat.size(); j++){
                     val[j] = (Integer)dat.get(j);
                  
                  
                 }
                 res.add(val, nxt.getName());
             }
             res.sort();
             result.add(res);
         }
         return result;
     }*/
     /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getBlockBoundaries()
     */
    public List<Integer> getBlockBoundaries(){
              SortedSet<Integer> res = new TreeSet<Integer>();
              res.add(0);
         for(Iterator<String> keys = this.viterbi.keySet().iterator(); keys.hasNext();){
             getBlockBoundaries(keys.next(), res);
         }
         return new ArrayList<Integer>(res);
     }
     
     /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getBlockBoundaries(java.lang.String, java.util.Set)
     */
    public void getBlockBoundaries(String key, Set<Integer> res){
  
        PIGData dat = this.viterbi.get(key);
        Comparable prev =  dat.getElement(0);
        for(int i=1; i<dat.length(); i++){
            Comparable comp =  dat.getElement(i);
           if(!comp.equals(prev)){
               res.add(i);
           }
            prev = comp;
        }
       //  return res;
     }
     
     /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#restricToAlias(java.util.Collection)
     */
  //  public abstract void restricToAlias(Collection<String> alias);
     
     /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getL(java.lang.String)
     */
     
     
      
   /*  public void restrictSites(int i) {
         if(i<this.length()){
             
             if(loc.size()>0) loc = loc.subList(0,i);
             if(names.size()>0) names = names.subList(0,i);
             if(maf!=null) maf.restrict(i);
             if(this.homo_count.size()>0) homo_count = homo_count.subList(0, i);
         }
             
     }*/
    
     
    
    
     
     public EmissionState getL(String st) {
         
         EmissionState dat = dataL.get(st);
         if(dat==null) return null;
         if(!dat.name.equals(st)) throw new RuntimeException("!!");
         return dat;
      }
     
     public Map<String, boolean[]> getUncertainPositions(){
         Map<String, boolean[] > m = new HashMap<String, boolean[]>();
         for(Iterator<String> it = this.getKeyIterator(); it.hasNext();){
             String key = it.next();
             m.put(key, getUncertainGenotypePositions(key));
         }
         return m;
     }
     /** returns for each position true if the position is uncertain */
     public boolean[] getUncertainGenotypePositions(String key){
         EmissionState emst = this.dataL(key);
         boolean[] res = new boolean[emst.length()];
         for(int i=0; i<emst.length(); i++){
             if(emst.getFixedInteger(i)!=null) {
               res[i] = false;  
             }
            
             else{
                 double[] emiss = emst.getEmiss(i);
                 double sum = Constants.sum(emiss);
                 int max = Constants.getMax(emiss);
                 if(emiss[max] >0.999999){
                     res[i] = false;
                 }
                 else
                 res[i] = true;
             }
          //   res[i] = !res[i];
         }
         return res;
     }
     protected DataCollection(DataCollection dat){
         this.index = dat.index;
         this.pheno =dat.pheno;
        this.name = dat.name;
        this.chrom = dat.chrom;
         this.loc = new ArrayList<Integer>(dat.loc);
         this.length = dat.length;
       //  this.cc = dat.cc;
         this.dataL = new TreeMap<String, EmissionState>();
         for(Iterator<Entry<String, EmissionState>> it = dat.dataL.entrySet().iterator(); it.hasNext();){
             Entry<String, EmissionState> nxt = it.next();
             dataL.put(nxt.getKey(), (EmissionState) nxt.getValue().clone());
         }
       
         this.length = dat.length;
     //    this.maf = dat.maf!=null ? (EmissionState) dat.maf.clone() : null;
      //   this.names = dat.names;
       
         this.ped = dat.ped;
         this.snpid= new ArrayList<String>(dat.snpid);
      //   this.homo_count = dat.homo_count;
         this.data = new TreeMap<String, PIGData>();
         for(Iterator<PIGData> it = dat.data.values().iterator(); it.hasNext();){
             PIGData data_i = (PhasedDataState) it.next().clone();
             data.put(data_i.getName(), data_i);
         }
         this.alleleA =dat.alleleA!=null ?  new ArrayList(dat.alleleA) : null;
         this.alleleB = dat.alleleB!=null  ? new ArrayList(dat.alleleB) : null;
        // this.count_null_sites = dat.count_null_sites;
     }
  
       
   
     
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getFirstIndexAbove(long)
     */
    public int getFirstIndexAbove(long pos){
        for(int i=0; i<this.loc.size(); i++){
            if(this.loc.get(i)>pos) return i;
        }
        return this.loc.size();
    }
     
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getIndex(java.lang.Integer)
     */
    public int getIndex(Integer name){
            return loc.indexOf(name);
           
        }
    
    
    
    String getPrintString(int no, String st){
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<no; i++){
            sb.append(st);
        }
        return sb.toString();
    }
    
   
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#length()
     */
    public final int length() {
        return length;
    }

public Comparator<PIGData> comp = new Comparator<PIGData>(){

    public int compare(PIGData arg0,
            PIGData arg1) {
        String name1 = arg0.getName();
        String name2 = arg1.getName();
      
           
           return name1.compareTo(name2);
        //}
    }
    
};

    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#printHapMapFormat(java.io.File, java.lang.String)
     */
    public void printHapMapFormat(File f, String chr) throws Exception{
        String len = "%"+this.getKeys().iterator().next().length()+"s";
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
        StringBuffer header = new StringBuffer("%8s\t%8s");
        String[] toPrint = new String[2 + data.size()];
        toPrint[0] = "chrom"; toPrint[1] = "pos"; 
        List<PIGData> list = new ArrayList(data.values());
        Collections.sort(list, comp);
        for(int i=2; i<toPrint.length; i++){
            toPrint[i] = list.get(i-2).getName();
            header.append("\t"+len);
        }
        String headerSt = header.toString();
       // int pos_index =s_i==null ? 0 : getFirstIndexAbove(s_i[0]);
        pw.println(String.format(headerSt, toPrint));
       
        for( int pos_index =0; pos_index < loc.size(); pos_index++){
            toPrint[0] = chr; toPrint[1] = this.loc.get(pos_index).toString(); 
            
            for(int i=2; i<toPrint.length; i++ ){
              ComparableArray comp = ((ComparableArray)list.get(i-2).getElement(pos_index));
               toPrint[i] = comp.toStringPrint();
             
              
            }
            pw.println(String.format(headerSt, toPrint));
        }
        pw.close();
    }
  
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#printHapMapFormat(java.io.PrintWriter, java.util.Collection, boolean)
     */
    public void printHapMapFormat(PrintWriter pw, Collection<Integer> toD, boolean style) throws Exception{
       // this.printHapMapFormat(pw, toD, style, this.data, );
        this.printHapMapFormat(pw, this.indiv(),toD, style, new String[] {"snpid", "loc"}, new String[0], new String[] {"data"}, "%7s");
    }
    
    public void printHapMapFormat(PrintWriter pw, 
    		List<String> indiv,
            Collection<Integer> toD, boolean style,
            String[] initTags,
            String[] phenoTags,
            String[] tags, 
            String len, int start, int end) throws Exception{
        Integer[] phenIndex = new Integer[pheno==null ? 0 : pheno.size()];
        for(int i=0; i<phenIndex.length; i++){
            phenIndex[i] = i;
        }
        printHapMapFormat(pw, indiv, toD, style, initTags, phenoTags, tags, len,start, end, new int[] {0,1,2}, phenIndex, true);
        
    }
    
    public void printHapMapFormat(PrintWriter pw, 
    		List<String> indiv,
            Collection<Integer> toD, boolean style,
            String[] initTags,
            String[] phenoTags,
            String[] tags, 
            String len) throws Exception{
        Integer[] phenIndex = new Integer[pheno==null ? 0 : pheno.size()];
        for(int i=0; i<phenIndex.length; i++){
            phenIndex[i] = i;
        }
        printHapMapFormat(pw, indiv, toD, style, initTags, phenoTags, tags, len, loc.get(0), loc.get(this.loc.size()-1)+1, new int[] {0,1,2}, phenIndex, true);
        
    }
  
    /*class for printing compressed file */
    class WriteCompressed {
    	File dir1;
    	PrintWriter indiv;
    	final int start;
    	final int end;
    	WriteCompressed(File dir2, int start, int end, boolean state) throws Exception{
    		
    		
    		File dir4 = new File(dir2,name);
    		dir4.mkdir();
    		File dir3 = new File(dir4,state ? "states" : "geno");
    		dir3.mkdir();
    		
    		File dir = new File(dir3,Constants.chrom0());
    		dir.mkdir();
    		this.state = state;
    		this.start = start;
    		this.end = end;
    		this.dir1 = dir;
    		PrintWriter osw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "SNPS"))));
    		DataCollection dc = DataCollection.this;
    		for(int i=0; i<dc.snpid().size(); i++){
    			  if(loc.get(i)<start ) continue;
                  if(loc.get(i)>=end) break;
                    osw.write("chr"+dc.chrom()+"\t");
                    osw.write(dc.loc().get(i)+"\t");
                    osw.write((dc.loc().get(i)+40)+"\t");
                    osw.write(dc.snpid().get(i)+"\t");
                    if(dc.alleleA()!=null && dc.alleleA().size()>0 && dc.alleleB()!=null && dc.alleleB().size()>0){
                        Character maj = dc.alleleA().get(i);
                        Character min = dc.alleleB().get(i);
                        if(maj==null) maj = 'N';
                        if(min==null) min = 'N';
                        osw.write(maj); osw.write("\t");
                        osw.write(min);
                   
                    }//
                    osw.write("\n");
                }
    		osw.close();
    		osw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Name"))));
    		try{
    		 String headSNP = dc.headSNP();
    	        osw.write(headSNP+"\n");
    	        
    	        osw.write(dc.head_snp()+"\n");
    	        osw.write("id\tdata_index");
    	     
    	       osw.write("\n");
    		osw.close();
    		}catch(Exception exc){
    			exc.printStackTrace();
    		}
    		indiv = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Samples"))));
    	}
    	final boolean state;
    	
    	  public void print(String key) throws Exception{
          	indiv.println(key);
          	  for( int pos_index =0; pos_index < loc.size(); pos_index++){
          		  if(loc.get(pos_index)<start ) continue;
                  if(loc.get(pos_index)>=end) break;
          		PrintWriter osw = new PrintWriter(new FileWriter(new File(dir1, snpid.get(pos_index)), true));
         // 	  PseudoDistribution hes_ =  ((HaplotypeEmissionState)dataL.get(key)).emissions[pos_index];
	              String hes = getCompressedString(key, pos_index,true,state);
	              String hes1 = getCompressedString(key, pos_index,false,state);
	          
	              osw.write( hes1+
	                      (hes==null  ? "" :"\t"+hes));
	              osw.write("\n");
	          		osw.close();
          	  }
    	  }
    	  public void print(List<String> keys) throws Exception{
    		  for(int k=0; k<keys.size(); k++){
    			  indiv.println(keys.get(k));
    		  }
            	
            	  for( int pos_index =0; pos_index < loc.size(); pos_index++){
            		  if(loc.get(pos_index)<start ) continue;
                    if(loc.get(pos_index)>=end) break;
            		PrintWriter osw = new PrintWriter(new FileWriter(new File(dir1, snpid.get(pos_index)), true));
           // 	  PseudoDistribution hes_ =  ((HaplotypeEmissionState)dataL.get(key)).emissions[pos_index];
            		 for(int k=0; k<keys.size(); k++){
  	              String hes = getCompressedString(keys.get(k), pos_index,Constants.noSamples()>1,state);
  	              String hes1 = null;//getCompressedString(keys.get(k), pos_index,false,state);
            		
  	              osw.write( hes1+
  	                      (hes==null  ? "" :"\t"+hes));
  	              osw.write("\n");
            		 }
  	          		osw.close();
            	  }
      	  }
    }
    
    
    /*class for printing averages */
     class PrintHapMapFormat {
    	
        Collection<Integer> toD;boolean style;
        String[] initTags;
      //  String[] phenoTags;
        String[] tags;
       // int[] tags_length;
        //String len;
        int start; int end; 
        //int[] type; 
      StringBuffer[][] sb;
       PrintWriter pw;
        
        public void print(String key) throws Exception{
        	
        	outer:for(int i=0; i<tags.length; i++){
        		int len = sb[i].length;
        		for(int j=0; j<len; j++){
        			this.sb[i][j] = new StringBuffer();
        			if(len==1)sb[i][j].append(key+"_"+tags[i]);
        			else sb[i][j].append(key+"_"+tags[i]+"."+j);
        		}
        		
        	  for( int pos_index =0; pos_index < loc.size(); pos_index++){
        		  if(loc.get(pos_index)<start ) continue;
                  if(loc.get(pos_index)>=end) break;
                  if(toD!=null&& toD.contains(pos_index)) continue;
                 
                  String info = getInfo(tags[i], key, pos_index, style);
                  if(info==null) continue outer;
                  else info = info.trim();
                  if(len==1)sb[i][0].append("\t"+info);
                  else{
                	  String[] info1 = info.trim().split("\\s+");
                	  	if(pos_index==0 && info1.length>sb[i].length){
                		  sb[i] = new StringBuffer[info1.length];
                		  for(int j=0; j<info1.length; j++){
                			  this.sb[i][j] = new StringBuffer();
                  			if(len==1)sb[i][j].append(key+"_"+tags[i]);
                  			else sb[i][j].append(key+"_"+tags[i]+"."+j);
                		  }
                	  	}
                  for(int j=0; j<info1.length; j++){
                	 
                	 String inf =info1[j] ;
                            sb[i][j].append("\t"+inf);
                  }
                  }
               }
        	  for(int j=0; j<len; j++){
                 pw.println(sb[i][j].toString());
        	  }
               
              }
        }
       
        public PrintHapMapFormat(PrintWriter pw,
        	
                Collection<Integer> toD, boolean style,
                String[] initTags,
               
                String[] tags, 
                int start, int end,   boolean printHeader
        ){
        	this.pw = pw;
        	this.toD = toD;
        	this.style = style;
        	this.initTags = initTags;
        	this.tags = tags;
        //	this.tags_length = new int[tags.length];
        	sb = new StringBuffer[tags.length][];
        	for(int i=0; i<tags.length; i++){
        		//tags_length[i] = tags[i].equals("state")? Constants.modify0.length : 1;
        		sb[i] = new StringBuffer[tags[i].equals("state") ? Constants.getFinalModelLength() : 1];
        	}
        	this.start = start; this.end = end;  
        	if(printHeader){
             for(int i=0; i<initTags.length; i++){
                pw.print(initTags[i]);
                for(int k=0; k<loc.size(); k++){
                	  if(loc.get(k)<start ) continue;
                      if(loc.get(k)>=end) break;
                      if(toD!=null&& toD.contains(k)) continue;
                	pw.print("\t"+getInfo(initTags[i], k));
                }
                pw.println();
             }
        	}
           
             //pw.close();
        }
        
    }
    
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#printHapMapFormat(java.io.PrintWriter, java.util.Collection, boolean, java.lang.String[], java.lang.String[], java.lang.String)
     */
    public void printHapMapFormat(PrintWriter pw,
    		List<String> indiv,
            Collection<Integer> toD, boolean style,
            String[] initTags,
            String[] phenoTags,
            String[] tags, 
            String len, int start, int end, int[] type, Integer[] phenIndex, boolean printHeader
//            Map<String, PIGData> data, 
  //          Map<String, Double[]> uncertainty,
    //        Map<String, Double[]> uncertaintyPhase,
      //      Map<String, EmissionState> emState
    ) throws Exception{
       PrintHapMapFormat phf = new PrintHapMapFormat(pw, toD,  printHeader, initTags,   tags, start,  end, printHeader);
       for(int i=0; i<indiv.size(); i++){
    	   phf.print(indiv.get(i));
       }
       phf.pw.close();
    }
	public void updateIndex(int i) {
		
	}
  

public int numClasses(int phenIndex){
      if( pheno.phenotypeDistribution[phenIndex] instanceof SkewNormal) return 0;
      else{
          return ((SimpleExtendedDistribution)this.pheno.phenotypeDistribution[phenIndex]).probs.length;
      }
  }
   int currentPhenScIndex = -1;
   int currentPosScIndex = -1;
   int currentType = -1;
   //public static  boolean scoreRegression = false;
 public synchronized String getPhenInfo(String string , int pos_index, int phenIndex, int type){
     int numCl = numClasses(phenIndex);
    return "";
    
 }
 

  

public String getInfo(String string, int pos_index) {
		if(string.startsWith("chisq")){
		    throw new RuntimeException("!!");
		}
		else if(string.startsWith("armitage")){
		    throw new RuntimeException("!!");
		}
		/*else if(string.equals("maf")){
		double[] res =  this.maf.getEmiss(pos_index);
			Double[]r = new Double[res.length];
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<res.length; i++){
				r[i] = res[i];
				sb.append("%5.3g/");
			}
			
			return String.format(sb.toString(), r);
		}*/
		else if(string.equals("snpid")){
			return snpid==null ? "" : snpid.get(pos_index);
		}
		else if(string.equals("loc")){
			return loc.get(pos_index)+"";
		}
		else if(string.equals("index")){
            return pos_index+"";
        }
		else if(string.startsWith("hwe")){
            Double[]r = new Double[hwe.length];
            for(int i=0; i<r.length; i++){
                r[i] = hwe[i].get(pos_index);
            };
            StringBuffer sb = new StringBuffer();
            for(int i=0; i<r.length; i++){
            //    r[i] = res[i];
                sb.append("%5.3f/");
            }
            return String.format(sb.toString(), r);
		}
		else{
		    try{
		    List l =(List) this.getClass().getField(string).get(this);
		    if(l.size()==0) return null;
		    if(l.get(pos_index)==null) return "null";
		    else return l.get(pos_index).toString();
		    }catch(Exception exc){
		        exc.printStackTrace();
		        return null;
		    }
		}
		//else throw new RuntimeException("!!");
	}
	/*private double getHWE(int pos_index, boolean state, boolean ignoreMissing, int data_index) {
	//    EmissionStateSpace emstsp1 =this.stSp[0];
	    CompoundEmissionStateSpace emstsp;
	    if(this.dataL.size()>0) emstsp = (CompoundEmissionStateSpace) this.dataL.values().iterator().next().getEmissionStateSpace();
	    else emstsp
	    =Emiss.getEmissionStateSpace(1);
      double[] count1 = new double[emstsp.getMembers()[0].size()];
      double[] count = new double[emstsp.size()];  
      int cnt =0;
      outer: for(Iterator<String> it = this.getKeyIterator(); it.hasNext();){
          String key = it.next();
         
          PhasedDataState dat = (PhasedDataState) data.get(key);
          if(dat.noCopies()!=2) continue;
          int d_index = dat.emissions[pos_index].getDataIndex();
          if(d_index ==-1) throw new RuntimeException("!!");
          if(d_index!=data_index) continue;
          
          cnt++;
          if(!state){
                  ComparableArray compa = (ComparableArray) dat.getElement(pos_index);
                  if(ignoreMissing &&  compa.countNull()==2) continue outer;
                  int ind = emstsp.get(compa);
                //  maf.addCount(ind, 1.0, i);
                  int[] indices = emstsp.getMemberIndices(ind);
                      for(int k=0; k<indices.length; k++){
                          count1[indices[k]]+= 1.0;
                      }
                      count[ind]+=1.0;
          }
          else{
                  EmissionState nxt = dataL.get(key);
                  Integer fixed = nxt.getFixedInteger(pos_index);
                  if(fixed!=null){
                      int ind = nxt.getFixedInteger(pos_index);
                      if(ignoreMissing &&  ((ComparableArray)emstsp.get(ind)).countNull()==2) continue outer;
                    //  maf.addCount(ind, 1.0, i);
                      int[] indices = emstsp.getMemberIndices(ind);
                          for(int k=0; k<indices.length; k++){
                              count1[indices[k]]+= 1.0;
                          }
                          count[ind]+=1.0;
                  }
                  else{
                      double[] prob = nxt.getEmiss(pos_index);
                    double sum =  Constants.sum(prob);
                  
                      inner: for(int ind =0; ind<prob.length; ind++){
                          if(ignoreMissing &&  ((ComparableArray)emstsp.get(ind)).countNull()==2) continue inner;
                          int[] indices =((CompoundEmissionStateSpace) nxt.getEmissionStateSpace()).getMemberIndices(ind);
                          for(int k=0; k<indices.length; k++){
                              count1[indices[k]]+= prob[ind];
                          }
                          count[ind]+=prob[ind];
                      }
                  }
                  
          }
         
          
      }
     
      if(cnt==0 || Constants.sum(count1)<0.001) return 1.0;
//      double[] exp_count = new double[count.length];
      SimpleExtendedDistribution.normalise(count1);
      double sum = Constants.sum(count);
      double statistic =0;
      double sum1 = 0;
      int cntNonZero =0;
      for(int i=0; i< count1.length; i++){
          if(count1[i]>0.001) cntNonZero++;
      }
     for(int i=0; i<emstsp.size(); i++){
         int[] memb = emstsp.getMemberIndices(i);
         double exp_count = 
             (count1[memb[0]] * count1[memb[1]])*sum;
         if(memb[0]!=memb[1]) exp_count = 2*exp_count;
         sum1+=exp_count;
         double diff =exp_count - count[i]; 
         if(exp_count<1e-10 && count[i]>0.01) throw new RuntimeException("!! "+exp_count+" "+count[i]+" "+pos_index);
         if(Math.abs(diff)>0.001 && exp_count>0.001){
             statistic += Math.pow(diff,2) / exp_count;
         }
        
     }
     if(Math.abs(sum- sum1) > 0.001) throw new RuntimeException("!! "+sum+" "+sum1);
     return ChiSq.chi2prob(degF(cntNonZero), statistic);
     
      
    }*/
    private int degF(int cntNonZero) {
      return (cntNonZero*(cntNonZero-1))/2;
    }

    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#size()
     */
	//public abstract int size() ;
	 public int size() {
	      return dataL.size();
	    }

    
	 

    
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getDeletedPositions()
     */
    public final List<Aberation> getDeletedPositions(int[] cn, String key, PrintWriter pw, double thresh){
        List<Aberation> obj= new ArrayList<Aberation>();
        EmissionStateSpace emstsp = this.dataL.get(key).getEmissionStateSpace();
    	 List<Integer> geno = new ArrayList<Integer>();
    		for(int i=0; i<cn.length; i++){
    			int[] res = emstsp.getGenoForCopyNo(cn[i]);
    			
    			if(res!=null){
    			for(int k=0; k<res.length; k++){
    				int[] res2 = emstsp.getHaplopairForGeno(res[k]);
    				for(int jj=0; jj<res2.length; jj++){
    				int[] res1 = emstsp.getHaploFromHaploPair(res2[jj]);//emstsp.getHaplopairForGeno(res[k]);
    				for(int j=0; j<res1.length; j++){
    					geno.add(res1[j]);
    				}
    				}
    			}
    			}
    		}
        //for( Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
    		HaplotypeEmissionState dat1 = (HaplotypeEmissionState)dataL.get(key);
            HaplotypeEmissionState dat = Constants.noSamples()>1 ?  dat1: (HaplotypeEmissionState)data.get(key);
         //List l = Arrays.asList(dat.emissions).subList(807, 840);
        // System.err.println(cn+": "+geno);
            obj.addAll(dat.getDeletedPositions(geno, this.loc, pw, cn[0], thresh, Constants.continueThresh()));
       // }
        
        Collections.sort(obj);
        return obj;
    }
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#getDeletedPositions(double, java.util.SortedSet, java.util.SortedSet)
 */
public  void getDeletedPositions(double mafThresh, SortedSet<Integer> ampl, SortedSet<Integer> del){
    int[] count = new int[loc.size()];
    int[] countAmpl = new int[loc.size()];
    double cntThresh = mafThresh * this.data.size()*2;
    System.err.println("count thresh is "+cntThresh+" "+this.data.size()+" "+mafThresh);
    Arrays.fill(count, 0);
      for( Iterator<PIGData> it = this.iterator(); it.hasNext();){
          PIGData dat = it.next();
          for(int i=0; i<dat.length(); i++){
              ComparableArray compa = (ComparableArray)dat.getElement(i);
              for(int k=0; k<compa.size(); k++){
                  Emiss em = (Emiss) compa.get(k);
                  if(em.noCopies()==0) count[i]++;
                  else if(em.noCopies()==2) countAmpl[i]++;
              }
          }
    }
      for(int i=0; i<count.length; i++){
          if(count[i]>cntThresh)del.add(i);
          if(countAmpl[i]>cntThresh) ampl.add(i);
      }
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#printLocations(java.io.PrintWriter)
 */
public void printLocations(PrintWriter pw){
    for(int i=0; i<this.loc.size(); i++){
        pw.println(loc.get(i));
    }
    pw.close();
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#getDeletedPositions()
 */
public List<Aberation> getDeletedPositions(boolean deletion){
    List<Aberation> obj= new ArrayList<Aberation>();
    for( Iterator<PIGData> it = this.iterator(); it.hasNext();){
        PhasedDataState dat = (PhasedDataState)it.next();
        obj.addAll(dat.getDeletedPositions(dataL.get(dat.getName()), dat.noCop(), deletion));
    }
    
    Collections.sort(obj);
    return obj;
}

public void fix(Locreader loc, int thresh){
    List<Integer> fixed = new ArrayList<Integer>();
    List<Integer> nonFixed = new ArrayList<Integer>();
    for(int i=0; i< this.loc.size(); i++){
        int pos = this.loc.get(i).intValue();
        Location overl = loc.contains(pos,thresh);
       if(overl==null){
          fixed.add(this.loc.get(i));
        for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
            HaplotypeEmissionState nxt = ((HaplotypeEmissionState)it.next());
            short di = nxt.emissions[i].getDataIndex();
           EmissionStateSpace emstsp =  nxt.getEmissionStateSpace();
            double[] emiss = nxt.getEmiss(i);
            double[] newd = new double[emiss.length];
            for(int k=0; k<newd.length; k++){
                int cn = emstsp.getCN(k);
                if(cn!=2) newd[k] =0;
                else newd[k] = emiss[k];
            }
            SimpleExtendedDistribution.normalise(newd);
            int max = Constants.getMax(newd);
            if(newd[max]>0.999){
            	
                nxt.emissions[i] = new IntegerDistribution(max, nxt.getEmissionStateSpace());
                nxt.emissions[i].setDataIndex(di);
            }
            else nxt.emissions[i] = new SimpleExtendedDistribution(newd, Double.POSITIVE_INFINITY);
        }
       }
       else{
           nonFixed.add(this.loc.get(i));
       }
    }
  //  this.calculateMaf(true);
    System.err.println("finished fixing \nfixed:"+fixed+"\n not fixed"+nonFixed);
}
public Locreader getMergedDeletions(boolean extend, boolean deletion){
    List<Aberation> l= getDeletedPositions(deletion);
  
    Locreader locr = new Locreader((long)Long.MAX_VALUE-1,"");
    for(int i=0; i<l.size(); i++){
        Aberation ab = l.get(i);
        Integer start = this.loc.get(ab.start);
        Integer finish = this.loc.get(ab.end);
        if(extend){
            if(ab.start>0){
                start = this.loc.get(ab.start-1)+1;
            }
            if(ab.end <loc.size()-1){
               finish = this.loc.get(ab.end+1)-1;
            }
        }
        Location loc = new Location(start, finish, "");
        loc.setName("");
        loc.setNoCop(0);
        locr.add(loc);
    }
    locr.merge(0.00);
    return locr;
}


class PrintAberations{
	private PrintWriter pw;
	private PrintWriter pw1;
	//final int rep;
	  String pst = new String("%-7s  %7s %7s %7s %7s  %7s %7s %7s");
	    String pst1 = new String("%-7s  %7s %7d %7d %7d  %5.3g %7s %5.3g");
	 
	    public PrintAberations(File parentFile) {
	    	File cnvFile = new File(parentFile, "cnv");
	    	File pointwise = new File(parentFile, "cnv1");
	   	   // this.rep =Emiss.getSpaceForNoCopies(Constants.backgroundCount()).copyNumber.size();
	    	cnvFile.mkdir();
	    	pointwise.mkdir();
	    	try{
	    	this.pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(cnvFile, name+".txt"))));
	    	  pw.println(String.format(pst, new Object[] { "Sample", "Chr", "FirstProbe", "LastProbe" ,"NbSnp" ,"length_min", "Type", "Avg_certainty"}));
		   		
	    	this.pw1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(pointwise, name+".txt"))));
	    	}catch(Exception exc){
	    		exc.printStackTrace();
	    	}
			// TODO Auto-generated constructor stub
		}
		public void print(String key){
			double[] scaleLoc = Constants.scaleLoc();
			List<Aberation> obj;
			if(true){
	    	  obj = getDeletedPositions(new int[] {0}, key,pw1, Constants.longThresh());
	    	 HaplotypeEmissionState data_ = ((HaplotypeEmissionState)dataL.get(key));
	    	  Integer rep = data_.noCop();
	    	   for(int i=1; i<=rep*Constants.maxCopies(); i++){
	    		   if(i!=rep){
	    			   obj.addAll(getDeletedPositions(new int[] {i}, key,pw1, Constants.longThresh()));
	    		   }
	    	   }
			}else{
				 HaplotypeEmissionState data_ = ((HaplotypeEmissionState)dataL.get(key));
				  Integer rep = data_.noCop();
				int[] lt = new int[rep];
				for(int k =0; k<rep; k++){
					lt[k] = k;
				}
				int[] gt = new int[rep*Constants.maxCopies()-(rep)];
				for(int k=0; k<gt.length; k++){
					gt[k] = rep+k+1;
				}
				 obj = getDeletedPositions(lt, key,pw1, Constants.longThresh());
		    			   obj.addAll(getDeletedPositions(gt, key,pw1, Constants.longThresh()));
		    		
				
			}
	    	   for(int i=0; i<obj.size(); i++){
	    	        Aberation ab = obj.get(i);
	    	        Object[] obj1 = new Object[8];
	    	      
	    	        obj1[0] = ab.name;
	    	       
	    	        File dir = new File(Constants.outputDir());
	    	        String chrom = Constants.chrom(0);
	    	        int end_index = ab.end;
	    	        int start_index = ab.start;
	    	        double start = loc.get(start_index);
	    	        double end, nextstart;
	    	      
	    	        if(scaleLoc!=null){
	    	        	double a1 = ((double)start/scaleLoc[1]);
	    	            final double chr = Math.floor(a1);
	    	        	start = (a1 -chr)*scaleLoc[0];
	    	        	chrom = ((int)chr)+"";
	    	        	int maxindex = chrToMaxIndex.get((int)chr-1);
	    	        	if(end_index>=maxindex){
	    	        		if(end_index>maxindex){
		    	        		   ab.start = maxindex+1; 
		    	        		   i = i-1;
		    	        		}
	    	        		end_index = maxindex;
	    	        		nextstart =  loc.get(maxindex) + (loc.get(maxindex) - loc.get(maxindex-1));
	    	        	}else{
	    	     
	    	        		nextstart = loc.get(end_index+1);
	    	        	}
	    	        	if(end_index<start_index) throw new RuntimeException("end less than start index");
	    	        	
	    	        	end = loc.get(end_index);
	    	        	
	    	        	end = (((double)end/scaleLoc[1])-chr)*scaleLoc[0];
	    	        	nextstart = (((double)nextstart/scaleLoc[1])-chr)*scaleLoc[0];
	    	         }else{
	    	        	end= loc.get(end_index);
	  	    	        nextstart = (end_index+1 < loc.size()) ? loc.get(end_index+1) : loc.get(end_index) + (loc.get(end_index) - loc.get(end_index-1));
	    	        }
	    	        end = Constants.proportionOfDistanceToNextProbeToExtend()*nextstart  + (1-Constants.proportionOfDistanceToNextProbeToExtend())*end;
		    	     
	    	        obj1[1] = chrom;
	    	        obj1[2] = (int) start;
	    	        obj1[3] =(int) end;
	    	        
	    	        obj1[4] = (end_index-start_index+1);
	    	        obj1[5] = ((Number)obj1[3]).doubleValue() - ((Number)obj1[2]).doubleValue();
	    	        if(((Number)obj1[5]).doubleValue()<0){
	    	        
	    	        	throw new RuntimeException("end less than start");
	    	        }
	    	      //  obj1[5] = "0";
	    	        obj1[6] = Constants.inversion() ? (ab.copy==0 ? "no_recomb" : ab.copy==1 ? "normal" : "inversion") : ab.copy+"";
	    	        obj1[7] = ab.certainty;
	    	       String st = String.format(pst1, obj1);
	    	      // System.err.println(st);
	    	        pw.println(st);
	    	        
	    	    }
	    	   if(DataCollection.this.dc instanceof MatchedDistributionCollection){
					pw.println("##Estimated cellularity:"+((MatchedDistributionCollection)dc).cellularity[0]);
					pw.println("##Estimated ratio:"+((MatchedDistributionCollection)dc).ratio[0]);
				    pw.println("##Version:"+2.0);
				}
	    }
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#printDeletedPositions(java.io.PrintWriter)
 */
public final void printDeletedPositions(File f){
   PrintAberations pa= new PrintAberations(f);
   for(Iterator<String> it =this.dataL.keySet().iterator(); it.hasNext();){
	   pa.print(it.next());
   }
   pa.pw.close();
   pa.pw1.close();
}


/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#writeLocation(java.io.PrintWriter, java.util.Collection)
 */
public void writeLocation(PrintWriter pw, Collection<Integer> toDrop){
   // if(toDrop.size()==0){
     //   throw new RuntimeException("size is 0!");
    //}
    //also print locs
   // if(toDrop.size()==loc.size()) throw new RuntimeException("!!");
    List<Integer> loc1 = new ArrayList<Integer>(loc.size());
    for(int i = 0 ; i<loc.size(); i++){
        if(toDrop==null || !toDrop.contains(i)){
            loc1.add(loc.get(i));
        }
    }
 //int sz = toDrop.size();
    
    SimpleScorableObject.printIdLine("", pw, loc1.size());
    pw.println("# locs");
    Integer maxLoc = loc1.get(loc1.size()-1);
    int numberPos= maxLoc.toString().length();
    char[][] sb = new char[numberPos+1][loc1.size()+numberPos+10];
    char sp = ' ';
    for(int i=0; i<sb.length; i++){
       Arrays.fill(sb[i], sp);
    }
    int j=0;
   
    for(int i=0; i<loc1.size(); i++){
        char[] st = loc1.get(i).toString().toCharArray();
       System.arraycopy(st, 0, sb[j],i,  st.length);
            j++;
            if(j==sb.length) j=0;
    }
    for( j=0; j<sb.length; j++){
       pw.println(new String(sb[j]));
    }
    pw.println("# end locs");
}
public static void printUncertainty(double[] unc, PrintWriter pw, Collection<Integer> toDrop){
    StringBuffer[] sb = new StringBuffer[]{new StringBuffer(), new StringBuffer()};
    sb[0] = new StringBuffer();
    sb[1] = new StringBuffer();
    for(int j=0; j<unc.length; j++){
        if(toDrop==null || ! toDrop.contains(j)){
        double d = Math.round(100*unc[j]);
         if(d>=100){
             sb[0].append("*");sb[1].append("*");
         }
         else{
             char[] st = (d+"").toCharArray();
             sb[0].append(st[0]); sb[1].append(st[1]);
         }
        }
    }
    pw.println(sb[0].toString());
    pw.println(sb[1].toString());
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#writeFastphase(java.util.Map, java.util.Map, java.util.Map, java.io.PrintWriter, boolean, boolean, java.util.Collection)
 */
public void writeFastphase(Map<String, PIGData> data1,
          Map<String, double[]>uncertaintyPhase,
        PrintWriter[] pw,  
    //    boolean restrictToPairs, 
        boolean printUncertainty, boolean mark, Collection<Integer> toDrop) throws Exception{
    writeFastphase(data1,  uncertaintyPhase, pw, printUncertainty, mark, true, toDrop);

}

public class WriteMostLikely{
	PrintWriter[] pw;
	public WriteMostLikely(PrintWriter[] pw){
		for(int i=0; i<pw.length; i++){
		pw[i].println(size());
	        
	        pw[i].println(loc.size());
		}
		this.pw = pw;
		
	}
	public void write(String key, List<Character> a_all, List<Character> b_all){
		  //PIGData dat =  DataCollection.this.viterbi.get(key) ;
		  HaplotypeEmissionState dat = DataCollection.this.viterbiL.get(key);
        if(key.startsWith("!bg")) return;
     List<String> nme = ((SimpleEmissionStateSpace)((CompoundEmissionStateSpace)  ((PhasedDataState) dat).getEmissionStateSpace()).getMembers()[0]).nme;
     PseudoDistribution[] ems = ((PhasedDataState)dat).emissions;
     PseudoDistribution dist = ems[(int)Math.floor((double)ems.length/2.0)];
    
     int ind = DataCollection.this.getIndex(key);
     pw[ind].println(key+"\t"+dist.getProbString(nme));
    		
       System.err.println("T");
     // writeFastphase(key,dat,expand,mark,toDrop,pw[ind], a_all, b_all);
	}
	public void close() {
	for(int i=0; i<pw.length; i++)pw[i].close();
		
	}
}
/** class for writing phased results.*/
public class WritePhased{
	
	PrintWriter[] pw;
	boolean printUncertainty;
	boolean mark,expand;
	Collection<Integer> toDrop;
	boolean states;
	public WritePhased(PrintWriter[] pw, boolean printUncertainty, boolean mark,
			boolean expand, Collection<Integer> toDrop, boolean states){
		for(int i=0; i<pw.length; i++){
		pw[i].println(size());
	        
	        pw[i].println(loc.size());
		}
		this.pw = pw;
		this.printUncertainty = printUncertainty;
		this.mark = mark;
		this.toDrop = toDrop;
		this.states = states;
	}
	public void write(String key, List<Character> a_all, List<Character> b_all){
		  PIGData dat = states ? DataCollection.this.viterbi.get(key) : data.get(key);
		  EmissionState datL = states ? DataCollection.this.viterbiL.get(key) : dataL.get(key);
          if(key.startsWith("!bg")) return;
         int ind = DataCollection.this.getIndex(key);
        writeFastphase(key,dat,(HaplotypeEmissionState)datL,expand,mark,toDrop,pw[ind], a_all, b_all);
	}
	public void close() {
		for(int k=0; k<this.pw.length; k++){
			pw[k].close();
		}
		
	}
}

/* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#writeFastphase(java.util.Map, java.util.Map, java.util.Map, java.io.PrintWriter, boolean, boolean, boolean, java.util.Collection)
     */
    public void writeFastphase(Map<String, PIGData> data1,
              Map<String, double[]>uncertaintyPhase,
            PrintWriter[] pw,  
        //    boolean restrictToPairs, 
            boolean printUncertainty, boolean mark, boolean expand, Collection<Integer> toDrop
            
    ) throws Exception{
        
    	WritePhased wp = new WritePhased(pw, printUncertainty, mark, expand, toDrop, false);
        for(Iterator<String> it =data1.keySet().iterator(); it.hasNext();){
        
            String key = it.next();
          wp.write(key, this.alleleA, this.alleleB);
       
          
        }
        //pw.close();
    }
  
    private void writeFastphase(String key, PIGData dat, HaplotypeEmissionState datL, boolean expand,
			boolean mark, Collection<Integer> toDrop, PrintWriter pw,
			List<Character> alleleA, List<Character> alleleB
			) {
    	
    	//int ind = getIndex(key);
    	  if(dat==null){
          	System.err.println("warning is null "+key);
          	return;
          }
           dat.print(pw, expand, mark , toDrop, alleleA, alleleB, datL.emissions) ;
          double[] unc = null;//dat.getUncertainty(dataL.get(dat.getName()));
      
         if(unc!=null){
             double[] uncP = uncertaintyPhase.get(dat.getName());
             pw.println("certainty");
             printUncertainty(unc, pw, toDrop);
             pw.println("certainty phase");
            if(uncP!=null) printUncertainty(uncP, pw, toDrop);
         }
		
	}
    protected int getIndex(String key){
    	return 0;
    }
	public Iterator<String> getKeyIterator() {
       /* if(cc!=null && cc.size()>0){
            final Set<String> cases = new HashSet<String>();
            final Set<String> controls = new HashSet<String>();
            final Set<String> other = new HashSet<String>();
            for(Iterator<String> it = this.getKeys().iterator(); it.hasNext();){
                String key = it.next();
                Boolean b = cc.get(key);
                if(b==null) other.add(key);
                else if(b) cases.add(key);
                else controls.add(key);
            }
            return new Iterator<String>(){
                Iterator<String>[] it = new Iterator[] {other.iterator(), cases.iterator(), controls.iterator()};
                int ind =0;
                {
                    while(!it[ind].hasNext()){
                        ind++;
                    }
                }
                public boolean hasNext() {
                    return ind < it.length && it[ind].hasNext();
                }
                public String next() {
                   String nxt = it[ind].next();
                   if(!it[ind].hasNext()){
                       ind++;
                   }
                   return nxt;
                }
                public void remove() {
                    // TODO Auto-generated method stub
                    
                }
            };
        }
        else */return this.getKeys().iterator();
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getNames()
     */
    public List<String> getNames(){
       return new ArrayList(this.getKeys());
    }

    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#trim(int)
     */
    public void trim(int i) {
    	 List<String> toKeep = new ArrayList<String>();
    	 List<String> names = this.getNames();
    	 toKeep.addAll( Arrays.asList(Constants.toInclude1));
        if(names.size()>i ){
             toKeep.addAll(names.subList(0, Math.max(1, i-toKeep.size())));
    //   if(bg!=null && !toKeep.contains(this.bg.name)) toKeep.add(bg.name);
    //    		this.restricToAlias(toKeep);
    //    }
       
        	this.restricToAlias(toKeep);
        	
        }
    }

   
   /* public int restrict(int st1, int end1) {
        int il = this.firstGreaterThan(this.loc, st1)-1;
        int ir = this.firstGreaterThan(loc, end1);
        
      
        this.loc = this.loc.subList(il, ir);
        this.snpid = this.snpid.subList(il, ir);
        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
            it.next().restrictSites(il, ir);
        }
        this.length = loc.size();
        return il;
    }*/
   

    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#getLocations()
     */
    public List<Integer> getLocations() {
        // TODO Auto-generated method stub
        return loc;
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#replace(java.util.Map)
     */
    public void replace(Map<String, PIGData> newData){
        this.data = newData;
       // this.makeDataIndex();
    }

    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#arrangeDataAccordingToPedigree()
    
    public void arrangeDataAccordingToPedigree(){
        System.err.println("arranging according to pedigree");
        data = this.arrangeDataAccordingToPedigree(this.data);
        viterbi = this.arrangeDataAccordingToPedigree(this.viterbi);
        if(this.dataL.size()>0){
        Map<String, String> mother = this.ped.mother;
        Map<String, String> father = this.ped.father;
        Map<String, EmissionState> toAdd= new HashMap<String, EmissionState>();
          for(Iterator<String> it = mother.keySet().iterator(); it.hasNext();){
              String chi = it.next();
         
              EmissionState moth =getL( mother.get(chi));
              EmissionState child = getL( chi);
              String fathKey = father.get(chi);
              EmissionState trio;
              if(fathKey!=null){
                  EmissionState fath = getL(fathKey);
                  trio =  EmissionState.getEmissionState( moth,fath, child);
              }
              else{
                  trio =   EmissionState.getEmissionState( moth, child);
              }
                  toAdd.put(trio.getName(), trio);
            
          }
         replaceL(toAdd);
        }
    } */
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#replaceL(java.util.Map)
     */
    public void replaceL(Map<String, EmissionState> newData){
        this.dataL = newData;
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.collection.DataC#arrangeDataAccordingToPedigree(java.util.Map)
    
    public final Map<String, PIGData> arrangeDataAccordingToPedigree(Map<String, PIGData> datai) {
       Map<String, PIGData> toAdd= new HashMap<String, PIGData>();
        for(Iterator<String> it = ped.mother.keySet().iterator(); it.hasNext();){
            String chi = it.next();
            PhasedDataState fath = (PhasedDataState)datai.get(ped.father.get(chi));
                PhasedDataState moth =(PhasedDataState) datai.get( ped.mother.get(chi));
                PhasedDataState child = (PhasedDataState) datai.get( chi);
                if( moth==null || child==null) continue;
               // String[] chi_nam = child.getName().split("\\.");
         //       if(!chi_nam[0].startsWith(moth.getName())) throw new RuntimeException("!!");
           //     if(!chi_nam[1].startsWith(fath.getName())) throw new RuntimeException("!!");
                if(fath==null){
                    PIGData moth1 = SimpleScorableObject.make(new PhasedDataState[] { moth, child}, false,";", this.stSp[1]);
                    toAdd.put(moth1.getName(), moth1);
                }
                else{
                    PIGData moth1 = SimpleScorableObject.make(new PhasedDataState[] { moth,fath, child}, false, ";", this.stSp[1]);
                    toAdd.put(moth1.getName(), moth1);
                }
          
        }
        return toAdd;
      // replace(toAdd);
    //  this.makeDataIndex();
    } */

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#iterator()
 */
public Iterator<PIGData> iterator() {
  return data.values().iterator();
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#set(int, lc1.dp.data.representation.PIGData)
 */
public void set(int ij, PIGData newDat) {
    this.data.put(newDat.getName(), newDat);
    
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#setPedigree(lc1.ensj.PedigreeDataCollection)
 */
public void setPedigree(PedigreeDataCollection pedData) {
   this.ped = pedData;
    
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#getKeys()
 */
public Set<String> getKeys() {
	
    if(dataL.size()!=0)
   return this.dataL.keySet();
    else return this.data.keySet();
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#checkConsistency()
 
public final double[][] checkConsistency(){
   double[] cons1 = this.checkConsistency(this.data);
   double[] cons2 =  this.checkConsistency(this.viterbi);
   return new double[][] {cons1, cons2};
}*/

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#checkConsistency(java.util.Map)
 
public final double[] checkConsistency(Map<String, PIGData> datai) {
    Map<String, PIGData> toAdd= new HashMap<String, PIGData>();
    int incons =0;
    int cons =0;
    int total =0;
     for(Iterator<String> it = ped.mother.keySet().iterator(); it.hasNext();){
         String chi = it.next();
             PIGData fath = datai.get(ped.father.get(chi));
             PIGData moth =datai.get( ped.mother.get(chi));
             PIGData child = datai.get( chi);
             if(fath==null || moth==null || child==null) continue;
             total+=fath.length();
            outer: for(int i=0; i<fath.length(); i++){
                ComparableArray c = (ComparableArray) child.getElement(i);
                ComparableArray m = (ComparableArray) moth.getElement(i);
                ComparableArray f = (ComparableArray) fath.getElement(i);
                if(c.size()==1){
                    if(m.contains(c.get(0)))  cons++;
                    else incons++;
                }
                else if(c.size()==2){
                    if(m.contains(c.get(0)) && f.contains(c.get(1))) cons++;
                    else if(m.contains(c.get(1)) && f.contains(c.get(0))) cons++;
                    else incons++;
                }
                else throw new RuntimeException("!!");
            }
     }
     if(incons+cons!=total) throw new RuntimeException("!!");
     return new double[] {incons,cons};
  
 }*/
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#containsKey(java.lang.String)
 */
public boolean containsKey(String name2) {
   return dataL.containsKey(name2);
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#summarise()
 
public void summarise() {
    this.split();
   double[] ld_av =  this.calcLDAverage();
   System.err.println("ld_av is "+ld_av[0]+" "+ld_av[1]);
 int no = this.loc.size();
 int length = this.loc.get(no-1)-loc.get(0);
 List<Integer> l = new ArrayList<Integer>();
 double desn1 =0;
 for(int i=1; i<loc.size(); i++){
     l.add(loc.get(i)-loc.get(i-1));
     desn1 +=loc.get(i)-loc.get(i-1);
 }
 Collections.sort(l);
 int mid =(int)Math.round((double) l.size()/2.0);
 double desn = l.get(mid);

 //  double desn = (double) length/(double)no;
   System.err.println("length is "+length);
   System.err.println("density is "+desn+" "+((double)desn1/ (double)loc.size()));
   System.err.println("num is "+no);
   System.exit(0);
    
}
*/
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#calcLDAverage()
 */
public double[] calcLDAverage(){
    double[][]res = calcLD();
    double[] sum = new double[] {0,0};
    for(int i=0; i<res[0].length; i++){
        sum[0]+=res[0][i];
        sum[1]+=res[1][i];
    }
    sum[0] = sum[0] /(double) (res[0].length);
    sum[1] = sum[1] /(double) (res[0].length);
    return sum;
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#calcLD()
 */
public double[][] calcLD(){
    return calcLD(Emiss.a().toString(), Emiss.b().toString());
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#calcLD(java.lang.String, java.lang.String)
 */
public double[][] calcLD(String st1, String st2){
    LDCalculator lower = new R2Calculator();
    LDCalculator upper = new DPrimeCalculator();
    double[][] res = new double[2][this.length()-1];
    for(int i=1; i<this.length(); i++){
        int j = i-1;
        String [] st_i = new String[] {st1, st2};
        String [] st_j = new String[] {st1, st2};
        double[][] d = getLD(i,j,st_i, st_j);
        res[0][i-1] = upper.calculate(d);
        res[1][i-1] = lower.calculate(d);
        if(Double.isNaN(res[0][i-1])) throw new RuntimeException("!! "+i+" "+d[0][0]+" "+
                d[0][1]+" "+d[1][0]+" "+d[1][1]);
    }
    return res;
    
}
/*public void split(){
    split(1.0, false);
}*/
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#split()
 
public void split(double prob, boolean repackage){
    Set<PhasedDataState >set =new HashSet<PhasedDataState>();
    for(Iterator<PIGData> it = this.iterator(); it.hasNext();){
        PhasedDataState da = (PhasedDataState) it.next();
        if(Constants.rand.nextDouble()<=prob){
         
            PhasedDataState[] dat = da.split();
            for(int j=0; j<dat.length; j++){
                dat[j].setName(da.getName()+"_"+j);
                PhasedDataState toadd = repackage ? 
                        (PhasedDataState)  SimpleScorableObject.make(new PhasedDataState[]{dat[j]}, false, "", this.stSp[0]):
                            dat[j];
                        Double[] phenV = da.phenValue();
                 toadd.setPhenotype(phenV);
                set.add(
                        toadd);
            }
        }
        else{
            set.add(da);
        }
    }
    this.data.clear();
    for(Iterator <PhasedDataState> it = set.iterator(); it.hasNext();){
        PhasedDataState da = it.next();
        this.data.put(da.getName(), da);
    }
    this.indiv = indiv();
}*/
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#calcLD(calc.LDCalculator, calc.LDCalculator, lc1.dp.data.representation.PIGData, double[][])
 */
public void calcLD(
        LDCalculator upper, 
        LDCalculator lower,
        PIGData poss,
        double[][]res
        ){
    for(int i=0; i<this.length(); i++){
        res[i][i]=1.0;
        for(int j=i+1; j<this.length(); j++){
                ComparableArray c_i = (ComparableArray)poss.getElement(i);
                ComparableArray c_j = (ComparableArray)poss.getElement(j);
                String [] st_i = new String[] {c_i.get(0).toString(), c_i.get(1).toString()};
                String [] st_j = new String[] {c_j.get(0).toString(), c_j.get(1).toString()};
                double[][] d = getLD(i,j,st_i, st_j);
                res[i][j] = lower.calculate(d);
                res[j][i] = upper.calculate(d);
        }
    }
}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#calcLD(calc.LDCalculator, lc1.dp.data.representation.PIGData, boolean, int, int)
 */
public List<Double>  calcLD(
        LDCalculator upper, 
        PIGData poss,
         boolean left,
        int pos, int lenthresh
        ){
    List<Double> res = new ArrayList<Double>();
  //  res[pos] = 1.0;
    res.add(0.0);
    for(int ik=1; ; ik++){
       int i = left ? pos-ik : pos+ik; 
       if(i<0 || i >= loc.size() ||
               Math.abs(this.loc.get(i) - this.loc.get(pos))>lenthresh) break;
        ComparableArray c_i = (ComparableArray)poss.getElement(i);
        ComparableArray c_j = (ComparableArray)poss.getElement(pos);
                String [] st_i = new String[] {c_i.get(0).toString(), c_i.get(1).toString()};
                String [] st_j = new String[] {c_j.get(0).toString(), c_j.get(1).toString()};
                double[][] d = getLD(i,pos,st_i, st_j);
                try{
                res.add( upper.calculate(d));
                }catch(Exception exc){
                    System.err.println(d[0][0]+" "+d[0][1]+" "+d[1][0]+" "+d[1][1]);
                    exc.printStackTrace();
                }
    }
    return res;
}

/*
 * st1 are the different possible haplotypes at pos1
 * st2 are the different possible haplotypes at pos2
 */
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#getLD(int, int, java.lang.String[], java.lang.String[])
 */
public double[][] getLD(int pos1, int pos2, String[] st1, String[] st2){
    double[][] res = new double[][] {{0.0,0.0},{0.0,0.0}};
    int[] start = new int[] {pos1, pos2};
    List<Emiss>[] hapl = this.getHaplotypes(new int[] {pos1, pos2});
   // String[][] potential = new String[st1.length][st2.length];
    for(int i=0; i<st1.length; i++){
        for(int j=0; j<st2.length; j++){
       
       //     String[] str = new String[] {st1[i], st2[j]};
          ///** res[i][j] = countHaplotypes(hapl, st1[i]+st2[j]);
         // Logger.global.info("count "+res[i][j]);
        }
    }
    if(allzero(res[0]) && allzero(res[1])) throw new RuntimeException("!!");
    return res;
}

private boolean allzero(double[] ds) {
   for(int i=0; i<ds.length; i++){
       if(ds[i]!=0) return false;
   }
   return true;
}

private double countHaplotypes(List<String> hapl, String str) {
    int cont =0;
    outer: for(int i=0; i<hapl.size(); i++){
    //	System.err.println("comp "+hapl[i])
           if(hapl.get(i).equals(str)) cont++;
    }
    return cont;
}

public List<Emiss>[] getHaplotypes(int[] pos){
	List<Emiss>[] res = new List[pos.length];
	for(int i=0; i<res.length; i++){
		res[i] = getHaplotypes(pos[i]);
	}
	
	return res;
}



public List<Emiss> getHaplotypes(int pos){
   List<Emiss> res = new ArrayList<Emiss>();
 //  Set<String> keys = data.keySet();
    Iterator<PIGData> it = this.data.values().iterator();
   
    outer: for(int i=0; it.hasNext(); i++){
        PIGData nxt = it.next();
       // nxt.getStringRep(start, end);
        ComparableArray compa = (ComparableArray) nxt.getElement(pos);
        ;
        for(int k=0; k<compa.size(); k++){
        	res.add(((Emiss)compa.get(k)));
        }
        
    }
    return res;
}
public List<Comparable> getHaplotypes1(int pos){
	   List<Comparable> res = new ArrayList<Comparable>();
	 //  Set<String> keys = data.keySet();
	    Iterator<PIGData> it = this.data.values().iterator();
	   
	    outer: for(int i=0; it.hasNext(); i++){
	        PIGData nxt = it.next();
	        if(nxt.getName().equals("!bg")) continue;
	        double[] certainty = this.uncertaintyPhase(nxt.getName());
	        double[] probs = ((SimpleExtendedDistribution)this.dataL.get(nxt.getName()).emissions(pos)).probs;
	      //  if()
	    //    HaplotypeEmissionState hes = (HaplotypeEmissionState) this.dataL.get(nxt.getName());
	        ComparableArray compa = (ComparableArray) nxt.getElement(pos);
	        double maxp = probs[Constants.getMax(probs)];
	        if(certainty[pos]>Constants.imputedThresh(0) && maxp > Constants.imputedThresh(0)){
		       // nxt.getStringRep(start, end);
		       
		        for(int k=0; k<compa.size(); k++){
		    //    	if(compa.get(k).toString().equals("B")){
		      //  		System.err.println("B");
		        //	}
		        	res.add(ComparableArray.make(((AbstractEmiss)compa.get(k))));
		        }
	        }
	        else{
	        	for(int k=0; k<compa.size(); k++){
	        		res.add(null);
	        	}
	        }
	        
	        
	    }
	    return res;
	}

public List<Comparable> getValues(int pos, int type, boolean avg){
	 List<Comparable> res = new ArrayList<Comparable>();
	 //  Set<String> keys = data.keySet();
	    Iterator<EmissionState> it = this.dataL.values().iterator();
	   
	    outer: for(int i=0; it.hasNext(); i++){
	       HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
	       EmissionStateSpace em = Emiss.getSpaceForNoCopies(nxt.noCop());
	       // nxt.getStringRep(start, end);
	       if(nxt.emissions[pos]==null){
	    	   System.err.println(this.name);
	    	   res.add(Double.NaN);
	       }
	       else{
	       Integer fixed = nxt.emissions[pos].fixedInteger();
	      
	       if(fixed!=null){
	    	   double v = getVal(fixed.intValue(),em,type);// type==0 ? em.getCN(fixed.intValue()) :
	    		 //em.getBCount(fixed.intValue(), type) ;
	    			  
	          /* if(cn_restrict!=null && 
	        		   em.getCN(fixed.intValue())!=cn_restrict)
	        		  {
	        	   v = Double.NaN;
	           }*/
	        
	        res.add(v);
	       }else if(nxt.emissions[pos] instanceof SimpleExtendedDistribution){
	    	   double[] v = PairEmissionState.pool.getObj(em.genoListSize());
	    	   Arrays.fill(v, 0);
	    	   double[] probs = ((SimpleExtendedDistribution)nxt.emissions[pos]).probs;
	    	   getAvg(probs,em,  type,v);
	    	   int max =Constants.getMax(v);
	    	   if(v[max]<Constants.imputedThresh(0)){
	    		   res.add(Double.NaN);
	    	   }
	    	   else{
	    	   res.add(
	    			   
	    			   avg ? getAvg(v)
	    					   
	    					   : max
	    	   );
	    	   }
	    	   PairEmissionState.pool.returnObj(v);
	       }
	       }
	        
	    }
	    return res;
}

private double getAvg(double[] v) {
	
	double res = 0;
	for(int i=0; i<v.length; i++){
		res+=v[i]*(double)i ;
	}
	return res;
}
public List<Comparable> getGenotypes(int pos){
	   List<Comparable> res = new ArrayList<Comparable>();
	 //  Set<String> keys = data.keySet();
	    Iterator<EmissionState> it = this.dataL.values().iterator();
	   
	    outer: for(int i=0; it.hasNext(); i++){
	       HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
	       EmissionStateSpace em = Emiss.getSpaceForNoCopies(nxt.noCop());
	       // nxt.getStringRep(start, end);
	       Integer fixed = nxt.emissions[pos].fixedInteger();
	      
	       if(fixed!=null){
	       ComparableArray compa = (ComparableArray) em.get(fixed.intValue())
	       
	        ;
	        res.add(compa);
	       }else if(nxt.emissions[pos] instanceof SimpleExtendedDistribution){
	    	   res.add(em.get(Constants.getMax(((SimpleExtendedDistribution)nxt.emissions[pos]).probs)));
	       }
	        
	    }
	    return res;
	}

public List<Comparable> getCN(int pos){
	   List<Comparable> res = new ArrayList<Comparable>();
	 //  Set<String> keys = data.keySet();
	    Iterator<EmissionState> it = this.dataL.values().iterator();
	   
	    outer: for(int i=0; it.hasNext(); i++){
	       HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
	       EmissionStateSpace em = Emiss.getSpaceForNoCopies(nxt.noCop());
	       // nxt.getStringRep(start, end);
	       Integer fixed = nxt.emissions[pos].fixedInteger();
	      
	       if(fixed!=null){
	       ComparableArray compa = (ComparableArray) em.get(fixed.intValue())
	       
	        ;
	        res.add(compa);
	       }else if(nxt.emissions[pos] instanceof SimpleExtendedDistribution){
	    	   double[] pr = ((SimpleExtendedDistribution)nxt.emissions[pos]).probs;
	    	   double sum=0;
	    	   for(int k=0; k<pr.length; k++){
	    		   sum+=pr[k]*em.getCN(k);
	    	   }
	    	   res.add(em.get(sum));
	       }
	        
	    }
	    return res;
	}

Double NaN =null;
public List<Comparable> getIntensity(int pos, boolean r){
	   List<Comparable> res = new ArrayList<Comparable>();
	 //  Set<String> keys = data.keySet();
	    Iterator<EmissionState> it = this.dataL.values().iterator();
	   
	    outer: for(int i=0; it.hasNext(); i++){
	        PseudoDistribution nxt = ((HaplotypeEmissionState) it.next()).emissions[pos];
	        res.add(nxt.getIntensity(r));
	  
	    }
	    return res;
	}

/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#countHaplotypes(int[], int[], java.lang.String[])
 */
public int countHaplotypes(int[]st,int[] en, String[] str){
 int count=0;
    outer: for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
        PIGData nxt = it.next();
        for(int i=0; i<st.length; i++){
            int start = st[i];
            int end = en[i];
            String st1 = nxt.getStringRep(start, end);
          //System.err.println("cf "+st1+" "+str[i]);
            if(!st1.equals(str[i])){
                continue outer;
            }
        }
        count++;
        
    }
    return count;
}
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#printIndiv(java.io.PrintWriter)
 */
public void printIndiv(PrintWriter pw_indiv) {
	for(Iterator<String> it = this.getKeys().iterator(); it.hasNext();){
		pw_indiv.println(it.next());
	}
}

public  void writeFastphase(File pw_hap2, boolean states) throws Exception{
	//pw_hap2.mkdir();
	if(this.loc.size()==0) return;
	PrintWriter[] pw = 
		getPrintWriter(pw_hap2);
		   if(states){
        
    }
//this.writeFastphase(this.,uncertainty,uncertaintyPhase, pw,  false, true, false,calculateMaf1().getConstantPos());
    else{
        this.writeFastphase(data,uncertaintyPhase, pw,  false, true, false,null);
    }
    for(int i=0; i<pw.length; i++){
    	pw[i].close();
    }
   
}  
public PrintWriter[] getPrintWriter(File file) throws Exception{
	file.mkdir();
	return new PrintWriter[] {
			new PrintWriter(new BufferedWriter(new FileWriter(new File(file, 
					(Constants.writeMergedAverageFiles() ? Constants.experiment() : this.name)+
					".txt"))))
	};
}
   

public  double[] get(Map<String, double[]> uncertainty1, String key){
    double[] res =  uncertainty1.get(key);
    if(res ==null){
        uncertainty1.put(key, res = new double[this.length]);
    }
    return res;
}


public  final double[] uncertaintyPhase(String key){
    return get(uncertaintyPhase, key);
    
}

public  double[] uncertaintyVitPhase(String key){
    return get(uncertaintyVitPhase, key);
}

public   void setData(String key, PIGData dat, HaplotypeEmissionState datL, PIGData datvit, HaplotypeEmissionState datvitL, double[] certainty){
    dat.setName(key);
    datL.setName(key);
    this.data.put(key, dat);
    this.uncertaintyPhase.put(key, certainty); //NB storing certainty in genotype here, not phase
   // System.err.println(dataL.size());at
   // EmissionState st1 = this.dataL.get(key);
   // st1.
if(!Constants.useOriginalLikelihoodAsUncertainty())  this.dataL.put(key, datL);
  //  System.err.println(dataL.size());
    if(datvit!=null){
        datvit.setName(key);
        datvitL.setName(key);
        this.viterbi.put(key, datvit);
        this.viterbiL.put(key, datvitL);
    }
    try{
  if(true){
	  this.writeResults(key, this.alleleA, this.alleleB);
  }
    }catch(Exception exc){
    	exc.printStackTrace();
    }
  if(!Constants.savePhasedConfiguration()){
	 // if(key.equals("95111")){
	//	  System.err.println("here");
	 // }
	 this.removeKey(key, datvit!=null);
  }
  //this.writeFastphase(key);
  
   
}
public void removeKey(String key,boolean states){
	  if(states){
	       
	        this.viterbi.put(key,null );
	        this.viterbiL.put(key,null);
	    }
}
PrintAberations printAb = null;

private WriteMostLikely wp_mostlikely;
public void initialisePrinting(File parentFile) throws Exception{
	//if(Constants.prin)
//	StringBuffer suff = new StringBuffer();
	//String suffix = suff.toString();
	File phasedFile = new File(parentFile, "phased2");
	File phasedFileStates = new File(parentFile, "phasedStates");
	File mostLikelyStates = new File(parentFile, "mostLikelyState");
	File avgFile = new File(parentFile, "avg");
	File compressedFile  = new File(parentFile, "res");
	File samplesFile = new File(parentFile,"sample");
	samplesFile.mkdir();
	this.samplespw = new PrintWriter(new FileWriter(new File(samplesFile,this.name+".txt")));
	compressedFile.mkdir();
	phasedFile.mkdir();
	phasedFileStates.mkdir();
	avgFile.mkdir();
	
	int[] core = Constants.core();
	PrintWriter[] pw = 
		this.getPrintWriter(phasedFile);
		//new PrintWriter(new BufferedWriter(new FileWriter(new File(phasedFile, this.name+".txt"))));
	PrintWriter[] pw_states = 
		this.getPrintWriter(phasedFileStates);
//		new PrintWriter(new BufferedWriter(new FileWriter(new File(phasedFileStates, this.name+".txt"))));
	PrintWriter[] pw_statesMostlikely = 
		this.getPrintWriter(mostLikelyStates);
	
	PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(avgFile, 
			(Constants.writeMergedAverageFiles() ? Constants.experiment() : this.name)+".txt"))));
			
	if(Constants.writePhased()){
		boolean mark = false;
	this.wp = new WritePhased(pw, true, mark, false, null,false);
	}
	if(Constants.saveStates()){
		this.wp_states = new WritePhased(pw_states, false, true, false, null,true);
		if(Constants.mostLikely()) this.wp_mostlikely = new WriteMostLikely(pw_statesMostlikely);

	}
	if(Constants.printAb()){
	 this.printAb = new PrintAberations(parentFile);
	}
	
	String[] strings = Constants.writeAverages;
	
		int start = core[0];
		int end = core[1];
		if(Constants.reverse()){
			end = -core[0];
			start = -core[1];
		}
		if(strings!=null && strings.length>0 && !strings[0].equals("null")){
	this.phm = new PrintHapMapFormat(pw1,null, true,  new String[] {"snpid", "loc"}, 
            strings,   start, end,  true);
		}
		
	if(Constants.writeRes()!=null && Constants.writeRes()){
		this.writeCompressGeno = new WriteCompressed(compressedFile,start,end, false);;
		if(Constants.saveStates()){
			this.writeCompressStates = new WriteCompressed(compressedFile,start,end, true);;
		}
	}
}
PrintHapMapFormat phm;
public WriteCompressed writeCompressGeno, writeCompressStates;
WritePhased wp, wp_states;
PrintWriter samplespw;;
public void writeResults(String key, List<Character> allA, List<Character> allB) throws Exception{
	if(dataL.containsKey(key)){
		this.samplespw.println(key);
	if(phm!=null){
		phm.print(key);
	}
	if(!Constants.bufferCompress()){
	if(this.writeCompressGeno!=null){
		writeCompressGeno.print(key);
	}
	if(this.writeCompressStates!=null){
		writeCompressStates.print(key);
	}
	}
	if(this.wp!=null ){
		wp.write(key,allA,allB);
	}
	if(this.wp_states!=null){
		wp_states.write(key,null,null);
	}
	if(this.wp_mostlikely!=null){
		wp_mostlikely.write(key,null,null);
	}
	if(this.printAb!=null){
		printAb.print(key);
	}
	
	}
}
public void finishedPrinting(){
	if(phm!=null) phm.pw.close();
	if(wp!=null) wp.close();
	if(this.wp_states!=null)wp_states.close();
	if(this.wp_mostlikely!=null)this.wp_mostlikely.close();
	if(printAb!=null){
		printAb.pw.close();
		printAb.pw1.close();
	}
	if(this.writeCompressGeno!=null){
		writeCompressGeno.indiv.close();
		CompressDir.compress(writeCompressGeno.dir1);
	}
	if(this.writeCompressStates!=null){
		writeCompressStates.indiv.close();
		CompressDir.compress(writeCompressStates.dir1);
	}
	if(samplespw!=null) samplespw.close();
}



public  void setViterbi(String key, PIGData datvit){
    this.viterbi.put(key, datvit);
}

public  void setRecSites(String name,
        SortedMap<Integer, Integer>[] sampleRecSites){
    this.recSites.put(name, sampleRecSites);
}

public  void clearViterbi(){
    this.viterbi.clear();
}

public  Iterator<EmissionState> dataLvalues(){
    return this.dataL.values().iterator();
}

public  void clearData(){
    this.data.clear();
}

public  SortedMap<Integer, Integer>[] recSites(String j){
    return this.recSites.get(j);
}




/*public  EmissionState maf(){
    return this.maf;
}*/

public  Object ped(){
    return this.ped;
}

public EmissionState dataL(String string) {
    EmissionState res = this.dataL.get(string);
    return res;
}

public void makeDistributions(int index){}

public Object headerObject(String st,Object prev){
    String str = st.toLowerCase();
    if(str.startsWith("index") || str.startsWith("sample")) return null;
    else if(str.startsWith("name") || str.startsWith("snp")) return this.snpid;
    else if(str.startsWith("pos") || str.startsWith("loc")) return this.loc;
    else if(str.startsWith("chr")) return null;
    else if(str.startsWith("address")) return null;
    else{
        try{
                return this.getClass().getField(st).get(this);
               // if(f!=null) return null;
        }catch(Exception exc){}
    }
    return null;
//        return SimpleScorableObject.make(st, this.loc.size());
}



public boolean process(HaplotypeEmissionState data,HaplotypeEmissionState dataL, String header, String geno, int i, double[] missing){
	/*if(header.toLowerCase().indexOf("geno")>=0){
        //if(data.emissions[i]==null){
        int ind = trans(geno);
             data.emissions[i] = new IntegerDistribution(ind);
             data.emissions[i].setDataIndex(this.index);
        //}
             return true;
    }*/
	return false;
    
}

public   Boolean   process(String indiv,  String[] header,  String[] geno, int i, int ploidy,  double[] missing){
    try{
      //  boolean doneGeno = false;
        PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
        HaplotypeEmissionState dataL = (HaplotypeEmissionState) this.dataL.get(indiv);
        for(int k=0; k<header.length ; k++){
        	process(data,dataL,  header[k], k<geno.length  ? geno[k] : null, i,missing);
        	
        }
       
      if(Constants.likelihoodInput()){
    	  ((SimpleExtendedDistribution)dataL.emissions[i]).normalise();
      }
      return dataL.emissions[i].probeOnly();
    }catch(Exception exc){
        exc.printStackTrace();
    }
    return null;
   
}


public static void fixDupl(List<String> indiv){
	Set<String> m = new HashSet<String>();
	for(int i=0; i<indiv.size(); i++){
    	String ind = indiv.get(i);
    	if(m.contains(ind)){
    		//if(Constants.duplicates(i).equals("separate")){
    			String ind1 = ind+Constants.duplSep+"1";
    			for(int k=0; m.contains(ind1); k++){
    				ind1 = ind+Constants.duplSep+k;
    			}
    			indiv.set(i, ind1);
    			
    		//}
    	}
    	else{
        m.add(ind);
    	}// data.put(indiv.get(i), SimpleScorableObject.make(indiv.get(i),loc.size(), stSp[nocopies-1], this.index));
    }
}
public void createDataStructure(List<String> indiv, List<Integer> ploidy, List<Integer> dToIncl){
    for(int i=0; i<dToIncl.size(); i++){
    	String ind = indiv.get(dToIncl.get(i));
    	
        dataL.put(ind, null);
       // data.put(indiv.get(i), SimpleScorableObject.make(indiv.get(i),loc.size(), stSp[nocopies-1], this.index));
    }
}
public EmissionStateSpace[] getNoCopies(){
    
    SortedMap<Integer, EmissionStateSpace> ss = new TreeMap<Integer, EmissionStateSpace>();
    
    for(Iterator<PIGData> it = this.iterator(); it.hasNext();){
        PIGData data = it.next();
        int noCop = data.noCopies();
        if(!ss.containsKey(noCop-1)){
            String name = data.getName();
            EmissionState stat = this.getL(name);
            EmissionStateSpace emstsp =    stat==null ?
                Emiss.getEmissionStateSpace(noCop-1) : stat.getEmissionStateSpace();
            
            ss.put(noCop-1, emstsp);
        }
    }
    EmissionStateSpace[] res  = new EmissionStateSpace[ss.lastKey()+1];
    for(int i=0 ;i <res.length; i++){
        res[i] = ss.get(i);
    }
    return res;
}
public void readSNPInfo(InputStream is, String chr) throws Exception{
    BufferedReader br = new BufferedReader(new InputStreamReader(is));
    String st = "";
    
    while((st = br.readLine())!=null){
        String[] str = st.split("\t");
        if(str[1].equals(chr)){
            int id = snpid.indexOf(str[0]);
            if(id>=0){
                this.alleleA.set(id,str[3].charAt(0));
                this.alleleB.set(id, str[4].charAt(0));
            }
        }
        System.err.println(st);
    }
}
public String[]  restrictToIds(Set<Integer> name2) {
	
	List<String> toDrop = new ArrayList<String>();
	List<String> toKeep = new ArrayList<String>();
	if(name2!=null) {
	outer: for(Iterator<String> it = this.getKeys().iterator(); it.hasNext();){
		String key = it.next();
		HaplotypeEmissionState hes = (HaplotypeEmissionState) ((DataCollection)this).data.get(key);
		Set<Short> s = new HashSet<Short>();
		hes.dataIndices(s);
		for(Iterator<Integer> it1= name2.iterator(); it1.hasNext();){
			if(s.contains(it1.next())){
				toKeep.add(key);
				continue outer;
			}
		}
		System.err.println("reporting "+hes.getName());
		toDrop.add(hes.getName());
	}
	}
	return toDrop.toArray(new String[0]);
	
}
public static int firstGreaterThanOrEqual(List<Integer>loc, int readPos){
    
    int read =0;
    for(read=0; read <loc.size(); read++){
        if(loc.get(read)>readPos) break;
    }
    return read;
}

String[] header;
int lrr_index=0;
List<Integer> chrToMaxIndex = null;  //if using 'all' this keeps track of maximum co-ord for each chrom

List<String> headsLowerCase;
List<String> header_snp;
List<String> header_sample;
//int[] numLevels;
public ProbabilityDistribution[] numLevels(){
    if(pheno==null) return new ProbabilityDistribution[0];
    return pheno.phenotypeDistribution;
}
static private final String REAL_NUMBER = "^[-+]?\\d+(\\.\\d+)?$";
static private final String NONINTEGER = "\\D";


public Phenotypes pheno;
public Set<String> readPhenotypes(BufferedReader br,List<String> header , File incl, List<String> indiv) throws Exception{
	   if(incl.exists() && incl.length()>0 ){
	       pheno = new Phenotypes(incl);
	       for(int i=0; i<header.size(); i++){
	           header.set(i, header.get(i).trim());
	       }
	  //   Logger.global.info("phenotypes are "+pheno.phen);
	       Set<String> todo = new HashSet<String>(indiv);
	        
	           String[][] res = new String[indiv.size()][pheno.phen.size()];
	        
	         int[] alias = new int[pheno.size()];
	         for(int i=0; i<pheno.size(); i++){
	             alias[i] = header.indexOf(pheno.phen.get(i));
	          /*  if(alias[i]<0){
	                 return null;
	             }*/
	         }
	         String str = "";
	         while(todo.size()>0 && (str = br.readLine())!=null){
	            
	            String[] st =  str.split("\\s+");
	        //    System.err.println(st.length);
	          //  if(st.length!=header.size()){
	            //    throw new RuntimeException("!!");
	            //}
	            String specialSplit = Constants.specialCode(this.index);
	            String nme = specialSplit==null ? st[0] : st[0].split(specialSplit)[0];
	            if(todo.contains(nme)){
	                int index = indiv.indexOf(nme);
	               
	              
	                for(int j=0; j<alias.length; j++){
	                	if(alias[j]>=0){
	                    res[index][j] = st[alias[j]];
	                    if(res[index][j].equals("NA") ||res[index][j].equals("null")) res[index][j]=null;
	                    
	                	}
	                		
	                }
	              
	                todo.remove(nme);
	            }
	         }
	         Logger.global.info("did not find pheno for "+((double)todo.size()/(double)indiv.size())+"\n"+todo);
	         int[]   numLevels = new int[alias.length];
	         int[] type = pheno.type(); //0 is real, 1 is intenger, 2 is text
	         Double[][] res1 = new Double[indiv.size()][pheno.size()];
	         double[] min = new double[alias.length];
	         double[] max = new double[alias.length];
	         Arrays.fill(max, Double.NEGATIVE_INFINITY);
	         Arrays.fill(min, Double.POSITIVE_INFINITY);
	       
	         for(int i=0; i<alias.length; i++){
	           /*  int nonNull = 0;
	             try{
	            	 if(alias[i] <0) continue;
	                 System.err.println("doing "+alias[i]+" "+res[0][i]);
	              while(res[nonNull][i] == null || res[nonNull][i].equals("NA") || res[nonNull][i].length()==0){
	                 System.err.println("is null for "+nonNull+" "+i+" "+res[nonNull][i]);
	                 nonNull++;
	             }
	            
	             }catch(Exception exc){
	                 exc.printStackTrace();
	             }*/
	             if(type[i]==0 || type[i]==1){
	                 if(type[i]==0){
	                     numLevels[i] = 0;
	                     for(int k=0; k<res.length; k++){
	                         if(res[k][i]!=null &&  !res[k][i].equals("NA") && res[k][i].length()!=0){
	                             res1[k][i]= Double.parseDouble(res[k][i]);
	                             if(res1[k][i] < min[i]) min[i] = res1[k][i];
	                             if(res1[k][i] < max[i]) max[i] = res1[k][i];
	                         }
	                     }
	                 }
	                 else{
	                     pheno.phenVals[i] = new HashMap<String, Integer>();
	                     for(int k=0; k<res.length; k++){
	                         if(res[k][i]!=null &&  !res[k][i].equals("NA") && res[k][i].length()!=0){
	                           /*  Integer val;
	                             if(res[k][i].equals("1.0")) res1[k][i] = 1;
	                             else if(res[k][i].equals("0.0")) val =0;
	                             else if(res[k][i].equals("NaN")) val = null;
	                             else val = Integer.parseInt(res[k][i]);*/
	                             res1[k][i] = Double.parseDouble(res[k][i]);
	                             if(res1[k][i]!=null && !Double.isNaN(res1[k][i])){
	                            	 pheno.phenVals[i].put(res1[k][i]+"", (int) res1[k][i].doubleValue());
	                             }
	                         }
	                     }
	                     numLevels[i] = pheno.phenVals[i].size();
	                     if(numLevels[i]>10) numLevels[i] = 0;
	                     pheno.phenVals[i] = null;
	                 }
	               
	             }
	             else{
	                 if(pheno.phenVals[i]!=null){
	                     for(int k=0; k<res.length; k++){
	                         if(res[k][i]!=null &&  !res[k][i].equals("NA") && res[k][i].length()!=0){
	                             Integer lev = pheno.phenVals[i].get(res[k][i]);
	                             res1[k][i] =lev==null ? null :  lev.doubleValue();
	                             numLevels[i] = pheno.phenVals[i].size();
	                         }
	                     }
	                 }
	                 else{
	                     pheno.phenVals[i] = new HashMap<String, Integer>();
	                     int curr_lev =0;
	                     for(int k=0; k<res.length; k++){
	                         if(res[k][i]!=null &&  !res[k][i].equals("NA") && res[k][i].length()!=0){
	                             Integer lev = pheno.phenVals[i].get(res[k][i]);
	                             if(lev==null){
	                                 pheno.phenVals[i].put(res[k][i], lev = curr_lev );
	                                 curr_lev++;
	                             }
	                             res1[k][i] = lev.doubleValue();
	                          
	                             numLevels[i] = pheno.phenVals[i].size();
	                         }
	                     }
	                 }
	                 numLevels[i] = pheno.phenVals[i].size();
	                 if(numLevels[i]>10) numLevels[i] = 0;
	             }
	         }
	         pheno.phenotypeDistribution =    new ProbabilityDistribution[pheno.size()];
	         for(int k=0; k<pheno.phen.size(); k++){
	        	//if(numLevels[k]==0) throw new RuntimeException("is continuous trait");
	             pheno. phenotypeDistribution[k] = numLevels[k] >0 ? 
	                 (ProbabilityDistribution)   new SimpleExtendedDistribution1(numLevels[k]) :
	                     (ProbabilityDistribution) new TrainableNormal("", 0.0, 1.0, 1000, 1.0);
	                             //0.0,1.0,0.001, min[k], max[k], 1000.0, 10.0);
	                 //phenotypeDistribution[k].setName(p.get(k));
	         }
	         for(int i=0; i<alias.length; i++){
	        	 if(alias[i]<0) continue;
	             for(int k=0; k<res1.length; k++){
	                 Double v = res1[k][i];
	                 if(v!=null){
	                 pheno.phenotypeDistribution[i].addCount(v,1.0);
	                 }
	             }
	             pheno.phenotypeDistribution[i].transfer(0.0);
	         }
	         for(int k=0; k<indiv.size(); k++){
	             ((HaplotypeEmissionState)this.dataL.get(indiv.get(k))).setPhenotype(res1[k]);
	         }
	         return todo;
	   }
	   return new HashSet<String>();
	  
	}

public void printPheno(PrintWriter pw) throws Exception{
    if(pheno!=null){
        String formatString = this.pheno.print(pw);//.replaceAll("7s", "5.3g");
        String[] toPrint = new String[pheno.phen.size()];
        for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
            HaplotypeEmissionState nxt = (HaplotypeEmissionState) it.next();
            pw.print(nxt.getName()+"\t");
            Double[] phenV = nxt.phenValue();
            if(phenV==null){
                Arrays.fill(toPrint, "null");
            }
            for(int k=0; k<phenV.length; k++){
                toPrint[k] = phenV[k]==null ? "null" : phenV[k].toString();
            }
            pw.println(String.format(formatString,toPrint));
        }
    }
    pw.close();
}
public void writeCompressed(File dir, boolean writeDistr){
	  if(this.loc.size()==0) return;
    dir.mkdir();
    int[] core = Constants.core();
    File out = new File(dir, name.split("_")[0]);
    		
   
    try{
        CompressDC cc = new CompressDC( out, this);
        cc.writeDistribution = writeDistr;
        List<List<String>> l = new ArrayList<List<String>>();
        l.add(this.indiv());
        cc.run(l, core[0], core[1]);
    }catch(Exception exc){
        exc.printStackTrace();
    }
}
protected List<String>indiv;

protected List avgDepth;  //second one is for matched normal
//is called by constructor in case subclasses have to initialise
public void initialise(){
}

public Boolean[] probeOnly;
public  boolean hasABGenos(ZipFile zf){
	return true;
}
boolean abGenos = false;

protected int geno_id, baf_ind, baf_ind1;
protected boolean baf_geno;



public File dir;

public ParseTree tree = null;

/*
 * index is the data index
 * no_copies is the number of chromosomes being modelled
 * mid = [start, end]
 * bf is the build file
 * kb is the distance left and right
 * snp_ids is an optional set of snps to restrict to
 */
DataProjection dp;
public  DataCollection (File f, short index, int no_copies, final int[][] mid,File bf, Collection<String> snp_ids_to_restrict) throws Exception{
//	System.err.println(Constants.print(mid[0]));
	 this.chrom = f.getName().split("\\.")[0];//.split("_")[0];
	
	if(f.exists()){
		
		System.err.println("specified file exists: opening "+f.getName());
	}
	else{
		System.err.println("specified file does not exist "+f.getAbsolutePath());
		String nme = f.getName();
		int len = nme.length();
		
		
		int ind = getFirstNonNumeric(nme);
		
		if(nme.charAt(ind)=='p' || nme.charAt(ind)=='q'){
		
			f = new File(f.getParentFile(), nme.substring(0,ind)+".zip");//nme.substring(ind+1));
			 this.chrom = f.getName().split("\\.")[0].split("_")[0];
		}
		else{
			System.err.println("looking for p/q files");
			File f1 = getKaryoFile(f,mid[0][0]);
			File f2 = getKaryoFile(f,mid[mid.length-1][1]);
			if(!f1.equals(f2)) {
				throw new RuntimeException("going across centromere at ");
			}
			f = f1;
			if(f1.exists()) this.chrom = f1.getName().split("\\.")[0].split("_")[0];
		}
		
		if(!f.exists() || f.length()==0){
			
			File[] fs = f.getParentFile().listFiles(new FileFilter(){
		
				@Override
				public boolean accept(File arg0) {
					if( arg0.getName().startsWith(chrom)){
						final	String[] strs = arg0.getName().split("\\.")[0].split("_");
						int i1 = Integer.parseInt(strs[1])-20000;
						int i2 = Integer.parseInt(strs[2])+20000;
						int m1 = mid[0][0];
						int m2 = mid[0][1];
						int o1 = m2 - i1;
						int o2 = i2 - m1;
						double overl = Math.min(o1,o2);
					
					return overl>0;
					}else return false;
				}
				
			});
			
			if(fs.length==0) f = new File(f.getParentFile(),"all.zip");
			else f = fs[0];
			
		}
	}
	if(!f.exists() || f.length()==0){
		throw new RuntimeException("could not find "+f);
		//f = findFile(f,mid);
	}
	System.err.println("opening "+f+" "+mid[0][0]+"-"+mid[0][1]+" index "+index);
    this.index = index;
   
  this.dir = f.getParentFile();
   File plate = new File(f.getParentFile(), "plate.txt");
    //this.no_copies = no_copies;
   String name1 = Constants.inputDir[index].replaceAll("/", "__"); 
	   //f.getParentFile().getName();
   if(bf!=null && Constants.build.length>1){
	
	  int ind = name1.indexOf('_');
	  name = ind>=0 ? name1.substring(0,ind)+"-"+bf.getName()+name1.substring(ind) :
		  name1+"-"+bf.getName();
	 
	  
   }
   else{
	   name = name1;
   }
   name = name+Constants.suffix(this.index);
   
   
    //chrom = Constants.chrom0();//.substring(2);//chrom.substring(0,chrom.indexOf(".zip")).split("_")[0];
    //List<Integer >loc1; 
    String[] head;  int read=0; int end;
     indiv = new ArrayList<String>();
    
     avgDepth = new ArrayList();
 
    List<Integer> ploidy = new ArrayList<Integer>();

    ZipFile zf = new ZipFile(f);
    String ent_name = ((ZipEntry)zf.getEntries().nextElement()).getName();
   String f_name = f.getName().split("\\.")[0];
   String pref = "";// ent_name.startsWith(f_name) ? f_name+"/" : "";
   //String pref =  ent_name.startsWith(f_name) ? f_name+"/" : "";
   pref = "";
//   pref = "";
        this.abGenos = hasABGenos(zf);
      //  loc1 = new ArrayList<Integer>();
        BufferedReader name_bf = ApacheCompressor.getBufferedReader(zf, pref+"Name");
        if(name_bf==null){
        	name_bf = new BufferedReader(new FileReader(new File(f.getParentFile(),pref+"Name")));
        }
         List<String> headers=  ApacheCompressor.getIndiv(name_bf,  null);
         header = headers.get(0).split("\t");
          headsLowerCase = Arrays.asList(header);
         for(int j=0; j<headsLowerCase.size(); j++){
        	 headsLowerCase.set(j, headsLowerCase.get(j).toLowerCase());
         }
         header_snp = Arrays.asList(headers.get(1).split("\\s+"));
         header_sample = Arrays.asList(headers.get(2).split("\t"));
         for(int i=0; i<header_snp.size(); i++){
        	 header_snp.set(i, header_snp.get(i).trim());
         }
         this.lrr_index = headsLowerCase.indexOf("pc0");
         if(this.lrr_index <0)  this.lrr_index=headsLowerCase.indexOf("log r ratio");
         if(this.lrr_index <0)  this.lrr_index=headsLowerCase.indexOf("depth");
      geno_id =-1;
      int b_ind = -1;
         inner: for(int k=0; k<this.header.length; k++){
 			if(geno_id<0 && header[k].toLowerCase().indexOf("geno")>=0 ||
 					header[k].toLowerCase().indexOf("plus_allele")>=0
 					){
 				geno_id = k;
 			//	break inner;
 			}
 			if(b_ind<0 && header[k].toLowerCase().indexOf("b allele freq")>=0){
 				b_ind = k;
 			}
 		}
      int sample_id = header_sample.indexOf("id");
      int depth_id = header_sample.indexOf("depth");
      int median_id = header_sample.indexOf("median");
      if(sample_id<0) sample_id = header_sample.indexOf("sampleID");
      if(sample_id<0) sample_id = header_sample.indexOf("sample");
      if(sample_id<0) sample_id = header_sample.indexOf("Sample");
      File samplesFile = new File(f.getParentFile(), "Samples");
         if(samplesFile.exists()){
         indiv = ApacheCompressor.getIndiv(samplesFile, sample_id);
         }  else{
        	 indiv = ApacheCompressor.getIndiv(zf, pref+"Samples", sample_id);
         }
         if(indiv.size()==0){
        	 File f1 = new File(f.getParentFile(), "all.zip");
        	 if(f1.exists()){
        		 ZipFile zf1 = new ZipFile(f1);
        		  indiv = ApacheCompressor.getIndiv(zf1, "Samples", sample_id);
        		 zf1.close();
        	 }
         }
         File dpfile = new File(f.getParentFile(),"pcs_in.zip");
         if( Constants.numPcs()>=0 &&  dpfile.exists()){
        	 int column = depth_id < 0 ? median_id : depth_id;
         dp = new DataProjection(dpfile, Constants.numPcs(), indiv, 
        			ApacheCompressor.getIndiv(zf, pref+"Samples",  column)
        		, depth_id >= 0);
         }
         File tree = new File(f.getParentFile(),"tree.txt");
         if(tree.exists()){
        	this.tree = new ParseTree(tree, indiv); 
         }
         List<Integer> dToInc = getSamps();
     
        
         
         int plo_id = header_sample.indexOf("ploidy");
         if(plo_id>0){
        	 List<String> pl = ApacheCompressor.getIndiv(zf,pref+"Samples", plo_id);
        	 for(int k_=0; k_<dToInc.size(); k_++){
        		 int k = dToInc.get(k_);
        		 ploidy.add(Integer.parseInt(pl.get(k)));
        	 }
         }
         
         else{
        	 for(int i=0; i<indiv.size(); i++){
        		 ploidy.add(no_copies);
        	 }
         }
         //code to deal with average depth for each subject which
         //is read from the samples file (column named avgdepth)
        int avgDepthCol = header_sample.indexOf("avgdepth");
       
        if(avgDepthCol<0){
        	avgDepthCol = header_sample.indexOf("info_1");
        }
       
      
        {
        	List<String> ad ;
        	//File samplesFile = new File(f.getParentFile(), "Samples");
        	if(samplesFile.exists()){
        		ad  = ApacheCompressor.getIndiv(samplesFile);//,pref+"Samples");
        	}else{
        		ad = ApacheCompressor.getIndiv(zf,pref+"Samples");
        	}
        	 for(int k_=0; k_<dToInc.size(); k_++){
        		 int k = dToInc.get(k_);
	       		String[] str = ad.get(k).split("\\s+");//split(":");
	       		int avgDepthCol1 = avgDepthCol;
	       		if(avgDepthCol1>=0){
	       			avgDepth.add(Double.parseDouble(str[avgDepthCol1]));
	       		}
	       		else if(avgDepthCol1<0 && str.length>header_sample.size()){ avgDepthCol1 = header_sample.size();
	       		 	avgDepth.add(Double.parseDouble(str[avgDepthCol1]));
	       		}else{
	       			avgDepth.add(ploidy.get(k).doubleValue());
	       		}
	       		
	    	 }
        }
       
       
        /* this is to deal with repeats
         * for(int i=0; i<indiv.size(); i++){
        	 String ik =  indiv.get(i).split("#")[0];
        	 int ind = indiv.indexOf(ik);
        	 if(ind>=0 && ind!=i){
        		 ik = ik+"_"+i;
        	 }
             indiv.set(i,ik);
         }
         */
         Map<String, String> platem = new HashMap<String, String>();
         Map<String, String> alternative = new HashMap<String, String>();
         BufferedReader plate_br = null;
         List plate_file_header = null;
         int plat_head=-1;
         int pat_head= header_sample.indexOf("id");
         
         if(this.header_sample.indexOf("Plate")>=0 || this.header_sample.indexOf("PLATE")>=0){
        	 plate_br = ApacheCompressor.getBufferedReader(zf, "Samples");
        	 plat_head = Math.max(header_sample.indexOf("Plate"),header_sample.indexOf("PLATE"));
        	 plate_file_header = header_sample;
        	 String st = "";
        //int alt_head = 
      		//   Constants.convertIds() ? 
      		 // plate_file_header.indexOf("ALTERNATIVE") : -1;
      
      
      		
    
      	   while((st = plate_br.readLine())!=null){
      		   String[] str = st.split("\t");
      		  
      		//   int ind1 = indiv.indexOf(str[pat_head]);
      		 //  if(ind1<0) ind1 = indiv.indexOf(str[pat_head].split("_")[0]);
      		//   if(ind1>=0){
      			   
      			//   indiv.set(ind1, str[pat_head]);
      		  if(plat_head>=0) platem.put(str[pat_head], str[plat_head]);
      		
      		   
      	   }
      	   plate_br.close();
      	   //following code is for decoding ids, not needed for normal operation
      	  /* if(Constants.convertIds() ){
      		 //  PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("decoding.txt"))));
      		  List<String> res = Compressor.readZipFrom(zf, "Samples");
      		  for(int i=0; i<res.size(); i++){
      			  String str = res.get(i);
      			 // int ind = str.indexOf('\t');
      			//pw.println(alternative.get(str.substring(0,ind))+str.substring(ind));
      		  }
      		  pw.close();
      	   }*/
         }
         if(plate.exists() && plate.length()>0){
        	 BufferedReader plbr = new BufferedReader(new FileReader(plate));
        	 String st = plbr.readLine();
        	 List<String> header = Arrays.asList(st.split("\t"));
        	 int alt_head = header.indexOf("ALTERNATIVE");
        	 int gend_head = header.indexOf("Gender");
        	 pat_head = header.indexOf("PATIENT");
        	 if(pat_head <0){
        		 pat_head = header.indexOf("id");
        	 }
        	 while((st = plbr.readLine())!=null){
        		 String[]str = st.split("\t");
        		  if(alt_head>=0  && Constants.convertIds()){
           			   alternative.put(str[pat_head], str[alt_head]);
           		   }
           		   if((Constants.haploidMale) && gend_head>=0){
           			   String gend = str[gend_head];
           			   
           			   if(gend.equals("M")){
           				   int index2 = indiv.indexOf(str[pat_head]);
           				   if(index2>=0){
           				   ploidy.set(index2, 1);
           				   }
           			   }
           			   
           		 //  }
           		   }
           		   
           		   
        	 }
         }
        if(platem.size()==0){
      	   for(int i=0; i<indiv.size(); i++){
      		   platem.put(indiv.get(i), "");
      	   }
         }
     
      List<String> chr = new ArrayList<String>();
      // if(bf!=null) {
        	  
        	   BufferedReader  buildF = bf==null ? null : getBuildReader(f.getParentFile(), bf,zf);
        	   int[] allel = new int[] {indexOf(header_snp,new String[] {"A","AlleleA","REF", "alleleA"}),
        			   indexOf(header_snp, new String[] {"B", "AlleleB","ALT","alleleB"})};
        	   int strand_id = Constants.strand(index)==null ? indexOf(header_snp,new String[] {"strand","anc_eq_ref"}) : -1;
        	  
        	  if(strand_id>=0 && header_snp.get(strand_id).equals("anc_eq_ref"))  {
        		  strand_represents_ref = true;
        	  }
        	   if(buildF==null){
   	            System.err.println("WARNING - USING INTERNAL BUILD FILE "+(bf==null ? "" : bf.getAbsolutePath() +" does not exist"));
   	         	buildF = ApacheCompressor.getBufferedReader(zf, pref+"SNPS");  
   	         
   	        }
         int bin_id  = header_snp.indexOf("bins");
         int snpid_ = header_snp.indexOf("id");
         if(snpid_<0) snpid_ = header_snp.indexOf("snpid");
         if(snpid_ <0) throw new RuntimeException("no index of id in "+header_snp);
            String chr_ = chrom;//Constants.chrom0().split("_")[0];
            if(chr_.indexOf("chr")<0) chr_ = "chr"+chr_;
             readBuildFile(zf, pref,buildF,chr_,mid,  
                     this.loc,   chr, this.snpid,this.alleleA, this.alleleB,this.strand,
                     1, header_snp.indexOf("chr"), 
                     snpid_,strand_id,bin_id,
                           allel , snp_ids_to_restrict);
             
             String str = Constants.strand(index);
             Boolean str_ = str==null ? null : str.equals("+");
           if(strand.size()==0){
        	   for(int i=0; i<loc.size(); i++){
        		   strand.add(str_);
        	   }
           }
           int[] alias = checkOrder();
           if(alias!=null){
        	   snpid =reorder(snpid,alias); loc=reorder(loc,alias);
        	   alleleA = reorder(alleleA,alias); alleleB = reorder(alleleB,alias);
        	   
           }
           if(Constants.reverse()){
            	 Collections.reverse(snpid);
            	 Collections.reverse(this.alleleA);
            	 Collections.reverse(this.alleleB);
            	 Collections.reverse(this.loc);
            	 for(int i=0; i<loc.size(); i++){
            		 loc.set(i, -loc.get(i));
            	 }
           }
             buildF.close();
             
    
        this.length = loc.size();
   
      this.initialise();
    
    
    String stri = "";
    this.length = loc.size();
    this.fixDupl(indiv);
    this.createDataStructure(indiv, ploidy,dToInc);

    Logger.global.info("reading pheno");
    File parent = f.getParentFile();
   /* File inclFile =  new File(parent, "include.txt");
    while(!exists(inclFile)  && parent!=null){
        parent = parent.getParentFile();
        if(parent!=null){
        	inclFile = new File(parent, "include.txt");
        }
    }*/
    BufferedReader br = DataCollection.getBufferedReader(new File(f.getParentFile(), pref+"Samples.txt"));
    List<String> headerS = this.header_sample;
    if(br==null){
    	br = DataCollection.getBufferedReader(new File(f.getParentFile().getParentFile(), pref+"Samples.txt"));
        if(br==null){
        	ZipArchiveEntry entry  = zf.getEntry(pref+"Samples");
        	if(entry!=null)
        	br = new BufferedReader(new InputStreamReader(zf.getInputStream(entry)));
        	
        }
        else{
        	 headerS =  Arrays.asList(br.readLine().split("\t"));
        }
    }
    else{
       headerS =  Arrays.asList(br.readLine().split("\t"));
  }
    //this.readPhenotypes(br, headerS, inclFile,indiv);
   
  //  EmissionState  es = dataL.values().iterator().next();
    Logger.global.info("processing snps... ");
   
    probeOnly = new Boolean[snpid.size()];
   
    baf_ind = this.indexOf1(headsLowerCase, "gtype:geno".split(":"));
    baf_ind1 = this.indexOf1(headsLowerCase, "b allele:b_allele".split(":"));
    if(baf_ind>=0 && headsLowerCase.get(baf_ind).equals("genotype")) baf_geno = false;
    else baf_geno = baf_ind>=0;
    if(!baf_geno) {
    	baf_ind = baf_ind1;
    }
    
   /* if(Constants.CHECK){
		   checkStatesForNull();
	  }*/
    double[] missing = new double[4]; //first two are missing/non missing; next two are allele freqs
    Boolean cnvP = Constants.cnvP(this.index);
   Set<Integer> todrop = new TreeSet<Integer>();
    double missThresh = Constants.NAthresh();
    double[] lrrVals = new double[dToInc.size()];
    List<Integer> failed = new ArrayList<Integer>();
    int cumR = Constants.cumulativeR(index);
    int max = (int) Math.floor((float)this.length/(float)cumR)*cumR;
    if(max==0) throw new RuntimeException("for this size region set cumR to one");
    for(int i=0; i<max; i++){
    	
    	 if(geno_id>=0 && alleleA.size()<=i && ! this.abGenos) {
    		 this.addAlleles(i, ApacheCompressor.getIndiv(zf, snpid.get(i), geno_id),
    				 b_ind<0 ? null : ApacheCompressor.getIndiv(zf, snpid.get(i), b_ind)
    		 );
    	 }
    	
    	try{
    		String sno = this.snpid.get(i);
    		if(sno==null) throw new RuntimeException("is null");
    		if(Constants.removeIndels() && i< this.alleleA.size() && (alleleA.get(i).equals('I') || alleleA.get(i).equals('D'))){
    			System.err.println("excluded "+this.snpid().get(i)+" "+alleleA.get(i));
    			failed.add(i);
    			baf.add(Double.NaN);
    		}else{
    		try{
    		probeOnly[i] = process(pref+sno,(int) Math.floor((float)i/(float)cumR),zf, ploidy,dToInc,missing, lrrVals);
    		}catch(Exception exc){
    			failed.add(i);
    			baf.add(Double.NaN);
    			System.err.println(exc.getMessage());
			exc.printStackTrace();
    		}
    		}
    		if(Constants.standardiseLRR()){
    			double mean = Constants.calcMean(lrrVals);
    			double sd = Constants.calcSd(lrrVals, mean);
    			
    			for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
    				((HaplotypeEmissionState)it.next()).emissions[i].standardise(mean, sd);
    			}
    		}
           if(missing[0]/(double)dToInc.size()>missThresh){
        	   todrop.add(i);
           }
           if(cnvP!=null && probeOnly[i]!=null){
        	   if(cnvP != probeOnly[i] ){
            	   todrop.add(i);
               }
           }
           
    		if( Constants.fillLikelihood(this.index)>=1.0) probeOnly[i] = true;
    		
    		
    		
    	}catch(Exception exc){
    		System.err.println("could not process "+snpid.get(i));
    		exc.printStackTrace();
    	}
        
     }
    
    if(cumR>1){
    	//need to collapse snpid, loc, etc
    	int rem = this.snpid.size() % cumR;
       //Constants.minB[this.index]=0;
    	snpid = thin(this.snpid,cumR); loc = thin(this.loc,cumR);
    	probeOnly = (Boolean[]) thin(Arrays.asList(this.probeOnly),cumR).toArray(new Boolean[snpid.size()]); 
    	alleleA = thin(this.alleleA,cumR); alleleB = thin(this.alleleB,cumR); strand = thin(this.strand,cumR);
    	baf = thin(this.baf,cumR); 
    	this.length = loc.size();
    	System.err.println("CHECKING");
    	if(!(DataCollection.this instanceof MatchedSequenceDataCollection)){
    	for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
    		HaplotypeEmissionState nxt = (HaplotypeEmissionState) it.next();
    		PseudoDistribution[] ems = nxt.emissions;
    		for(int i=0; i<ems.length; i++){
    			((IlluminaRDistribution)ems[i]).divide(cumR);
    		}
//    		((IlluminaRDistribution)ems[ems.length-1]).divide(rem);
    	}
    	}
    	todrop = new TreeSet<Integer>();
    }
    this.indiv =sublist( indiv,dToInc);
    todrop.addAll(findLowDepth());
    todrop.addAll(this.findLowBAF());
    Constants.parentObj(this.indiv);
  //  todrop.addAll(this.findHets(Constants.parentObj.keySet().toArray()));
    {
    	 this.drop(new ArrayList(todrop), false);
    }
    if(Constants.CHECK){
		   checkStatesForNull();
	  }
     double[] scaleLoc = Constants.scaleLoc();
     if(scaleLoc!=null){
    	 chrToMaxIndex = new ArrayList<Integer>();
 		
 		
 		
    	 for(int i=0; i<loc.size(); i++){
    		 double start = loc.get(i);
    			double a1 = ((double)start/scaleLoc[1]);
	            int chr1 = (int)Math.floor(a1)-1; // problem if negative 
	            while(chrToMaxIndex.size()<=chr1){
	     			chrToMaxIndex.add(0);
	     		}
	          chrToMaxIndex.set(chr1, i);
    	 }
     }
    
    Logger.global.info(" ..done ");
   // if(loc.size()==0) return;
    this.makeDistributions(index);
    Set<String> toDelete = new HashSet<String>();
     if(Constants.loess(index) || Constants.gc(index)){
         toDelete.addAll(this.applyLoess(zf, platem, Constants.loess(index), Constants.gc(index),pref));
     }
    /* if(platem.size()>0){
    	List<String> plateToInclude = Arrays.asList(Constants.platesToInclude(index));
    	List<String> plateToExclude = Arrays.asList(Constants.platesToExclude(index));
    	Set<String> plates = new HashSet<String>();
    	 if(plateToInclude.contains("all"))
    	 {
    		plates.addAll(platem.values()); 
    	 }
    	 else{
    		 plates.addAll(plateToInclude);
    	 }
    	 if(!plateToExclude.contains("null")){
    		 plates.removeAll(plateToExclude);
    	 }
    	 for(Iterator<String> it = this.indiv.iterator(); it.hasNext();){
    		 String nxt = it.next().split(Constants.duplSep)[0];
    		 
    		 if(!plates.contains(platem.get(nxt))){
    			 toDelete.add(nxt);
    			 System.err.println("dropping "+nxt+" "+platem.get(nxt));
    		 }
    	 }
     }*/
     if(Constants.median_correction(index) && !(chrom.startsWith("X") || chrom.startsWith("Y"))){
         this.applyMedianCorrection(zf,dToInc);
     }
     this.applyVarianceThreshold(zf, 
    		 Constants.standardiseVariance(index) 
    		 && !(chrom.startsWith("X") || chrom.startsWith("Y")), dToInc,sample_id);
  if(length!=loc.size()) {
      throw new RuntimeException("!!");
  }
  
  //  this.calculateMLGenotypeData(true);
  this.setMinMaxRValues();
 
   //this.calculateMaf(true);
   zf.close();
  
    Logger.global.info("free memory "+Runtime.getRuntime().freeMemory());
    this.dropIndiv(toDelete.toArray(new String[0]));
   // Collections.sort(indiv);
    try{
    	if(strand.size()==0  &&  loc.size()>0
    			&& Constants.strand(index).equals("-")){
    		 for(int i=0; i<this.strand.size(); i++){
    			 flipStrand(i);
    		 }
    	}else{
   for(int i=0; i<this.strand.size(); i++){
	   if(strand.get(i)!=null && !strand.get(i)){
		   if(strand_represents_ref ){
			   //if(true) throw new RuntimeException("!!");
			   this.swapAlleles(i);
		   }
		   else{
			   this.flipStrand(i);
		   }
	   }
   }
    	}
    }catch(Exception exc){
    	exc.printStackTrace();
    }
  // EmissionState state = this.dataL.values().iterator().next();
   if(alternative.size()>0){
   	 this.rename(alternative);
  	/* for(int i=0; i<indiv.size(); i++){
  		 String alt = alternative.get(indiv.get(i));
  		 if(alt!=null){
  			 indiv.set(i, alt);
  			 String pl = platem.remove(indiv.get(i));
  			 platem.put(alt, pl);
  		 }
  		 
  	 }*/
   }
 
   if(!Constants.duplicates(index).startsWith("separate")){
	   List<EmissionState> toRemove = new ArrayList<EmissionState>();
	   String end1 = Constants.duplSep+"1";
	   for(int k=0; k<indiv.size(); k++){
		   String indv =indiv.get(k); 
		   if(indv.endsWith(end1)){
			   String first = indv.substring(0,indv.length() -Constants.duplSep.length()-1);
			   int orig = indiv.indexOf(first);
			   List<EmissionState> l = new ArrayList<EmissionState>();
			   add(l,dataL.get(first) );
			   add(l, dataL.get(indv));
			   for(int k1=0; (orig = indiv.indexOf(first+Constants.duplSep+k1))>=0; k1++){
				   add(l, dataL.get(first+Constants.duplSep+k1));
			   }
			   if(l.size()>1){
				   if(Constants.duplicates(index).startsWith("first")||Constants.duplicates(index).startsWith("merge")){
					 toRemove.addAll(l.subList(1, l.size()));
				   }
				   else if(Constants.duplicates(index).startsWith("last")  ){
					   toRemove.addAll(l.subList(0,l.size()-1));
				   }
				   if(Constants.duplicates(index).equals("merge")){
					   HaplotypeEmissionState fir =(HaplotypeEmissionState) l.get(0);
					   if(fir!=null){
						   PseudoDistribution[] st = (fir).emissions;
						   PseudoDistribution[][] l1 = new PseudoDistribution[l.size()][];
						   for(int j=0; j<l.size(); j++){
							  l1[j] =  ((HaplotypeEmissionState)(l.get(j))).emissions; 
						   }
						   for(int j=0; j<st.length; j++){
							   EmissionStateSpace emstsp =  Emiss.getSpaceForNoCopies(ploidy.get(k));
							  st[j] =  new lc1.stats.CompoundDistribution(l1[0][j], l1[1][j],emstsp );
							  for(int j1=2; j1<l.size(); j1++){
								 ((lc1.stats.CompoundDistribution) st[j]).addDist(l1[j1][j], emstsp);
							  }
						   }
					   }
				   }
			   }
		   }else if(indv.indexOf("###")>=0){
			   System.err.println(indv);
		   }
	   }
	   this.dropIndiv(getNames(toRemove.toArray(new EmissionState[0])));
	  this.subsample(Constants.thin(index));
	  if(failed.size()>0) this.drop(failed, false);
	  if(Constants.removeIndels()){
		  this.drop(Constants.indexOf(this.alleleA, 'I'), false);
		  this.drop(Constants.indexOf(this.alleleB, 'I'), false);
		  this.drop(Constants.indexOf(this.alleleA, 'N'), false);
		  this.drop(Constants.indexOf(this.alleleB, 'N'), false);
		  this.drop(Constants.indexOf(this.alleleA, '.'), false);
		  this.drop(Constants.indexOf(this.alleleA, null), false);
		  
	  }
	  /*if(false & Constants.removeAmbiguousStrand()){
		  List stra = this.strand;
		  this.drop(Constants.indexOf(stra,null), false);
	  }*/
	//	EmissionState st = getData("1000412");
	  if(Constants.dropMonomorphic(index)){
		 // if(index>0 || this.name.toLowerCase().indexOf("haplo")<0) throw new RuntimeException("!! "+this.name);
		  //removing identical sites
		  this.drop(Constants.indexOf(this.baf, 0.0, 1e-7), false);
		  this.drop(Constants.indexOf(this.baf, 1.0,1e-7), false);
		  this.drop(findMonomorphic(), false);
	  }
	  if(this.length==0){
		  Logger.global.warning("REMOVED ALL SITES!!!!!!!");
	  }
	  if(this.dc!=null) this.addMixture();
   }
}

private Collection<? extends Integer> findHets(Object[] array) {
	HaplotypeEmissionState[] hes = new HaplotypeEmissionState[array.length];
	for(int k=0; k<hes.length; k++){
		hes[k] = (HaplotypeEmissionState)this.dataL.get(array[k].toString());
	}
	List l = new ArrayList<Integer>();
	EmissionStateSpace emstsp = hes[0].emissions[0].getEmissionStateSpace();
	outer: for(int i=0; i<loc.size(); i++){
		for(int k=0; k<hes.length; k++){
			ComparableArray comp = (ComparableArray) emstsp.get(hes[k].getBestIndex(i));
			if(comp.het()) {
				l.add(i);
				continue outer;
			}
		}
	}
	
	return l;
}
protected Collection<? extends Integer> findLowDepth() {
	return new ArrayList<Integer>();
}
protected Collection<? extends Integer> findLowBAF() {
	List l = new ArrayList<Integer>();
	if(Constants.excludeBafThresh()>0){
	
	for(int k=0; k<this.baf.size(); k++){
		if(baf.get(k) < Constants.excludeBafThresh() || 1-baf.get(k) < Constants.excludeBafThresh()){
			l.add(k);
		}
	}
	}
	return l;
}
private int[]  checkOrder() {
	List<Integer> loc1 = new ArrayList<Integer>(loc);
  Collections.sort(loc1);
  if(!loc1.equals(loc)){
	  int[] alias = new int[loc.size()];
	  for(int i=0; i<loc1.size(); i++){
		  alias[i] = loc.indexOf(loc1.get(i));
	  }
	  return alias;
  }
  return null;
	
}
private List reorder(List l, int[] alias){
   if(l.size()>0){
   List l1 = new ArrayList();
   for(int k=0; k<alias.length; k++){
	   l1.add(l.get(alias[k]));
   }
   return l1;
   }else return l;
}
private int getFirstNonNumeric(String nme) {
	 Pattern p = Pattern.compile("\\D");
	 Matcher mat  = p.matcher(nme);
	 mat.useAnchoringBounds(false);
	 boolean matc = mat.find(nme.startsWith("chr") ? 3 : 0);
	 if(matc) return  mat.start();
	 else return -1;
}
private List thin(List baf, int cumR) {
	int l = (int)Math.floor((float)baf.size()/(float) cumR);
	List res = new ArrayList(l);
	for(int i=0; i < l; i++){
		res.add(baf.get(i*cumR));// + (int)Math.floor(((float)cumR)/2.0)));
	}
	return res;
}
private HaplotypeEmissionState getData(String string) {
	for(Iterator<String> it = dataL.keySet().iterator(); it.hasNext();){
		String key = it.next();
		if(key.indexOf(string)>=0) return (HaplotypeEmissionState) dataL.get(key);
	}
	return null;
}
private List<Integer> findMonomorphic() {
	int[] bcount = new int[this.loc.size()];
	int[] count =  new int[this.loc.size()];
	for(Iterator<EmissionState> ems = this.dataLvalues();ems.hasNext();){
		HaplotypeEmissionState nxt = ((HaplotypeEmissionState)ems.next());
		for(int i=0; i<bcount.length; i++){
			PseudoDistribution dist = nxt.emissions[i];
			if(dist instanceof IntegerDistribution){
				bcount[i]+=((IntegerDistribution)dist).noB;
				count[i]+=((IntegerDistribution)dist).noCop;
			}
		}
	}
	List<Integer> res = new ArrayList<Integer>();
	for(int k=0; k<count.length; k++){
		if(count[k]>0 && (bcount[k]==0 || bcount[k]==count[k] )){
			res.add(k);
		}
		//else res.add(k);
	}
	return res;
}
protected void setMinMaxRValues() {
	 if(dc!=null) {
		 ((DistributionCollection)dc).setMinMax(Constants.minR(index), Constants.maxR(index),Constants.minB(index), Constants.maxB(index));
	 }
	
}
private void checkStatesForNull() {
	for(Iterator<EmissionState> states = this.dataLvalues(); states.hasNext();){
		   HaplotypeEmissionState state = (HaplotypeEmissionState)states.next();	
		for(int i=0; i<state.emissions.length; i++){
		  if(state.emissions[i] instanceof IlluminaRDistribution && ((IlluminaRDistribution)state.emissions[i]).r()==null){
	  throw new RuntimeException("!!");
		  }
		}
	}
	
}
private void modifySampList(List<String[]> phenoBoth1, File pheno, List<String> toinc, boolean add ) throws Exception{
	List<String[]> phenoBoth = Constants.read1(phenoBoth1);
    if(phenoBoth.size()>0 && pheno.exists() && pheno.length()>0){
   	 BufferedReader br = new BufferedReader(new FileReader(pheno));
   	 String st = br.readLine();
   	 List<String> header = Arrays.asList(st.split("\t"));
   	 int pat_id = header.indexOf("PATIENT");
   	 int[] incl_inds =  new int[phenoBoth.size()] ;
   	 List<String>[] incl_val = new List[phenoBoth.size()];
   	 boolean[] all = new boolean[phenoBoth.size()];
   	 boolean[] none = new boolean[phenoBoth.size()];
   	 boolean[] lt = new boolean[phenoBoth.size()];
   	 boolean[] gt = new boolean[phenoBoth.size()];
   	 boolean[] perc = new boolean[phenoBoth.size()];
   	 boolean anyperc = false;
   	 for(int i=0; i<incl_inds.length; i++){
   		 String[] strs = phenoBoth.get(i);
   		 if(strs==null || strs[0]==null || strs[0].equals("null")){
   			 none[i] = true;
   		 }
   		 else if(strs[0].equals("all")){
  			 all[i] = true;
  		 }else{
  			 incl_inds[i] = header.indexOf(strs[0]);
  			 if(incl_inds[i]<0) incl_inds[i] = header.indexOf(strs[0].toLowerCase());
  			 if(incl_inds[i]<0){
  				 throw new RuntimeException("did not find "+strs[0]);
  			 }
  			 incl_val[i] = Arrays.asList(strs).subList(1, strs.length);
  			 if(incl_val[i].size()==1){
  				 if(incl_val[i].get(0).startsWith("<")){
  					 lt[i] = true;
  					 incl_val[i].set(0, incl_val[i].get(0).substring(1));
  				 }
  				 else if(incl_val[i].get(0).startsWith(">")) {
  					 gt[i] = true;
  					 incl_val[i].set(0, incl_val[i].get(0).substring(1));
  				 }
  				 if(incl_val[i].get(0).endsWith("%")){
  					 perc[i] = true;
  					 anyperc = true;
  					 String st1 = incl_val[i].get(0);
  					 incl_val[i].set(0, st1.substring(0,st1.length()-1));
  				 }
  			 }
  		 }
   	 }
   	 if(anyperc){
   		 List<Double>[]vals = new ArrayList[incl_inds.length];
   		 for(int k=0; k<incl_inds.length; k++){
   			 vals[k]  = new ArrayList<Double>();
   		 }
   		 while((st = br.readLine())!=null){
   	   		String[] str = st.split("\t");
   	   		for(int j=0;j<incl_inds.length; j++){
   	   			vals[j].add(Double.parseDouble(str[incl_inds[j]]));
   	   		}
   		 }
   		 br.close();
   		for(int k=0; k<incl_inds.length; k++){
  			List<String> l =  incl_val[k];
  			List<Double> val = vals[k];
  			Collections.sort(val);
  			for(int j=0; j<l.size(); j++){
  				
  				int pos =(int) Math.round((Double.parseDouble(l.get(j))/100.0)*(double) val.size());
  				l.set(j, ""+val.get(Math.min(val.size()-1, pos)));
  			}
  		 }
   		br = new BufferedReader(new FileReader(pheno));
        st = br.readLine();
   	 }
   	 while((st = br.readLine())!=null){
   	//	 st = st.replaceAll("NA", "NaN");
   		
   		String[] str = st.split("\t");
   		for(int j=0;j<incl_inds.length; j++){
   			if(!none[j] && 
   					(all[j] || 
   						(!lt[j] && !gt[j] && incl_val[j].contains(str[incl_inds[j]])) ||
   						(lt[j] && !Double.isNaN(Double.parseDouble(str[incl_inds[j]])) && Double.parseDouble(str[incl_inds[j]]) <  Double.parseDouble(incl_val[j].get(0))) || 
   						(gt[j] &&  !Double.isNaN(Double.parseDouble(str[incl_inds[j]])) && Double.parseDouble(str[incl_inds[j]]) >  Double.parseDouble(incl_val[j].get(0)))
   						)){
   					if(add ){
   						if(!toinc.contains(str[pat_id]))toinc.add(str[pat_id]);
   					}
   					else{
   						toinc.remove(str[pat_id]);
   					}
   			}
   		}
   		
   	 }
    }
}

private List<String> expand(List<String> st){
	List<String> indiv = this.indiv;
	Set<String>st1 = new HashSet<String>();
	for(int k=0; k<st.size();k++){
		String str = st.get(k);
		int ind = str.indexOf('*');
		if(ind>=0){
			String pref = str.substring(0,ind);
			String suff = str.substring(ind+1);
			for(int j=0; j<this.indiv.size(); j++){
				String str1 = indiv.get(j);
				if((str1.startsWith(pref) && pref.length()>0 )||  (suff.length()>0 && str1.endsWith(suff))) st1.add(str1);
			}
		}else{
			st1.add(str);
		}
	}
	return new ArrayList<String>(st1);
}
protected List<Integer> getSamps() throws Exception{
	List<Integer> dToInc  = new ArrayList();  

//<<<<<<< .mine
    List<String> todel = modify(Arrays.asList(Constants.toDel(this.index).split(";")));
     List<String> toinc = modify(new ArrayList<String>(Arrays.asList(Constants.toInclude(this.index).split(";"))));
   
/*=======
    List<String> todel = expand(Arrays.asList(Constants.toDel(this.index).split(";")));
     List<String> toinc = expand(new ArrayList<String>(Arrays.asList(Constants.toInclude(this.index).split(";"))));

>>>>>>> .r224*/
    if(toinc.size()==1 && toinc.get(0).equals("all")){
    	 toinc.remove(0);
    	 toinc.addAll(indiv);
     }
     
     if(toinc.size()==1 && toinc.get(0).equals("null")){
    	 toinc.remove(0);
     }
     File pheno = new File(this.dir,"pheno.txt");
     File plate = new File(this.dir,"plate.txt");
     modifySampList(Arrays.asList(Constants.phenoToInclude(this.index)), pheno, toinc, true);
     modifySampList(Arrays.asList(Constants.platesToInclude(this.index)), plate, toinc, true);
     if(!todel.get(0).equals("null")){
    	 toinc.removeAll(todel);
     }
     modifySampList(Arrays.asList(Constants.phenoToExclude(this.index)), pheno, toinc, false);
     modifySampList(Arrays.asList(Constants.platesToExclude(this.index)), plate, toinc, false);
     for(int k=0; k<toinc.size(); k++){
    	int ind = indiv.indexOf(toinc.get(k));
    	if(ind>=0) dToInc.add(ind);
     }
    
   
   if(dToInc.size()>indiv.size()  || dToInc.size()==0){
	   
	   throw new RuntimeException("!!\nd\n"+dToInc+"\nindiv\n"+indiv);  
   }
   System.err.println("SAMPLES "+dToInc.size());
     return dToInc;
}
private void modify(String st, List<String> n){
	if(st.toLowerCase().startsWith("usedef")){
		String[] str = Constants.useDataAsModel;
		for(int k=0; k<str.length; k++){
			tree.expand(n, str[k]);
		}
	}
	if(st.endsWith("*")){
		String str = st.substring(0,st.length()-1);
		if(this.tree!=null)	this.tree.expand(n,str );
		else{
			for(Iterator<String> it = this.indiv.iterator();it.hasNext();){
				String nxt = it.next();
				if(nxt.startsWith(str)){
					n.add(nxt);
				}
			}
		}
	}else if(st.startsWith("*")){
		String str = st.substring(1,st.length());
		for(Iterator<String> it = this.indiv.iterator();it.hasNext();){
			String nxt = it.next();
			if(nxt.endsWith(str)){
				n.add(nxt);
			}
		}
	}
		else{
	
		n.add(st);
	}
}
private List<String> modify(List<String> l1) {
	List<String> l = Constants.read(l1, this.index);
	List<String> n = new ArrayList<String>();
	File f = null;
	for(int i=0; i<l.size(); i++){
		String st = l.get(i);
		modify(st, n);
	}
	return n;
}
private List<String> sublist(List<String> indiv2, List<Integer> toinc) {
	List<String> l = new ArrayList<String>();
	for(int k=0; k<toinc.size(); k++){
		int j1 = toinc.get(k);
		l.add(indiv2.get(j1));
	}
	return l;
}
private File findFile(File f, final int[][] mid) {
	 String chr1 = f.getName().split("\\.")[0];
	 final String chr = chr1.endsWith("p")  || chr1.endsWith("q")  ?  chr1.substring(0, chr1.length()-1) : chr1;
	File[] f1 = f.getParentFile().listFiles(new FileFilter(){

		@Override
		public boolean accept(File arg0) {
			boolean eq =  arg0.getName().startsWith("chr"+chr+"_");
			if(eq){
				String[] str = arg0.getName().substring(3).split("_");
				int start = Integer.parseInt(str[1]);
				int end = Integer.parseInt(str[2].split("\\.")[0]);
				
				for(int k=0; k<mid.length; k++){
					if(Math.min(end - mid[k][0],mid[k][1] - start)>=0){
						return true;
					}
				}
			}
			return false;
		}
		
	});
	if(f1.length>0){
		return f1[0];
	}
	else return f;
		//f.getName().substring(3).split("_")[0];
	
}
private static int indexOf(List<String> header_snp2, String[] strings) {
	for(int k=0; k<strings.length; k++){
		int i = header_snp2.indexOf(strings[k]);
		if(i>=0) return i;
	}
	return -1;
}
protected static int indexOf1(List<String> header_snp2, String[] strings) {
	for(int k=0; k<strings.length; k++){
		for(int j=0; j<header_snp2.size(); j++){
		if(header_snp2.get(j).indexOf(strings[k])>=0) return j;
		}
	}
	return -1;
}
private String[] getNames(EmissionState[] array) {
	String[] res = new String[array.length];
	for(int k=0; k<res.length; k++){
		res[k] = array[k].getName();
	}
	return res;
}
private static void add(List l, Object o){
	   if(o!=null){
		   l.add(o);
	   }
}
Map<String, Double> karyo = null;

public  Map<String,Double> readKaryotypes()  {
	if(karyo!=null) return karyo;
	karyo = new HashMap<String, Double>();
	double[] pos = new double[2];
	int pos1 = (int)Math.floor(Constants.decode(loc.get(loc.size()-1), pos, false));
	int pos2 = (int) Math.floor(Constants.decode(loc.get(0), pos, false));
	boolean singlechrom  =false;
	if( pos1==pos2){
		singlechrom = true;
	
	}
	try{
	{
	File dirp = new File(Constants.baseDir);
	File[] karyofiles = dirp.listFiles(new FileFilter(){

		@Override
		public boolean accept(File pathname) {
			return pathname.getName().indexOf("karyo")>=0 && pathname.getName().indexOf(Constants.build(0).split("\\.")[0])>0;
		}
		
	});
	if(singlechrom && Constants.annotateMB()>0){
		for(int k=0; k<300; k+=Constants.annotateMB()){
			karyo.put(k+"", Constants.recode(new double[] {pos1, (k*1e6)}));
		}
	}
	else if(karyofiles.length>0){
		
	BufferedReader br = new BufferedReader(new FileReader(karyofiles[0]));
	String st = "";
	while((st = br.readLine())!=null){
		String[] str = st.split("\t");
		String txt = str[0];
		int chrind = Integer.parseInt(str[0].toLowerCase().replaceAll("x", "23").replaceAll("y", "24").replaceAll("mt", "24"));
		
		if(!singlechrom || str.length==2){
			int posi = Integer.parseInt(str[1]);
			double v = Constants.scaleLoc==null ?posi : Constants.recode(new double[] {chrind,posi});
			karyo.put(txt, v);
		}else{
		String[] strs = str[2].substring(1,str[2].length()-1).split(",");
		for(int k=0;k<strs.length; k++){
			String[] strs1 = strs[k].split("=");
			String band = txt+strs1[0];
			String[] stend = strs1[1].split("-");
			double mid = (Double.parseDouble(stend[0])+Double.parseDouble(stend[1]))/2.0;
			double v1= Constants.scaleLoc==null ? mid : Constants.recode(new double[] {chrind,mid});

			karyo.put(band,v1);
		}
		}
	}
	br.close();
	}
	}
	}catch(Exception exc){
		exc.printStackTrace();
	}

return karyo;
}


public static File getKaryoFile(File f, int pos)  throws Exception{
	String buildf = Constants.build(0);
	Pattern p  = Pattern.compile("[0-9]");
	Pattern p1  = Pattern.compile("[a-zA-Z]");
	Matcher m = p.matcher(buildf);
	m.useAnchoringBounds(false);
	//boolean m1 = mat.find();
	if(m.find()){
		buildf = buildf.substring(m.start());
	}
	m = p1.matcher(buildf);

	if(	m.find()){
		buildf = buildf.substring(0,m.start());
	}
	final String buildf1 = buildf;
	String nme = f.getName();
	int cent = Constants.getKaryoFile(f.getParentFile());
	String pq = pos < cent ? "p" : "q";
	 int ind = nme.indexOf('.');
		f = new File(f.getParentFile(), nme.substring(0,ind)+pq+nme.substring(ind));

	
		 return f;
	
	
}
public static BufferedReader getBuildReader(File dir, File bf, ZipFile zf) throws Exception{
	BufferedReader buildF=null;
	 if(bf!=null){
         if(!exists(bf)){
    	       bf = new File(dir, bf.getName());
    	     //  System.err.println("h");
    	   }
    	   }
        if(bf==null || !exists(bf)){
        	return null;
        }
        else{
        	buildF =    DataCollection.getBufferedReader(bf);
        }
        return buildF;
}

public void addAlleles(int i, List<String> l, List<String>b){
//	if(true) return ;
	Character A=null;
	Character B= null;
	double av=0;
	outer: for(int k=0; k<l.size(); k++){
		String st  = l.get(k);
		if(st.indexOf('N')>=0 || st.indexOf('-')>=0 || st.indexOf('X')>=0 || st.indexOf('_')>=0 || st.indexOf('Y')>=0|| st.indexOf('Z')>=0)  continue;
		char[] ch = l.get(k).toCharArray();
		for(int j=0; j<ch.length; j++){
			if(A==null){
				A = ch[j];
				if(b!=null) av = Double.parseDouble(b.get(k));
			}
			else if(B==null && ch[j]!=A){
				B = ch[j];
				break outer;
			}
		}
	}
	if(B==null) B = 'N';
	if(A==null) A = 'N';
	if(av>0.5){
		Character tmp = A;
		A=B;
		B=tmp;
	}
	if(i<alleleA.size()){
		this.alleleA.set(i, A);
		this.alleleB.set(i,B);
	}else{
		this.alleleA.add(A);
		this.alleleB.add(B);
	}

	
}
public void process(String[] str, int i,int no,  int loc_index, int[] maf_index,int chr_index, int strand_index,int snp_index,
		List l, List chr, List majorAllele, List alleleB, List forward, int bin_index){
	  l.add(loc_index<0 ? i : str[loc_index]);
      loc.add(no);
      if(chr_index>=0) chr.add(str[chr_index].startsWith("chr") ? str[chr_index].substring(3) : str[chr_index]);
      snpid.add(str[snp_index]);
      if(maf_index!=null &&maf_index[0]>=0 && maf_index[1]>=0   && str.length>maf_index[0]){
    	  char ch0 = str[maf_index[0]].charAt(0);
    	  char ch1 = str[maf_index[1]].charAt(0);
    	  if(str[maf_index[0]].length()>1 ||str[maf_index[1]].length()>1 ){
    		  int comp = (new Integer(str[maf_index[0]].length())).compareTo(str[maf_index[1]].length());
    	//	  if(comp==0) throw new RuntimeException("prob with alleles "+Arrays.asList(str));
    		  majorAllele.add(comp>0 ? 'D' : 'I');  alleleB.add(comp>0 ? 'I' : 'D');
    	  }
    	  else if(ch0=='-'){
    		majorAllele.add('D'); alleleB.add('I'); //note, we are using C,G as complementary alleles to represent I
    		// and D.
    	}else if(ch1=='-'){
    		majorAllele.add('I'); alleleB.add('D');
    	}
    	else{
    		majorAllele.add(ch0); alleleB.add(ch1);
    	}
      }else if(maf_index!=null && maf_index[0]>=0 && str.length > maf_index[0]){
    	  String str_ = str[maf_index[0]];
    	  if(str_.length()==0) majorAllele.add(null);
    	  else{char ch0 = str_.charAt(0);
    	  majorAllele.add(ch0=='-' ? 'D' : ch0);
    	  }
    	  alleleB.add(null);
      }
      else if(maf_index!=null && maf_index[1]>=0 && str.length > maf_index[1]){
    	  String str_ = str[maf_index[1]];
    	  char ch0 = str_.charAt(0);
    	  alleleB.add(ch0=='-' ? 'D' : ch0);
    	  majorAllele.add(null);
      }
      if(strand_index>=0){// && strand_index < str.length){
    	  char ch = str[strand_index].charAt(0);
    	  if(ch=='+'|| str[strand_index].toLowerCase().equals("true")
    		//  || ch=='D'
    			  ) forward.add(true);
    	  else if(ch=='-'|| str[strand_index].toLowerCase().equals("false")
    		//  || ch=='I'
    			  ) forward.add(false);
    	  else forward.add(null);
      	
      }
}


public  void readBuildFile(ZipFile zf, String prefix, BufferedReader br, String chrom1, final int [][] fromTo,List<Integer> loc, List<String> chr,  List<String> snpid, 
        List<Character> majorAllele, List<Character> alleleB,
        List<Boolean> forward,
        int loc_index, int chr_index, int snp_index, int strand_index, int bin_index, int[] maf_index,
        Collection<String> snp_ids_to_restrict) throws Exception{
    List<String> l = new ArrayList<String>();
  
    
    	
   
    String[] todrop = Constants.toDropPrefix();
    for(int i=0; i<fromTo.length; i++){
    	if(fromTo[i][0]> fromTo[i][1]) throw new RuntimeException("!! from is greater than end "+fromTo[i][0]+" "+fromTo[i][1]);
    }
    boolean drop = todrop.length>0;
    String st = "";
    Set<String> done = new HashSet<String>();
    Set<Integer> done1 = new HashSet<Integer>();
    String chrom = chrom1;
   
    int st_ind = this.getFirstNonNumeric(chrom);
    String chrom_nop = st_ind>=0 ? chrom.substring(0,st_ind) : chrom;
    
   /* if(chrom.endsWith("p") || chrom.endsWith("q")){
    	chrom_nop = chrom.substring(0,chrom.length()-1);
    }*/
   outer: for(int i=0;(st = br.readLine())!=null; i++){
        String[] str = st.split("\t+");
        if(i==0 && chr_index>=0 && ! str[chr_index].startsWith("chr") && chrom1.toLowerCase().startsWith("chr")){
        	chrom = chrom1.substring(3);
        }
        if(snp_ids_to_restrict!=null && !snp_ids_to_restrict.contains(str[snp_index])){
        	continue outer;
        }
       if(drop){
    	   for(int k=0; k<todrop.length; k++){
    		   if(snp_index>=0 && str[snp_index].startsWith(todrop[k])){
    			   continue outer;
    		   }
    	   }
       }
       if(maf_index[1]>=0){
    	   if(str[maf_index[1]].indexOf(",")>=0){
    		   continue outer;
    	   }
       }
       String id = str[snp_index];
      String chr_id = str[chr_index];
     int pind =chr_id.startsWith("chrpcs") ? -1 :  Math.max(chr_id.indexOf('p'), chr_id.indexOf('q'));
     if(pind>=0) chr_id = chr_id.substring(0,pind);
//      if(chr_id.endsWith("p")) chr_id = chr_id.substring(0)
        if( (chr_id.equals(chrom_nop) || chr_id.equals(chrom)) 		
   ){
            int no = loc_index<0 ? i*Constants.lengthMod() : Integer.parseInt(str[loc_index]);
            for(int k=0; k<fromTo.length; k++){
             
                if(no >= fromTo[k][0] ){
                	if( no <= fromTo[k][1]){
                	
                		  if(zf.getEntry(prefix+id)!=null ){
                	 if(!Constants.excludeMultiAllelicSites() ||maf_index==null ||  ((maf_index[0]<0 || str[maf_index[0]].indexOf(',')<0)  
                			 && (maf_index[1] <0 || str[maf_index[1]].indexOf(',')<0 ))){       	
                	String snp = str[snp_index];
                	  
                	if(done.contains(snp) || done1.contains(no)) {
                		Logger.global.warning("duplicate SNPS in "+name+" "+snp);
                		continue outer;
                	}
                	else 
                	{done.add(snp);
                	done1.add(no);
                	}
                	
                	this.process(str, i, no, loc_index, maf_index, chr_index, strand_index, snp_index, l, chr, majorAllele, alleleB, forward, bin_index);
                	// process(str,i);
                  
                    	
                  
                 //   System.err.println(no);
                    continue outer;
                		  }
                	}
                		 
                	
                    else{
                    	//int ind = l1.indexOf(id);
                 	   Logger.global.warning("did not find " +id+" in file");
                    }
                	 }
                	}
                }
                //else break;
            }
        }
   
   // }
   br.close();
   // return readZip(zf, l);
    
}

private boolean excludeFromTo(int[][] fromTo, int no) {
	  for(int k=0; k<fromTo.length; k++){
          
          if(no >= fromTo[k][0]){
          	if( no <= fromTo[k][1]){
          		return true;
          	}
          }
	  }
          return false;
}
/** returns if it is a probeOnly probe */

protected  Boolean process(String snpid, int i,ZipFile zf,
		List<Integer> ploidy, 
		List<Integer> sampToInc, double[] missing, double[] lrr) throws Exception {

	 
     //    System.err.println("reading "+snp_id);
         List<String> l = null;
         String snp_id = process(snpid);
       l = ApacheCompressor.getIndiv(zf, snp_id, null);
       if(l==null){
    	   l = ApacheCompressor.getIndiv(zf, snp_id+".txt", null);
    	   if(l==null){
        	   l = ApacheCompressor.getIndiv(zf, snp_id+".txt.gz", null);
           }
       }
        
        if(l!=null && l.size()!=indiv.size()){
        	List<String> l1 = new ArrayList<String>();
        	for(int ii=0; ii<l.size(); ii++){
        		String[]str = l.get(ii).split("\\s+");
        		for(int kk=0; kk<str.length; kk++){
        			l1.add(str[kk]);
        		}
        	}
        	 l = l1;
        	//l = l.subList(0, indiv.size());
        //	if(Constants.CHECK) 
        		//throw new RuntimeException("!!" +l.size()+" "+indiv.size() + " "+snp_id);
         }
       if(dp!=null) this.dp.convert(l,this.lrr_index);

      boolean  probeOnly = true;
       boolean allNull = true;
      
      
       Arrays.fill(missing,0);
         for(int j_=0;j_<sampToInc.size(); j_++){
        	 int j = sampToInc.get(j_);

        	 if(l==null){
        		((HaplotypeEmissionState) dataL.get(indiv.get(j))).emissions[i] =  
        			Emiss.getSpaceForNoCopies(ploidy.get(j)).getHWEDist1(0.0);
        	 }
        	 else{
	             String stri =l.get(j);
	             String[] st =stri.trim().split("\\s+");
	             Boolean po = this.process( indiv.get(j), header, st, i, ploidy.get(j),missing);
	         if(this.lrr_index>=0){
	        	 lrr[j_] = Double.parseDouble(st[lrr_index]);///this.avgDepth.get(j);
	         }
	             if(po!=null){
	            	 allNull = false;
	                 probeOnly = probeOnly && po;
	             }
        	 }
         }
         double no = missing[2]+missing[3];
         if(allNull || l==null){
        //	 PseudoDistribution dist = this.dataL.get(indiv.get(0)).emissions(i);
        	 this.baf.add(0.0);
        	 return null;
         }else if(no>0){
        	 this.baf.add(missing[3]/no);
        	 return probeOnly;
         }
         
        
         double baf_ = defaultBAF();
         if(baf_geno){
        	 baf_ = getBaf(l,baf_ind,'B');
         }
         if(baf_ind1>=0 && Double.isNaN(baf_)){
        	
        	 baf_ = getBaf(l,baf_ind1);
        	
         }
         if(Double.isNaN(baf_)){
        	 baf_  = calcBaf(l);
         }
         //}
        // double missingrate = this.getMissing(l,);
         this.baf.add(baf_);
         //}
       //  EmissionState ems = dataL.get("chrY.1000506.dedup.bam");
        // System.err.println("h");
       //  if(header[k].toLowerCase().indexOf("geno")>=0){
         return probeOnly;
}
protected double defaultBAF() {
	// TODO Auto-generated method stub
	return Double.NaN;
}
protected double calcBaf(List<String> l) {
	//if(true) throw new RuntimeException("calculating baf from sample");
	return Double.NaN;
}
protected double getBaf(List<String> l, int ind, char c) {
	double cnt=0;
	double bcnt =0;
	for(int k=0; k<l.size(); k++){
		String st = l.get(k).trim().split("\\s+")[ind];
		
		if(st.indexOf('N')<0){
			
			bcnt+=this.count(st, c);
			cnt+=2;
		}
	}
	return (bcnt+5)/(cnt+10);
}


private double getBaf(List<String> l, int ind) {
	double cnt=0;
	double bcnt =0;
	for(int k=0; k<l.size(); k++){
		String[] str = l.get(k).trim().split("\\s+");
		String st = str[ind];
		if(st.toLowerCase().indexOf('n')<0){
			double d = Double.parseDouble(st);
			if(!Double.isNaN(d)){
				bcnt+=d;
				cnt+=1;
			}
		}
	}
	if(cnt==0 || Double.isNaN(bcnt)){
		return 0.5;//throw new RuntimeException(" problem calculating b allele freq");
	}
	return bcnt/cnt;
}
protected String process(String snp_id) {
	// TODO Auto-generated method stub
	return snp_id;
}
public List<String> indiv(){
	if(indiv==null || indiv.size()!=dataL.keySet().size()){
		List<EmissionState> l = new ArrayList<EmissionState>(dataL.values());
		Collections.sort(l,new Comparator<EmissionState>(){

			public int compare(EmissionState o1, EmissionState o2) {
				int i1 = ((HaplotypeEmissionState)o1).dataIndex();
				int i2 = ((HaplotypeEmissionState)o1).dataIndex();
				if(i1==i2) return o1.getName().compareTo(o2.getName());
				else return i1<i2 ? -1 :1;
			}
			
		});
		indiv = new ArrayList<String>(l.size());
		for(int i=0; i<l.size(); i++){
			indiv.add(l.get(i).getName());
		}
	}
	return indiv;
}
public static boolean exists(File inclFile) {
	if(inclFile==null){
		return false;
	}
   if(inclFile.exists()){
	   return true;
   }
   File f1 = new File(inclFile.getParentFile(),inclFile.getName()+".gz");
   if((f1).exists()){
       return true;
   }
   else {
	   return false;
   }
}

public void calculateMLGenotypeData(boolean change){
    if(change) this.data.clear();
    for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
        EmissionState ld = it.next();
    
         PIGData res  = ld.getGenotypeData();
        if(change)  this.data.put(res.getName(), res);
   }
}
private Integer whichIndex(String[] strings, String string) {
    String st1 = string.toLowerCase();
   for(int i =0; i<strings.length; i++){
       if(strings[i].toLowerCase().indexOf(st1)>=0) return i;
   }
   return null;
}
public static int firstGreaterThan(List<Integer>loc, int readPos){
    
    int read =0;
    for(read=0; read <loc.size(); read++){
        if(loc.get(read)>readPos) break;
    }
    return read;
}
public static int getEndIndex(List<Integer>loc){
    int readPos = Constants.end();
    int read =0;
    for(read=loc.size()-1; read >=0; read--){
        if(loc.get(read)<readPos) break;
    }
   return read+1;
}
public  List<Integer> randPos(double d){
    List<Integer> l = new ArrayList<Integer>();
   // l.add(0);
    for(int i=0; i<this.length(); i++){
        double nxt = Constants.rand.nextDouble();
        if(nxt<d){
            l.add(i);
        }
    }
    return l;
}
// if setAsMissing is true, then we set to missing instead
public void dropRandom( 
        double d, boolean setAsMissing) {
 
   this.drop(randPos(d), setAsMissing);
    
}

public final String toString(){
	return this.name+"_"+this.size()+"_"+this.loc.size();
//   return this.loc.toString();
    
}
/** drop individuals with missing greater than thres 
public void dropMissing(double d) {
    List<String> toKeep = new ArrayList<String>();
    for(Iterator<String> it = this.getKeyIterator(); it.hasNext();){
        String st = it.next();
        
        PIGData dat = this.get(st);
        int nll = dat.countNull();
        if((double )nll < d * (double)dat.length()) toKeep.add(st);
    }
    this.restricToAlias(toKeep);
    
} */

public void dropIndiv(String[] toDel) {
	System.err.println("before "+this.dataL.size());
    System.err.println("dropping from "+this.name+" "+Arrays.asList(toDel));
   // if(indiv==null) indiv = new ArrayList<String>(this.dataL.keySet());
    int size1 = this.indiv.size();
    for(int i=0; i<toDel.length; i++){
    	if(!toDel[i].equals("null")){
    	if(indiv!=null) indiv.remove(toDel[i]);
        dataL.remove(toDel[i]);
        data.remove(toDel[i]);
    	}
    }
    int size2 = this.indiv.size();
    System.err.println("removed "+size1+" +size2 "+this.data.size()+" "+this.dataL.size());
    
    
}



public void setIndex(int i) {
    this.index = (short)i;
    for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
    	it.next().setIndex(i);
    }
    if(dc!=null) ((DistributionCollection)this.dc).setIndex((short)i);
    
}

public void calcIndiv() {
	this.indiv = new ArrayList<String>(this.dataL.keySet());
	
}

/** takes this data set (original) and maps it to coords of inferred dataset 
 * can be either maximal or minimal (largest possible aberations, or smallest possible aberations 
 * */
public void mapTo(DataCollection inferred, boolean maximal) {
	this.getMergedDeletions(maximal, true);
	
}



public int getMinDist(int pos, boolean right) {
	int min_dist = Integer.MAX_VALUE;
	if(right){
		for(int i=this.loc.size()-1;i>=0; i--){
			if(loc.get(i)<pos) break;
			int dist = Math.abs(loc.get(i) - pos);
			if(dist < min_dist){
				min_dist = dist;
			}
		}
	}
	else{
		for(int i=0; i<this.loc.size(); i++){
			if(loc.get(i)>pos) break;
			int dist = Math.abs(loc.get(i) - pos);
			if(dist < min_dist){
				min_dist = dist;
			}
		}
	}
	return min_dist;
}

public void restrictToMarkers(String[] strings) {
	if(strings==null || strings.length==1 && strings[0].equals("null")) return ;
	List<Integer> toDrop = new ArrayList<Integer>();
	Set<String> s = new HashSet<String>(Arrays.asList(strings));
	for(int i=0; i<loc.size(); i++){
		if(!s.contains(snpid.get(i))){
			toDrop.add(i);
		}
	}
this.drop(toDrop, false);
	
}

	public void dropROutliers(double[] ds) {
		
		if(ds==null) return;
	//	if(true) throw new RuntimeException("!!");
		for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
			HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
			for(int k=0; k<hes.emissions.length; k++){
				PseudoDistribution dist = hes.emissions[k];
				if(dist instanceof IlluminaRDistribution){
					((IlluminaRDistribution)dist).dropROutlier(ds);
				}
			}
		}
		
	
}
public double weight=1.0;
	public double[] weights() {
		return new double[] {weight};
	}

	public void writeAverages(File file, String[] strings) throws Exception {
		file.mkdir();
		int[] core = Constants.core();;
		PrintWriter pw_hap2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(file,Constants.experiment()+".txt"))));
		//PrintWriter pw_hap2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(outdir, chrom+".txt"))));
		// dc.writeFastphase(pw_hap2, false);
		printHapMapFormat( 
                 pw_hap2, this.indiv(), null, true, 
                 new String[] {"snpid", "loc"}, 
                 new String[0],
                  strings, "%7s", core[0], core[1]);
          //pw_hap1.close();
         pw_hap2.close();
		
	}


public DataCollection[] split(File f,String type1, int maxG){
	try{
		String[] type = type1.split(";");
		BufferedReader br = new BufferedReader(new FileReader(f));
		String st = br.readLine();
		String[]str = st.split("\t");
		int[] col = new int[type.length];
		int ind = Arrays.asList(str).indexOf("PATIENT");
		for(int i=0; i<type.length; i++){
			col[i] = Arrays.asList(str).indexOf(type[i]);
			if(col[i]<0) throw new RuntimeException("can't split with "+type[i]+"\n Should be one of "+Arrays.asList(str));
		}
		System.err.println("HEADER "+Arrays.asList(str));
		System.err.println("type is "+type);
		Map<String, String>[] m = new Map[col.length];//new HashMap<String, String>();
		List<String>[] vals = new List[col.length];
		for(int k=0; k<m.length; k++){
			m[k] = new HashMap<String, String>();
		}
		while((st = br.readLine())!=null){
			if(st.length()==0) continue;
			str = st.split("\t");
			
			for(int k=0; k<col.length; k++){
			
			try{
			//	if(!str[col[k]].equals("NA") && !str[col[k]].equals("NaN")){
				
				String v; 
				if(col[k] < str.length) v =   str[col[k]];
				else{
					System.err.println("WARNING string too short "+Arrays.asList(str));
					v = "NaN";
				}
				
			m[k].put(str[ind], v);
			//	}
			}catch(Exception exc){
				exc.printStackTrace();
				System.err.println("problem with "+Arrays.asList(str));
				System.exit(0);
			}
			}
		}
		Map<String, Integer> m2 = null;
		List<String> types2 = null;
		for(int k=0; k<m.length; k++){
			
		 vals[k] = new ArrayList<String>(new HashSet<String>(m[k].values()));
		List[] vals_ = new List[vals[k].size()];
		for(int i=0; i<vals_.length; i++){
			vals_[i] = new ArrayList<String>();
		}
		for(Iterator<String> it = m[k].keySet().iterator(); it.hasNext();){
			String key = it.next();
			vals_[vals[k].indexOf(m[k].get(key))].add(key);
		}
		
		boolean isnumeric = true;
		try{
			Iterator<String> it = vals[k].iterator();//
			String st1 = "";
			while((st1 = it.next()).length()==0){}
			if(!st1.equals("NA"))
			Double.parseDouble(st1);
		}catch(Exception exc){
			isnumeric  = false;
		}
		Map<String, Integer> m1;
		List<String> types;
		PhenoGroup[] quantiles = Constants.quantiles(index);
		if(isnumeric && quantiles!=null && vals[k].size()>quantiles.length){
			
				 types = new ArrayList<String>();
				 m1 =   ArmitageCalculator.convertToQuantiles(m[k], types, quantiles);
				List[] l = new ArrayList[types.size()];
				for(int i=0; i<l.length; i++){
					l[i] = new ArrayList<String>();
				}
				for(Iterator <Entry<String,Integer>> it = m1.entrySet().iterator(); it.hasNext();){
					Entry<String, Integer> entry = it.next();
					l[entry.getValue()].add(entry.getKey());
				}
				
				
			
		}
		else if(vals[k].size()>maxG){
			
			//else{
					System.err.println("AGGREGATING!! "+vals+" "+maxG);
					List[] n = new List[maxG];
				    StringBuffer[] str1 = new StringBuffer[maxG];
					for(int i=0; i<n.length; i++){
						n[i] = new ArrayList();
						str1[i] = new StringBuffer();
					}
					for(Iterator<String> it = vals[k].iterator(); it.hasNext();){
						String ntxt =it.next();
						int i = Constants.nextInt(n.length);
						n[i].add(ntxt);
						str1[i].append(";"+ntxt);
					}
					m1 = new HashMap<String, Integer>();
					types = new ArrayList<String>(str1.length);
					for(int i=0; i<str.length; i++){
						String strn = str1[i].toString();
						int len = Math.min(15, strn.length());//
						types.add(strn.substring(0,len));
					}
					for(Iterator<String> it = m[k].keySet().iterator(); it.hasNext();){
						String key = it.next();
						String val = m[k].get(key);
						int i=0;
						for(;!n[i].contains(val);i++){}
						m1.put(key, i);
					}
					//return split(m1, null);
			//}
		}
		else{
			m1 = new HashMap<String, Integer>();
			types = new ArrayList<String>(new HashSet<String>(m[k].values()));
			for(Iterator<String> it = m[k].keySet().iterator(); it.hasNext();){
				String key = it.next();
				m1.put(key,types.indexOf(m[k].get(key)));
			}
		}
		if(m2==null){
			m2 =m1;
			types2 = types;
		}
		else{
			List<String> vals3 = new ArrayList<String>();
			//HashMap<String, Integer> m3 = new HashMap<String, Integer>();
			for(Iterator<String> it = types2.iterator(); it.hasNext();){
				String st1 = it.next();
				for(Iterator<String> it2 = types.iterator(); it2.hasNext();){
					String st2 = it2.next();
					vals3.add(st1+"_"+st2);
				}
			}
			
			for(Iterator<String> it = m2.keySet().iterator(); it.hasNext();){
				String key = it.next();
				m2.put(key, vals3.indexOf(types2.get(m2.get(key))+"_"+types.get(m1.get(key))));
			}
			types2 = vals3;
		}
		}
		 return split(m2, types2,type1);
	}catch(Exception exc){
		exc.printStackTrace();
		return null;
	}
	
}

public List<AssociationCalculator>[][] ac;

public HWECalculator[] getHWECalculator(){
	if(hwe_calc==null){
		SortedSet<Integer> ploidy = getPloidy();
		hwe_calc = new HWECalculator[ploidy.last()];
		for(Iterator<Integer> it = ploidy.iterator(); it.hasNext();){
			int ploi = it.next();
		    hwe_calc[ploi-1] = new HWECalculator(this,ploi);
		}
	}
	return hwe_calc;
}


public static String[] inputDirLoc = null;

public static String inputDirLoc(int i){
	if(inputDirLoc!=null){
		return inputDirLoc[i];
	}
	else return Constants.inputDir(i);
}

public File[] getPhenoFiles(){
	String inpDir = Constants.inputDirLoc[this.index];
	File f;
	int index;
	while(!(f=new File(Constants.baseDir, inpDir+"/pheno.txt")).exists() && (index = inpDir.lastIndexOf('/'))>=0){
		inpDir = inpDir.substring(0,index);
		
	}
	if(f.exists())
	return new File[] {f};
	else {
		return null;
	}
}

public List<AssociationCalculator>[][] getArmitageCalculator() {
	
	// TODO Auto-generated method stubt
	if(ac == null){
		SortedSet<Integer> ploidy = getPloidy();
		ac = new List[ploidy.last()][1];
		for(int i=0; i<ac.length; i++){
			ac[i][0] = new ArrayList<AssociationCalculator>();
		}
		 //String dirF = Constants.baseDir;
        // if(dirF.equals(".")){
         	String dirF = System.getProperty("user.dir");
        // }
        File dir = new File(dirF);
        File[] phenoF;
        File f1 = new File(dir,Constants.experimentPhenoFile());
        if(f1.exists()){
        	phenoF = new File[]{f1};
        }
        else{
        	phenoF = getPhenoFiles();
        }
       
      if(phenoF==null) return null;
		try{
		
				String stratify = Constants.stratify();
				Set<String> vals = null;
				if(stratify!=null){
					vals = new HashSet<String>();
				
				for(int i=0; i<phenoF.length; i++){
					if(phenoF[i].exists()){
					BufferedReader br = new BufferedReader(new FileReader(phenoF[i]));
					String st = br.readLine();
					int ind = Arrays.asList(st.split("\t")).indexOf(stratify);
					while((st = br.readLine())!=null){
						vals.add(st.split("\t")[ind]);
					}
					br.close();
					}
				}
				}
				for(Iterator<Integer> it = ploidy.iterator(); it.hasNext();){
					int ploi = it.next();
					if(Constants.phenoToAssoc(this.index)!=null ){
					
					if(Constants.assocTest()==AssociationCalculator.linear){
						if(stratify==null){
					      ac[ploi-1][0].add(new LinearRegressionCalculator(this, phenoF, Constants.phenoToAssoc(this.index),ploi,null, null ));
						}
						else{
							Iterator<String> it1 = vals.iterator();
							for(int k=0; k<vals.size(); k++){
								   ac[ploi-1][0].add(new LinearRegressionCalculator(this, phenoF, Constants.phenoToAssoc(this.index),ploi,stratify, it1.next() ));
							}
						}
					}
					else{
						if(stratify==null){
						ac[ploi-1][0].add(new ArmitageCalculator(this, phenoF, Constants.phenoToAssoc(this.index)[0],ploi,null, null));
						}
						else{
						Iterator<String> it1 = vals.iterator();
						for(int k=0; k<vals.size(); k++){
							ac[ploi-1][0].add(new ArmitageCalculator(this, phenoF, Constants.phenoToAssoc(this.index)[0],ploi, stratify, it1.next()));
						}
						}
					}
					for(Iterator<AssociationCalculator> it1 = ac[ploi-1][0].iterator(); it1.hasNext();){//; i<ac[ploi-1].size(); i++){
						AssociationCalculator ac1 = it1.next();
						if(ac1.keysToIndex.size()==0 || (ac1.width==0 && ac1 instanceof LinearRegressionCalculator)) {
							try{
								System.err.println(this.getKeys());
								throw new RuntimeException("no pheno ids matched "+this.name+" "+phenoF[0]);
							}catch(Exception exc){
								exc.printStackTrace();
								System.exit(0);
							}
							it1.remove();
						}
					}
					}
				//	System.err.println("h");
				//	if(Constants.calcAssoc && this instanceof MergedDataCollection){
				//		ac[ploi-1].add(new ArmitageCalculator((MergedDataCollection)this, ploi));
					//}
				
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	return ac;
}
public SortedSet<Integer> getPloidy(){
	SortedSet<Integer> ploidy = new TreeSet<Integer>();
	for(Iterator<EmissionState> it =this.dataL.values().iterator(); it.hasNext(); ){
		ploidy.add(it.next().noCop());
	}
	return ploidy;
}
HWECalculator[] hwe_calc;



public void add(DataCollection right){
	List<String> snpid = new ArrayList<String>();
	List<Integer> loc = new ArrayList<Integer>();
	if(this.loc.get(loc.size()-1)>= right.loc.get(0)) throw new RuntimeException("!!");
	this.loc.addAll(right.loc);
	this.snpid.addAll(right.snpid);
	this.alleleA.addAll(right.alleleA);
	this.alleleB.addAll(right.alleleB);
	this.strand.addAll(right.strand);

	this.probeOnly =(Boolean[]) Constants.append(probeOnly, right.probeOnly, Boolean.class);
	this.length+=right.length;
	List<String> toD = new ArrayList<String>();
	for(Iterator<String> it = this.dataL.keySet().iterator(); it.hasNext();){
		String key = it.next();
		HaplotypeEmissionState nxt = (HaplotypeEmissionState)dataL.get(key);
		HaplotypeEmissionState nxt1 =(HaplotypeEmissionState) right.dataL.get(key);
		if(nxt1!=null) {
			nxt.append(nxt1);
		}
		else{toD.add(nxt.name);}
	}
	//((DistributionCollection)this.dc).append((DistributionCollection)right.dc);
	this.dropIndiv(toD.toArray(new String[0]));

}

public void transfer(DataCollection dc2) {
	this.loc = dc2.loc;
	this.name = dc2.name;
	this.snpid = dc2.snpid;
	this.dataL = dc2.dataL;
	if(alleleA.size()>0)this.alleleA = dc2.alleleA;
	if(alleleB.size()>0)this.alleleB = dc2.alleleB;
	//if(strand.size()>0) 
		this.strand = dc2.strand;
	this.probeOnly = dc2.probeOnly;
	this.length = dc2.length;
//	this.dc = dc2.dc;
}


public DataCollection[]  split(Map<String, Integer> pheno, List<String> types, String type1){
	
	if(true){
		//this adds a null category if not included
			//Constants.quantiles(index)==null || Arrays.asList(Constants.quantiles(index)).contains("null")){
		adjust(pheno,types);
	}
	//List vals = new ArrayList(new HashSet(pheno.values()));
	DataCollection[] res = new DataCollection[types.size()];
//	Object st = pheno.get("5055302_1671121232_A");
	//int ind_ = indiv.indexOf("5055302_1671121232_A");
String restrPheno = Constants.restrictPheno(this.index);
	for(int i=0; i<res.length; i++){
		if(restrPheno!=null && !types.get(i).equals(restrPheno)) continue;
		res[i] = new DataCollection();
		res[i].loc = loc;
		res[i].baf = baf;
		res[i].snpid = snpid;
		res[i].chrom = chrom;
		res[i].probeOnly = probeOnly;
	//	res[i].maf = maf;
		res[i].alleleA = new ArrayList<Character>(alleleA);
		res[i].alleleB = new ArrayList<Character>(alleleB);
		res[i].indiv = new ArrayList<String>();
		res[i].strand = new ArrayList<Boolean>(this.strand);
		res[i].length = this.length;
		//res[i].stSp = stSp;
		//res[i].stSp1 = stSp1;
		
		res[i].name = this.name+"_"+type1+"="+types.get(i);
	if(dc!=null){
		res[i].dc =i==0 ||(Constants.sepModels()==false && false) ? dc :  new DistributionCollection((DistributionCollection)dc, true) ;//dc==null ? null : this.dc.clone();
	}
		if(Constants.modelbg()){
			res[i].bg = bg;//i==0 ? bg : (BackgroundEmissionState) bg.clone();
			res[i].bg1 = bg1;
			res[i].dataL.put(bg.name, bg);
			res[i].data.put(bg1.name, bg1);
		}
		
		//res[i].indiv.add(bg.name);
	}
	for(Iterator<String> it = this.dataL.keySet().iterator(); it.hasNext();){
		String key = it.next();
		Integer phenot = pheno.get(key);
		
		int ind =phenot;
		DataCollection resi = res[ind];
		HaplotypeEmissionState em = (HaplotypeEmissionState) this.dataL.get(key);
		resi.dataL.put(key, em);
		resi.indiv.add(key);
	}
/*	for(Iterator<Map.Entry<String, Object> > it = pheno.entrySet().iterator(); it.hasNext();){
		Map.Entry<String, Object> m = it.next();
		String key = m.getKey();
		int ind = vals.indexOf(m.getValue());
		DataCollection resi = res[ind];
		HaplotypeEmissionState em = (HaplotypeEmissionState) this.dataL.get(key);
		if(em!=null){
			// em.setIndex(ind);
			 resi.dataL.put(key, em);
		}
		else if(this.name.indexOf("SLEGEN")>=0){
			
			for(Iterator<String> it1 = this.dataL.keySet().iterator(); it1.hasNext();){
			 String key1 = it1.next();
			 inner: if(key1.startsWith(key)){
				 em = (HaplotypeEmissionState) this.dataL.get(key1);
				 resi.dataL.put(key1, em);
				 break inner;
			 }
			 
			}
			if(em==null){
				System.err.println("missing "+key);
			}
		}
		PIGData em1 =  this.data.get(key);
		if(em1!=null){
			resi.data.put(key,em1);
			resi.indiv.add(key);
		}
		
		
	}
//	if(res[0].dataL.get("65111")!=null && res[1].dataL.get("65111")!=null){
	//	throw new RuntimeException("!!");
//	}*/
	return res;

}
private void adjust(Map pheno2, List<String> types) {
	
	int null_inde = types.size();
	boolean addedNull = false;
	outer: for(Iterator<String> it1 = this.dataL.keySet().iterator(); it1.hasNext();){
		 String key1 = it1.next();
		 if(!pheno2.containsKey(key1)){
		 if(this.name.indexOf("SLEGEN")>=0){
			 if(true) throw new RuntimeException("!!");
				 for(Iterator<String> it2 = pheno2.keySet().iterator(); it2.hasNext();){
					 String key = it2.next();
					 if(key1.startsWith(key) || key.startsWith(key1)){
						pheno2.put(key1, pheno2.get(key));
						System.err.println("substituing "+key1+" "+key);
						continue outer;
					 }
				 }
			 }
			 pheno2.put(key1, null_inde);
			 Logger.global.warning("no pheno found for "+key1);
			 addedNull=true;
		 }
		 
		
		 
		}
	if(addedNull){
		types.add("null");
	}
		
	
}
/** first is array of R dists, second B dists */
public Object[] probDists() {
	// TODO Auto-generated method stub
	return null;
}

public void initialisationStep() {
	// TODO Auto-generated method stub
if(bg!=null)
		this.bg.initialiseCounts();
	
}
public EmissionState bg() {
	return bg;
}
public double[] calcDist(Integer noCop, int i_hmm, PseudoDistribution dist) {
	return this.bg.calcCNDist(noCop, i_hmm, dist);
}
public Boolean probeOnly(int i) {
	if(probeOnly==null || probeOnly.length==0) return null;
	else return probeOnly[i];

}
public boolean hasIntensity(int i) {
	return true;
}
public List<Integer> intensityProbes() {
	List<Integer> res = new ArrayList<Integer>();
	for(int i=0; i<this.loc.size(); i++){
		if(this.hasIntensity(i)){
			res.add(i);
		}
	}
	return res;
}
public void rename(Set<String> non_unique) {
//if(true) throw new RuntimeException("!!");
	if(rename==null) rename = new HashMap<String, String>();
	for(Iterator<String> it = non_unique.iterator(); it.hasNext();){
		String key = it.next();
		if(dataL.containsKey(key)){
			EmissionState st = dataL.remove(key);
			st.name = key+"###"+this.name;
			rename.put(key, st.name);
			dataL.put(st.name, st);
			if(data.containsKey(key)){
				
				data.put(st.name, data.remove(key));
			}
		
		}
	}
	this.indiv = null;
	this.indiv();
	
}

public void rename(Map<String, String> convert) {
	for(Iterator<String> it = convert.keySet().iterator(); it.hasNext();){
		String key = it.next();
		String val = convert.get(key);
		if(dataL.containsKey(key)){
			EmissionState st = dataL.remove(key);
			st.name = val;
			dataL.put(st.name, st);
			if(data.containsKey(key)){
				
				data.put(st.name, data.remove(key));
			}
		
		}
	}
	this.indiv = null;
	this.indiv();
	
}
public static DataCollection append(DataCollection dataCollection,
		DataCollection resu) {
	DataCollection mdc = MergedDataCollection.getMergedDataCollection(new DataCollection[] {dataCollection, resu}, 
		dataCollection.name+"_"+resu.name	,null, false
	);
	return mdc;
}


public Map<String, String>[] getRenaming() {
	return new Map[]{this.rename};
}
public double getMaf(int j) {
	double[] bcount = new double[] {0,0};
	
	for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
		PseudoDistribution dist = ((HaplotypeEmissionState)it.next()).emissions[j];
		dist.getB(bcount,1);
	}
	return bcount[0]/bcount[1];//(double) dataL.size();
}
public void setDC(AbstractDistributionCollection res) {
	this.dc = res;
	
}
public boolean updateIndex(String rs) {
	// TODO Auto-generated method stub
	
	return false;
}
public boolean canProcess(String rs){
	return false;
}
public void printcompressed(List<String> keys) throws Exception {
	if(writeCompressGeno!=null)writeCompressGeno.print(keys);
	if(writeCompressStates!=null) writeCompressStates.print(keys);
}
public double baf(int i) {
	if(this.probeOnly[i]!=null && probeOnly[i]) return 0.0;
	return this.baf.get(i);
}

   
}
