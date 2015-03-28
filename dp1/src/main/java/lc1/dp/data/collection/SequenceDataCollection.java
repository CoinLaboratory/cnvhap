package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.CGH.Location;
import lc1.dp.data.representation.Emiss;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.PhasedDataState;
import lc1.stats.DepthDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.util.ApacheCompressor;
import lc1.util.Constants;

import org.apache.commons.compress.archivers.zip.ZipFile;

public class SequenceDataCollection extends LikelihoodDataCollection{

	
	
	
	public SequenceDataCollection() {
		super();
		// TODO Auto-generated constructor stub
	}

	public SequenceDataCollection(DataCollection dat) {
		super(dat);
		// TODO Auto-generated constructor stub
	}

	public List<Integer> getBackgroundCN(){
		/*	if(Constants.modelbg()){
			 return this.stSp[1].copyNumber;	
			}
			else*/
		    return Arrays.asList(new Integer[] {1});
		}
	

	@Override    
	public final void makeDistributions(int index) {
		this.dc = new DistributionCollection(getBackgroundCN(),index, this.loc.size(),this.dir);
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

	public SequenceDataCollection(File f, short index, int noCopies,
			int[][] mid, File bf,
			Collection<String> snpidrest) throws Exception {
		super(f, index, noCopies, mid, bf,snpidrest);
		
		
		// TODO Auto-generated constructor stub
	}
	
	

	

	public SequenceDataCollection(List<Integer> locs) {
		super(locs);
		// TODO Auto-generated constructor stub
	}

	public SequenceDataCollection(String[] categories) throws Exception {
		super(categories);
		// TODO Auto-generated constructor stub
	}

	@Override
	public DataC clone() {
		// TODO Auto-generated method stub
		return super.clone();
	}

	

	@Override
	public HaplotypeEmissionState createEmissionState(String key, int noCopies) {
		// TODO Auto-generated method stub
		return super.createEmissionState(key, noCopies);
	}

	@Override
	public void extractFromTrioData() {
		// TODO Auto-generated method stub
		super.extractFromTrioData();
	}

	@Override
	public double[] getR(Location loc, List<Double> l, List<Integer> l1,
			List<Integer> l2, List<Double> b, String name) {
		// TODO Auto-generated method stub
		return super.getR(loc, l, l1, l2, b, name);
	}

	@Override
	public boolean hasIntensity(int i) {
		// TODO Auto-generated method stub
		return super.hasIntensity(i);
	}

	@Override
	public void maximisationStep(double[] pseudo, int i) {
		// TODO Auto-generated method stub
		super.maximisationStep(pseudo, i);
	}

	@Override
	public void print(PrintWriter pw) {
		// TODO Auto-generated method stub
		super.print(pw);
	}

	@Override
	public boolean process(HaplotypeEmissionState data,
			HaplotypeEmissionState dataL, String header, String geno, int i, double[] miss) {
		// TODO Auto-generated method stub
		return super.process(data, dataL, header, geno, i,miss);
	}
	
	//List<String> snps_all = new ArrayList<String>();
	
	
	@Override

	protected  Boolean process(String snpid, int i,ZipFile zf,
			List<Integer> ploidy, List<Integer> sampToInc, double[] miss, double[] lrr) throws Exception {
		int counta_ind = this.indexOf1(headsLowerCase, "counta".split(":"));
		int    countb_ind = this.indexOf1(headsLowerCase, "countb".split(":"));
		if(counta_ind<0 && countb_ind<0){
			counta_ind = this.indexOf1(headsLowerCase, "ad".split(":"));
			countb_ind = counta_ind;
		}
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
	        	//l = l.subList(0, indiv.size());
	        //	if(Constants.CHECK) 
	        		throw new RuntimeException("!!" +l.size()+" "+indiv.size());
	         }
	    if(dp!=null)    this.dp.convert(l, this.lrr_index);
	      boolean  probeOnly = true;
	       boolean allNull = true;
	       this.minR.add(Double.MAX_VALUE);
	       this.maxR.add(Double.MIN_VALUE);
	       for(int j_=0;j_<sampToInc.size(); j_++){
	        	 int j = sampToInc.get(j_);

	        	 if(l==null){
	        		((HaplotypeEmissionState) dataL.get(indiv.get(j))).emissions[i] =  
	        			Emiss.getSpaceForNoCopies(ploidy.get(j)).getHWEDist1(0.0);
	        	 }
	        	 else{
		             String stri =l.get(j);
		             String[] st =stri.trim().split("\\s+");
		             Boolean po = this.process( indiv.get(j), header, st, i, ploidy.get(j),avgDepth.get(j),j);
		             if(this.lrr_index>=0){
			        	 lrr[j_] = Double.parseDouble(st[lrr_index]);
			         }
		             if(po!=null){
		            	 allNull = false;
		                 probeOnly = probeOnly && po;
		             }
		             
	        	 }
	                
	         }
	        
	         if(allNull || l==null){
	             //	 PseudoDistribution dist = this.dataL.get(indiv.get(0)).emissions(i);
	             	 this.baf.add(0.0);
	             	 return null;
	              }
	              
	             
	              double baf_ = Double.NaN; 
	              
	             try{
	             	
	             	 baf_ = getBaf(l,counta_ind, countb_ind);
	             	 
	             }catch(Exception exc){
	            	 Exception exc1 = new RuntimeException("problem with "+snpid+" "+exc.getMessage());
	            	 exc1.printStackTrace();
	            	 baf_ = 0.5;
	             }
	              //}
	             // double missingrate = this.getMissing(l,);
	              this.baf.add(baf_);
	         return probeOnly;
	}
	

	private double getBaf(List<String> l, int inda, int indb) {
		double cnt=Constants.bafPseudo();
		double bcnt =0.5*Constants.bafPseudo();
		for(int k=0; k<l.size(); k++){
			String[] str = l.get(k).trim().split("\\s+");
			if(inda>str.length){
				
				int cnta,cntb;
				if(indb==inda){
					String[] st1 = str[inda].split(",");
					 cnta = Integer.parseInt(st1[0]);
					 cntb = Integer.parseInt(st1[1]);
				}else{
			 cnta = Integer.parseInt(str[inda]);
			 cntb = Integer.parseInt(str[indb]);
			bcnt+=cntb;
			cnt+=cnta+cntb;
			}
			}
		}
		if(cnt==0 || Double.isNaN(bcnt)){
			System.err.println("WARNING TOTAL CNT IS ZERO");
			return 0.5;
//			throw new RuntimeException(" problem calculating b allele freq "+cnt+ " "+bcnt);
		}
		return bcnt/cnt;
	}

	 @Override
	    public  Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy, double[] miss){
	        boolean doneB = false;
	        boolean nullB = false;
	        boolean nullR = false;
	        IlluminaNoBg state = (IlluminaNoBg)this.dataL.get(indiv);
	      
	       for(int k=0; k<header.length; k++){
	           if(geno[k].equals("null")) geno[k] = "NaN";
	             String hk = header[k].toLowerCase();
	             if(hk.indexOf("y raw")>=0   || hk.indexOf("y_raw")>=0){
	                 Double b = Double.parseDouble(geno[k])/1000.0;
	                 if(b==null || Double.isNaN(b)) nullB  = true;
	                ((IlluminaNoBg)state).setB(i,b);
	                doneB =true;
	                if(b<Constants.minB[index]) b = Constants.minB[index];
	                 if(b>Constants.maxB[index]) b = Constants.maxB[index];
	                 
	            }
	             else if((hk.indexOf("x raw")>=0 || hk.indexOf("x_raw")>=0)){
	                 double r = Double.parseDouble(geno[k])/1000.0;
	                 if( Double.isNaN(r)) nullR  = true;
	               
	                // if(this.snpid.get(i).startsWith("cnv1950")) r+=0.6;
	                 state.setR(i,r);
	                 if(r<Constants.minR[index]) r = Constants.minR[index];
	                 if(r>Constants.maxR[index]) r = Constants.maxR[index];
	                
	                 
	             }
	             else if(false && header[k].indexOf("AgilentPred")>=0){
	              /*   EmissionStateSpace emStSp = state.getEmissionStateSpace();
	                 PhasedDataState pig = (PhasedDataState)this.data.get(indiv);
	                 if(geno[k].equals("0")){
	                     pig.emissions[i] = new IntegerDistribution(super.trans("AA"));
	                 }
	                 else{
	                     IlluminaRDistribution dist = 
	                    	 Constants.suppressB() ?
	                    			 new IlluminaRDistribution(this.index):
	                    	 new IlluminaDistribution(this.index);
	                     double r = Double.parseDouble(geno[k]);
	                     dist.setR(r);
	               //      IlluminaProbR probR = this.getR();
	                //     IlluminaProbB probB = this.getB();
	                    // double[] d = dist.calcDistribution( state.distribution, emStSp);
	                  //   int max = Constants.getMax(d);
	                   
	                   
	                  //  pig.emissions[i] = new IntegerDistribution(max);
	                 }*/
	             }
	                    }
	       if(nullR ){
	    	   ((IlluminaNoBg)state).emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);
	    	   
	       }
	       else if(!doneB
	    		   || snpid.get(i).toLowerCase().indexOf("cnv")>=0
	    		   || snpid.get(i).toLowerCase().indexOf("A_")>=0
	    		 || nullB
	    		   ){
	         Number r =  ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).r();
	        ((IlluminaNoBg)state).emissions[i] = new IlluminaRDistribution(this.index);
	        ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).setR(r.doubleValue());
	         //  throw new RuntimeException("!!");
	       }
	       try{
	    	      //  boolean doneGeno = false;
	    	        PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
	    	        for(int k=0; k<header.length; k++){
	    	        /*if(header[k].toLowerCase().indexOf("geno")>=0){
	    	            if(data!=null && data.emissions[i]==null){
	    	            int ind = trans(geno[k]);
	    	        
	    	                 data.emissions[i] = new IntegerDistribution(ind);
	    	            }
	    	        }*/
	    	        
	    	    }
//	    	        if(data.emissions[i]==null){
	  //  	            EmissionState sta = this.dataL(indiv);
	    //	             data.emissions[i] = new IntegerDistribution(sta.getBestIndex(i));
	    	//         }
	    	  //      data.emissions[i].setDataIndex(this.index);
	    	       
	    	    }catch(Exception exc){
	    	        exc.printStackTrace();
	    	    }
	    	    Boolean res =  state.emissions[i].probeOnly();
	    	 /*  if(res==null){
	    		   int pos = this.loc.get(i);
	    	   	Logger.global.info("all NC at "+i+" "+this.index+" "+pos);
	    	   }*/
	    	    return res;
	    }

	
	 //@Override
	 public Boolean process(String indiv, String[] header, String[] geno, int i,
			 int ploidy,double averageDepth, int ind_indiv) {

		 boolean doneDepth = false;
		 boolean nullDepth = false;
		 int numPoints=0;
		 HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);

		 Integer depth=null;

		 for(int k=0; k<header.length; k++){	   		
			 String hk = header[k].toLowerCase();
			 if(hk.indexOf("depth")>=0 ||hk.indexOf("dp")>=0 ){
				 if(geno[k].startsWith("nul")) geno[k] = "NaN";
				 depth = Integer.parseInt(geno[k]);
				 // if(depth>0){
				 //	 System.err.println(depth);
				 // }
				 if(depth==null || Double.isNaN(depth)) nullDepth  = true;
				 //            ((IlluminaNoBg)state).setB(i,depth);
				 doneDepth =true;
			 }
		 }
		 if(nullDepth){
			 state.emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);

		 }
		 else{
			 state.emissions[i] = new DepthDistribution(this.index,depth,averageDepth/(double)ploidy, ind_indiv);

		 }
		 return true;		




		 //		return super.process(indiv, header, geno, i, ploidy);
	 }
	
/*	//@Override
	public Boolean process(String indiv, String[] header, String[] geno, int i,
			int ploidy,double averageDepth) {

        boolean doneDepth = false;
        boolean nullDepth = false;
        int numPoints=0;
        HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);

        int positionInGeno = i-this.snpid.indexOf(this.snpid.get(i));
	   		
         String hk = header[0].toLowerCase();
         Integer depth = null;
         if(hk.indexOf("depth")==0){
        	 if(geno[positionInGeno].startsWith("nul")) geno[positionInGeno] = "NaN";
        	 depth = Integer.parseInt(geno[positionInGeno]);
        	// if(depth>0){
        	//	 System.err.println(depth);
        	// }
             if(depth==null || Double.isNaN(depth)) nullDepth  = true;
//            ((IlluminaNoBg)state).setB(i,depth);
            doneDepth =true;
        }
       
       if(nullDepth){
    	   state.emissions[i] = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);
    	   
       }
       else{
    	   state.emissions[i] = new DepthDistribution(this.index,depth,averageDepth/(double)ploidy);
    	  
       }
       return true;		
		
		
		
		
//		return super.process(indiv, header, geno, i, ploidy);
	}
	*/

	 
	 
	/*@Override
	public void process(String[] str, int i,int no,  int loc_index, int[] maf_index,int chr_index, int strand_index,int snp_index,
			List l, List chr, List majorAllele, List alleleB, List forward, int bin_index){
		// TODO Auto-generated method stub
//		super.process();
		
	  l.add(loc_index<0 ? i : str[loc_index]);
      loc.add(no);
      if(chr_index>=0) chr.add(str[chr_index].startsWith("chr") ? str[chr_index].substring(3) : str[chr_index]);
      //to keep the size of snpid right add the location to it.
      snpid.add(str[snp_index]);
      if(maf_index!=null &&maf_index[0]>=0 && maf_index[1]>=0   && str.length>maf_index[0]){
          majorAllele.add(str[maf_index[0]].charAt(0));
          alleleB.add(str[maf_index[1]].charAt(0));
          }
      if(strand_index>=0 && strand_index < str.length){
      	forward.add(str[strand_index].charAt(0)=='+');
      	
      }
	}*/

	@Override
	public void readBuildFile(ZipFile zf, String prefix, BufferedReader br, String chrom1, final int [][] fromTo,List<Integer> loc, List<String> chr,  List<String> snpid, 
	        List<Character> majorAllele, List<Character> alleleB, List<Boolean> forward, int loc_index, int chr_index, int snp_index, 
	        int strand_index, int bin_index, int[] maf_index,
	        Collection<String> snpidrest
	) throws Exception{
		
		// int[][] fromTo= getExpanded(fromTo1,50);
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
	    String chrom_nop = chrom;
	    if(chrom.endsWith("p") || chrom.endsWith("q")){
	    	chrom_nop = chrom.substring(0,chrom.length()-1);
	    }
	  
	    //int firstPos = -1;
	   if (br==null)
		   throw new RuntimeException("!! error with build file");
	   else
	   {
		   outer: for(int i=0;(st = br.readLine())!=null; i++){
			   String[] str = st.split("\\s+");
			   if(i==0 && chr_index>=0 && ! str[chr_index].startsWith("chr")){
				   chrom = chrom1.substring(3);
			   }
			   if(drop){
		    	   for(int k=0; k<todrop.length; k++){
		    		   if(snp_index>=0 && str[snp_index].startsWith(todrop[k])){
		    			   continue outer;
		    		   }
		    	   }
		       }
			   
			   String id = str[snp_index];
		        if(chr_index<0 ||  str[chr_index].equals(chrom) || str[chr_index].equals(chrom_nop))
		        {
		            int no = loc_index<0 ? i*Constants.lengthMod() : Integer.parseInt(str[loc_index]);
		            for(int k=0; k<fromTo.length; k++)
		            {
		    //        	if(no>=fromTo1[k][0]){
		      //      		if( no <= fromTo1[k][1]){
		            		//	this.snps_all.add(id);
		            	if(no >= fromTo[k][0])	{
		                	if( no <= fromTo[k][1]){
		                		if(zf.getEntry(prefix+id)!=null ){

				                	String snp = str[snp_index];
				                	if(/*done.contains(snp) || */done1.contains(no)) {
				                		Logger.global.warning("duplicate SNPS in "+name+" "+snp);
				                		continue outer;
				                	}
				                	else 
				                	{
				                		done.add(snp);
				                		done1.add(no);
				                	}
				                	this.process(str, i, no, loc_index, maf_index, chr_index, strand_index, snp_index, l, chr, majorAllele, alleleB, forward, bin_index);
				                	// process(str,i);
				                  
				                 //   System.err.println(no);
				                    continue outer;
		                		}	
		                	}
		//                }
		  //          		}
		            	}
		                //else break;
		            }
		        }
		    }
	   	br.close();
		   		}
	   }

	/*private int[][] getExpanded(int[][] fromTo, int spacing) {
		int[][] res = new int[fromTo.length][2];
		for(int k=0; k<res.length; k++){
			res[k][0] = fromTo[k][0] - Constants.extra()*spacing;
			res[k][1] = fromTo[k][1]+Constants.extra()*spacing;
		}
		return res;
	}*/

	}


	

//}
