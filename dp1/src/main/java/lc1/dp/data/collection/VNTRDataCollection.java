package lc1.dp.data.collection;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.Constants;

public class VNTRDataCollection extends SimpleDataCollection {

	
	
	static String[] header = new String[] {"genotype"};
	public VNTRDataCollection(File file, short i2, int i, int[][] mid,
			File buildF,Collection<String> snpidrest) throws Exception{
		super(file,i2,i,mid, buildF,snpidrest);
	}
	
	//Map<Integer, int[]> indices;// = new HashMap<Integer, int[]>();
	String getInd(String geno, int k){
		//String[] res = new String[Constants.alleles[index]];
		//for(int k=1; k<res.length; k++){
			char ch = Emiss.conv.get(k);
			int count = count(geno, ch);
			if(count==0){
				return  "AA";
			}
			else if(count ==1){
				return "AB";
			}
			else{
				return "BB";
			}

	//	}
		//return res;
	}
	protected String process(String snp_id) {
		// TODO Auto-generated method stub
		return snp_id.split("_")[1];
	}
	
	
	Map<String, List<String>> bins;
	//@Override
	 public Boolean  process(String indiv,  String[] header,  String[] geno, int ii, int ploidy){
    	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	 //  int[] ind = indices.get(ii);
    	   int i = ii;
    //	 snpid.get(ii));
    	   String[] snp_ = this.snpid.get(i).split("_");
    	   int allel=bins.get(snp_[1]).indexOf(snp_[0])+1;
    	 
    	try{
        	
        	
       
            for(int k=0; k<header.length; k++){
            if(header[k].toLowerCase().indexOf("geno")>=0 || header[k].toLowerCase().indexOf("plus_allele")>=0 ||  header[k].toLowerCase().indexOf("gtype")>=0){
            	if(geno[k].indexOf("null")>=0 ){
            		//for(int kk=0; kk<ind.length; kk++)
                		dataSt.emissions[i] =stsp.getHWEDist1(null);
            		//}
            		return false;
            	}
                if(dataSt.emissions[i]==null){
                	String geno1 = getInd(geno[k], allel);
                	//for(int kk=0; kk<geno1.length; kk++){
                		 int ind1 = 
                         	!this.abGenos || alleleA.size()==0 && alleleB.size()==0 ? 
                         	trans(geno1, stsp, null, null,index) :
                         		trans(geno1, stsp,alleleA.get(i), alleleB.get(i),index)
                         		;
                        
                         
                             	 dataSt.emissions[i] = stsp.getIntDist(ind1);
                	//}
                   
                        
                        
                    
                }
            }
            
        }
           
           // data.emissions[i].setDataIndex(this.index);
        }catch(Exception exc){
            exc.printStackTrace();
        }
       return false;
    }
	/*@Override
	 public Boolean  process(String indiv,  String[] header1,  String[] geno1, int i, int ploidy){
		
    	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	   String[] geno;
    	  
    	   if(geno1[0].equals("null") || geno1[1].equals("null")){
    		   geno = new String[] {"null"};
    	   }
    	   else{
    		   geno = new String[] {geno1[0]+geno1[1]};
    	   }
    	try{
        	
        	
          //  boolean doneGeno = false;
        //    EmissionStateSpace stsp = this.stSp[this.no_copies-1];
        	   boolean deletionIsNa = Constants.useDeletion(i);
        	
           
           // PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
         
            for(int k=0; k<header.length; k++){
            if(header[k].toLowerCase().indexOf("geno")>=0 || header[k].toLowerCase().indexOf("plus_allele")>=0){
            	if(geno[k].equals("null")){
            	
                		dataSt.emissions[i] =stsp.getHWEDist1(0.0);
                	
            		return false;
            	}
                if(dataSt.emissions[i]==null){
                	
                	 int ind = 
                     	!this.abGenos ? 
                     	trans(geno[k], stsp, null, null) :
                     		trans(geno[k], stsp,alleleA.get(i), alleleB.get(i))
                     		;
                    if(geno[k].equals("null")){
                    	Logger.global.info("h");
                    }
                   double soften = Constants.softenHapMap();
                   // dataSt.emissions[i] = new IntegerDistribution(ind,stsp);
                    if(deletionIsNa & stsp.getCN(stsp.getHaploPairFromHaplo(ind))==0){
                    	if(Constants.format.length==1){
                    		dataSt.emissions[i] =stsp.getHWEDist1(0.0);
                    	}
                    	else{
                    		dataSt.emissions[i] = null;
                    	}
                    }
                    else{
                        if(soften>0){
                        	
                        	 dataSt.emissions[i] = stsp.getSoftenedDist(ind);
                        }
                        else{
                        	 dataSt.emissions[i] = stsp.getIntDist(ind);
                        }
                        
                    }
                    
                }
            }
            
        }
           
           // data.emissions[i].setDataIndex(this.index);
        }catch(Exception exc){
            exc.printStackTrace();
        }
       return false;
    }*/
	//List<List<String>> bins = new ArrayList<List>();
	@Override
	public void process(String[] str, int i,int no,  int loc_index, int[] maf_index,int chr_index, int strand_index,int snp_index,
			List l, List chr, List majorAllele, List alleleB, List forward, int bin_index){
	//	int alleles = Constants.alleles(index)-1;
		//if(indices==null) indices = new HashMap<Integer, int[]>();
	if(this.bins==null) this.bins = new HashMap();
		//this.indices.put(i,ind);
		String[] bins = str[bin_index].split(":");
		this.bins.put(str[snp_index], Arrays.asList(bins));
		int[] ind = new int[bins.length];
		if(Constants.alleles==null) Constants.alleles = new int[Constants.format.length];
		if(Constants.alleles.length<=index){
			int[] all = new int[Constants.format.length];
			System.arraycopy(Constants.alleles,0,all,0,Constants.alleles.length);
		}
	//	Constants.alleles[index] = bins.length;
		for(int k=0; k<bins.length; k++){
			ind[k] = loc.size();
		  l.add(str[loc_index]+k+1);
	      loc.add(no+k+1);
	      chr.add(str[chr_index].substring(3));
	      snpid.add(bins[k]+"_"+str[snp_index]);//+"_"+));
	      majorAllele.add('A');
	      alleleB.add('B');
/*	      if(maf_index!=null &&maf_index[0]>=0 && maf_index[1]>=0   && str.length>maf_index[0]){
	          majorAllele.add(str[maf_index[0]].charAt(0));
	          alleleB.add(str[maf_index[1]].charAt(0));
	      }*/
	      if(strand_index>=0 && strand_index < str.length){
	      	forward.add(str[strand_index].charAt(0)=='+');
	      }
		}
	}
}
