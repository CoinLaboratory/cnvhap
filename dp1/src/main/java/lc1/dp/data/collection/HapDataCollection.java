package lc1.dp.data.collection;

import java.io.File;
import java.util.Collection;
import java.util.logging.Logger;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.Constants;

public class HapDataCollection extends SimpleDataCollection {

	static final double lowerB = 1032;
	static final double  upperB = 1044;
	static String[] header = new String[] {"genotype"};
	public HapDataCollection(File file, short i2, int i, int[][] mid,
			File buildF,Collection<String> snpidrest) throws Exception{
		super(file,i2,i,mid, buildF, snpidrest);
	}
//	@Override
	 public Boolean  process(String indiv,  String[] header1,  String[] geno1, int i, int ploidy){
		
    	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	   String[] geno;
    //	  System.err.println("geno is "+Arrays.asList(geno1)+" hap is "+this.ch1);
    	   if(geno1[0].equals("null")){
    		   geno = new String[] {"null"};
    	   }
    	   else{
    		   geno = new String[] {trans(geno1[0])};
    	   }
    	try{
        	
        	
          //  boolean doneGeno = false;
        //    EmissionStateSpace stsp = this.stSp[this.no_copies-1];
        	   boolean deletionIsNa = Constants.useDeletion(i);
        	
           
           // PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
         
            for(int k=0; k<header.length; k++){
            if(header[k].toLowerCase().indexOf("geno")>=0 || header[k].toLowerCase().indexOf("plus_allele")>=0){
            	if(geno[k].equals("null")){
            	
                		dataSt.emissions[i] =stsp.getHWEDist1(null);
                	
            		return false;
            	}
                if(dataSt.emissions[i]==null){
                	
                	 int ind = 
                     	!this.abGenos ? 
                     	trans(geno[k], stsp, null, null,index) :
                     		trans(geno[k], stsp,alleleA.get(i), alleleB.get(i),index)
                     		;
                    if(geno[k].equals("null")){
                    	Logger.global.info("h");
                    }
                   double soften = Constants.softenHapMap(index);
                   // dataSt.emissions[i] = new IntegerDistribution(ind,stsp);
                    if(deletionIsNa & stsp.getCN(stsp.getHaploPairFromHaplo(ind))==0){
                    	if(Constants.format.length==1){
                    		dataSt.emissions[i] =stsp.getHWEDist1(null);
                    	}
                    	else{
                    		dataSt.emissions[i] = null;
                    	}
                    }
                    else{
                        if(soften>0){
                        	
                        	 dataSt.emissions[i] = stsp.getSoftenedDist(ind, 1-soften);
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
    }
	private String trans(String string) {
		char[] ch = new char[] {string.charAt(0)==ch1 ? 'B' :'A',
				string.charAt(1)==ch1 ? 'B' :'A'
				};
		return new String(ch);
	}
	public static char ch1 = Constants.hapAllele();
	private String getAllele(String st) {
		double i =Double.parseDouble(st);
		if(i<lowerB) {
			return "A";
		}
		else if(i<=upperB){
			return "B";
		}
		else {
			return "A";
		}
/*		if(i>=lowerB && i<=upperB){
			return "A";
		}
		else{
			return "_";
		}*/
	}
	
}
