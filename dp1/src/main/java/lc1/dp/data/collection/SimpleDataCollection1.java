package lc1.dp.data.collection;

import java.io.File;
import java.util.Collection;
import java.util.List;

import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PhasedDataState;

public class SimpleDataCollection1 extends SimpleDataCollection {
    public SimpleDataCollection1(){
        
    }
    protected SimpleDataCollection1(DataCollection dat){
      super(dat);
      
      
       
    }
    
  
   
    /*  private List<PIGData> data_bu;
    private List<PIGData> trio_bu;*/
    
    public SimpleDataCollection1(List<PIGData> name) {
       super(name);
    }
    
   
   
    
  /*  public void resetToTrio(){
        if(trio_bu!=null){
            this.data = trio_bu;
        }
    }*/
    
 
   public SimpleDataCollection1(File file,short index, 
		   int no_copies, int[][] mid,  File bf,Collection<String> snpidrest) throws Exception{
       super(file, index, no_copies, mid, bf,snpidrest);
       //if((Constants.soften(0))>0.01) modify();
      // this.cnv_region.clear();
    }
    
    public  SimpleDataCollection1 clone(){
        return new SimpleDataCollection1(this);
    }
    
   /* List<Boolean> cnv_region;
    double[] a_dist ;
    double[] non_cnv_region;
    double[] non_cnv_in_cnv_region;
    double[][] cnv; //goes from 0 to 4 copies
    static double[][] cnv_pqr = new double[][] {
    		new double[] {0.5, 0.5,0.0},
    		new double[] {0.4, 0.6,0.0},
    		new double[] {0.0,1.0,0.0},
    		new double[] {0,0.6,0.4},
    		new double[] {0,0.5,0.5}
    		
    };*/
    @Override
    public void initialise(){
    }
    
    /*public  double[] get(double[] pqr, double fracB, EmissionStateSpace stsp){
    	 double[] d =getDistOverCn(pqr);
    	 return HaplotypeEmissionState.fillCN(stsp, d, a_dist);
    }*/
   
    
   public static double[] getDistOverCn(double[] pqr){
    	double p = pqr[0];
    	double q = pqr[1];
    	double r = pqr[2];
    	return  new double[] {p*p, 2*p*q,q*q+2*p*r, 2*q*r, r*r}; //0,1,2,3,4 
    }
    
  /*  public void modify(){
    	   EmissionStateSpace stsp = Emiss.getSpaceForNoCopies(Constants.backgroundCount());
    	for(Iterator<String> it = this.dataL.keySet().iterator(); it.hasNext();){
    		String key = it.next();
    		HaplotypeEmissionState hes = (HaplotypeEmissionState) dataL.get(key);
    		for(int i=0; i<cnv_region.size(); i++){
    			if(cnv_region.get(i)){
    				int ind = hes.getFixedInteger(i);
    				if(ind!=this.no_copies){
    				//	if(ind==0 || ind==4){
    					//double[] dist_cnv = new double[5];
    					hes.emissions[i] = new SimpleExtendedDistribution(cnv[ind], Double.POSITIVE_INFINITY);
    				//	}
    					//else{
    						
    					//}
    				}
    				else{
    					//hes.emissions[i] = new SimpleExtendedDistribution(non_cnv_in_cnv_region, Double.POSITIVE_INFINITY);
    				}
    			}
    			else{
    				//hes.emissions[i] = new SimpleExtendedDistribution(non_cnv_region, Double.POSITIVE_INFINITY);
    			}
    		}
    	}
    }*/
   
   
    public Boolean  process(String indiv,  int ploidy, String[] header,  String[] geno, int i){
        try{
       
            PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
            HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
            EmissionStateSpace stsp = Emiss.getSpaceForNoCopies(dataSt.noCop());
            for(int k=0; k<header.length; k++){
            if(header[k].toLowerCase().indexOf("geno")>=0){
            	
              /*  if(data.emissions[i]==null){
                    int ind = trans(geno[k]);
                   int cn = stsp.getCN(stsp.getGenoForHaplopair(ind));
                    data.emissions[i] = new IntegerDistribution(ind);
                   if(cn!=this.no_copies) {
                    	if(cnv_region.size()<=i)cnv_region.add(true);
                    	else cnv_region.set(i, true);
                    }
                    else{
                    	if(cnv_region.size()<=i)cnv_region.add(false);
                    	//else cnv_region.set(i, false);
                    }
                   
                        
                         dataSt.emissions[i] = new IntegerDistribution(ind);
                   // }
                }*/
            }
            
        }
           
            data.emissions[i].setDataIndex(this.index);
        }catch(Exception exc){
            exc.printStackTrace();
        }
       return null;
    }
    
  
  
  /* if(cn!=this.no_copies){
                    	double[] dist = new double[stsp.size()];
                    	Arrays.fill(dist, 0.0);
                    	dist[ind] = 1-Constants.soften();
                    	String[] copies = cn<no_copies ? new String[] {"__","A_"} : 
                    		(Constants.modelCNP==6 ?new String[] {"AX","XX"} :  new String[] {"AX","XX","TX", "TT"});
                    	
                    	
                        double res = ((1.0-dist[ind])/(double) (copies.length));
                    	for(int kk=0; kk<copies.length; kk++){
                    		
                    		int pos = trans(copies[kk]);
//	                    	if(pos!=ind){
	                    		dist[pos] += res;
	//                    	}
	                    	
	                    //	throw new RuntimeException("!!");
	                      
                    	}
                    	double sum = Constants.sum(dist);
                    	if(Math.abs(sum-1.0)>0.001){
                    		throw new RuntimeException("!! "+sum);
                    	}
                    	 dataSt.emissions[i] = new SimpleExtendedDistribution(dist, Double.POSITIVE_INFINITY);
                    }
                    else{
   */
}
