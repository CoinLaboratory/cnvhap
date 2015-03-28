package lc1.dp.data.collection;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.PhasedDataState;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRThetaDistribution;
import lc1.util.Constants;

public class ASCNDataCollection extends LikelihoodDataCollection{
     // IlluminaProbB distB;
   //   IlluminaProbR distR;
  
  //public IlluminaProbB getB(String key){
   //   return distB;
 // }
  
  //    public Object[] probDists() {
  	//	return new Object[] {new IlluminaProbR[] {this.distR}};
  	//}
//  public IlluminaProbR getR(String key){
  //    return distR;
 // }
  /*  @Override
 public double[]  getR(Location loc, List<Double> l, List<Integer > l1, List<Integer> l2, List<Double> b, String name){
      ASCNEmissionState emiss = (ASCNEmissionState) this.dataL.get(name);
      if(emiss==null) return null;
      EmissionStateSpace emstsp = emiss.getEmissionStateSpace();
      int start = DataCollection.firstGreaterThanOrEqual(this.loc, (int)loc.min);
      int end = IlluminaRDataCollection.firstGreaterThan(this.loc, (int)loc.max)-1;
      if(end-start >=1){
       //   System.err.println(this.loc.get(start)+"-"+this.loc.get(end) +" cf "+loc);
          double sum = 0;
          double geom = 0.0;
          double max =0;
          for(int i=start; i<=end; i++){
              l.add(emiss.r(i));
              b.add(emiss.b(i));
              if(loc.noCop()==0 && emiss.r(i) < max){
                  max = emiss.r(i);
              }
              else if(loc.noCop()==2 && emiss.r(i) > max){
                  max = emiss.r(i);
              }
              sum+=emiss.r(i);
              geom+=Math.log(emiss.r(i));
              l1.add(emstsp.getCN(emiss.mostLikely(i)));
              l2.add(this.loc.get(i));
          }
          geom = Math.exp(geom / (double)(end-start+1));
          sum = sum/(double)(end-start+1);
        //  System.err.println(l);
        //  System.err.println(l1);
        //  System.err.println(l2);
          Boolean val = null;
         return new double[] {sum, max};
      }
      return null;
  }*/
@Override
  public HaplotypeEmissionState createEmissionState(String key, int no_copies){
      //if(stSp[1].size()==stSp1[1].size()) 
      if(index <0) throw new RuntimeException("index should be positive");
      return   new IlluminaNoBg(key, null,  this.snpid, Constants.cumulativeR(index),index);
         
      //return   new IlluminaNoBg(key, stSp[no_copies-1],  this.snpid, index);
     // return  new Illumina1NoBg(key, stSp[1],stSp1[1], trans[1], getR(key),getB(key), this.length, index) ;
     
  }
    
    
  //public static double minR=0;
  //public static double maxR=0; 
@Override
public boolean hasIntensity(int i) {
	return true;
}
  
  @Override
    public final Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy, double[] missing){
      //  boolean doneB = false;
      //  boolean nullB = false;
        HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);
        IlluminaDistribution dist =  new IlluminaRThetaDistribution(this.index);
      state.emissions[i] =dist;
       for(int k=0; k<header.length; k++){
           if(geno[k].equals("null")) geno[k] = "NaN";
             String hk = header[k].toLowerCase();
             if(hk.indexOf("b")>=0){
            	 // double r = Double.parseDouble(geno[k]);
                 double b = Double.parseDouble(geno[k]);
                
                dist.setB(b);
                //dist.setR(r);
                /*if(r<Constants.minR){
                	Constants.minR = r;
                	Constants.minB = r;
                }
                if(r>Constants.maxR){
                	Constants.maxR = r;
                	Constants.maxB = r;
                }*/
                 
            }
             else if(hk.indexOf("a")>=0){
                 double r = Double.parseDouble(geno[k]);
                // if(this.snpid.get(i).startsWith("cnv1950")) r+=0.6;
                 dist.setR(r);
                 /*if(r<Constants.minR){
                	 Constants.minR = r;
                		Constants.minB = r;
                 }
                 if(r>Constants.maxR){
                	 Constants.maxR = r;
                		Constants.maxB = r;
                 }*/
                 
             }
            
       }
       try{
       ((IlluminaRThetaDistribution)state.emissions[i]).transformToThetaR();
       }catch(Exception exc){
    	   state.emissions[i] = 	Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(0.0);
       }
       try{
    	      //  boolean doneGeno = false;
    	        PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
    	        for(int k=0; k<header.length; k++){
    	      /*  if(header[k].toLowerCase().indexOf("geno")>=0){
    	            if(data.emissions[i]==null){
    	            int ind = trans(geno[k]);
    	        
    	                 data.emissions[i] = new IntegerDistribution(ind);
    	            }
    	        }*/
    	        
    	    }
    	       /* if(data.emissions[i]==null){
    	           // EmissionState sta = this.dataL(indiv);
    	             data.emissions[i] = new IntegerDistribution(0);
    	         }
    	        data.emissions[i].setDataIndex(this.index);*/
    	     //   if( Constants.format()[index].toLowerCase().equals("ascn") ){
   	        	
   	       // }
    	    }catch(Exception exc){
    	        exc.printStackTrace();
    	    }
    	    return state.emissions[i].probeOnly();
    }
  
    
    public  ASCNDataCollection (File f, short index, int no_copies, 
    		int[][] mid,  File bf,Collection<String> snpidrest) throws Exception{
        super(f, index, no_copies, mid, bf,snpidrest);
        
    }
 
    public ASCNDataCollection(ASCNDataCollection collection) {
      super(collection);
      this.dc = collection.dc;
     
    }
    public ASCNDataCollection(DataCollection obj) {
        super(obj);
   //     this.make();
       // this.stSp = obj.stSp;
        //this.stSp1 = obj.stSp1;
        //this.trans = obj.trans;
        if(true) throw new RuntimeException("!!");
       // this.data = obj.data;
//        CompoundEmissionStateSpace stSp = Emiss.getEmissionStateSpace(1, 2);
 //       CompoundEmissionStateSpace stSp1 = Emiss.getEmissionStateSpace(1)  ;
   //     EmissionStateSpaceTranslation trans = new EmissionStateSpaceTranslation(stSp, stSp1, true);
        makeDistributions(obj.index);
        
        StringBuffer sb = new StringBuffer();
        for(int i=0; i<loc.size(); i++){
            sb.append("%5.3f ");
        }
     
        for(Iterator<String> it = this.data.keySet().iterator(); it.hasNext();){
            String key = it.next();
            PIGData d = this.data.get(key);
            data.put(key, d);
            IlluminaNoBg ldc = 
               // stSp[1].size()==stSp1[1].size() ? 
                       // new ASCNEmissionState(d.getName(), stSp[1], this.snpid, index);
            	 new IlluminaNoBg(d.getName(), null,  this.snpid, Constants.cumulativeR(index), index);
              //  new Illumina1NoBg(d.getName(), stSp[1],stSp1[1], trans[1], getR(key),getB(key),this.loc.size(), index);
                 
          /*  for(int i=0; i<d.length(); i++){
                ComparableArray arr = (ComparableArray) d.getElement(i);
                int obj_index = stSp[1].getGenotype((arr));
            //    int ba = bg.sample(i);
               // System.err.println(ba);
//                        Constants.sample(bg.getEmiss(i)));
                int cn = stSp[1].getCN(obj_index);
                ldc.set(i, getR(key).sampleR(2,cn), getB(key).sampleB(obj_index));
            }*/
            this.dataL.put(key, ldc);
        }
        System.err.println("done");
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
    	this.dc = new DistributionCollection(this.getBackgroundCN(),index, this.loc.size(),this.dir);
        
    }
    public Integer getBGCount(int data_index, int pos){
    	return 1;
    	/*if(Constants.modelbg())
    	return bg.getFixedInteger(pos);
    	else return Constants.backgroundCount();*/
    	/*Integer val = bg.getFixedInteger(pos);
    	if(val!=null) return bg.getEmissionStateSpace().getCN(val);
    	else return null;
    	*/
    }

    @Override
    public void initialisationStep(){
     //  if(true)  throw new RuntimeException("!!");
     //  this.distR.initialiseRCounts();
     //  this.distB.initialiseBCounts();
    }
    @Override
    public void maximisationStep(double[] pseudo, int i){
        if(true) throw new RuntimeException("!!");
        super.maximisationStep(pseudo, i);
      /*  if( pseudo[3] < 1000){
           this.distR.rMaximistation(pseudo[3], i);
           this.distB.bMaximistation(pseudo[3], i);
        }*/
    }
    
    public DataC clone(){
        return new ASCNDataCollection(this);
    }
    public void printDist(PrintWriter pw) {
   //     super.printDist(pw);
   //   this.distR.print(pw);
    //  this.distB.print(pw);
    
        
    }
    public void print(String[] print, PrintWriter pw){
     for(int i=0; i<print.length; i++){
         pw.print(print[i]);
         pw.print(i<print.length-1  ? "\t": "\n");
     }
    }
   /* public void printHapMapFormat(File f, String chr) throws Exception{
        String len = "%7s";
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
        StringBuffer header = new StringBuffer("%-7s");
        List<String> list = new ArrayList(data.keySet());
        String[] toPrint = new String[5 + 3*list.size()];
        toPrint[0] = "Index";
        toPrint[1] = "Name";
        toPrint[2] = "Address";
        toPrint[3] = "Chr";
        toPrint[4] = "Position";
        for(int i=1; i<5; i++){
            header.append(" "+len);
        }
        for(int i=0; i<list.size(); i++){
            toPrint[i*3+5] = list.get(i)+".GType";
            toPrint[i*3+6] = list.get(i)+".B Allele Freq";
            toPrint[i*3+7] = list.get(i)+".Log R Ratio";
            header.append(" "+len);
        }
        String headerSt = header.toString();
       // int pos_index =s_i==null ? 0 : getFirstIndexAbove(s_i[0]);
       // pw.println(Format.sprintf(headerSt, toPrint));
       print(toPrint, pw);
        for( int pos_index =0; pos_index < loc.size(); pos_index++){
            toPrint[0] = "-"; toPrint[1]="-"; toPrint[2] =  "-";
            toPrint[3] = chr; toPrint[4] =this.loc.get(pos_index).toString(); 
            
            for(int i=0; i<list.size(); i++){
              String key = list.get(i);
              ComparableArray comp = (ComparableArray) this.data.get(key).getElement(pos_index);
              ASCNEmissionState ill = (ASCNEmissionState) this.dataL.get(list.get(i));
              toPrint[i*3+5] = comp.toStringPrint();
              ACSNDistribution dist = (ACSNDistribution) ill.emissions[i];
              toPrint[i*3+6]  = Format.sprintf("%5.3f", new Object[] {dist.dist1.r()});
              toPrint[i*3+7]  = Format.sprintf("%5.3f", new Object[] {dist.dist2.r()});
            }
            print(toPrint, pw);
//            pw.println(Format.sprintf(headerSt, toPrint));
        }
        pw.close();
    }*/


    public void simulate() {
      for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
          HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
          PhasedDataState dat =(PhasedDataState) this.data.get(hes.getName());
          EmissionStateSpace emstsp = dat.getEmissionStateSpace();
          for(int i=0; i<hes.noSnps(); i++){
            
              int intvalue = dat.emissions[i].fixedInteger();
              Comparable comp = emstsp.get(intvalue);
              int haploP = emstsp.getHaploPairFromHaplo(intvalue);
              int cn = emstsp.getCN(haploP);
          
             double rval = this.dc.sampleR(0,2,cn,i);
              //((IlluminaDistribution) hes.emissions[i]).setB(bval);
              ((IlluminaDistribution) hes.emissions[i]).setR(rval);
          }
      }
        
    }


	
}
