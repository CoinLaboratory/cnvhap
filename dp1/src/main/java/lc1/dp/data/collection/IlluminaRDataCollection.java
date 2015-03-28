package lc1.dp.data.collection;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;
import lc1.dp.states.PhasedDataState;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.util.Constants;

public class IlluminaRDataCollection extends LikelihoodDataCollection{

	
  public HaplotypeEmissionState createEmissionState(String key, int no_copies){
      //if(stSp[1].size()==stSp1[1].size()) 
      if(index <0) throw new RuntimeException("index should be positive");
     
          return   new IlluminaNoBg(key, Emiss.getSpaceForNoCopies(no_copies),  this.snpid, Constants.cumulativeR(index),index);
     // return  new Illumina1NoBg(key, stSp[1],stSp1[1], trans[1], getR(key),getB(key), this.length, index) ;
     
  }
    
    @Override
  public boolean hasIntensity(int i) {
		return true;
	}
  
  
 
  
  @Override
    public Boolean process(String indiv, String[] header,  String[] geno, int i, int ploidy,  double[] missing){
        boolean doneB = false;
        boolean doneR = false;
        boolean nullB = true;
        boolean nullR = true;
   
        IlluminaNoBg state = (IlluminaNoBg)this.dataL.get(indiv);
     /*   if(header[0].equals("left")){
        	double r = ( Double.parseDouble(geno[0]) - Double.parseDouble(geno[1]))/Double.parseDouble(geno[0]);
      	  state.setR(i,r);
      	  doneR = true;
      	  return true;
      }*/
        
       for(int k=0; k<geno.length; k++){
           if(geno[k].startsWith("nul") || geno[k].equals("nan") )geno[k] = "NaN";
             String hk = header[k].toLowerCase();
             if(Constants.transformTheta() && (hk.indexOf("theta")>=0 || hk.indexOf("countb")>=0)){
            	  Double b = Double.parseDouble(geno[k]);
                  nullB  = (b==null || Double.isNaN(b));
                 ((IlluminaNoBg)state).setB(i,b);
            
                 doneB =true;
             }
           /*  else if(false & hk.indexOf("right")>=0){
            	  Double b = Double.parseDouble(geno[k]);
            	  Double r = Double.parseDouble(geno[k-1]);
            	  double v = Math.log(r/b);
            	  if(Double.isNaN(v)) v = 0;
            	  if(v<Constants.minR[index]){
                  	v = Constants.minR[index];
                   }
                   if(v>Constants.maxR[index]) {
                  	v =  Constants.maxR[index];
                   }
            	 
            	  ((IlluminaNoBg)state).setR(i,v);
            			  //(r-b)/Math.max(b,r));
            	  doneR = true;
            	  nullR = Double.isNaN(v);
             }*/
             else if(!doneB && (hk.indexOf("b allele")>=0   || hk.indexOf("b_allele")>=0)){
                 Double b = Double.parseDouble(geno[k]);
                 nullB  = (b==null || Double.isNaN(b));
                ((IlluminaNoBg)state).setB(i,b);
                doneB =true;
                 
            }
             else if(Constants.transformTheta() && (hk.equals("r") || hk.indexOf("counta")>=0)){
            	  double r =Math.log( Double.parseDouble(geno[k]));
            	 nullR  = Double.isNaN(r);
                doneR = true;
                 // if(this.snpid.get(i).startsWith("cnv1950")) r+=0.6;
                 if(!nullR){
                  if(r<Constants.minR[index]){
                 	r = Constants.minR[index];
                  }
                  if(r>Constants.maxR[index]) {
                 	r =  Constants.maxR[index];
                  }
                  state.setR(i,r);
                 
                 }
             }
             else if(
            		 !doneR && 
(            		 (header[k].indexOf("pc0")>=0 || header[k].indexOf("Log")>=0 || header[k].indexOf("log")>=0) 
		&& header[k].indexOf("Log RU")<0 || header[k].indexOf("Depth")>=0 || header[k].indexOf("depth")>=0)){
                 double r = Double.parseDouble(geno[k]);///avgDepth;
                 nullR  = Double.isNaN(r);
               
                // if(this.snpid.get(i).startsWith("cnv1950")) r+=0.6;
                if(!nullR){
                 if(r<Constants.minR[index]){
                	r = Constants.minR[index];
                 }
                 if(r>Constants.maxR[index]) {
                	r =  Constants.maxR[index];
                 }
                 state.setR(i,r);
                }
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
    	   ((IlluminaNoBg)state).emissions[i] = null;
    	   
       }
       else if(!doneB
    		   || snpid.get(i).toLowerCase().indexOf("cnv")>=0
    		   || snpid.get(i).toLowerCase().indexOf("A_")>=0
    		 || nullB
    		   ){
        Number r =  ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).r();
        ((IlluminaNoBg)state).emissions[i] = new IlluminaRDistribution(this.index);
       if(r!=null) ((IlluminaRDistribution)((IlluminaNoBg)state).emissions[i]).setR(r.doubleValue());
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
//    	        if(data.emissions[i]==null){
  //  	            EmissionState sta = this.dataL(indiv);
    //	             data.emissions[i] = new IntegerDistribution(sta.getBestIndex(i));
    	//         }
    	  //      data.emissions[i].setDataIndex(this.index);
    	       
    	    }catch(Exception exc){
    	        exc.printStackTrace();
    	    }
    	    Boolean res =  state.emissions[i]==null ? null : state.emissions[i].probeOnly();
    	 /*  if(res==null){
    		   int pos = this.loc.get(i);
    	   	Logger.global.info("all NC at "+i+" "+this.index+" "+pos);
    	   }*/
    	 if(nullB) missing[1]++;
    	 if(nullR) missing[0]++;
    	    return res;
    }
   /* @Override
    public void fix(Locreader loc, int thresh){
        List<Integer> fixed = new ArrayList<Integer>();
        List<Integer> nonFixed = new ArrayList<Integer>();
        EmissionStateSpace emstsp = this.getEmStSpace();
        double[] d = new double[emstsp.size()];
        for(int i=0; i< this.loc.size(); i++){
            int pos = this.loc.get(i).intValue();
            Location overl = loc.contains(pos,thresh);
           if(overl==null){
              fixed.add(this.loc.get(i));
            for(Iterator<EmissionState> it = this.dataLvalues(); it.hasNext();){
                HaplotypeEmissionState nxt = ((HaplotypeEmissionState)it.next());
                nxt.emissions[i].calcDistribution(this.getB(nxt.getName()), d, emstsp, this.no_copies);
                short di = nxt.emissions[i].getDataIndex();
                SimpleExtendedDistribution.normalise(d);
                int max = Constants.getMax(d);
                if(d[max]>0.999){
                    nxt.emissions[i] = new IntegerDistribution(max);
                }
                else{
                	nxt.emissions[i] = new SimpleExtendedDistribution(d, Double.POSITIVE_INFINITY);
                }
                nxt.emissions[i].setDataIndex(di);
            }
           }
           else{
               nonFixed.add(this.loc.get(i));
           }
        }
        this.calculateMaf(true);
        System.err.println("finished fixing \nfixed:"+fixed+"\n not fixed"+nonFixed);
    }*/
    
    public  IlluminaRDataCollection (File f, short index, 
    		int no_copies, int[][] mid,  File bf,Collection<String> snpidrest) throws Exception{
        super(f, index, no_copies, mid, bf,snpidrest);
        
    }
 
    public IlluminaRDataCollection(IlluminaRDataCollection collection) {
      super(collection);
      this.dc = collection.dc;
   //   this.distR = collection.distR;
    //  this.distB = collection.distB;
    }
    public IlluminaRDataCollection(DataCollection obj) {
        super(obj);
   //     this.make();
       // this.stSp = obj.stSp;
       // this.stSp1 = obj.stSp1;
       // this.trans = obj.trans;
        if(true) throw new RuntimeException("!!");
      
    }
public List<Integer> getBackgroundCN(){
	if(Constants.modelbg()){
	 throw new RuntimeException("!!");
	}
	else
    return Arrays.asList(new Integer[] {Constants.backgroundCount(this.index)});
}
@Override    
public void makeDistributions(int index) {
	
    	this.dc = new DistributionCollection(getBackgroundCN(),index, this.loc.size(),this.dir);
          }

    @Override
    public void initialisationStep(){
    //   if(true)  throw new RuntimeException("!!");
     //  this.distR.initialiseRCounts();
      // this.distB.initialiseBCounts();
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
        return new IlluminaRDataCollection(this);
    }
    public void printDist(File pw) {
   //     super.printDist(pw);
     // this.distR.print(pw);
    //  this.distB.print(pw);
    
        
    }
    public void print(String[] print, PrintWriter pw){
     for(int i=0; i<print.length; i++){
         pw.print(print[i]);
         pw.print(i<print.length-1  ? "\t": "\n");
     }
    }
    public void printHapMapFormat(File f, String chr) throws Exception{
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
              IlluminaNoBg ill = (IlluminaNoBg) this.dataL.get(list.get(i));
              toPrint[i*3+5] = comp.toStringPrint();
              IlluminaDistribution dist = (IlluminaDistribution) ill.emissions[i];
              toPrint[i*3+6]  = String.format("%5.3f", new Object[] {dist.b()});
              toPrint[i*3+7]  = String.format("%5.3f", new Object[] {dist.r()});
            }
            print(toPrint, pw);
//            pw.println(Format.sprintf(headerSt, toPrint));
        }
        pw.close();
    }


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
             double bval =  this.dc.sampleB(0,haploP, i);
             double rval = this.dc.sampleR(0,2,cn,i);
              ((IlluminaDistribution) hes.emissions[i]).setB(bval);
              ((IlluminaDistribution) hes.emissions[i]).setR(rval);
          }
      }
        
    }


	
}
