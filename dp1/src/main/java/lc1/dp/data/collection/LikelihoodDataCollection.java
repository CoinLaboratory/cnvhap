package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import lc1.CGH.Location;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.stats.IntegerDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import pal.statistics.NormalDistribution;

public class LikelihoodDataCollection extends DataCollection{
    
    public static String[] cat =  new String[] {
        "",
        "A",
        "B",
        "AA",
        "AB",
        "BB",
        "AAA",
        "AAB",
        "ABB",
        "BBB",
        "AAAA",
        "AAAB",
        "AABB",
        "ABBB",
        "BBBB"
      } ;
    
    public boolean hasIntensity(int i) {
    	return false;
    }
  
    public static DataCollection read(File f){
        try{
         return new LikelihoodDataCollection(cat);
        }catch(Exception exc){
            exc.printStackTrace();
            System.exit(0);
            return null;
        }
    }
    
    public double[]  getR(Location loc, List<Double> l, List<Integer > l1, List<Integer> l2, List<Double> b, String name){
    throw new RuntimeException("!!");
    }
  
  //  int noIndiv;
    
    public DataC clone(){
        return new LikelihoodDataCollection(this);
    }
  
   
    protected LikelihoodDataCollection(DataCollection dat){
        super(dat);
       // this.dataL = new HashMap<String, EmissionState>(dat.dataL);
    }
    
    public void maximisationStep(double[] pseudo, int i){
        if(Constants.trainData()){
            for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
               EmissionState dataLi =  it.next();
               dataLi.transferCountsToProbs(pseudo[4]);
            }
        }
    }
   
    
    static boolean readLine(BufferedReader[] br, List<String>[] str) throws Exception{
        int lengthRest = Constants.restrict()[0];
        for(int i=0; i<str.length; i++){
            if(br[i]!=null){
                String stri = br[i].readLine();
                if(stri==null) return false;
                else{
                    List<String> l  = Arrays.asList(stri.split("\\s+"));
                    int min = Math.min(l.size(), lengthRest+1);
                    str[i] = l.subList(1, min);
                }
            }
        }
        return true;
    }
    
    static List<String> readNames(File id_file){
        if(!id_file.exists() || id_file.length()==0) return null;
        List<String > l = new ArrayList<String>();
        try{
          if(!id_file.exists()|| id_file.length()==0) throw new RuntimeException("!!!");
           BufferedReader br = new BufferedReader(new FileReader(id_file));
           String st = "";
           while((st = br.readLine())!=null){
              l.add(st);
           }
         
        }catch(Exception exc){
            exc.printStackTrace();
        }
        return l;
    }
    
    public void extractFromTrioData(){
        super.extractFromTrioData();
        Map<String, EmissionState> l = new HashMap<String, EmissionState>();
        for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
            
            EmissionState dat_i = it.next();
            if(dat_i instanceof CompoundState){
                EmissionState[] data_i = dat_i.split();
                for(int j=0; j<data_i.length; j++){
                    l.put(data_i[j].getName(), data_i[j]);
                  
                }
            }
            else{
                l.put(dat_i.getName(), dat_i);
              
            }
        }
        this.dataL = l;
      
    }
    
    public HaplotypeEmissionState createEmissionState(String key, int no_copies){
        //if(stSp[1].size()==stSp1[1].size()) 
    
            return   SimpleScorableObject.make(key, loc.size(), null, this.index);
       // return  new Illumina1NoBg(key, stSp[1],stSp1[1], trans[1], getR(key),getB(key), this.length, index) ;
       
    }
    
    @Override
    public  void createDataStructure(List<String> indiv,List<Integer> ploidy, List<Integer> sampIdToIncl){
    	Iterator<Integer> it = 
    		sampIdToIncl==null ? Constants.iterator(indiv.size()) :
    		sampIdToIncl.iterator();
    	
        while(it.hasNext()){
        	int i =it.next();
            String key = indiv.get(i);
            int pl =ploidy.get(i);
              HaplotypeEmissionState  value = createEmissionState(key, pl);
           
            value.setNoCop(pl);
          if(dataL.containsKey(key)) Logger.global.warning("repeat samples found "+key);
          else dataL.put(key, value);
          //  data.put(key, SimpleScorableObject.make(key,loc.size(), value.getEmissionStateSpace(),this.index));
          
        }
    }
    
   @Override
    public   Boolean   process(String indiv,  String[] header,  String[] geno, int i, int ploidy, 
    		double[]missing){
        try{
          //  boolean doneGeno = false;
        	if(true) throw new RuntimeException("!!");
            PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
            HaplotypeEmissionState dataL = (HaplotypeEmissionState) this.dataL.get(indiv);
           dataL.emissions[i] = null;
            CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(dataL.noCop());
            for(int k=0; k<header.length ; k++){
            	process(data,dataL,  header[k], k<geno.length  ? geno[k] : null, i, missing);
            	
            }
            if(dataL.emissions[i]==null){
            	boolean norm = false;
            	boolean na = false;
        	inner:	for(int k=0; k<header.length; k++){
        			if(header[k].toLowerCase().startsWith("geno")){
        				//System.err.println(geno[k]);
        				if( geno[k].indexOf('N')<0){
        					 int ind = 
        	                    	!this.abGenos  || alleleA.size()==0 || alleleB.size()==0? 
        	                    	SimpleDataCollection.trans(geno[k], stsp, null, null,index) :
        	                    		SimpleDataCollection.trans(geno[k], stsp,alleleA.get(i), alleleB.get(i),index)
        	                    		;
        				 dataL.emissions[i] = new IntegerDistribution(
   								ind,
   								stsp);
        				}
        				else{
        					 dataL.emissions[i] = null; 
        					 
        				}
        				break inner;
        			}
        			/*else if(header[k].equals("Log R")){
        				  double d =  Double.parseDouble(geno[k]);
        				  if(d<-0.1 || d> 0.1){
        					  na = true;
        				  }
        			}*/
        			else if(header[k].equals("B allele")){// || header[k].equals("Log R")){
        				norm = true;
        				boolean b = header[k].equals("B allele");
   					  double d =  Double.parseDouble(geno[k]);
   					  if(dataL.emissions[i]==null){
   						double[] pr = new double[stsp.defaultList.size()];
   						Arrays.fill(pr,1);
						dataL.emissions[i] = new SimpleExtendedDistribution(pr, Double.POSITIVE_INFINITY, stsp);
   					  }
   					  double[] pr = dataL.emissions[i].probs();
   						  
   						  int bg = Constants.backgroundCount(index);
   						
   						inner1: for(int k1=0; k1<pr.length; k1++){
   							int cn = stsp.getCN(k1);
   							if(cn!=bg) pr[k1] =0;
   							else{
   							if(b && k1==0) continue inner1;
   							double mean = header[k].equals("B allele") ? 
   									(double)stsp.getBCount(k1)/(double)cn:Double.NaN;
   	   							//Constants.r_mean[0][stsp.getCN(k1)];
						
   						  pr[k1] = NormalDistribution.pdf(d,  mean, 0.05);
   							}
   						
   					  }
   				   }
   				   
   				 //  return false;
   			   }
            	if(na){
            		dataL.emissions[i] = null;
            	}
            	else if(norm  ) {
            		Constants.normalise(dataL.emissions[i].probs());
            	//	System.err.println(stsp.get(Constants.getMax(dataL.emissions[i].probs()))+" "+Arrays.asList(geno));
            	}
            	
   			  
        	}
            if(dataL.emissions[i]==null){
 				  String p_id =  this.snpid.get(i).toLowerCase();
 				  boolean probeOnly = p_id.startsWith("cnv") || p_id.startsWith("a");
 				   dataL.emissions[i] = ((CompoundEmissionStateSpace)stsp).getHWEDist1(probeOnly ? (new Double(0)): null);
 			   }
          if(Constants.likelihoodInput() ){
        	  ((SimpleExtendedDistribution)dataL.emissions[i]).normalise();
          }
          return dataL.emissions[i].probeOnly();
        }catch(Exception exc){
            exc.printStackTrace();
        }
        return null;
       
    }
    
    @Override
    public boolean process(HaplotypeEmissionState data, HaplotypeEmissionState dataL, String header, String geno, int i, double[] missing){
	   if(geno!=null && !geno.equals("null") && !super.process(data,dataL,  header,geno,i, missing) && !header.toLowerCase().startsWith("geno") ){
			PseudoDistribution dist = dataL.emissions[i];
			CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(dataL.noCop());
			
			if(header.equals("distr")){
				if(geno==null){
				if(dist==null){
					double[] pr = new double[stsp.defaultList.size()];
					dist = new SimpleExtendedDistribution(pr, Double.POSITIVE_INFINITY, stsp);
					dataL.emissions[i] = dist;
				}
				   double[] probs = dataL.emissions[i].probs();
				Arrays.fill(probs, 1.0/ (double) probs.length);
				}
			  // double[] probs = dataL.emissions[i].probs();
				
				String[] str = geno.split("=");
				String sta;
			//	double Double.parseDouble(str[1]) = ;
				if(str.length==1){
					dist = new IntegerDistribution(stsp.getFromString(str[0].replaceAll("_","")), stsp);
					dataL.emissions[i] = dist;
				}
				else if(str.length==2 && Double.parseDouble(str[1].split(";")[0])>=0.99 ){
					String st1 = str[0];
					st1 = st1.replace(';',',');
					String[] str_ = st1.split(":");
					sta = str_[0];
					
					dist = new IntegerDistribution(stsp.getHapl(sta), stsp);
					dataL.emissions[i] = dist;
				}
				else{
					double[] pr = new double[stsp.defaultList.size()];
					Arrays.fill(pr, 0);
					
					dist = new SimpleExtendedDistribution(pr, Double.POSITIVE_INFINITY, stsp);
					dataL.emissions[i] = dist;
					double[] probs = dist.probs();
				//if(str.length>2){
				//	Logger.global.info("h");
				//}
				for(int k=0; k<str.length-1; k++){
					
					if(k==0){
						String st1 = str[k];
						st1 = st1.replace(';',',');
						String[] str_ = st1.split(":");
						sta = str_[0];
					}
					else{
						String st1 = k==0 ? str[k] : str[k].replaceFirst(";", ":");
						st1 = st1.replace(';',',');
						String[] str_ = st1.split(":");
						sta = str_[1];
					}
					
					
					
					 Integer index = stsp.getHapl(sta);
					 if(index==null){
						 throw new RuntimeException("prob with "+sta+"\n"+geno);
					 }
					 Integer ind1 = stsp.getHaploPairFromHaplo(index);
					
					 probs[ind1] = Double.parseDouble(str[k+1].split(";")[0]);
				}
				double sum = Constants.sum(probs);
				if(sum<1.0){
					double diff = 1-sum;
					double[] p1 = stsp.getHWEDist1(null).probs();
					for(int i1=0; i1<probs.length; i1++){
						probs[i1]+=p1[i1]*diff;
					}
				}
				//Constants.normalise1(probs);
				}
			}
			else{
				header = header.replace("X", "AA").replace("Y", "AB").replace("Z", "BB").replace("0", "").replace("1", "A").replace("2", "AA");
			   Integer index =stsp.getHapl(header);

	    
	    	
	        if(index!=null){
				   if(dist==null){
						double[] pr = new double[stsp.defaultList.size()];
						dist = new SimpleExtendedDistribution(pr, Double.POSITIVE_INFINITY, stsp);
						dataL.emissions[i] = dist;
					}
				   double[] probs = dataL.emissions[i].probs();
	             int genot = stsp.getHaploPairFromHaplo(index);
	             //System.err.println(header[k]+" "+emstsp.get(genot));
	             double	v = Double.parseDouble(geno);
	           //  int cn = stsp.getCN(genot);
	            // if(cn>2 && v>0){
	            //	 Logger.global.info("h");
	            // }
	            probs[genot] =  v;
			}
			}
         }
	  // if(Constants.su)
	   return true;
   }
    
    
   /* @Override
    public  void process(String indiv, String[] header,  String[] geno, int i){
        HaplotypeEmissionState state = (HaplotypeEmissionState)this.dataL.get(indiv);
        EmissionStateSpace emstsp=  state.getEmissionStateSpace();
        if(state.emissions[i]==null){
            state.emissions[i] = new SimpleExtendedDistribution(state.getEmissionStateSpace().size());
        }
        SimpleExtendedDistribution dst = (SimpleExtendedDistribution) state.emissions[i];
        double[] probs = dst.probs;
        
       for(int k=0; k<header.length; k++){
           if(header[k].indexOf("Genot")<0){
              
           }
           else super.process(indiv, header, geno,i);
       }
       if(Constants.CHECK){
    	   double sum = Constants.sum(probs);
    	   if(Math.abs(1.0-sum)>0.01) throw new RuntimeException(" "+sum);
       }
      
        //dst.normalise(probs);
     
    }*/
    
   
  //  public LikelihoodDataCollection
   

    
   
  public LikelihoodDataCollection(){
        
    }
  
 // EmissionStateSpace[] emStSp;
 
 
   public LikelihoodDataCollection(String[] categories) throws Exception{
       if(Constants.onlyCopyNo()){
           for(int i=0; i<categories.length; i++){
               char[] ch = categories[i].toCharArray();
               Arrays.fill(ch, 'A');
               categories[i] = new String(ch);
           }
       }
       File dir1 = new File(Constants.baseDir());
       File fi = new File(dir1,"L.50K_Merged.CEPH.Affy.txt");
       String[] indiv = readIndiv(new File(dir1, "indiv.txt"));
       BufferedReader br = new BufferedReader(new FileReader(fi));
       String st = br.readLine();
       boolean header = st.startsWith("chrom");
       if(header){
           st = br.readLine();
       }
       this.loc = readPosInfo(fi,1, header);
       int read = IlluminaRDataCollection.firstGreaterThan(loc, Constants.offset());
       loc = loc.subList(read, loc.size());
       for(int i=0; i<read; i++){
           st = br.readLine();
       }
       CompoundEmissionStateSpace stSp = Emiss.getEmissionStateSpace(1);
       EmissionStateSpace stSp1 = Emiss.getEmissionStateSpace(1);
       EmissionStateSpaceTranslation trans = new EmissionStateSpaceTranslation(stSp, stSp1, true);
       this.length = Math.min(Constants.restrict()[0], loc.size());
         int[] categoryToGenotypeIndex = new int[categories.length];
         Arrays.fill(categoryToGenotypeIndex, -1);
         outer: for(int i=0; i<categories.length; i++){
             String str = categories[i];
             for(int k=0; k<stSp.genoListSize(); k++){
                 String genString = stSp1.getGenotypeString(stSp.getGenotype(k));
                 if(str.equals(genString)){
                     categoryToGenotypeIndex[i]=k;
                     continue outer;
                 }
             }
         }
        
        EmissionState[] ldata=null;
        String[][] str = null;
     
        for(int k=0; k<length; k++){
            //try{
                boolean NA = st.indexOf("NA")>=0;
            List<String> string  = Arrays.asList(st.split("\\s+"));
             if(k==0){
                 if(Math.abs(Math.IEEEremainder(string.size()-2, categories.length)) > 0.001) throw new RuntimeException("!!");
                 double len = ((double)string.size()-2.0)/(double)categories.length;
                
                 ldata = new EmissionState[(int)len ];
                if(ldata.length!=indiv.length) throw new RuntimeException("!!");
                 for(int i=0; i<ldata.length; i++){
                     ldata[i] = EmissionState.getEmissionState(indiv[i], stSp1, length);
                 }
                str = new String[ldata.length][];
             }
           
             for(int i=0; i<ldata.length; i++){
                 str[i] = string.subList(2+i*categories.length, 2+(i+1)*categories.length).toArray(new String[0]);
             }
             for(int i=0; i<ldata.length; i++){
                 double[] d = new double[str[i].length];
                 for(int k1=0; k1<d.length; k1++){
                     d[k1] = NA ? -1 :  -1*Double.parseDouble(str[i][k1]);
                 }
                // ((AffyEmissionState)ldata[i]).addDataPoint(d, categoryToGenotypeIndex,null, k, trans, stSp);
             }
             st = br.readLine();
            
         }
         for(int i=0; i<ldata.length; i++){
           //  ldata[i].updateBestIndex();
             dataL.put(ldata[i].getName(), ldata[i]);
             }
      this.calculateMLGenotypeData(true);
    //  initialiseMaf();
       br.close();
 
    }
   
  
public LikelihoodDataCollection(List<Integer> locs) {
    this.length = locs.size();
    this.loc = locs;
}



public LikelihoodDataCollection(File f, short index, int no_copies,int[][] mid,  File bf,
		Collection<String> snpidrest) throws Exception {
    super(f, index, no_copies, mid,  bf, snpidrest);
}




   
   /*public void initialiseMaf() {
       CompoundEmissionStateSpace stSp = (CompoundEmissionStateSpace) Emiss.getEmissionStateSpace(1);
       EmissionStateSpace stSp1 = stSp.getMembers()[0];
       this.maf = makeMafState(stSp1);
     for(int i=0; i<length; i++){
         for(Iterator<EmissionState> it =this.dataL.values().iterator(); it.hasNext();){
             EmissionState nxt = it.next();
            double[] probs = nxt.getEmiss(i);
            for(int k=0; k<probs.length; k++){
                int[] indices = stSp.getMemberIndices(k);
                for(int k1=0; k1 < indices.length; k1++){
                    maf.addCount(indices[k1], probs[k], i);
                }
            }
         }
     }
     maf.transferCountsToProbs(0);
     maf.initialiseCounts();
    
}*/



/*public LikelihoodDataCollection(List<EmissionState> name) {
       EmissionStateSpace stSp = Emiss.getEmissionStateSpace(1);
       this.length = name.get(0).noSnps();
       for(Iterator<EmissionState> it = name.iterator(); it.hasNext();){
           EmissionState nxt = it.next();
           this.dataL.put(nxt.getName(), nxt);
           PhasedIntegerGenotypeData res  = nxt.getGenotypeData(stSp);
           this.data.put(res.getName(), res);
       }
    // TODO Auto-generated constructor stub
}*/



public void print(PrintWriter pw){
       List<Integer> posi = new ArrayList<Integer>();
       StringBuffer sb2= new StringBuffer();
    
       for(int i=0; i<this.length(); i++){
           posi.add(i);
           sb2.append("%8i ");
       }
       pw.println(" "+String.format(sb2.toString(), posi.toArray()));
       for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
           
           it.next().print(pw, "", null);
           pw.println();
       }
    }














}
