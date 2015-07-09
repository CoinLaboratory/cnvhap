package lc1.dp.emissionspace;

import java.awt.Color;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;

import lc1.dp.data.representation.AbstractEmiss;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.IntegerEmiss;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.illumina.IlluminaProbR;
import lc1.stats.IntegerDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;


//Note default list behaviour is over haplo pairs - ie unorder lists of haplotypes!!!!
public abstract class  EmissionStateSpace implements List<Comparable>, Serializable {
    
	static Color getMix(double[] c0, double[] c1, double mix){
		return new Color((int) Math.floor (c0[0] *mix + c1[0]*(1-mix)),
						(int) Math.floor (c0[1] *mix + c1[1]*(1-mix)),
						(int) Math.floor (c0[2] *mix + c1[2]*(1-mix)));
	}
	
	 public static 	Color[] col;
	 static{
		int[] compArray = new int[4];
		float[] compArray1 = new float[4];
		 String[] str = Constants.color();
		 try{
			if(str.length==1 && str[0].startsWith("null") || str[0]==null){
				  double bg = Constants.backgroundCount(0);
				  double noStates = Constants.modify0[0].length;
				   double max = noStates * bg;
				  col = new Color[(int) (max+1)];
				
				 double[] c0 = new double[] {255,0,0};//r
				 double[] c1 = new double[] {0,0,255};//b
				 double[] c2 = new double[] {0,255,0};//g
			//	 double[] c3 = new double[] {125,0,255};
				col[0] = Color.PINK;
				 for(double k=1; k<bg; k++){
					col[(int)k] = getMix(c1,c0,(k-1)/(bg));
					 
				 }
				 double len = (double) col.length - bg;
				 double len_2 = Math.ceil(len/2.0)+bg;
				 for(double k=bg+1.0; k<len_2; k++){
						col[(int)k] = getMix(c2,c1,(k-(bg+1.0))/((double)len_2-1.0-(bg+1.0)));
						
						 
					 }
				 for(double k=len_2; k<col.length; k++){
					 
					 col[(int)k] = getMix(c0,c2,(k-(len_2))/((double)col.length-len_2));
				 }
				 col[(int)bg] = Color.GRAY;
				}
			else{
		 col= new Color[str.length];
		
		 for(int k=0; k<str.length; k++){
			 if(str[k].startsWith("#")){
				 
				 for(int jj =0; jj<compArray.length; jj+=2){
					
					 compArray[jj] =  Integer.parseInt(str[k].substring(1+2*jj, 1+2*jj+2),16) ;
				 }
				
				 col[k] = new Color(compArray[0], compArray[1],compArray[2],compArray[3]);
			 }else{
			 String[] str1 = str[k].split("~");
			 col[k] =(Color) Color.class.getField(str1[0]).get(null);
			 if(str1.length>1){
				 col[k] = new Color(col[k].getColorSpace(), col[k].getRGBComponents(compArray1), Float.parseFloat(str1[1]));
				// System.err.println(col[k]);
			 }
			 }
		 }
		 
			}
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
		// System.err.println("set cols ");
	 }
	 
	
	public String toString(){
    	return this.haploList.toString();
    }
  interface StringMethod extends Serializable {
        public abstract String getString(Comparable comp);
    }
  private PseudoDistribution[] intDist;
  public PseudoDistribution getIntDist(int ind) {
		
		if(intDist==null){
			intDist = new PseudoDistribution[this.size()]; 
		}
		if(intDist[ind]==null){
		 intDist[ind] = new IntegerDistribution(ind, this);
		}
		return intDist[ind];
	}
  
  Color[] color,color_muted ;
  
  public Color[] getColor(boolean muted){
	  if(color==null){
		  float[] compArray = new float[4];
		  color  = new Color[this.size()];
		  color_muted  = new Color[this.size()];
		 // SortedSet<Double> bafs = new TreeSet<Double>();
		  try{
		  for(int i=0; i<this.size(); i++){
			  //  shape[i] = IndividualPlot.shapes[emstsp.getBCount(i)];
			   // b_alias[i] = emstsp.getGenoForHaplopair(i);
			    color[i] = col[Math.min(col.length-1, getCN(i))];
			    color_muted[i] = 
			    	new Color(color[i].getColorSpace(), color[i].getRGBComponents(compArray),
			    			(float) Constants.muteAlpha());
				
			    		
		  }
		  }catch(ArrayIndexOutOfBoundsException exc){
			  System.err.println("too few colors "+col.length+" "+this.size());
			 exc.printStackTrace();
			  
			  System.exit(0);
		  }
	  }
	  return muted ? color_muted : color;
  }
 
  
  
  public int[] haploPairToGeno() {
		return this.haploPairToGeno;
	}
//public static EmissionStateSpace emStSp;;
    public  StringMethod gsm = new StringMethod(){

        public String getString(Comparable comp) {
         return getGenotypeString(comp);
        }
        
    };
    public  StringMethod hsm = new StringMethod(){

        public String getString(Comparable comp) {
         return getHaploPairString(comp);
        }
        
    };
   //  private int noCopies;
public String getHeader(){
	StringBuffer sb = new StringBuffer();
    for(int i=0; i<size(); i++){
        if(i>0) sb.append("\t");
        sb.append( getHaploPairString(get(i)));
       
    }
    return sb.toString();
}
     public List defaultList;

     
     public int getCNIndex(int i){
         for(int j=0;j<copyNumber.size() ; j++){
         if(this.copyNumber.get(j)==i) return j;
         }
         throw new RuntimeException("!!");
     }
     public List<Comparable> getGenotypes(){
    	 return genotypeList;
     }
//     public int noGenotypes;  //elements [0 ... noGenotypes] enumerate possible genotypes elements[noGenotypes to end] enumerate other haplolists
     private List<Comparable> genotypeList = new ArrayList<Comparable>();
     private List<Comparable> haplopairList = new ArrayList<Comparable>();
     protected List<Comparable> haploList = new ArrayList<Comparable>();
     public List<Integer> copyNumber = new ArrayList<Integer>();

     
     public int haploSize(){
    	 return this.haploList.size();
     }
     
    private  Map<Comparable, Integer> stateSpaceToCopyNumber = new TreeMap<Comparable, Integer>();//getGenoComparator());
    private  Map<Comparable, Integer> stateSpaceToGenotypeIndex = new TreeMap<Comparable, Integer>();//getGenoComparator());
    private  Map<Comparable, Integer> stateSpaceToHaplolistIndex = new HashMap<Comparable, Integer>();
    private  Map<Comparable, Integer> stateSpaceToHaploPairIndex = new TreeMap<Comparable, Integer>();//getHaploPairComparator());
public int numHaplotypes(){
    return this.haploList.size();
}

public Integer getHaploPairIndex(String st){
	return this.stateSpaceToHaploPairIndex.get(st);
}
public List<Comparable> getGenotypeList(){
    return genotypeList;
}
     double[][] genoToHaploW;
     Comparable[][] genoToHaplos;
     Comparable[][] haploPairToHaplos;

     private int[][] copyNoIndexToGeno;
     private int[][] genoToHaploPair; 
     private int[][] genoToHaplo;
   protected int[][] haploPairToHaplo; //haploPair is unordered, haplo is ordered
     protected int[] haploToHaploPair;
     private double[] haploPairIndexToWeight;
     private int[] haploPairToCN;
     private int[] haploPairToBCount;
     public int getBCount(int index){
    	 return haploPairToBCount[index];
     }
     
     public int getBCount(int index, int allele){
    	 ((AbstractEmiss)this.defaultList.get(index)).noB(allele);
    	 return haploPairToBCount[index];
     }
     public int getHaploPairFromHaplo(int hapl){
         return this.haploToHaploPair[hapl];
     }
     public int getCN(int index){
         return haploPairToCN[index];
     }
     
     private int[] haploPairToGeno;
    public int[] getGenoForCopyNo(int i){
        int j = this.copyNumber.indexOf(i);
        if(j<0){
        	return null;
//        	j = copyNumber.size()-1;
  //      	Logger.global.info("warning copy number "+i+" greater than max "+copyNumber.get(j));
        }
        return copyNoIndexToGeno[j];
    }
    public int[] getGenoForCopyNo(int[] is){
    	Set<Integer> res = new TreeSet<Integer>();
    	for(int i=0; i<is.length;i++){
        int j = this.copyNumber.indexOf(is[i]);
        if(j<0){
        	j = copyNumber.size()-1;
        	
        //	
        	Logger.global.info("warning2 copy number "+i+" greater than max "+copyNumber.get(j));
        }
        int[] gen = copyNoIndexToGeno[j];
        for(int k=0; k<gen.length; k++){
        	res.add(gen[k]);
        }
    	}
    	int[] res1 = new int[res.size()];
    	Iterator<Integer> it = res.iterator();
    	for(int k=0; k<res1.length; k++){
    		res1[k] = it.next();
    	}
        return res1;
    }
    
    public int[] getGenoExcludingCopyNo(int i){
    	List<Integer> l = new ArrayList<Integer>();
    	for(int i1=0; i1<this.cnLength(); i1++){
    		if(i1!=i){
    			 //int j = this.copyNumber.indexOf(i);
    	         int[]res = this.getGenoForCopyNo(i1);
    	         for(int k=0; k<res.length; k++){
    	        	 l.add(res[k]);
    	         }
    		}
    	}
    	int[] res = new int[l.size()];
    	for(int k=0; k<res.length; k++){
    		res[k] = l.get(k);
    	}
    	return res;
       
    }
    public int[] getHaplopairForGeno(int i){
        return genoToHaploPair[i];
    }

    public Comparable[] getHaploForHaploPair(Comparable comp){
        return haploPairToHaplos[this.getHaploPair(comp).intValue()];
    }
public Integer getHapl(String st){
    return this.stateSpaceToHaplolistIndex.get(st);
}
    public Integer getHapl(Comparable omp){
    	String str = getHaploString(omp);
       Integer res =  this.stateSpaceToHaplolistIndex.get(str);
       if(res ==null){
    	   String str1 = getHaploString(omp);
           throw new RuntimeException("nothing found for "+omp +" in \n"+this.stateSpaceToHaplolistIndex.keySet());
       }
       return res;
    }
    public Integer getHaploPair(Object omp){
        Integer res = this.stateSpaceToHaploPairIndex.get(this.getHaploPairString((Comparable)omp));
        return res;
       
     }
    public Comparable getHapl(Integer i){
        return ((ComparableArray)this.haploList.get(i)).copy();
    }
    public Comparable getHaploPair(Integer i){
        return this.haplopairList.get(i);
    }

//    following depend on what we decide is the state space
    /** what is the emission state index for comp */
    public Integer get(Object comp){
        return this.getHaploPair(comp);
    }
    /** what are the haplotype indices for the ith element in the state space */
    public int[] getHaps(int i){
        return this.haploPairToHaplo[i];
    }
  /*  double[][] calcArry = null;
    public double[] getCalcArray(int len) {
        return calcArray(len)
        // TODO Auto-generated method stub
        return null;
    }*/
    
    public boolean exclude(ComparableArray list) {
        if(list.get(0) instanceof Emiss && list.size()>1){
            if(Math.abs(((AbstractEmiss)list.get(0)).noCopies() -((AbstractEmiss)list.get(1)).noCopies()) >1) return true;
            
        }
        return false;
    }
    public Comparable getGenotype(int comp){
        return this.genotypeList.get(comp);
    }
    public Integer getGenotype(Object comp){
    	String str = getGenotypeString((Comparable)comp);
        return stateSpaceToGenotypeIndex.get(str);
    }
public int getRealCN(ComparableArray compA){
    return getCN(compA);
}
    public int getCN(ComparableArray compA) {
      return this.stateSpaceToCopyNumber.get(getGenotypeString(compA));
    }

    /** what is the i^th emission */
    public Comparable get(int i){
        return haploList.get(i);
    }
    public int[] getGenotypeConsistent(int i){
        return this.genoToHaploPair[i];
    }
    public double[] getWeights(int stSpIndex) {
     return this.genoToHaploW[stSpIndex];
    }
    protected Comparator<Comparable> getHaploPairComparator(){
        return new Comparator<Comparable>(){ //DEFAULT_COMPARATOR
            public int compare(Comparable o1, Comparable o2) {
            	String st1 = getHaploPairString(o1);
            	String st2 = getHaploPairString(o2);
            	int i1 = st1.indexOf(',');
            	int i2 = st2.indexOf(',');
            	if(i1==i2)
                return st1.compareTo(st2);
            	else return i1<i2 ? 1 :-1;
            }
            
        };
    }

    protected Comparator<Comparable> getGenoComparator(){
        return new Comparator<Comparable>(){ //DEFAULT_COMPARATOR
            public int compare(Comparable o1, Comparable o2) {
                String st2 = getGenotypeString(o2);
                String st1 = getGenotypeString(o1);
                return st1.compareTo(st2);
            }
            
        };
    }
    
    private  static void map1(int[][] copyNumberToGeno, List<Comparable> genotypeList, 
            Map<Comparable, Integer> stateSpaceToCNIndex, StringMethod sm, List<Integer> copyNumber
           ){
        for(int i=0; i<copyNumberToGeno.length; i++){
            int cn = copyNumber.get(i);
            List<Integer> hapl = new ArrayList<Integer>();
            for(int j=0; j<genotypeList.size(); j++){
                Comparable st = genotypeList.get(j);
                String str  = sm.getString(st);
                int ind = stateSpaceToCNIndex.get(
                        str);
                if(ind==cn){
                    hapl.add(j);
                }
              //  else{
              //      System.err.println("excl ");
              //  }
            }
            copyNumberToGeno[i] = new int[hapl.size()];
            for(int j=0; j<hapl.size(); j++){
                copyNumberToGeno[i][j] = hapl.get(j);
                //if(genoToHaplos!=null) genoToHaplos[i][j] =haploList.get(hapl.get(j));
            }
           
        }
        
    }
  
    /**fills in genoToHaplo
     * fills in genosToHaplos  (equivalent but with Comparable rather than index) can be null 
     *  */
    private static void map(int[][] genoToHaplo, Comparable[][] genoToHaplos, List<Comparable> haploList, 
            Map<Comparable, Integer> stateSpaceToGenotypeIndex, StringMethod sm
           ){
        for(int i=0; i<genoToHaplo.length; i++){
            
            List<Integer> hapl = new ArrayList<Integer>();
            for(int j=0; j<haploList.size(); j++){
                Comparable st = haploList.get(j);
                String str  = sm.getString(st);
                int ind = stateSpaceToGenotypeIndex.get(
                        str);
                if(ind==i){
                    hapl.add(j);
                }
             //   else{
              //      System.err.println("excl ");
             //   }
            }
            genoToHaplo[i] = new int[hapl.size()];
            if(genoToHaplos!=null) genoToHaplos[i] = new Comparable[hapl.size()];
            for(int j=0; j<hapl.size(); j++){
                genoToHaplo[i][j] = hapl.get(j);
                if(genoToHaplos!=null) genoToHaplos[i][j] =haploList.get(hapl.get(j));
            }
           
        }
        
    }
    public EmissionStateSpace(){
       // this.noCopies = noCopies;
    }
        public void init(List<Comparable> list){
            this.defaultList = this.haplopairList;
            Collections.sort(list, new Comparator<Comparable>(){

				@Override
				public int compare(Comparable arg0, Comparable arg1) {
					if(arg0 instanceof ComparableArray){
						if(arg1 instanceof ComparableArray){
							int n1 = ((ComparableArray)arg0).numLevels();
							int n2 =((ComparableArray)arg1).numLevels();
							if(n1==n2) return arg0.compareTo(arg1);
							else return n1<n2 ? -1 :1;
							
						}
						else return 1;
					}
					else if(arg1 instanceof ComparableArray) return -1;
					else return arg0.compareTo(arg1);
				}
            	
            });
            Set<Comparable> genotypeSet = new TreeSet<Comparable>(getGenoComparator());
            Set<Comparable> haplopairSet = new TreeSet<Comparable>(getHaploPairComparator());
            Set<Comparable> haplotypeSet = new HashSet<Comparable>();
            Set<Integer> cn = new TreeSet<Integer>();
            for(int i=0; i<list.size(); i++){
                Comparable el = list.get(i);
                genotypeSet.add(el);
                haplopairSet.add(el);
                haplotypeSet.add(el);
                cn.add(((AbstractEmiss)el).noCopies());
               
            }
            haplotypeSet.removeAll(haplopairSet);
            haplopairSet.removeAll(genotypeSet);  //so the sets are mutually exclusive
            this.copyNumber.addAll(cn);
            this.genotypeList.addAll(genotypeSet);
            this.haplopairList.addAll(genotypeSet);
            haploList.addAll(genotypeSet);
            this.haplopairList.addAll(haplopairSet);
            haploList.addAll(haplopairSet);
            
            haploList.addAll(haplotypeSet);
           
           // this.noCopies = noCopies;
            this.initialise();
            this.invCNVSize = 1.0/copyNumber.size();
            this.invSize = 1.0/defaultList.size();
            this.genoToHaploW = new double[this.genotypeList.size()][];
            this.haploPairIndexToWeight = new double[this.haplopairList.size()];
          //  Arrays.fill(genoToHaploW, 1.0);
         
            this.copyNoIndexToGeno = new int[this.copyNumber.size()][];
            this.genoToHaplo = new int[this.genoListSize()][];
            this.genoToHaploPair = new int[this.genoListSize()][];
            //this.genoToHaploW = new double[this.size()][];
            this.genoToHaplos = new Comparable[this.genoListSize()][];
            this.haploPairToHaplo = new int[this.haplopairList.size()][];
            this.haploToHaploPair= new int[this.haploList.size()];
            this.haploPairToHaplos = new Comparable[this.haplopairList.size()][];
            this.haploPairToGeno = new int[this.haplopairList.size()];
            this.haploPairToCN = new int[this.haploSize()];
            this.haploPairToBCount = new int[this.haploSize()];
            map(genoToHaplo, genoToHaplos, haploList, stateSpaceToGenotypeIndex, gsm);
            map1(copyNoIndexToGeno,  genotypeList,stateSpaceToCopyNumber, gsm, this.copyNumber);
         
        
            map(this.genoToHaploPair,null,  haplopairList, stateSpaceToGenotypeIndex, gsm);
            map(haploPairToHaplo,haploPairToHaplos, haploList, this.stateSpaceToHaploPairIndex, hsm);
            for(int i=0; i<this.haplopairList.size(); i++){
                Comparable hpli = this.haplopairList.get(i);
                String hplis = this.getGenotypeString(hpli);
                haploPairToGeno[i] = this.stateSpaceToGenotypeIndex.get( hplis);
                       
            }
            for(int i=0; i<this.haploPairToHaplo.length; i++){
                int[] hPToHi = this.haploPairToHaplo[i];
                for(int j=0; j<hPToHi.length; j++){
                    this.haploToHaploPair[hPToHi[j]] = i;
                }
            }
            //map(this.haploPairToGeno, null, this.genotypeList, this.stateSpaceToHaploPairIndex);
            
            for(int i=0; i<this.genoToHaploPair.length; i++){
                int[] haploPairIndices = genoToHaploPair[i];
               genoToHaploW[i] = new double[haploPairIndices.length];
                for(int j=0; j<haploPairIndices.length; j++){
                    Comparable comp1 =  haplopairList.get(haploPairIndices[j]);
                    if(comp1 instanceof Emiss || comp1 instanceof IntegerEmiss){
                        genoToHaploW[i][j] = 1.0;
                    }
                    else{
                        boolean exclude = exclude((ComparableArray)comp1);
                        if(exclude){
                            genoToHaploW[i][j]  = Constants.exclude();
                        }
                        else genoToHaploW[i][j]  = 1.0;
                      
                    }
                    this.haploPairIndexToWeight[haploPairIndices[j]] = genoToHaploW[i][j];
                }
                SimpleExtendedDistribution.normalise(genoToHaploW[i]);
               
            }
            for(int i=0; i<this.haploSize(); i++){
                Comparable comp = this.haploList.get(i);
                if(comp instanceof ComparableArray){
                    haploPairToCN[i] = ((ComparableArray)comp).noCopies(true);
                   this.haploPairToBCount[i] =(int)  ((ComparableArray)comp).noB();
                }
                else if(comp instanceof Emiss){
                    haploPairToCN[i] = ((Emiss)comp).noCopies();
                    this.haploPairToBCount[i] =(int)  ((Emiss)comp).noB();
                }
                else{
                	haploPairToCN[i] = 2;
                	  this.haploPairToBCount[i] = -100;
                }
              
            }
          
            if(true || Constants.annotate()){
                StringBuffer sw = new StringBuffer();
                sw.append("emission state space");
                
                for(int i=0; i<this.haploList.size(); i++){
                    Comparable comp = haploList.get(i);
                    sw.append(comp instanceof ComparableArray ? ((ComparableArray)comp).elements().toString() : comp.toString());
                    sw.append("  ");
                    if(i==this.genotypeList.size()-1 || i==this.haplopairList.size()-1) sw.append("\n");
                }
              for(int i=0; i<this.defaultList.size(); i++){
                  int[] hapL = this.getHaps(i);
                  for(int i1=1; i1<hapL.length; i1++){
                      if(this.haploList.get(i1).toString().equals(this.haploList.get(0).toString())) throw new RuntimeException("!!");
                  }
              }
                Logger.global.info(sw.toString());
            }
        }
       
      
        
      



        
      
        public void initialise(){
            for(int i=0; i<genotypeList.size(); i++){
                Comparable com = genotypeList.get(i);
                String genString = getGenotypeString(com);
                stateSpaceToGenotypeIndex.put( genString, i);
                if(com instanceof Emiss){
                    this.stateSpaceToCopyNumber.put(genString, ((Emiss)com).noCopies());
                }
                else if(com instanceof IntegerEmiss){
                    this.stateSpaceToCopyNumber.put(genString, 1);
                }
                else{
                    this.stateSpaceToCopyNumber.put(genString,((ComparableArray)com).noCopies(true));
                }
           //     this.stateSpaceToCopyNumber.put(genString, this.getCNIndex(com instanceof Emiss ?  : ));
            }
            for(int i=0; i<haploList.size(); i++){
                stateSpaceToHaplolistIndex.put(getHaploString(haploList.get(i)), i);
            }
            Comparator comp = this.getHaploPairComparator();
            for(int i=0; i<this.haploList.size(); i++){
                int j=0;
                Comparable o2 = haploList.get(i);
                for(;j<this.haplopairList.size(); j++){
                    if(comp.compare(haplopairList.get(j), o2)==0) break;
                }
                if(j==haplopairList.size()) throw new RuntimeException("!!");
                String haploPairString =  this.getHaploPairString(haploList.get(i));
                stateSpaceToHaploPairIndex.put(
                       haploPairString, j);
            }
          System.err.println(stateSpaceToGenotypeIndex.size()+" "+stateSpaceToHaplolistIndex.size());
        }
        
        
        
        //implementing list interface
        public int size() {
           return this.defaultList.size();
        }
        public boolean isEmpty() {
            return this.defaultList.isEmpty();
        }
        public boolean contains(Object o) {
          return defaultList.contains(o);
        }
        public Iterator<Comparable> iterator() {
           return defaultList.iterator();
        }
        public Object[] toArray() {
            return defaultList.toArray();
        }
        public <T> T[] toArray(T[] a) {
           return (T[])this.defaultList.toArray((T[]) a);
        }
        public boolean add(Comparable o) {
          return  defaultList.add(o);
        }
        public boolean remove(Object o) {
          return defaultList.remove(o);
        }
        public boolean containsAll(Collection<?> c) {
          return defaultList.containsAll(c);
        }
        public boolean addAll(Collection<? extends Comparable> c) {
           return defaultList.addAll(c);
        }
        public boolean addAll(int index, Collection<? extends Comparable> c) {
           return defaultList.addAll(index, c);
        }
        public boolean removeAll(Collection<?> c) {
            return defaultList.removeAll( c);
        }
        public boolean retainAll(Collection<?> c) {
            return defaultList.retainAll( c);
        }
        public void clear() {
           defaultList.clear();
            
        }
        public Comparable set(int index, Comparable element) {
          return (Comparable) defaultList.set(index, element);
        }
        public void add(int index, Comparable element) {
            defaultList.add(index, element);
            
        }
        public Comparable remove(int index) {
            return (Comparable)  defaultList.remove(index);
        }
        public int indexOf(Object o) {
          return defaultList.indexOf(o);
        }
        public int lastIndexOf(Object o) {
            return defaultList.lastIndexOf(o);
        }
        public ListIterator<Comparable> listIterator() {
            return defaultList.listIterator();
        }
        public ListIterator<Comparable> listIterator(int index) {
            return defaultList.listIterator(index);
        }
        public List<Comparable> subList(int fromIndex, int toIndex) {
            return defaultList.subList(fromIndex, toIndex);
        }
        public int haplopairListSize() {
            return haplopairList.size();
        }
       /* public int noCopies() {
           return noCopies;
        }*/
        public Integer genoListSize() {
           return this.genotypeList.size();
        }
        
        
        public abstract String getHaploPairString(Comparable comp);
         
        public abstract String getGenotypeString(Comparable comp);
        public abstract String getHaploString(Comparable comp);
          
       
        public int getGenoForHaplopair(int obj_index) {
           return haploPairToGeno[obj_index];
        }
        
        public Integer getFromString(String sti) {
         // System.err.println(sti);
        	Integer res =  this.stateSpaceToGenotypeIndex.get(sti);
        // if(res==null){
          //   throw new RuntimeException("  !!"+sti);
         //}
          return res;
        }
        public double[] getArray(String c) {
          double[] res = new double[this.defaultList.size()];
          Arrays.fill(res, 0.0);
          if(c.equals("e") || c.equals("a")){
              Arrays.fill(res, 1.0/(double)res.length);
              return res;
          }
          else if(c.equals("f")){
              char[] ch = new char[] {'A', '_'};
              for(int j=0; j<ch.length; j++){
                  res[this.getFromString(ch[j]+"")] = 1.0 / (double)ch.length;
              }
              return res;
          }
          else{
             Integer i =  this.getFromString(c);
             if(i!=null){
                 res[i] = 1.0;
                 return res;
             }
             else{
            	 try{
            		 int nocop = Integer.parseInt(c);
            		 int[] ind = this.getGenoForCopyNo(nocop);
            		 for(int k=0; k<ind.length; k++){
            			 int[] haploPair = this.getHaplopairForGeno(ind[k]);
            			 for(int k1 = 0; k1<haploPair.length; k1++){
            				 res[haploPair[k1]] = 1.0;
            			 }
            		 }
            		
            	 }catch(Exception exc){
            		exc.printStackTrace(); 
            	 }
            	 Constants.normalise(res);
                return res;
             }
              /*EmissionStateSpace[] membs = ( (CompoundEmissionStateSpace)this).getMembers();
              for(int j=0; j<membs.length; j++){
                  Integer i1  = membs[j].getFromString(c);
                  if(i1!=null){
                    
                      for(int ik=0; ik<this.size(); ik++){
                          if(( (CompoundEmissionStateSpace)this).getMemberIndices(ik)[j]==i1.intValue()){
                              res[ik] = 1.0;
                          }
                      }
                      SimpleExtendedDistribution.normalise(res);
                  return res;
                  }
              }
              throw new RuntimeException("!!not found");*/
             //}
          }
        }
        private String str(char c) {
			if(c=='_') return "";
			else return c+"";
		}
		public double getWeight(int haploPairIndex) {
           return haploPairIndexToWeight[haploPairIndex];
        }
       int[] swtch = null;
		public int[] getSwitchTranslation() {
			if(swtch==null){
          swtch = new int[this.haplopairList.size()];
            for(int i=0; i<swtch.length; i++){
            	Comparable comp = this.haplopairList.get(i);
            	Comparable comp1 = SimpleScorableObject.switchAlleles(comp);
               swtch[i] = this.getHaploPair( comp1);
           }
			}
            return swtch;
        }
        public double[] mix(double[] probs, double[] array, double[] ds) {
          double[] res = new double[probs.length];
          for(int i=0; i<res.length; i++){
              res[i] = probs[i] *ds[0]+array[i]*ds[1];
          }
          return res;
        }
        public int[] getHaploFromHaploPair(int hapl) {
          return this.haploPairToHaplo[hapl];
        }
        public int flip(int hap1_j) {
           int[] res = this.getHaploFromHaploPair(this.getHaploPairFromHaplo(hap1_j));
           if(res.length==1) return res[0];
           if(res.length>2) throw new RuntimeException("!!");
           if(res[1]==hap1_j) return res[0];
           else return res[1];
        }
        double invCNVSize;
        double invSize;
        public double invCNVSize() {
            // TODO Auto-generated method stub
            return invCNVSize;
        }
        public double bSpaceSize(int j){
            return this.getGenoForCopyNo(this.getCN(j)).length;
        }
        public double invSize() {
           return invSize;
        }
        int[] noBIndex;
		public int[] noBIndex() {
			if(noBIndex==null){
				List<Integer> nob = new ArrayList<Integer>();
				for(int i=0; i<this.defaultList.size(); i++){
					if(((ComparableArray)defaultList.get(i)).noB()==0){
						nob.add(i);
					}
				}
				noBIndex = new int[nob.size()];
				for(int i=0; i<nob.size(); i++){
					noBIndex[i] = nob.get(i);
				}
			}
			return noBIndex;
		}
		public int cnLength() {
			return this.copyNumber.size();
		}
		public int[][] aliasNB;
		public int[][] aliasNB() {
			if(aliasNB==null){
				int cns = this.copyNumber.get(copyNumber.size()-1)+1;
			aliasNB = new int[cns][cns+1];
			for(int i=0; i<aliasNB.length; i++){
				int[] forcn = this.getGenoForCopyNo(i);
			//	aliasNB[i] = new int[i+1];
				Arrays.fill(aliasNB[i],-1);
				if(forcn!=null){
				for(int j=0; j<forcn.length; j++){
					int bc =  this.getBCount(forcn[j]);
					aliasNB[i][bc] = forcn[j];
					
				}
				}
			}
			
			}
			return aliasNB;
		}
		
		public int getByAlias(int noCop, int cntB) {
			if(aliasNB==null){
				aliasNB();
			}
			return aliasNB[noCop][cntB];
		}
		public boolean allOneLength() {
			 boolean allOneLength = true;
		      for(int i=0; i< this.haploPairToHaplo.length; i++){
		          int[] equiv = haploPairToHaplo[i];
		          if(equiv.length>1) allOneLength = false;
		      }
		      return allOneLength;
		}



		
      
    

    
    
}
