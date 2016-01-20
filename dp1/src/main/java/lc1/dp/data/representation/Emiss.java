package lc1.dp.data.representation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.CompoundEmissionStateSpace2;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.SimpleEmissionStateSpace;
import lc1.util.Constants;

public class Emiss implements AbstractEmiss{
  // final int i;
	
	
	public static  EmissionStateSpace[] spaceByCN ;
    public static  EmissionStateSpace mergedSpace;
    public static  CompoundEmissionStateSpace getEmissionStateSpace(int i){
        return emiss.spaceByPloidy[i];
      }
    
    public static void switchHWE() {
    	for(int k=0; k<emiss.spaceByPloidy.length; k++){
    if(emiss.spaceByPloidy[k]!=null)emiss.spaceByPloidy[k].recalcHWEDist();
    	}
		// TODO Auto-generated method stub
		
	}
    
    public static CompoundEmissionStateSpace getSpaceForNoCopies(int no_copies){
      	return emiss.getSpaceForNoCopies(no_copies);
      }
    public static Emiss switchElement(Emiss comp) {
      return emiss.switchElement(comp);
      }
	public static Convert conv = new Convert();
	public static CompoundEmissionStateSpace getEmissionStateSpace(
			int[] numCopies) {
		return emiss.getEmissionStateSpace(numCopies);
	}
	
	public static int[] alleles;
	
	public static Emiss[] ems, ems1;
	public static Emiss A,N;
	
	 static Map<String, Emiss> m= new HashMap<String, Emiss>();
     public static Emiss get(String string) throws Exception{
         Emiss res =  m.get(string);
         if(res==null) throw new Exception("nothing for "+string);
         return res;
        
      }
	 public  static Emiss a(){return ems[1];}
     public  static Emiss N(){return ems[0];}
     public  static Emiss b(){if(ems.length>2) return ems[2]; else return null;}  
     public static EmissionStateSpaceForNoAlleles emiss;// = new EmissionStateSpaceForNoAlleles[0];
   //  public static EmissionStateSpace[] mergedSpace = new EmissionStateSpace[Constants.alleles.length];
    // new EmissionStateSpace[datac instanceof MergedDataCollection? ((MergedDataCollection)datac).ldl.length : 1];
     static{
    		// List<Integer> allelesL = new ArrayList<Integer>();
    		 int alleles1 = Constants.maxAlleles();
    		 Logger.global.info("alleles "+alleles1);
         	ems = new Emiss[alleles1+1];
           	ems1 = new Emiss[alleles1];
           	ems[0] = new Emiss('_',0, "_");
           	for(int i=1; i<ems.length; i++){
           		char ch =conv.get(i);
           		ems[i]=new Emiss(ch,i, ch+"");
           		ems1[i-1] = ems[i];
           	}
           	A = ems[1];
           	N = ems[0];
           	        for(int i=0; i<ems.length; i++){
           	        System.err.println("putting "+ems[i].toStringShort()+" : "+ems[i]);
           	        m.put(ems[i].toStringShort(), ems[i]);
           	        }
           	//for(int i=0; i<emiss.length; i++){
           	    /*  int index = allelesL.indexOf(Constants.alleles[i]);
           	      if(index>=0){
           	    	  emiss[i] = emiss[index];
           	    	 // mergedSpace[i] = mergedSpace[index];
           	      }
           	      else{*/
           	    	//  allelesL.add(Constants.alleles[i]);
           	    	  emiss = new EmissionStateSpaceForNoAlleles(alleles1);
          Emiss.ems = emiss.ems;
          Emiss.ems1 = emiss.ems1;
          Emiss.spaceByCN = emiss.spaceByCN;
          Emiss.mergedSpace = emiss.mergedSpace;
          
           	    	 // mergedSpace[i] = emiss[i].mergedSpace;
           	      //}
           	//}
     }
	
    final char c;
  //  public static Emiss[] ems,ems1;
   // final int i; //used for ordering
    public Emiss(char c, int i, String shortString){
     this.index=(short)i;
        this.c = c;
        this.shortString = shortString;
    }
    final String shortString;
    
    
    public String toStringShort(){
        return this.shortString;

   }

  final short index;
    public int noCopies(){
        return index==(short)0 ? 0 :1;
    }
    public String toStringPrint(){
        return c+"";
    }
    
    public char getChar(){
        return c;
    }
    public int noB() {
        if(index==(short)2 )
        	return 1;
        else 
        	return 0;
    }
    
    public int noB(int i) {
        if(index ==(short)i) return 1;
        else return 0;
    }
   
   
  private static Comparable[] getNumF(int numFounders){
      Comparable[] founders = new Comparable[numFounders];
      for(int i=0; i<founders.length; i++){
          founders[i] =new IntegerEmiss( i+1);
      }
      return founders;
  }
    
  private static Map<Integer, EmissionStateSpace>[] stateEmissionStateSpace = 
  new Map[] {
      new HashMap<Integer, EmissionStateSpace>(), new HashMap<Integer, EmissionStateSpace>() };
  
  public static EmissionStateSpace getStateEmissionStateSpace(int numF){
      EmissionStateSpace res = stateEmissionStateSpace[0].get(numF);
      if(res==null){
          stateEmissionStateSpace[0].put(numF, res = new SimpleEmissionStateSpace(getNumF(numF)));
      }
      return res;
  }
  
  public static EmissionStateSpace getStateEmissionStateSpace(int[] emStSp) {
      int numF = Constants.product(emStSp);
      EmissionStateSpace res = stateEmissionStateSpace[1].get(numF);
      if(res==null){
          EmissionStateSpace[] em = new EmissionStateSpace[emStSp.length];
          for(int i=0; i<em.length; i++){
              em[i] = getStateEmissionStateSpace(emStSp[i]);
          }
       //   if(emStSp.length==2)
          stateEmissionStateSpace[1].put(numF, res = new CompoundEmissionStateSpace(em,false, false));
      }
      return res;
   }
  
    
    
    
       
    
   
    public static class EmissionStateSpaceForNoAlleles{
    	int noAlleles;
    	Emiss[] ems1,ems;
    	EmissionStateSpace alleleSpace ;//=  new SimpleEmissionStateSpace( ems1);
    	
    	public  EmissionStateSpace[] spaceByCN ;
        public  EmissionStateSpace mergedSpace;
        public CompoundEmissionStateSpace[] spaceByPloidy = new CompoundEmissionStateSpace[Constants.maxPloidy()];
        NCube spaceByCN2;
       
        
        public CompoundEmissionStateSpace getEmissionStateSpace(
    			int[] numCopies) {
    		Arrays.sort(numCopies);
    		Object val = spaceByCN2.get(numCopies);
    		if(val==null){
    		boolean allAlleleSpace = true;
    			boolean allEqual = true;
    			for(int i=0; i<numCopies.length; i++){
    				if(numCopies[i]!=1) allAlleleSpace = false;
    				if(numCopies[i]!=numCopies[0]) allEqual = false;
    			}
    			/*if(allAlleleSpace && numCopies[0]>1){
    				val = spaceByCN[numCopies[0]];
    			}*/
    			//else{
    				EmissionStateSpace[] stsp = new EmissionStateSpace[numCopies.length];
    			//	SimpleEmissionStateSpace stsp1 = new SimpleEmissionStateSpace(stsp);
    			//	Arrays.fill(stsp, stsp1);
    				for(int i=0; i<stsp.length; i++){
    					stsp[i] = spaceByCN[numCopies[i]];
    				}
    				val = allEqual ?new CompoundEmissionStateSpace(stsp,false, false) : new CompoundEmissionStateSpace2(stsp,false);
    			//}
    			spaceByCN2.set(numCopies, val);
    		}
    		
    			return (CompoundEmissionStateSpace)val;
    		
    	}
        
        public  Emiss switchElement(Emiss comp) {
            if(comp==ems[1]) return ems[2];
            else if(comp==ems[2]) return ems[1];
            else if(comp==ems[0]) return ems[0];
            else throw new RuntimeException("!!");
          }
   
       
        
        
        public EmissionStateSpaceForNoAlleles(int alleles){
        	this.noAlleles = alleles;
        	ems = new Emiss[alleles+1];
           	ems1 = new Emiss[alleles];
           	ems[0] = Emiss.ems[0];
           	for(int i=1; i<ems.length; i++){
           		ems[i]=Emiss.ems[i];
           		ems1[i-1] = Emiss.ems1[i-1];
           	}
           	alleleSpace =  new SimpleEmissionStateSpace( ems1);
         	
           	        
        	
        	spaceByCN = new EmissionStateSpace[Constants.maxCopies()+1];
        	spaceByCN2 = new NCube(Constants.maxCopies(), Constants.maxCopies()+1);
        	spaceByCN[0] = new SimpleEmissionStateSpace(new Emiss[] {ems[0]});
        	spaceByCN[1] = alleleSpace;
        	for(int i=2; i<spaceByCN.length; i++){
        		EmissionStateSpace[] tmp = new EmissionStateSpace[i];
        		Arrays.fill(tmp, alleleSpace);
        		spaceByCN[i] = new CompoundEmissionStateSpace(tmp,false, false);
        	}
        String[] mod = Constants.modify(0);

        	boolean includeCopyZero = false;
        	for(int i=0; i<mod.length; i++){
        		if(mod[i].equals("0")) includeCopyZero=true;
        	}

        	
        	mergedSpace = new SimpleEmissionStateSpace(includeCopyZero ? spaceByCN :
        		Arrays.asList(spaceByCN).subList(1, spaceByCN.length).toArray(new SimpleEmissionStateSpace[0]));
        	for(int i=0; i<spaceByPloidy.length; i++){
        		EmissionStateSpace[] tmp = new EmissionStateSpace[i+1];
        		Arrays.fill(tmp, mergedSpace);
        		spaceByPloidy[i] = new CompoundEmissionStateSpace(tmp,false,false);
        	}
        }
        public  CompoundEmissionStateSpace getEmissionStateSpace(int i){
            return spaceByPloidy[i];
          }
          
          public  CompoundEmissionStateSpace getSpaceForNoCopies(int no_copies){
          	return spaceByPloidy[no_copies-1];
          }
          public  Emiss A,N;
          public Emiss translate(char c){
//            if(gaps.indexOf(c)>=0) return ems[0];
          //  if(c=='*') return star;
            for(int i=0; i<ems.length; i++){
                if(ems[i].c==c) return ems[i];
            }
         //   if(c=='0') return A;
          //  else if(c=='1') return B;
          
                throw new RuntimeException(new String(new char[]{c})+" "+(int)c);
            
        }	
          
    }
  
   
  
   
    /** can improve this class!!! */
  public static class NCube{
	  
	  SortedMap<int[], Object> m = new TreeMap<int[], Object>(CompoundEmissionStateSpace.comp1);
	  
	  
	  
	  //Object[] obj;
	  //int[] mult;
	  //int d; //dimension
	  //int n; //number perdimension
	  public NCube(int d, int n){
		//  this.obj=  new Object[(int)Math.pow(n, d)];
		 // mult = new int[d];
		  //for(int i=1; i<=mult.length; i++){
		//	  mult[mult.length-i] = (int) Math.pow(n, i-1);
		 // }
	  }
	  public void set(int[] ind, Object o){
		  /*Arrays.sort(ind);
		  int index =0;
		  for(int i=0; i<ind.length; i++){
			  index += mult[i]*ind[i];
		  }
		  obj[index] =o; */
		  m.put(ind, o);
	  }
	  public Object get(int[] ind){
		  return m.get(ind);
/*		  Arrays.sort(ind);
		  int index =0;
		  for(int i=0; i<ind.length; i++){
			  index += mult[i]*ind[i];
		  }
		  return obj[index]; */
	  }
  }
    	
   
    
 //   private static EmissionStateSpace stateSpace = new CompoundEmissionStateSpace(spaceByCN);
       
        
    
                
    
    public static String gaps =  new String(new char[] {'_', (char)0});
 //   public static sSpace = new Em
    
  
    
    
    
   
   
    /** i =  noCopies - 1 
    public static CompoundEmissionStateSpace getEmissionStateSpace(int i, EmissionStateSpace base, EmissionStateSpace[] store){
        if(store[i]==null){
            if(i==0 || i==1){
                EmissionStateSpace[] sp = new EmissionStateSpace[i+1];
                Arrays.fill(sp,base);
                store[i] = new CompoundEmissionStateSpace(sp,false);
            }
            else throw new RuntimeException("!!");
        }
        return (CompoundEmissionStateSpace) store[i];
    }*/
    
    public boolean equals(Object obj){
        if(!(obj instanceof Emiss)) return false;
        return c== ((Emiss)obj).c;
    }

   
    
    public int compareTo(Object o) {
    	if(! (o instanceof Emiss)) return -1;
   // 	if(!(o instanceof Emiss)) return false;
      char i1 = ((Emiss)o).c;
       if(c==i1) return 0;
       else return c<i1 ? -1 :1;
    }
    public String toString(){
        return c+"";
      
    }
    
  
    
   
    

    public ComparableArray expand() {
        return ComparableArray.make(this);
    }

  /*  public static int indexOf(Emiss n2) {
      for(int i=0; i<stateSpace.size(); i++){
          if(stateSpace.get(i)==n2){
              return i;
          }
      }
      throw new RuntimeException("!");
    }*/

   

	

   

	public static class Convert{
		static int A = (int)'A';
		static int Z = (int) 'Z';
		static int _ = (int) '_';
		public int get(char ch){
			return l.indexOf(ch);
		}
		public char get(int i){
			return l.get(i);
		}
		List<Character> l = new ArrayList<Character>();
			Convert(){
				l.add('_');
				for(int i = A; i<=Z; i++){
					l.add((char) i);
				}
			 for(int i=31; i<126; i++){
				 if(!l.contains((char) i)){
					 l.add((char)i);
				 }
			 }
			}

		}

	
	
   
}
