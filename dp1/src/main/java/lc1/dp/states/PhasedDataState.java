package lc1.dp.states;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import lc1.CGH.Aberation;
import lc1.dp.data.representation.CSOData;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.IntegerDistribution;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;

public class PhasedDataState extends HaplotypeEmissionState implements PIGData{

    public Comparable[] phased;
	/** u controls the sampling.  A high value results in a closer copy
     * use st_to_init to initialise and st_to_pseudo to set pseudo counts
     *  */
    public PhasedDataState(PhasedDataState st_to_init,  String name) {
       super(st_to_init, name);
      
    }
   
 public PhasedDataState(PhasedDataState state_j) {
     super(state_j);
     this.phased = new Comparable[state_j.phased.length];
     for(int i=0; i<phased.length; i++){
    	 phased[i] = state_j.phased[i];
     }
 }
 public int noCopies(){
		return this.noCop();
	}
 public PhasedDataState(String name, int noSnps, EmissionStateSpace emStSp, short data_index){
     super(name, noSnps, emStSp, data_index);
     this.phased = new Comparable[noSnps];
 }
 
 

 public static String getName(CSOData[] unit, String sep){
     StringBuffer buf = new StringBuffer(unit[0].getName());
     for(int i=1; i<unit.length; i++){
         buf.append(unit[0].getElement(0) instanceof ComparableArray && ((ComparableArray)unit[0].getElement(0)).size()>1 ? sep : sep);
         buf.append(unit[i].getName());
     }
    return buf.toString();
 }
    
    //following is for PIGData
    
 public PhasedDataState(PIGData[] statesD, List<int[]> list) {
    this(statesD[0].getName(), list.size(), ((PhasedDataState)statesD[0]).getEmissionStateSpace(),(short) -1);
    for(int i=0; i<list.size(); i++){
        int[] ind = list.get(i);
        this.emissions[i] = ((PhasedDataState)statesD[ind[0]]).emissions[ind[1]];
    }
}
 
 
 public PhasedDataState(PhasedDataState[] unit, boolean merge, String join, EmissionStateSpace emstsp) {
	 super(getName(unit, join), unit[0].length(), emstsp, unit[0].data_index);
     int noCopies = 0;//unit.length;
   //  if(!merge){
         for(int i=0; i<unit.length; i++){
             noCopies+=unit[i].noCop();
         }
  //   }
     outer: for(int i=0; i<unit[0].length(); i++){
         List<Comparable> l = new ArrayList<Comparable>();
         for(int j=0; j<unit.length; j++){
             Comparable comp =(Comparable) unit[j].getElement(i);
             if(comp instanceof ComparableArray){
            ComparableArray comp1 =  (ComparableArray) comp;
            if(merge)l.addAll(comp1.elements());
            else l.add(comp);
             }
             else{
                 l.add(comp);
             }
         }
         this.addPoint(i,new ComparableArray(l));
     }
	}
 
public PhasedDataState clone(){
    return new PhasedDataState(this);
}
public PhasedDataState(String string, List<String> st,
        EmissionStateSpace emStSp, short s) {
    this(string, st.get(0).length(), emStSp, s);
    for(int i=0; i<st.get(0).length(); i++){
        if(st.get(0).charAt(i)==' ') continue;
        char[] point = new char[2*st.size()-1];
        for(int j=0; j<st.size(); j++){
            point[2*j]= st.get(j).charAt(i);
            if(j<st.size()-1) point[2*j+1] = ',';
        }
        String str = (new String(point)).replaceAll("_", "");
        
        int ind = emStSp.getHapl(str);
        this.emissions[i] = new IntegerDistribution(ind,emStSp);
    }
}



/* (non-Javadoc)
  * @see lc1.dp.data.representation.CSOData#getUncertainty(lc1.dp.states.EmissionState)
  */
 public double[] getUncertainty(EmissionState st){
     double[] cert = new double[this.length()];
     for(int i=0; i<length(); i++){
         ComparableArray compa = (ComparableArray)getElement(i);
         cert[i] = st.score(
                 this.emStSp.getHaploPairFromHaplo(((IntegerDistribution)emissions[i]).fixedInteger())
                 , i);
     }
     return cert;
 }
   
 /* (non-Javadoc)
  * @see lc1.dp.data.representation.PIGData#getDeletedPositions(lc1.dp.states.EmissionState)
  */
 public Collection<Aberation> getDeletedPositions(EmissionState st, int noCop, Boolean deletion) {
     int[][] cnts = new int[noCop][this.length()];
     double[] uncertainty = getUncertainty(st);
     for(int i=0; i<length(); i++){
         ComparableArray compa = (ComparableArray)getElement(i);
         for(int k=0; k<compa.size(); k++){
            Comparable em =  compa.get(k);
                cnts[k][i] = (em instanceof Emiss) ? ((Emiss)em).noCopies() :((ComparableArray)em).getGenotypeString().length();
         }
       
     }
     List<Aberation> res = new ArrayList<Aberation>();
     for(int k=0; k<cnts.length; k++){
        res.addAll( Aberation.getAberation(this.getName()+"_"+k, cnts[k], uncertainty, deletion));
     }
     return res;
 }
 
 
 
public String getStringRep(int start, int end) {
 StringBuffer sb = new StringBuffer();
 for(int i= start; i<=end; i++){
     sb.append(this.getStringRep(i));
 }
 return sb.toString();
    
}



public Set<Integer>[] getSwitches() {
   Set<Integer>[] switches = new TreeSet[this.noCop()];
   
   Comparable[] prev = new Comparable[switches.length];
   for(int j=0; j<prev.length; j++){
       prev[j] = ((ComparableArray)this.getElement(0)).get(j);
       switches[j] = new TreeSet();
   }
   for(int i=1; i<this.length(); i++){
       for(int j=0; j<prev.length; j++){
           Comparable nxt  = ((ComparableArray)this.getElement(i)).get(j);
           if(!prev[j].equals(nxt)){
               switches[j].add(i);
           }
           prev[j] = nxt;
       }
   }
   return  switches;
}
   
public PIGData recombine(Map<Integer, Integer> recSites, int starting_index) {
   if(this.noCop()!=2 && this.noCop()!=1) throw new RuntimeException("!!");
  
   PhasedDataState rec = new PhasedDataState(this.getName()+"r", this.length(),  ((CompoundEmissionStateSpace)emStSp).getMembers()[0], this.data_index);
   int index = starting_index;
   for(int i=0; i<this.length(); i++ ){
       if(recSites.containsKey(i)){
           index = 1-index;
       }
       ComparableArray comp = (ComparableArray) this.getElement(i);
      // Integer index = (Integer) ((ComparableArray)recSites.getElement(i)).get(0);
       rec.addPoint(i, ComparableArray.make(((ComparableArray)this.getElement(i)).get(index)));
   }
   return rec;
}

private String conc(List<String> names){
    StringBuffer sb = new StringBuffer();
    for(Iterator<String> it = names.iterator(); it.hasNext();){
        sb.append(it.next());
        if(it.hasNext())sb.append("_");
    }
    return sb.toString();
}
    
public PhasedDataState[] split(int[] is) {
    List<String> names = Arrays.asList(this.getName().split("_"));
   PhasedDataState[] res = new PhasedDataState[is.length];
   int[] from = new int[is.length];
   int[] to = new int[is.length];
   for(int j=0; j<res.length; j++){
       from[j] = j==0 ? 0 : to[j-1];
       to[j] = from[j]+is[j];
       res[j] = new PhasedDataState(conc(names.subList(from[j], to[j])), this.length(), this.emStSp, this.data_index);
   }
   for(int i=0; i<this.length(); i++){
       ComparableArray compa = (ComparableArray) this.getElement(i);
       for(int j=0; j<is.length; j++){
           ComparableArray compa_1 = new ComparableArray(compa.elements().subList(from[j], to[j]));
           res[j].addPoint(i,compa_1);
       }
   }
   return res;
}
   
   
//if elements are TrComparableArrays, it returns the index rather than the state!
/* (non-Javadoc)
 * @see lc1.dp.data.representation.CSOData#split()
 */
public PhasedDataState[] split() {
    //List<String> names = Arrays.asList(this.getName().split(";"));
 CompoundEmissionStateSpace ces = (CompoundEmissionStateSpace)this.getEmissionStateSpace();
 EmissionStateSpace[] em = ces.getMembers();
  PhasedDataState[] res = new PhasedDataState[((ComparableArray)this.getElement(0)).size()];
//  int[] from = new int[res.length];
 // int[] to = new int[res.length];
  for(int j=0; j<res.length; j++){
   //   from[j] = j==0 ? 0 : to[j-1];
    // Comparable compa =  ((ComparableArray)this.getElement(0)).get(j);
     // to[j] = from[j]+
     // (compa instanceof ComparableArray ? 
      //((ComparableArray)compa).size() : 1);
      //String name = this.getName()+"_"+j;
      res[j] = new PhasedDataState(name, this.length(), em[j], this.data_index);
      Double[] phenV = this.phenValue;
      res[j].setPhenotype(phenV);
  }
  for(int i=0; i<this.length(); i++){
      int index = this.getFixedInteger(i);
    int[] ind =  ces.getMemberIndices(index);
     
      for(int j=0; j<ind.length; j++){
          res[j].emissions[i] = new IntegerDistribution(ind[j], ces.getMembers()[j]);
//          res[j].setFixedIndex(i, ind[j]);
//          res[j].addPoint(i,compa.getReal(j));
      }
  }
  return res;
}
   
    
   
   
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.CSOData#phaseCorrect(lc1.dp.data.representation.CompoundScorableObject, java.util.Collection, java.util.Collection)
     */
        public int[] phaseCorrect(CSOData ph, Collection<Integer> noCopiesL,
        Collection<Integer> noCopiesR , PrintWriter out, List<Integer> loc, double[] cert  , double thresh     
        ){
            PhasedDataState original = (PhasedDataState)ph;
            Set<Integer> poss =new HashSet<Integer>();//compar1);
            int correct = 0;
            int wrong = 0;
            Integer[] phasePos = this.updatePhasePositions();//, includeGaps);
            outer: for(int jk=0; jk<phasePos.length;jk++){
                int i = phasePos[jk];
                if(cert!=null && cert[i]<thresh) continue;
                int hapl = this.getFixedInteger(i);
                int haploPair = emStSp.getHaploPairFromHaplo(hapl);
                if(haploPair!=emStSp.getHaploPairFromHaplo( original.getFixedInteger(i))
                        || !noCopiesR.contains((Integer) emStSp.getCN(haploPair))){
                 //   Logger.global.info("excluded "+emiss+" "+noCopiesR);
                    continue; 
                }
                poss.clear();
                int[] possib = emStSp.getHaploFromHaploPair(haploPair);
                addAll(poss, possib);
                inner: for(int jk1 = jk-1; jk1>=0 ; jk1--){
                    int i1 = phasePos[jk1];
                    int hapl_1 = this.getFixedInteger(i1);
                    int haploPair_1 = emStSp.getHaploPairFromHaplo(hapl_1);
                    if(haploPair_1!= emStSp.getHaploPairFromHaplo( original.getFixedInteger(i1))
                    ||  !noCopiesL.contains((Integer) emStSp.getCN(haploPair_1))   
                    ) {
                     //   Logger.global.info("excluded "+emiss1+" "+noCopiesL);
                        continue outer; //note this means we do not worry about relative phase if the call was wrong, we could
                        //also do continue inner, if we want to phase with the previous position
                    }
                    
                          
                    phaseRelative(original, i, i1, poss);
                   if(!poss.contains(hapl)){
                       out.println(this.getName()+" "+loc.get(i)+" "+loc.get(i1)+" wrong "+this.getElement(i1)+"->"+this.getElement(i)+" cf "+original.getElement(i1)+"->"+original.getElement(i)+" "+poss);
                       wrong++;
                       break inner;
                   }
                   else{
                	   out.println(this.getName()+" "+loc.get(i)+" "+loc.get(i1)+" correct"+this.getElement(i1)+"->"+this.getElement(i)+" cf "+original.getElement(i1)+"->"+original.getElement(i)+" "+poss);
                           correct++;
                       //if(poss.size()!=1) throw new RuntimeException("!! check");
                      // break inner;
                   }
                  if(poss.size()==1)  break inner;
                }
              
            }
         //  System.err.println("comp "+wrong+" "+correct);
           
            
            return new int[] {wrong, correct+wrong};
        }
        
        public static int countDiff(Comparable o1, Comparable o2){
            ComparableArray obj1 = (ComparableArray)o1;
            ComparableArray obj2 = (ComparableArray)o2;
            int cnt = 0;
            for(int i=0; i<obj1.size(); i++){
                 if(!obj1.get(i).equals(obj2.get(i))) cnt++;
            }
            return cnt;
        }
        
        
        public static void printIdLine(String idLine, PrintWriter pw, int len, char end){
            char[] ch = new char[Math.max(idLine.length(),len)];
            Arrays.fill(ch, ' ');
            System.arraycopy(idLine.toCharArray(), 0, ch, 0, idLine.length());
          int jmp =10;
          int jmp2 = 100;
            for(int i=0; i<ch.length; i+=jmp){
                if(Math.IEEEremainder(i, jmp2)==0 && idLine.length()+i < ch.length){
                    System.arraycopy(idLine.toCharArray(), 0, ch, i, idLine.length());
                }
                if(i>=idLine.length()){
                    ch[i] = '|';
                }
            }
            pw.print((new String(ch))+end);
        }
        
        /* (non-Javadoc)
         * @see lc1.dp.data.representation.CSOData#isNested()
         */
        public boolean isNested() {
            if(this.getElement(0) instanceof ComparableArray){
         return (((ComparableArray)this.getElement(0)).get(0)) instanceof ComparableArray;
            }
            return false;
        }
        /* (non-Javadoc)
         * @see lc1.dp.data.representation.CSOData#print(java.io.PrintWriter, boolean, boolean, boolean, java.util.Collection)
         */
        public void print(PrintWriter pw, boolean idline, boolean expand, boolean mark, Collection<Integer> toDrop
        		,List<Character> alleleA, 
                List<Character> alleleB, PseudoDistribution[] unc
        ){
           char eol = '\t';
          
           int ploidy = ((ComparableArray)this.emStSp.defaultList.get(0)).size();
          for(int i=0; i<ploidy; i++){
        	  print1(pw,  idline, expand, mark,  toDrop, alleleA, alleleB, i, eol, unc);
          }
         
        }
        
        /* (non-Javadoc)
         * @see lc1.dp.data.representation.CSOData#print(java.io.PrintWriter, boolean, boolean, boolean, java.util.Collection)
         */
        public void print_old(PrintWriter pw, boolean idline, boolean expand, boolean mark, Collection<Integer> toDrop
        		,List<Character> alleleA, 
                List<Character> alleleB		
        ){
       char eol = '\t';
           if(idline) {
               if(mark) printIdLine("# id "+this.getName(), pw, this.length(), eol);
               else pw.println(("# id "+this.getName()));
           }
           if(toDrop!=null) throw new RuntimeException("!!");
         
        }
        
        
        
        public Comparable getElement(int i) {
            return this.phased[i];
        }
        
        public void print1(PrintWriter pw, boolean idline, boolean expand, boolean mark, Collection<Integer> toDrop
        		,List<Character> alleleA, 
                List<Character> alleleB, int chrom, char eol, PseudoDistribution[] unc
        ){
        	//pw.print("testing");
        //	if(emissions.length!=alleleA.size()) throw new RuntimeException("diff sizes");
        	 if(idline) {
                 if(mark) printIdLine("# id "+this.getName()+":"+chrom, pw, this.length(), eol);
                 else pw.println(("# id "+this.getName()));
             }
            Object[][] l =new Object[this.maxCopies(chrom)][];
            for(int i=0; i<l.length; i++){
                l[i] = new Object[this.length()];
            }
            for(int i=0; i<length(); i++){
                List<Comparable> list = new ArrayList<Comparable>();
                Comparable compa = ((ComparableArray)this.getElement(i)).get(chrom);
               if(compa instanceof ComparableArray) ((ComparableArray) compa).addObjectsRecursive(list);
               else list.add(compa);
                int k=0;
                for(Iterator<Comparable> it = list.iterator(); it.hasNext();k++){
                    Object nxt = it.next();
                    l[k][i] = nxt;
                }
            }
            double th =Constants.imputedThreshGraph(0); 
            for(int j=0; j<l.length; j++){
                Object[][] res = expand ? expand(l[j]) : new Object[][] {l[j]};
                for(int k=0; k<res.length; k++){
                    for(int i=0; i<res[k].length; i++){
                        if(toDrop==null || !toDrop.contains(i)){
                        	double p =1.0;
                        	String top;
                        	if(unc!=null && unc[i]!=null){
                        		PseudoDistribution dist = unc[i];
                        		int max = dist.getMax();
                        	   p = dist.probs(max);
                        	}
                        	if(p>th){
                        		top  =printElement( res[k][i], expand);
                        	
	                        	if(alleleA!=null && alleleA.size()>0 && Constants.replaceAB()){
	                        		Character torep = alleleA.get(i);
	                        		Character torep1 = alleleB.get(i);
	                        		if(torep!=null && torep1 !=null){
	                        			top =top.replace('A', torep).replace('B', torep1); 
	                        		}
	                        	}
                        	}else{
                        		top = "N";
                        	}
                            
                        	pw.print(top);
                        }
                    }
                    pw.println();
                }
            }
        }
       private  int maxCopies() {
    	   int max=0;
    	   for(int i=0; i<length(); i++){
               List<Comparable> list = new ArrayList<Comparable>();
               Comparable compa = this.getElement(i);
               if(compa instanceof ComparableArray){
            	   ((ComparableArray) this.getElement(i)).addObjectsRecursive(list);
               }
               else{
            	   list.add(compa);
               }
               int k=0;
               for(Iterator<Comparable> it = list.iterator(); it.hasNext();k++){
                   Object nxt = it.next();
                 //  l[k][i] = nxt;
               }
               if(k>max) max = k;
           }
    	   return max;
			
		}
       
       private  int maxCopies(int chrom) {
    	   int max=0;
    	   for(int i=0; i<length(); i++){
               List<Comparable> list = new ArrayList<Comparable>();
               Comparable compa = ((ComparableArray)this.getElement(i)).get(chrom);
               if(compa instanceof ComparableArray){
            	   ((ComparableArray)compa).addObjectsRecursive(list);
               }
               else{
            	   list.add(compa);
               }
               int k=0;
               for(Iterator<Comparable> it = list.iterator(); it.hasNext();k++){
                   Object nxt = it.next();
                 //  l[k][i] = nxt;
               }
               if(k>max) {
            	   max = k;
               }
           }
    	   return max;
			
		}

	int maxNoCopies() {
			// TODO Auto-generated method stub
			return 0;
		}

	private static Object[][] transpose(Object[][]obj){
           int max = 0;
           for(int i=0; i<obj.length; i++){
               if(obj[i].length>max){
                   max = obj[i].length;
               }
           }
           Object[][] res = new Object[max][obj.length];
           for(int i=0; i<obj.length; i++){
               for(int j=0; j<obj[i].length; j++){
                   res[j][i] = obj[i][j];
               }
           }
           return res;
       }
        
        static Object[][] expand(Object[] comp){
           List<Object[]> l = new ArrayList<Object[]>();
           for(int i=0; i<comp.length; i++){
               if(comp[i] instanceof Emiss){
                   ComparableArray compa = ((Emiss)comp[i]).expand();
                   l.add(compa.elements().toArray());
              
               }
               
               else l.add(new Object[] {comp[i]});
           }
           return transpose(l.toArray(new Object[0][]));
        }
        /* (non-Javadoc)
         * @see lc1.dp.data.representation.CSOData#print(java.io.PrintWriter, boolean, boolean, java.util.Collection)
         */
        public void print(PrintWriter pw, boolean expand, boolean mark, Collection<Integer> toDrop,
        		List<Character> alleleA, 
                List<Character> alleleB,PseudoDistribution[] unc){
              // if(!expand) super.print(false, pw);
              print(pw, true, expand, mark, toDrop, alleleA, alleleB, unc);
           }
   
        

      
          
           
           
   
            /* (non-Javadoc)
             * @see lc1.dp.data.representation.CSOData#samplePhase(lc1.dp.emissionspace.EmissionStateSpace, java.util.List, int, int, java.util.Set)
             */
            public  double samplePhase(EmissionStateSpace emStSp, List<CSOData> spList, int i1, int i2, Set<Integer> poss) {
                if(poss.size()>2) throw new RuntimeException("!!");
            //    int len = this.length();
               
                int hap1 = this.getFixedInteger(i1);
               int hap2 = this.getFixedInteger(i2);
               int haplopairIndex1 =  emStSp.getHaploPairFromHaplo(hap1);
               int haplopairIndex2 =  emStSp.getHaploPairFromHaplo(hap2);
              Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
              for(Iterator<Integer> it = poss.iterator(); it.hasNext();){
                  counts.put(it.next(), 0);
              }
                int sumMatch=0;
                 for(int j=0; j<spList.size(); j++){
                     PhasedDataState data = (PhasedDataState)spList.get(j); 
                     int hap1_j = data.getFixedInteger(i1);
                     int haplopair_j = emStSp.getHaploPairFromHaplo(hap1_j);
                     if(haplopair_j==haplopairIndex1 ){
                         int hap2_j = data.getFixedInteger(i2);
                         int haplopair2_j = emStSp.getHaploPairFromHaplo(hap2_j);
                          if(haplopair2_j==haplopairIndex2){
                              if(hap2_j!=hap2){
                                  hap1_j = emStSp.flip(hap1_j);
                              }
                              counts.put(hap1_j, counts.get(hap1_j)+1);
                              sumMatch++;
                          }
                     }
                 
                 }
             
              int max_index = getMax(counts);
              double max = counts.get(max_index);
              double certainty = sumMatch==0 ? 0 : max / (double)sumMatch;
                     List<Integer> l1 = new ArrayList<Integer>();
                     for(Iterator<Integer> it = poss.iterator(); it.hasNext();){
                         int next = it.next();
                        if(counts.get(next) < max) it.remove();
                     }
                    if(poss.size()==0){
                        poss.addAll(l1);
                     }
                     return certainty;
             }
            
            /* checking phase of i relative to i1, so poss returns possible phases at i*/
            private void phaseRelative(CSOData orig, int i, int i1, Collection<Integer> poss){
                PhasedDataState original = (PhasedDataState)orig;
               // List<ComparableArray> comp = new ArrayList<ComparableArray>(poss);
              int hapl = this.getFixedInteger(i);
               int hapl_1 = this.getFixedInteger(i1);
               int o_hapl = original.getFixedInteger(i);
                int o_hapl_1 = original.getFixedInteger(i1);
                if(o_hapl_1!=hapl_1){
                    o_hapl = this.emStSp.flip(o_hapl);
                }
                    //int noDiff = countDiff(emStSp.getHapl(o_hapl), emStSp.getHapl(o_hapl_1));
                    for(Iterator<Integer> it = poss.iterator(); it.hasNext();){
                        if(it.next().intValue()!=o_hapl) it.remove();
                    }
            }
            
            private int getMax(Map<Integer, Integer> counts) {
                Iterator<Integer> it = counts.keySet().iterator();
                Integer key = it.next();
              while( it.hasNext()){
                  Integer key1 = it.next();
                 if(counts.get(key1)> counts.get(key)) key = key1; 
              }
              return key;
            }

            private int getMax(Integer[][] diff, int i) {
                int max =0;
               for(int k=1; k<diff.length; k++){
                   if(diff[k][i]>diff[max][i]){
                       max = k;
                   }
               }
               return max;
            }

            /* (non-Javadoc)
             * @see lc1.dp.data.representation.CSOData#samplePhase(java.util.List, double[], lc1.dp.emissionspace.EmissionStateSpace)
             */
            public void samplePhase(List<CSOData> spList, double[] uncertainty, EmissionStateSpace emStSp1_){
                if(((ComparableArray)this.getElement(0)).isNested()){
                
                    throw new RuntimeException("!!");
                }
                else{
                 //   EmissionStateSpace emStSp = Emiss.getEmissionStateSpace(spList.get(0).noCopies()-1);
                    Set<Integer> poss =new HashSet<Integer>();
                    Integer[] phasePos = this.updatePhasePositions();//, true);
                    inner: for(int jk=1; jk<phasePos.length;jk++){
                        poss.clear();
                        int i = phasePos[jk];
                        int hapl = this.getFixedInteger(i);
                        ComparableArray compa =(ComparableArray) this.phased[i];
                       EmissionStateSpace emStSp1 = emStSp1_==null ? Emiss.getSpaceForNoCopies(compa.size()) :
                    	   emStSp1_;
                        int haplPair = emStSp1.getHaploPairFromHaplo(hapl);
                        int[] possib = emStSp1.getHaploFromHaploPair(haplPair);
                        addAll(poss, possib);
                        for(int jk1 = jk-1; jk1>=0 && poss.size()>1; jk1--){
                            int i1 = phasePos[jk1];
                           // if(this.length()!=len)  throw new RuntimeException("!!");
                            uncertainty[i]*=samplePhase(emStSp1,  spList, i, i1, poss);
                          //  if(this.length()!=len)  throw new RuntimeException("!!");
                        }
                        if(poss.size()>0){
                            Integer res = poss.iterator().next();
                            ((IntegerDistribution)this.emissions[i]).setFixedIndex(res);
                        }
                       
                    }
                }
            }
           
            
    private void addAll(Set<Integer> poss, int[] possib) {
             for(int i=0; i<possib.length; i++){
                 poss.add(possib[i]);
             }
           //     if(poss.size()>2) throw new RuntimeException("!!");
            }

    public void set(int i, Comparable obj) {
       this.addPoint(i, obj);
        
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.CSOData#updatePhasePositions(lc1.dp.emissionspace.EmissionStateSpace)
     */
    public Integer[]  updatePhasePositions(){
        List<Integer> l = new ArrayList<Integer>();
        for(int i=0; i<this.length(); i++){
            int hapl = this.getFixedInteger(i);
            int haploPair = emStSp.getHaploPairFromHaplo(hapl);
            int[] poss = emStSp.getHaploFromHaploPair(haploPair);
            if(poss.length>1) l.add(i);
           }
       // }
       return l.toArray(new Integer[0]);
    }
   /*****NEEED TO CHECK THIS IS THE RIGHT THING TO GET**/
    public void addPoint(int i, Comparable i1) {
    	try{
        int res = this.emStSp.getHapl(i1);
    	this.phased[i] = i1;
        short d_i = emissions[i]!=null ? emissions[i].getDataIndex() : -1;
        this.emissions[i] = new IntegerDistribution(res, emStSp);
        emissions[i].setDataIndex(d_i);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
      //  System.err.println(i1);
       // if(Constants.CHECK && !i1.equals(emStSp.getHapl(res))) throw new RuntimeException("!!");
    }
    
   
    /*public boolean allNull() {
        // TODO Auto-generated method stub
        return false;
    }*/
    public Class clazz() {
        throw new RuntimeException("!!");
    }
   
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#getStringRep(int)
     */
    public  String getStringRep(int i){
        Comparable compa = this.getElement(i);
        if(compa instanceof ComparableArray  && ((ComparableArray)compa).size()==1){
            return ((ComparableArray)compa).get(0).toString();
        }
        return compa.toString();
    }
   
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#printElement(java.io.PrintWriter, java.lang.Object, boolean)
     */
    public String printElement(Object el, boolean expand){
        if(el==null) return (" ");
        else{ 
         return el instanceof Emiss ? 
                    (expand ? ((Emiss)el).toStringShort() : ((Emiss)el).toStringPrint()) : 
            el instanceof Integer ? ((Integer)el).toString((Integer)el, Constants.radix()) : el.toString();
        }
    }

   
    
    
   
   
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.CSOData#getUncertainty(lc1.dp.states.EmissionState)
     */
    public double getUncertainty(EmissionState st, int i){
        EmissionStateSpace emstsp = st.getEmissionStateSpace();
      
        return st.score(
                emstsp.getHaploPairFromHaplo(
                ((IntegerDistribution)emissions[i]).fixedInteger()), i);
    }

   
   
    public void mkTrCompArray() {
        throw new RuntimeException("!!");
        
    }
    public void print(PrintWriter pw, String prefix){
        StringBuffer sb1= new StringBuffer(prefix);
       
        for(int i=0; i<this.length(); i++){
          sb1.append(this.getElement(i));
        }
       pw.println(prefix+" "+sb1.toString());
       
    }
    public String toString(){
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        this.print(pw,this.getName());
//        StringBuffer sb = new StringBuffer(this.getName()+"\n");
        
        return sw.getBuffer().toString();//this.getName()+":"+this.l1;
    }
}
