package lc1.dp.emissionspace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import lc1.dp.data.representation.ComparableArray;
import lc1.util.CopyEnumerator;

public class CompoundEmissionStateSpace2 extends CompoundEmissionStateSpace {

	
	int[][] haploToTypeIndex;
	public int[] getMemberTypeIndices(int key) {
		return haploToTypeIndex[key];
	}
	public CompoundEmissionStateSpace2(EmissionStateSpace[] stateSpaces, boolean onlyRepeats
			) {
		  this.members = stateSpaces;
	        init(initStateSpace(stateSpaces, onlyRepeats, false));
	      
	        haploToMember = new int[haploList.size()][];
	        haploToTypeIndex = new int[haploList.size()][];
	        for(int i=0; i<this.haploList.size(); i++){
	            ComparableArray compa = (ComparableArray) backTranslate(this.haploList.get(i));
	            haploToMember[i] = new int[compa.size()];
	            haploToTypeIndex[i] = new int[compa.size()];
	            for(int ij=0; ij<compa.size(); ij++){
	            	
	            	//inner: for(int k=0; k<members.length; k++){
	            		String va = null;//
	            		Integer val = null;
	            		int k=0;
	            		Comparable comp_ = compa.get(ij);
	            	    for(; val==null; k++){
	            	    	va = members[k].getHaploString(comp_);
	            	    	val = members[k].getHapl(va);
	            	    }
	            	    k--;
	            		if(val!=null){
	            			haploToMember[i][ij] = val;
	            			haploToTypeIndex[i][ij] = k;
	            			//if(k==ij) break inner;
	            		}
	            		else{
	            			throw new RuntimeException("!!");
	            			/* 
	            			k = 1-ij;
	            			 va = members[k].getHaploString(compa.get(ij));
		            		 val = members[k].getHapl(va);
		            			haploToMember[i][ij] = val;
		            			haploToTypeIndex[i][ij] = k;
		            			*/
	            			
	            		}
	            		if(k>=members.length) 
	            			throw new RuntimeException("!!");
	            	//}
	            }
	        }
	        for(int i=0; i<this.haploToMember.length; i++){
	            int[] members = haploToMember[i];
	            membersIndexToIndex.put(members, i);
	        }
	        for(int i=0; i<this.defaultList.size(); i++){
	            int[] hapL = this.getHaps(i);
	            for(int i1=1; i1<hapL.length; i1++){
	                if(this.haploList.get(i1).toString().equals(this.haploList.get(0).toString())) throw new RuntimeException("!!");
	            }
	        }
	        SortedSet<Integer> s = new TreeSet<Integer>();
	        for(int i=0; i< this.haploPairToHaplo.length; i++){
	            s.add(haploPairToHaplo[i].length);
	        }
	}
	


	@Override
	 List<Comparable> initStateSpace(final List<Comparable>[] stateSpaces1, final boolean onlyRepeats, final boolean limitParents){
        final List<Comparable> set = new ArrayList<Comparable>();
       final  List<Comparable>[] stateSpaces = new List[stateSpaces1.length];
        for(Iterator<int[]> it = getReOrderings(stateSpaces.length).iterator(); it.hasNext();){
        	  int[] nxt= it.next();
        	  for(int j=0; j<nxt.length; j++){
        		  stateSpaces[j] = stateSpaces1[nxt[j]];
        	  }
        	
              CopyEnumerator posm = new CopyEnumerator(stateSpaces.length){
                  public Iterator<Comparable> getPossibilities(int depth) {
                    return stateSpaces[depth].iterator();
                  }
                  @Override
                  public void doInner(Comparable[] list) {
                     // if(!Constants.exclude() || !exclude(list)){ 
                      set.add(
                              translate(new ComparableArray(Arrays.asList(list))));//PairEmissionState.perm.get(res).get(0));
                      //}
                  }
                
                  @Override
                  public boolean exclude(Comparable[] list) {
                	if(!onlyRepeats)
                      return false;
                	else{
                		for(int i=1; i<list.length; i++){
                			if(!list[i].equals(list[0])) return true;
                		}
                		return false;
                	}
                  }
                @Override
                  public boolean exclude(Object obj, Object previous) {
                      return false;
                  }
             };
             posm.run();
        }
          //  List<ComparableArray> ss =  new ArrayList<Comparable>(set);
             return set;
    }
	 /** returns all possible orderings if arrays of length i - e.g. [0,1] [1,0] */
    private List<int[]> getReOrderings(int length) {
    	if(length==1){
    		List<int[]> res = new ArrayList<int[]>();
    		res.add(new int[] {0});
    		return res;
    	}
    	List<int[]> l = new ArrayList<int[]>();
		for(int i=0; i<length; i++){
			int[] left = new int[length-1];
			int k1=0;
			for(int k=0; k<length; k++){
				if(k!=i){
					left[k1] = k;
					k1++;
				}
			}
			List<int[]> l1= getReOrderings(length-1);
			for(int j=0; j<l1.size(); j++){
				int[] res_1 = l1.get(j);
				int[] res_new = new int[length];
				
				res_new[0] = i;
				for(int k=1; k<res_new.length; k++){
					res_new[k] = left[res_1[k-1]];
				}
				l.add(res_new);
			}
		}
		return l;
	}
}
