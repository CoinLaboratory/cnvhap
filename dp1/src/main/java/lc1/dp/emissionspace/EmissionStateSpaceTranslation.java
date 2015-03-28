package lc1.dp.emissionspace;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;


public class EmissionStateSpaceTranslation {
   // EmissionStateSpace big; ///bigger space!
   // EmissionStateSpace small; //smaller space;
    
    Map<Integer, Integer> bigToSmall = new HashMap<Integer, Integer>();
    Map<Integer, int[]> smallToBig = new HashMap<Integer, int[]>();
    
    public EmissionStateSpaceTranslation(EmissionStateSpace big, EmissionStateSpace small, boolean replaceB){
       
        Map<Integer, Set<Integer>> smallToBig1 = new HashMap<Integer, Set<Integer>>();
        String[] stri = getHapStrings(big);
        String[] strj = getHapStrings(small);
        if(stri.length!=strj.length){
        outer: for(int i=0; i<stri.length; i++){
            String st  = stri[i].replaceAll("X", "AA").replaceAll("Y", "BB").replaceAll("Z", "AB")
            .replaceAll("T", "AAA").replaceAll("U", "AAB").replaceAll("V", "ABB").replaceAll("W", "BBB");
            if(replaceB) st = st
            .replaceAll("B", "A")
            ;
             // big.getHaploPairString(big.get(i));
            for(int j=0; j<strj.length; j++){
                //String st1  = small.getGenotypeString(small.getGenotype(j));
                String st1  = strj[j].replaceAll("0", "_").replaceAll("1", "A").replaceAll("2", "AA").replaceAll("3", "AAA");
                if(st.equals(st1)){
                    bigToSmall.put(i, j);
                    Set<Integer> l = smallToBig1.get(j);
                    if(l==null) smallToBig1.put(j, l = new HashSet<Integer>());
                    l.add(i);
                    continue outer;
                }
            }
       //   throw new RuntimeException("nothing found for "+st+" in "+Arrays.asList(strj));
        }
        for(Iterator<Entry<Integer, Set<Integer>>> it = smallToBig1.entrySet().iterator(); it.hasNext();){
            Entry<Integer, Set<Integer>> nxt = it.next();
            int[] res = new int[nxt.getValue().size()];
            int i=0;
            for(Iterator<Integer> it1 = nxt.getValue().iterator();it1.hasNext(); i++){
                res[i] = it1.next();
            }
            smallToBig.put(nxt.getKey(), res);
        }
        }
        else{
            for(int i=0; i<stri.length; i++){
                bigToSmall.put(i, i);
                smallToBig.put(i, new int[]{i});
            }
        }
      
      //  Logger.global.info(bigToSmall+"");
    }
    
    private static String[] getHapStrings(EmissionStateSpace small) {
        String[] res = new String[small.size()];
        for(int j=0; j<small.size(); j++){
            //String st1  = small.getGenotypeString(small.getGenotype(j));
           res[j] = small.getHaploPairString(small.get(j));
         
        }
        return res;
    }

    public Integer bigToSmall(int i){
        return bigToSmall.get(i);
    }
    public int[] smallToBig(int i){
        return smallToBig.get(i);
    }
}
