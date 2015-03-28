/**
 * 
 */
package lc1.util;

import java.util.Comparator;

class ArrayComparator implements Comparator<Comparable[]> {
    public int compare(Comparable[] o1, Comparable[] o2) {
        for(int i=0; i<o1.length; i++){
            if(o1[i]==null && o2[i]==null) continue;
            else if(o1[i]==null) return -1;
            else if(o2[i]==null) return 1;
            else{
               int res = o1[i].compareTo(o2[i]);
               if(res!=0) return res;
           }
        }
        return 0;
     }
}