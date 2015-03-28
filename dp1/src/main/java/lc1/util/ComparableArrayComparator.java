/**
 * 
 */
package lc1.util;

import java.util.Comparator;

public class ComparableArrayComparator implements Comparator{
    public ComparableArrayComparator(){
        
    }
    public int compare(Object obj1, Object obj2) {
        Comparable[] c1 = (Comparable[])obj1;
        Comparable[] c2 = (Comparable[])obj2;
        if(c1.length!=c2.length){
            throw new RuntimeException("sizes not equal!!!");
        }
      for(int i=0; i<c1.length; i++){
          Comparable n1 = c1[i];
          Comparable n2 = c2[i];
           if(n1!=n2){
               if(n1==null) return 1;
               else if(n2==null) return -1;
               else return (n1.compareTo(n2));
           }
       }
       return 0;
    }
    
}