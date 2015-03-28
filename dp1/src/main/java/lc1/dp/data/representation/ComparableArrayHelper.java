package lc1.dp.data.representation;

import java.util.Comparator;
import java.util.Iterator;


public class ComparableArrayHelper {
   // public final  static Permutations perm = new Permutations();
    public static class OrderComparator implements Comparator<ComparableArray>{
        public int compare(ComparableArray c1, ComparableArray c2) {
            if(c1.size()!=c2.size()){
                throw new RuntimeException("sizes not equal!!! "+c1+ " cf "+ c2);
            }
            Iterator<Comparable> o1 = (c1).iterator();
            Iterator<Comparable>  o2 = c2.iterator();
          while(o1.hasNext()){
              Comparable n1 = o1.next();
              Comparable n2 = o2.next(); 
         //     if(n1 instanceof ComparableArray) throw new RuntimeException("!!");
              int res = n1.compareTo(n2);
              if(res!=0) return res;
           }
           return 0;
        }
        
    };
    public static Comparator ORDER = new OrderComparator();
}
