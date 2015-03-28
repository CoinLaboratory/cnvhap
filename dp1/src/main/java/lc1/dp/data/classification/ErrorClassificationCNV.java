/**
 * 
 */
package lc1.dp.data.classification;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.emissionspace.EmissionStateSpace;

public class ErrorClassificationCNV extends ErrorClassificationAbstract{
     
    EmissionStateSpace ss;
    public ErrorClassificationCNV(EmissionStateSpace ss, double thresh, 
             PrintWriter log){
          super(getCNV(ss), thresh, null, null, log);
        }
    private static List<Comparable> getCNV(EmissionStateSpace ss){
        Set<Integer> s = new TreeSet<Integer>();
        s.add(0);
        s.add(1);
        for(int i=0; i<ss.size(); i++){
            s.add(((ComparableArray)ss.get(i)).noCopies(true));
        }
        return new ArrayList<Comparable>(s);
    }
    @Override
    protected void compare(ComparableArray orig, ComparableArray pred,  int[] fromTo) {
        fromTo[0] = orig.noCopies(true);
        fromTo[1] = pred.noCopies(true);
    }
    @Override
    protected int wobbleRoom(ComparableArray compA) {
       return 4;
    }
   
    public String getString(int j){
        return this.emSt.get(j).toString();
      
    }

   
          
}
