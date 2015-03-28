/**
 * 
 */
package lc1.dp.data.classification;

import java.io.PrintWriter;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.emissionspace.EmissionStateSpace;

public class ErrorClassificationGenotype extends ErrorClassificationAbstract{
    
    
    public ErrorClassificationGenotype(EmissionStateSpace ss, double thresh, PrintWriter log){
          super(ss, thresh,  
                  new ErrorClassificationCNV(ss, thresh,null),
                  ss.genoListSize(),  log);
        }
    @Override
    protected void compare(ComparableArray orig, ComparableArray pred,  int[] fromTo) {
       fromTo[0] = ((EmissionStateSpace)this.emSt).getGenotype(orig);
       fromTo[1] =  ((EmissionStateSpace)this.emSt).getGenotype(pred);
    }
    @Override
    protected int wobbleRoom(ComparableArray compA) {
        EmissionStateSpace ss = ((EmissionStateSpace)this.emSt);
       return ss.getGenoForCopyNo(ss.getCN(compA)).length;
    }
   
    public String getString(int j){
        return ((EmissionStateSpace)emSt).getGenotypeString(((EmissionStateSpace)emSt).get(j));
      
    }
    
          
}
