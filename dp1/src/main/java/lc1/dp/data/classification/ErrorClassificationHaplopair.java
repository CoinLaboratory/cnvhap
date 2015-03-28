/**
 * 
 */
package lc1.dp.data.classification;

import java.io.PrintWriter;
import java.util.Comparator;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;

public class ErrorClassificationHaplopair extends ErrorClassificationAbstract{
   
       public ErrorClassificationHaplopair(EmissionStateSpace ss, double thresh, PrintWriter log){
          super(ss, thresh, new ErrorClassificationGenotype(ss, thresh, log), ss.haplopairListSize(), log);
           
        }
       
        static Comparator comp = new Comparator<Comparable>(){
            public int compare(Comparable o1, Comparable o2) {
               if(o1==o2) return 0;
               else if(o1==null) return -1;
               else if(o2==null) return 1;
               else return o1.compareTo(o2);
            }
        };
      

        @Override
        protected int wobbleRoom(ComparableArray compA) {
            EmissionStateSpace ss = ((EmissionStateSpace)this.emSt);
           return ss.getHaplopairForGeno(ss.getGenotype(compA)).length;
        }
       
       @Override
       protected  void compare(ComparableArray orig, ComparableArray pred,  int[] fromTo){
         //  List<Comparable> orig1 = new ArrayList<Comparable>(orig.elements());
          // int[] swap = Emiss.swap;
         //  List<Comparable> pred1 = new ArrayList<Comparable>(pred.elements());
       //   Collections.sort(orig1);
       //   Collections.sort(pred1);
          //  for(int i=0; i<orig1.size(); i++){
              //  int cas = orig1.get(i)==orig_mod1.get(i)? 0:1;
               
                Integer from = ((EmissionStateSpace)emSt).getHaploPair(orig);
                Integer to = ((EmissionStateSpace)emSt).getHaploPair(pred);
                if(from==null){
                    throw new RuntimeException("no haplo pair for "+orig+" "+emSt.size()+" "+emSt);
                }
                /*if(prob[1]>prob[0]){
                    from = swap[from];
                    to = swap[to];
                }*/
                fromTo[0] = from;
                fromTo[1] = to;
           // }
        }



       public String getString(int j){
           Comparable comp  = this.emSt.get(j);
           if(comp instanceof ComparableArray){
               return ((ComparableArray)comp).toString();
           }
           else if(comp instanceof Emiss){
               return ((Emiss)comp).toStringShort();
           }
           else return comp.toString();
       }
     
     
    }