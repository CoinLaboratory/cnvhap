package lc1.simulation;

import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

import JSci.maths.statistics.ProbabilityDistribution;
import cern.jet.random.Uniform;

public class ExtNode {
  /* private double bl;
   public double bl(){
       return bl;
   }*/
   public String toString(){
       return loc.toString();
   }
   // BigInteger pos;
   ExtNode(double bl, int start, int end){
       this.bl = new TreeMap<Double,Location>();
       
       this.loc = new RangeLocation(start, end);
       this.bl.put(bl,loc );
   }
   Location [] mutations = new Location[3]; //ordered by position within eachof 3 categories
 
  Location loc;
  Map<Double, Location> bl;
   
 /*   public int compareTo(Object o) {
       return Double.compare(bl, ((ExtNode)o).bl);
    }*/
    
 /*   public List<Integer[]>[] extract( int start, int end){
        List<Integer[]>[] res = new List[3];
        for(int category =0; category < 3; category++){
        List<Integer[]> l = mutations[category];
        if(l==null || l.size()==0 ){
            res[category] = null;
        }
        else{
            
            int st =0;
            for(; st < l.size(); st++){
                if(l.get(st)[0]+l.get(st)[1]>=start) break;
            }
            int en =st;
            for(; en < l.size(); en++){
                if(l.get(en)[0]>end) break;
            }
            res[category] =  l.subList(st, en);
           // System.err.println("got one!");
        }
        }
        return res;
    }*/
    //public List<Integer>
static Comparator comp = new Comparator<Integer[]>(){

    public int compare(Integer[] o1, Integer[] o2) {
       return o1[0].compareTo(o2[0]);
    }

    
};
    public void simulate(int[] noEvents, Uniform uniform, ProbabilityDistribution[] lengthDist,double length) {
        for(int j=0; j<noEvents.length; j++){
        if(noEvents[j]>0){
//            mutations[j] = 
         //   Integer [][] n = new Integer[noEvents[j]][];
            for(int i=0; i<noEvents.length; i++){
                //Integer[] n = new Integer[2];
                for(int k=0; k<noEvents[i]; k++){
                    int st =(int) Math.round( uniform.nextDouble()* ( length)); //location
                    int len = (int) Math.round(lengthDist[j].inverse(uniform.nextDouble())); //length
                    Location mut = new RangeLocation(st, st+len-1);
                    if(mut.overlaps(this.loc)){
                        Location intersection = mut.intersection(loc);
                        if(mutations[j]==null) mutations[j] = intersection;
                        else mutations[j] = mutations[j].union(intersection);
                    }
                }
               // if(n[i][1]>1) System.err.println(n[i][0]+" "+n[i][1]+" "+j);
            }
         //   Arrays.sort(n, comp);
          //  mutations[j] = Arrays.asList(n);
          //  System.err.println(this.mutations[j]); 
        }
        }
        
    }
 /*   public void setbl(double bl2) {
       this.bl = bl2;
        
    }*/
    public void extend(double bl, int start, int end) {
    	Location l1 =new RangeLocation(start, end); 
        this.loc 
          = this.loc.union(l1);
        Location l2 = this.bl.get(bl);
        if(l2!=null){
        	this.bl.put(bl, l2.union(l1));
        }
        else{
        	this.bl.put(bl,l1);
        }
        //Logger.global.info("new loc "+loc);
        
    }
  
   
   
}
