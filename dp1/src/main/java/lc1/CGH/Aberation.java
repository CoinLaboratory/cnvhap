package lc1.CGH;

import java.util.ArrayList;
import java.util.List;

public class Aberation implements Comparable{
    public int start;
    public int end;
    final public int copy;
    final public String name;
   public  double certainty;
    public Aberation(String name, int start, int copy){
       this.name = name;
        this.start = start;
        this.copy = copy;
    }
    public String toString(){
        return this.name+":"+this.start+":"+this.end+":"+this.copy+" "+this.certainty;
    }
    public static List<Aberation> getAberation(String name, int[] noCopies, double[] cert, Boolean deletion){
        int len = noCopies.length;
        List<Aberation> l = new ArrayList<Aberation>();
        for(int i=0; i<len; i++){
        	   int noCop =noCopies[i];
            if((noCop!=1 && deletion==null) ||
            		(deletion!=null && noCopies[i]<1 && deletion)
            			|| (deletion!=null && noCopies[i]>1 && !deletion)
            ){
             
                double cert_i = 0;
                double cnt =0;
                Aberation ab = new Aberation(name, i, noCop);
                for(;  i<len &&noCopies[i]==noCop ;i++ ){
                    cert_i+=cert[i];
                    cnt++;
                }
                ab.certainty = cert_i/cnt;
                    //Math.pow(cert_i, 1.0/(double) (ab.end-ab.start+1));
                i--;
                ab.end = i;
                l.add(ab);
            }
        }
        return l;
    }
    public int compareTo(Object o) {
       Aberation ab = (Aberation) o;
       if(start!=ab.start) return start < ab.start ? -1 : 1;
       else if(copy!=ab.copy) return copy < ab.copy ? -1 :1;
       else if(end!=ab.end) return end < ab.end ? -1 :1;
       else  return name.compareTo(ab.name);
    }
}
