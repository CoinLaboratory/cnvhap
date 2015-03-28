/**
 * 
 */
package lc1.sequenced;

import java.util.HashMap;
import java.util.Map;

public class Deletion implements Comparable{
    public int chr;
    public int start;
    public int end;
    public int id;
    public int mid(){
        return (int) Math.round(((double)(end + start))/2.0);
    }
    public String toString(){
        return id+"_"+chr+"_"+start+"_"+end;//+"_"+m;
    }
    Deletion(String st){
        String[] str = st.split("_");
        try{
        id = Integer.parseInt(str[2]);
        }catch(Exception exc){
            id = Integer.parseInt(str[1]);
        }
        end = Integer.parseInt(str[str.length-1]);
        start = Integer.parseInt(str[str.length-2]);
        chr = Integer.parseInt(str[str.length-3]);
    }
    public void add(String st1, String st2){
        if(st1.length()>1){
            m.put(st1,st2.indexOf("homo")>=0);
        }
    }
    public Map<String, Boolean> m = new HashMap<String, Boolean>();
    public int compareTo(Object o) {
       Deletion o1 = (Deletion)o;
       if(chr !=o1.chr) return chr < o1.chr ? -1 : 1;
       else if(start!=o1.start) return start <o1.start ? -1 :1;
       else if(end !=o1.end) return  end < o1.end ? -1 :1;
       else return 0;
    }
    public String toStringShort() {
        return id+"_"+chr+"_"+start+"_"+end;
    }
}