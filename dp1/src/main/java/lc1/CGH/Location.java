package lc1.CGH;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

public class Location implements Comparable{
    
    /** returns a sublist of locations with start betweent min and max inclusive
     * perc is percentage of loc
     *  */
    public List<Location> getStartOverlaps(Iterator<Location> it,  double perc){
        List<Location> res = new ArrayList<Location>();
        while(it.hasNext()){
            Location nxt = it.next();
            double overl = this.overlaps(nxt) ;
           // System.err.println(nxt+" "+overl);
            if(overl>perc * Math.min(nxt.size(), this.size())){
                res.add(nxt);
                
             
            }
        }
       return res;
    }
 //   List<String> names
    public String chr;
 public  long min;
 public   long max;
 public double logL=0;
 public String probeId = "";;
 public void setProbeId(String id){
     this.probeId = id;
 }
 public void setLogL(double logL){
     this.logL = logL;
 }
    private int noCop=-1;
   private String name = "";
public List<String >noObs = new ArrayList<String>();
   public void incrObs(Location i){
       noObs.addAll(i.noObs);
   }
    Properties ann;
   // Map<String, String> properties = new HashMap<String, String>();
    public long size(){
        return max - min+1;
    }
    public String toString(){
        StringBuffer sb = new StringBuffer();
        if(chr!="") sb.append(chr+":");
        sb.append(min);
        if(max!=min) sb.append("-"+max);
        if(noCop>=0) sb.append("/"+noCop);
        return sb.toString();
//        if(min==max) return "Chr"+chr+"_"+name+":"+min+"-"+max+"/"+noCop+";";//+noObs.toString();
  //      return "Chr"+chr+"_"+name+":"+min+"-"+max+"/"+noCop+";";//+noObs.toString();
    }
    
  /*  public Location subsquent(){
        return new Location(chr, this.max+1, this.max+1, 0);
    }*/
    
    /** what fraction of this noObs are contained in loc */
    public int getIndivOverlaps(Location loc) {
        int res = 0;
        for(Iterator<String> it = noObs.iterator(); it.hasNext();){
            if(loc.noObs.contains(it.next())) res++;
        }
        return res;
     }
    public Location(String st){
        String[] str = st.split("_");
        try{
        this.name = Integer.parseInt(str[2])+"";
        }catch(Exception exc){
            name = ""+Integer.parseInt(str[1]);
        }
        max = Integer.parseInt(str[str.length-1]);
        min = Integer.parseInt(str[str.length-2]);
        chr = ""+Integer.parseInt(str[str.length-3]);
    }
    public Location(String str, long min, long max){
        this.chr = str;
        this.min = min;
        this.max = max;
       
    }
   
    public Location(Collection<Location> overlap) {
        Iterator<Location> it = overlap.iterator();
        Location first = it.next();
        noCop = first.noCop;
        min = first.min;
        max = first.max;
        chr = first.chr;
        while(it.hasNext()){
            Location nxt = it.next();
            if(nxt.noCop!=noCop){
                throw new RuntimeException("!!");
            }
            if(nxt.min<this.min) min = nxt.min;
            if(nxt.max >this.max)max = nxt.max;
        }
    }
public Location(int start, int end, String chr){
    this.min = (long)start;
    this.max = (long)end;
    this.chr = chr;
}
    public Location(Location loc) {
      this.chr = loc.chr;
      this.name = loc.name;
      this.min = loc.min;
      this.max = loc.max;
      this.noCop = loc.noCop;
      this.noObs = new ArrayList(loc.noObs);
      this.noSnps = loc.noSnps;
    }
    public boolean contains(int pos){
        return pos>=min && pos <=max;
    }
    /*returns overlap size */
    public double overlaps(Location loc){
        if(loc.chr.equals(chr)){
            double min_interval = Math.max(min, loc.min);
            double max_interval = Math.min(max, loc.max);
            return   max_interval-min_interval+1;
        }
        else return Double.NEGATIVE_INFINITY;
    }
    
    public double overlaps(int pos) {
        double min_interval = Math.max(min, pos);
        double max_interval = Math.min(max, pos);
        return   max_interval-min_interval+1;
    }
    
    public Location getOverlap(Location loc) {
        if(loc.chr.equals(chr)){
           Location newLoc = new Location(chr,(int) Math.max(min, loc.min),
           (int)Math.min(max, loc.max));
           return newLoc;
        }
        else return null;
        
    }
    
   /* public int compareTo(Object o1) {
      
     }
~*/public boolean equals(Object o){
    return this.compareTo(o)==0;
}

    public int compareTo(Object o) {
        Location o2 = (Location) o;
        int res1 = this.chr.compareTo(o2.chr);
        if(res1!=0) return res1;
        if(this.min != o2.min) return this.min < o2.min ? -1 :1;
        else    if(this.max != o2.max) return this.max < o2.max ? -1 :1;
        else if(this.noCop!=o2.noCop) return this.noCop < o2.noCop ? -1 :1;
        else{
            return  this.name.compareTo(o2.name);
        }
    }
    public void plotCode(PrintWriter pw, boolean b, int mino, int maxo,int i) {
       pw.print((b ? "plot" :"lines" ));
       double off = i*0.2;
       pw.print(     "(x = c("+min+","+max+"), y = c("+(off)+","+(off)+"), col =  "+(i+1)+" ,lwd ="+2 );
       if(b) pw.print(     " ,xlim = c("+mino+","+maxo+"), ylim = c("+0+","+0.5+") " );
       if(b) pw.print(",xlab = \"Location\", ylab = \"Copy Number\"");
       if(b) pw.println(             ",type = \"l\")");
       else pw.println(")");
        
    }
    public void setNoCop(int noCop){
       // if(AberationFinder.includeCN()){
            this.noCop = noCop;
        //}
    }
    public void setName(String nme1) {
        if(nme1==null){
            this.name = "";
            this.noObs.add(name);
            return;
        }
        else{
            name = 
//                AberationFinder.includeName() ? 
                        nme1 
                       // : ""
                            ;
            this.noObs.add(nme1);
        }
       
    }
    public int noCop() {
       return noCop;
    }
    public String name() {
        return name;
    }
/*[number_of_individuals_in_this_loc_detected_by_other_method,
  total_num_of_indiv_in_this_loc,
  num_of_indiv_in_other_loc detected by this
  total_num of indiv in other method */  
    public int[] overl1; 
//overl1 is as a fraction of this total individuals
    public void setOverl(int[] is) {
       overl1 = is;
     
        
    }
    public String toStringPrint() {
      return this.chr+"_"+this.min+"_"+this.max;
    }
   public int noSnps =0;
    public void setNoSnps(int noSnps) {
        this.noSnps = noSnps;
    }
    public int mid() {
        return (int) Math.floor(((double)min+(double)max)/2.0);
    }
    public void setNeg() {
       this.min = -min;
       this.max = -max;
        
    }
	public void print(PrintWriter regions) {
		regions.println(this.min+"\t"+this.max+"\t"+this.name+"\t"+this.noCop+"\t"+this.noSnps);
	}
   
   
}
