package lc1.CGH;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import lc1.dp.data.collection.IlluminaRDataCollection;

public class Locreader {
  public final long lengthLimit;
  
  IlluminaRDataCollection dc;
  final String name;
    public SortedSet<String> keys = new TreeSet<String>();
    
    public SortedSet<Integer> probes;
    private Map<String, List<Location>> deletions = new HashMap<String, List<Location>>();
    private Map<String, List<Location>> amplifications = new HashMap<String, List<Location>>();
    public void sort(){
     for(Iterator<List<Location>> it = deletions.values().iterator(); it.hasNext();){
         Collections.sort(it.next());
     }
     for(Iterator<List<Location>> it = amplifications.values().iterator(); it.hasNext();){
         Collections.sort(it.next());
     }
    }
    
    
    public Integer[] detected(Locreader loc1, int overlp, PrintWriter out1){
        out1.println("how many does "+name+" find of "+loc1.name+" predictions");
        boolean[] del  = new boolean[] {false, true};
        Integer[] res = new Integer[] {0,0, 0,0};
        for(int i=0;i<del.length; i++){
            for(Iterator<String> it = loc1.keys.iterator(); it.hasNext();){
                detected(loc1, it.next(), del[i], res, overlp, out1);
            }
        }
        return res;
    }
    public void printDetails(PrintWriter out1, Location loc, String name){
        if(dc!=null){
            if(name.equals("") && false){
                for(int i=0; i<loc.noObs.size(); i++){
                    printDetails(out1, loc, loc.noObs.get(i));
                }
            }
            else{
            out1.println(this.name);
            List<Double> l = new ArrayList<Double>();
            List<Integer> l1 = new ArrayList<Integer>();
            List<Integer> l2 = new ArrayList<Integer>();
            List<Double> b = new ArrayList<Double>();
            double[] max = null;//dc.getR(loc, l, l1, l2, b, name);
            if(max==null) return;
            StringBuffer sb = new StringBuffer();
            StringBuffer sb1 = new StringBuffer();
            for(int i=0; i<l.size(); i++){
                sb.append("%7.3f ");
            }
            for(int i=0; i<l.size(); i++){
                sb1.append("%7i ");
            }
            out1.println("max "+max[0]+" "+max[1]);
            out1.println(String.format(sb.toString(), l.toArray()));
           if(b.size()>0) out1.println(String.format(sb.toString(), b.toArray()));
            out1.println(String.format(sb1.toString(), l1.toArray()));
            out1.println(String.format(sb1.toString(), l2.toArray()));
            }
        }
    }
    public void detected(Locreader loc1, String key, boolean del, Integer[] res, int overlapNoProbes, PrintWriter out1){
        List<Location> map1 =  del ? deletions.get(key): amplifications.get(key);
        List<Location> dectable = new ArrayList<Location>();
        SortedSet<Integer> detectableProbes = new TreeSet<Integer>();
        SortedSet<Integer> detectedProbes = new TreeSet<Integer>();
        detectable(loc1, key, del, dectable, detectableProbes, overlapNoProbes);
        List<Location> detected = new ArrayList<Location>();
        outer: for(int i=0; i<dectable.size(); i++){
            Location loc = dectable.get(i);
            if(map1!=null){
            for(int j=0; j<map1.size(); j++){
                Location loc_1 = map1.get(j);
                
                if(loc.overlaps(loc_1)>=0){
                    detectedProbes.addAll(noPosIn(loc.getOverlap(loc_1)));
                    detected.add(loc);
                    out1.println("detected "+loc+" "+loc_1);
                  printDetails(out1, loc, loc.name());
                  loc1.printDetails(out1, loc, loc.name());
                    continue outer;
                }
            }
            
            }
            out1.println("not detected "+loc);
            printDetails(out1, loc, loc.name());
            loc1.printDetails(out1, loc, loc.name());
        }
       
       res[0]+=detected.size();
       res[1]+=dectable.size();
       res[2]+=detectedProbes.size();
       res[3]+=detectableProbes.size();
     //  out1.println("detected  "+detected);
     //  out1.println("detectable "+dectable);
     //  out1.println(key+" "+del+" :  "+Arrays.asList(res));
       out1.flush();
      // if(this.name.equals("snp") || this.name.equals("french")){
      //     System.err.println("detecable probes "+name+" "+probes.size()+" "+probes);
      // }
    }
    
    public void detectable(Locreader loc1, String key, boolean del, List<Location> detecable, SortedSet<Integer> detectableProbes, int overlapNoProbes){
        List<Location> map2 =  del ? loc1.deletions.get(key): loc1.amplifications.get(key);
      if(map2==null) return;
       for(int i=0; i<map2.size(); i++){
           Location loc = map2.get(i);
           SortedSet<Integer> noPos = noPosIn(map2.get(i));
           if(noPos!=null&& noPos.size()>=overlapNoProbes){
               detecable.add(loc);
               detectableProbes.addAll(noPos);
           }
       }
    }
    
    public SortedSet<Integer> noPosIn(Location loc){
        SortedSet<Integer> tail = probes.tailSet((int)loc.min);
        if(tail==null) return null;
        else return tail.headSet((int)loc.max+1);
    }
    
    public int number(){
        int num=0;
        for(Iterator<String> it = keys.iterator(); it.hasNext();){
            num+=number(it.next());
        }
        return num;
    }
    public int number(String st) {
        
        List<Location> loc = deletions.get(st);
        List<Location> loc2 = amplifications.get(st);
        int num = loc==null ? 0: loc.size();
        num+=loc2==null ? 0 :loc2.size();
        return num;
     }
public int number(String st, int noCop) {
        List<Location> loc = noCop==0 ? deletions.get(st) : amplifications.get(st) ;
        return loc==null ? 0: loc.size();
       
     }
    public long totalSize(String key){
       long size =0;
        for(Iterator<Location> it = iterator(key); it.hasNext();){
            size+=(it.next().size());
        }
        return size;
    }
    public long totalSize(String key, int nocop){
        long size =0;
         for(Iterator<Location> it = iterator(key, nocop); it.hasNext();){
             size+=(it.next().size());
         }
         return size;
     }
    public long[] getTotal(){
        long[] res = new long[]{0,0};
        for(Iterator<Location> it = this.iterator(); it.hasNext();){
            Location loc = it.next();
            if(loc.noCop()==0) res[0]+= loc.max-loc.min;
            else res[1]+= loc.max-loc.min;
        }
        return res;
    }
    public Location getFirst(String name){
        Location loc1 = getFirst(0, name);
        Location loc2 = getFirst(2, name);
        if(loc1==null) return loc2;
        else if(loc2==null) return loc1;
        return loc1.compareTo(loc2) < 0 ? loc1 : loc2;
    }
    public Location getFirst(){
        SortedSet<Location > l = new TreeSet<Location>();
        for(Iterator<String> it = this.keys.iterator(); it.hasNext();){
            Location loc = getFirst(it.next());
            if(loc!=null) l.add(loc);
        }
        return l.first();
    }
    public Location getLast(){
        SortedSet<Location > l = new TreeSet<Location>();
        for(Iterator<String> it = this.keys.iterator(); it.hasNext();){
            Location loc = getLast(it.next());
            if(loc!=null) l.add(loc);
        }
        return l.last();
    }
    public Location getLast(String name){
        Location loc1 = getLast(0, name);
        Location loc2 = getLast(2, name);
        if(loc1==null) return loc2;
        else if(loc2==null) return loc1;
        return loc1.compareTo(loc2) > 0 ? loc1 : loc2;
    }
    
    public Location getFirst(int noCop, String nm){
        if(noCop==0) return deletions.get(nm).size()>0 ? deletions.get(nm).get(0) : null;
        else return amplifications.get(nm).size()> 0?  amplifications.get(nm).get(0): null;
    }
    public Location getLast(int noCop, String name){
        if(noCop==0) return deletions.get(name).size()>0 ? deletions.get(name).get(deletions.get(name).size()-1) : null;
        else return amplifications.get(name).size()> 0? amplifications.get(name).get(amplifications.get(name).size()-1) : null;
    }
    public void addAll(Locreader loc1){
        for(Iterator<Location> it = loc1.iterator(); it.hasNext();){
            add(new Location(it.next()));
        }
    }
    public void addAll(Iterator<Location> it){
        for(; it.hasNext();){
            add(it.next());
        }
    }
    public Locreader(long lengthLim, String name){
        this.lengthLimit = lengthLim;
        this.name = name;
    }
   /* public void addAll(Locreader loc1, int num){
        for(Iterator<Location> it = loc1.iterator(); it.hasNext();){
            Location loc = it.next();
            if(loc.noSnps>=num){
                add(loc);
            }
        }
    }*/
    
    int getMax() {
       return (int) getLast().max;
     }
    
    int getMin() {
      
        return (int) getFirst().min;
     }
     void printPlotCode(PrintWriter pw, int i) {
        int maxo = getMax();
        int mino = getMin();
      //  pw.println("pdf(file = \"rplot.pdf\")");
        Iterator<Location> it = iterator();
        if(it.hasNext()){
        it.next().plotCode(pw, true && i==0, mino, maxo, i);
        while(it.hasNext()){
            it.next().plotCode(pw, false, 0, maxo, i);
        }
        }
    }
     
     
    public void merge(double frac){
        Logger.global.info("merging ");
       // Logger.global.info("before merging "+this.number()+" "+this.getClass());
         for(Iterator<String> it = this.keys.iterator(); it.hasNext();){
             String nm = it.next();
             merge(0, nm, frac);
             merge(2, nm, frac);
         }
         Logger.global.info("finished merging ");
       //  Logger.global.info("after merging "+this.number()+" "+this.getClass());
     }
    
    public void mergeNames(){
        int[] ind = new int[] {0,2};
        for(int j=0; j<ind.length;  j++){
            int i = ind[j];
            Map<String, List<Location>> m = i==2 ? amplifications : deletions;
            List<Location> list = new ArrayList<Location>();
            
            for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
               String key = it.next();
                   List<Location> loc = m.get(key);
                   for(Iterator<Location> it1 = loc.iterator(); it1.hasNext();){
                       Location loc1 = it1.next();
                       loc1.setName("");
                       list.add(loc1);
                       it1.remove();
                   }
                   if(loc.size()==0) it.remove();
                  // it.remove();
            }  
//            if(m.keySet().size()>0) throw new RuntimeException("!!");
            m.put("", list);
            Collections.sort(list);
       }
        keys.clear();
        keys.add("");
       //this.sort();
    }
public void merge(int noCop, String name, double frac){
    List<Location> ss = noCop==0 ? this.deletions.get(name) : this.amplifications.get(name);
    if(ss.size()==0) return;
    Location loc = ss.get(0);
    long st = Integer.MIN_VALUE;
    while(loc!=null){
        if(loc.min<st) throw new RuntimeException("!!");
        st = loc.min;
       // System.err.println("growing loc "+loc+" "+this.l.size());
       Location nxt =  growRight(loc, name, frac);
       loc = nxt;
     
    }
}

public void removeWithObsLessThan(int thresh){
    for(Iterator<Location> it = this.iterator(); it.hasNext();){
        if(it.next().noObs.size()<thresh) it.remove();
    }
}
/* frac is fraction of minimum */
public Location growRight(Location loc, String name, double frac){
    if(loc==null) throw new RuntimeException("!!");
   List<Location> abs = loc.noCop()==0 ? deletions.get(name) : amplifications.get(name);
   // if(!this.deletions.contains(loc)) throw new RuntimeException("!!");
   int index = abs.indexOf(loc);
   Iterator<Location>  it =abs.subList(index+1, abs.size()).iterator();
   Location loc1 = null;
  // it.next();
   while(it.hasNext()){
      loc1 = it.next();
      double overl = loc.overlaps(loc1);
       if(overl>=0){
           if(overl >= frac * Math.min(loc.size(), loc1.size())){
         //  System.err.println("merging "+loc+" with "+loc1);
               loc.min = Math.min(loc.min, loc1.min);
               loc.max = Math.max(loc.max, loc1.max);
               loc.incrObs(loc1);
               it.remove();
           }
       }
       else return abs.get(index+1);
   }
   return null;
}


/** merges with overlaps with start greater than or equal */
    public void add(Location loc){
        if(loc.size() > this.lengthLimit
             
                ){
           throw new RuntimeException(""+loc);
        }
        this.keys.add(loc.name());
        Map<String, List<Location>> m1 = loc.noCop()==0 ? deletions : amplifications;
        if(!m1.containsKey(loc.name())){
            deletions.put(loc.name(), new ArrayList<Location>());
           amplifications.put(loc.name(), new ArrayList<Location>());
        }
        List<Location> abs = m1.get(loc.name());
      
        abs.add(loc);
    } 

    
    public Iterator<Location> iterator() {
        final Iterator<String> it = this.keys.iterator();
        return new Iterator<Location>(){
            Iterator<Location> current;
            Iterator<Location> prev;
            {
                getNextIt();
                prev = current;
            }
            private void getNextIt(){
                while(it.hasNext() && (current==null || !current.hasNext())){
                    current= iterator(it.next());
                }
            }
            public boolean hasNext() {
                return current!=null && current.hasNext();
            }

            public Location next() {
               Location res = current.next();
               prev = current;
               if(!current.hasNext()) getNextIt();
               return res;
            }

            public void remove() {
               prev.remove();
            }
            
        };
    }
    
   
   public Iterator<Location> iterator(String name) {
       Collection<Location> del = this.deletions.get(name);
       if( del==null) del = Arrays.asList(new Location[0]);
       Collection<Location> ampl = this.amplifications.get(name);
       if( ampl==null) ampl = Arrays.asList(new Location[0]);
        final Iterator<Location > it1 = del.iterator();
        final Iterator<Location> it2 =ampl.iterator();
      return new Iterator<Location>(){
            boolean use1 = true;
            public boolean hasNext() {
               return it1.hasNext() || it2.hasNext();
            }
    
            public Location next() {
              if(!use1) return it2.next();
              else if(it1.hasNext()) return it1.next();
              else {
                  use1 = false;
                  return it2.next();
              }
            }
    
            public void remove() {
                if(use1) it1.remove();
                else it2.remove();
            }
          
      };
    }
   public List<Location> get(String name, int noCop){
       return   noCop==0 ? this.deletions.get(name) : this.amplifications.get(name);
   }
   public Iterator<Location> iterator(String name, int noCop) {
       Collection<Location> del = get(name, noCop);
         
        if(del==null) return Arrays.asList(new Location[0]).iterator();
        else return del.iterator();
    }
    public void thin(int thres){
        for(Iterator<Location> it = iterator(); it.hasNext();){
            Location nxt = it.next();
            if(nxt.noObs.size()<thres) it.remove();
        }
    }
 
    public List<Integer> getLocs() {
        List<Integer> l = new ArrayList<Integer>();
        for(Iterator<Location> it = this.iterator(); it.hasNext();){
            Location lo = it.next();
            if(lo.min!=lo.max) throw new RuntimeException("!!");
            l.add((int)lo.min);
        }
        return l;
    }
    public Location  overlaps(Location li, int thresh) {
      for(Iterator<Location> it = this.iterator(); it.hasNext(); ){
          Location loc = it.next();
          if(loc.overlaps(li)>thresh) return loc;
      }
      return null;
      
    }

    
    public Location  contains(int pos, int thresh) {
        for(Iterator<Location> it = this.iterator(); it.hasNext(); ){
            Location loc = it.next();
            if(loc.overlaps(pos)>thresh) {
            	return loc;
            }
        }
        return null;
        
      }

    public void setChr(String chr) {
        for(Iterator<Location> it = iterator(); it.hasNext();){
            Location loc = it.next();
            loc.chr = chr;
        }
    }


    public void restrict() {
       for(Iterator<Location> it = this.iterator(); it.hasNext();){
           Location loc = it.next();
           SortedSet<Integer> ss = this.noPosIn(loc);
           if(ss.size()>0){
           loc.min = ss.first();
           loc.max = ss.last();
           }
       }
        
    }
}
