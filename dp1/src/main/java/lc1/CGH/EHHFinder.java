/**
 * 
 */
package lc1.CGH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;

import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.sequenced.Convert;
import lc1.util.Constants;

public class EHHFinder{
    Collection<String> set;
    String base;
    public static void main(String[] args){
        try{
            Constants.parse(args);
            EHHFinder.thresh = 0;
            EHHFinder.lim=  Integer.MAX_VALUE;
            EHHFinder.thresh1 = 0.0;
            File user = new File(System.getProperty("user.dir"));
            Location[] del = Convert.getDeletions(new File(user, "../CGH/deletion_samples.txt"), 
                    Integer.parseInt(user.getName()),
                    Constants.mid()[0]);
            Collection<String> l = del[0].noObs;
        SimpleDataCollection sdt = 
            SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader(user, "phased1.txt"),Emiss.class, Emiss.getEmissionStateSpace(1));
        sdt.removeKeyIfStartsWith("NA");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(user, "hap_summary.txt"))));
        EHHFinder ehh = new EHHFinder(sdt);
        
   int st1 = sdt.loc.indexOf(del[0].min);
   int end1 = sdt.loc.indexOf(del[0].max);
           Set<String> coreHaps =  sdt.getAllHaplotypes(st1, end1).keySet();
          for(Iterator<String> it = coreHaps.iterator(); it.hasNext();){
              String st = it.next();
              pw.println("HAPLOTYPES FOR "+st);
              ehh.base = st;
              System.err.println(st);
             if(st.startsWith("_") && st.endsWith("_")) ehh.set = l;
             else {
              //  continue;
                ehh.set = null;
             }
            ehh.extend(st, st1, end1, pw);
     
            pw.println("##############################################");
            pw.flush();
          }
//            List<String> indiv = sdt.getNames();
             
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   Comparator comp= new Comparator<String>(){
       public int compare(String o1, String o2) {
          int l1 = o1.length();
          int l2 = o2.length();
          if(l1!=l2) return l1> l2 ? -1 : 1;
          else return o1.compareTo(o2);
       }
         
     };
     static int thresh =3 ;
     static double thresh1 = 0.2;
    final SimpleDataCollection sdt;
    public EHHFinder(SimpleDataCollection sdt2) {
        this.sdt = sdt2;
    }
    private Set<String> getIndividuals(List<PIGData> value) {
        Set<String> res = new HashSet<String>();
        for(int i=0; i<value.size(); i++){
            res.add(value.get(i).getName());
        }
        return res;
      }
    int orig_size;
    public void extend(Aberation ab, PrintWriter pw){
         String[] nme = ab.name.split("_");
         PIGData haplotype = sdt.get(nme[0]).split()[Integer.parseInt(nme[1])];
          String string = haplotype.getStringRep(ab.start, ab.end);
          extend(string, ab.start, ab.end, pw);
    }
    
    public void extend(String string, int start_, int end_,  PrintWriter pw){
     
        List<PIGData> dat = sdt.getAllHaplotypes(start_, end_).get(string);
       /// check(dat);
        if(dat==null) return;
        if(set!=null){
            Set<String> hapl_names = new HashSet<String>();
            for(Iterator<PIGData> it = dat.iterator(); it.hasNext();){
                PIGData dat_ = it.next();
                String name = dat_.getName().split("_")[0].split("\\|")[0];
                if(name.startsWith("NA")) it.remove();
                
                else if(!set.contains(name)){
                 //   CompoundScorableObject dat1_ = (CompoundScorableObject) dat_.clone();
                  //  dat1_.restrictSites(start_-5, end_+5);
                    pw.println("deletion_samples.txt does not contain "+name+" "+dat_.getStringRep(start_-5, end_+5));
                   // it.remove();
                 }
                else hapl_names.add(name);
            }
            for(Iterator<String>  it = set.iterator(); it.hasNext();){
                String key = it.next();
                if(!hapl_names.contains(key)){
                    pw.println("phased deletion does not contain "+key);
                }
            }
        }
      orig_size = dat.size();
      pw.println("original size is "+orig_size);
      SortedMap<Integer, SortedMap<String, Object[]>> map = new TreeMap<Integer, SortedMap<String, Object[]>>();
   //   System.err.println("HEEEERE");
   //  check(dat);
      extend(start_, end_, string, dat, map, 0, string.length());
      pw.println("aberation "+string+" "+sdt.loc.get(start_)+" "+sdt.loc.get(end_));
   
      int cnt =0;
      for(Iterator<Integer> it = map.keySet().iterator(); it.hasNext() ;cnt++){
          Integer key1 = it.next();
          SortedMap<String, Object[]> m1 = map.get(key1);
          List<PIGData> vals = null;
          for(Iterator<String> it1 = m1.keySet().iterator(); it1.hasNext();){
          String key = it1.next();
          Object[] val = map.get(key1).get(key);
          List<PIGData> val0 = (List<PIGData>) val[0];
          if(//val0.size() < 0.5*(double)orig_size && 
                  vals!=null && val0.equals(vals)) continue;
          else vals = val0;
          int start = (Integer) val[1];
          int end = (Integer) val[2];
        //  String st = nxt.getKey();
          PIGData da = SimpleScorableObject.make("", 
        		  Arrays.asList(new String[] {key}),  
        		  Emiss.getEmissionStateSpace(1),
        		  (short)-1);
          StringWriter sw = new StringWriter();
          PrintWriter pw1 = new PrintWriter(sw);
          da.print(pw1, false, true, true, null, null,null);
          String st1 = sw.toString();
          Set<String> indv = getIndividuals(val0);
       //  if(cnt<5 || val0.size()>=4) {
             pw.println(base+" : "+da.length()+"\tr"+val0.size()+" of "+orig_size+"\tl "+(sdt.loc.get(end) - sdt.loc.get(start))+"\t"+indv);
             pw.println(st1.substring(0, st1.indexOf("\n")-1)+"|\t\t\t"+base+"\t\t\t"+((double)val0.size() / (double)orig_size));
       //  }
            /* if(key.length()>1){
                 String keyL = key.substring(1);
                 String keyR = key.substring(0, key.length()-1);
                inner: for(int key11 = key1; ; key11++){
                     SortedMap<String, Object[] > m11 = map.get(key11);
                     if(m11==null) break inner;
                     if(m11.containsKey(keyL)){
                         Set<String> valL = getIndividuals((List<PIGData>) m11.get(keyL)[0]);
                         valL.removeAll(indv); 
                         pw.println("\t\tleft lo "+valL);
                         break inner;
                     }
                }
                 inner: for(int key11 = key1; ; key11++){
                     SortedMap<String, Object[] > m11 = map.get(key11);
                     if(m11==null) break inner;
                     if(m11.containsKey(keyR)){
                         Set<String> valR = getIndividuals((List<PIGData>) m11.get(keyR)[0]);
                         valR.removeAll(indv); 
                         pw.println("\t\tright lo "+valR);
                         break inner;
                     }
                } 
                  
             }*/
          }
          //+map+"\n");
      }
    //    System.err.println(string+" "+dat.size()+" "+sdt.getKeys().size()+" "+left.getValue().size()+" "+right.getValue().size());

    }
   private void print(SortedMap<Integer, SortedMap<String, Object[]>> map, Integer key1, String key) {
        // TODO Auto-generated method stub
        
    }
/* public static void check(Collection<PIGData> dat){
        for(Iterator<PIGData> it = dat.iterator(); it.hasNext();){
            if(it.next().length()<35) throw new RuntimeException("!!");
        }
    }*/
   static Integer lim = 1;
    public void extend(int start, int end, String key, List<PIGData> dat,
            SortedMap<Integer, SortedMap<String, Object[]>>map, int base_index0, int base_index1){
         if(dat.size() < thresh || dat.size() < thresh1 * orig_size) return;
     //   else if(key.length()>10){
            SortedMap<String, Object[]> map1  = map.get(dat.size());
             if(map1==null) map.put(dat.size(), map1 = new TreeMap<String, Object[]>(comp));
             char[] ch = key.toCharArray();
             Arrays.fill(ch, base_index0, base_index1, '*');
             String key1 =new String(ch);
            if( map1.containsKey(key1)) return;
          
           // if(base_index1 - base_index0>1){
             //   throw new RuntimeException("!!");
          //  }
            map1.put(key1, new Object[] {dat, start, end});
       // }
        System.err.println(key+" "+dat.size());
        if(end < dat.get(0).length()-1){
            int k=0;
            for(Iterator<Entry<String, List<PIGData>>> it =
                sdt.getHaplotypes(sdt.getHaplotypes(start,end+1, dat)).iterator();
            it.hasNext() && k<lim;k++){
               Entry<String, List<PIGData>> right= it.next();
               extend(start, end+1, right.getKey(), right.getValue(), map, base_index0, base_index1);
            }
        }
        if(start-1>0){
            int k=0;
            for(Iterator<Entry<String, List<PIGData>>> it =sdt.getHaplotypes(sdt.getHaplotypes(start-1,end, dat)).iterator();
            it.hasNext()&& k<lim;k++){
            Entry<String, List<PIGData>> left = it.next();
            extend(start-1, end, left.getKey(), left.getValue(), map, base_index0+1, base_index1+1);
            }
           
        }
    
    }
}