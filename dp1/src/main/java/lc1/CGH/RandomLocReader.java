package lc1.CGH;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;

import lc1.util.Constants;

public class RandomLocReader extends Locreader {
    List<String> keys;
    
    boolean randomLoc;
  public RandomLocReader(  Locreader data, List<Integer> locs)
  {
      super(Long.MAX_VALUE,"");
      randomLoc = AberationFinder.randomLoc;
      Logger.global.info("making random");
      this.keys = new ArrayList<String>(data.keys);
      for(Iterator<String> it1 = data.keys.iterator(); it1.hasNext();){
          String key = it1.next();
      for(Iterator<Location> it = data.iterator(key); it.hasNext();){
          Location loc = it.next();
        //  Location loca = null;
      //    while(loca==null){
          Location    loca = 
              randomLoc ?
                      getRandom1(loc.chr, locs,  loc):      
              getRandom(loc.chr, locs,  loc);
        //  }

          this.add(loca);
      }
      }
     
      this.sort();
      Logger.global.info("done making random");
  }
static Random rand = Constants.rand;
private Location getRandom(String chr,  List<Integer> locs, Location templ) {
    int noCop = templ.noCop();
    int len = (int) templ.size();
    int[] len1 = getStartEnd(templ, locs);
 //  int len= len1[0];
   // if(len<1) throw new RuntimeException(" "+templ.noSnps);
    int min =  0; //locs.get(0);
    int max = locs.size();//locs.get(locs.size()-1) - (int)len;
    if(max<min) throw new RuntimeException("!!");
   int r = rand.nextInt(max-min);
   int st = locs.get(r+min) - len1[1]; //start index
   int end = st +len-1; //end_index
 int noSnps =0;
 for(int i=0; i<locs.size(); i++){
     if(locs.get(i) >=st && locs.get(i)<=end) noSnps++;
 }
           Location loc = new Location(chr, st,end);
           loc.setNoCop(noCop);
//                  (locs.get(i)-st==len || rand.nextBoolean() || i==0)? locs.get(i) : locs.get(i-1));
           loc.setName(this.keys.get(rand.nextInt(this.keys.size())));
           loc.setNoSnps(noSnps);
           return loc;
//   return null;
}

private int[] getStartEnd(Location templ, List<Integer> locs) {
    int cnt =0;
    int top=0;
    int bottom=0;
    boolean first = true;
   for(int i=0; i<locs.size(); i++){
      
       if(locs.get(i) >= templ.min){
           if(first){
               bottom = locs.get(i)-(int)templ.min  ;
               first = false;
           }
           cnt++;
       }
       if(i < locs.size()-1 && locs.get(i+1)>templ.max){
           top = (int) templ.max - locs.get(i);
           break;
       }
    
   }
   return new int[] {cnt, bottom, top};
}

private Location getRandom1(String chr,  List<Integer> locs, Location templ) {
    int noCop = templ.noCop();
    int len = (int) templ.size();
    int min =locs.get(0);
    int max = locs.get(locs.size()-1) - (int)len;
   int r = rand.nextInt(max-min);
   int st = r+min; //start index
   int end = st +len-1; //end_index
 
           Location loc = new Location(chr, st,end);
           loc.setNoCop(noCop);
//                  (locs.get(i)-st==len || rand.nextBoolean() || i==0)? locs.get(i) : locs.get(i-1));
           loc.setName(this.keys.get(rand.nextInt(this.keys.size())));
           return loc;
//   return null;
}

}
