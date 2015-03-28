package lc1.CGH;

import java.io.BufferedReader;
import java.io.File;

public class SNPLocationReader1 extends Locreader {

  
    public SNPLocationReader1 (File dir, String chromosome, Location region)throws Exception{
      super(Long.MAX_VALUE,"");
        //  List<Integer> locs = AberationFinder.readPosInfo(new File(dir.getParentFile(), "data.txt"), 4, true);
        BufferedReader br = AberationFinder.getBufferedReader(dir, "snp.txt");
       br.readLine();
        String st;
   //    br.readLine();
   //    boolean first = false;
        while((st = br.readLine())!=null){
            if(st.startsWith("|")) st = br.readLine();
           String[] str =  st.trim().split("\\s+");
         //  for(int i=0; i<str.length; i++){
            long min = Integer.parseInt(str[2]);
           // if(first){
           //     if(!locs.contains((int)min)) throw new RuntimeException("!!!! "+chromosome);
           //     first = false;
           // }
            Location location =  new Location(chromosome,  min, min) ;
         //  System.err.println(location);
            if(region==null || region.overlaps(location)>=0){
                add(location);
            }
        }
        br.close();
        sort();
    }
}
