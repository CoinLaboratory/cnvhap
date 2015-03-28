package lc1.CGH;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

public abstract class AbstractAberatiionReader extends Locreader {
    public AbstractAberatiionReader(long lengthLim, String name) {
        super(lengthLim, name);
        // TODO Auto-generated constructor stub
    }
    protected int[] col; 

    public static Iterator<String[]> getIterator(final BufferedReader br) throws Exception{
        return new Iterator<String[] >(){
            String[] nxt = br.readLine().split("\\t");
            public boolean hasNext() {
              return nxt!=null && nxt.length>0 && (nxt.length>1 || nxt[0].length()>1) ;
            }
            public String[] next() {
               String[] res = nxt;
               try{
               String str = br.readLine();
               if(str!=null){
                   nxt = str.trim().split("\\t");
               }
               else {
                   nxt = null;
               }
               }catch(Exception exc){ exc.printStackTrace(); System.exit(0);}
               return res;
            }
            public void remove() {}
        };
    }
    public abstract String getName(String[] str);
    public abstract String getChr(String[] str);
    public abstract String getStart(String[] str);
    public abstract String getEnd(String[] str);
    public abstract int getNoProbes(String[] str);
    public  String getProbeId(String[] str){
        return "";
    }
    public abstract double getNoCopy(String[] str);
   
    public boolean exclude(String[] str){
        return Math.abs(getNoCopy(str))<AberationFinder.lowIntens;
    }
   
    public abstract String[]  getCols();
    
    public void initialise( BufferedReader br,String chromosome, Location region,int first, int second,  String name, Set<String > indiv)throws Exception{
        initialise(getIterator(br), chromosome, region, first,second, name, indiv);
    }
    public void initialise( File dir,String chromosome, Location region, int first, int second, String name,
    Set<String> indiv       
    )throws Exception{
        File f = new File(dir, name);
        BufferedReader br;
        if(f.exists()){
            br = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
        }
        else{
            
            f = new File(dir, name+".gz");
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
        }
        initialise(br, chromosome, region, first,second, name, indiv);
        br.close();
    }
    public void initialise( Iterator<String[]> it ,String chromosome, Location region,   int firstCount, int secondCount, String name,
            Set<String> indiv)throws Exception{
     
        for(int i=0; i<firstCount; i++){
            
            System.err.println("skipping "+Arrays.asList(it.next()));
        }
         fillMap(it.next(), getCols());   
        for(int i=0; i<secondCount; i++){
            System.err.println("skipping "+Arrays.asList(it.next()));
        }
        while(it.hasNext()){
            String[] str = it.next();
            
            String chrom = getChr(str);
            if(chromosome!=null && !chrom.equals(chromosome)){
                continue;
            }
            if(exclude(str)){
                Logger.global.info("excluding "+Arrays.asList(str));
                 continue;
             }
            long min = Long.parseLong(getStart(str));//Integer.parseInt(pos[0]);
            long max = Long.parseLong(getEnd(str));//Integer.parseInt(pos[1]);
           
           double noCop = 
              getNoCopy(str);
           
            Location location = max > min ? new Location(chrom, min,max) :new Location(chrom, max,min) ;
           if(indiv!=null && !indiv.contains((getName(str)))) continue;;
            location.setName(getName(str));
            if(noCop==0) location.setNoCop(-1);
            else location.setNoCop(noCop>0 ? 2 : 0);
            location.setLogL(noCop);
            location.setNoSnps(getNoProbes(str));
            location.setProbeId(getProbeId(str));
            if(region==null || region.overlaps(location)>0){
                add(location);
            }
            else{
            //    Logger.global.info("excluded "+location+" does not match "+region);
            }
        }
        sort();
    }

    protected void fillMap(String[] header, String[] cols) {
        col = new int[cols.length];
            Arrays.fill(col, -1);
          for(int j=0; j<cols.length; j++){
              inner: for(int i=0; i<header.length; i++){
                  if(header[i].startsWith(cols[j])){
                      col[j] = i;
                      break inner;
                  }
              }
          }
        
      }
    public void addAll(Collection<Location> probesinMultiRegion) {
       for(Iterator<Location> it = probesinMultiRegion.iterator(); it.hasNext();){
           this.add(it.next());
       }
        
    }
    public void restrictEnds(Locreader p244k) {
        for(Iterator<Location> it = this.iterator(); it.hasNext();){
            Location loc = it.next();
           // System.err.println("bef "+loc);
            List<Location > l = loc.getStartOverlaps(p244k.iterator(), 0);
            Collections.sort(l);
            if(l.size()==0) throw new RuntimeException("!!");
            else{
                loc.min = l.get(0).min;
                loc.max = l.get(l.size()-1).max;
            }
         //   System.err.println("aft "+loc);
        }
        
    }
    public void getNoProbes(Locreader locsL) {
        for(Iterator<Location> it = this.iterator(); it.hasNext();){
            Location loc = it.next();
            List<Location > l = loc.getStartOverlaps(locsL.iterator(), 0);
            
            loc.noSnps = l.size();
            if(l.size()==0){
                System.err.println("removed "+loc);
                it.remove();
            }
        }
        
    }

}





