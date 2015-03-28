package lc1.sequenced;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.CGH.Location;
import lc1.dp.data.collection.LikelihoodDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.SSOData;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.SimpleExtendedDistribution;

public class Convert {

    public static void main(String[] args){
        try{
       File f = new File("deletion_samples.txt");
    Location[] l = getDeletions(f, 6,  null);
     LikelihoodDataCollection ldl = getLikelihoodDataCollection(l);
     EmissionState emSt = ldl.dataL.get("21938");
       for(int i=0; i<l.length; i++){
           System.err.println(l[i]);
       }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
 
    
    public static Location[] getDeletions(File f, Integer chr, int[] mid) throws Exception{
        Location midLoc = 
            mid!=null ? 
            new Location(""+ chr, (long) mid[0], (long) mid[1]) : null;
        BufferedReader br = new BufferedReader(new FileReader(f));
        String st = br.readLine();
        String[] str = st.split("\t");
        Location[] l = new Location[(int) Math.ceil((double) str.length/2.0)] ;
        for(int i=0; i<l.length; i++){
            l[i] = new Location(str[i*2]);
            
        }
        while((st = br.readLine())!=null){
           str = st.split("\t");
           for(int i=0; i<l.length; i++){
               if(str.length>i*2){
                   if(str[i*2].length()>0){
                   l[i].noObs.add(str[i*2]+"_"+ str[i*2+1].trim());
                   }
               }
           }
        }
        List<Location> res = new ArrayList<Location>();
        for(int i=0; i<l.length; i++){
            double  overl = midLoc==null  ? 10 : l[i].overlaps(midLoc);
            System.err.println("overl "+midLoc+" "+l[i]+" "+overl+" "+chr);
            if(overl>0){
               // if(mid==null || l[i].start<= mid[0]){
                 //   if(mid==null || l[i].end >=mid[0]){
                        res.add(l[i]);
                    }
                    else{
                        Logger.global.info("excluded "+l[i]+" cf "+mid);
                    }
               
        }
        Collections.sort(res);
        return res.toArray(new Location[0]);
    }
   
    public static LikelihoodDataCollection getLikelihoodDataCollection(Location[] d){
        if(d.length==0) throw new RuntimeException("zero length");
        Set<String> keys = new HashSet<String>();
        List<Integer> locs = new ArrayList<Integer>();
        for(int i=0; i<d.length; i++){
            
            keys.addAll(d[i].noObs);
            locs.add((int)d[i].min);
            locs.add((int) d[i].max);
        }
      
   //   keys.add("22300_homoA");
        LikelihoodDataCollection ldl = new LikelihoodDataCollection(locs);
   //     keys.clear();
        EmissionStateSpace emStSp = Emiss.getEmissionStateSpace(1);
        if(emStSp.defaultList.size()<3) throw new RuntimeException("!!");
        for(Iterator<String> it = keys.iterator();it.hasNext(); ){
            String[] key = it.next().split("_");
            SSOData sso =null;// new HaplotypeData(key[0],d.length);
           
            HaplotypeEmissionState emSt = new HaplotypeEmissionState(key[0], 2*d.length, emStSp.size(), emStSp, -1, null);
          
            for(int j=0; j<d.length; j++){
                Boolean homo = key[1].equals("homo");
                    double[] prob =new double[emStSp.defaultList.size()];
                    Arrays.fill(prob, 0);
                    for(int k=0; k<prob.length; k++){
                        int cn = emStSp.getCN(k);
                        if(key[1].equals("homo")){
                            if(cn==0) prob[k]=1.0;
                        }
                        else if(key[1].equals("hetero")){
                            if(cn==1) prob[k]=1.0;
                        }
                        else if(key[1].equals("homoA")){
                            if(cn==4) prob[k] = 1.0;
                        }
                        else if(key[1].equals("heteroA")){
                            if(cn==3) prob[k] = 1.0;
                        }
                        else throw new RuntimeException("!!"+key[1]);
                    }
                    SimpleExtendedDistribution.normalise(prob);
                   emSt.setTheta(prob, j*2);
                   emSt.setTheta(prob, j*2+1);
              
            }
          //  HaplotypeEmissionState emSt = EmissionState.getEmissionState(sso, Emiss.getEmissionStateSpace(1, 0), 0.0);
           // emSt.train_j = false;
            ldl.dataL.put(key[0], emSt);
        }
       ldl.calculateMLGenotypeData(true);
     
       if(ldl.snpid==null){
           ldl.snpid = new ArrayList();
       for(int i=0; i<ldl.loc.size(); i++){
           ldl.snpid.add("del"+i);
           
       }
       }
        return ldl;
        
        
        
    }
    
    
}
