/**
 * 
 */
package lc1.dp.data.collection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.SimpleExtendedDistribution;

public class MergedPhenotypes extends Phenotypes{
    Phenotypes[] phens;
    
    
    static Comparator cc = new Comparator<Map.Entry<String, Integer>>(){

        public int compare(Entry<String, Integer> o1, Entry<String, Integer> o2) {
           return o1.getValue().compareTo(o2.getValue());
        }
        
    };
    class Info{
        public Info(String trait){
            this.trait = trait;
            this.phenAlias = new int[phens.length];
            Arrays.fill(phenAlias, -1);
            this.vals=new Map[phens.length];
            this.dists = new ProbabilityDistribution[phens.length];
            recoding = new Map[phens.length];
        }
        /** ie we find this phenotype in the jth position on the ith dataset*/
        public void setAlias(int dataset_index, int pos_ind_dataset){
            phenAlias[dataset_index] = pos_ind_dataset;
            this.dists[dataset_index] = phens[dataset_index].phenotypeDistribution[pos_ind_dataset];
            this.vals[dataset_index] = phens[dataset_index].phenVals[pos_ind_dataset];
        }
        String trait;
        int[] phenAlias;
        Map<String, Integer>[] vals;
        ProbabilityDistribution[] dists;
        
        Map<Integer, Integer>[] recoding;
        
        //with all datasets considered
        Map<String, Integer> coding =null;//
        ProbabilityDistribution dist = null;
        
        public void calcRecoding(){
          
            int counter = 0;
            for(int i=0; i<this.vals.length; i++){
                if(phenAlias[i]>=0){
                    if(vals[i]!=null){
                        if(dist!=null) throw new RuntimeException("!!");
                        if(coding==null) coding =  new HashMap<String, Integer>();
                        recoding[i] = new HashMap<Integer, Integer>();
                        List<Map.Entry<String, Integer>> l = new ArrayList<Map.Entry<String, Integer>>(vals[i].entrySet());
                        Collections.sort(l, cc);
                        for(Iterator<Map.Entry<String, Integer>> it = l.iterator(); it.hasNext();){
                            Map.Entry<String, Integer> nxt = it.next();
                            Integer ind = coding.get(nxt.getKey());
                            if(ind==null){
                                coding.put(nxt.getKey(), ind  = counter);
                                counter++;
                            }
                            recoding[i].put(nxt.getValue(), ind);
                        }
                    }
                    else{
                        if(coding!=null) throw new RuntimeException("!!");
                        if(dist==null){
                            dist = dists[i].clone();
                          //  dist.initialise();
                        }
                        else{
                            dist.addCounts(dists[i]);
                        }
                        
                    }
                }
            }
            if(dist==null){
                dist = new SimpleExtendedDistribution(coding.keySet().size());
                for(int i=0; i<this.vals.length; i++){
                    if(phenAlias[i]>=0){
                        ((SimpleExtendedDistribution)dist).addCounts((SimpleExtendedDistribution)dists[i], recoding[i]);
                    }
                }
            }
            dist.transfer(0.0);
            
        }
        public Double transform(int i, Double phenValue) {
          if(this.recoding[i]==null|| phenValue==null) return phenValue;
          else{
              
              return (double) recoding[i].get(phenValue.intValue());
          }
        }
    }
    Map<String, Info> traits = new HashMap<String, Info>();
    
   
  MergedPhenotypes(DataCollection[] ldl){
      this.phens = new Phenotypes[ldl.length];
      for(int i=0; i<phens.length; i++){
          phens[i] = ldl[i]==null ? null : ldl[i].pheno;
      }
       getPhenAlias();
       this.phen = new ArrayList<String>(traits.keySet());
       this.phenotypeDistribution = new ProbabilityDistribution[phen.size()];
       this.phenVals = new Map[phen.size()];
       for(int i=0; i<phen.size(); i++){
           String pheno = phen.get(i);
           Info inf = traits.get(pheno);
           inf.calcRecoding();
           phenotypeDistribution[i] = inf.dist;
           phenVals[i] = inf.coding;
       }
       
  }
    
    private void getPhenAlias() {
      for(int i=0; i<phens.length; i++){
          if(phens[i]!=null){
              for(int j=0; j<phens[i].size(); j++){
                  String ph = phens[i].phen.get(j);
                  Info res = traits.get(ph);
                  if(res==null){
                      traits.put(ph, res = new Info(ph));
                      
                  }
                  res.setAlias(i, j);
              }
          }
      }
    }

    public Double[] getRecodedValues(HaplotypeEmissionState[] tmp) {
       Double[] phenD = new Double[this.phen.size()];
       for(int k=0; k<phen.size(); k++){
           Info inf = this.traits.get(phen.get(k));
           for(int i=0; i<tmp.length; i++){
               int alias = inf.phenAlias[i];
               if(tmp[i]!=null &&  alias>=0){
                   Double sc = inf.transform(i, tmp[i].phenValue()[alias]);
                   if(phenD[k]!=null){
                       if(Math.abs(phenD[k] - sc)>0.001) throw new RuntimeException("!!");
                   }
                   else phenD[k] = sc;
               }
           }
           
       }
       return phenD;
       
        
    }
   
}