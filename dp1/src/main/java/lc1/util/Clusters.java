/**
 * 
 */
package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

public class Clusters{
   int noStates;
   
   public void compare(Clusters clusters){
      double[] overlap = new double[] {0,0};
      int numpos = this.getNumPos();
      if(clusters.getNumPos()!=numpos) throw new RuntimeException("!!");
      
       for(int i=0;i<=numpos ; i++){
           SortedSet<Group> g1 = this.getClustersAt(i);
           compare(clusters.getClustersAt(i), g1, overlap);
           if(g1.size()==0) break;
       }
       System.err.println("total : overlap "+overlap[0] / (double)(numpos+1)+" : "+overlap[1]/ (double)(numpos+1));
   }
   
   private int getNumPos() {
       int max = 0;
       for(Iterator<Group> it = this.clusters.iterator(); it.hasNext();){
           List<Integer> pos = it.next().positions;
           int max_i = pos.get(pos.size()-1);
           if(max_i > max){
               max = max_i;
           }
       }
       return max;
}

public static void main(String[] args){
       try{
          Clusters cl1 = new Clusters(new File("cluster_output6_true.txt"));
          Clusters cl2 = new Clusters(new File("cluster_output6_false.txt"));
          System.err.println("train emissions as background");
           cl1.compare(cl2);
           System.err.println("fixed emissions as background");
           cl2.compare(cl1);
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
   
   public Clusters(File f) throws Exception{
       BufferedReader br = new BufferedReader(new FileReader(f));
       String st = "";
       while((st = br.readLine())!=null){
           String[] str = st.split(":");
           String[] membr = str[1].trim().substring(1, str[1].length()-2).split(",");
           Group g = new Group();
           for(int i=0; i<membr.length; i++){
               g.add(Integer.parseInt(membr[i].trim()));
           }
           String[] pos = str[0].trim().substring(0, str[0].length()-1).split(",");
           for(int i=0; i<pos.length; i++){
               String[] pos_i = pos[i].split("-");
               int start = Integer.parseInt(pos_i[0]);
               int end = Integer.parseInt(pos_i[1]);
               for(int j=start; j<=end; j++){
                   g.positions.add(j);
               }
           }
           clusters.add(g);
           
       }
       br.close();
   }
   
   private void compare(SortedSet<Group> clustersAt, SortedSet<Group> c1, double[] overlap) {
       Iterator<Group> it1 = clustersAt.iterator(); 
       Iterator<Group> it2 = c1.iterator();
       while(it2.hasNext()){
           Group g1 = it1.hasNext() ? it1.next() : null;
           Group g2 = it2.next();
           overlap[0]+=g2.size();
           overlap[1] += g2.overlap(g1);
       }
       
   }

public SortedSet<Group> getClustersAt(Integer i){
       SortedSet<Group> sub_clust = new TreeSet<Group>(clustcomp);
    //   int sum=0;
       for(Iterator<Group> it = clusters.iterator(); it.hasNext();){
           Group g = it.next();
           if(g.positions.contains(i)){
               sub_clust.add(g);
              // sum+=g.size();
           }
       }
      // System.err.println(sum);
       return sub_clust;
   }
   SortedSet<Group> clusters = new TreeSet<Group>(clustcomp);
   //Integer[] startPos;
    public Clusters(ArrayList<Integer[][]> res, int noStates){
        this.noStates = noStates;
        int len = res.get(0)[0].length;
        for(int i=0; i<len; i++){
            Cluster clust = new Cluster(res, i);
            for(Iterator<Group> it = clust.values.iterator(); it.hasNext();){
                Group group = it.next();
                if(group.size()>0) {
                    boolean contains = clusters.contains(group);
                    if(contains){
                        SortedSet<Group> tail = clusters.tailSet(group);
                        Group first = tail.first();
                        if(!first.equals(group)) throw new RuntimeException("!!");
                        first.positions.add(i);
                    }
                    else{
                        clusters.add(group);
                    }
                }
            }
        }
    }
    
   static  class ClustComp implements Comparator<Group> {
        public int compare(Group l1,Group l2) {
          int minSize = (int) Math.min(l1.size(), l2.size());
          for(int i=0;i<minSize; i++){
              int comp = l1.get(i).compareTo(l2.get(i));
              if(comp!=0) return comp;
          }
          if(l1.size()==l2.size()) return 0;
          else return l1.size() < l2.size() ? -1 : 1;
        }
    }
   static  Comparator clustcomp = new ClustComp();
    
    public void print(PrintWriter pw){
        for(Iterator<Group> it = this.clusters.iterator(); it.hasNext();){
            it.next().print(pw);
        }
    }
    
    class Group{
        List<Integer> positions = new ArrayList<Integer>();
        List<Integer> members = new ArrayList<Integer>();
        Group(){
            
        }
        Group(int pos){
            positions.add(pos);
        }
        public int overlap(Group g1) {
            if(g1==null) return 0;
           List<Integer> members1 = new ArrayList<Integer>(this.members);
           members1.retainAll(g1.members);
           return members1.size();
        }
        public void print(PrintWriter pw){
            StringBuffer pos = new StringBuffer();
            for(int i=0; i<positions.size(); i++){
                if(i==0 || positions.get(i)-1 != positions.get(i-1)){
                    pos.append(positions.get(i)+"-");
                }
                 if(i==positions.size()-1 || positions.get(i)+1!=positions.get(i+1)){
                    pos.append(positions.get(i)+",");
                }
            }
            pw.println(pos.toString()+": "+members);
           // pw.println(this.positions);
        }
        public String toString(){
            return members.toString();
        }
        public Integer get(int i) {
            return members.get(i);
        }
        public int size() {
         return members.size();
        }
      
        public boolean equals(Object obj){
            return members.equals(((Group) obj).members);
        }
        public void add(int k) {
           members.add(k);
            
        }
    }
    
    class Cluster{
        /** maps the cluster number to a membership list of inidividuals */
        List<Group> values;
       
        
        /*start inclusive, end exclusive */
        Cluster(ArrayList<Integer[][]> res, int pos){
            SortedMap<Integer, Group> membership = new TreeMap<Integer, Group>();
            for(int j=1; j<noStates; j++){
                    membership.put(j,  new Group(pos));
            }
            for(int k=0; k<res.size(); k++){
                Integer[][] vit = res.get(k);
                for(int j=0; j<vit.length; j++){
                        Group cluster = membership.get(vit[j][pos]);
                        cluster.add(k);
                }
            }
            values = new ArrayList<Group>(membership.values());
            Collections.sort(values, clustcomp);
        }
    }
}