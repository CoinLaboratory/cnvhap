/**
 * 
 */
package summary;

class Entry implements Comparable{
    double sc;
   // double sc1;
  //  double sc2;
 //  String  maf;
    String name;
    int loc;
    
    Entry(String st){
         String[] str = st.trim().split("\\t");//
         String[] str1 =str[3].split("/"); 
        // maf = str
         sc = Double.parseDouble(str1[0]);
     //    sc1 = -1*Math.log10(Double.parseDouble(str1[1]));
     //    sc2 = -1*Math.log10(Double.parseDouble(str1[2]));
         name = str[0];
         loc = Integer.parseInt(str[1].trim());
    }
    public String toString(){
        return name+"\t"+loc+"\t"+Math.log(sc);//+"\t"+sc1+"\t"+sc2;
    }
    public int compareTo(Object o) {
       Entry e1 = (Entry) o;
       if(e1.sc !=sc) return sc< e1.sc ? -1:1;
       else return name.compareTo(e1.name);
    }
}