/**
 * 
 */
package lc1.dp.states;

class PosTrans1 implements Comparable{
   Integer pos;
   int[] ind_loc;
   String snp_id;
   public String toString(){
       return pos+"_"+ind_loc[0]+"_"+ind_loc[1];
   }
   public void replace(PosTrans1 pt) {
   this.ind_loc = pt.ind_loc;
    
}
PosTrans1(int pos, int ind, int loc, String snp_id){
       this.pos = pos;
       this.snp_id = snp_id;
       this.ind_loc = new int[] {ind, loc};
   }
public int compareTo(Object o) {
    return  pos.compareTo(((PosTrans1)o).pos);
}
}