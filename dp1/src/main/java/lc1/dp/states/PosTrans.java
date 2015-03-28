/**
 * 
 */
package lc1.dp.states;

class PosTrans extends PosTrans1{
  
PosTrans(int pos, int ind, int loc, String snp_id){
     super(pos, ind, loc, snp_id);
   }
public int compareTo(Object o) {
     int possc =  pos.compareTo(((PosTrans)o).pos);
     if(possc==0) return snp_id.compareTo(((PosTrans)o).snp_id);
     else return possc;
}
}