/**
 * 
 */
package lc1.simulation;

import pal.tree.Node;

public  class NodeNum implements Comparable{
  public NodeNum(Node n, int num, int length){
      this.n = n;
      this.num = num;
      this.length = length;
  }
  Node n;
  int num; //position on (0,1) of mutation.
    int length; //length in bp of mutation event
    public int compareTo(Object o) {
       int d1 = ((NodeNum)o).num;
       int d = num;
       if(d==d1) return 0;
       else if(d<d1) return -1;
       else return 1;
    }
    public String toString(){
        return 
        String.format("%2s_%5i %5i", new Object[] {
                n.getIdentifier(), num, length});
    }
    public int start() {
        // TODO Auto-generated method stub
        return num;
    }
    public int end() {
        // TODO Auto-generated method stub
        return num+length;
    }
}