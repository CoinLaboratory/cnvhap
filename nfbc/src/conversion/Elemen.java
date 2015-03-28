/**
 * 
 */
package conversion;

class Elemen{
	  int start;
	  int end;
	  String id;
	  int cn;
	  public boolean contains(int pos){
		  return pos>=start-5 && pos<=end+5;
	  }
	  public Elemen(int start, int end,  int cn,String id){
		  this.start = start;
		  this.end = end;
		  this.id = id;
		  this.cn = cn;
	  }
	public String toString(){
		return start+":"+end+"_"+id+"_"+cn;
	}
	public void set(int pos) {
		start = pos+2;
		end = pos-2;
		
	}
	public int overlap(Elemen next) {
	return Math.min(next.end - start, end - next.start);
	}
  }