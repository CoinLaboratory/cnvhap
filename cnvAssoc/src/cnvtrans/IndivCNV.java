package cnvtrans;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

public class IndivCNV {
     String id;
   
	 SortedSet<CNV> cnvs;
	 int index;
	 IndivCNV(String id1, int index){
		 this.id = id1;
		 cnvs = new TreeSet<CNV>();
		 this.index = index;
	
	 }

	public IndivCNV(IndivCNV indivCNV) {
		this.id = indivCNV.id;
		this.index = indivCNV.index;
		 this.cnvs = new TreeSet<CNV>();
	}

	public void add(CNV cnv) {
	   cnvs.add(cnv);
		
	}
	
	//CNV tmp = new CNV();
	/* first row is count all, then count 0, then count ampl
	public boolean nocopies(int position, double[][] res,int base){
		Iterator<CNV> it =cnvs.iterator();
		
		while(it.hasNext()){
			CNV  cnv = it.next();
			if(cnv.contains(position)){
				res[0][cnv.type]+=cnv.cert;
				res[1][Math.max(0,base-cnv.type)]+=cnv.cert;
				res[2][Math.max(0,cnv.type-base)]+=cnv.cert;
				if(cnv.type<base){
					res[0][cnv.type+1]+=(1-cnv.cert);
					res[1][Math.max(0,base-(cnv.type+1))]+=(1-cnv.cert);
					res[2][0]+=(1-cnv.cert);
				}
				else if(cnv.type>base){
					res[0][cnv.type-1]+=(1-cnv.cert);
					res[1][0]+=(1-cnv.cert);
					res[2][Math.max(0,(cnv.type-1)-base)]+=1-cnv.cert;
				}
				return true;
			}
		}
		
			res[0][base]=1.0;
			res[1][0] =1; res[2][0] = 1;
			return false;
	} */
	
	/* this does not take into account certainty */
	public CNV nocopies(int position){
		Iterator<CNV> it =cnvs.iterator();
		
		while(it.hasNext()){
			CNV  cnv = it.next();
			if(cnv.contains(position)){
				return cnv;
			}
		}
		return null;
			
	}

	
	
	public CNV hasOverlap(CNV cnv) {
		int maxOverl = -1000000000;
		CNV max = null;
	for(Iterator<CNV> it = this.cnvs.iterator(); it.hasNext();){
		CNV cnv2 = it.next();
		int overl = cnv.overlaps(cnv2);
		int len2 = cnv2.length();
		int len1 = cnv.length();
		if(overl>maxOverl){
			max = cnv2;
			maxOverl = overl;
		
			//	System.err.println("overlap "+cnv+" "+cnv2+" "+overl+" "+len1+" "+len2);
			//if(overl > 0.5*len1 && overl > 0.5*len2) return cnv2;
			if(overl>= CalcOverlaps.overlThresh)return cnv2;
		}
	
	}
	//if(maxOverl<=0){
	//	System.err.println("h");
	//}
//	if(maxOverl>=0) return max;
	 return null;
		
	}

	public IndivCNV annontateNum(SortedSet<Integer> list,int num) {
		IndivCNV res = new IndivCNV(this);
		for(Iterator<CNV> it = this.cnvs.iterator(); it.hasNext();){
			CNV cnv = it.next();
			if(list!=null){
				int nosnp = cnv.getNoSnp(list);
			if(nosnp>=num){
				CNV cnv_new =	new CNV(cnv);
			      cnv_new.nosnp = nosnp;
				res.add(cnv_new);
			}
		}
			//else{
			//	System.err.println("h");
			//}
			
		}
		return res;
		
	}

	public int[] getVals(SortedSet<Integer> list, int v, boolean left){
		SortedSet<Integer> ts = list.tailSet(v-1);
		int v1 = v;
		int v2 = v;
		int v3 =list.tailSet(v).first(); 
		if(v3==v){
			if(left){
				 v1 = v>list.first() ? list.headSet(v-1).last(): v-500;
				 v2 = v;
			}else{
				v2 = v<list.last() ? list.tailSet(v+1).first() : v+500;
				v1 = v;
			}
			
		}else{
		v2 = v<list.last() ? list.tailSet(v-1).first() : v+500;
		 v1 = v>list.first() ? list.headSet(v+1).last(): v-500;

	
		}
		
		if(v1==v2) {
		throw new RuntimeException("!!");
		}
		return new int[] {v1,v2};
	
	}
	
	public IndivCNV getFlanking(SortedSet<Integer> list, double expand) {
		IndivCNV res = new IndivCNV(this.id,this.index);
		for(Iterator<CNV> it = this.cnvs.iterator(); it.hasNext();){
			CNV cnv = it.next();
			
			if(list!=null){
				CNV cnv_left = new CNV(cnv);
				CNV cnv_right = new CNV(cnv);
				cnv_left.end = cnv.start;
				cnv_right.start = cnv.end;
				cnv_right.nosnp = 1;
				cnv_left.nosnp = 1;
			int[] vals = getVals(list,cnv.start,true );
			int[] valsR = getVals(list,cnv.end ,false);
				cnv_left.start =(int)  Math.round(vals[0]*expand + vals[1]*(1-expand));
				SortedSet<Integer> right =list.tailSet(cnv.end+1);
				cnv_right.end = (int)  Math.round(valsR[1]*expand + valsR[0]*(1-expand));
			//	System.err.println(cnv);
			//	System.err.println(cnv_left);
			//	System.err.println(cnv_right);
				res.add(cnv_right);
				res.add(cnv_left);
				//if(expand){
					cnv.start = cnv_left.start;
					cnv.end = cnv_right.end;
					
				//}
		   }
			
			
		}
		return res;
	}
}
