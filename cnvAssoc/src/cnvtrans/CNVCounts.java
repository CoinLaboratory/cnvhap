<<<<<<< .mine
package cnvtrans;

import java.util.Set;
import java.util.SortedSet;

public abstract class CNVCounts {
  public abstract Integer fixedPos();
  public abstract int get(int i);  
  public abstract void addto(SortedSet<Integer> cnvTypes) ;
  
  public static class IntC extends CNVCounts{
	  Integer type;
	  final static int cnt=1;
	  double cert;

	public IntC(int i, double d) {
		this.type = i;
		//this.cnt = 1;
		cert = d;
	}

	public double cert(int i){
		if(i==type) return cert;
		else return Double.NaN;
	}
	
	@Override
	public Integer fixedPos() {
		return type;
	}

	@Override
	public int get(int i) {
		if(i==type) return cnt;
		else return 0;
	}

	@Override
	public void addto(SortedSet<Integer> cnvTypes) {
	cnvTypes.add(type);
		
	}

	@Override
	public void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
			int[] certCount) {
		cnvCounts[type]+=cnt;
		sb.append(type);
		certCount[type]+=cnt;
		avgCert[type]+=cert*cnt;
	}

	@Override
	public int max() {
	return type;
	}

	@Override
	public boolean containedIn(Set<Integer> cnvTypesToInclude) {
	return cnvTypesToInclude.contains(this.type);
	}

	@Override
	public boolean passCertThresh(double cert_thresh) {
		return this.cert>=cert_thresh;
	}
  }
  public static class MultC extends CNVCounts{
	 int[] cnts;
	 double[] certs;
	 
	public MultC(String typeStr, String certStr) {
		String[] str= typeStr.split(":");
		String[] strC= certStr.split(":");
		cnts = new int[str.length];
		certs = new double[str.length];
		for(int k=0; k<cnts.length; k++){
			cnts[k] = Integer.parseInt(str[k]);
		//	System.err.println("len"+strC.length);
		//	System.err.println(strC[0]);
		//	System.err.println(k);
			if(k>=strC.length) certs[k] = Double.NaN;
			else certs[k] = (strC[k].length()==0 )? 
					Double.NaN : 
						Double.parseDouble(strC[k]);
		}
	}
	@Override
	public boolean containedIn(Set<Integer> cnvTypesToInclude) {
	for(int i=0; i<cnts.length; i++){
		if(cnts[i]>0 && cnvTypesToInclude.contains(i)) return true;
	}
	return false;
	}

	@Override
	public Integer fixedPos() {
		return null;
	}

	@Override
	public int get(int i) {
		return i<cnts.length ? cnts[i] : 0;
	}
	public double cert(int i) {
		return i<cnts.length ? certs[i] : 0;
	}
	@Override
	public void addto(SortedSet<Integer> cnvTypes) {
		for(int i=0; i<cnts.length; i++){
			if(cnts[i]>0) cnvTypes.add(i);
		}
		
	}

	@Override
	public void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
			int[] certCount) {
		for(int i=0; i<this.cnts.length; i++){
			
			sb.append(cnts[i]);
			if(i<cnts.length-1) sb.append(":");
			
			if(cnts[i]>0){
				cnvCounts[i]+=cnts[i];
				certCount[i]+=cnts[i];
				avgCert[i]+=certs[i]*cnts[i];
			}
		}
	}

	@Override
	public int max() {
		return cnts.length;
	}
	@Override
	public boolean passCertThresh(double cert_thresh) {
		// TODO Auto-generated method stub
		return true;
	}
  }
public abstract void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
		int[] certCount);
	// TODO Auto-generated method stub
public abstract int max() ;
public abstract  boolean containedIn(Set<Integer> cnvTypesToInclude) ;
public abstract boolean passCertThresh(double cert_thresh);
	


}


=======
package cnvtrans;

import java.util.Set;
import java.util.SortedSet;

public abstract class CNVCounts {
  public abstract Integer fixedPos();
  public abstract int get(int i);  
  public abstract void addto(SortedSet<Integer> cnvTypes) ;
  
  public static class IntC extends CNVCounts{
	  Integer type;
	  final static int cnt=1;
	  double cert;

	public IntC(int i, double d) {
		this.type = i;
		//this.cnt = 1;
		cert = d;
	}

	public double cert(int i){
		if(i==type) return cert;
		else return Double.NaN;
	}
	
	@Override
	public Integer fixedPos() {
		return type;
	}

	@Override
	public int get(int i) {
		if(i==type) return cnt;
		else return 0;
	}

	@Override
	public void addto(SortedSet<Integer> cnvTypes) {
	cnvTypes.add(type);
		
	}

	@Override
	public void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
			int[] certCount) {
		cnvCounts[type]+=cnt;
		sb.append(type);
		certCount[type]+=cnt;
		avgCert[type]+=cert*cnt;
	}

	@Override
	public int max() {
	return type;
	}

	@Override
	public boolean containedIn(Set<Integer> cnvTypesToInclude) {
	return cnvTypesToInclude.contains(this.type);
	}
  }
  public static class MultC extends CNVCounts{
	 int[] cnts;
	 double[] certs;
	 
	public MultC(String typeStr, String certStr) {
		String[] str= typeStr.split(":");
		String[] strC= certStr.split(":");
		cnts = new int[str.length];
		certs = new double[str.length];
		for(int k=0; k<cnts.length; k++){
			cnts[k] = Integer.parseInt(str[k]);
		//	System.err.println("len"+strC.length);
		//	System.err.println(strC[0]);
		//	System.err.println(k);
			if(k>=strC.length) certs[k] = Double.NaN;
			else certs[k] = (strC[k].length()==0 )? 
					Double.NaN : 
						Double.parseDouble(strC[k]);
		}
	}
	@Override
	public boolean containedIn(Set<Integer> cnvTypesToInclude) {
	for(int i=0; i<cnts.length; i++){
		if(cnts[i]>0 && cnvTypesToInclude.contains(i)) return true;
	}
	return false;
	}

	@Override
	public Integer fixedPos() {
		return null;
	}

	@Override
	public int get(int i) {
		return i<cnts.length ? cnts[i] : 0;
	}
	public double cert(int i) {
		return i<cnts.length ? certs[i] : 0;
	}
	@Override
	public void addto(SortedSet<Integer> cnvTypes) {
		for(int i=0; i<cnts.length; i++){
			if(cnts[i]>0) cnvTypes.add(i);
		}
		
	}

	@Override
	public void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
			int[] certCount) {
		for(int i=0; i<this.cnts.length; i++){
			cnvCounts[i]+=cnts[i];
			sb.append(cnts[i]);
			if(i<cnts.length-1) sb.append(":");
			certCount[i]+=cnts[i];
			if(cnts[i]>0) avgCert[i]+=certs[i]*cnts[i];
		}
	}

	@Override
	public int max() {
		return cnts.length;
	}
  }
public abstract void append(StringBuffer sb, int[] cnvCounts, double[] avgCert,
		int[] certCount);
	// TODO Auto-generated method stub
public abstract int max() ;
public abstract  boolean containedIn(Set<Integer> cnvTypesToInclude) ;
	


}


>>>>>>> .r259
