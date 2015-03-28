/**
 * 
 */
package lc1.util;

public class PhenoGroup{
	double min;
	double max;
	String str;
	public PhenoGroup(String string) {
		this.str = string;
		if(string.startsWith("<")){
			min = Double.NEGATIVE_INFINITY;
			max = Double.parseDouble(string.substring(1));
		}
		else if(string.startsWith(">")){
			max = Double.POSITIVE_INFINITY;
			min = Double.parseDouble(string.substring(1));
		}
		else{
			String[] str = string.split("\\^-");
			min = Double.parseDouble(str[0]);
			max = Double.parseDouble(str[1]);
		}
	}
	public String toString(){
		return str;
	}
	public boolean  within(String v){
		if(v.equals("NA")) return false;
		if(v.equals(str)) return true;
		if(v.length()==0) return false;
		try{
		double v1 = Double.parseDouble(v);
		if(v1>min && v1<max) return true;
		else {
			return false;
		}
		}catch(Exception exc){
			exc.printStackTrace();
			return false;
		}
	}
}