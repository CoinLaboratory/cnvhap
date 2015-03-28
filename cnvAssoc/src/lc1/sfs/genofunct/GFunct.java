/**
 * 
 */
package lc1.sfs.genofunct;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;

// we assume ploidy of 2 at this stage
public abstract class GFunct{
	//abstract double calc(int[] i, int[] pop1);

	public String toString(){
		return this.getName();
	}
	
	//pops is the number of individuals, whereas for a Funct (alleles) pops is number of alleles
	public abstract double calc(double[] v, double[] vhet, int[] pops);

public static GFunct[] getFuncts(String[] str, List<String> nme)  {
	List<GFunct> l = new ArrayList<GFunct>();
	try{
	
	for(int i=0; i<str.length; i++){
		String[] str1 = str[i].split("_");
		Class clazz = Class.forName("lc1.sfs.genofunct."+str1[0]);
		if(str1.length>1){
			String[] str2 = str1[1].split(";");
			for(int k=0; k<str2.length; k++){
				String[] str3 = str2[k].split("-");
				Double[] arg = new Double[2];
				arg[0] = Double.parseDouble(str3[0]);
				arg[1] = str3.length==1 ? 1.0 :  Double.parseDouble(str3[1]);
				Constructor c = clazz.getConstructor(new Class[] {Double.class, Double.class});
				GFunct fun = (GFunct)c.newInstance(arg);
				l.add(fun);
				nme.add(fun.getName());
			}
		}else{
			GFunct fun = (GFunct)clazz.newInstance();
			nme.add(fun.getName());
		   l.add(fun);	
		}
	}
	}catch(Exception exc){
		exc.printStackTrace();
	}
	return l.toArray(new GFunct[0]);
}



 public abstract String getName() ;

}
