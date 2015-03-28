/**
 * 
 */
package lc1.sfs.freqfunct;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;

public abstract class Funct{
	//abstract double calc(int[] i, int[] pop1);
 
	
	
	public abstract double calc(double[] v, int[] pops);

public static Funct[] getFuncts(String[] str, List<String> nme)  {
	List<Funct> l = new ArrayList<Funct>();
	try{
	
	for(int i=0; i<str.length; i++){
		String[] str1 = str[i].split("_");
		
		Class clazz = Class.forName("lc1.sfs.freqfunct."+str1[0]);
		if(str1.length>1){
			String[] str2 = str1[1].split(";");
			for(int k=0; k<str2.length; k++){
				String[] str3 = str2[k].split("-");
				Double[] arg = new Double[2];
				arg[0] = Double.parseDouble(str3[0]);
				arg[1] = str3.length==1 ? 1.0 :  Double.parseDouble(str3[1]);
				Constructor c = clazz.getConstructor(new Class[] {Double.class, Double.class});
				Funct fun = (Funct)c.newInstance(arg);
				l.add(fun);
				nme.add(fun.getName());
			}
		}else{
			Funct fun = (Funct)clazz.newInstance();
			nme.add(fun.getName());
		   l.add(fun);	
		}
	}
	}catch(Exception exc){
		exc.printStackTrace();
	}
	return l.toArray(new Funct[0]);
}



 public abstract String getName() ;

}
