package lc1.sfs;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;


public class ProbMap extends TreeMap<Double, Double>{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	public void add(double key, double prob){
		this.put(key, this.get(key)+prob);
	}
	public void init(){
		for(Iterator<Map.Entry<Double, Double>> it = this.entrySet().iterator(); it.hasNext();){
			it.next().setValue(0.0);
		}
	}
	
	/**  */
	public void getQuantiles(double[] d, double[] quants, double[] v){
		double sum =0;
		int index =0;
		for(Iterator<Map.Entry<Double, Double>> it = this.entrySet().iterator(); it.hasNext() && index < d.length;){
			Map.Entry<Double, Double> nxt = it.next();
			sum+=nxt.getValue();
			while(index < d.length && sum>=d[index] ){
				quants[index] = nxt.getKey();
				v[index] = sum;
			}
			//it.next().setValue(sum);
		}
	}
	
	
}
