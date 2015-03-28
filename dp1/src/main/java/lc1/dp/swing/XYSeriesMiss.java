package lc1.dp.swing;

import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;

public class XYSeriesMiss extends XYSeries{

	//int nan_both =0;
	double nan_x = 0;
	//int nan_y = 0;
	public String toString(){
		return this.getKey().toString();
	}
	public XYSeriesMiss(Comparable key) {
		super(key);
		// TODO Auto-generated constructor stub
	}
	public void add(XYDataItem item, boolean notify){
		if(Double.isNaN(item.getX().doubleValue()) || 
			Double.isNaN(item.getY().doubleValue())){
			
				nan_x++;
			
		}
		else{
			super.add(item, notify);
		}
	}

}
