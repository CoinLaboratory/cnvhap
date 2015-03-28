package lc1.dp.swing;

import java.util.ArrayList;
import java.util.List;

import org.jfree.data.xy.DefaultXYZDataset;

public class DefaultXYZDataset1 extends DefaultXYZDataset {
	List< double[][]>m = new ArrayList<double[][]>();
@Override
public void addSeries(Comparable seriesKey, double[][] data) {
	super.addSeries(seriesKey, data);
	this.m.add( data);
}
public void updateSeries(int i1, int i, Integer integer, Double double1,
		Double double2) {
	double[][] d = m.get(i1);
	d[0][i] = integer;
	d[1][i] = double1;
	d[2][i] = double2 ;
	
}
public void clear(){
	m.clear();
}
}
