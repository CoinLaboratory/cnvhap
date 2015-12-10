package lc1.dp.swing;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

public class Rainbow {

	public List<Color> colors = new ArrayList<Color>();
	
	public Color get(int i){
		return colors.get(i);
	}
	/*public void addCols(int[] start, int [] end, int steps){
		for(double k=0; k<steps; k++){
			int[] mid = new int[3];
			double prop = k/steps;
			for(int j=0; j<3; j++){
					mid[j] = (int) Math.floor((1-prop) * start[j] + prop * end[j]);
			}
		}
		
	}*/
	
	public Rainbow(int ndel1, int nampl1, int nampl2, int nampl3){
	//    colors.add(Color.orange);
		//int ndel = ndel1-1;
	//	int nampl1 = nampl11-1;
		
		for(int b=0; b<ndel1; b++)colors.add(new Color(255, 0,b*255/ndel1));//  [red -> purple)
		colors.add(Color.gray);
		for(int r=nampl1; r>0; r--)colors.add(new Color(r*255/nampl1, 0, 255));//    [purple ->   blue)
		//colors.add(new Color(0, 0*255/100,  255)); //blue
		for(int g=0; g<nampl2; g++)colors.add(new Color(0, g*255/nampl2,  255)); // [blue -> cyan )
		for(int b=nampl3;b>0;  b--)colors.add(new Color(0, 255,b*255/nampl3)); // [cyan -> green)
		//for(int r=0; r<n; r++)colors.add(new Color(r*255/100, 255, 0));
		//for(int g=100; g>0; g--)colors.add(new Color(255, g*255/100,  0));
		//
		//colors.add(new Color(0,255,0));
	}
	public static Color[] getColors(int maxcn,  int bgcount) {
		//int ndel=bgcount;
		//int nampl1 = bgcount;
		//int nampl2 = Math.max(0,2*bgcount -(nampl1+bgcount)); 
		int nampl3 = Math.max(0,maxcn -(3*bgcount)); 
		Rainbow r = new Rainbow(bgcount, bgcount,bgcount, nampl3);
		return r.colors.toArray(new Color[0]);
	}
	
}
