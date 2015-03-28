package lc1.dp.swing;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.freehep.graphics2d.VectorGraphics;

public class SummaryLogoPanel extends JPanel {
	public static void main(String[] args){
		try{
			Map<String, Double > m = new HashMap<String, Double>();
			m.put("AABAA", 1.0);
			m.put("BBBBB", 1.2);
			showTables(new SummaryLogoPanel(m), "test");
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	private static void showTables(JPanel jp, String title) {
		Dimension dim = new Dimension(500,500);
		jp.setSize(dim);
		jp.setMinimumSize(dim);
		jp.setPreferredSize(dim);
        //   JTable table = new JTable(model);
           JFrame fr = new JFrame(title);
           fr.setDefaultCloseOperation(fr.DISPOSE_ON_CLOSE);
           JScrollPane jB = new JScrollPane(jp, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
          // fr.setSize(new Dimension(500,800));
           fr.getContentPane().add(jB);
           fr.pack();
           fr.setVisible(true);
            
        }
	
	SummaryLogo left;
	SummaryLogo right;
	double noSnps;
	double total =0;
	int midSNP;
	public SummaryLogoPanel(Map<String , Double> m){
		String str = m.keySet().iterator().next();
		noSnps = str.length();
		 midSNP = (int) Math.floor(((double)str.length()-1.0)/2.0);
		Map<String, Double> leftM  = new HashMap<String, Double>();
		Map<String, Double> rightM = new HashMap<String, Double>();
		for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String nxt = it.next();
			Double value = m.get(nxt);
			total+=value;
		}
		for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String nxt = it.next();
			Double value = m.get(nxt);
			String rights = nxt.substring(midSNP);
			append(rightM, rights, value/total);
//			rightM.put(rights, value/total);
			StringBuffer lefts = new StringBuffer();
			for(int i=midSNP; i>=0; i--){
				lefts.append(nxt.charAt(i));
			}
			append(leftM, lefts.toString(), value/total);
//			leftM.put(lefts.toString(), value/total);
		}
		
		right = new SummaryLogo(rightM);
		left = new SummaryLogo(leftM);
		
	}
	private void append(Map<String, Double> m, String key, Double v){
		if(m.containsKey(key)){
			m.put(key, m.get(key)+v);
		}
		else{
			m.put(key, v);
		}
	}
	
 public void paint( Graphics g ) {
	  g.setColor(Color.white);
      g.fillRect(0,0, this.getWidth(), this.getHeight());
		 double height = this.getHeight();
		 double width = this.getWidth()/noSnps;
		 double mid = width*midSNP;
		 VectorGraphics vg = VectorGraphics.create(g);
		 vg.setColor(Color.BLACK);
		 vg.drawRect(mid, 0, width, height);
		 left.paint(vg, mid, 0, width,  height, false);
		 right.paint(vg, mid, 0, width,  height, true);
	 }
	
}
