/**
 * 
 */
package lc1.dp.swing;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JPanel;

import org.freehep.graphics2d.VectorGraphics;

public class NamePanel extends JPanel{
	 String name;
	 NamePanel(String name){
		 this.name = name;
	
		 this.setSize(100, 100);
	 }
		@Override
		public void paint(Graphics g){
			  g.setColor(Color.white);
	          g.fillRect(0,0, this.getWidth(), this.getHeight());
	          VectorGraphics vg = VectorGraphics.create(g);
	          vg.setFont(HaplotypePanel1.font12);
	          vg.setColor(Color.black);
	         vg.drawString(name, 10, 10);
	     
			
	
		}
 }