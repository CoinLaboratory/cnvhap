package lc1.dp.swing;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.font.GlyphVector;
import java.util.Arrays;

import javax.swing.JPanel;

import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;

public class PanelName extends JPanel {
	public static  Font font16 = new Font("Arial", Font.PLAIN, 16);//(int)(30*ROC.resolution));
	 double[] frac = null;
	public void setFrac(double[] frac1){
		double sum = Constants.sum(frac1);
		for(int k=0; k<frac1.length; k++){
			
			frac[k] = frac1[k] /sum;
		}
	 }
	public PanelName(String[] names, Dimension dim, int[] offset){
		this.names = names;
		this.setSize(dim);
		this.setMinimumSize(dim);
		this.offset = offset;
		this.frac = new double[offset.length];
		Arrays.fill(frac, 1.0/(double)frac.length);
		int width = this.getWidth();
		System.err.println("width");
	}
	int[] offset;
	final String[] names;
	@Override
	public void paint(Graphics g){
		VectorGraphics vg = VectorGraphics.create(g);
		double height = this.getHeight();
		double width = this.getWidth();
		   g.setColor(Color.WHITE);
	       g.fillRect(0,0, getWidth(), getHeight());
		vg.setColor(Color.black);
		vg.setFont(font16);
		double h_p =0;
		//double h = height/(double) names.length;
		for(int i=0; i<names.length; i++){
			double h = height* frac[i];
			GlyphVector gv = font16.createGlyphVector(((Graphics2D)g).getFontRenderContext(), names[i].toCharArray());
			double h1 = gv.getOutline().getBounds().getHeight();
			vg.drawString(names[i], 0, h_p+h1+offset[i]);
			h_p+=h;
		}
	}
	
}
