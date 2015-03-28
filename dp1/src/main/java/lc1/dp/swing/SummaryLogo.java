package lc1.dp.swing;

import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.freehep.graphics2d.VectorGraphics;

public class SummaryLogo  {

	
	//Rectangle rect;
	final Map<Character, Double> heights = new HashMap<Character,Double>();
	
	SummaryLogo(Map<String, Double> m){
		Map<Character, Map<String, Double>> m1 = new HashMap<Character, Map<String, Double>>();
		double total = 0;
		for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String st = it.next();
			Double value = m.get(st);
			total+=value;
			Character ch = st.charAt(0);
			Double v = heights.get(ch);
			heights.put(ch, v==null ? value : value+v);
			if(st.length()==1){
			    rightMap.put(ch, null);
			}
			else{
				Map<String,Double> m2 =  m1.get(ch);
				if(m2==null){
					m1.put(ch, m2 = new HashMap<String, Double>());
				}
				
				String substring = new String(st.substring(1));
				m2.put(substring, value);
			}
		}
		if(m1.size()>0){
			for(Iterator<Character> it = m1.keySet().iterator(); it.hasNext();){
				Character ch = it.next();
				rightMap.put(ch, new SummaryLogo(m1.get(ch)));
			}
		}
		
	}
	
	Map<Character, SummaryLogo> rightMap = new TreeMap<Character, SummaryLogo>();
	static Map<Character, Color> m = new HashMap<Character, Color>();
	static{
		m.put('-', Color.red);
		m.put('X', Color.green);
		m.put('Y', Color.green);
		m.put('Z', Color.green);
		m.put('Z', Color.green);
		m.put('A', Color.BLUE);
		m.put('B', Color.BLUE);
	}
	//double total;
	public static  Font font= new Font("Arial", Font.PLAIN, 10);
	  
	public void paint(VectorGraphics vg, double offset_w, double offset_h, double width, double tot_height, boolean right){
		  FontRenderContext frc = vg.getFontRenderContext();
		  double cum_ht =offset_h;
		for(Iterator<Character> it = rightMap.keySet().iterator(); it.hasNext();){
			Character ch = it.next();
			
			
			GlyphVector gv = font.createGlyphVector(frc, new char[] {ch});
            Shape outline = gv.getOutline();
            Rectangle2D obounds = outline.getBounds2D();
            AffineTransform at = new AffineTransform();
            double height = this.heights.get(ch)*tot_height;
           
            at.setToTranslation(0.0+offset_w, cum_ht ); 
            //1
           
            //sc = 1;
            at.scale((width/obounds.getWidth()), (height/obounds.getHeight()));
            at.translate(-obounds.getMinX(), -obounds.getMinY());
            Color col = m.get(ch);
            if(col==null) col = Color.black;
            vg.setColor(col);
            outline  = at.createTransformedShape(outline);
            vg.fill(outline);
            vg.draw(outline);
            SummaryLogo obj = rightMap.get(ch);
            if(obj!=null){
	           if(right) {
	        	   obj.paint(vg, offset_w+width, cum_ht, width, tot_height,  right);
	           }
	           else{
	        	   obj.paint(vg, offset_w-width, cum_ht, width, tot_height, right);
	           }
	          
            }
            cum_ht+=height;
		}
	}
	
}

