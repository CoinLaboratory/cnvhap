/**
 * 
 */
package lc1.dp.swing;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import javax.swing.JPanel;

import lc1.ensj.GenesForRegion;
import lc1.util.Constants;

import org.ensembl.datamodel.KaryotypeBand;
import org.ensembl.datamodel.Location;
import org.freehep.graphics2d.VectorGraphics;
import org.jfree.chart.ChartPanel;
import org.jfree.data.Range;

public class GenePanel extends JPanel implements MouseListener{
	double y_min = 2;
      double y_0 = 10;
      double y_1 =  20;
      double y_11 =  15;
      double y_111 =  10;
      double y_2 = 30;
      
      static final String[] chrs_ = "1:10:11:12:13:14:15:16:17:18:19:2:20:21:22:3:4:5:6:7:8:9".split(":");
    static final	String[] length_ =   "247249719:135374737:134452384:132349534:114142980:106368585:100338915:88827254:78774742:76117153:63811651:242951149:62435964:46944323:49691432:199501827:191273063:180857866:170899992:158821424:146274826:140273252".split(":");

      
    //  double x_offset = 0;
      Dimension dim5 = new Dimension(500,100);
ChartPanel cp;
      //List<Integer> location;
      static GenesForRegion grf = new GenesForRegion();
      
      SortedMap<String, List<int[]>> exons;
      double reg_length;
     // int noSnps;
      final Font font;
      ColorAdapter ca = new ColorAdapter();
      double sp=0;
     double ep=0;  
     Map<String, String> ids = new HashMap<String, String>();
     public void setStartEnd1(){
    	 l.clear();
    	 int j_up =0;
    	 boolean up = true;
         int j_down =0;
         if(exons!=null){
  	 for(Iterator<Map.Entry<String, List<int[]>>> it = exons.entrySet().iterator(); it.hasNext();){
        Map.Entry<String, List<int[]>>  nxt = it.next();
        String nme = nxt.getKey();
        if(nme ==null) nme = "null";
        boolean drawnName = false;
    
     for( int i=0; i<nxt.getValue().size(); i++){
       int[] pos = nxt.getValue().get(i);
       
        // double x_loc = getPos(pos[0]);
        // double x_loc1 = getPos(pos[1]);
         Entity1 ent = new Entity1(nme, pos, ids.get(nme), !drawnName, up, j_up, j_down);
         if(!drawnName ){
      	   drawnName = true;
             if(up){
             	   j_up++;
             	   if(j_up>=numIncr){
             		   j_up =0;
             	   }
                }
                if(!up){
             	   j_down++;
             	   if(j_down>=numIncr){
             		   j_down =0;
             	   }
                }
                up =!up;
         }
        
         this.l.add(ent);
     }
  	 }
         }
     }
   
      public void setStartEnd(String chr, String chr_end, double min, double max){
    	  this.sp = min;
    	  this.ep = max;
    	  this.chr =chr;
    	  this.chr_end = chr_end;
    	  if(chr.equals(chr_end) && (max - min) < Constants.showEnsThresh()){
    	  reg_length = ep-sp;
          try{
          Location loc = new Location("chromosome:"+chr+":"+(int)Math.round(min)+"-"+(int)Math.round(max));
        ids.clear();
          exons = grf.fetchGenesByLocation(loc, ids);
          System.err.println("exons are"+exons);
         
         
         if(exons!=null){ 
        	//  setStartEnd1();
         }
          }catch(Exception exc){
              exc.printStackTrace();
          }}
      }
      
     public  static double[] chr_length;
     public static List[] karyotypes;
      static{
    	  try{
	     chr_length = new double[24];
	     karyotypes = grf.fetchKaryotypes();
	      grf.fetchSequenceRegionsSuchAsChromosomeOrContig(chr_length);
	    	  }catch(Exception exc){
	    		  exc.printStackTrace();
	    		  for(int i=0; i<chrs_.length; i++){
	    			  String name = chrs_[i];
	    			  if(name.equals("X")){
	    				   chr_length[22] =Integer.parseInt(length_[22]);
	    			   }
	    			  else if(name.equals("Y")){
	    				   chr_length[23] =Integer.parseInt(length_[23]);
	    			   }
	    			   else{
	    				   try{
	    					   int pos = Integer.parseInt(name);
	    					 
	    					   chr_length[pos-1] = Integer.parseInt(length_[i]);
	    				   }catch(Exception exc1){
	    					   System.err.println(name);
	    				   }
	    			   }
	    		  }
	    	  }
      }
      String chr, chr_end;
      int chr_, chr_end_;
      /* start end are as fraction of chromosome length */
      GenePanel(double min, double max,  Font f, ChartPanel cp, double width) {
         //this.location = location;
    		int chr =(int)  Math.floor(min);
    		int chr1 =(int)  Math.floor(max);
    	
    		this.setStartEnd(chr, chr1, min-chr, max-chr1);
         this.cp = cp;
         this.font = f;
         this.width = width;
         //this.chr_length = chr_l;
          this.setPreferredSize(dim5);
          this.setSize(dim5);
          this.setMinimumSize(dim5);
       //  this.setStartEnd(chr, chr_end, start, end);
         this.addMouseListener(this);
      }
     
      final double width;
      public GenePanel(String chr, Integer start, Integer end, Font font10, ChartPanel  cp, int width) {
    	  this.setStartEnd(chr, chr, start, end);
    	  this.cp = cp;
    	  this.width = width;
          this.font = font10;
          //this.chr_length = chr_l;
           this.setPreferredSize(dim5);
           this.setSize(dim5);
           this.setMinimumSize(dim5);
        //  this.setStartEnd(chr, chr_end, start, end);
          this.addMouseListener(this);
	}

public Color color(String st){
	
	if (st.equals("gneg")) return Color.white;
	else if(st.startsWith("gpos")){
		int i =(int) Math.floor(255*(Double.parseDouble(st.substring(4))/100.0));
		return new Color(i,i,i);
	}
	else return Color.gray;
	
}
	double prev_loc =0;
     MiniBrowser browser, browser1;
      int numIncr = 10;
      
      List<Entity1> l = new ArrayList<Entity1>();
      
      class Entity1{
    	  public Entity1(String nme, int[] pos, String id, boolean drawName, boolean up, int j_up,int  j_down) {
			this.name = nme;
			this.drawName = drawName;
			this.up = up;
			this.pos0 = pos[0];
			pos1 = pos[1];
			this.j_down = j_down;
			this.j_up = j_up;
		//	this.x_loc = x_loc;
		//	this.x_loc1 = x_loc1;
			url = "http://www.ensembl.org/Homo_sapiens/geneview?gene="+id;
		//	url1 = "http://www.genenames.org/data/hgnc_data.php?hgnc_id="+id2;
		}
    	  String url, url1;
    	 int pos0, pos1;
	//	double x_loc,x_loc1;
		boolean drawName;
		int j_up, j_down;
		boolean up;
    	  String name;
    	  //URL url;
      }
      
      public void paint(Graphics g){
    	   rec = new Rectangle2D.Double();
    	 if(cp!=null){
    		
    	        rec.setRect(cp.getScreenDataArea());
    	       
    		//cp.getChart().getPadding().trim(rec);
    		 double x = rec.getX();
    		 double y = rec.getY();
    		 double w = rec.getWidth();
    		 double h = rec.getHeight();
 	        rec.setRect(x+10, y, w-10, h);
    	 }
    	 else{
    		  rec.setRect(this.getBounds());
    	 }
    	 this.paint(g,rec);
//          armitage[k].getChart().getPadding().trim(r);
      }
      public Double getLocation(double frac){
    	  double pos = frac*(this.max - this.min)+min;
    	  int chrP =(int) Math.floor(pos);
    	  if(chrP>=1){	
    	  Double l = this.chr_length[chrP-1]*(pos-(double) chrP)/1000000;
    	  return l;
    	  }
    	  else return null;
      }
      public void paint(Graphics g, Rectangle2D r){
    	  this.setStartEnd1();
          g.setColor(Color.white);
     //     this.x_offset = cp!=null ? getWidth() -( cp.getScreenDataArea().getWidth() ): 0;
          double h = getHeight();
          this.y_0 = 0.4*h;
          this.y_1 = 0.9*h;
          this.y_11 = 0.5*h;
          this.y_2 =h;
         this.y_111 = 0.2*h;
         double incr = ((1.0/3.0)*h)/(double) numIncr; 
          g.fillRect(0,0, getWidth(), getHeight());
          
          VectorGraphics vg = VectorGraphics.create(g);
          vg.setFont(this.font);
      
          	vg.setColor(Color.BLACK);
          vg.setLineWidth(0.5);
         /* {
              double x_0 = getX(0, true);
              double x_1 = getX(noSnps-1, true);
            
              
              g.setColor(background_color);
              g.fillRect(0,0, getWidth(), getHeight());
               vg = VectorGraphics.create(g);
               vg.drawRect(x_0, y_0, x_1 - x_0, y_1-y_0);//(arg0, arg1, arg2, arg3)
          }*/
       
          FontRenderContext frc = vg.getFontRenderContext();
          for(double ii=0; ii<10; ii++){
        	  double frac = ii/10.0;
        	  Double loc = getLocation(frac);
        	  if(loc!=null){
        	  vg.setColor(l.size()==0 ? Color.black : Color.LIGHT_GRAY);
        	  double x_pos = frac*r.getWidth()+r.getMinX();
        	  vg.drawString(String.format("%5.3f mb", new Object[] {loc}), x_pos, 8.0);
              vg.drawLine(x_pos, 8,x_pos,8+y_111);
        	  }
        	  
          }
          if(l.size()==0){
        	  double tot_len = this.max - this.min;
        	  double wid = r.getWidth()/tot_len; //width for each chrom
        	  for(int ii= chr_; ii<chr_end_; ii++){
        		  double offset = (ii - min)*wid+r.getMinX();
        		  paint(vg, offset, wid, ii);
        	  }
        	 
        	  
        	//  for(int chr_i = this.chr_; chr_i < )
          }
         for(int ij=0; ij<l.size(); ij++){
              Entity1 ent = l.get(ij);
             
         //     boolean drawnName = false;
              Color c = ca.getColor(ent.name);
              vg.setColor(c);
            
          
             
               double x_loc = this.getPos(ent.pos0);
               double x_loc1 = this.getPos(ent.pos1);
               if(ent.drawName){
            	   GlyphVector gv = font.createGlyphVector(frc,ent.name.toCharArray());
            	   Shape outline = gv.getOutline();
                   Rectangle2D obounds = outline.getBounds2D();
                /*   AffineTransform at = new AffineTransform();
                   at.setToTranslation(x_loc-obounds.getMinX(), (up ? (y_min+j_up*incr) : y_11+j_down*incr) -obounds.getMinY());
                
                   outline  = at.createTransformedShape(outline);*/
                 //  j_up =0;
                  vg.drawString(ent.name, x_loc-obounds.getMinX(), (ent.up ? (y_min+ent.j_up*incr)-obounds.getMinY() : y_2-ent.j_down*incr) );
                  
                   //vg.fill(outline);
                //   vg.drawString( x_loc, up ? y_0 : y_2);
                 
               }
         
               vg.drawRect(x_loc, y_0, x_loc1 - x_loc, y_111);
          }
      }
    
	private void paint(VectorGraphics vg, double offset, double wid, int chr_2) {
		if(chr_2<1) return;
		
		List kary = karyotypes[chr_2-1];
		double sz = chr_length[chr_2-1];
		vg.setFont(this.font);
		double prev_loc =offset;
		for(int i=0; i<kary.size(); i++){
			 KaryotypeBand kb = (KaryotypeBand)kary.get(i);
			 Location loc = kb.getLocation();
		//	System.err.println("stain "+kb.getStain());
			 vg.setColor(this.color(kb.getStain()));
			 double st = (loc.getStart() / sz)*wid;
			 double end = (loc.getEnd() / sz)*wid;
				boolean centLeft = i>=1 && kb.getBand().startsWith("q") && ((KaryotypeBand)kary.get(i-1)).getBand().startsWith("p");
				boolean centRight = i<kary.size()-1 && kb.getBand().startsWith("p") && ((KaryotypeBand)kary.get(i+1)).getBand().startsWith("q");
				if(centLeft){
					double[] x = new double[] {end+offset, end+offset, st+offset};
					double[] y = new double[] {y_0+y_11, y_0, y_0+y_11/2.0};
					vg.fillPolygon(x, y, 3);
					vg.setColor(Color.BLACK);
					vg.drawPolygon(x, y, 3);
				}
				else if(centRight){
					double[] x = new double[] {st+offset, st+offset, end+offset};
					double[] y = new double[] {y_0+y_11, y_0, y_0+y_11/2.0};
					vg.fillPolygon(x, y, 3);
					vg.setColor(Color.BLACK);
					vg.drawPolygon(x, y, 3);
					
				}
				else{
					double rnd = i==0 ? 20.0 : 0.0;
					 vg.fillRoundRect(st+offset, y_0, end-st, y_11,rnd, rnd);
					vg.setColor(Color.BLACK);
					vg.drawRoundRect(st+offset, y_0, end-st, y_11,rnd, rnd);
					if(st+offset - prev_loc > 30){
					  vg.drawString(kb.getBand(),st+offset, y_0-2);
					  prev_loc = st+offset;
					}
		                
				}
		}
	}
	// final double[] chr_length;
    private double getPos(double pos) {
     //  double  wid =  ((double)rec.getWidth()  ) / reg_length;
       double  wid =  (width ) / reg_length;
      //  int pos1 = location.get(sp);
        return wid * (pos - sp)+rec.getMinX();
    }
	public void setRange(Range range) {
		 min = range.getLowerBound();
		 max = range.getUpperBound();
		int chr =(int)  Math.floor(min);
		int chr1 =(int)  Math.floor(max);
		this.setStartEnd(chr, chr1, min-chr, max-chr1);
		
		
	}
	double min, max;
	
	public void setStartEnd(int chr, int chr1, double min, double max){
		String chrom = chr==23 ? "X": (chr)+"";
		String chrom1 = chr1==23 ? "X": (chr1)+"";
		this.chr_end_ = chr1;
		this.min = chr+min;
		this.max = chr1+max;
		if(chr>=1)
		this.setStartEnd(chrom,chrom1, (min)*chr_length[chr-1], (max)*chr_length[chr1-1]);
		
	}
	Rectangle2D rec;
	public void mouseClicked(MouseEvent e) {
		double x = e.getX();
		Entity1 ent = null;
		double dist = Double.POSITIVE_INFINITY;
		for(int i=0; i<this.l.size(); i++){
			double dist1 = Math.abs(this.getPos(this.l.get(i).pos0) - x);
			if(dist1<dist){
				ent = l.get(i);
				dist = dist1;
			}
		}
		if(ent!=null){
			if(this.browser==null){
				browser = new MiniBrowser();
				browser.setAddress(ent.url);
				browser.actionGo();
				browser.show();
				
			}
		}
		
	}
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
  }