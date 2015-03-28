package lc1.dp.swing;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Arrays;
import java.util.Properties;
import java.util.logging.Logger;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.SwingUtilities;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.FreeHaplotypeHMM;
import lc1.dp.model.FreeSiteTrans1;
import lc1.dp.states.EmissionState;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;

/*@Author Lachlan Coin*/
public class RateMatrixPanel extends JPanel implements PropertyChangeListener{
 
    Color  background_color = Color.WHITE,
            line_color = Color.BLACK,
            font_color = Color.BLACK;
   // int x_start =0;
  //  int x_end;
    
    final double[][] mat;
 
   final  int mult = 500;
//final File hmmF;
    
   /* public void write(){
        try{
        for(int i=0; i<this.getComponentCount(); i++){
            InnerPanel p = (InnerPanel) this.getComponent(i);
            this.writeToZipFile(p,hmmF,  p.getName());
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    public void writeToZipFile(Component charts, File dir, String id) throws Exception{
        File out = new File(dir, (id)+".png");
        ImageGraphics2D g = new ImageGraphics2D(out,charts, ImageConstants.PNG); 
        g.startExport();
        charts.print(g);
        g.endExport();
    }*/
//final double noIndiv;
    Font small_font,
        large_font,
        small_italic_font,
        large_italic_font;
    FontMetrics fm_small,
                fm_large,
                fm_small_italic,
                fm_large_italic;

    Logger logger = Logger.global.getAnonymousLogger();
   // int[] y_loc;
  //  int domain_thickness = 30;
   // double ratio;
   
 //   double[][] mat;
    //double[][] hittingProb;
    public void update(){
    	if(plot==2|| Constants.plot()>=2){
    		 
    	this.overallRate = hmm.trans.updateRates(mat, pi);
    	//	this.hmm.getCounts(this.hmm.noSnps, mat);
    		
    	
   
    	SwingUtilities.invokeLater(new Runnable(){
    		public void run(){
    			updateUI();
    		}
    	});
    	}
     
    }
    
   // Phenotypes pheno;
    

    InnerPanel ip;
   
   
   
   
 double[] pi;
 double overallRate;  
 public RateMatrixPanel(FreeHaplotypeHMM hmm,   File outdir){
 
	 {
		
		ca = ColorAdapter.get(hmm);
	    
	  }
	 
	 this.maxLineWidth = 100.0 * (4.0/(double)hmm.modelLength());
	
	
	//this.location = location;
    //this.major = major;
	 
    //this.minor = minor;
    //this.noIndiv = noIndiv;
    //this.pheno = pheno;
   this.chartBF = IndividualPlot.makeOrDelete(outdir, "hmmRate");
  // this.hittingProb = hmm.getHittingProb(hmm.noSnps);
   // String[] stToP = Constants.phenToPlot();
    this.setLayout(new BorderLayout());
   // this.locP = new LocPanel(location,location.size()*mult, 100, this.offset_x, "hmm");
    ip = new InnerPanel();
   JScrollPane jp =  new JScrollPane(
   		ip,
   		  JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
  
 this.add(jp);
   mat = new double[hmm.rateLength()-1][hmm.rateLength()-1];
   pi = new double[hmm.rateLength()-1];
   overallRate = (hmm.trans).updateRates(mat, pi);
	//hmm.getCounts(hmm.noSnps, mat);
     //  hmm.getHittingProb(hmm.noSnps, mat);
   /* for(int i=0; i<mat.length; i++){
        for(int j=1; j<mat[i].length; j++){
            mat[i][j] = transform(mat[i][j]); 
        }
    }*/
    this.numFounders = hmm.modelLength()-1;
//    this.noSnps = hmm.noSnps;
     this.hmm = hmm;
  //   this.ratePanel.add(this.getChart());
     int locsize = numFounders+1;
    dim = new Dimension(locsize*mult, len);
     this.setMinimumSize(dim);
     this.setPreferredSize(dim);
    ip.setMinimumSize(dim);
     ip.setPreferredSize(dim);
     noSnps = this.numFounders;
   //  this.minloc = location.get(0);
    // this.reg_length = location.get(noSnps-1)-minloc;
 }
final int noSnps;

Dimension dim;
     FreeHaplotypeHMM hmm;
   
   // final int index=0;
  
   //  final double textOffset;
  
  public double getX(int i, boolean useLoc){
	  double w = this.getWidth();
	  double offx = offset_x;
      //if(useLoc) return offset_x +(((double)location.get(i)-minloc)/this.reg_length) *(this.getWidth()-(offset_x+offset_x)) ;
      double res =  offset_x +((((double)i)/(double)this.noSnps)) *(w-(offset_x+offset_x)) ;
      return res;
  }
  
  public double getY(int i){
      return ((double)i) *height +offset_y;
  }
 
  
  
  public double getStartX(){
      return 5.0;
  }
 
  public double getStartY(){
      return (double)  (this.getHeight())/2.0;
  }

  public static  Font font16 = new Font("SansSerif", Font.PLAIN, Constants.hmmFontSize());
   double wid;
   double wid1;
   double height;
  double shape;
  double max_shape;
 
  
   //double  prev_loc =0;
   
 //  final int minloc;
  // final double reg_length;
    double maxLineWidth = 100;
   static int len = 300;
 //  static int width = 1000;
   //final int noSnps;
   final int numFounders;
   final double offset_x = 100;//2*mult;
  // final double offset_xR = mult;
   double offset_y = 0;
   
   static Color LIGHTGREEN = new Color(0,100,0);
  static Color DARKGREEN  = new Color(189, 183, 107);
 
 
  public  ColorAdapter ca;
  
  
 // static ColorAdapter ca1 = new ColorAdapter(new Color[] {Color.BLUE, Color.RED, Color.RED, Color.GREEN, Color.GREEN, Color.GREEN});
     
  
  
 // XYSeriesCollection ratesSeries = new XYSeriesCollection();
  
 
 /* public ChartPanel getChart(){
	// double upperBound =  this.updateRates();
	  final JFreeChart chart = ChartFactory.createXYLineChart(
             "Log Rate",
               "Relative position ", // domain axis label
            null, // range axis label
               rates, // data
               PlotOrientation.VERTICAL, false, // include legend
               true, // tooltips?
               false // URL generator? Not required...
               );
	  final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
	  for(int i=0; i<rates.getSeriesCount(); i++){
		  Color c = ca.getColor(i);
		  renderer.setSeriesPaint(i, c);
		  Integer noC = ((EmissionState)this.hmm.getState(i+1)).noCop();
		  if(noC!=null){
		  renderer.setSeriesShapesVisible(i, noC.intValue()==1);
		  }
	  }
	  chart.getXYPlot().setRenderer(renderer);

	//  chart.getXYPlot().getRangeAxis().get
	  //chart.getXYPlot().setRangeAxis(new LogarithmicAxis("rates"));
	  final ChartPanel cp =   new ChartPanel(
			  chart,
					 this.getWidth(), //width
		  		 200, //height
		  		this.getWidth(), //mindrawWidth
		           100, //mindrawHeight
		           this.getWidth(), //maxDrawWith
		           400,//maxDrawHeight
		           ChartPanel.DEFAULT_BUFFER_USED,
		           true,  // properties
		           true,  // save
		           true,  // print
		           true,  // zoom
		           true   // tooltips		
		  
		  );
	  {
	 Range r =  chart.getXYPlot().getDomainAxis().getRange();
	 chart.getXYPlot().getDomainAxis().setAutoRange(false);
	 chart.getXYPlot().getDomainAxis().setRange(1,r.getUpperBound());
	  }
	  {
	 Range r =  chart.getXYPlot().getRangeAxis().getRange();
	if(r.getLowerBound() < upperBound){
	 chart.getXYPlot().getRangeAxis().setAutoRange(false);
	 chart.getXYPlot().getRangeAxis().setRange(r.getLowerBound(), upperBound);
	}
	  }
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
	 cp.setBorder(null);
	  return cp;
  }
 */ 
  class InnerPanel extends JPanel{
      VectorGraphics vg;
      final int phenIndexToPaint;// =0;
      final double[] angle;
      final double[] angle1;
      
      public InnerPanel(){
          this.phenIndexToPaint = -1;
          this.angle = null;
          this.angle1 = null;
          
      }
 /*  public InnerPanel(String stToP) {
       this.phenIndexToPaint =pheno.phen.indexOf(stToP);
       ProbabilityDistribution dist = pheno.phenotypeDistribution[phenIndexToPaint];
       if(dist instanceof SkewNormal){
           this.angle = new double[Constants.segments()];
       }
       else this.angle = new double[((SimpleExtendedDistribution)dist).probs.length];
       double incr = 1.0/(double) angle.length;
       for(int i=0; i<angle.length; i++){
           angle[i] = i*incr*360;
       }
       if(dist instanceof SkewNormal){
           this.angle1 = ((TrainableNormal)dist).getAngle(angle.length);
              
       }
       else this.angle1 = angle;
    }*/
      
      public void paint0(){
          double x_start = RateMatrixPanel.this.getX(0, false);
        
       
          for(int j=0; j< numFounders;
          j++){
       	
              double y_start = RateMatrixPanel.this.getY(j);
             Color c = ca.getColor(j+"");
         
            
//             double[] pi1 = pi;
              double shape_j =pi[j]*10;
            	  //max_shape* transform2( mat[i][j]);
          
              vg.setColor(c);
           
                  fillOval(x_start, y_start, shape_j, 0, 360);
                  vg.setColor(Color.black);
                  vg.drawString(
           			   String.format("%5.3g",pi[j]),
           			   
                           x_start, y_start);
          
          }
        
      }

   void drawOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
       vg.drawArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
   }
   void fillOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
       vg.fillArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
   }
   public void write(int i){
       double x_start = RateMatrixPanel.this.getX(i, false);
     // double textOffset =  ((((double)0.5)/(double)noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
       for(int j=0; j< numFounders;
       j++){
           double y_start = RateMatrixPanel.this.getY(j);
           vg.setColor(Color.BLACK);
           Integer fixed = ((EmissionState)hmm.getState(j+1)).emissions(i).fixedInteger();
           EmissionStateSpace emStSp=  ((EmissionState)hmm.getState(j+1)).getEmissionStateSpace();
           
           if(fixed!=null){
        	   int bcount = emStSp.getBCount(fixed);
        	   vg.drawString(
        			   1.0+"",
                        x_start, y_start+bcount*Constants.hmmFontSize);
        	   
           }
           else{
        	   EmissionState stat = (EmissionState)hmm.getState(j+1);
           double[] probs= (stat).emissions(i).probs();
           double rem =1.0;
           for(int k=0; k<probs.length; k++){
        	 
        	   {
        		   double p =  Math.round((k==probs.length-1 ? rem :probs[k])*100.0)/100.0;
        		   rem -=p;
        		   if(p>0){
        		   int k1 = emStSp.getBCount(k);
        		   vg.drawString(
            			  //emStSp.getGenotypeString(emStSp.get(k))+
            			p+"",
                            x_start, y_start+k1*Constants.hmmFontSize);
        		   }
        		  // sb.append(emStSp.get(k)+":"+probs[k]+"\n");
        	   }
           }
           }
         
          
       }
   }
   
   public void write0(){
       double x_start = 0;
     // double textOffset =  ((((double)0.5)/(double)noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
       for(int j=0; j< numFounders;
       j++){
           double y_start = RateMatrixPanel.this.getY(j);
           vg.setColor(Color.BLACK);
          
           EmissionStateSpace emStSp=  ((EmissionState)hmm.getState(j+1)).getEmissionStateSpace();
           
          
           Integer cn = ((EmissionState)hmm.getState(j+1)).noCop();
          
           if(cn!=null){
        	   int[] geno = emStSp.getGenoForCopyNo(cn);
	           for(int k=0; k<geno.length; k++){
	        	   
	        		  
	        		   vg.drawString(
	            			  emStSp.getGenotypeString(emStSp.get(geno[k]))
	            		,
	                            x_start, y_start+k*Constants.hmmFontSize);
	        		   
	        		  // sb.append(emStSp.get(k)+":"+probs[k]+"\n");
	        	   
	           }
           }
           else{
        	   for(int k=0; k<emStSp.genoListSize(); k++){
	        	   
	        		  
        		   vg.drawString(
            			  emStSp.getGenotypeString(emStSp.get(k))
            		,
                            x_start, y_start+k*Constants.hmmFontSize);
        		   
        		  // sb.append(emStSp.get(k)+":"+probs[k]+"\n");
        	   
           }
           }
           }
         
          
       }
  
   public void paint( int i){
       double x_start = RateMatrixPanel.this.getX(i+1, false);
       vg.setColor(Color.BLACK);
    
       for(int j=0; j< numFounders;
       j++){
    	
           double y_start = RateMatrixPanel.this.getY(j);
          Color c = ca.getColor(j+"");
      
         // double[][]mat1 = mat;
         // double[] pi1 = pi;
         // double overallR = overallRate;
           double shape_j =max_shape* transform2( Math.abs(mat[i][j])/overallRate);
       
           vg.setColor(c);
        if(i!=j){
               fillOval(x_start, y_start,shape_j, 0, 360);
        }else{
        	
        	  drawOval(x_start, y_start, shape_j, 0, 360);
        }
               vg.setColor(Color.BLACK);
              
            	   vg.drawString(
            			   String.format("%5.3g; %5.3g", new Double[] {mat[i][j]/overallRate,
            					   overallRate}),
            			   
                            x_start, y_start);
            	   
             
       
       }
     
   }
   float[] dash = new float[2];
   private Stroke getStroke(double d) {
	   dash[0] = Math.max(0,(float) (d*10.0));
	   dash[1] = Math.max(0,(float) ((1-d)*10.0));
	return new BasicStroke(1.0f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);
}
public void paint( Graphics g ) {
     //  double x_0 = HMMPanel.this.getX(0, true);
    //   double x_1 = HMMPanel.this.getX(noSnps-1, true);
     offset_y = ((double)getHeight())*0.1;
       
       g.setColor(RateMatrixPanel.this.background_color);
       g.fillRect(0,0, getWidth(), getHeight());
        vg = VectorGraphics.create(g);
    //    vg.drawRect(x_0, y_0, x_1 - x_0, y_1-y_0);//(arg0, arg1, arg2, arg3)
   //    wid =  ((double)getWidth() - (2*offset_x) ) / ((double) reg_length);
       wid1 =  ((double)getWidth() - (2*offset_x) ) / ((double) noSnps);
       height  = ((double)getHeight()-2*offset_y ) / ((double) numFounders-1.0);
       shape =Math.min(wid1, height);
       max_shape = Math.min(shape*5, Math.max(wid1, height));
     
       vg.setFont(font16);
       paint0();
       for(int i=0; i<noSnps; i++){
           paint(i);
         //  write(i);
       } 
      /* write0();
       for(int i=0; i<noSnps; i++){
         
         //  write(i);
       }*/ 
       
      
   }
  }
   static final double log2 = Math.log(2);
  // double[][] mat;
   
   public double mod( double[] probs, int max, int[] comp ){
       double sum = 0;
       for(int i=0; i<comp.length; i++){
           sum+= -probs[comp[i]]*(Math.log(probs[comp[i]])/log2);
       }
       double info = -probs[max]*(Math.log(probs[max])/log2);
      return  transform1((sum-info)/sum);
   }
 public static Color modify(Color c, double frac){
   
     try{
     return new Color(c.getRed(), c.getGreen(), c.getBlue(), (int) Math.floor(c.getAlpha()*(frac)));
     }catch(Exception exc){
    	 exc.printStackTrace();
    	 System.err.println(frac);
    	 return Color.BLACK;
        // System.exit(0);
     }
   //  return null;
 }
  
  
  // final EmissionStateSpace emStSp;
   
   
   Character uncertain = '-';
   
  
   
 //   double y_0 = 10;
  //  double y_1 = 20;
    
    
    
  
   
   private double transform2(double d) {
	   //System.err.println(" t "+d);
	  // if(d<0.0001) return 0;
	  // else return d;
   return Math.log(d)/100;
 //  Math.pow(d, Constants.hmmBubblePow());//Math.pow(d,2);//Math.pow(10, Math.log(d)/Math.log(2));
    }
   private double transform1(double d) {
       return Math.pow(10, Math.log10(d)/Math.log(1.5));
    }
private Color getLineColor(double d, Color c) {
       try{
      return new Color(c.getRed(),  c.getGreen(), c.getBlue(),(int)Math.round(d*255));
       }catch(Exception exc){
           System.err.println(d);
           exc.printStackTrace();
           System.exit(0);
//           return null;
       }
       return null;
    }


int plot = 0;
public void setToPlot(int i){
   plot= i;
}
   final File chartBF;

   
   public synchronized void propertyChange(PropertyChangeEvent arg0) {
    String nme = arg0.getPropertyName();
  /*if(nme.equals("emiss")){
        if(plot<2 && Constants.plot()<2) return;
        
	   Object[] obj = (Object[]) arg0.getNewValue();
	 
	   StateDistribution dist = (StateDistribution) obj[0];
	   Integer l = (Integer) obj[1];
	   if(l==0) clear();
	   Integer i = (Integer)obj[2];
	   this.addInformation(dist,i);
   }*/    
    if(nme.equals("setToPlot")){
        int level = (Integer) arg0.getNewValue();
        setToPlot(level);
        return;
    }
    if(nme.equals("done") ){
    	// if(plot<2 && Constants.plot()<2) return;
          //this.update();
    	try{
    		
    		ip.setMinimumSize(dim);
            
             
             ip.setSize(dim);
      //       ratePanel.setSize(dim);
       //    int lp_height = locP.getSize().height;
             if(Constants.printPlots() && plot==2){
            	 {
            	ip.updateUI();
            	 System.err.println("printing hmm");
            	// int maxSNPS =Constants.maxSNPS();
               //  int n = (int) Math.ceil((double) this.noSnps/ maxSNPS);
            	//  for(int i=0; i<n; i++){
            		//  for(int j=0; j<this.getTabCount(); j++){
			            	  Properties p = new Properties();
			           			p.setProperty("PageSize", "A4");
			           		
			           			//int len1 = maxSNPS;//Math.min(maxSNPS, location.size()  - i*maxSNPS);
			           	Dimension dim1 = new Dimension(dim.width,dim.height);	
			         //  	AffineTransform at = new AffineTransform();
			             ImageGraphics2D g = new ImageGraphics2D(new File(this.chartBF, "hmm.png")
			             ,ip, ImageConstants.PNG) ;
			             g.setProperties(p);
			             g.startExport();
			         //    at.setToTranslation(0, 0);
			         //    g.setTransform(at);
			             this.ip.print(g);
			          
			             g.endExport();
			             System.err.println("done printing hmm");
            		 // }
            	//  }
            	 }
         /*   if(false)	{
            		  Properties p = new Properties();
	           			p.setProperty("PageSize", "A4");
	           			ChartPanel cp = (ChartPanel) this.ratePanel.getComponent(0);
	           			cp.setSize(dim);
	           			ImageGraphics2D g = new ImageGraphics2D(new File(this.chartBF, "rates.png"), cp,
	        					ImageConstants.PNG);
	        		
	        			g.setProperties(p);
	        			g.startExport();
	        			cp.print(g);
	        			g.endExport();
            		
            	 }*/
            	 
             }
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
     }

  
   /* if(nme.equals("init")){
       
    }
    else if(nme.equals("emiss")){
       
    }
  
    else 
    else if(nme.equals("expec_i")){
        
    }
    else if(nme.equals("expectation1")){
       
      //  this.updateUI();
    }
    else if(nme.equals("dist_maximisation")){
       
    }
    else*/ 
        if(nme.equals("hmm_maximisation")){
        update();
    }
   // else throw new RuntimeException("!! "+nme);
    
}

private void clear() {
	for(int i=0; i<this.mat.length; i++){
		Arrays.fill(mat[i], 0.0);
	}
	
}
   

}
