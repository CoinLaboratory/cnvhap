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
import java.util.List;
import java.util.Properties;
import java.util.logging.Logger;

import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.SwingUtilities;

import lc1.dp.data.collection.Phenotypes;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.FreeHaplotypeHMM;
import lc1.dp.states.EmissionState;
import lc1.dp.transition.AbstractTransitionProbs;
import lc1.stats.IntegerDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SkewNormal;
import lc1.stats.TrainableNormal;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/*@Author Lachlan Coin*/
public class HMMPanel extends JPanel implements PropertyChangeListener{
    final List<Integer> location;
    final List<Character> major;
    final List<Character> minor;
    Color  background_color = Color.WHITE,
            line_color = Color.BLACK,
            font_color = Color.BLACK;
   // int x_start =0;
  //  int x_end;
    
   
 
    int mult = 25;
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
final double noIndiv;
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
    	this.mat =//this.hmm.getHittingProb(this.hmm.noSnps, mat);
    		this.hmm.getCounts(this.hmm.noSnps, mat);
    		
    	
    	this.ratePanel.removeAll();
    	this.ratePanel.add(this.getChart());
    	SwingUtilities.invokeLater(new Runnable(){
    		public void run(){
    			updateUI();
    		}
    	});
    	}
     
    }
    
    Phenotypes pheno;
    

    InnerPanel ip;
   
    LocPanel locP;
   
   
    public HMMPanel(FreeHaplotypeHMM hmm, List<Integer> location, List<Character> major, List<Character> minor, int noIndiv, Phenotypes pheno, File outdir){
 
	 {
		
		ca = ColorAdapter.get(hmm);
	    
	  }
	 this.maxLineWidth = 100.0 * (4.0/(double)hmm.modelLength());
	 this.ratePanel = new JPanel();
	 
	 for(int k=0; k<hmm.modelLength()-1;k++){
		 this.rates.addSeries(new XYSeries(""+k));
	 }
	this.location = location;
    this.major = major;
    this.minor = minor;
    this.noIndiv = noIndiv;
    this.pheno = pheno;
   this.chartBF = IndividualPlot.makeOrDelete(outdir, "hmmF");
  // this.hittingProb = hmm.getHittingProb(hmm.noSnps);
    String[] stToP = Constants.phenToPlot();
    this.setLayout(new BorderLayout());
    this.locP = new LocPanel(location,location.size()*mult, 100, this.offset_x, "hmm");
    JSplitPane sjp =   new JSplitPane1(JSplitPane.VERTICAL_SPLIT,  
       		ip = new InnerPanel(),new JSplitPane(JSplitPane.VERTICAL_SPLIT,locP, ratePanel)
     );
    sjp.setDividerLocation(0.8); 
    this.add(sjp);
  
   /* if(stToP.length==0){
        this.addTab("no phen", new InnerPanel());
      
    }
    else{
        for(int i=0; i<stToP.length; i++){
            this.addTab(stToP[i], new InnerPanel(stToP[i]));
        }
        
    }*/
  //  emStSp =Emiss.getSpaceForNoCopies(Constants.backgroundCount());
    	// ((EmissionState)hmm.getState(1)).getEmissionStateSpace();
  //  List<Integer> cn = emStSp.copyNumber;
   /* for(int i=0; i<cn.size(); i++){
        int[] geno = emStSp.getGenoForCopyNo(cn.get(i));
        for(int k=0; k<geno.length; k++){
            ca1.getColor(geno[k]+"");
        }
    }*/
   mat = new double[hmm.noSnps][hmm.modelLength()];
	//hmm.getCounts(hmm.noSnps, mat);
       hmm.getHittingProb(hmm.noSnps, mat);
   /* for(int i=0; i<mat.length; i++){
        for(int j=1; j<mat[i].length; j++){
            mat[i][j] = transform(mat[i][j]); 
        }
    }*/
    this.numFounders = hmm.modelLength()-1;
    this.noSnps = hmm.noSnps;
     this.hmm = hmm;
     this.ratePanel.add(this.getChart());
    dim = new Dimension(location.size()*mult, len);
     this.setMinimumSize(dim);
     this.setPreferredSize(dim);
     this.minloc = location.get(0);
     this.reg_length = location.get(noSnps-1)-minloc;
 }


Dimension dim;
     FreeHaplotypeHMM hmm;
   
   // final int index=0;
  
   //  final double textOffset;
  
  public double getX(int i, boolean useLoc){
      if(useLoc) return offset_x +(((double)location.get(i)-minloc)/this.reg_length) *(this.getWidth()-(offset_x+offset_x)) ;
      else return offset_x +((((double)i)/(double)this.noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
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
 
  
   double  prev_loc =0;
   
   final int minloc;
   final double reg_length;
    double maxLineWidth = 100;
   static int len = 300;
   static int width = 1000;
   final int noSnps;
   final int numFounders;
   final double offset_x = mult;
  // final double offset_xR = mult;
   double offset_y = 0;
   
   static Color LIGHTGREEN = new Color(0,100,0);
  static Color DARKGREEN  = new Color(189, 183, 107);
 
 
  public  ColorAdapter ca;
  
  
 // static ColorAdapter ca1 = new ColorAdapter(new Color[] {Color.BLUE, Color.RED, Color.RED, Color.GREEN, Color.GREEN, Color.GREEN});
     
  final JPanel ratePanel;
  XYSeriesCollection rates = new XYSeriesCollection();
 // XYSeriesCollection ratesSeries = new XYSeriesCollection();
  
  public double  updateRates(){
	for(int k=0; k<rates.getSeriesCount(); k++){
		rates.getSeries(k).clear();
	}
	  double upper = Double.NEGATIVE_INFINITY;
	  AbstractTransitionProbs[] trans = hmm.trans.transProbs;
	  for(int i=1; i<trans.length; i++){
	
		 AbstractTransitionProbs rates = trans[i];
		
		 for(int j=0; j<this.rates.getSeriesCount(); j++){
			 double r = Math.log10(rates.getRate(j));
			// if(Double.isNaN(r)) throw new RuntimeException("!");
			 if(r>upper){
				 upper = r;
			 }
			 this.rates.getSeries(j).add((double)i, r);
			 }
		
	  }
	  return upper;
  }
  public ChartPanel getChart(){
	 double upperBound =  this.updateRates();
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
   public InnerPanel(String stToP) {
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
    }
public void paint0(){
    /*   vg.setColor(Color.BLACK);
      drawOval(getStartX(), getStartY(),max_shape, 0, 360);
       if( hmm!=null){
          // double[] trans  = new double[numFounders];
           for(int k=0; k<numFounders; k++){
            //   double shape_k = max_shape* mat[1][k+1];
               double val  = transform2(hmm.getTransitionScore(0, k+1, 0));
               vg.setColor(Color.black);//getLineColor(val,  Color.BLACK));
               vg.setLineWidth(maxLineWidth *val);
               
//               vg.drawLine(getStartX(),getStartY(), HMMPanel.this.getX(0, false),  HMMPanel.this.getY(k));
           }
       }*/
   }
   void drawOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
       vg.drawArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
   }
   void fillOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
       vg.fillArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
   }
   public void write(int i){
       double x_start = HMMPanel.this.getX(i, false);
     // double textOffset =  ((((double)0.5)/(double)noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
       for(int j=0; j< numFounders;
       j++){
           double y_start = HMMPanel.this.getY(j);
           vg.setColor(Color.BLACK);
           PseudoDistribution ems =((EmissionState)hmm.getState(j+1)).emissions(i);
           Integer fixed = ems.fixedInteger();
           EmissionStateSpace emStSp=  ((EmissionState)hmm.getState(j+1)).getEmissionStateSpace();
           
           if(fixed!=null){
        	   emStSp = ((IntegerDistribution)ems).emstsp;
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
           double y_start = HMMPanel.this.getY(j);
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
       double x_start = HMMPanel.this.getX(i, false);
       vg.setColor(Color.BLACK);
     //  vg.setLineWidth(0.5);
       double x_loc = HMMPanel.this.getX(i, true);
     //  vg.drawLine(x_loc, y_0, x_loc, y_1);
      // vg.drawLine(x_loc, y_1, x_start, HMMPanel.this.getY(0));
     /*  if(i==0 || x_loc> 100+prev_loc){
           String str = HMMPanel.this.location.get(i)+"";
           vg.setColor(Color.BLACK);
           vg.drawString(str, x_loc, 10);
           prev_loc = x_loc;
          
       }*/
       double x_start1 = i<noSnps-1 ? HMMPanel.this.getX(i+1,false) : HMMPanel.this.getX(i,false);
       for(int j=0; j< numFounders;
       j++){
    	   if(Constants.CHECK){
    		   double sum = Constants.sum(mat[i]);
    		   try{
    		   if(Math.abs(1.0-sum)>0.001) throw new RuntimeException("!! "+sum);
    		   }catch(Exception exc){
    			   exc.printStackTrace();
    			   System.exit(0);
    		   }
    	   }
    	   
    	   
           double y_start = HMMPanel.this.getY(j);
          Color c = ca.getColor(j+"");
        //   vg.setColor(c);
           double shape_j =max_shape* transform2( mat[i][j+1]);
         //  vg.setLineWidth(line_width);
         //  drawOval(x_start, y_start, shape_j);
          EmissionState st = (EmissionState)hmm.getState(j+1);
           PseudoDistribution dist = (st).emissions(i);
         
           Integer fixed = dist.fixedInteger();
           int max;
           Color c1 ;
           if(fixed!=null){
               max = fixed;
               c1 = c;//this.ca1.getColor(fixed+"");
             
              
           }
           else{
              double[] probs =  ((SimpleExtendedDistribution)dist).probs();
               max = Constants.getMax(probs);
               EmissionStateSpace emStSp = st.getEmissionStateSpace();
              int[] comp = emStSp.getGenoForCopyNo(emStSp.getCN(max));
              double frac = mod(probs, max, comp);
             //System.err.println(emStSp.get(max)+" "+comp.length+" "+frac);
              c1 = HMMPanel.this.modify(c, frac);
             // c1 = Color.black;
           }
           vg.setColor(c1);
           ProbabilityDistribution[][] emDT = 
               ((EmissionState)hmm.getState(j+1)).emissionsDatatype();
           ProbabilityDistribution distDT = phenIndexToPaint >=0 && emDT!=null && emDT.length>phenIndexToPaint ? emDT[phenIndexToPaint][i] : null;
           if(distDT==null ){
               fillOval(x_start, y_start, shape_j, 0, 360);
           }
           else{
               double[] counts = distDT.getCount(angle1);
               for(int k=0; k<angle.length; k++){
//                   double[] angle_k = angle[k];
                   double shape_jk =shape_j*counts[k]*(double)angle.length;
                   double t = k ==angle.length-1 ? 360 : angle[k+1];
                   fillOval(x_start, y_start, shape_jk, angle[k], (t-angle[k]));
               }
           }
          
              
         // vg.fillOval(x_start, y_start, shape_j,shape_j);
         // vg.setColor(getColor());
         
         // vg.setLineWidth(HMMPanel.this.maxLineWidth);
        //  double[][] mat1 = mat;
          double mult  = Constants.plotFlux() ? mat[i][j+1] : 0.2;
          if(i<noSnps-1 && hmm!=null){
              for(int k=0; k<numFounders; k++){
            	//  hmm.getHittingProb(hmm.noSnps);
            	  double pr = hmm.getTransitionScore(j+1, k+1, i+1);
                 double val = transform2(pr*(mult));
                 
                  vg.setColor(c);//getLineColor(val, c));
                  double wid = HMMPanel.this.maxLineWidth*val;
              
                //  vg.setStroke(getStroke(Math.pow(pr,1.0)));
                  
                 if(wid>0) {
                	 vg.setLineWidth(wid);
                	 vg.drawLine(x_start,HMMPanel.this.getY(j), x_start1,  HMMPanel.this.getY(k));
                 }
              }
          }
         
       
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
       
       g.setColor(HMMPanel.this.background_color);
       g.fillRect(0,0, getWidth(), getHeight());
        vg = VectorGraphics.create(g);
    //    vg.drawRect(x_0, y_0, x_1 - x_0, y_1-y_0);//(arg0, arg1, arg2, arg3)
       wid =  ((double)getWidth() - (2*offset_x) ) / ((double) reg_length);
       wid1 =  ((double)getWidth() - (2*offset_x) ) / ((double) noSnps);
       height  = ((double)getHeight()-2*offset_y ) / ((double) numFounders-1.0);
       shape =Math.min(wid1, height);
       max_shape = Math.min(shape*5, Math.max(wid1, height));
     
       vg.setFont(font16);
       paint0();
       for(int i=0; i<noSnps; i++){
           paint(i);
           write(i);
       } 
       write0();
       for(int i=0; i<noSnps; i++){
         
           write(i);
       } 
       
      
   }
  }
   static final double log2 = Math.log(2);
   double[][] mat;
   
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
    
    
    
   private double transform(double d) {
	   //System.err.println(" t "+d);
	   if(Constants.showAll(1))return d;
	   if(Constants.plotFlux && d<1e-4 || ! Constants.plotFlux && d <0.05) return 0;
	   else return d;
//       return Math.pow(d,2);//Math.pow(10, Math.log(d)/Math.log(2));
    }
   
   private double transform2(double d) {
	   //System.err.println(" t "+d);
	  // if(d<0.0001) return 0;
	  // else return d;
   return Math.pow(d, Constants.hmmBubblePow());//Math.pow(d,2);//Math.pow(10, Math.log(d)/Math.log(2));
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
             ratePanel.setSize(dim);
           int lp_height = locP.getSize().height;
             if(Constants.printPlots() && plot==2){
            	 {
            	
            	 System.err.println("printing hmm");
            	// int maxSNPS =Constants.maxSNPS();
               //  int n = (int) Math.ceil((double) this.noSnps/ maxSNPS);
            	//  for(int i=0; i<n; i++){
            		//  for(int j=0; j<this.getTabCount(); j++){
			            	  Properties p = new Properties();
			           			p.setProperty("PageSize", "A4");
			           			boolean inclLocP = false;
			           			//int len1 = maxSNPS;//Math.min(maxSNPS, location.size()  - i*maxSNPS);
			           	Dimension dim1 = new Dimension(dim.width,dim.height+this.ratePanel.getHeight()+(inclLocP? lp_height : 0));	
			           	AffineTransform at = new AffineTransform();
			             ImageGraphics2D g = new ImageGraphics2D(new File(this.chartBF, "hmm.png")
			             ,dim1, ImageConstants.PNG) ;
			             g.setProperties(p);
			             g.startExport();
			             at.setToTranslation(0, 0);
			             g.setTransform(at);
			             this.ip.print(g);
			             at.setToTranslation(0,dim.height );
			              g.setTransform(at);
			             this.ratePanel.print(g);
			           if(inclLocP){ 
			        	   at.setToTranslation(0,dim.height +ratePanel.getHeight());
			              g.setTransform(at);
			           
			              this.locP.print(g);
			           }
			             g.endExport();
			             System.err.println("done printing hmm");
            		 // }
            	//  }
            	 }
            if(false)	{
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
            		
            	 }
            	 
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
