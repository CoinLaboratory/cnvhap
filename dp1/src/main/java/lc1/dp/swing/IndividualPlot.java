package lc1.dp.swing;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.Sampler;
import lc1.dp.data.collection.AssociationCalculator;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LinearRegressionCalculator;
import lc1.dp.data.collection.MatchedDistributionCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.stats.CompoundDistribution;
import lc1.stats.IlluminaDistribution;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.IntegerDistribution;
import lc1.stats.Mixture;
import lc1.stats.Mixture2;
import lc1.stats.MixtureDistribution;
import lc1.stats.OrthogonalProbabilityDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.PseudoDistribution;
import lc1.stats.PseudoMixture;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SkewNormal;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.freehep.graphicsio.svg.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.AxisChangeEvent;
import org.jfree.chart.event.AxisChangeListener;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;
import org.jfree.util.ShapeList;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;


public class IndividualPlot extends JTabbedPane implements PropertyChangeListener{
    
    final List<Integer> index;
    final String[] names1;
    String toPlotString = "expectation2";
    final boolean plotMerged;
    double minrv = Constants.logplot() ? 1e-3: Double.NEGATIVE_INFINITY;
      //"hmm_maximisation"
    final boolean useVals;
 final   Dimension dim5 = new Dimension(Constants.scatterWidth(), Constants.scatterWidth());
    
   static  int width = 	 Constants.r_panel_width();
final static int type_index =0; //which type 0 is CN, 1 is B
   static int height = Constants.plotHeight();
   static Dimension dim2 = new Dimension(width, height);
   static final  Dimension dim = new Dimension(width, height);
   static final Dimension dimBoth =  new Dimension(width, height*2);
   static final  Dimension dim_ = new Dimension(width, height*2);
   static final  Dimension dimG = new Dimension(800, 800);
   static EmissionStateSpace[] stateEm=  new EmissionStateSpace[2];//null;
   CompoundEmissionStateSpace emStSp;
 EmissionStateSpace emStSp1;
 final MarkovModel hmm;
//   final int[] b_alias;
  boolean singlechrom = false;
   
    JComponent jpB;
    JComponent[] jB, jR;
    JComponent[] chartB;
	JComponent[] jG = new JComponent[2];  
//	List<ColoredLine> []lines;
    JComponent[] chartR;
 //   JComponent[][]chartDist;
    JComponent[][][]chartG = new JPanel[2][][]; //first dimension is normal,global
    ChartPanel[] chartBNew, chartRNew;
    ChartPanel[][][]chartGNew = new ChartPanel[2][][]; //first dimension is normal, global
    JComponent jpBS, jpRS;
    JComponent[] jGS = new JComponent[2]; //normal, global;
    public static  Font font16 = new Font("SansSerif", Font.PLAIN, (int) Math.round(4* Constants.shapeSize1()));
    public static  Font font20 = new Font("SansSerif", Font.BOLD, (int) Math.round(4*Constants.shapeSize1()));
    public static  Font font10 = new Font("SansSerif", Font.PLAIN,(int) Math.round(4*Constants.shapeSize1()));
    public static  Font font101 = new Font("SansSerif", Font.PLAIN,(int) Math.round(3*Constants.shapeSize1()));

    public static  Font font4 = new Font("SansSerif", Font.PLAIN, (int) Math.round(4*Constants.shapeSize1()));
    final int noSnps;
    BaumWelchTrainer bwt;
    List<Integer> location;
    
    List<String> rbSeriesNames = new ArrayList<String>();
    public XYSeriesCollection[] getRSeriesCollection(int k){
    	if(Constants.r_panel()){
    	    int[] b_alias = this.emStSp.haploPairToGeno();
      XYSeriesCollection[] r = rdc.get(k);
      if(r==null){
           r = new XYSeriesCollection[this.index.size()];
           for(int i=0; i<r.length; i++){
        	   r[i] = new XYSeriesCollection();
          
          
        	   for(int ij=0; ij<b_alias.length; ij++){
                   r[i].addSeries(new XYSeries(emStSp.get(ij)+""));
               }
           
           }
            rdc.put(k, r);
            
      }
           return r;
    	}
    	else return null;
      
   }
  
public XYSeriesCollection[] getBSeriesCollection(int k){
	if(Constants.b_panel()){
           int[] b_alias = this.emStSp.haploPairToGeno();
           XYSeriesCollection[] b =  bdc.get(k);
           if(b==null){
                b =new XYSeriesCollection[this.index.size()];
                for(int i=0; i<b.length; i++){
                	b[i] = new XYSeriesCollection();
               for(int ij=0; ij<b_alias.length; ij++){
                   b[i].addSeries(new XYSeries(emStSp.get(ij)+""));
               }
                }
               bdc.put(k, b);
           }
           return b;
	}
	else return null;
   }
    
   final Map<Integer, XYSeriesCollection[]> rdc,bdc;  // rdc ordered by distribution, bdc ordered by emissionstatespace order
   
  // XYSeriesCollection current_b, current_r;
   public File chartDistF, chartBF, chartRF;
   public static File makeOrDelete(File outdir, String st){
       File fi = new File(outdir, st);
       if(!fi.exists()) fi.mkdir();
       else{
       File[] f = fi.listFiles();
       for(int i=0; i<f.length; i++){
           f[i].delete();
       }
       }
       return fi;
   }
 List<Integer> indicesToInclude = new ArrayList<Integer>();
 //String[][] probesToPlot;
 //double[][] pvalue;
 //int[][] snp_alias;
 static int divsize = 10;
 
 public LocPanel updateLocPanel(int[] alias){
	 int start =SignificancePlot.start_;
	
	 int end = SignificancePlot.end_+1;
	 if(alias[start]>=0 && alias[end-1]>=0)
	 return  new LocPanel(this.location.subList(alias[start], alias[end-1]+1), Constants.scatterWidth()*end-start, 50,(double)Constants.scatterWidth()/2.0, "scatter");
	 else return new LocPanel(this.location, Constants.scatterWidth()*end-start, 50,(double)Constants.scatterWidth()/2.0, "scatter");
 }
 
 
 public LocPanel updateLocPanel(int[] alias, int start, int max_len){
	 int end = Math.min(alias.length-1, start+max_len);
	 return  new LocPanel(this.location.subList(alias[start], alias[end]), Constants.scatterWidth()*end-start, 50, (double)Constants.scatterWidth()/2.0, "scatter");
 }
 
//public  static boolean updateProbes=false;
 //AssociationCalculator acalc;
 final int phenoIndex;
    public IndividualPlot(BaumWelchTrainer bwt, String[] names, List<Integer> loc,List<String> snpid, List<Integer> intensityPr,
    		String name, File outdir, int ploidy, short phenoIndex, boolean plotMerged){
        super();
        this.plotMerged = plotMerged;
        this.useVals = phenoIndex>=0;
        this.phenoIndex = Math.max(0,phenoIndex);
        this.ploidy = ploidy;
       /* String[] toann = Constants.rsToAnnotate();
        if(toann!=null){
        	int ind = snpid.indexOf(toann[0]);
        	if(ind>=0)
        	SignificancePlot.minId[0] = ind;
        }*/
        String[] samp = Constants.annotateSamples;
        this.samplesToAnnotate =
        	samp==null ? new HashSet<String>() :
        	new HashSet<String>(
        		Arrays.asList(samp));
        this.snpid = snpid;
        this.names1 = plotMerged ? new String[] {"merged"} : names;
        this.setName(name);
        this.hmm = bwt.hmm;
       // this.lines = new List[snpid.size()];
        String[]nme = DataCollection.datC.getUnderlyingDataSets();
       Integer[] ind = new Integer[nme.length];
       
        for(int i=0; i<ind.length; i++){
        	ind[i] = i;
        }
        if(plotMerged) ind = new Integer[]{0};
        this.index=  Arrays.asList(ind);
        this.bwt = bwt;
       // this.b_frac = DistributionCollection.dc.b_frac();
     
       // else b_alias =null;
        this.location = loc;
        this.chrom = Constants.chrom0();
        {
        	int pos1 = (int)Math.floor(Constants.decode(loc.get(loc.size()-1), this.pos, false));
        	int pos2 = (int) Math.floor(Constants.decode(loc.get(0), this.pos, false));
        if( pos1==pos2){
        	this.singlechrom = true;
        	if(Constants.scaleLoc!=null) this.chrom = this.pos[0]+"";
        	
        }
        }
     //   this.ca_b = new ColorAdapter[index.size()];
      //  this.ca_r = new ColorAdapter[index.size()];
        this.emStSp =Emiss.getSpaceForNoCopies(ploidy);
        this.genoToAnnotate = new HashSet<Integer>();
        int[] cnToAnnotate = Constants.cnToAnnotate();
       
        if(cnToAnnotate!=null){
        	 int[] bafToAnnotate = Constants.bafToAnnotate();
        for(int i=0; i<cnToAnnotate.length; i++){
        	if(cnToAnnotate[i]<0){
        		genoToAnnotate.add(emStSp.genoListSize());
        	}
        	else{
        		int[] toadd = emStSp.getGenoForCopyNo(cnToAnnotate[i]);
        		if(toadd!=null){
        			for(int k=0; k<toadd.length; k++){
        				if(bafToAnnotate==null || emStSp.getBCount(toadd[k]) == bafToAnnotate[0] )
        					genoToAnnotate.add(toadd[k]);
        			}
        		}
        	}
        }
        }
     //   this.b_alias = emStSp.haploPairToGeno();
        int maxi =0;
        this.strokes = new Map[index.size()];
        for(int i=0; i<index.size(); i++){
        //	ca_r[i] = DistributionCollection.dc.getCAr();//[index.get(i)].ca_r;
        //	ca_b[i] = DistributionCollection.dc.getCAb();
        	if(index.get(i)>maxi){
        		maxi = index.get(i);
        	}
        	strokes[i] = new HashMap<String, Stroke>();
        	
        }
        for(int ij=0; ij<this.emStSp.genoListSize(); ij++){
            
          
         	  rbSeriesNames.add(emStSp.get(ij)+"");
         	 
          
        }
        this.data_index_alias = new int[maxi+1];
        Arrays.fill(data_index_alias, -1);
        for(int i=0; i<index.size(); i++){
        	data_index_alias[index.get(i)] = i;
        }
        if(SignificancePlot.snp_alias==null ) SignificancePlot.init(snpid, this.names1.length);
       this.rGlobal = new XYSeriesCollection[this.index.size()][2];
       for(int i=0; i<index.size(); i++){
    	   for(int k=0; k<rGlobal[i].length; k++){
    		   rGlobal[i][k] = this.newRBSeries();
    	   }
    	  // this.rGlobal[i] = this.newRBSeries();
    	   //this.rGlobalProbeOnly[i] = this.newRBSeries();
       }
        this.noSnps = loc.size()-1;
      //  chartDist = new JPanel[index.size()][3];
      //  JComponent[] jpDist = getBottomPane(chartDist);
        chartBF = makeOrDelete(outdir, "plot_"+this.getName()+"_predictionsB");
        chartRF = makeOrDelete(outdir, "plot_"+this.getName()+"_predictionsR");
        chartDistF = makeOrDelete(outdir, "plot_"+this.getName()+"_distributions");
       this.emStSp1 = emStSp.getMembers()[0];
        rdc = new HashMap<Integer, XYSeriesCollection[]>();
        bdc = new HashMap<Integer, XYSeriesCollection[]>();
       this.rb = new Map[index.size()];
       this.annotation = new Map[index.size()];
       for(int i=0; i<rb.length; i++){
    	   rb[i] =  new HashMap<Integer, XYSeriesCollection>();
    	   annotation[i] = new HashMap<Integer,  AnnotationSeries>();
       }
      
      if(Constants.showScatter()) this.updateClusterSeries();
      
        this.include = new Boolean[bwt.data.length];
       
        if(Constants.plot()>=2){
        	 jB = new JPanel[this.index.size()];
             jR = new JPanel[this.index.size()];
          if(Constants.b_panel){
        	  chartB = new JPanel[bwt.data.length];
        	   chartBNew = new ChartPanel[bwt.data.length];
        	    getJFrame(chartB, jB, false,dim2, BoxLayout.Y_AXIS);
              
               
          }
           if(Constants.r_panel) {
        	   chartR= new JPanel[bwt.data.length];
        	   chartRNew= new ChartPanel[bwt.data.length];
        	   getJFrame(chartR, jR,false,dim2, BoxLayout.Y_AXIS);
           }
            
          
            String[] pr = SignificancePlot.probesToPlot;
            if(Constants.showScatter()){
            for(int k=0; k<chartG.length; k++){
             chartG[k] = new JPanel[index.size()][];
             chartGNew[k] = new ChartPanel[index.size()][];
            }
            }
             try{
             locP = updateLocPanel(SignificancePlot.snp_alias);
             }catch(Exception exc){
            	 System.err.println(exc.getMessage());
             }
             if(Constants.showScatter()){
                 for(int ii=0; ii<chartGNew[0].length; ii++){
                	 chartGNew[0][ii] = new ChartPanel[pr==null ? loc.size() : pr.length];
                	 chartG[0][ii] = new JPanel[pr==null ? loc.size() : pr.length];
                 }
                 for(int ii=0; ii<chartGNew[1].length; ii++){
                	 chartGNew[1][ii] = new ChartPanel[2];
                	 chartG[1][ii] = new JPanel[2];
                 }
                 jG[0] = getJFrame(chartG[0],  true,dim5);
                 jG[1] = getJFrame(chartG[1],  true,dim5);
        		//Dime
        	//jG.setMinimumSize(dim5);
        	//jG.setSize(dim5);
        
          //   for(int k=0; k<jGS.length; k++){
             jGS[0] = new JScrollPane(
            		 new JSplitPane(JSplitPane.VERTICAL_SPLIT, locP, jG[0]),
            		  JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
             jGS[0].setMinimumSize(dimG);
             jGS[1] =// new JScrollPane(
            //	 new JScrollPane(
            			 jG[1]
            	//                    ,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER)
             ;
             jGS[1].setMinimumSize(dimG);
             }
            // }
             //jG = new JPanel[this.index.size()];
             
             jpBS = //new JSplitPane(JSplitPane.VERTICAL_SPLIT,  
            	 Constants.b_panel() ? 
            		 getSplitPane(jB,this.names1,JSplitPane.VERTICAL_SPLIT) : null; 
            		// jpDist[1]);
           
             jpRS =
            	 Constants.r_panel() ? 
            	 // new JSplitPane(JSplitPane.VERTICAL_SPLIT, 
            		 getSplitPane(jR,names1,JSplitPane.VERTICAL_SPLIT) : null;//,jpDist[0]);
          //   jB.setMinimumSize(dim2);
            if(jGS[0]!=null) jGS[0].setMinimumSize(dim2);
            if(jGS[1]!=null) jGS[1].setMinimumSize(dim2);
             //  jR.setMinimumSize(dim2);
             if(jpBS!=null) jpBS.setMinimumSize(dim2);
            
            // jpBS.setDividerLocation(400);
            // jpRS.setDividerLocation(400);
             jpB = jpBS;//,"b allele", 2*width);
             //if(probesToPlot.length>0){
        	 this.addTab("Scatter", jGS[0]);
        	 this.addTab("Scatter", jGS[1]);
        // }
             if(jpRS!=null) this.addTab("log R", jpRS);
             if(jpBS!=null) this.addTab("B allele", jpBS);
            
           
           
         
             this.setSelectedIndex(0);
            // setLayout(new BorderLayout());
          //   add(BorderLayout.CENTER, jpBS);
          //   add(BorderLayout.AFTER_LAST_LINE, jpDist);
        }
       /* else{
        	for(int i=0; i<jB.length; i++){
        	if(jpBS!=null) {
        		jB[i] = new JPanel();
        		 jB[i].setLayout(new BorderLayout());
        	}
        	if(jpRS!=null){
        		jR [i]= new JPanel();
       	  		jR[i].setLayout(new BorderLayout());
        	}
        	}
        	jG = new JPanel();
        	 
        	  jG.setLayout(new BorderLayout());
        }*/
        for(int i=0; i<include.length; i++){
        	if(include(i)) this.indicesToInclude.add(i);
        //	if(include(i))numToInclude++;
        }
        if(Constants.plot()==2) this.setSelectedIndex(Constants.dataPanelToPlot());
        if(Constants.plot()>=1 && Constants.printPlots()){
        	int len = (int) Math.floor((double)indicesToInclude.size()/(double)this.numPerPage);
        //	this.gB_ = jpBS==null ?  null : new AbstractVectorGraphicsIO[len+1];
        //	this.gR_ = jpRS==null ? null : new AbstractVectorGraphicsIO[len+1];
             	 try{
             	 
             		 File out = new File( "scatter_"+ploidy+".png");
             	     chartDistF.mkdir();
             	 //   Properties p = new Properties();
         		//	p.setProperty("PageSize", "A5");
         			//Dimension d = new Dimension(
         				//	Constants.scatterWidth()*snp_alias.length,
         					//Constants.scatterWidth()*this.snp_alias[0].length);
        // 			 g_G =new PDFGraphics2D(new File(this.chartDistF, "scatter.pdf"),d);
         			//	new ImageGraphics2D(new File(this.chartDistF, "scat_"+i2+"_ "+ii+".png"),dim5, ImageConstants.PNG); 
         //			g_G.setProperties(p);
         		//	AffineTransform at1 = new AffineTransform();
         		//	at1.setToTranslation(0, 300);
         			
         			//g_G1.setDeviceIndependent(Constants.plot()==1);
         	//		g_G.startExport();
             		
             	 }catch(Exception exc){
             		 exc.printStackTrace();
             	 }
        }
       
    }
    
   

   

	public static JComponent getFixedPane(JComponent[] jg2, String [] names, int split, Dimension dim1, Dimension dim2){
		JPanel res = new JPanel();
		res.setLayout(new BoxLayout(res, split==0 ? BoxLayout.Y_AXIS : BoxLayout.X_AXIS));
		for(int i=0; i<names.length; i++){
			JPanel jres = new JPanel();
    		jres.setToolTipText(names[0]);
    		jres.setName(names[0]);
    	/*	JComponent jsp = split==JSplitPane.VERTICAL_SPLIT ? new JScrollPane(
    				jg2[0],JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED) : jg2[0];
        				
        			*/	
    		jg2[i].setSize(dim1);
    		jg2[i].setPreferredSize(dim1);
    		jg2[i].setMinimumSize(dim1);
    		jres.add(BorderLayout.CENTER,jg2[i]);
    		res.add(jres);
		}
		res.setSize(dim2);
		res.setMinimumSize(dim2);
		return res;
	}
    
    public static  JComponent getSplitPane(JComponent[] jg2, String[] names, int split) {
   
		if(jg2.length==1) {
    		JPanel jres = new JPanel();
    		jres.setLayout(new BorderLayout());
    		jres.setToolTipText(names[0]);
    		jres.setName(names[0]);
    		JComponent jsp = split==JSplitPane.VERTICAL_SPLIT ? new JScrollPane(
    				jg2[0],JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED) : jg2[0];
        				
        				
    		jres.add(BorderLayout.CENTER,jsp
    				
    		);
    		JPanel jtp = new NamePanel(names[0]);
    		jtp.setSize(jg2[0].size());
    	//	jtp.setText(names[0]);
    		jres.add(BorderLayout.BEFORE_FIRST_LINE,jtp);
    		return jres;
    	 	}
    	else{
    		int mid =(int) Math.ceil(((double)jg2.length)/2.0);
    		JComponent[] left = new JComponent[mid];
    		JComponent[] right = new JComponent[jg2.length-mid];
    		System.arraycopy(jg2, 0, left, 0, mid);
    		System.arraycopy(jg2, mid, right, 0, right.length);
    		String[] leftnames = new String[left.length];
    		String[] rightnames = new String[right.length];
    		System.arraycopy(names, 0, leftnames, 0, mid);
    		System.arraycopy(names, mid, rightnames, 0, right.length);
    		JSplitPane jpBS =    new JSplitPane1(split,  
    				getSplitPane(left, leftnames,split), getSplitPane(right,rightnames,split));
    		
    		jpBS.setOneTouchExpandable(true);
    		  jpBS.setDividerSize(divsize);
    		  jpBS.setDividerLocation((double)left.length / (double)(left.length+right.length));
    		  return jpBS;
    	}
    }
    
  /*  public JComponent[] getBottomPane(JComponent[][] chartDist){
    	  JTabbedPane pane = new JTabbedPane();
    	  JTabbedPane pane1 = new JTabbedPane();
        for(int i=0; i<chartDist.length; i++){
           pane.add( Constants.format()[index.get(i)],getBottomPane(chartDist[i]));
           pane1.add( Constants.format()[index.get(i)],chartDist[i][2]);
        }
    
        return new JComponent[]{pane, pane1};
    }*/
   /* public JComponent getBottomPane(JComponent[] chartDist){
        for(int i=0; i<chartDist.length; i++){
            chartDist[i] = new JPanel();
           chartDist[i].setMinimumSize(dim);
        }
        
        JComponent tabs =  getTabbedPane(
                new JScrollPane(chartDist[0], JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED),
                new JScrollPane(chartDist[1], JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED)
                          );
     
        tabs.setMinimumSize(dim);
        return tabs;
    }
    public JComponent getTabbedPane(JComponent r_pair, JComponent r_single){
        JTabbedPane pane = new JTabbedPane();
        pane.addTab("r_pair", r_pair);
        pane.addTab("r_single", r_single);
        return pane;
    }*/
   
    
    
  
    int plot = 0;
    public void setToPlot(int i){
       plot= i;
    }
   public void writeCurrentCharts(int i) {
       try{
           if(!chartBF.exists()) chartBF.mkdir();
           if(!chartBF.exists()) chartRF.mkdir();
           JComponent  currentRPanel = this.getCurrentRPanel(i);
           JComponent  currentBPanel = this.getCurrentBPanel(i);
           if(currentRPanel==null) return ;
         if(currentRPanel!=null)  this.writeToZipFile(currentBPanel, chartBF, this.bwt.data[i].getName());
         if(currentBPanel!=null)  this.writeToZipFile(currentRPanel, chartRF, this.bwt.data[i].getName());
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
   
   public void writeToZipFile(Component charts, File dir, String id) throws Exception{
        File out = new File(dir, (id)+".png");
        SVGGraphics2D g = new SVGGraphics2D(out,charts);//ImageConstants.PNG); 
        g.setDeviceIndependent(Constants.plot()==1);
        g.startExport();
        charts.print(g);
        g.endExport();
   }
    
    public void writeToZipFile(Component[] charts, File dir, String[] id) throws Exception{
       for(int i=0; i<charts.length; i++){
           File out = new File(dir, (id[i])+".png");
           charts[i].setSize(500,500);
           ImageGraphics2D g  = new ImageGraphics2D(out,charts[i], ImageConstants.PNG); 
           g.setDeviceIndependent(Constants.plot()==1);
           //GIFGraphics2D g = new GIFGraphics2D(out,charts[i]); 
           g.startExport();
           charts[i].print(g);
           g.endExport();
       }
    }
   
  
  /* static JFrame getJFrame(JComponent jpBS, String name, int width){
       JFrame jf = new JFrame(name);
       jf.getContentPane().add(jpBS);
       jf.setSize(width,800);
       jf.setDefaultCloseOperation(jf.EXIT_ON_CLOSE);
       return jf;
   }*/
    final int ploidy;
    public boolean include(int k){
    	if(true) return true;
        if(include!=null){
            if(include[k]==null){
            	include[k] = false;
            	Set<Short> s = new HashSet();
            	((HaplotypeEmissionState)this.bwt.data[k]).getDataIndices(s);
            	for(int i=0; i<index.size(); i++){
            	
            	include[k] =include[k] || s.contains(new Short((short)this.index.get(i).intValue()));
            	}
            	include[k] = include[k]&& this.bwt.data[k].noCop()==ploidy;
            	if(Constants.annotateSamples!=null && !Constants.annotateSamples()[0].equals("all") && ! Arrays.asList(Constants.annotateSamples()).contains(bwt.data[k].name)){
            		include[k] = false;
            	}
            }
            return include[k];
        }
        else{
        	Set<Short> s = new HashSet();
        	((HaplotypeEmissionState)this.bwt.data[k]).getDataIndices(s);
        	boolean include = false;
        	for(int i=0; i<index.size(); i++){
        		include = include || ((HaplotypeEmissionState)this.bwt.data[k]).getDataIndices(s).contains(new Short((short)this.index.get(i).intValue()));
        	}
        	include = include && this.bwt.data[k].noCop()==ploidy;
        	if(Constants.annotateSamples!=null && !Constants.annotateSamples()[0].equals("all") && ! Arrays.asList(Constants.annotateSamples()).contains(bwt.data[k].name)){
        		include = false;
        	}
        		
        	return include;
        }
    	
        
    }
    Boolean[] include;
    JComponent[] jp;
    JComponent getJFrame( JComponent[][] cp, boolean all, Dimension dim){
    	 jp= new JComponent[cp.length];
    	 Dimension dim1 = new Dimension(cp[0].length*dim.width,dim.height);
    	for(int i=0; i<cp.length; i++){
    		jp[i] =getJFrame(cp[i], all,dim,dim1);
    	}
    	Dimension dim2 = new Dimension(cp[0].length*dim.width,(dim.height)*cp.length);
    	JComponent jcomp =   this.getFixedPane(jp,names1, JSplitPane.VERTICAL_SPLIT,dim1,dim2);
    	jcomp.setSize(dim2);
    	return jcomp;
    }
    JComponent[] getJFrame( JComponent[] cp, JComponent[] jpBS, boolean all,
    		Dimension dim1, int axis
    ){
      for(int i=0; i<jpBS.length; i++){
    	  jpBS[i] = new JPanel();
    	  jpBS[i].setLayout(new BoxLayout(jpBS[i], axis));
    	  jpBS[i].setSize(dim1);
    	  jpBS[i].setMinimumSize(dim1);
      }
       
        for(int i=0; i<cp.length; i++){
            if(all || include(i)){
            cp[i] = new JPanel();
            jpBS[((HaplotypeEmissionState)this.bwt.data[i]).dataIndex()].add(cp[i]);
            }
        }
        return jpBS;
    }
    
    JComponent getJFrame( JComponent[] cp, boolean all, Dimension dim, Dimension dim1){
        JPanel jpBS = new JPanel();
    
     	
     	  jpBS.setLayout(new BoxLayout(jpBS, BoxLayout.X_AXIS));
       
        
    //    int height =  dim.height;
         for(int i=0; i<cp.length; i++){
             if(all || include(i)){
             cp[i] = new JPanel();
             jpBS.add(cp[i]);
             cp[i].setSize(dim);
             cp[i].setMinimumSize(dim);
             cp[i].setPreferredSize(dim);
           // cp[i].setMinimumSize(new Dimension(dim.width,height*cp.length ));
                // jpBS.add(cp[i]);
             }
         }
         jpBS.setSize(dim1);
        jpBS.setPreferredSize(dim1);
        jpBS.setMinimumSize(dim1);
        
    //    this.pack();
         
         return jpBS;
     }
    
  
    public double[] noSamples;
    
    
   /* public ChartPanel[] getDistCharts(int index){
    	if(noSamples==null){
    		noSamples = new double[bwt.probDists[index].mvf.single.size()];
    	}
        ChartPanel[] cp = new ChartPanel[3];
        bwt.probDists[index].mvf.updateSampleSize(noSamples);
          cp[0] = 
        	  Constants.format()[index].startsWith("geno") ? null : 
        	  plotR(bwt.probDists[index].s1, bwt.probDists[index].mvf, null,index,  0, this.lu_rpair, "Conditional distributions over LRR", "LRR", bwt.probDists[index].ca_r);
        
          cp[1] = plotR(bwt.probDists[index].mvf.single, null, 
        		 noSamples,index, 1, lu_r, "Conditional distributions over LRR", "LRR",  bwt.probDists[index].ca_r);
          cp[2] = plotR(bwt.probDists[index].s2, null, null,index,2, lu_b, "Conditional distributions over BAF","BAF", bwt.probDists[index].ca_b);
          
          return cp;
    }*/
    
  
    public XYSeriesCollection refresh(XYSeriesCollection xys){
    	XYSeriesCollection series = new XYSeriesCollection();
    	 for(int k=0; k<xys.getSeriesCount(); k++){
			 XYSeries series_  =   xys.getSeries(k);
			 series.addSeries(new XYSeriesMiss(series_.getKey()));
	            
	        }
    	 return series;
    }
    public void reinitialise(int l) {
      //  this.currentIndividual = l;
    	if(l==0){
    		for(int i=0; i<rb.length; i++){
    		
    		for(Iterator<Integer> it = this.rb[i].keySet().iterator(); it.hasNext();){
    			Integer v = it.next();
    			 rb[i].put(v, this.refresh(rb[i].get(v)));
    			
    		}
    		for(Iterator<AnnotationSeries> it = this.annotation[i].values().iterator(); it.hasNext();){
    			AnnotationSeries xys = it.next();
    			xys.clear();
    		}
    		}
    		for(int i=0; i<this.rGlobal.length; i++){
    			for(int k=0; k<rGlobal[i].length; k++){
    				this.rGlobal[i][k] = this.refresh(rGlobal[i][k]);
    			}
    		}
    	}
       XYSeriesCollection[] current_b = this.getBSeriesCollection(l);
       XYSeriesCollection[] current_r  = this.getRSeriesCollection(l);
       if(current_r!=null){
       for(int j=0; j<current_r.length; j++){
        for(int k=0; k<current_r[j].getSeriesCount(); k++){
            current_r[j].getSeries(k).clear();
        }
       }
       }
       if(current_b!=null){
       for(int j=0; j<current_b.length; j++){
        for(int k=0; k<current_b[j].getSeriesCount(); k++){
            current_b[j].getSeries(k).clear();
        
        }
       }
       }
        
    }
    
    private XYSeriesCollection newRBSeries(){
    	XYSeriesCollection coll = new XYSeriesCollection();
    	for(int k=0; k<this.rbSeriesNames.size(); k++){
        	
        	
    		coll.addSeries(new XYSeriesMiss(this.rbSeriesNames.get(k)));
    	
    	}
    	coll.addSeries(new XYSeriesMiss("unclassified"));
    	return coll;
    }
 
    private XYSeriesCollection getRBSeries(int ind, int i) {
    	int alias =  SignificancePlot.snp_alias(i);
    	XYSeriesCollection coll = rb[ind].get(alias);
    	if(coll!=null) return coll;
    	else{
    		rb[ind].put(i, coll =  newRBSeries());
    		this.annotation[ind].put(i, new AnnotationSeries(this.emStSp.getColor(true)));
    	
    	return coll;
    	}
    
    }
    final int[] data_index_alias; //converts from index to position in index
   final Map<Integer, XYSeriesCollection>[] rb;
   final Map<Integer,AnnotationSeries>[] annotation;
   
   static class AnnotationSeries{
	   Map<String, XYTextAnnotation> l = new HashMap<String, XYTextAnnotation>();
	   final Color[] col;
	   AnnotationSeries(Color[] col){
		   this.col = col;
	   }
	 //  Color[] col =  DistributionCollection.dc.getColR();
	public void clear() {
		this.l.clear();
		
	}

	public void put(String st, XYDataItem val, int i) {
		XYTextAnnotation annot = new XYTextAnnotation(st,val.getXValue(), val.getYValue());
		l.put(st, annot);
		annot.setFont(font5);
	//	if(i<0) annot.setPaint(Color.LIGHT_GRAY);
	//	else annot.setPaint(col[i]);
		annot.setPaint(Color.RED);
		
	}

	public void addAnnotation(XYPlot xyp) {
		for(Iterator<XYTextAnnotation> it = l.values().iterator(); it.hasNext();){
			XYTextAnnotation xyt = it.next();
		xyt.setPaint(Color.RED);
			xyp.addAnnotation(xyt);	
		}
		
	}
   }
 
   final Set<Integer> genoToAnnotate;
   final Set<String> samplesToAnnotate;
   // short[] ind_ = new short[1];
   
   XYSeriesCollection[][] rGlobal;// rGlobalProbeOnly;
   double[] pos = new double[2];
   /** i is index of individual
     * stateDistribution has the distribution over genotypes
     *  */
    public void addedInformation(StateDistribution emissionC, int ll, int i, PseudoDistribution dist, int k, HaplotypeEmissionState sta){
        double x =location.get(i).doubleValue();
        if(singlechrom) x = Constants.decode(x, pos, this.singlechrom)/1e6;
        
        short ind =  plotMerged ? 0 : dist.getDataIndex() ;
//        	(dist instanceof CompoundDistribution) ? (short) ((CompoundDistribution)dist).getDataIndex(k):
       
        //ind =0;
        if(ind<0) ind=(short)(data_index_alias.length-1);
        String name = this.bwt.data[ll].name;
        if(ind>=data_index_alias.length || ind<0) return;
      int d_i =  data_index_alias[ind];
      if(d_i ==-1)return;
        XYSeriesCollection current_b = Constants.b_panel() ? this.getBSeriesCollection(ll)[d_i] : null;
        XYSeriesCollection current_r  = Constants.r_panel() ? this.getRSeriesCollection(ll)[d_i] : null;
      int i1 =  indexOf(SignificancePlot.snp_alias, i);
      Boolean probeOnly = dist.probeOnly();
      
      XYSeriesCollection rbig = probeOnly==null ? null : this.rGlobal[ind][probeOnly ? 1:0];// rGlobalProbeOnly[ind] : rGlobal[ind]; 
    //  if(probeOnly){
    //  }
        XYSeriesCollection rbi = i1>=0 ? rb[ind].get(i1) : null;
    //    XYSeriesCollection rbi1 = rb[ind].get(-1) ;
        AnnotationSeries annoti = i1>=0 ? this.annotation[ind].get(i1) : null;
   //     AnnotationSeries annoti1 = i1>=0 ? this.annotation[ind].get(-1) : null;
        double[] prob =  PairEmissionState.pool.getObj(Emiss.getSpaceForNoCopies(sta.noCop()).genoListSize());
        
        double sumNonMix = Constants.allowComponent() ?  Sampler.getProbOverStates(emissionC, bwt.hmm, sta,i,prob,0) : 0;
        
       
        double sum = 	Sampler.getProbOverStates(emissionC, bwt.hmm, sta, i,prob, Constants.isLogProbs());
        if(!Constants.allowComponent) sumNonMix = sum;
//        double sumNonMix = Constants.allowComponent() ? Sampler.getProbOverStates(emissionC, bwt.hmm, sta, i,prob) : sum;//,0);	
     //    double[] prob_b = new double[current_b.getSeriesCount()];
       /* if(dist instanceof IlluminaDistribution && Constants.getMax(prob)==5){
        	  IlluminaDistribution distr1 = (IlluminaDistribution) dist;
        	  if(distr1.r().doubleValue()<-3 && distr1.b() < 0.01){
        		  System.err.println("h");
        		  Sampler.getProbOverStates(emissionC, bwt.hmm, sta, i,prob);
        	  }
          }*/
       
         double[] prior_st = new double[emStSp.copyNumber.size()];
       //  double[] prob_r =new double[emStSp.copyNumber.size()];
           /* if(dist instanceof CompoundDistribution){// && ! (dist instanceof ACSNDistribution)){
                dist = ((CompoundDistribution)dist).getForIndex(ind);
            }*/
            if(dist!=null && true){//dist.containsIndex(ind)){
            	
            
             //  getProbOverRDists(emissionC,prob_r);
             /*   Arrays.fill(prob_b, 0);
              //  prob_b[Constants.getMax(prob)]=1.0;
                for(int k=0; k<prob.length; k++){
                   int alias = b_alias[k];
                   if(alias < prob_b.length){ //relies on uniform being at end
                       prob_b[alias]+=prob[k];
                   }
                }*/
                int tot_s = Constants.getMax(prob);
            //    int r_s = Constants.getMax(prob_r);
                int b_s = tot_s;//Constants.getMax(prob);
               if( dist instanceof IlluminaDistribution){
            	   Number rv = ((IlluminaRDistribution)dist).r1();
            	   if(useVals){
                   	 List<AssociationCalculator>[][]ac = DataCollection.datC.ac;
                     rv = ((LinearRegressionCalculator) this.sp.get(ind).ac
                    		 ).getValue(sta.getName(), phenoIndex);
                   }
            	   Double bv =  ((IlluminaDistribution)dist).b(i);
            	//   bv = rv.doubleValue()==0 ? Math.random() : bv/(rv.doubleValue());
            	
            	   
                	if(rv!=null && current_r!=null) current_r.getSeries(b_s).addOrUpdate(x, Math.max(rv.doubleValue(), minrv));
                	if(bv!=null && current_b!=null) current_b.getSeries(b_s).addOrUpdate(x, bv.doubleValue());
                	 if(rv!=null && bv!=null){
                		
                		XYDataItem item = 
                			new XYDataItem(rv==null ? Double.NaN : rv,  bv);
                		if(prob[tot_s] >= Constants.imputedThreshGraph(0) && sumNonMix>=Constants.imputedThreshGraph(1)*sum){
                			/* if(bv>0.95){
                    			 System.err.println("h");
                    			 double sumNonMix1 = Constants.allowComponent() ?  Sampler.getProbOverStates(emissionC, bwt.hmm, sta,i,prob,0) : 0;
                    		        
                    		       
                    		        double sum1 = 	Sampler.getProbOverStates(emissionC, bwt.hmm, sta, i,prob);
                    		 }*/
                			if(rbi!=null){
                				rbi.getSeries(tot_s).add(item);
                				if((genoToAnnotate.contains(tot_s)|| samplesToAnnotate.contains(name))){
                    				annoti.put(name, item, tot_s);
                    		//		annoti1.put(name, val, i1)
                    			}
                			}
                			if(rbig!=null){
                			//	Comparable c = this.emStSp.get(tot_s);
                			//	System.err.println(c+" "+rv+" "+bv);
                				rbig.getSeries(tot_s).add(item);
                			}
                		//	rbi1.getSeries(tot_s).add(item);
                			
                		}
                		else{
                			if(rbi!=null){
                				rbi.getSeries(prob.length).add(item);
                				if(genoToAnnotate.contains(prob.length)|| samplesToAnnotate.contains(name)){
                    				annoti.put(name, item, -1);
                    			}
                			}
                			if(rbig!=null) rbig.getSeries(prob.length).add(item);
                			
                		}
                	}
                	 /*else{
                		 for(int j=0; j<prob.length; j++){
                			 if(prob[j]>0){
                				((XYSeriesMiss) rbi.getSeries(j)).nan_x+=prob[j]/(double) prob.length;
                			 }
                		 }
                	 }*/
            		}
              /* else if( dist instanceof ACSNDistribution){
            	   Double rv = ((ACSNDistribution)dist).dist1.r();
            	   
            	   double bv =  ((ACSNDistribution)dist).dist2.r();
                	if(rv!=null) current_r.getSeries(r_s).addOrUpdate(x, rv);
                	current_b.getSeries(b_s).addOrUpdate(x, bv);
                	if(rbi!=null && rv!=null) rbi.getSeries(r_s).add(rv.doubleValue(), bv);
            		}*/
                else 
                if(dist instanceof IlluminaRDistribution){
                	   Number rv = ((IlluminaRDistribution)dist).r1();
                	   if(useVals){
                         	 List<AssociationCalculator>[][]ac = DataCollection.datC.ac;
                           rv = ((LinearRegressionCalculator) this.sp.get(ind).ac
                          		 ).getValue(sta.getName(), phenoIndex);
                         }
                	   if(rv!=null)
                		   if(current_r!=null) current_r.getSeries(b_s).addOrUpdate(x, rv.doubleValue());
                	if(rv!=null ){
                		XYDataItem item = new XYDataItem( rv.doubleValue(), Constants.rand.nextDouble()*Constants.bRandom());
                		if(rbi!=null){
                			rbi.getSeries(tot_s).add(item);
                			if(genoToAnnotate.contains(tot_s)|| samplesToAnnotate.contains(name)){
                				annoti.put(name, item, tot_s);
                		//		annoti1.put(name, val, i1)
                			}
                		}
                		if(rbig!=null) rbig.getSeries(tot_s).add(item);
                		
                	   }
                	
                		else{
                   		 for(int j=0; j<prob.length; j++){
                   			 if(prob[j]>0){
                   				((XYSeriesMiss) rbi.getSeries(j)).nan_x+=prob[j];
                   			 }
                   		 }
                   	
                	}
                	//current_b.getSeries(b_s).addOrUpdate(x,  ((IlluminaDistribution)dist).b());
                }	
                else if(Constants.b_panel() && dist instanceof IntegerDistribution){ //r_s!=bwt.probDists[index].normal_index
                	//double[] probr1 = 
                	int pos =( (IntegerDistribution)dist).fixedInteger();
                	//int cn = this.emStSp.getCN(pos); 
              	//	sta.getEmissionStateSpace().getCN(pos);
//                	 if(current_r!=null) current_r.getSeries(b_s).addOrUpdate(x,cn);
                    if(current_b!=null) current_b.getSeries(b_s).addOrUpdate(x,  emStSp.getBCount(pos)/2.0);
                }
                else if(Constants.b_panel() && dist instanceof SimpleExtendedDistribution){ //r_s!=bwt.probDists[index].normal_index
                	//double[] probr1 = 
                	
                	sum( dist.probs(), prior_st);
                	int bcount =  emStSp.getBCount(Constants.getMax(dist.probs()));
                //	if(current_r!=null)  current_r.getSeries(b_s).addOrUpdate(x,prob[cn]);
//                	 if(current_b!=null) current_b.getSeries(b_s).addOrUpdate(x,  prob[cn]);
//                     current_b.getSeries(b_s).addOrUpdate(x,  bcount/2.0);
                }
                else if(dist instanceof MatchedDistributionCollection.BackgroundDistribution){
                	//int cn = emStSp.getCN(tot_s);
                 //	current_r.getSeries(b_s).addOrUpdate(x,  ((MatchedDistributionCollection.BackgroundDistribution)dist).r1(cn));
                	current_b.getSeries(b_s).addOrUpdate(x,  ((MatchedDistributionCollection.BackgroundDistribution)dist).b1());
                }
            }
            PairEmissionState.pool.returnObj(prob);
      //  }
    }
    private int indexOf(int[] is, int p) {
    	if(is==null) return p;
	for(int i=0; i<is.length; i++){
		if(is[i]==p) return i;
		
	}
	return -1;
}

	private void getProbOverRDists(StateDistribution emissionC, double[] prob) {
	//double[] prob = new double[this.emStSp.copyNumber.size()];
	Arrays.fill(prob, 0.0);
	for(int i=1; i<emissionC.dist.length; i++){
		double v = emissionC.dist[i];
		prob[((EmissionState)this.hmm.getState(i)).noCop()]+=v;
	}
	//return prob;
}

	protected void sum(double[] emiss, double[] prior_st) {
		Arrays.fill(prior_st, 0.0);
		for(int i=0; i<emiss.length; i++){
			prior_st[emStSp.getCN(i)]+=emiss[i];
		}
		
	}

	double[] res;
   /* protected synchronized double[] getProbOverStates(StateDistribution emissionC,
            MarkovModel hmm, HaplotypeEmissionState obj, int i) {
    	  EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
      if(res==null) res = new double[emstsp.defaultList.size()];
      Arrays.fill(res,  0.0);
      double sum = 0;
      for(int j=1; j<emissionC.dist.length; j++){
          double p = emissionC.dist[j];
          EmissionState state_j = (EmissionState) hmm.getState(j);
          EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
          if(p>0){
              sum+=obj.calcDistribution((EmissionState) state_j, i,  p);
              double[] dist = state_j.distribution;
              for(int k=0; k<dist.length; k++){
            	  res[emstsp.get(emstsp1.get(k))]+=dist[k];
              }
          }
      }
      for(int k=0; k<res.length; k++){
          res[k] = res[k] / sum;
      }
        return res;
    }*/
    
    public JComponent getCurrentBPanel(int l){
        if(this.chartB!=null) return chartB[l];
        else return new JPanel();
    }
    Stroke stroke = new BasicStroke(1.0f);
    Stroke stroke1 =      new BasicStroke(Constants.strokeWidth(), BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, new float[] {1.0f, 1.0f, 1.0f}, 1.0f);
    
  
    public JComponent getCurrentGPanel(int l, int k, int ind){
        if(this.chartG[ind]!=null && chartG[ind][l][k]!=null) {
        	
        	return chartG[ind][l][k];
        }
        else {
        	JPanel jp = new JPanel();
        
        	return jp;
        }
    }
  
    public JComponent getCurrentRPanel(int l){
        if(this.chartR!=null) return chartR[l];
        else return new JPanel();
    }
    
    public static  Font font= new Font("Tahoma", Font.PLAIN, (int) Math.round(4*Constants.shapeSize1()));
    
    public void update(){
       if(chartR!=null){ for(int i=0; i<chartR.length; i++){
            if(chartRNew[i]==null) continue;
            JComponent currentRPanel = this.getCurrentRPanel(i);
            JComponent currentBPanel = this.getCurrentBPanel(i);
         
            if(currentRPanel==null) continue;
            currentBPanel.setMinimumSize(dim);
            currentRPanel.setMinimumSize(dim);
           currentBPanel.setSize(dim);
            currentRPanel.setSize(dim);
              currentRPanel.removeAll(); 
             currentBPanel.removeAll();
             
             currentRPanel.add(chartRNew[i]);
             currentBPanel.add(chartBNew[i]);
        }
       }
      
       if(chartG!=null){
    	   for(int kk=0; kk<chartG.length; kk++){
    		   if(chartG[kk]!=null){
        for(int i=0; i<chartG[kk].length; i++){
        	for(int ik=0; ik<chartG[kk][i].length; ik++){
	        	if(chartGNew[kk][i][ik]!=null){
		        	JComponent currentGPanel = this.getCurrentGPanel(i,ik, kk);
		        	
		        	if(currentGPanel!=null) {
		        		
		        		
		        		currentGPanel.removeAll();
		        		currentGPanel.add(chartGNew[kk][i][ik]);
		        		currentGPanel.setMinimumSize(dim5);
		        		currentGPanel.setSize(dim5);
		        	}
	        	}
        	}
        }
    		   }
    	   }
       }
    }
    public static Stroke dashed =
        new BasicStroke(0.2f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, new float[] {10.0f, 1.0f, 10.0f}, 1.0f);
    public static Stroke dotted =
        new BasicStroke(0.8f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, new float[] {2.0f, 20.0f}, 0.0f);
 
   /* public void check(XYSeriesCollection datas){
        Map<Double, Double> s = new HashMap<Double, Double>();
        for(int i=0; i<datas.getSeriesCount(); i++){
            XYSeries series_i = datas.getSeries(i);
            for(int k=0; k<series_i.getItemCount(); k++){
                double x = series_i.getX(index).doubleValue();
                double y = series_i.getY(index).doubleValue();
                Double y_1 = s.get(x);
                if(y_1!=null) throw new RuntimeException("already in !"+x+" "+y+" "+y_1);
               s.put(x, y);
               
            }
        }
    }*/
    
    public Shape getShape(char ch, double sz){
    	Graphics2D gr = (Graphics2D)this.getGraphics();
     Shape sh = this.font10.createGlyphVector(gr.getFontRenderContext(), new char[] {ch}).getOutline();//new Ellipse2D.Double(sz,sz,sz,sz);//new Rectangle(1,1);
     AffineTransform at = new AffineTransform();
     double sz1 = sz/sh.getBounds().getHeight();
      at.setToScale(sz1, sz1);
      Shape sh1 = at.createTransformedShape(sh);
     /* double diffx = sh.getBounds2D().getMaxX() - sh1.getBounds2D().getMaxX();
      double diffy = sh.getBounds2D().getMaxY() - sh1.getBounds2D().getMaxY();
      at.setToTranslation(diffx, diffy);
      return at.createTransformedShape(sh1);*/
      return sh1;
       // 
       // return 
    }
    
    Shape ellipse = new Ellipse2D.Double(1,1,1,1);
	private String chrom;
    public static Shape getShape(Shape sh, double sz){
    //  Shape sh = //new Rectangle(1,1);
     AffineTransform at = new AffineTransform();
    double sz1 = sz/ sh.getBounds().getHeight();
    if(Math.abs(sz1-1)<0.01) return sh; 
    at.setToScale(sz1, sz1);
    Shape sh1 = at.createTransformedShape(sh);
  /*  double diffx = sh.getBounds2D().getMaxX() - sh1.getBounds2D().getMaxX();
    double diffy = sh.getBounds2D().getMaxY() - sh1.getBounds2D().getMaxY();
    at.setToTranslation(diffx, diffy);
    
  
      return at.createTransformedShape(sh1);*/
    return sh1;
       // Graphics2D gr = (Graphics2D)this.getGraphics();
       // return this.font10.createGlyphVector(gr.getFontRenderContext(), new char[] {ch}).getOutline();
    }
  // static char[] shapes1 = new char[] {'x', 'o', '+', '#','!','=','+','=', '%','*','a','b','c','d','e','f'};
   public static Shape[] shapes = new Shape[Math.max(7, Constants.maxPloidy()*Constants.maxCopies()+2)];
   public static Shape shape_null;
 
   static{
	   for(int k=0;k<shapes.length; k++){
		   int len =  DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE.length;
		shapes[k] = getShape(   DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE[k%len],Constants.shapeSize());
	   }
	   shape_null = getShape(new Ellipse2D.Double(1,1,1,1), Constants.shapeSize()/3.0);
		  
	   AffineTransform rot = new AffineTransform();
	   rot.setToRotation(Math.PI/4);
	   shapes[0] = getShape(new Rectangle2D.Double(0,0,2,2), Constants.shapeSize());
	   shapes[4] =getShape( new Polygon(new int[] {-1,1,0}, new int[] {-1,-1,1},3), Constants.shapeSize());
	   shapes[2] = getShape(new Polygon(new int[] {-1,1,0}, new int[] {1,1,-1},3), Constants.shapeSize());
	   shapes[3] = getShape(rot.createTransformedShape(shapes[0]), Constants.shapeSize());
	   shapes[1] = getShape(new Ellipse2D.Double(1,1,1,1), Constants.shapeSize());
	   shapes[5] = getShape(new Polygon(new int[] {-2,0,2,0}, new int[] {0,-1,0,1},4), Constants.shapeSize());
	   shapes[6] = getShape(rot.createTransformedShape(shapes[5]), Constants.shapeSize());
	   Shape sh = getShape(new Polygon(new int[] {-2,-2,0,2,2,2,0,-2}, new int[] {0,-1,-1,-1, 0,1,1,1},8), Constants.shapeSize());
	  for(int i=7; i<shapes.length; i++){
	     if(i==7) shapes[i] = sh;
	     else{
	    	 shapes[i] = getShape(rot.createTransformedShape(shapes[i-1]), Constants.shapeSize());
	     }
	  }
	  // shapes[5] = 
   }
/*<<<<<<< .mine
    public static  Font font8 = new Font("SansSerif", Font.PLAIN, 30);
    public static  Font font5 = new Font("SansSerif", Font.PLAIN, 30);
    public static  Font font6 = new Font("SansSerif", Font.PLAIN, 24);
=======*/
    public static  Font font8 = new Font("SansSerif", Font.PLAIN, Constants.shapeMult()*(int) Constants.shapeSize());
    public static  Font font5 = new Font("SansSerif", Font.PLAIN, Constants.shapeMult()*(int) Constants.shapeSize());
    public static  Font font6 = new Font("SansSerif", Font.PLAIN, Constants.shapeMult()*(int) Constants.shapeSize());
//>>>>>>> .r286
    public  JFreeChart graph(XYSeriesCollection[] datas,  String title, Color[] ca, Shape[]sha, boolean b ) {
    	AbstractXYItemRenderer[] renderer = new AbstractXYItemRenderer[datas.length];
    	Shape[] shape = new Shape[datas.length];
    	Double[] shapeSize = new Double[datas.length];
    	//Arrays.fill(shapeSize,2);
    	for(int i=0; i<datas.length; i++){
    		String[] str = Constants.plotType(i);
    		try{
    		renderer[i] = (AbstractXYItemRenderer)Class.forName("org.jfree.chart.renderer.xy."+str[0]).getConstructor(new Class[0]).newInstance(new Object[0]);
    		}catch(Exception exc){
    			exc.printStackTrace();
    			
    		}
    		if(renderer[i] instanceof XYLineAndShapeRenderer){
    			XYLineAndShapeRenderer rend = (XYLineAndShapeRenderer)renderer[i];
	    		if(str[1].equals("nodot")){
	    			rend.setBaseShapesVisible(false);
	    		}
	    		else{
	    			rend.setBaseShapesVisible(true);
	    			if(str[1].length()==1){
	    				shape[i] = getShape(ellipse,Constants.shapeSize());
	    			}
	    			else {
	    			//	ShapeUtilities.
	    				shape[i] = getShape(new Rectangle(2,2),Constants.shapeSize());
	    			}
	    			
	    		}
	    		if(str[2].equals("noline")){
	    			rend.setBaseLinesVisible(false);
	    		}
	    		else{
	    			rend.setBaseLinesVisible(true);
	    		}
	    		if(str[1].startsWith("circ")){
	    	//		shape[i] = new Ellipse2D.Float(20,20,20,20);
	    			rend.setBaseShape(shape[i]);
	    			rend.setShapesFilled(false);
	    			
	    		}
	    		if(str.length>3){
	    			shapeSize[i] = Double.parseDouble(str[3]);
	    		}
    		}
	   else{
    			XYBarRenderer rend = (XYBarRenderer) renderer[i];
    			rend.setDrawBarOutline(true);
    			rend.setBase(str.length>1 ? Double.parseDouble(str[1]): 0);
    		}
    	}
        final JFreeChart chart = ChartFactory.createXYLineChart(
        		title, //Constants.experiment(),
//               b && Constants.plot!=2? "": title,
               !b ? "": "Position on chromosome "+chrom + (singlechrom ? " (mb)" :""), // domain axis label
               b ?  "Proportional read depth" : "Total read depth", //"HD statistic", // range axis label
                datas[0], // data
                PlotOrientation.VERTICAL, Constants.includeLegend(), // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
      
        XYPlot plot = (XYPlot) chart.getPlot();
        for(int i=0; i<datas.length; i++){
        	  plot.setDataset(i,datas[i]);
        }
      
      
     
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderPaint(Color.white);
        plot.setBackgroundPaint(Color.white);
            NumberAxis yAxis = (NumberAxis)plot.getRangeAxis();
            NumberAxis xAxis = (NumberAxis)plot.getDomainAxis();
            chart.getTitle().setFont(font8);
            yAxis.setTickLabelFont(font6);
            yAxis.setLabelFont(font6);
            xAxis.setTickLabelFont(font6);
            xAxis.setLabelFont(font6);
            if(location.size()>0 && false){
	            xAxis.setAutoRange(false);
	           xAxis.setLowerBound(location.get(0)-1000);//, pos, this.singlechrom));
	            xAxis.setUpperBound(location.get(location.size()-1)+1000);//, pos, this.singlechrom));
            }
         
            
            
         
            plot.setDomainGridlinePaint(Color.WHITE);
           plot.setRangeGridlinePaint(Color.WHITE);
        
            for(int kk=0; kk<datas.length; kk++){
           // int ik=0;
            plot.setRenderer(kk, renderer[kk]);
            for(int i=0; i<datas[kk].getSeriesCount(); i++){
            	int noCop =emStSp.getCN(i);
                Color colors = ca[noCop];
               int noB = emStSp.getBCount(i);
               // Color colors = ca[i];
                
                renderer[kk].setSeriesPaint(i, colors);
                renderer[kk].setSeriesFillPaint(i,colors);
              if(!(renderer[kk] instanceof XYBarRenderer) ) renderer[kk].setSeriesShape(i,
            		 shapeSize[kk]==null   ? sha[noB]: new Ellipse2D.Float(-10f,-10f,20f,20f));
            			 //this.getShape(sha[noB],shapeSize[kk] ));
              
              //  ik++;
            }
            }
        
  
              chart.setBackgroundPaint(Color.white);
             return chart;
}
    
  ShapeList shapeList = new ShapeList();
   

  
  public XYSeriesCollection getDistr(XYSeriesCollection datas, boolean x){
	 List<Double>[] li = new List[emStSp.copyNumber.size()];
	 for(int i=0; i<li.length; i++){
		 li[i] = new ArrayList<Double>();
	 }
	  for(int i=0; i<datas.getSeriesCount(); i++){
		 int cn =  this.emStSp.getCN(i);
		 int j = this.emStSp.getGenoForCopyNo(cn)[0];
		List<Double> d = li[j];
		  XYSeries series = datas.getSeries(i);
		  for(int k=0; k<series.getItemCount(); k++){
			  XYDataItem xy = series.getDataItem(k);
			  if(x){
				 d.add( x? xy.getXValue() : xy.getYValue());
			  }
		  }
	  }
	  XYSeriesCollection datas_new = new XYSeriesCollection();
	  for(int j=0; j<li.length; j++){
		  List<Double> d = li[j];
		   Collections.sort(d);
		   XYSeries series_new = new XYSeries(emStSp.copyNumber.get(j)+"");
		   for(int k=0; k<d.size(); k++){
			   series_new.add(d.get(k),new Double(((double)k+0.5)/(double)d.size()));
		   }
		  datas_new.addSeries(series_new);
	  }
	  return datas_new;
  }
  ///this
  public  JFreeChart graph(XYSeriesCollection datas_,  String title, Color[] ca , int data_index, Double fracR, Double fracB) {
      //  this.shapeList.setShape(index, shape)
	  XYSeriesCollection datas = this.useVals ? getDistr(datas_,true) : datas_;
    	if(Constants.scatterWidth<300){
    		font8 = new Font("SansSerif", Font.PLAIN, 1*(int) Constants.shapeSize1());
    		font6 = new Font("SansSerif", Font.PLAIN, 1*(int) Constants.shapeSize1());
    		font4 = new Font("SansSerif", Font.PLAIN, 1*(int) Constants.shapeSize1());
    	}
        final JFreeChart chart = ChartFactory.createXYLineChart(
                title+" chromosome "+chrom+
                (Constants.collapseScatterInd() ? "" : 
                (fracR==null ? "" : 
                	String.format(" Q= %5.3g", new Object[] {-Math.log10(1-fracR)}))),
                useVals? "Phenotype quantile" : "LRR ", // domain axis label
                useVals? "Cumulative distribution within CN band" : "BAF ", // range axis label
                datas, // data
                PlotOrientation.VERTICAL, Constants.includeLegend(), // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
        XYPlot plot = (XYPlot) chart.getPlot();
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderPaint(Color.white);
        plot.setBackgroundPaint(Color.white);
     
        if(fracR!=null){
        	Color c = HMMPanel.modify(Color.white, 
             		Math.max(Math.pow(fracR*fracB,Constants.backgroundPower()), 0.05));
        	chart.setBackgroundPaint( c);
        	plot.setBackgroundPaint(c);
        	plot.setForegroundAlpha((float)	Math.max(Math.pow(fracR*fracB,Constants.foregroundPower()), 0.05));
        }
            NumberAxis yAxis = (NumberAxis)plot.getRangeAxis();
          
            NumberAxis xAxis = (NumberAxis)plot.getDomainAxis();
           
            chart.getTitle().setFont(font8);
            yAxis.setTickLabelFont(font6);
            yAxis.setLabelFont(font6);
            xAxis.setTickLabelFont(font6);
            xAxis.setLabelFont(font6);
         
            plot.setDomainGridlinePaint(Color.WHITE);
           plot.setRangeGridlinePaint(Color.WHITE);
            final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
            renderer.setShapesFilled(Boolean.FALSE);
            renderer.setLinesVisible(false);
            final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
            final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
            renderer1.setBaseShape(new Rectangle(1,1));
          
           if(dotted!=null) renderer1.setBaseStroke(dotted);
            
            int ik=0;
            int len = datas.getSeriesCount();
           
            for(int i=0; i<len; i++){
            
            	
               // if(ca[noCop]==Color.cyan &&! datas.getSeries(i).isEmpty()){
             //   }
            	 Color colors;
            	 int noB;
                if(i==len-1){
                	colors = Color.LIGHT_GRAY;
                	  noB = 0;
                }
                else{
                	int noCop =emStSp.getCN(ik);
                    colors = ca[noCop];
                    noB = emStSp.getBCount(i);
                }
             
                renderer.setSeriesPaint(i, colors);
              //  Sh
                renderer.setSeriesShape(i, i==len-1 ? shape_null : shapes[noB]);//getShape(shapes[noB],noCop==2 ? 2.0 : 4.0));
                renderer.setSeriesShapesFilled(i, false);
                renderer1.setSeriesShapesFilled(i, false);
                renderer1.setSeriesFillPaint(i, this.noColor);
                renderer2.setSeriesShapesFilled(i, false);
                renderer1.setSeriesPaint(i, colors);
                renderer2.setSeriesPaint(i, colors);
                ik++;
                
               
                
            }
          //renderer1.setBaseShape(shape)
            renderer1.setShapesVisible(false);
            renderer2.setShapesVisible(false);
            renderer2.setLinesVisible(true);
            renderer2.setLinesVisible(true);
            //renderer1.setStroke(new BasicStroke(BasicStroke.JOIN_MITER));
            renderer1.setStroke(dashed);
          if(dotted!=null)  renderer1.setStroke(dotted);
            plot.setRenderer(0,renderer);
          if(false){
            try{
                NumberAxis yAxis1 = (NumberAxis) yAxis.clone();
                NumberAxis xAxis1 = (NumberAxis) xAxis.clone();
                plot.setDomainAxis(1, xAxis1);
                yAxis1.setTickLabelFont(null);
                //xAxis.
                plot.setRangeAxis(1, yAxis1);
                }catch(CloneNotSupportedException exc){exc.printStackTrace();}
          }
           //   chart.setBackgroundPaint(Color.white);
             return chart;
}
    Color noColor = new Color(1,1,1,0);
    
  /*  static  JFreeChart getChart(XYSeriesCollection exp, 
             XYSeriesCollection obs, String name, String yaxis, ColorAdapter ca,
             Map<Double[], String[]> sn_params
            ) {
     
         final JFreeChart chart = ChartFactory.createXYLineChart(
                 name,
                 yaxis, // domain axis label
                 "Frequency expected", // range axis label
                 obs, // data
                 PlotOrientation.HORIZONTAL, false, // include legend
                 true, // tooltips?
                 false // URL generator? Not required...
                 );
         
         ((NumberAxis) ((XYPlot) chart.getPlot()).getRangeAxis())
                 .setAutoRangeIncludesZero(false);
         chart.setBackgroundPaint(Color.white);
//         chart.getLegend().setAnchor(org.jfree.chart.Legend.SOUTH);
         final XYPlot plot = chart.getXYPlot();
         // plot.setRenderer(new XYLineAndShapeRenderer());
         // ( (XYLineAndShapeRenderer)plot.getRenderer()).set.s
         plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
         // plot.getDomainAxis().setFixedAutoRange()
         plot.setDataset(1, exp);
         plot.mapDatasetToRangeAxis(1, 1);
       //  domainAxis.setRange(min[0], min[1]);
         final ValueAxis axis1 =Constants.logplot() ?  new LogarithmicAxis("Conditional probability distribution") : new NumberAxis("Conditional probability distribution");
         final ValueAxis axis2 =Constants.logplot() ? new LogarithmicAxis("Cumulative observed count") : new NumberAxis("Cumulative observed count");
         axis1.setTickLabelsVisible(true);
         axis2.setTickLabelsVisible(true);
         axis1.setTickLabelFont(font4);
         axis2.setTickLabelFont(font4);
         if(Constants.logplot()){
        	 axis1.setAutoTickUnitSelection(false);
        	 axis2.setAutoTickUnitSelection(false);
         }
       
        // axis1.setStandardTickUnits(new NumberTickUnit());
         plot.setRangeAxis(0,axis1);
         plot.setRangeAxis(1, axis2);
         final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
         final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
         
         int ik=0;
         for(int i=0; i<obs.getSeriesCount(); i++){
            // if(ik==color.length) ik=0;
          Color colors = ca.getColor(obs.getSeries(i).getKey().toString());
          renderer1.setSeriesPaint(i, colors);
          renderer2.setSeriesPaint(i, colors);
         // renderer.setSeriesShape(i font16.createGlyphVector(arg0, arg1));
          ik++;
         }
         // renderer2.setToolTipGenerator(new
         // StandardCategoryToolTipGenerator());
         if(dotted!=null)renderer2.setBaseStroke(dotted);
         renderer1.setShapesVisible(false);
         renderer2.setShapesVisible(false);
         plot.setRenderer(0, renderer1);
         plot.setRenderer(1, renderer2);
         plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
      //  LegendTitle legend =(LegendTitle) chart.getLegend();
        
     // legend.setItemFont(font8);
         if(sn_params!=null){
         for(Iterator<Map.Entry<Double[], String[]>> para = sn_params.entrySet().iterator(); para.hasNext();){
        	 Map.Entry<Double[], String[]> ent = para.next();
        	 String[] val = ent.getValue();
        	 Double[] x = ent.getKey();
        	 if(x!=null && Constants.logplot() && Constants.addAnnotationToDistGraphs()){
	     	XYTextAnnotation annot = new XYTextAnnotation(val[0]
					+ "", x[0], x[1]);
			annot.setFont(font10);
			annot.setPaint(ca.getColor(val[1]));
			plot.addAnnotation(annot);
        	 }
         }
         }
         return chart;
     }*/
    
    
   // XYSeriesCollection[] lu_rpair =new XYSeriesCollection[2];
   // XYSeriesCollection[] lu_r= new XYSeriesCollection[2];
  //  XYSeriesCollection[] lu_b = new XYSeriesCollection[2];
   /* public ChartPanel plotR(List<ProbabilityDistribution> ss1,
            ProbMultivariate mvf, double[] sum, int data_index, int index1, XYSeriesCollection[] lowerupper,
            String name, String namey, ColorAdapter ca) {
        List<Double[]> cis = new ArrayList<Double[]>();
        int len = ss1==null ? mvf.pairs.size() : ss1.size();
        for(int i=0; i<len; i++){
            cis.add(new Double[] {null, null});
        }
        double[] range = new double[] {0.05, 0.95};
        List<String> names = new ArrayList<String>();
        if(mvf!=null){
            for(int i=0; i<mvf.pairs.size(); i++){
                names.add(mvf.getCompoundName(i));
            }
        }
        else{
            for(int i=0; i<ss1.size(); i++){
                names.add(ss1.get(i).name());
            }
        }
        Map<Double[], String[]> det = new HashMap<Double[], String[]>();
        XYSeriesCollection[] obs = mvf==null ?  plot(ss1, sum, cis, range, names, det) : plot(mvf.pairs, sum, cis, range, names, det);
         XYSeriesCollection theor = obs[1];
       
            lowerupper[0] = new XYSeriesCollection();
            lowerupper[1] = new XYSeriesCollection();
          
            for(int i=0; i<theor.getSeriesCount(); i++){
                XYSeries lowers = new XYSeries(theor.getSeries(i).getKey());
                XYSeries uppers = new XYSeries(theor.getSeries(i).getKey());
               
                Double[] ci =cis.get(i);
                lowers.add(this.location.get(0), ci[0]); lowers.add(this.location.get(noSnps), ci[0]);
                uppers.add(this.location.get(0), ci[1]); uppers.add(this.location.get(noSnps),  ci[1]);
                lowerupper[0].addSeries(lowers);
                lowerupper[1].addSeries(uppers);
            }
       
       // JFreeChart chart = getChart(obs[0], obs[1],name,namey, ca, det);
     
            chartDist[data_index][index1].removeAll();
            final ChartPanel cp  = new ChartPanel(chart,
            		 1200, //width
            		 400, //height
            		 1200, //mindrawWidth
                     200, //mindrawHeight
                     1200, //maxDrawWith
                     400,//maxDrawHeight
                     ChartPanel.DEFAULT_BUFFER_USED,
                     true,  // properties
                     true,  // save
                     true,  // print
                     true,  // zoom
                     true   // tooltips		
            
            );

          //  chartDist[index].setMinimumSize(cp.getMinimumSize());
            chartDist[data_index][index1].add(cp);
            cp.setMinimumSize(dim);
            cp.setSize(dim);
            chartDist[data_index][index1].setMinimumSize(dim);
            chartDist[data_index][index1].setSize(dim);
         return cp;
//            Rectangle r = chartDist[index].getBounds();
          //  chartDist[index].repaint(1000);//(jpDist.getGraphics());
    }*/




    private double getYMax(XYSeriesCollection theor) {
        double max = Double.NEGATIVE_INFINITY;
        for(int i=0; i<theor.getSeriesCount(); i++){
            XYSeries series = theor.getSeries(i);
            for(int j=0; j<series.getItemCount(); j++){
                double y = series.getY(j).doubleValue();
                if(y>max){
                    max = y;
                }
            }
        }
        return max;
    }
    private double[] getCI(XYSeries ser, double pl, double pm) {
        double max = ser.getY(ser.getItemCount()-1).doubleValue();
        double[] res = new double[2];
        for(int i=0; i<ser.getItemCount(); i++){
            if(ser.getY(i+1).doubleValue() / max> pl ){
                res[0] = i;
                break;
            }
        }
        for(int i=ser.getItemCount()-1; i>=0;i--){
            System.err.println(ser.getY(i-1).doubleValue()/max);
            if(ser.getY(i-1).doubleValue() / max< pm ){
                res[1] = i;
                break;
            }
        }
        return res;
    }


    private static double sum(Collection<ProbabilityDistribution> s1) {
        double sum=0;
        for(Iterator<ProbabilityDistribution> it = s1.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
                sum+=pdist.sum();
        }
        return sum;
    }




    public static final List<ProbabilityDistribution> extract(Map<String, List<ProbabilityDistribution>> m){
        List<ProbabilityDistribution> res = new ArrayList<ProbabilityDistribution>();
        for(Iterator<List<ProbabilityDistribution>> it = m.values().iterator(); it.hasNext();){
            res.add(it.next().get(0));
        }
        return res;
    }
    public static final Map<String, List<ProbabilityDistribution>> transform(Collection<ProbabilityDistribution> s2){
        Map<String, List<ProbabilityDistribution>> m = new HashMap<String, List<ProbabilityDistribution>>();
        
        for(Iterator<ProbabilityDistribution> it = s2.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
            String st = pdist.toString();
            List<ProbabilityDistribution > l = m.get(st);
            if(l==null){
                m.put(pdist.toString(), l = new ArrayList<ProbabilityDistribution>());
            }
            l.add(pdist);
        }
        return m;
    }
    
   static boolean rescale = false;
public static double plotThresh = 1e-3;
    public static XYSeriesCollection[] plot(List<ProbabilityDistribution> coll, double[] sums, List<Double[]> cis,double[] mm, List<String> names,
         Map<Double[], String[]> det){
        XYSeriesCollection datas1 = new XYSeriesCollection();
        XYSeriesCollection datas2 = new XYSeriesCollection();
    //    double tot =0;
        for(int ij=0; ij<coll.size(); ij++){
        	ProbabilityDistribution pdist1 =  coll.get(ij);
        	SkewNormal pdist = pdist1 instanceof SkewNormal? (SkewNormal) pdist1 :  (SkewNormal) ((Mixture)pdist1).dist[0];
            if(pdist !=null ){
               double sum = 
                   sums==null ? 
                  pdist.sum() : sums[ij];
         //      tot+=sum;
               //   if(sum==0) sum=1;
                  XYSeries obs = new XYSeries(names.get(ij));
                  String key = names.get(ij).toString();
                 boolean swtch =   key.startsWith("1.0");
                  double cum = pdist.plotObservations(names.get(ij), true, obs, swtch);
                  
                   XYSeries theor = new XYSeries(names.get(ij));
                   pdist.plotTheoretical(names.get(ij), false, theor);
                  
                    double max = 0;
                    double theoSum=0;
                  for(int ik=0; ik<theor.getItemCount(); ik++){
                      theoSum+=theor.getY(ik).doubleValue();
                      if(theor.getY(ik).doubleValue()>max) max = theor.getY(ik).doubleValue();
                  }
                  double cumSum =0;
                  Double[] ci = cis.get(ij);
                  if(Math.abs(pdist.shape())<0.001){
                      pdist.ci(ci, mm);
                  }
                  else{
                  for(int ik=0; ik<theor.getItemCount(); ik++){
                      cumSum+=theor.getY(ik).doubleValue();
                      if(ci[0]==null && cumSum/theoSum>mm[0]){
                          ci[0] = theor.getX(ik).doubleValue();
                      }
                      if(ci[1]==null && cumSum/theoSum>=mm[1]){
                          ci[1] = theor.getX(ik).doubleValue();
                          break;
                      }
                  }
                  }
                  Double[] mode = new Double[] {pdist.location,0.5*sum};//( !Constants.logplot() ? 1.0 : sum)};//0.1*(sum/max) };
                  det.put(mode, new String[] {pdist.toString(), names.get(ij)});
                
               //   int obsCount =0;
                  System.err.println("key "+key);
                  /*if(swtch){
                	  for(int ik=obs.getItemCount()-1; ik>=0; ik--){
                		  double v = obs.getY(ik).doubleValue();
                		  double newV =cum - v; 
                		  if(newV > plotThresh){
                		  obs.update(ik, new Double( newV));
                		//  obsCount++;
                		  }
                		  else{
                			  obs.remove(ik);
                		  }
                	  }
                  }*/
                  for(int ik=theor.getItemCount()-1; ik>=0; ik--){
                	  double frac = (theor.getY(ik).doubleValue()/max);
                		/*  if(swtch){
                			  frac = 1.0 - frac;
                          }*/
                	  if(cum<plotThresh && ! Constants.logplot())   {
                      	theor.remove(ik);
                      }
                	  else if(rescale){
                        double the =  frac* sum;//( !Constants.logplot() ? 1.0 : sum);
                        if(the>plotThresh ||  Constants.logplot()){
                          theor.update(ik, new Double(the));
                        }
                	  }
                    
                    }
                /*  if(normalise){
                  for(int ik=0; ik<obs.getItemCount(); ik++){
                     obs.update(ik, new Double((obs.getY(ik).doubleValue()/sum)));
                  }
                  }*/
                 // if(theor.getItemCount()>0)
                    datas2.addSeries(theor);
                  //if(obs.getItemCount()>0)  
                  datas1.addSeries(obs);
            }
        }
        
        return new XYSeriesCollection [] {datas1, datas2};
//      
        /*if(true){
            try{
                Thread.currentThread().wait(10000);
            }catch(Exception exc){
                exc.printStackTrace();
            }
        }*/
    }
    
    long lastupdate = -1;
    
 public void changeAxis(Range range, boolean dox,  int k1){
		
		for(int kk=0; kk<chartRNew.length; kk++){
			if(chartRNew[kk]==null) continue;
			//if(kk!=k1){
				ValueAxis axisr, axisb;
				if(dox){
				
					axisr = chartRNew[kk].getChart().getXYPlot().getDomainAxis();
					axisb = chartBNew[kk].getChart().getXYPlot().getDomainAxis();
					
				}
				else{
					axisr = chartRNew[kk].getChart().getXYPlot().getRangeAxis();
					axisb = chartBNew[kk].getChart().getXYPlot().getRangeAxis();
				}
				if(!axisr.getRange().equals(range)){
					axisr.setRange(range, true, false);
				}
				if(!axisb.getRange().equals(range)){
					axisb.setRange(range, true, false);
				}
				chartRNew[kk].updateUI();
				chartBNew[kk].updateUI();
			//}
		}
    }
//	AffineTransform at = new AffineTransform();
 public void changeAxis1(Range range, boolean dox){
		for(int kkk=0; kkk<chartGNew.length; kkk++){
		for(int kk=0; kk<this.chartGNew[kkk].length; kk++){
			for(int kk1 = 0; kk1< this.chartGNew[kkk][kk].length; kk1++){
			if(chartGNew[kkk][kk][kk1]==null) continue;
			//if(kk!=k1){
				ValueAxis axisr;
				if(dox){
				
					axisr = chartGNew[kkk][kk][kk1].getChart().getXYPlot().getDomainAxis();
					
					
				}
				else{
					axisr = chartGNew[kkk][kk][kk1].getChart().getXYPlot().getRangeAxis();
					
				}
				if(!axisr.getRange().equals(range)){
					axisr.setRange(range, true, false);
				}
				
				chartGNew[kkk][kk][kk1].updateUI();
				//chartBNew[kk].updateUI();
			//}
		}
		}
		}
 }
 //final ColorAdapter[] ca_r, ca_b;
 final List<String> snpid;
  //int cnt1 =0;int cnt=0;
 final double[]annot_inter = Constants.annotationInterval();
 final double[] res_inR = new double[3];
 final double[] mean = new double[2];
 final DoubleMatrix2D matr = new DenseDoubleMatrix2D(2,2);
 final double[] res_inB = new double[3];
 
 //final double[] b_frac;
 //final boolean plotSig ;
public void updateClusterSeries(){
//	System.err.println("updating");
	for(int ii=0; ii<SignificancePlot.probesToPlot.length; ii++){
		for(int i=0; i<this.index.size(); i++){
			this.getRBSeries( i, ii);
		}
		
	}
	for(int ii=0; ii<rGlobal.length; ii++){
		for(int j=0; j<rGlobal[ii].length; j++){
			rGlobal[ii][j] = this.newRBSeries();
		}
	}
//	this.rGlobal
}
     
 
List<SignificancePlot> sp = new ArrayList<SignificancePlot>();

public void setSignficancePlot(SignificancePlot sp){
	this.sp.add(sp);
}
public void setSignficancePlot(List<SignificancePlot> sp){
	this.sp = sp;
}
LocPanel locP;

//JPanel sigPanel;
public synchronized void propertyChange(PropertyChangeEvent arg0) {
    	try{
        String nme = arg0.getPropertyName();
     //   Logger.global.info("prop change "+nme);
        if(nme.equals("pval") && SignificancePlot.snp_alias!=null ){
	        		
	        		
	      if(Constants.showScatter())      	updateClusterSeries();
	            if(plot==2)	this.locP = this.updateLocPanel(SignificancePlot.snp_alias);
	            /*	if(Constants.plot>=2 || plot==2){
	            		this.sigPanel.removeAll();
	            		ChartPanel currentGPanel = sp.getChartPanel(
lowIndex, true, jG.getSize().width)[0];
	            		currentGPanel.setMinimumSize(dim5);
		        		currentGPanel.setSize(dim5);
	            	   sigPanel.add(currentGPanel);
	            	}*/
        }
        if(nme.equals("setToPlot")){
        	if(Constants.printPlots()){
            int level = (Integer) arg0.getNewValue();
            setToPlot(level);
        	}
            return;
        }
        if(Constants.plot()<=1 && plot==0) return;
     
        if(nme.equals("init")){
            
            if(Constants.plot()>1){
                Object[] obj = (Object[]) arg0.getNewValue();
                Integer l = (Integer)obj[1];
                reinitialise(l);
              
            }
         
        }
         
        else if(nme.equals("pre_exp")){
          /*if(plot==2 && Constants.printPlots()){
        	setup(this.cnt);
          }*/
        }
        else if(nme.equals("emiss")){
            if(plot<2 && Constants.plot()<2) return;
            Object[] obj = (Object[]) arg0.getNewValue();
            Integer l = (Integer)obj[1];
            StateDistribution dist = (StateDistribution) obj[0];
            Integer i = (Integer)obj[2];
            HaplotypeEmissionState sta = (HaplotypeEmissionState) bwt.data[l];
            PseudoDistribution dist1 = sta.emissions[i];
            PseudoDistribution distr = dist1 instanceof MixtureDistribution ? ((MixtureDistribution)dist1).dist[0] : dist1;
          
      //     if(sta.noCop()!=ploidy) return;
           if(distr instanceof CompoundDistribution){
        	   List<PseudoDistribution> l1 =  ((CompoundDistribution)distr).l;
        	   for(int ii=0; ii<l1.size(); ii++){
        		   PseudoDistribution distr1 = l1.get(ii);
        		   if(distr1 instanceof PseudoMixture){
                	   addedInformation(dist, l, i, ((PseudoMixture)distr1).dist[0], 0, sta);
        		   }
        		   else addedInformation(dist, l, i, distr, ii, sta);
        	   }
           }else if(distr instanceof PseudoMixture){
        	   addedInformation(dist, l, i, ((PseudoMixture)distr).dist[0], 0, sta);
           }else{
        	   addedInformation(dist, l, i, distr, 0, sta);
           }
        }
        else if(nme.equals("finished")){
        	if( !Constants.r_panel() || !Constants.b_panel()) return;
            if(plot<2 && Constants.plot()<2 ) return;
         
            Object[] obj = (Object[]) arg0.getNewValue();
            final  Integer l = (Integer)obj[1];
             if(!include(l)) return ;
             int ind1 = indicesToInclude.indexOf(l);
             int ind = (int) Math.floor((double)ind1/(double) this.numPerPage);
             int cnt1 = ind1-this.numPerPage*ind;
             String st = bwt.data[l].name;
        int datind = bwt.data[l].dataIndex(0);
          String datname = (datind>=0) ?  this.names1[datind]+"_" : "";
             try{
            	final  ChartPanel cr = this.plotR(datname, st, l);
            	
            cr.getChart().getXYPlot().getRangeAxis().setAutoRange(true);
            cr.getChart().getXYPlot().getRangeAxis().setAutoRangeMinimumSize(3.0);
           final ChartPanel cb = this.plotB(datname, st, l,0);
           cr.setMinimumSize(dim);
            cb.setMinimumSize(dim);
            cr.setSize(dim);
            cb.setSize(dim);
       
            if(Constants.plot()>1){
                chartRNew[l] = cr;
                chartBNew[l] = cb;
                
               
            }
            if(plot==2){
                JComponent currentBPanel = new JPanel();
                JComponent currentRPanel = new JPanel();
                currentBPanel.setMinimumSize(dim);
                currentRPanel.setMinimumSize(dim);
               currentBPanel.setSize(dim);
                currentRPanel.setSize(dim);
                  currentRPanel.removeAll(); 
                 currentBPanel.removeAll();
                 currentRPanel.add(cr);
                 currentBPanel.add(cb);
               
               
                    try{
                    	if(Constants.printPlots())
                    	{
                    
                    		if(true){
                    			boolean joinPlots = false;
                    			boolean png = true;
                    			Dimension dim_ = joinPlots ? dimBoth: this.dim;
                    			   //ImageGraphics2D g  = 
                    			   VectorGraphics g = 
                    				   png?
                    						   new ImageGraphics2D(new File(chartBF, "_"+st+".png"),dim_, ImageConstants.PNG):
                    				   new PDFGraphics2D(new File(chartBF, "_"+st+".pdf"),dim_);
            	              //     ImageGraphics2D g1  = new ImageGraphics2D(new File(chartRF, "_"+st+".png"),currentRPanel, ImageConstants.PNG); 
            	                   g.setDeviceIndependent(Constants.plot()==1);
            	                   g.startExport();
            	                   g.setDeviceIndependent(Constants.plot()==1);
            	                 //  g1.startExport();
            	                   cr.print(g);
            	                   if(joinPlots){
            	                	   AffineTransform at1 = new AffineTransform();
            	                	   at1.setToTranslation(0,this.height);
            	                	   g.setTransform(at1);
            	                   }else{
            	                	   g.endExport();
            	                	   g = 
                        				   png?
                        						   new ImageGraphics2D(new File(chartBF, "B_"+st+".png"),dim_, ImageConstants.PNG):
                        				   new PDFGraphics2D(new File(chartBF, "_"+st+".pdf"),dim_);
                        						   g.setDeviceIndependent(Constants.plot()==1);
                            	                   g.startExport();
                            	                   g.setDeviceIndependent(Constants.plot()==1);
                            	                 //  g1.startExport();
                            	                   cb.print(g);
                        						   
            	                   }
            	                   g.endExport();
            	                //   g1.endExport();
                    		}
		                   
                    	/*	else{
	                    			AffineTransform at = new AffineTransform();
	                    		
			                   at.setToTranslation(0,this.height*cnt1);
			                   AffineTransform at1 = new AffineTransform();
			                   at1.setToTranslation(0,this.height*cnt1);
			               
			                  
			                   gB.setTransform(at);
			                   gR.setTransform(at1);
	             
	                   cb.print(gB);
	                   cr.print(gR);
                    		}*/
             
                  
                 
                    	}
                                             }catch(Exception exc){
                                                
                                                
                           exc.printStackTrace();
                       }
                 
                  
              
                
         
               
          }
             }catch(Exception exc){
                // System.err.println(this.index+" "+r_name);
                 exc.printStackTrace();
                 
             }
        }
        else if(nme.equals("expec_i")){
        
            
        }
        else if(false && nme.equals("hmm_maximisation")  ){
        	AbstractDistributionCollection dc = DistributionCollection.dc;
        	if(dc!=null){
		        for(int i2=0; i2<chartGNew[0].length; i2++){
		    		for(int ii=0; ii<chartGNew[0][i2].length; ii++){
		    			
		    			addAnnotation(chartGNew[0][i2][ii].getChart().getXYPlot(),i2,ii,false, stroke1, this.emStSp.getColor(true), dc);
		    		}
		    		
		    	}
        	}
        }
        else if(nme.equals(toPlotString) ){
        	double minx, maxx, miny, maxy;
        	Range xrange, yrange;
            if(plot<1 && Constants.plot()<2) return;
            int[] snp_alias = SignificancePlot.snp_alias;
            
            if(Constants.plot>=2 ||( plot ==2 && Constants.plot()==1)){
            	try{
        this.locP = this.updateLocPanel(SignificancePlot.snp_alias);
            	}catch(Exception exc){
            		System.err.println(exc.getMessage());
//            		exc.printStackTrace();
            	}
	           	 Properties p = new Properties();
	  			p.setProperty("PageSize", "A4");	
	  			int exp = Constants.scatterPlotsPerPage();
	  			int max_len = Math.min(exp, SignificancePlot.snp_alias.length);
	  			//Dimension d= new Dimension(Constants.scatterWidth()*this.snp_alias.length
	  				//	, Constants.scatterWidth()*Math.min(exp, this.snp_alias[0].length));
	  			Dimension d= 
	  				new Dimension(Constants.scatterWidth()*max_len, Constants.scatterWidth()* index.size()+locP.getSize().height);
	  					
	  			
	  			Dimension d2= new Dimension(Constants.scatterWidth()*max_len
	  					, locP.getSize().height);
	  			Dimension d3 = new Dimension(Constants.scatterWidth()*2
	  					, Constants.scatterWidth()* index.size());
	  			// jG.getSize(d);
	  			if(jG!=null){
	  				if(jG[0]!=null){
	  				jG[0].setMinimumSize(d);
	  				jG[0].setPreferredSize(d);
	  				}
	  				if(jG[1]!=null){
	  				jG[1].setMinimumSize(d);
	  				jG[1].setPreferredSize(d);}
	  			}
	  			
	  			
	  			//JFrame frame = new JFrame();
	  			//frame.setContentPane(jG);
	  			//frame.pack();
	  			//frame.setVisible(true);
	  			 int cnt2=0;
	  		    int start =SignificancePlot.start_;
	  		    int end = SignificancePlot.end_;
	  		AffineTransform at = new AffineTransform();
	  		AffineTransform at1 = new AffineTransform();
	  			VectorGraphics g_G=null;// new PNGGraphics2D(new File(this.chartDistF, "scatter.png"),d);
	  			LocPanel locP1 = null;
	  		if(plot==2 && Constants.showScatter()){ 
	  			Dimension d_ =flip(d); 
	  			g_G =	new ImageGraphics2D(new File(this.chartDistF, this.location.get(SignificancePlot.snp_alias[start])+".png"),
	  					d_,ImageConstants.PNG);
	  			g_G.setProperties(p);
	  			g_G.startExport();
	  			g_G.setColor(Color.white);
	  			g_G.fillRect(0, 0, d_.getWidth(), d_.getHeight());
	  			g_G.setBackground(Color.white);
	  		  locP1 = this.updateLocPanel(snp_alias, start,max_len );
  	  		  locP1.setSize(d2);
  	  		  if(Constants.flip){
  	  			  at.setToIdentity();
  	  			 at.setToRotation(Math.PI/2,d2.getHeight()/2,d2.getHeight()/2);//,d2.getWidth()/2,d2.getHeight()/2);
  	  			at1.setToTranslation(Constants.scatterWidth()* index.size(), 0);
  	  			
  	  			 at1.concatenate(at);
		  		//  at.setToRotation(-Math.PI/2);// .setToTranslation(0,0);// Constants.scatterWidth()*index.size());
  	  			 //at.setToRotation(Math.PI/2,d2.getWidth()/2,d2.getHeight()/2);
	  			// at.setToTranslation(-d2.height, -d2.width);
  	  		  }
  	  	 g_G.setTransform(at1);
  	  		  locP1.print(g_G);
//	  			g_G.drawRect(0, 0, d.width, d.height);
	  			
	  		
	  		}
	  		
	  		ChartPanel[] cr_ = new ChartPanel[index.size()];
	  		ChartPanel[] cr1_ = new ChartPanel[index.size()];
	  		int[] prev_index = new int[index.size()];
	  		Arrays.fill(prev_index,-1);
	  		 for(int ii=0; ii<2; ii++){
	  			Arrays.fill(cr_, null);
		  		
	  		 for(int i2=0; i2<index.size(); i2++){
	  			if(!Constants.showScatter()) continue;
         	 XYSeriesCollection current_r =  this.rGlobal[i2] [ii];//: this.rGlobalProbeOnly[i2];//  this.rb[i2].get(ii);
         	 String name =  this.names1[index.get(i2)]+"_" +(ii==0 ? "global" : "global probe only");
         	 Double fracR = DistributionCollection.getFracGlobS(i2, ii>0, true);
         	 Double fracB = DistributionCollection.getFracGlobS(i2, ii>0, false);
         	 String[] form = DistributionCollection.getFormGlobS(i2);
        	// String[] formB = DistributionCollection.getFormGlobS(i2, false);
         	 if(current_r!=null) {
         		final ChartPanel cr =  new ChartPanel(graph(current_r, 
         			name, emStSp.getColor(false),i2,
         			fracR, fracB),
         			Constants.scatterWidth(),Constants.scatterWidth(),Constants.scatterWidth(),Constants.scatterWidth(),Constants.scatterWidth(),
         			Constants.scatterWidth(),
         	          false,  true,      true,  true,    true,   true );  // tooltips	
         	  //  JFrame fr = new JFrame();
         	//	fr.getContentPane().add(cr);
         	//	fr.pack();
         	//	fr.setVisible(true);
         		cr_[i2] = getSize(current_r)==0 ? null : cr;
      //   		atG.setToTranslation(Constants.scatterWidth*ii, Constants.scatterWidth()*i2);
  	  		
  	  		 
  	  		  //locP1.print(g_G);
  	  		
         		XYPlot xyp = cr.getChart().getXYPlot();
         		AbstractDistributionCollection dc = DistributionCollection.dc;
         		if(dc!=null){
         		addAnnotation(xyp,i2,ii, true,stroke, emStSp.getColor(true), dc);
         		}
         		if((plot==2 || Constants.plot>=2 )&& Constants.globalRange()){
            		double[] range = Constants.range();
            		minx = 0;
            		maxx = 0;
            		miny = range[2];
            		maxy = range[3];
            	
            		
	    				if(cr!=null){
	    					Range rx = cr.getChart().getXYPlot().getDomainAxis().getRange();
	    					Range ry = cr.getChart().getXYPlot().getRangeAxis().getRange();
	    					if(rx.getLowerBound() < minx) minx = rx.getLowerBound()-1;
	    					if(rx.getUpperBound() > maxx) maxx = rx.getUpperBound()+1;
	    					if(ry.getLowerBound() < miny) miny = ry.getLowerBound();
	    					if(ry.getUpperBound() > maxy) maxy = ry.getUpperBound();
	    				}
            		
            		xrange = new Range(minx, maxx);
            		yrange = new Range(miny - (0.1*(maxy-miny)), maxy);
            		if(form!=null){
            		xyp.getDomainAxis().setLabel(form[0]+"\t"+form[1]);
            	//	xyp.getDomainAxis(1).setLabel(form[1]);
            		if(form.length>2){
            			xyp.getRangeAxis().setLabel(form[2]+"\t"+form[3]);	
            			//xyp.getRangeAxis(1).setLabel(form[3]);	
            		}
            		}
	    			//xyp.getDom
	    				if(cr!=null){
	    					cr.getChart().getXYPlot().getDomainAxis().setAutoRange(false);
	    					cr.getChart().getXYPlot().getRangeAxis().setAutoRange(false);
	    					cr.getChart().getXYPlot().getDomainAxis().setRange(xrange);
	    					cr.getChart().getXYPlot().getRangeAxis().setRange(yrange);
	    				}
            	
            		
    		
	    			
    			
            	}
         		
         		 if(plot==2){
      	  			 cr.setSize(d3);
      	  			 cr.setMinimumSize(d3);
      	  			VectorGraphics g_Glob =		new ImageGraphics2D(new File(this.chartDistF, name+".png"),cr,ImageConstants.PNG);
    	  			g_Glob.setProperties(p);
    	  			g_Glob.startExport();
      	  			 cr.print(g_Glob);
      	  			if(g_Glob!=null) g_Glob.endExport();
      	  		 }
         		cr.setSize(dim5);
         	if(Constants.plot>=2)	this.chartGNew[1][i2][ii] =  cr;
         	
         	     }
         	 
         		}
	  		 }
	  		 
	  		
	  		
//	  		ChartPanel[] cr_ = new ChartPanel[index.size()];
	  	for(int ii=start; ii<=end; ii++){ 
	  		if(g_G!=null){
	  		if(((ii-start)-cnt2*exp)>=exp){
	  		 
	  			g_G.endExport();
	  			cnt2++;
	  			Dimension d_ = flip(d);
	  			g_G =// new PNGGraphics2D(new File(this.chartDistF, "scatter.png"),d);
	  				new ImageGraphics2D(new File(this.chartDistF, this.location.get(snp_alias[ii])+"_"+this.ploidy+"_"+".png"),d_,ImageConstants.PNG);
	  			g_G.setProperties(p);
	  			g_G.startExport();
	  			 locP1 = this.updateLocPanel(snp_alias, cnt2*exp+start,exp);
		  		  locP1.setSize(d2);
		  		g_G.setColor(Color.white);
	  			g_G.fillRect(0, 0, d_.getWidth(), d_.getHeight());
		  		 // at.setToTranslation(0, Constants.scatterWidth()*index.size());
		  		at.setToIdentity();
		  		 if(Constants.flip){
		  			 //2
		  			  at.setToIdentity();
		  	  			 at.setToRotation(Math.PI/2,d2.getHeight()/2,d2.getHeight()/2);//,d2.getWidth()/2,d2.getHeight()/2);
		  	  			at1.setToTranslation(Constants.scatterWidth()* index.size(), 0);
		  	  			
		  	  			 at1.concatenate(at);
		  			 
		  			//  at.setToIdentity();
		  	  		//	 at.setToRotation(Math.PI/2,d2.getHeight()/2,d2.getHeight()/2);//,d2.getWidth()/2,d2.getHeight()/2);
		  	  		//	at1.setToTranslation(Constants.scatterWidth()* index.size(), 0);
		  	  		//	at.concatenate(at1);
		  			// at.setToRotation(Math.PI/2);
		  			// at.setToTranslation(-d2.height, -d2.width);
		  			 
		  		 }
		  		
		  		  g_G.setTransform(at);
		  		 
		  		  //locP1.print(g_G);
	  		
	  			
	  		}
	  		}
	  		//at.setToTranslation(-cnt2*exp, 0);
	  		//g_G.setTransform(at);
	  		//this.sigPanel.print(g_G);
	  		
	  		
	  		
	  		
            for(int i2=0; i2<index.size(); i2++){
            	if(cr_[i2]!=null){
            		cr1_[i2] = cr_[i2];
            		cr_[i2]=null;
            		prev_index[i2] = ii-1;
            	}
            	
            	if(snp_alias[ii]>=0 &&  Constants.allowLocalDist()){
            			
            			
            		int i3 = Constants.plotMerged() ? 0 :i2;	
            	 XYSeriesCollection current_r = this.rb[i2].get(ii);
            	 double pv = 	sp.size()==0 ? 1.0 : this.sp.get(i3).getRSeriesCollection(type_index).getSeries(this.phenoIndex).getY(ii).doubleValue();
            	 String name =  this.names1[index.get(i2)]+"_"+
     			SignificancePlot.probesToPlot[ii]+"_"+this.location.get(snp_alias[ii])+"_"+
     		(	Constants.collapseScatterInd()  || !Constants.calcAssoc? "" :
     			String.format("pval= %5.3g",
     				//	-Math.log10(
     						pv));
     				//	SignificancePlot.pvalueToP[this.index.get(i2)][ii]));
            
            //	 DataCollection.datC.getBGCount(i2, snp_alias[i2][ii]);
            	 Double fracR = DistributionCollection.getFracS(i2,snp_alias[ii], true);
            	 Double fracB = DistributionCollection.getFracS(i2,snp_alias[ii], false);
            	 String[] form = useVals ? null : DistributionCollection.getFormS(i2,snp_alias[ii]);
            	// String formB = DistributionCollection.getFormS(i2,snp_alias[ii], true);
            //	 System.err.println("frac "+frac);
            //	 if(current_r.getI)
            	 if(current_r!=null) {
            		final ChartPanel cr =  new ChartPanel(graph(current_r, 
            			name, emStSp.getColor(false),i2,
            			fracR, fracB),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            	          false,
            	            true,  // properties
            	            true,  // save
            	            true,  // print
            	            true,  // zoom
            	            true   // tooltips		
            	);
            		//cr.translateJava2DToScreen(current_r.getSeries(0).getDataItem(0));
            		//cr
            		//int i3 = Constants.collapseScatterInd() ? 0  : i2;
            		cr_[i2] = getSize(current_r)==0 ? null : cr;
            		//double minPv = this.sp.get(index.get(i2)).minPv;
            		if(pv>1.301){
            		//	cr.getChart().setBorderVisible(true);
            		//	cr.getChart().setBackgroundPaint(Color.yellow);
            		//	cr.getChart().setBorderPaint(Color.LIGHT_GRAY);
            		//	cr.getChart().getTitle().setFont(font20);
            		}
            	//	current_r.g
            		XYPlot xyp = cr.getChart().getXYPlot();
            		if(!this.useVals && Constants.annotateClusterPosition()){
            			AbstractDistributionCollection dc = DistributionCollection.dc;
            			if(dc!=null){
            				addAnnotation(xyp,i2,ii,false, stroke, emStSp.getColor(true), dc);
            			}
            		}
            		if(Constants.annotateName()){
                		addAnnotation1(xyp,i2,ii);
                		}
            		
            		cr.setSize(dim5);
            		if(form!=null){
            		xyp.getDomainAxis().setLabel(form[0]+" // "+form[1]);
            		//xyp.getDomainAxis(1).setLabel(form[1]);
            		if(form.length>2){
            			xyp.getRangeAxis().setLabel(form[2]+" // "+form[3]);
            		//	xyp.getRangeAxis(1).setLabel(form[3]);
            		}
            		}
            	if(Constants.plot>=2)	{
            		this.chartGNew[0][i2][ii] =  cr;
            		
            	}
            	
            		// if(Constants.plot==2){
            			// addList(cr, false); addList(cr, true);
                 //		}
            	     }
            	 
            		}
            		
            	}
            	if(plot==2 || Constants.plot>=2){
            		double[] range = Constants.range();
            		minx = range[0];
            		maxx = range[1];
            		miny = range[2];
            		maxy = range[3];
            	
            		for(int jk=0; jk<cr_.length; jk++){
            			ChartPanel cr = cr_[jk];
            			
	    				if(cr!=null){
	    					Range rx = cr.getChart().getXYPlot().getDomainAxis().getRange();
	    					Range ry = cr.getChart().getXYPlot().getRangeAxis().getRange();
	    					if(rx.getLowerBound() < minx) minx = rx.getLowerBound()-1;
	    					if(rx.getUpperBound() > maxx) maxx = rx.getUpperBound()+1;
	    					if(ry.getLowerBound() < miny) miny = ry.getLowerBound();
	    					if(ry.getUpperBound() > maxy) maxy = ry.getUpperBound();
	    				}
            		}
            		xrange = new Range(minx, maxx);
            		yrange = new Range(miny, maxy);
            		
	    			for(int jk=0; jk<cr_.length; jk++){
	    				ChartPanel cr = cr_[jk];
	    				if(cr!=null){
	    					if(Constants.globalRange()){
	    					cr.getChart().getXYPlot().getDomainAxis().setAutoRange(false);
	    					cr.getChart().getXYPlot().getRangeAxis().setAutoRange(false);
	    					cr.getChart().getXYPlot().getDomainAxis().setRange(xrange);
	    					cr.getChart().getXYPlot().getRangeAxis().setRange(yrange);
	    					}
	    					if(plot==2){
	    						int jk1 = jk;
	    						if(Constants.collapseScatterInd() && jk>0 && ((MergedDataCollection)DataCollection.datC).map[snp_alias[ii]][jk]!=null){
	    	            			
	    	            				jk1 = 0;
	    	            			
	    	            		}
	    					if(Constants.flip) {
	    						at.setToTranslation(jk1*Constants.scatterWidth(), ((ii-start) - exp*cnt2)*Constants.scatterWidth());
	    					}
	    					else{
	    						at.setToTranslation( ((ii-start) - exp*cnt2)*Constants.scatterWidth(),
	    							
	    							d2.getHeight()+jk1*Constants.scatterWidth());
	    					}
	    	    			g_G.setTransform(at);
	    					cr.print(g_G);
	    					if(Constants.drawAnnotationLines()&& prev_index[jk]>=0 && !this.useVals){
	                			//ChartPanel cr_1 = chartGNew[0][jk][ii-1];
	                			addAnnotation1(cr,cr1_[jk], jk, ii, prev_index[jk], g_G);
	                		}
	    					
	    					}
	    				}
            		}
            		
    		
	    			
    			
            	}
            	
              
            }
           
           
            if(Constants.plot()>1 ||( plot==2 && Constants.plot()==1)){
                update();
                this.updateUI();
            }
            if(plot==2 && g_G!=null){
                try{
                     /*   if(gB_!=null){
                      for(int i=0; i<gB_.length; i++){
                    	  if(gB_[i]!=null){
                    	  gB_[i].endExport();
                    	  gR_[i].endExport();
                    	  }
                      }
                        }*/
                  //  	jG.setMinimumSize(d);
        	  	//		jG.setPreferredSize(d);
                 //   	jG.print(g_G);
                      
             			g_G.endExport();
             			//chartGNew[0][0].print(g_G);
             	//		g_G1.setTransform(at1);
             		//	g_G1.setDeviceIndependent(Constants.plot()==1);
             		//	cr.print(g_G1);
             		//	g_G1.endExport();
                    	// this.g_G.endExport();
                }catch(Exception exc){
                    exc.printStackTrace();
                }
                this.plot = 0;
                
            }
            }
            
        }
        else if(nme.equals("dist_maximisation")){
           
        }
        else if(nme.equals("hmm_maximisation")){
            
        }
    //    else throw new RuntimeException("!! "+nme);
     //   Logger.global.info("exiting "+nme);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    	
    }
   private Dimension flip(Dimension d) {
if(Constants.flip){
	return new Dimension(d.height, d.width);
}
else{
	return d;
}
}

private int getSize(XYSeriesCollection current_r) {
	int cnt=0;
	for(int k=0; k<current_r.getSeriesCount(); k++){
		cnt+=current_r.getItemCount(k);
	}
	return cnt;
}

private void addList(final ChartPanel cr, boolean domain) {
	   Axis ax = domain ? cr.getChart().getXYPlot().getDomainAxis() :cr.getChart().getXYPlot().getRangeAxis()
		   ;
		ax.addChangeListener(new AxisChangeListener(){
				

				public void axisChanged(AxisChangeEvent event) {
					
					ValueAxis axis = ((ValueAxis)event.getAxis());
					boolean dox = axis==cr.getChart().getXYPlot().getDomainAxis();
					Range range = axis.getRange();
					if(dox){
					//	System.err.println("fired property change in "+getName());
						IndividualPlot.this.firePropertyChange("axis", null, range);
					}
					changeAxis1(range, dox);
				}
         	
         });
	
}

ChartPanel plotR(String datname, String st, final  Integer l) {
	  String r_name =  
        	datname+
        	 st
        	// "LRR on Agilent 185k array colored by cnvHap CN prediction";
     	 //" LRR on Illumina 317k chip colored by cnvHap CN prediction";
        	 //+"LRR on Illumina 1M chip (dots) and on Agilent 244k aCGH array (bars) colored by CN genotype prediction"
        	 ;
        	// "LRR on 244k aCGH array colored by cnvHap CN prediction";
        	//"LRR on Illumina Human1M chip colored by cnvHap CN prediction";
        	// "Predicted CN using ADM2 on 244k aCGH array (bars) mapped to 1M chip (dots)";
        	 //"Predicted CN using PennCNV on Illumina 1M chip (dots) mapped to 244k aCGH array (bars)";//"Copy number "+st;
  XYSeriesCollection[] current_r = Constants.plot()>1 ?  this.getRSeriesCollection(l) : rdc.remove(l);
 
  if(current_r==null) return null;
 
  final ChartPanel cp =   new ChartPanel(graph(current_r , r_name, emStSp.getColor(false),shapes, false),
  		 Constants.r_panel_width(), //width
  		 Constants.r_panel_height(), //height
  		 Constants.r_panel_width(), //width
  		 Constants.r_panel_height(), //height
  		 Constants.r_panel_width(), //width
  		 Constants.r_panel_height(), //height
           ChartPanel.DEFAULT_BUFFER_USED,
           true,  // properties
           true,  // save
           true,  // print
           true,  // zoom
           true   // tooltips		
  
  );
  if(Constants.logplot()){
	 
 // final ValueAxis axis1 =Constants.logplot() ?  new LogarithmicAxis("Conditional probability distribution") : new NumberAxis("Conditional probability distribution");
  final ValueAxis axis2 =Constants.logplot() ? new LogarithmicAxis("Log Read Depth") : new NumberAxis("Log Read Depth");
  //axis1.setTickLabelsVisible(true);
  axis2.setTickLabelsVisible(true);
  //axis1.setTickLabelFont(font4);
  axis2.setTickLabelFont(font4);
  axis2.setLabelFont(font4);
 // if(Constants.logplot()){
 	// axis1.setAutoTickUnitSelection(false);
 	 axis2.setAutoTickUnitSelection(false);
  
  cp.getChart().getXYPlot().setRangeAxis(axis2);
  }
  cp.getChart().getXYPlot().getRangeAxis().setTickLabelInsets(new RectangleInsets(5,5,5,5));
	//cp.getChart().getXYPlot().getDomainAxis()
  if(Constants.plot>1){
  	cp.getChart().getXYPlot().getDomainAxis().addChangeListener(new AxisChangeListener(){
		

		public void axisChanged(AxisChangeEvent event) {
			
			ValueAxis axis = ((ValueAxis)event.getAxis());
			boolean dox = axis==cp.getChart().getXYPlot().getDomainAxis();
			Range range = axis.getRange();
			if(dox){
			//	System.err.println("fired property change in "+getName());
				IndividualPlot.this.firePropertyChange("axis", null, range);
			}
			changeAxis(range, dox,  l);
		}
 	
 });

  }
  return cp;
	}

public ChartPanel plotB(String datname, String st, final int l, int ik) {
	 String b_name = datname+"_"+st
  //+	 "BAF on Illumina 317k chip colored by cnvHap CN genotype prediction "
	 ;
       XYSeriesCollection[] current_b =  Constants.plot()>1 ? this.getBSeriesCollection(l) : bdc.remove(l);
   final  ChartPanel cp =   new ChartPanel(graph(current_b, st, emStSp.getColor(false),
    		shapes, true),
    		 Constants.r_panel_width, //width
    		 Constants.plotHeight, //height
    		 Constants.r_panel_width, //width
    		 Constants.plotHeight, //height
    		 Constants.r_panel_width, //width
    		 Constants.plotHeight, //height200
             ChartPanel.DEFAULT_BUFFER_USED,
             true,  // properties
             true,  // save
             true,  // print
             true,  // zoom
             true   // tooltips		
    );
    ValueAxis va = cp.getChart().getXYPlot().getRangeAxis();
    if(Constants.cumulativeR(0)<=1){
    va.setAutoRange(false);
    double minB = Constants.minB(ik);
    double maxB = Constants.maxB(ik);
    double extra = 0.0 *(maxB - minB);
    minB = minB - extra;
    maxB = maxB +extra;
    va.setRange(new Range(minB, maxB));
    }
    JFreeChart chart = cp.getChart();
    chart.getXYPlot().getRangeAxis().setTickLabelInsets(new RectangleInsets(5,5,5,5));
/*    final NumberFormat nf =new DecimalFormat("0.00");
    //nf.setMinimumFractionDigits(2);
    //nf.setMaximumFractionDigits(2);
    //nf.setMaximumIntegerDigits(5);
    //nf.setMinimumFractionDigits(5);
   //nf.setRoundingMode(RoundingMode.)
  // chart.getXYPlot().getRangeAxis().setTickMarkInsideLength(2);
   TickUnitSource tus =  NumberAxis.createStandardTickUnits();
   
    chart.getXYPlot().getRangeAxis().setStandardTickUnits(tus);
    /*new StandardTickUnitSource(){
  	   public TickUnit getCeilingTickUnit(double size) {
  		 NumberTickUnit un = (NumberTickUnit)super.getCeilingTickUnit(size);
  		   return   new NumberTickUnit(un.getSize(), 
  	               nf);
  	   }
  	   public TickUnit getLargerTickUnit(TickUnit unit) {
  		   NumberTickUnit un = (NumberTickUnit)super.getLargerTickUnit(unit);
  	        return new NumberTickUnit(un.getSize(), 
  	                nf);
  	    }
   // });
    chart.getXYPlot().getRangeAxis().setAutoTickUnitSelection(false);*/
    if(Constants.hideAxis() && !singlechrom){
    	
    	   XYPlot plot= cp.getChart().getXYPlot();
    	    ValueAxis domaxis = plot.getDomainAxis(0);
    	    ValueAxis raxis = plot.getRangeAxis(0);
    	   Range r=  raxis.getRange();
    	  // double pos2 = r.getUpperBound()-0.2*(r.getUpperBound()-r.getLowerBound();
    	   double pos1 = r.getUpperBound()+0.05;
    	   raxis.setRange(r.getLowerBound(), r.getUpperBound()+0.2*(r.getUpperBound()-r.getLowerBound()));

       	Map<String,Double> karyo = DataCollection.datC.readKaryotypes();
    	   double u = r.getUpperBound();
    	for(Iterator<String> it = karyo.keySet().iterator(); it.hasNext();){
    		String nxt = it.next();
    		
    		XYTextAnnotation annot = new XYTextAnnotation(nxt, karyo.get(nxt),pos1);
    		if(!this.singlechrom) annot.setFont(this.font101);
    		//annot.setPaint(Color.GRAY);
    		plot.addAnnotation(annot);
   //plot.getDomainAxis(0).setLabelPaint(Color.WHITE);
   
   domaxis.setTickLabelPaint(Color.white);
    	}
    }
  // plot.getDomainAxis().
  //  plot.setDomainAxis(1, xAxis2);
  //  plot.mapDatasetToDomainAxis(2, 1);
//    cp.getChart().getXYPlot()plot.setDomainAxis(1, xAxis2 );
    
    if(Constants.plot>1){
      	cp.getChart().getXYPlot().getDomainAxis().addChangeListener(new AxisChangeListener(){
    		

    		public void axisChanged(AxisChangeEvent event) {
    			
    			ValueAxis axis = ((ValueAxis)event.getAxis());
    			boolean dox = axis==cp.getChart().getXYPlot().getDomainAxis();
    			Range range = axis.getRange();
    			if(dox){
    			//	System.err.println("fired property change in "+getName());
    				IndividualPlot.this.firePropertyChange("axis", null, range);
    			}
    			changeAxis(range, dox,  l);
    		}
     	
     });

      }
    return cp;
	}


Algebra lg = new Algebra();
void addAnnotation1(XYPlot xyp, int i2, int ii) {
	AnnotationSeries m = this.annotation[i2].get(ii);
	m.addAnnotation(xyp);
}

public double[]  translateJava2DToScreen(ChartPanel cr,XYTextAnnotation annot) {
    Rectangle2D area = cr.getChartRenderingInfo().getPlotInfo().getDataArea();
	ValueAxis ra = cr.getChart().getXYPlot().getRangeAxis();
	ValueAxis da = cr.getChart().getXYPlot().getDomainAxis();
	double widthx = area.getWidth();
	double widthy = area.getHeight();
	double minx = area.getMinX();
	double miny = area.getMinY();
   double x =  (((annot.getX() - da.getLowerBound())/da.getRange().getLength()) * widthx  + minx);
   double y =    (((ra.getUpperBound() - annot.getY())/ra.getRange().getLength()) * widthy  + miny);
   return new  double[] {x,y};
}


final Map<String, Stroke>[] strokes;
void addAnnotation1(ChartPanel cr, ChartPanel cr1,  int i2, int ii,int ii_prev, VectorGraphics g) {
	AnnotationSeries m = this.annotation[i2].get(ii);
	AnnotationSeries m1 = this.annotation[i2].get(ii_prev);
	int numb = ii- ii_prev;
	for(Iterator<String> it = m.l.keySet().iterator(); it.hasNext();){
		String key = it.next();
		XYTextAnnotation annot1 = m1.l.get(key);
		if(annot1!=null ){
			
			XYTextAnnotation annot = m.l.get(key);
			if(annot1.getPaint().equals(annot.getPaint())){
			double[] point = translateJava2DToScreen(cr, annot);
		
			double[] point1 = translateJava2DToScreen(cr1, annot1);
			g.setPaint(annot.getPaint());
			Map<String, Stroke> strokes_ = strokes[i2];
			Stroke stroke_ = strokes_.get(annot.getText());
			if(stroke_==null){
				//int order = 2;//strokes_.size();
				double stroke_width = 
					
					Math.pow(Constants.joinStrokeWidth, Constants.decayStrokeWidth);//Math.max(0.5, 3.0*Math.pow(order, -0.8));  //power of -0.5 for slower decay
				strokes_.put(annot.getText(),
						stroke_ =  new BasicStroke((float)stroke_width, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, 
								//new float[] {5,5,5},1.0f));
								new float[] {3,9},//21.0f, 9.0f, 3.0f, 9.0f },
								//Constants.randomFloat(5,10f), 
								1.0f));
			    
						
			}
			g.setStroke(stroke_);
			if(Constants.flip){
				g.drawLine(point[0], point[1], point1[0] , point1[1]- numb*cr.getWidth());
			}else{
				g.drawLine(point[0], point[1], point1[0] - numb*cr.getWidth(), point1[1]);
			}
			}
//			l.add( new ColoredLine(point, point1, annot.getPaint()));
			
		}
	}
	g.setStroke(stroke);
	g.setColor(Color.black);
	//m.addAnnotation(xyp);
}



XYAnnotation[] annot = new XYAnnotation[2];

void addAnnotation(XYPlot xyp, int i2, int ii, boolean global, Stroke stroke, Color[] col, AbstractDistributionCollection dc) {
	   int bg = DataCollection.datC.getBGCount(i2, ii);
//	   Color[] col =  DistributionCollection.dc.getColR(bg);
	 for(int j=0; j<this.rbSeriesNames.size(); j++){
			int nocop = emStSp.getCN(j);
			XYSeriesCollection coll = rb[i2].get(ii);
			if(coll==null) continue;
			XYSeriesMiss series_ =(XYSeriesMiss) coll.getSeries(j);
		
			
		if(!Constants.showAll(0) && series_.getItemCount() < 1 && series_.nan_x==0) continue;
			
		//	double bfrac = this.b_frac[b_alias[j]];
			
			if(Constants.joint){
			
				ProbabilityDistribution2 dist_;
				if(global){
					dist_ = dc.getDistributionRBGlob((short)i2, nocop,emStSp.getBCount(j)
							);
				}
				else{
				  dist_ = dc.getDistributionRB((short)i2, nocop,emStSp.getBCount(j), 
						SignificancePlot.snp_alias[ii]);
				}
				if(dist_==null) continue;
				ProbabilityDistribution2 dist2 = dist_ instanceof Mixture2 ? ((Mixture2)dist_).dist[0] : dist_;
				 getAnnotation1(dist2, annot, col[j], stroke);
				for(int i=0; i<annot.length; i++){
					if(annot[i]!=null) xyp.addAnnotation(annot[i]);
				}
				if(Constants.showNaN() && series_.nan_x>0){
					double c = Math.log10(20*((double)series_.nan_x/DataCollection.datC.indiv().size()));
					
					double x = ((OrthogonalProbabilityDistribution) dist2).distx.getMean();
					double y = ((OrthogonalProbabilityDistribution) dist2).disty.getMean();
					Shape sh = new Ellipse2D.Double(x-c/2,y-c/2,c,c);
					  
					XYShapeAnnotation annot = new XYShapeAnnotation(sh,stroke,col[j]);
					xyp.addAnnotation(annot);
				}
				
				
			}
			else{
			
				ProbabilityDistribution dist, distB;
			dist= dc.getDistribution((short)i2, bg, nocop, SignificancePlot.snp_alias[ii]);
			 distB = dc.getDistribution((short)i2, j, SignificancePlot.snp_alias[ii]);
			
			if(dist==null) continue;
		
		
			dist.getInterval(annot_inter, res_inR);
		
				distB.getInterval(annot_inter, res_inB);
				
			
			/*	if(res_inR[2]<res_inR[1]){
					throw new RuntimeException("!!");
				}*/
		XYAnnotation annotR = new XYLineAnnotation(res_inR[0], res_inB[1], res_inR[2], res_inB[1], stroke,col[j]);
			XYAnnotation annotB = new XYLineAnnotation(res_inR[1], res_inB[0], res_inR[1], res_inB[2],stroke,col[j]);
			xyp.addAnnotation(annotR);
			xyp.addAnnotation(annotB);
			}
		}
		
	}

private void getAnnotation1(ProbabilityDistribution2 dist2,
		XYAnnotation[] annot2, Color col, Stroke stroke) {
	if(dist2 instanceof OrthogonalProbabilityDistribution){
		((OrthogonalProbabilityDistribution) dist2).getIntervalX(annot_inter, res_inR);
		((OrthogonalProbabilityDistribution) dist2).getIntervalY(annot_inter, res_inB);
		
		annot2[0] = new XYLineAnnotation(
				res_inR[0], res_inB[1] , res_inR[2],res_inB[1],
				stroke,col);
		annot2[1] = new XYLineAnnotation(
				res_inR[1], res_inB[0] , res_inR[1],res_inB[2],
				stroke,col);
	}
	else if(false){
		if(true) throw new RuntimeException("!!");
		dist2.getInterval(annot_inter, this.matr, mean);
		EigenvalueDecomposition evd = new EigenvalueDecomposition(matr);
		DoubleMatrix2D V = evd.getV();
		
		DoubleMatrix1D eig = evd.getRealEigenvalues();
		for(int i=0; i<V.columns(); i++){
			double len  = Math.sqrt(eig.get(i));
		    double x = V.get(0, i)*len;
		    double y = V.get(1, i)*len;
		    double x1 = mean[0] - x;
		    double y1 = mean[1] -y;
		    double x2 = mean[0]+x;
		    double y2 = mean[1]+y;
		   // double theta = Math.atan2(y, x);
		 //   if(Math.abs(theta) > Math.PI/4.0  && Math.abs(theta - Math.PI/2.0)>Math.PI/4.0){
		//    }
			XYAnnotation annotR = new XYLineAnnotation(
					x1, y1 , x2,y2,
				//    x2, mean[1] + V.get(1, i)) , 
					//mean[0] - V.get(i, 0)*Math.sqrt(eig.get(i)), mean[1] - V.get(i, 1)*Math.sqrt(eig.get(i)) , 
					//mean[0] + V.get(i, 0)*Math.sqrt(eig.get(i)), mean[1] + V.get(i, 1)*Math.sqrt(eig.get(i)) , 
					
					stroke,col);
			annot2[i] = annotR;
		}
	}
	
}

/*  private ChartPanel[][] getDistCharts() {
	ChartPanel[][] res = new ChartPanel[index.size()][];
	for(int i=0; i<res.length; i++){
		res[i] = getDistCharts(i);
	}
	return res;
}*/


/*private AbstractVectorGraphicsIO getVectorGraphics(int ind) throws Exception {
	  
	   AbstractVectorGraphicsIO res;
	   int numToInclude = indicesToInclude.size();
	   this.setSize(1200,200*(Math.min(numPerPage,numToInclude)));
	   if(b){
		   res = gB_[ind];
		   if(res==null){
			   for(int i=0; i<jB.length; i++){
			   jB[i].setSize(1200,200*(Math.min(numPerPage,numToInclude)));
			   gB_[ind] = Constants.pdf(0).equals("pdf") ? new PDFGraphics2D(new File(this.chartBF, "b"+ind+".pdf"),jB[i]) :
	      		   new EMFGraphics2D(new File(this.chartBF, "b"+ind+".emf"),this.jB[i]) ;//ImageConstants.PNG);
		         gB_[ind].setDeviceIndependent(Constants.plot()==1);
		         gB_[ind].startExport();
		         res = gB_[ind];
			   }
		      
		   }
	   }
	   else{
		   res = gR_[ind];
		   if(res==null){
			   for(int i=0; i<jR.length; i++){
			   jR[i].setSize(1200,200*(Math.min(numPerPage,numToInclude)));
			   gR_[ind] = Constants.pdf(0).equals("pdf") ? new PDFGraphics2D(new File(this.chartRF, "r"+ind+".pdf"),jR[i]) :
	      		   new EMFGraphics2D(new File(this.chartRF, "r"+ind+".emf"),this.jR[i]) ;//ImageConstants.PNG); 
			   gR_[ind].setDeviceIndependent(Constants.plot()==1);
		        gR_[ind].startExport();
		        res = gR_[ind];
			   }
		   }
	   }
	   return res;
	}*/

// AbstractVectorGraphicsIO     gB,gR;
    final static int numPerPage = 10;
    
    //AbstractVectorGraphicsIO[] gB_; //= new AbstractVectorGraphicsIO[];
   // AbstractVectorGraphicsIO[] gR_;// = new HashMap<Integer, AbstractVectorGraphicsIO>();

    
    
}
