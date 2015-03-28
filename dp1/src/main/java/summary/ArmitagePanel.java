package summary;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Polygon;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.TableModel;

import lc1.assoc.OptionBuild;
import lc1.assoc.SNP;
import lc1.assoc.VirtualDataCollection1;
import lc1.dp.data.collection.ListTableModel;
import lc1.dp.swing.GenePanel;
import lc1.dp.swing.SummaryLogoPanel;
import lc1.stats.ChiSq;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.ChartEntity;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.entity.XYItemEntity;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYZDataset;

public class ArmitagePanel {
	
    public ChartPanel armitage;
    public ChartPanel qq;
   final  boolean doqq;
   String experiment;
   
   Color[][] color = 
	   new Color[][] {
		   new Color[] {Color.red.brighter(), Color.green.brighter(), Color.blue.brighter()},
		   new Color[] {Color.red.darker(), Color.green.darker(), Color.blue.darker()}
	   };
 //  double[] chrom_sizes;
    public ArmitagePanel(final VirtualDataCollection1 vdc, int expt_index,   boolean doqq, boolean probs) throws Exception{
    //	int[][]ind = new int[vdc.name.length][];
    	experiment = vdc.name[expt_index];
    	//for(int k=0; k<ind.length; k++){
    	//	ind[k] = vdc.getHeaderIndices(pheno, output, k);
    	//}
    	//chrom_sizes = chrom_length;
    	//final double xmin=1; final  double xmax=24;
        XYSeriesCollection[] coll = new XYSeriesCollection[2];  
        String[] fields = vdc.getFields(expt_index, vdc.PVALUE);
        XYSeries[][] series = new XYSeries[2][fields.length];
        SortedMap<Double, Integer> pvals[]=new SortedMap[fields.length];//new TreeMap<Double, Integer>();
        for(int i=0; i<series.length; i++){
        	coll[i] = new XYSeriesCollection();
        	for(int j=0; j<series[i].length; j++){
        		series[i][j] = new XYSeries(fields[j]);
        		coll[i].addSeries(series[i][j]);
        	}
        }
        for(int j=0; j<fields.length; j++){
        	pvals[j] = new TreeMap<Double, Integer>();
        }
        this.doqq = doqq;
        boolean pp = OptionBuild.probs;
       
        int tot=0;
        int coll_ind =0;
        for(int ij=0; ij<vdc.chroms.length; ij++){//
           outer: for(Iterator<SNP> it = vdc.chromIterator(ij); it.hasNext();){
            	try{
            		SNP snp = it.next();
            		String[] l = vdc.getPheno(snp.id, expt_index, ij);
            		String[] res = vdc.getPheno(l, expt_index, vdc.PVALUE);
            		for(int k=0; k<res.length; k++){
            			double p = Double.parseDouble(res[k]);
            			
            			series[coll_ind][k].add(ij+1+((double)snp.pos/(double)GenePanel.chr_length[ij]), p);	
            			if(doqq){
            				
            				Double p1 = probs ? p :Math.max(1e-20, ChiSq.chi2prob(1, p)); ;
            				Integer cnt = pvals[k].get(p1);
            				pvals[k].put(p1, cnt==null ? 1 : cnt+1);
            			}
            		}
            		tot++;
               }catch(Exception exc){
            	   exc.printStackTrace();
               }
            }
        coll_ind = 1-coll_ind;
            }
        
        //	for(int k=0; k<vdc[nonnull].name.length; k++){
        	//	if(ind[k].length>0){
        		//final int k1 = k;
              if(doqq) {
            	  qq= (new QQ(pvals, tot,vdc.name[expt_index], Math.log(OptionBuild.plotUpToPValue()))).chartPanel;
              }
      double  xmin =1;
      double        xmax = 23;
           //   if(coll[k]==null) continue;
                armitage =new ChartPanel(getChart(coll, "", expt_index==vdc.name.length-1,xmin, xmax ),
                		 1200, //width
                		 200, //height
                		 1200, //mindrawWidth
                         100, //mindrawHeight
                         1200, //maxDrawWith
                         400,//maxDrawHeight
                         ChartPanel.DEFAULT_BUFFER_USED,
                         true,  // properties
                         true,  // save
                         true,  // print
                         true,  // zoom
                         true   // tooltips		
                
                );
              
                armitage.getInsets().set(0, 0, 0, 0);
              
              armitage.addChartMouseListener(new ChartMouseListener(){

                public void chartMouseClicked(ChartMouseEvent arg0) {
                   JFreeChart chart =  arg0.getChart();
                  XYItemEntity entity =  ((XYItemEntity)arg0.getEntity());
                XYDataset ds =  chart.getXYPlot().getDataset();
                double x1 = arg0.getTrigger().getX();
                double y1 = arg0.getTrigger().getY();
              // Range range = chart.getXYPlot().getDomainAxis().getRange();
               if(entity==null){
	                EntityCollection entities = armitage.getChartRenderingInfo().getEntityCollection();
	              Iterator it =  entities.iterator();
	            
	             double dist = Double.POSITIVE_INFINITY;
	              while(it.hasNext()){
	            	  ChartEntity ent =(ChartEntity) it.next();
		              if(!(ent instanceof XYItemEntity)) continue;
		             XYItemEntity nxt = (XYItemEntity) ent;
	            	
	            	 double dist1 = getDistance(nxt, x1, y1);
	            	 if(dist1 < dist){
	            		 dist = dist1;
	            		 entity = nxt;
	            	 }
	              }
               }
             
                double x = ds.getX(entity.getSeriesIndex(), entity.getItem()).doubleValue();
                int fl = (int)Math.floor(x);
                int lo = (int)Math.round((x - fl)*(GenePanel.chr_length[fl-1]));
                for(int ki=0; ki<chart.getXYPlot().getDataset().getSeriesCount(); ki++){
                try{
                	double[] minmax = new double[] {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
                	if(vdc.dc[fl-1]!=null){
             //   XYZDataset xyz = vdc[fl-1].getBubbleChart(ki,lo, pheno, minmax);
              //  JFreeChart ch = getChart(xyz, pheno+" "+(ki==0 ? "no_cop" :(ki==1 ? "noA":"noB")), minmax);
              //  ChartPanel cp = new ChartPanel(ch);
              //  JFrame jf = new JFrame();
              //  jf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
              //  jf.getContentPane().add(cp);
              //  jf.pack();
              //  jf.setVisible(true);
                	}
                }catch(Exception exc){
                	exc.printStackTrace();
                }
                
                }
                 try{
                	 Map<String, Double> mco = new HashMap<String, Double>();
                	 Map<String, Double> mca = new HashMap<String, Double>();
                   TableModel[] tm = null;//vdc.getRS(lo, k1, mca, mco);
                   if(vdc.dc[fl-1]!=null){
                   SummaryLogoPanel slp = new SummaryLogoPanel(mco);
                   SummaryLogoPanel slp_ca = new SummaryLogoPanel(mca);
                  
                  Dimension dim = new Dimension(500,500);
                   showTables(slp, dim, "Controls");
                   showTables(slp_ca, dim, "Cases");
                   }
                   for(int i=0; i<tm.length; i++){
                       showTables(tm[i], null);
                   }
                   List<Object[]> genotypePicker = new ArrayList();
                   for(int i= -OptionBuild.contextsnps; i <=OptionBuild.contextsnps; i++ ){
                	   genotypePicker.add(new Object[] {"pos "+i, ""});
                   }
                   TableModel tm1 = new ListTableModel(genotypePicker, false);
                   ActionListener al = new ActionListener(){

						public void actionPerformed(ActionEvent e) {
							// TODO Auto-generated method stub
							
						}
                  
               	   
                  };
                   showTables(tm1, al);
                 }catch(Exception exc){
                     exc.printStackTrace();
                 }
                    
                }

                private double getDistance(XYItemEntity nxt, double x1, double y1) {
                	 Rectangle2D d =  nxt.getArea().getBounds2D();
                	 return Math.pow(d.getCenterX()-x1, 2)+Math.pow(d.getCenterY()-y1, 2);
				}

				private  void showTables(TableModel model, ActionListener al) {
                   JTable table = new JTable(model);
                   JFrame fr = new JFrame();
                   fr.getContentPane().setLayout(new BorderLayout());
                   fr.setDefaultCloseOperation(fr.DISPOSE_ON_CLOSE);
                   JScrollPane jB = new JScrollPane(table, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
                   fr.setSize(new Dimension(500,800));
                   fr.getContentPane().add(jB, BorderLayout.CENTER);
                   JButton jb = null;
                   if(al!=null){
                	   jb = new JButton("Apply");
                	   jb.addActionListener(al);
                	   fr.add(jb, BorderLayout.AFTER_LAST_LINE);
                   }
                   
                   fr.pack();
                   fr.setVisible(true);
                }
				private void showTables(JPanel jp, Dimension dim, String title) {
					
					jp.setSize(dim);
					jp.setMinimumSize(dim);
					jp.setPreferredSize(dim);
	                //   JTable table = new JTable(model);
	                   JFrame fr = new JFrame(title);
	                   fr.setDefaultCloseOperation(fr.DISPOSE_ON_CLOSE);
	                   JScrollPane jB = new JScrollPane(jp, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
	                   fr.setSize(new Dimension(500,800));
	                   fr.getContentPane().add(jB);
	                   fr.pack();
	                   fr.setVisible(true);
	                    
	                }

                public void chartMouseMoved(ChartMouseEvent arg0) {
                    // TODO Auto-generated method stub
                    
                }
                  
              });
            //    armitage[k].setPreferredSize(new Dimension(500,500));
        	//}
        	//}
    }
    
    JFreeChart getChart(XYZDataset coll_l, String name, double[] minmax) {
      //  ChartFactory.createBubbleChart(title, xAxisLabel, yAxisLabel, dataset, orientation, legend, tooltips, urls)
            final JFreeChart chart = ChartFactory.createBubbleChart(
                   name,
                 
                    "Number of copies ", // domain axis label
                    "phenotype", // range axis label
                    coll_l, // data
                    PlotOrientation.VERTICAL, false, // include legend
                    false, // tooltips?
                    false // URL generator? Not required...
                    );
            double len = minmax[1] - minmax[0];
            chart.getXYPlot().getRangeAxis().setRange(minmax[0]-len/5.0, minmax[1]+len/5.0);
            chart.getXYPlot().getDomainAxis().setRange(-1, 6);
           XYBubbleRenderer rend = new XYBubbleRenderer(XYBubbleRenderer.SCALE_ON_DOMAIN_AXIS);
        	chart.getXYPlot().setRenderer(rend);
           rend.setOutlinePaint(Color.white);
         //   chart.setBackgroundPaint(Color.white);
           // final XYPlot plot = chart.getXYPlot();
          //  plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
          //  final ValueAxis domainAxis = plot.getDomainAxis();
         //  if(coll_l.getSeries(0).getItemCount()<500){
          
            return chart;
            
    }

JFreeChart getChart(XYSeriesCollection[] coll_l, String name, boolean legend, double xmin, double xmax) {
    //  ChartFactory.createBubbleChart(title, xAxisLabel, yAxisLabel, dataset, orientation, legend, tooltips, urls)
        final JFreeChart chart = ChartFactory.createXYLineChart(
               name,
             
                "Location  ", // domain axis label
                "-log 10 pvalue", // range axis label
                coll_l[0], // data
                PlotOrientation.VERTICAL,true, // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
     
      XYPlot xy = chart.getXYPlot();
      for(int i=1; i<coll_l.length; i++){
    	  xy.setDataset(i, coll_l[i]);
      }
      if(!OptionBuild.probs){
    	 // XYSeriesCollection pvals = getPValLines(xmin, xmax, 1);
    	 // xy.setDataset(coll_l.length, pvals);
      }
        chart.getTitle().setFont(font10);
       // chart.getXYPlot().getDomainAxis().ge
        chart.getXYPlot().getRangeAxis().setLabelFont(font10);
       NumberAxis ra = (NumberAxis) chart.getXYPlot().getDomainAxis();
       ra.setAutoRange(false);
       ra.setRange(new Range(xmin, xmax));
        chart.getXYPlot().getDomainAxis().setLabelFont(font10);
        chart.getXYPlot().getRangeAxis().setTickLabelFont(font10);
        chart.getXYPlot().getDomainAxis().setTickLabelFont(font10);
        chart.getLegend().setItemFont(font10);
        chart.getLegend().setMargin(0,0	, 0, 0);
     
        AffineTransform at = new AffineTransform();
 	   at.setToScale(0.2, 0.2);
       final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
       final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
       final XYLineAndShapeRenderer renderer3 = new XYLineAndShapeRenderer();
       renderer3.setAutoPopulateSeriesFillPaint(false);
       renderer3.setAutoPopulateSeriesOutlinePaint(false);
       renderer1.setAutoPopulateSeriesFillPaint(false);
       renderer1.setAutoPopulateSeriesOutlinePaint(false);
       renderer3.setBaseFillPaint(Color.white);
       renderer3.setBaseOutlinePaint(Color.white);
       final Shape[] shapes = new Shape[3];
       int[] xpoints;
       int[] ypoints;

       // right-pointing triangle
       xpoints = new int[] {-3, 0, 3};
       ypoints = new int[] {-3, 3, -3};
       shapes[0] = new Polygon(xpoints, ypoints, 3);

       // vertical rectangle
       xpoints = new int[] {-3, 0, 3};
       ypoints = new int[] {3, -3, 3};
       shapes[1] = new Polygon(xpoints, ypoints, 3);
       // left-pointing triangle
       xpoints = new int[] {-3, -3, 3};
       ypoints = new int[] {-3, 3, 0};
       shapes[2] = new Polygon(xpoints, ypoints, 3);
   
       for(int k=0; k<coll_l[0].getSeriesCount() && k<3; k++){
    	   Shape shape =  renderer1.getSeriesShape(k);
    	  renderer1.setSeriesShape(k,shapes[k]);
    	  renderer1.setSeriesFillPaint(k, Color.white);// notify)
    	  renderer1.setSeriesOutlinePaint(k, color[0][k]);
    	  renderer3.setSeriesShape(k,shapes[k]);
    	  renderer3.setSeriesFillPaint(k, Color.white);// notify)
    	  renderer3.setSeriesOutlinePaint(k, color[1][k]);
    	//   renderer1.setSeriesShape(k, at.createTransformedShape(shape));
       }
       renderer2.setLinesVisible(true);
       renderer2.setShapesVisible(false);
       renderer2.setBasePaint(Color.black);
       for(int i=0; i<coll_l.length; i+=2){
    	   chart.getXYPlot().setRenderer(i, renderer1);
    	   chart.getXYPlot().setRenderer(i+1, renderer3);
       }
       chart.getXYPlot().setRenderer(coll_l.length, renderer2);
  /*  if(false && coll_l.getSeries(0).getItemCount()>3000){
    	renderer1.setShapesVisible(false);
    }
    else{*/
    	renderer1.setLinesVisible(false);
    	renderer3.setLinesVisible(false);
  //  }
      //  }
        return chart;
    }

double[] x = process(
		"1  2  3  4  5  " +
		"6  7  8  9 10 " +
		"11 12 13 14 15 " +
		"16");
double[] y = process(
		"2.705543  6.634897 10.827566 15.136705 19.511421 " +
		"23.928127 28.373987 32.841253 37.324893 41.821457" +
		" 46.328468 50.844169 55.365559 59.890642 64.101097" +
		" 68.412843");
private XYSeriesCollection getPValLines(double xmin, double xmax, int degf) {
	//double[] d = new double[] {1e-3, 1e-5};
	XYSeriesCollection xys = new XYSeriesCollection();
	for(int i=0; i<x.length; i++){
		//double y = ChiSquareDistribution.quantile(1.0-d[i], degf);
		XYSeries series = new XYSeries("1e-"+x[i]);
		series.add(xmin, y[i]);
		series.add(xmax, y[i]);
		xys.addSeries(series);
	}
	return xys;
}
private double[] process(String string) {
	String[] str = string.split("\\s+");
	double[] res = new double[str.length];
	for(int i=0; i<res.length; i++){
		res[i] = Double.parseDouble(str[i]);
	}
	return res;
}

public static  Font font10= new Font("Arial", Font.PLAIN, 8);
}
