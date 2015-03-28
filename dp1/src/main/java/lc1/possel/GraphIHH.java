package lc1.possel;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import javax.swing.JFrame;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class GraphIHH extends JFrame{

    XYSeries scoresDerived;
    XYSeries scoresAncestral;
    List<XYTextAnnotation> annotationD = new ArrayList<XYTextAnnotation>();
    List<XYTextAnnotation> annotationA = new ArrayList<XYTextAnnotation>();
    ChartPanel chartPanel;
    List<String >chrom;
   static  Dimension maxDim = new Dimension(500, 500);
   static double inc = 0.04;
   static double inc1 = 0.00;
    public GraphIHH(List<Double>[][] in,  Double res, Double p_val, int i1, int i2,  List<Integer> loc, File chr, Double maf_der, Double maf_anc, boolean show)
    throws Exception{
        super("IHH graph "+chr.getName().split("\\.")[0]+" "+loc.get(i1)+"-"+loc.get(i2));
        this.scoresDerived = new XYSeries("deletion allele");
        this.chrom = EHHMultiChromosome.split(chr);
       List<Double> derL = in[0][0];
       List<Double> ancL = in[0][1];
       List<Double> derR = in[1][0];
       List<Double> ancR = in[1][1];
        this.scoresAncestral = new XYSeries("non deletion allele");
        if(derL!=null){
            for(int k=0; k<derL.size(); k++){
                int pos = loc.get(i1 - k);
                double yD = (derL.get(k)/derL.get(0))*maf_der;
                double yA = (ancL.get(k)/ancL.get(0))*maf_anc;
                scoresDerived.add(pos, yD);
                scoresAncestral.add(pos, yA);
              //  this.annotationD.add(new XYTextAnnotation(""+derL.ch(k),pos, yD+inc1*maf_der));
               // this.annotationA.add(new XYTextAnnotation(""+ancL.ch(k),pos, yA+inc*maf_anc));
            }
        }
        if(derR!=null){
        for(int k=0; k<derR.size(); k++){
            int pos = loc.get(i2 + k);
            double yD = (derR.get(k)/derR.get(0))*maf_der;
            double yA = (ancR.get(k)/ancR.get(0))*maf_anc;
            scoresDerived.add(pos, yD);
            scoresAncestral.add(pos,yA);
            //this.annotationD.add(new XYTextAnnotation(""+derR.ch(k),pos, yD+inc1*maf_der));
            //this.annotationA.add(new XYTextAnnotation(""+ancR.ch(k),pos, yA+inc*maf_der));
        }
        }
        chartPanel = new ChartPanel(getChart(res, p_val, loc.get(i1), new int[] {loc.get(0), loc.get(loc.size()-1)}));
       // chartPanel.setName("Likelihood chart ");
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
       
        chartPanel.setPreferredSize(maxDim);
        this.getContentPane().add(chartPanel);
        this.pack();
        this.setVisible(show);
        
     
    }
    
    
    XYSeriesCollection coll_l;

    XYSeriesCollection coll_p;

    public String getInfo(){
        return chrom.get(0)+"-"+chrom.get(1)+"-"+chrom.get(2);
    }
    JFreeChart getChart(Double res, Double p_val, int pos, int min[]) {
        coll_l = new XYSeriesCollection();
        coll_p = new XYSeriesCollection();
        coll_l.addSeries(scoresDerived);
        coll_p.addSeries(scoresAncestral);
        Object[] toPrint = new Object[] {res, p_val};
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "IHH for deletion "+getInfo()+" at "+pos+
                String.format("\nIHH: %5.3f pval:  %5.3g", toPrint),
                "Location ", // domain axis label
                "Deleted allele  EHH * allele freq", // range axis label
                coll_l, // data
                PlotOrientation.VERTICAL, true, // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
        
        ((NumberAxis) ((XYPlot) chart.getPlot()).getRangeAxis())
                .setAutoRangeIncludesZero(false);
        chart.setBackgroundPaint(Color.white);
//        chart.getLegend().setAnchor(org.jfree.chart.Legend.SOUTH);
        final XYPlot plot = chart.getXYPlot();
        // plot.setRenderer(new XYLineAndShapeRenderer());
        // ( (XYLineAndShapeRenderer)plot.getRenderer()).set.s
        plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
        // plot.getDomainAxis().setFixedAutoRange()
        final ValueAxis domainAxis = plot.getDomainAxis();
        // domainAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
        plot.setDataset(1, coll_p);
        plot.mapDatasetToRangeAxis(1, 1);
      //  domainAxis.setRange(min[0], min[1]);
        final ValueAxis axis2 = new NumberAxis("Non deleted allele* allele freq");
        plot.setRangeAxis(1, axis2);
        final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
        final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
        // renderer2.setDrawShapes(true);
        // renderer2.setToolTipGenerator(new
        // StandardCategoryToolTipGenerator());
        plot.setRenderer(0, renderer1);
        plot.setRenderer(1, renderer2);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
        for(int i=0; i<this.annotationD.size(); i++){
            renderer1.addAnnotation(annotationD.get(i));
        }
        for(int i=0; i<this.annotationA.size(); i++){
            renderer2.addAnnotation(annotationA.get(i));
        }
        return chart;
    }
    public static void plot(List<Double>[][] in, double sc, double pval, int i, int i2, List<Integer> loc, File f, Double maf, boolean show)  throws Exception{
        GraphIHH ihhG = new GraphIHH(in, sc,pval, i, i2,loc, f, maf, 1.0-maf, show);
        Properties p = new Properties();
        p.setProperty("PageSize","A4");
        VectorGraphics g = new PDFGraphics2D(new File(ihhG.getInfo()+".pdf"), ihhG.maxDim); 
        g.setProperties(p); 
        g.startExport(); 
       ihhG.print(g); 
        g.endExport();
        
    }

    
}
