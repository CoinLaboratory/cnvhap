package summary;

import java.awt.Color;
import java.awt.Dimension;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;

import lc1.assoc.VirtualDataCollection1;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class QQ {
    XYSeries[] scores;
   public  ChartPanel chartPanel;
  //  ChartPanel chartPanel1;
    XYSeriesCollection coll_l;
 //   boolean log = false;
    XYSeries exp ;
    static  Dimension maxDim = new Dimension(500, 500);
    double mean;
    boolean log = true;
    public QQ(SortedMap<Double, Integer>[] entries, double tot, String name, double offset){
        exp = new XYSeries("expected");
        scores = new XYSeries[entries.length];
        double first = 1.0;
        double last = 0.0;
        for(int i=0; i<entries.length; i++){
            scores[i]= new XYSeries("scores_"+i);
            if(entries[i]==null || entries[i].size()==0) continue;
            double f = entries[i].firstKey();
            double f1 = entries[i].lastKey();
            if(f < first) first = f;
            if(f1>last) last = f1;
           double  cnt =0;
           for(Iterator<Map.Entry<Double, Integer>> it = entries[i].entrySet().iterator(); it.hasNext();){
              Map.Entry<Double, Integer> nxt = it.next();
              double nxtVal = nxt.getValue();
              double mid = (cnt+(nxtVal/2.0))/tot;
           double val = (double) nxt.getKey();
           if(log){
               val = -Math.log10(val);
               mid = -Math.log10(mid);
           }
           
              XYDataItem item =   new XYDataItem(mid, val);
                 
                    
             // System.err.println(item);
              scores[i].add(item);
              cnt+=nxt.getValue();
           }
        }
        if(log){
            first = - Math.log10(first);
            last = - Math.log10(last);
        }
        
        exp.add(first + offset, first+offset);
        exp.add(last,last);
       // double len = entries.size();
       
        chartPanel = new ChartPanel(getChart(!log, name));
        
         //chartPanel.setPreferredSize(maxDim);
       //  this.add(chartPanel);
    }
    
    JFreeChart getChart(boolean log, String name) {
        coll_l = new XYSeriesCollection();
        for(int i=0; i<scores.length; i++){
        coll_l.addSeries(scores[i]);
        }
        coll_l.addSeries(exp);
        final JFreeChart chart = ChartFactory.createXYLineChart(
               "QQ "+name,
             
                "Expected ", // domain axis label
                "Observed", // range axis label
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
        //plot
        // plot.setRenderer(new XYLineAndShapeRenderer());
        // ( (XYLineAndShapeRenderer)plot.getRenderer()).set.s
        plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
        // plot.getDomainAxis().setFixedAutoRange()
        final ValueAxis domainAxis = plot.getDomainAxis();
        if(log){
     plot.setDomainAxis(new LogarithmicAxis("expected"));
     plot.setRangeAxis(new LogarithmicAxis("observed"));
        }
      //  domainAxis.
        // domainAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
     //   plot.setDataset(1, coll_p);
    //    plot.mapDatasetToRangeAxis(1, 1);
      //  domainAxis.setRange(min[0], min[1]);
      //  final ValueAxis axis2 = new NumberAxis("Non deleted allele* allele freq");
     //   plot.setRangeAxis(1, axis2);
        final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
     //   final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
        // renderer2.setDrawShapes(true);
        // renderer2.setToolTipGenerator(new
        // StandardCategoryToolTipGenerator());
        plot.setRenderer(0, renderer1);
      //  plot.setRenderer(1, renderer2);
     
        return chart;
    }
    //QQ(List<Entry>){
        
    //}
}
