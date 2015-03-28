package lc1.possel;

import java.awt.Color;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class PlotScores  extends JFrame{
	
	public static void main(String[] args){
	    File user = new File(System.getProperty("user.dir"));
	    File norm = user.listFiles(new FilenameFilter1("score_1"))[0];
	    File dels = user.listFiles(new FilenameFilter1("score_0"))[0];
	    File ampl = user.listFiles(new FilenameFilter1("score_2"))[0];
	    
	    String[] str = new String[] {"1"};
	    for(int i=0; i<str.length; i++){
	    	File pritch = new File(user.getParent(), "ihs/ceu.ch"+str[i]);
	    	PlotScores ps = new  PlotScores(norm, dels, ampl, pritch, str[i], 158, 159);
	    }
	}
	static  Dimension maxDim = new Dimension(500, 500);
	
	XYSeries dels=new XYSeries("deletion allele");;
	XYSeries pritch=new XYSeries("pritchard");;
	XYSeries ampl=new XYSeries("amplified allele");;
	XYSeries norm=new XYSeries("normal allele");;
	
	String chr;
	 ChartPanel chartPanel;
	 final int start, end;
	PlotScores(File norm, File del, File ampl, File pritch, String chr, int start, int end){
		super(chr);
		this.chr = chr;
		this.start = start*1000000;
		this.end = end*1000000;
		int cnt0 = read(del, this.dels, chr);
		int cnt1 = read(norm, this.norm, chr);
		int cnt2 = read(ampl, this.ampl, chr);
		System.err.println("cnts are "+cnt0+" "+cnt1+" "+cnt2);
		readPritch(pritch, this.pritch);
		chartPanel = new ChartPanel(getChart());
	       
	        chartPanel.setPreferredSize(maxDim);
	        this.getContentPane().add(chartPanel);
	        this.pack();
	        this.setVisible(true);
	}

	private  int  read( File del, XYSeries dels2, String chr) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(del));
			String st = "";
			br.readLine();
			double prev_sc = -100;
			int prev_loc=0;
			int cnt=0;
			while((st = br.readLine())!=null){
				String[] str = st.split("\\s+");
//				if(str[0].equals(chr)){
					int loc = Integer.parseInt(str[3]);
					
					double ihs = Double.parseDouble(str[7]);
					if(loc >=start && loc<end && str[0].equals(chr)){
						dels2.add(loc, ihs);
					
					}
					if(Math.abs(ihs - prev_sc) > 0.1 || (loc - prev_loc)>1000*100) cnt++;
					prev_loc = loc;
					prev_sc = ihs;
	//			}
			}
			return cnt;
		}catch(Exception exc){
			exc.printStackTrace();
		}
		return 0;
	}
	private  void readPritch( File del, XYSeries dels2) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(del));
			String st = "";
			br.readLine();
			while((st = br.readLine())!=null){
				String[] str = st.split("\\s+");
				if(!str[3].equals("-")){
					int loc = Integer.parseInt(str[1]);
					if(loc >=start && loc<end){
					double ihs = Double.parseDouble(str[3]);
					dels2.add(loc, ihs);}
				}
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	
	  JFreeChart getChart() {
	        XYSeriesCollection coll_l = new XYSeriesCollection();
	        coll_l.addSeries(this.dels);
	       coll_l.addSeries(this.ampl);
	       coll_l.addSeries(this.pritch);
	     coll_l.addSeries(this.norm);
	 
	        final JFreeChart chart = ChartFactory.createXYLineChart(
	                "IHH for  "+chr,
	                "Location ", // domain axis label
	                "IHH", // range axis label
	                coll_l, // data
	                PlotOrientation.VERTICAL, true, // include legend
	                true, // tooltips?
	                false // URL generator? Not required...
	                );
	        
	        ((NumberAxis) ((XYPlot) chart.getPlot()).getRangeAxis())
	                .setAutoRangeIncludesZero(false);
	        chart.setBackgroundPaint(Color.white);
//	        chart.getLegend().setAnchor(org.jfree.chart.Legend.SOUTH);
	      //  final XYPlot plot = chart.getXYPlot();
	        // plot.setRenderer(new XYLineAndShapeRenderer());
	        // ( (XYLineAndShapeRenderer)plot.getRenderer()).set.s
	      //  plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
	        // plot.getDomainAxis().setFixedAutoRange()
	      //  final ValueAxis domainAxis = plot.getDomainAxis();
	        // domainAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
	     //   plot.setDataset(1,  coll_l);
	     //   plot.mapDatasetToRangeAxis(1, 1);
	    //  //  domainAxis.setRange(min[0], min[1]);
	      //  final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
	     //   final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
	        // renderer2.setDrawShapes(true);
	        // renderer2.setToolTipGenerator(new
	        // StandardCategoryToolTipGenerator());
	     //   plot.setRenderer(0, renderer1);
	        return chart;
	    }
}
