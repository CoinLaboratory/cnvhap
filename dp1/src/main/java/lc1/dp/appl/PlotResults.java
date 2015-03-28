package lc1.dp.appl;

import java.awt.Color;
import java.awt.Stroke;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import javax.swing.BoxLayout;
import javax.swing.JPanel;

import lc1.dp.swing.Headless;

import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.PageConstants;
import org.jfree.chart.ChartPanel;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class PlotResults {
	static boolean emf = false;
	
	public static void main(String[] args){
		run1(new File(System.getProperty("user.dir")));
	}
	public static void run1(File user){
	//	File user = 
		File[] fs = user.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.isDirectory() ;
			}
			
		});
		for(int i=0; i<fs.length; i++){
			if(user.getName().indexOf("244")>=0 || user.getName().indexOf("1M")>=0){
				if(fs[i].getName().indexOf("Cum")>=0) run(fs[i]);
				
			}
			else  {
				run1(fs[i]);
			}
			
		}
	}
	
	public static void run(File user){
		System.err.println("plotting for "+user.getAbsolutePath());
		try{
			
			File[] fs = user.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.isDirectory() && pathname.getName().indexOf("indiv")>=0;
				}
				
			});
			ChartPanel[] cp = new ChartPanel[fs.length];
			lc1.dp.swing.Headless jf; 
			  JPanel jp = new JPanel();
			  jf = new Headless(jp);//.setContentPane(jp);
			  jp.setLayout(new BoxLayout(jp, BoxLayout.X_AXIS));
			for(int i=0; i<fs.length; i++){
			PlotResults pr = new PlotResults(fs[i]);
			cp[i] = pr.cp;
			jf.getContentPane().add(cp[i]);
			}
			
			  jf.setVisible(true);
			  jf.pack();
				File out1 = new File(user, "graph2.png");
      		  Properties p = new Properties();
                p.setProperty("PageSize","A4");
            //    ExportDialog ep;
                p.setProperty(PageConstants.ORIENTATION, PageConstants.LANDSCAPE);
                ImageGraphics2D g  = new ImageGraphics2D(out1,jp, ImageConstants.PNG); 
                	
                g.setProperties(p); 
                g.startExport(); 
                jp.print(g);
                g.endExport();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static String[][] expt = new String[][]{ "1M_PennCNV 244k_ADM2 cnvPartion 244k_ADM2 185k_ADM2 Imputed244k".split("\\s+"),
		"1M_cnvHap	   317k_cnvHap	317k_1M_cnvHap".split("\\s+"),
		"244k_cnvHap	185k_cnvHap	185k_244k_cnvHap".split("\\s+"),
		"244k_1M_cnvHap 185_317_cnvHap 185_1M_cnvHap 244_317_cnvHap".split("\\s+")};
		
	public static String[][] expt1 = new String[][]{ "PennCNV_1M ADM2_244k cnvPartion_1M ADM2_244k ADM2_185k Imputed244k".split("\\s+"),
		"cnvHap_1M	   cnvHap_317k	cnvHap_317k_1M".split("\\s+"),
		"cnvHap_244k	cnvHap_185k	cnvHap_185k_244k".split("\\s+"),
		"cnvHap_244k_1M cnvHap_185k_317k cnvHap_185k_1M cnvHap_244k_317k".split("\\s+")};
		
	
	static Map<String, Color> colors = new HashMap<String, Color>();
	/*static {
		String[][] str = CompareBasic.expt;
		String[][] col = CompareBasic.col;
		try{
		for(int i=0; i<str.length; i++){
			for(int j=0; j<str[i].length; j++){
				Color colo = 	 (Color) Color.class.getField(col[i][j]).get(null);
				colors.put(str[i][j], colo);
			}
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}*/
	File f;
	XYSeries[] series;
	XYSeriesCollection collection = new XYSeriesCollection();
	ChartPanel cp;
	List<Color> result_set_color = new ArrayList<Color>();
	List<Stroke> result_set_stroke = new ArrayList<Stroke>();
	PlotResults(File f){
		this.f = f;
		for(int i=0; i<expt.length; i++){
			for(int j=0; j<expt[i].length; j++){
				m.put(expt[i][j].toLowerCase(), expt1[i][j]);
				System.err.println("put "+expt[i][j]+" "+expt1[i][j]);
			}
		}
		try{
			File[] fs = f.listFiles();
			series = new XYSeries[fs.length];
			for(int i=0; i<series.length; i++){
				series[i] = getSeries(fs[i]);
				result_set_color.add(colors.get(fs[i].getName()));
			//	result_set_stroke.add(
        				//str.length>3 ? 
        		//new BasicStroke(0.8f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,parse(str[3].split(":")) , 1.0f)
        	//	new BasicStroke(0.8f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,CompareBasic.generate(5) , 1.0f));
        	
				collection.addSeries(series[i]);
			}
			/*cp = ROC.getChartPanel(new XYSeriesCollection[] {collection}, result_set_color, result_set_stroke, "probeset", 
					"name",new Integer[]{f.getName().indexOf("deletion")>=0 ? 0 : 3}, 
					f.getName().indexOf("pop")>0 ? 1 : 0
					,0,
					"all", 
					"Minimum size of aberration", 
					"True positives", 
					"Cumulative true positive vs minimum aberration size", 
					true,
					true);*/
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	Map<String, String> m = new HashMap();
	private XYSeries getSeries(File file) throws Exception{
		String nme = file.getName();
		String nme1 = m.get(nme.toLowerCase());
		System.err.println(nme+" "+nme1);
		if(nme1==null) nme1 = nme;
		XYSeries series = new XYSeries(nme1);
		BufferedReader br = new BufferedReader(new FileReader(file));
		String st = "";
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			double x = Double.parseDouble(str[0]);
			double y = Double.parseDouble(str[1]);
			series.add(x, y);
		}
		br.close();
		return series;
	}
	
}
