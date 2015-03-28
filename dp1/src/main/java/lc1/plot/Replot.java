package lc1.plot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Stroke;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import lc1.dp.appl.CompareBasic;
import lc1.dp.data.classification.ROC;
import lc1.dp.data.classification.ROC.Line;
import lc1.dp.swing.Headless;
import lc1.util.Constants;

import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.PageConstants;
import org.jfree.chart.ChartPanel;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;



public class Replot {
static boolean useDirName =false;
static double[] overlap = new double[] {0.1, 0.1};
static List<String> replace = Arrays.asList(
		//""
		"all"
				.split(":"));
	static final Comparator fileC = new Comparator<File>(){

		@Override
		public int compare(File arg0, File arg1) {
			return arg0.getName().compareTo(arg1.getName());
		}
		
	};
	static String[][] titles = ROC.getSeriesTitles();
	static int index =0;
	public static boolean headless =!useDirName;
	static String[] types =
		
		("PennCNV(1M):cnvPartition(1M):QuantiSNP(1M):QuantiSNP(1M)(b):cnvHap(1M):cnvHap+snph(1M):" +
				"cnvHap(317k):cnvHap+snph(317k):"+
		"cnvHap(244k):cnvHap(185k):ADM(244k):ADM(185k):" +
	"244k+1M+185k+317k:244k:185k:1M:317k:244k+185k:1M+185k" +
	":185k+317k:244k+317k:1M+317k:244k+1M")
		.split(":");
	
	

	static String[] cols = 
		(
			"RED:GREEN:MAGENTA:PINK:DARK_GRAY:LIGHT_GRAY:DARK-BLUE:LIGHT-BLUE:" +
			"BLACK:CYAN:RED:GREEN: " +
			"MAGENTA:LIGHT_GRAY:CYAN:DARK_GRAY:DARK-BLUE:RED:ORANGE:LIGHT-GREEN:DARK-GREEN:LIGHT-BLUE:PINK"
		).split(":");
	
	static CompareBasic.ColorAdapt cola = new CompareBasic.ColorAdapt();
	static Map<String, String> m = new HashMap<String, String>();
	private static boolean sameAxis=index!=4;
	static{
		for(int i=0; i<types.length; i++){
			System.err.println(types[i]+" "+cols[i]);
			m.put(types[i],cols[i] );
		}
	}
	public static void main(String[] args){
		try{ 
			Constants.loess = null;
			Constants.gc = null;
			Constants.modify0 = new char[][] {{ '0', '1', '2' ,'3'}};
		
			
			
			File diro =new File(System.getProperty("user.dir"));
			File dir = new File(diro, args[0]+"/"+args[2]);
			File dir1 =  new File(diro, args[1]+"/"+args[2]);
			 File out  = new File(dir.getParentFile(), dir.getName()+"_new");
			out.mkdir();
		File 	out1   = new File(dir.getParentFile(), dir.getName()+"_new_noName");
			out1.mkdir();
			File[] f = dir.listFiles(new FileFilter(){

				@Override
				public boolean accept(File arg0) {
					return arg0.isDirectory();
				}
			});
			ROC.useSlope =true;
			Arrays.sort(f,fileC );
			String[][]names =
				index==0 ? 
				
				new String[][] {
					new String[] {" for detecting presence/absence of amplified probes in population",
								  " for detecting presence/absence of deleted probes in population"},
				 new String[] {" for detecting amplification breakpoints by sample",
								  " for detecting deletion breakpoints by sample"},
					new String[] {" for detecting amplified probes in individual sample",
					  " for detecting deleted probes in individual sample"}
					} : 
				new String[][] {{"",""},{"",""},{"",""}};
			ROC.topExtrax = 0.15;
			//ROC.bottomExtrax = 0.1;
			ROC.yperc = 0.95;
			ROC.useEnd = true;
			String todo = "graph1";
			for(int k=0; k<	f.length; k++){
				String nme = f[k].getName();
				if(index==0 && nme.equals("graph0")) ROC.topExtra =0.1;
				else ROC.topExtra = 0.2;
				if(nme.indexOf(todo)>=0){
					ROC.nameCurves = true;
					if(true){
						main(dir, dir1, out, f[k].getName(),names[k]);
					}
					else{
						main(dir, null, out, f[k].getName(), names[k]);
			        	main(dir1, null, out, f[k].getName(), names[k]);//"graph2");
					}
					if(Replot.headless){
						ROC.nameCurves = false;
						if(true){
							main(dir, dir1, out1, f[k].getName(),names[k]);
						}
						else{
							main(dir, null, out1, f[k].getName(), names[k]);
				        	main(dir1, null, out1, f[k].getName(), names[k]);//"graph2");
						}
					}
						
				}
			}
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static void main(File dir, File dir1,File out, String str, String[] title) throws Exception{
		File f1 = new File(dir, str);
		if( ! dir.exists() || ! f1.exists()){
			System.err.println("did not exist");
			return;
		}
		File out1 = new File(out, str+".png");
		
		plotDir(f1,dir1==null || !dir1.exists() ? null : new File(dir1,str),m, out1, title);
	}
	
	
	static void plotDir(File dir, File dir1, Map<String, String> m, File out, String[] title_) throws Exception{
		XYSeriesCollection[] xys  =readDir(dir,dir1==null || !dir1.exists() ? null :  dir1);
		
		List<Stroke> stroke = new ArrayList<Stroke>();
		List<Color> cols1 = new ArrayList<Color>();
		for(int i=0; i<xys[0].getSeriesCount(); i++){
			String key = xys[0].getSeries(i).getKey().toString();
			String col = m.get(key);
			System.err.println(key+" "+col);
			
			cols1.add(cola.getColor(col));
		stroke.add(new BasicStroke(0.01f));
//					BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,
	//				new float[] {5f,5f}, 1.0f));
		}
		
		Integer[] cn = new Integer[] {0,1};
	
		boolean[] log = ROC.getLogAxes()[index];
		boolean errorGraph = false;
		int max=1;
		int total = 1;
		double[] x_max1 = ROC.x_max()[index];
		boolean dots = true;
		String level="level";
		JPanel jp = new JPanel();
		String[] titles1 = titles[index];
		jp.setLayout(new BoxLayout(jp,
				BoxLayout.X_AXIS));
		for(int i=xys.length-1;i>=0; i--){
			ROC.overlap = overlap[i];
			 String title =  (useDirName ? dir.getAbsolutePath().split("ajhg1/")[1]+":":"")+titles1[0]+title_[i];
			                         ;//+(i==1 ? " (deletions)"2 : " (amplfications)");
		
			                      
		ChartPanel cp = ROC.getChartPanel(new XYSeriesCollection[] {xys[i]},
				cols1, stroke, "probeset", "result_set", cn, max, total, "chr",titles1[1],
				titles1[2], title, log[0], log[1], errorGraph, level, x_max1, dots);
		jp.add(cp);
		   analyse(xys[i],title, cp.getChart().getXYPlot().getDomainAxis().getRange().getUpperBound());
	//	cp.getChart().removeLegend();
		}
		if(sameAxis) CompareBasic.redrawAxes(new JPanel[]{jp});
		if(!headless){
			JFrame frame = new JFrame("");
			frame.setContentPane(jp);
			frame.setVisible(true);
			frame.pack();
			frame.show();
		}
		else{
			Headless jf;
			  jf = new Headless(jp);
				jf.pack();
				jf.setVisible(true);
				Properties p = new Properties();
				p.setProperty("PageSize", "A4");
				// ExportDialog ep;
				p.setProperty(PageConstants.ORIENTATION, PageConstants.LANDSCAPE);
				ImageGraphics2D g = new ImageGraphics2D(out, jp,
						ImageConstants.PNG);
				// ImageGraphics2D g = emf ? new EMFGraphics2D(out,cp_i.getSize()) :
				// new PDFGraphics2D(out1,
				// cp_i.getSize());
				g.setProperties(p);
				g.startExport();
				jp.print(g);
				g.endExport();
		}
		
	}
	
	private static void analyse(XYSeriesCollection seriesCollection,
			String title, double x) {
	double ref =Double.NaN;
	List<Line> l = new ArrayList<Line>();
		for(int i=0; i<seriesCollection.getSeriesCount(); i++){
			Line ser = new Line(seriesCollection.getSeries(i),false,false);
			if(ser.getKey().equals("cnvHap(1M)")){
				ref =  ser.getLt(x);
			}
			l.add(ser);
			
		}
		for(int i=0; i<seriesCollection.getSeriesCount(); i++){
		//	XYSeries ser = seriesCollection.getSeries(i);
			double d = l.get(i).getLt(x);
			double v = (ref-d)/d;
			String key = l.get(i).getKey().toString();
			System.err.println(x+ " COMPARE "+title+" "+key+" "+ref+" "+d+" "+((ref-d)/d));
		}
		
		
	}

	static int pos = 3;
	private static XYSeriesCollection[] readDir(File dir, File dir1) throws Exception{
		
		String[] pos1 = getPos(dir,pos);
		XYSeriesCollection[] ser = new XYSeriesCollection[pos1.length];
		for(int i=0; i<pos1.length; i++){
			XYSeries[] ser1 = readDir(dir, dir1, new PosFileFilter(pos, pos1[i], true));
			ser[i] = new XYSeriesCollection();
			for(int k=0; k<ser1.length; k++){
				ser[i].addSeries(ser1[k]);
			}
		}
		return ser;
	}
	

	private static XYSeries[] readDir(File dir, File dir1, FileFilter ff) throws Exception{
		SortedMap<String,  SortedMap<Double, Double>> m1 = new TreeMap<String,  SortedMap<Double, Double>> ();
		SortedMap<String,  SortedMap<Double, Double>> m2 = new TreeMap<String,  SortedMap<Double, Double>> ();
			double x = readDir(dir, ff,m1);
		  double x1 = dir1==null ? 0:readDir(dir1, ff,m2);
		  double scale = calcScale(m1,m2,x);
	if(m2!=null){
		for(Iterator<String> it = m2.keySet().iterator(); it.hasNext();){
			String key = it.next();
			if(m1.containsKey(key) && ! replace.contains(key) && ! replace.contains("all")){
		//		compare(key,m1.get(key),m2.get(key));
			}
			else{
				m1.put(key, m2.get(key));
			}
		}
	}
//		m1.putAll(m2);
		return getSeries(m1,Math.min(1.0, Math.max(x, x1)*5));
	}
	
	private static double calcScale(
			SortedMap<String, SortedMap<Double, Double>> targ,
			SortedMap<String, SortedMap<Double, Double>> src, double x) {
		// TODO Auto-generated method stub
		double sum=0;
		double cnt=0;
		for(Iterator<String> it = src.keySet().iterator();it.hasNext();){
			String key = it.next();
			if(targ.containsKey(key)){
				Line l1 = new Line(targ.get(key),logx, logy);
				Line l2 = new Line(src.get(key),logx, logy);
				sum+=l2.get(x)/l1.get(x);
				cnt++;
			}
		}
		return sum/cnt;
	}
	
	static boolean logx, logy;

	private static void compare(String str, SortedMap<Double, Double> map,
			SortedMap<Double, Double> map2) {
		ROC.Line  line1 = new ROC.Line(map,false,false);
		ROC.Line line2 = new ROC.Line(map2,false,false);
		double x = line1.mid();
		double y1 = line1.get(x);
		double y2 = line2.get(x);
		/*if(Math.abs(y1-y2)>0.01) {
			System.err.println("warning problem with "+str+" "+x+" "+y1+" "+y2);
		}*/
		double scale = y2/y1;
	/*	for(Iterator<Entry<Double,Double>> it = map.entrySet().iterator(); it.hasNext();){
			Entry<Double,Double> nxt = it.next();
			nxt.setValue(nxt.getValue()/scale);
		}*/
		System.err.println("scale is "+str+" "+scale);
	}


	private static XYSeries[] getSeries(Map<String, SortedMap<Double, Double>> m1, double x) {
		List<XYSeries> ser = new ArrayList<XYSeries>();
		for(Iterator<Entry<String,SortedMap<Double, Double>>> m = m1.entrySet().iterator(); m.hasNext();){
			Entry<String,SortedMap<Double, Double>> nxt = m.next();
			ser.add(getSeries(nxt.getValue(),x,nxt.getKey()));
		}
		return ser.toArray(new XYSeries[0]);
	}


	private static XYSeries getSeries(SortedMap<Double, Double> next,double x, String name) {
		ROC.Line line = new ROC.Line(next, false, false);
		XYSeries ser = new XYSeries(name);
		
		double lastKey = next.headMap(next.lastKey()).lastKey();
		
		
		for(Iterator<Entry<Double, Double>> ent = next.entrySet().iterator();ent.hasNext();){
			Entry<Double,Double> ent1 = ent.next();
			if(ent1.getKey()<=lastKey){
			ser.add(ent1.getKey(), ent1.getValue());
			}
			else{
				//double x = lastKey1+ (lastKey-firstKey);
				ser.add(x, line.get(x));
			}
			
		}
		return ser;
	}


	private static double  readDir(File dir, FileFilter ff,
			SortedMap<String, SortedMap<Double, Double>> m ) throws Exception{
		File[] in = dir.listFiles(ff);
		double x_ =0;
		//= new TreeMap<String, SortedMap<Double, Double>>();
		String pref = "cnvHap(";
		boolean allcnvhap = true;
		for(int i=0; i<in.length; i++){
			BufferedReader br = new BufferedReader(new FileReader(in[i]));
			String st = "";
			SortedMap<Double,Double> m1 = new TreeMap<Double,Double>();
			while((st = br.readLine())!=null){
				String[] str = st.split("\\s+");
				double x = Double.parseDouble(str[0]);
				double y = Double.parseDouble(str[1]);
				
				//if(x<0.5){
				m1.put(x,y);
				//}
			}
			br.close();
			if(m1.size()>2){
				double x1 = m1.headMap(m1.lastKey()).lastKey();
				if(x1>x_)x_ = x1;
			}
			String str = ROC.process(in[i].getName());
			if(m.containsKey(str)){
				m.put( str+"(b)", m.get(str));
				//str =;
			}
			m.put(str, m1);
			if(!str.startsWith(pref)){
				allcnvhap = false;
			}
		}
		if(allcnvhap){
			SortedMap<String, SortedMap<Double, Double>> m1 = new 
			TreeMap<String, SortedMap<Double,Double>>();
			
			for(Iterator<String> it =(new HashSet<String>( m.keySet())).iterator(); it.hasNext();){
				String nxt = it.next();
				//SortedMap<Double,Double> m1 = ;
				m1.put(nxt.substring(pref.length(),nxt.length()-1),
						m.remove(nxt));
				
			}
			m.putAll(m1);
		}	
		
		return x_;
	}
	static String[] getPos(File f, int pos){
		String[] lst = f.list();
		if(lst==null) {
			throw new RuntimeException(f.getName());
		}
		Set<String> s = new HashSet<String>();
		for(int i1=0; i1<lst.length; i1++){
			List<String> l = Arrays.asList(lst[i1].split("_"));
			List<String >l1 = l.subList(l.size()-pos, l.size());
			s.add(l1.get(0));
		}
		return s.toArray(new String[0]);
		}
	static class PosFileFilter implements FileFilter{
		final int pos;
		final String st;
		final boolean eq;
		
		public PosFileFilter(int pos, String st, boolean eq){
			this.pos = pos;
			this.st = st;
			this.eq = eq;
		}
		@Override
		public boolean accept(File arg0) {
		//	int ind1 = arg0.getName().indexOf("244k");
			//if(ind1>=0 && arg0.getName().indexOf("244k",ind1)>=0) return false;
			List<String> l = Arrays.asList(arg0.getName().split("_"));
			List<String >l1 = l.subList(l.size()-pos, l.size());
			if(l1.get(0).equals(st)) return eq;
			else return !eq;
		}
	}
	
}
