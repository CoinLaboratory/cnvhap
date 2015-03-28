package lc1.dp.data.classification;

import java.awt.Color;
import java.awt.Font;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.logging.Logger;

import javax.swing.JFrame;

import lc1.dp.appl.CompareBasic;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LightWeightDataCollection;
import lc1.dp.data.collection.LightWeightMergedDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.possel.HLALD;
import lc1.stats.HWECalc;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import lc1.util.PrintLatexTable;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;

public class ROC {
public	 static  boolean countRegs = false;
public static boolean combine =false;
public static boolean addOnlyLongestOverlap = true;
public static boolean useLimitFactor = false;
public static boolean drawNames = false;
public static String 	headless = System.getProperty("java.awt.headless");
public static boolean legend = false;//headless==null || !headless.equals("true");
public static boolean useSlope = true;
public static double resolution = 1.0;
	static boolean printDetails =!CompareBasic.combine;
	public final int minNumberTrue, minLengthTrue;
	int midB, midC;
	static public int min_len=0;
static public int[] len_bands = new int[] {Integer.MAX_VALUE-1};
	final int index = 0; // index of gold standard
	double prob_upper = 0.5;;
	double prob_lower = 0.5;
	double hwe_thresh = 0.0;
	final int centromere;
	
	
	static boolean  setShape = true;
	
	static Comparator len = new Comparator<int[]>(){

		@Override
		public int compare(int[] o1, int[] o2) {
			int len1 = o1[1]-o1[0];
			int len2 = o2[1]-o2[0];
			if(len1!=len2) return len1>len2 ? -1 : 1;
			else return 0;
					
		}
		
	};
	static Comparator len2 = new Comparator<double[]>(){

		@Override
		public int compare(double[] o1, double[] o2) {
			double len1 = o1[0];
			double len2 = o2[0];
			
			if(len1!=len2) return len1<len2 ? -1 : 1;
			else return 0;
					
		}
		
	};
	
	static Comparator len1 = new Comparator<List<Integer>>(){

		@Override
		public int compare(List<Integer> o1, List<Integer> o2) {
			int len1 = o1.get(1)-o1.get(0);
			int len2 = o2.get(1)-o2.get(0);
			if(len1!=len2) return len1>len2 ? -1 : 1;
			else return 0;
					
		}
		
	};
	
public static  List<List<Integer>> calcOverlap(List<List<Integer>> ov1, List<List<Integer>>ov2, List<List<Integer>>notfound){
		List<List<Integer>> map = new ArrayList<List<Integer>>();
		
		
		for(Iterator<List<Integer>> it = ov1.iterator(); it.hasNext();){
			List<Integer> v = it.next();
			int start = v.get(0);
			int end = v.get(1);
			boolean found = false;
			//SortedMap<Integer, Integer> head = ov2.headMap(end);
			List<List<Integer>> d = new ArrayList<List<Integer>>();
			for(Iterator<List<Integer>> it1 = ov2.iterator(); it1.hasNext();){
				List<Integer> v1 = it1.next();
				int start1 = v1.get(0);
				int end1 = v1.get(1);
				int new_end = Math.min(end1, end);
				int new_start = Math.max(start1, start);
				
				if(new_end>= new_start){
					found = true;
					System.err.println("LENGTHS "+(end-start)+" "+(end1-start1)+" "+(new_end - new_start));
					d.add(Arrays.asList(new Integer[]{new_start, new_end}));
				
				}
				
			}
			if(!found && notfound!=null) notfound.add(Arrays.asList(new Integer[] {start, end}));
			if(addOnlyLongestOverlap && d.size()>1){
				Collections.sort(d,len1);
			
				map.add( d.get(0));
			}
			else{
				map.addAll(d);
			}
		}
		return map;
	}

public static XYSeries getLengthHistogram(SortedMap<Integer, Integer> lengthCount, String name){
	XYSeries series = new XYSeries(name);
	int sum=0;
	for(Iterator<Entry<Integer, Integer>> it= lengthCount.entrySet().iterator();  it.hasNext();){
		Entry<Integer, Integer> ent = it.next();
		sum+=ent.getValue();
		series.add(ent.getKey().intValue(), sum);
	}
	return series;
}


public static XYSeries getLengthDist(Inner inner, int type, Inner limitFactor){
	SortedMap<Integer, Integer> length = new TreeMap<Integer, Integer>();
	if(type==0){
	for(Iterator<SortedMap<Integer, Integer>> it = inner.l_max.regions.values().iterator(); it.hasNext();){
		SortedMap<Integer, Integer> map = it.next();
		addToLengthDist(getList(map), length);
	}
	}
	else{
		for(Iterator<String> it = inner.l_max.regions.keySet().iterator(); it.hasNext();){
			String key = it.next();
			List<List<Integer>> limit_region = limitFactor==null ? null : getMergedLocations(limitFactor,key);
		for(int k=0; k<inner.l_indiv.length; k++){
		
		
				SortedMap<Integer, Integer> map = inner.l_indiv[k].regions.get(key);
				List<List<Integer>> l1 = getList(map);
				if(limit_region!=null){
					l1 =calcOverlap(l1, limit_region, null);
				}
				addToLengthDist(l1, length);
			}
		}
	}
	return getLengthHistogram(length, inner.name);
}

public static XYSeries getLengthDist(Inner inner1, Inner inner2, int type, Inner limitFactor){
	SortedMap<Integer, Integer> length = new TreeMap<Integer, Integer>();
	if(type==0){
		
		for(Iterator<String> it = inner1.l_max.regions.keySet().iterator(); it.hasNext();){
			String key = it.next();
			SortedMap<Integer, Integer> m1 = inner1.l_max.regions.get(key);
			SortedMap<Integer, Integer> m2 = inner2.l_max.regions.get(key);
			
			List<List<Integer>> notfound = new ArrayList<List<Integer>> ();
			if(ROC.useLimitFactor){
				addToLengthDist(calcOverlap(getMergedLocations(inner1,key), getMergedLocations(inner2,key), notfound), length);
			}
			else{
				addToLengthDist(calcOverlap(getList(m1), getList(m2), notfound), length);
			}
		
		}
	}
	else{
		
		
		
		
			
		for(Iterator<String> it = inner1.l_max.regions.keySet().iterator(); it.hasNext();){
			String key = it.next();
			List<List<Integer>> limit_region = limitFactor==null ? null : getMergedLocations(limitFactor,key);
			for(int k=0; k<inner1.l_indiv.length; k++){
				String ind = inner1.l_indiv[k].name;
				SortedMap<Integer, Integer> m1 = inner1.l_indiv[k].regions.get(key);
				SortedMap<Integer, Integer> m2 = inner2.l_indiv[k].regions.get(key);
				List<List<Integer>> notfound = new ArrayList<List<Integer>>();
				List<List<Integer>> l1 = getList(m1);
				List<List<Integer>> l2 = getList(m2);
				if(limit_region!=null){
					l1 = calcOverlap(l1, limit_region, null);
					l2 = calcOverlap(l2, limit_region, null);
				}
				addToLengthDist(calcOverlap(l1, l2, notfound), length);
		
			}
		}
	}
	return getLengthHistogram(length, inner1.name+"_"+inner2.name);
}

public static List<List<Integer>> getMergedLocations(Inner inner1, String key){
	List<List<Integer>> res = new ArrayList<List<Integer>>();
	for(int k=0; k<inner1.l_indiv.length; k++){
		res.addAll(getList(inner1.l_indiv[k].regions.get(key)));
	}
	return res;
}
static List<List<Integer>> getList(Map<Integer, Integer> m){
	List<List<Integer> > l = new ArrayList<List<Integer>>();
	for(Iterator<Entry<Integer, Integer>> it = m.entrySet().iterator();it.hasNext();){
		Entry<Integer, Integer> nxt = it.next();
		l.add(Arrays.asList( new Integer[] {nxt.getKey(),	nxt.getValue()}));
	
	}
	Collections.sort(l, len1);
	return l;
	
}

public static XYSeriesCollection getLenghtDists(Inner[] inn, int type){
	List<XYSeries> l = new ArrayList<XYSeries>();
	//Map m1 = inn[0].l_max.regions;
	
	//boolean eq = m1.equals(m2);
	Inner limitFactor = useLimitFactor && inn.length>1 ? inn[inn.length-1] : null;
	// System.err.println(eq);
	for(int i=0; i<inn.length; i++){
		l.add(getLengthDist(inn[i], type, limitFactor));
	}
	if(combine){
		for(int i=0; i<inn.length; i++){
			for(int j=i+1; j<inn.length; j++){
			l.add(getLengthDist(inn[i], inn[j], type, limitFactor));
			}
		}
	}
	XYSeriesCollection xys = new XYSeriesCollection();
	for(int i=0; i<l.size(); i++){
	xys.addSeries(l.get(i));
	}
	return xys;
}


public static  void addToLengthDist(List<List<Integer>> list, Map<Integer, Integer> lengthCount){
			for(Iterator<List<Integer>> it = list.iterator(); it.hasNext();){
				List<Integer> next = it.next();
				int len =next.get(1) - next.get(0);
				Integer cnt = lengthCount.get(len);
				lengthCount.put(len, cnt==null ? 1 : cnt+1);
			}
}
	
	//Constants.loess = null;;
	public static int thresh1 = Integer.MAX_VALUE;
	public static double[] fpr = new double[] {0.01,0.05, 0.25, };
	//static Map<Doutble, Double> replace = new HashMap<Double, Double>();
	public static boolean avgForR2 =true;
	/*static{
		replace.put(6.0, 4.0);
		replace.put(5.0, 4.0);  //for comparing to ref
	}*/
	public static String getCNString(Integer[] cn){
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < cn.length; i++) {
			sb.append(cn[i] + ".");
		}
		return sb.toString();
	}
	final static int ext = 5;
	public static final Comparator<? super Double> comp_back = new Comparator<Double>(){

		@Override
		public int compare(Double o1, Double o2) {
			return -1*o1.compareTo(o2);
		}
		
	};
	public static boolean useEnd = false;
	public static  double yperc = 0.98;
	public static String[][] getSeriesTitles() {
		String[][] res = new String[fpr.length + ext][];
		String rate = makeAsRate ? "rate" : "count";
		String freq = makeAsRate ? "frequency" : "count";
		String number = makeAsRate ? "Fraction" : "Number";
		res[0] = new String[] { "ROC ", "False postive "+rate,
				"True positive "+rate };
		res[1] = new String[] { "Error rate ", "Probability of aberration",
				"Fraction of all samples in probability bin with aberration" };
		res[2] = new String[] {"Cumulative Count", "Length(bp)", "Cumulative region count"};
		res[3] = new String[] {"Concordance of genotype calls", 
				number+" of calls"
				
				, "Discordance "+rate};
		res[4] = new String[] {"Correlation between real and predicted CN", "r^2", "Cumulative "+freq};
		
		for (int i = 0; i < fpr.length; i++) {

			res[i + ext] = new String[] {
					"Cumulative TP vs Size of Aberration at FPR of" + fpr[i]
							+ ". ", "Minimum size of aberration",
					"True Positive" };
		}
		// Arrays.fill(titles, );

		return res;// new String[][] {titles, x, y};
	}
	public static double[][] x_max(){
		double[][] res = new double[fpr.length + ext][];
		res[0] = new double[] {-0.001,0.05};
		res[1] = new double[] {-0.001,0.1};
		res[2] = new double[] {-0.002,0.1};
		res[3] = new double[] {-0.02,1.02};
		res[4] = new double[] {-0.02,1.02};
		for (int i = 0; i < fpr.length; i++) {
			res[i+ext] =new double[]{-0.001,0.1};
		}
		return res;
		
	}
	public static boolean makeAsRate = true;
public static boolean logconcy =true;
	public static boolean[][] getLogAxes() {
		boolean[][] res = new boolean[fpr.length + ext][];
		res[0] = new boolean[] {false, false };
		res[1] = new boolean[] { false, false};
		res[2] = new boolean[] {true, false};
		res[3] = new boolean[] {true, logconcy};
		res[4] = new boolean[] {false, false};
		res[5] = new boolean[] {false, false};
		for (int i = 0; i < fpr.length; i++) {
			res[i + ext] = new boolean[] { false, false };
		}
		return res;
	}
	public static boolean[] dots() {
		boolean[] res = new boolean[fpr.length + ext];
		res[0] =true;
		res[1] = true;
		res[2] = true;
		res[3] = true;
		res[4] = true;
		for (int i = 0; i < fpr.length; i++) {
			res[i + ext] = true;
		}
		return res;
	}

	public static int[] getTypes() {
		int[] res = new int[fpr.length + ext];
		Arrays.fill(res, 1);
		res[0] = len_bands.length;
		res[1] = 2;
		res[2] = 2;
		return res;
	}

	// boolean max;
	public String toString() {
		return this.base.name + ":" + this.toComp.name;
	}

	// static double thresh = 1e9; //anything further away from the original
	// probe locations is ignored
	final File outpDir;
	// List<Entr> l = new ArrayList<Entr>();
	// List<Entr> l_max = new ArrayList<Entr>();
	//String probeset;
	String source1, source2;
	// List<Integer> cn;
	// String cn_string;
	// //DataCollection[] ldl;
	DataCollection base;
	public DataCollection toComp;
	// String[] result_set_names;
	CompoundEmissionStateSpace emstsp;
	// SortedSet[] originalProbes;
	// final Integer[] indices;
	// Color[] cols;

	
static class Point {
		
		List<Boolean> toCount = new ArrayList<Boolean>();
		private 
		List<Double> l = new ArrayList<Double>();
		List<Integer> positions = new ArrayList<Integer>();
		

		int size() {
			return countRegs ?  positions.get(positions.size()-1)-positions.get(0) : sizeTrue();//return l.size();//last - first;
		}

		public void clear() {
			positions.clear();
			l.clear();
			toCount.clear();
		}

		public void add(double d, boolean canAdd, int pos) {
			this.l.add(d);
			this.positions.add(pos);
			this.toCount.add(canAdd);
			
		}

		public int sizeTrue() {
			int cnt =0;
			for(int i=0; i<this.toCount.size(); i++){
				if(toCount.get(i)) cnt++;
			}
			return cnt;
		}
	}
	
	HWECalc[] hwe;
	int norm_index;
	private String chr;

	public ROC(DataCollection base, File outp, String chr, int centromere, int minnum, int minLength) {
		this.lightweight = true;
		this.chr = chr;
		this.base =  base;
		this.minNumberTrue = minnum;
		this.minLengthTrue = minLength;
		if(base instanceof LightWeightDataCollection){
			midB =((LightWeightDataCollection)base).midpoint;
		}
		else{
			midB = ((LightWeightDataCollection)((LightWeightMergedDataCollection)base).ldl[0]).midpoint;
		}
	/*	if(!countRegs) minNumberTrue = 0;
		else if(base instanceof LightWeightDataCollection2) minNumberTrue = 0;
		else if(base instanceof LightWeightMergedDataCollection  
				&& ((LightWeightMergedDataCollection)base).ldl[0] instanceof LightWeightDataCollection2)
			minNumberTrue=0;
		else minNumberTrue=1;*/
	//	this.minNumberTrue = (!countRegs ||  base instanceof LightWeightDataCollection2 )? 0 : 1;
	this.centromere = centromere;
System.err.println("SET MIN NUMB TO "+minNumberTrue);
		// this.strokes = strokes.toArray(new Stroke[0]);
		// this.indices = indices_to_use_as_gold_standard.toArray(new
		// Integer[0]);
		this.outpDir = outp;
		// cols = result_set_color.toArray(new Color[0]);
		// this.originalProbes = originalProbes;
		// this.ldl = ldl;
		emstsp = Emiss.getSpaceForNoCopies(Constants
				.backgroundCount(0));
		norm_index = emstsp.getGenoForCopyNo(2)[0];
		this.hwe = new HWECalc[] { new HWECalc(emstsp), new HWECalc(emstsp) };
		indiv = new ArrayList<String>();
		for (Iterator<String> it = base.dataL.keySet().iterator(); it
				.hasNext();) {
			String key = it.next();
			// for(int i=0; i<hes.length; i++){
			if (
					 base.dataL.get(key) != null) {
				// hes[0].add((HaplotypeEmissionState) base.dataL.get(key));
				// hes[1].add((HaplotypeEmissionState)
				// toComp.dataL.get(key));
				indiv.add(key);
			}
		}
		// this.result_set_names = result_set_names;
		/*
		 * Set<String> indiv = new HashSet<String>(base.dataL.keySet()); for(int
		 * k=0; k<ldl.length; k++){ indiv.retainAll(ldl[k].dataL.keySet()); }
		 * for(int k=0; k<ldl.length; k++){ ldl[k].restricToAlias(indiv); }
		 */

	}

	/*
	 * public double area(boolean max){ List<Entr> l1 = max ? this.l_max :
	 * this.l; int cnt =0; int cnt_wrong=0; int area =0; int tot =0; for(int
	 * i=0; i<l1.size(); i++){ if(l1.get(i).prob[index]>this.prob_upper){ cnt++;
	 * area +=cnt; tot++; } else if(l1.get(i).prob[index] < this.prob_lower){
	 * cnt_wrong++; area +=cnt; tot++;
	 * 
	 * }
	 * 
	 * } return (double) area/(double) (tot*tot); }
	 * 
	 * public double OTT(boolean max, int errors){ List<Entr> l1 = max ?
	 * this.l_max : this.l; int cnt =0; int cnt_wrong=0; for(int i=0;
	 * i<l1.size(); i++){ if(l1.get(i).prob[index]>this.prob_upper){ cnt++; }
	 * else if(l1.get(i).prob[index] < this.prob_lower){ cnt_wrong++;
	 * if(cnt_wrong > errors) return cnt; } } return cnt; }
	 */
	final boolean lightweight;
	static Comparator<Double> c1 = new Comparator<Double>() {

		public int compare(Double o1, Double o2) {
			return -1 * o1.compareTo(o2);
		}

	};
public  static double topExtrax = 0.02;
public  static double bottomExtrax = 0.02;

	

	
	
	
	
	
	
	static void append(Map<Double, Integer> m1, Map<Double, Integer> m2){
		for(Iterator<Double> it =m2.keySet().iterator(); it.hasNext();){
			Double d = it.next();
			Integer v = m1.get(d);
			Integer v1 = m2.get(d);
			m1.put(d, v==null ? v1 : v+v1);
		}
		m2.clear();
	}
	
	static void append1(SortedMap<Double,int[]> m1, SortedMap<Double, int[]> m2){
		for(Iterator<Double> it =m2.keySet().iterator(); it.hasNext();){
			Double d = it.next();
			int[] v1 = m1.get(d);
			int[] v2 = m2.get(d);
			if(v1==null){
		     	m1.put(d, v2);
			}
			else{
				for(int k=0; k<v1.length; k++){
					v1[k]+=v2[k];
				}
			}
		}
	//	m2.clear();
	}
	

	final public List<String> indiv;
	public class Inner {
		
		
		
		class Mapper {
			public final  String name;

			//File outF;
	final String type;
	
	Map<String, SortedMap<Integer, Integer>> regions ;
	
	
	//private boolean added = false;
	public void add(Mapper m){
		//if(added ) throw new RuntimeException("!1");
		//else added = true;
		this.allvals.addAll(m.allvals);
		m.allvals.clear();
		append(l_False, m.l_False);
		this.totNeg+=m.totNeg;
		this.totPos += m.totPos;
		m.totNeg=0;
		m.totPos =0;
		for(int i=0; i<totalTrue.length; i++){
			totalTrue[i]+=m.totalTrue[i];
			m.totalTrue[i] = 0;
		}
		for(Iterator<Integer> it = m.l_True.keySet().iterator(); it.hasNext();){
			Integer v = it.next();
			SortedMap<Double, Integer> m1 = this.l_True.get(v);
			if(m1==null){
				this.l_True.put(v, m1 = new TreeMap<Double, Integer>());
				append(m1,  m.l_True.get(v));
			}
			else{
				append(m1,  m.l_True.get(v));
			
			}
		}
		for(Iterator<Integer> it = m.lengthCount.keySet().iterator(); it.hasNext();){
			Integer v = it.next();
			Integer m1 = this.lengthCount.get(v);
			if(m1==null){
				lengthCount.put(v, m.lengthCount.get(v));
			}
			else{
				Integer lc = lengthCount.get(v);
				Integer lc1 =  m.lengthCount.get(v);
				this.lengthCount.put(v,(lc1==null ? 0 :lc1.intValue())+ (lc==null ? 0 : lc.intValue()));
				
			}
		}
		if(regions!=null){
			this.regions.putAll(m.regions);
		}
		m.lengthCount.clear();
	}


			public Mapper(PrintWriter pw, String type, String name) throws Exception {
				if(type.equals("max") || true){
					corr = new ArrayList<Double>();
					regions = new HashMap();
					regions.put(base.chrom,new TreeMap<Integer, Integer>());
				}
				this.name = name;
				this.type = type;
				this.pw = pw;
				//pw = new PrintWriter(new BufferedWriter(new FileWriter(out_all)));
			//	if (printDetails)
			//		pw.println(Format.sprintf(form_h, header));
				// pw_max= new PrintWriter(new BufferedWriter(new
				// FileWriter(out_max)));
			}

			final PrintWriter pw;
			SortedMap<Integer, SortedMap<Double, Integer>> l_True = new TreeMap<Integer, SortedMap<Double, Integer>>();
			
			
			SortedMap<Double, Integer> l_False = new TreeMap<Double, Integer>(c1);
			// SortedMap<Double, Integer> l_True1 = new TreeMap<Double,
			// Integer>(c1);
			Map<String, Point> buffer = new HashMap<String, Point>();
			SortedSet<Double> allvals = new TreeSet<Double>(c1);
			SortedMap<Integer, Integer> lengthCount = new TreeMap<Integer, Integer>();

			SortedMap<Integer, Integer> cumulativeByLength = new TreeMap<Integer, Integer>();

			boolean cumulative = false;

			double[] fpThresh = new double[fpr.length];
			double[] fpValue = new double[fpr.length];

			SortedMap<Integer, Integer> makeCumulativeByLength(int k) {
				SortedMap<Integer, Integer> res = cumulativeByLength;
				res.clear();
				int cum = 0;
				for (Iterator<Integer> it = l_True.keySet().iterator(); it
						.hasNext();) {
					Integer key = it.next();
					SortedMap<Double, Integer> m = l_True.get(key);
					SortedMap<Double, Integer> tail = m.headMap(fpThresh[k]);
					cum += tail.size() == 0 ? 0 : tail.get(tail.lastKey());
					res.put(key, cum);
				}
				return res;
			}
			
			

			private void makeCumulative(SortedMap<Double, Integer> m, double[] fpr,
					double[] fpThresh,  double tot) {
				int cum = 0;
				int fpr_index = 0;
				int cum1=0;
				for (Iterator<Entry<Double, Integer>> it = m.entrySet().iterator(); it
						.hasNext();) {
					Entry<Double, Integer> v = it.next();
					cum += v.getValue();
					v.setValue(cum);
					
					// m.put(v, cum);
					if (fpr != null && fpr_index < fpr.length) {
						for (int i = fpr_index; i < fpr.length; i++) {
							if (cum > fpr[fpr_index] * tot) {
								fpThresh[fpr_index] = v.getKey();
								fpValue[fpr_index] = (cum1)/tot;
								fpr_index++;
							}
						}
					}
					cum1=cum;
				}

			}

			private void makeCumulative(SortedMap<Integer, Integer> m) {
				int cum = 0;
				for (Iterator<Entry<Integer, Integer>> it = m.entrySet().iterator(); it
						.hasNext();) {
					Entry<Integer, Integer> v = it.next();
					cum += v.getValue();
					v.setValue(cum);
				}

			}

			void makeCumulative() {
				if (!cumulative) {
					makeCumulative(l_False, fpr, fpThresh,  totNeg);
					makeCumulative(this.lengthCount);
					for (Iterator<SortedMap<Double, Integer>> it = l_True.values()
							.iterator(); it.hasNext();) {
						makeCumulative(it.next(), null,  null, 1);
					}
					cumulative = true;
				}
			}

			
			
			public Integer[] getLengthBands(int number) {
				return new Integer[] {1000,l_True.lastKey()+1};
			/*	if (number == 1)
					return new Integer[] { 0 };
				// if(!cumulative) throw new RuntimeException("!!");
				double max = (double) totalTrue / (double) number;
				Integer[] spl = new Integer[number];
				double[] valSpl = new double[number];
				for (int i = 0; i < valSpl.length; i++) {
					valSpl[i] = i * max;
				}

				double cum = 0;
				int index = 0;
				int prev_val = 0;
				for (Iterator<Entry<Integer, Integer>> it = lengthCount.entrySet()
						.iterator(); it.hasNext();) {
					if (cum >= valSpl[index]) {
						spl[index] = prev_val;
						index++;
						if (index == spl.length)
							break;
					}
					Entry<Integer, Integer> nxt = it.next();
					prev_val = nxt.getKey();
					cum += nxt.getValue().doubleValue();
				}

				return spl;*/
			}

			public double getNumRight(Double value, int maxLength) {
			
				double totPos = 0;
				for (Iterator<Integer> it = l_True.keySet().iterator();
				
					
						 it.hasNext();) {
					Integer key = it.next();
					if(key<maxLength){
					SortedMap<Double, Integer> m = l_True.get(key);
					
					Integer add = m.get(value);
					if (add != null) {
						totPos += add.doubleValue();
					}
					}
				}
				return totPos;

			}

			/*
			 * if(prob[0]>=prob_upper){
			 * 
			 * 
			 * } else if(prob[0]<prob_lower){ //transferBuffer(prob_buffer[k], l);
			 * 
			 * 
			 * }
			 */
			int[] totalTrue = new int[len_bands.length];
			int totNeg = 0;
			int totPos = 0;
			Object[] toP = new Object[5];
			public List<Double> corr;

			void putFalse(double[] d,String id,  int pos, int posl, int posr){
				if(canAdd(pos,posl, posr)){
					Integer cnt = l_False.get(d[1]);
					l_False.put(d[1], cnt == null ? 1 : cnt + 1);
					this.totNeg++;
					
				}
			}
			
		
			
			
			
			void put(double[] d, String id, int pos, int posl, int posr, boolean last, String rsid) {
				allvals.add(d[1]);
				if(printDetails){pw.print(String.format("%5.3g",d[1]).trim()+"\t");
				if(last){
					pw.println(d[0]+"\t"+(d[0]>prob_lower)+"\t"+String.format("%5.3g",d[0])+"\t"+chr.split("\\.")[0]+":"+pos+"\t"+id+"\t"+rsid);
					pw.flush();
				}
				}
				if (d[0] < prob_lower) {
					int len = clearBuffer(id);
					
					putFalse(d, id, pos, posl, posr);
				} else if (d[0] >= prob_upper) {

					Point p = buffer.get(id);

					if (p == null) {
						buffer.put(id, p = new Point());
					}
					/*else if (pos < p.last) {
						clearBuffer(id);
					}*/
				     
					if(p.positions.size()>0 && p.positions.get(0) < centromere && pos > centromere){
						clearBuffer(id);
						buffer.put(id, p=new Point());
					}
					
					
					p.add(d[1],canAdd(pos, posl, posr), pos);
					
				}
			}
			


			protected boolean canAdd(int pos, int posl, int posr) {
					if(thresh1> 1e6) return true;
				if(posl<0 || posr<0) return false;
				int left = pos-posl;
				int right = posr - pos;
				if(left > thresh1 && right > thresh1){
					return false;
				}
			
				return true;
			}

			void finish() {
				if (buffer != null) {
					for (Iterator<String> it = buffer.keySet().iterator(); it
							.hasNext();) {
						String id = it.next();
						int len = clearBuffer(id);
					
					}
					this.buffer = null;
				}
				//this.pw.close();
			}

			protected int  clearBuffer(String id) {
				// for(Iterator<Point> it =
				// buffer.values().iterator();it.hasNext();){
				Point p = buffer.get(id);
				// int sze = p.l.size();
				if (p != null && p.l.size() > 0) {
					int len = p.size();
					if(countRegs && len>6000000){
						System.err.println(base.name+" long CNV "+id+"\t"+p.positions);
					}
				//	System.err.println("len is "+len);
					if(len>min_len){
					SortedMap<Double, Integer> m1 = l_True.get(len);
					if (m1 == null)
						l_True.put(len, m1 = new TreeMap<Double, Integer>(c1));
					
					for (int i = 0; i < p.l.size(); i++) {
						Double d = p.l.get(i);
						Boolean b = p.toCount.get(i);
						if(b){
						Integer cnt = m1.get(d);
						m1.put(d, cnt == null ? 1 : cnt + 1);

						this.totPos++;
					
						}
					}
					}
					int incr = countRegs ? 1: p.sizeTrue();
					
					Integer l_c = lengthCount.get(len);
					if( p.l.size()>minNumberTrue){
						double len1 = p.positions.get(p.positions.size()-1) - p.positions.get(0);
						if(len1>=minLengthTrue){
							   lengthCount.put(len, l_c == null ? incr : l_c + incr);// value)
							   if(this.regions!=null){
								   regions.get(base.chrom).
								   	put(p.positions.get(0), p.positions.get(p.positions.size()-1));
							   }
						}
			        }
				
					int len_ind=-1;
					for(int i=len_bands.length-1; i>=0; i--){
						if(len<len_bands[i] && len>=min_len){
							totalTrue[i] += incr;
							len_ind=i;
						}
					}
					if(len_ind>=0 && this.corr!=null){
						for(Iterator<Double> it = corr.iterator(); it.hasNext();){
							Double res = it.next();
							
							Integer cnt = Inner.this.corr[len_ind].get(res);
							if(cnt==null){
								Inner.this.corr[len_ind].put(res, 1);
							}
							else{
								Inner.this.corr[len_ind].put(res, cnt+1);
							}
						}
				//	}
					}
					p.clear();
					return len;
				}
				return 0;
				// }
			}

			public Iterator<Double> values() {
				return allvals.iterator();
			}

			public double getNumFalse(Double value) {
				if (cumulative)
					throw new RuntimeException("!!");
				if (this.buffer != null)
					throw new RuntimeException("!!");
				Integer res = l_False.get(value);
				if (res == null)
					return 0.0;
				else
					return res.doubleValue();
			}

			public void check(Mapper l) {
				l.finish();
				this.finish();
				/*if (totNeg != l.totNeg)
					throw new RuntimeException("!! " + totNeg + " " + l.totNeg);
				if (totPos != l.totPos)
					throw new RuntimeException("!! " + totPos + " " + l.totPos);
	*/
			}


			

			/*
			 * public Set<Entry<Double, SortedMap<Double, Integer>>> entrySet() {
			 * return l.entrySet(); }
			 */
		}
		class MapperBreak extends Mapper{

			public MapperBreak(PrintWriter pw, String type, String name) throws Exception {
				super(pw, type, name);
				// TODO Auto-generated constructor stub
			}
			
			void put(double[] d, double actual, String id, int pos, int posl, int posr, boolean last, String rsid) {
					Point p = null;
					allvals.add(d[1]);
					if(printDetails){
					pw.print(String.format("%5.3g",d[1]).trim()+"\t");
					if(last){
						pw.println(d[0]+"\t"+(d[0]>prob_lower)+"\t"+String.format("%5.3g",d[0])+"\t"+chr.split("\\.")[0]+":"+pos+"\t"+id+"\t"+rsid);
					//	pw.println(d[0]+"\t"+(d[0]>prob_lower)+"\t"+id+"\t"+chr+":"+pos);
						pw.flush();
					}
					}
						if(actual<prob_lower){
						clearBuffer(id);
					}
					else if(actual>prob_upper){
						 p = buffer.get(id);

						if (p == null) {
							buffer.put(id, p = new Point());
						} 
						/*else if (pos < p.last) {
							clearBuffer(id);
						}*/
						
					}
					if (d[0] < prob_lower) {
						putFalse(d, id, pos, posl, posr);
					} else if (d[0] >= prob_upper && actual > prob_upper) {
						p.add(d[1], canAdd(pos,posl, posr), pos);
					}
				}
			
		}
		
		public void add(Inner i){
			if(!this.cn.equals(i.cn)) throw new RuntimeException("!!");
			this.total+=i.total;
			for(int k=0; k<this.confMatrix.length; k++){
				for(int k1=0; k1<this.confMatrix[k].length; k1++){
					confMatrix[k][k1]+=i.confMatrix[k][k1];
				}
			}
			this.cnt_all+=i.cnt_all;
			for(int j=0; j<l_left.length; j++){
				l_left[j].add(i.l_left[j]);
				if(l_right[j]!=null) l_right[j].add(i.l_right[j]);
				l_indiv[j].add(i.l_indiv[j]);
			}
			l_max.add(i.l_max);
			for(int k=0; k<corr.length; k++){
			append(corr[k], i.corr[k]);
			}
			for(int k=0;k<conc.length; k++){
				if(conc[k]!=null){
				append1(conc[k], i.conc[k]);
				}
			}
		}
		
		public Integer[] cn;
		public int total = 0;

	
		
		public String toString() {
			return ROC.this.toString() + "_" + Arrays.asList(cn);
		}

		// public PrintWriter aucpw;
		
		boolean hwe_correct;
		// Set<String> set;
		//StringBuffer sb;
		MapperBreak[] l_left;
		MapperBreak[] l_right;
		Mapper l_max;
		Mapper[] l_indiv;
	//	final Mapper[] l_cn;
		// public double totNeg, totNegMax, totPos, totPosMax;
		List<Integer>[] geno;
	
		PseudoDistribution[] dists = new PseudoDistribution[2];
		double[] sig = new double[hwe.length];
		double[] sig_min = new double[hwe.length];
		double[][] sig1 = new double[hwe.length][];
		// int[] i = new int[2];
	//	String cn_string;
		//File out_conf;
		final public String name;

	
	//	final PrintWriter pw_all, pw_break, pw_max, pw_ld, pw_conf;
		
		
		
		SortedMap<Double, int[]>[]conc;
		
	//	final boolean[] iscum;
	//int[] totalByCnBase;
		
		
		
		
	
		
		Inner(String name, Integer[] cn, boolean hwe_correct,
				PrintWriter[] pw_all, PrintWriter[] pw_breakL, PrintWriter[] pw_breakR, PrintWriter pw_max, PrintWriter pw_ld, PrintWriter pw_conf
		) {
			this.cn = cn;
			this.name = name;
		this.probs = new double[LightWeightDataCollection.last];
	this.corr = new SortedMap[len_bands.length];
	for(int i=0; i<corr.length; i++){
		corr[i] = new TreeMap<Double, Integer>(comp_back);
	}
			conc = new SortedMap[Emiss.getSpaceForNoCopies(Constants.backgroundCount(0)).copyNumber.size()];
			for(int i=0; i<conc.length; i++){
				conc[i] = new TreeMap<Double, int[]>(comp_back);
			}
		/*	this.iscum = new boolean[conc.length];
			Arrays.fill(iscum, false);*/
			this.pw_conf = pw_conf;
			this.pw_ld = pw_ld;
		/*	this.corr  = new SortedMap[len_bands.length];
			for(int k=0; k<corr.length; k++){
				corr[k] = new TreeMap<Double, Integer>(comp_back);
			}*/
			//this.name = toComp.name;
			this.hwe_correct = hwe_correct;
			geno = new List[cn.length];
			for (int i = 0; i < cn.length; i++) {
				
				int[] res = emstsp.getGenoForCopyNo(cn[i]);
				geno[i] = new ArrayList<Integer>();
				for (int k = 0; k < res.length; k++) {
					geno[i].add(res[k]);
				}
			}
			
			
			try {
			//	int sze = indiv.size();
				l_left = new MapperBreak[indiv.size()];
				l_right = new MapperBreak[indiv.size()];
				l_indiv = new Mapper[indiv.size()];
				for(int i=0; i<indiv.size(); i++){
					l_left[i] = new MapperBreak(pw_breakL[i], "left break", indiv.get(i));
					l_right[i] = new MapperBreak(pw_breakR[i], "right break", indiv.get(i));
					l_indiv[i] = new Mapper(pw_all[i], "indiv", indiv.get(i));
				}
				
				l_max = new Mapper(pw_max, "max","max");
		
				
			} catch (Exception exc) {
				exc.printStackTrace();
			}
			/*List<Mapper> l = new ArrayList<Mapper>();
			l.add(l_max);
			for(int i=0; i<l_left.length; i++){
				l.add(l_left[i]); l.add(l_right[i])
			}*/
			this.mapper = new Mapper[] {l_max, l_left[0], l_indiv[0]};
			if(mapper.length!=mapLength) throw new RuntimeException("!!");
			this.prob = new double[indiv.size()][2];
			this.prob_prev_indiv = new double[indiv.size()][2];
			this.max_cnv = new int[2];
			this.max_index= new int[2];
			this.prob_max_ind = new double[2];
			this.prob_max_cnv = new double[2];
			this.cnv_ind = new boolean[emstsp.copyNumber.size()];
			Arrays.fill(cnv_ind,false);
			for(int i=0;i< cn.length;i++){
				cnv_ind[cn[i]]=true;
			}
		//	this.prob_cn = new double[this.indiv.size()][cn.length][2];
			if (hwe_thresh > 0) {
				for (int i = 0; i < hwe.length; i++) {
					hwe[i].setList(indiv, i == 0 ? base : toComp);
					sig1[i] = new double[5];
				}
			}
			this.probByCN = new double[emstsp.copyNumber.size()];
			this.confMatrix = new int[Math.min(probByCN.length,5)+1][Math.min(probByCN.length,5)+1];
		start = new int[indiv.size()];
			end = new int[indiv.size()];
		}

		Object[] toP = new Object[4];
final public PrintWriter pw_conf, pw_ld;
		public void finish() {
			for(int k=0; k<l_indiv.length; k++){
				this.l_indiv[k].finish();
			this.l_left[k].finish();
			this.l_right[k].finish();
			this.l_left[k].add(l_right[k]);
			if(k>0){
				l_indiv[0].add(l_indiv[k]);
				l_left[0].add(l_left[k]);
			}
			}
			this.l_max.finish();
			
			//this.printConf(this.pw_conf);
			/*int cum=0;
			for(Iterator<Entry<Double, Integer>> it = this.corr.entrySet().iterator(); it.hasNext();){
				Entry<Double, Integer> nxt = it.next();
				double v = nxt.getValue();
				cum+=v;
				nxt.setValue(cum);
				
			}*/
		}
		final double[][]prob;
		//final double[][][]prob_cn;

		final double[][] prob_prev_indiv;
		
		double[] probByCN;
		
	 private void getProbByCN(PseudoDistribution dist){
		Arrays.fill(probByCN, 0.0);
		Integer fixed = dist.fixedInteger();
		if(fixed!=null){
			probByCN[emstsp.getCN(fixed.intValue())]+=1.0;
		}
		else{
		 double[] pr = dist.probs();
		 for(int k=0; k<pr.length; k++){
			probByCN[emstsp.getCN(k)]+=pr[k];
		 }}
	 }
		
		public void transfer(double[][] prob, double[][] prob1){
			for(int i=0; i<prob.length; i++){
				System.arraycopy(prob[i], 0, prob1[i], 0, prob[i].length);
			}
		}
	final 	public SortedMap<Double, Integer>[]corr ;
		
		int[] max_cnv, max_index;
		double[] prob_max_cnv, prob_max_ind;
		final boolean[] cnv_ind;
		int[][]confMatrix;
		public void printConf(PrintWriter out_conf){
			try{
				String[] col_head = new String[confMatrix.length];
				String[] row_head = new String[confMatrix.length];
				for(int i=0; i<col_head.length; i++){
					col_head[i] = i+"";
					row_head[i] = i+"";//+name;
				}
				col_head[col_head.length-1] = "NA";
				row_head[row_head.length-1] = "NA";
			//	out_conf.println(name.replaceAll("_", "\\_"));
				PrintLatexTable.print(out_conf,name.replaceAll("_", "\\_"),col_head, row_head, this.confMatrix);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		
		double[] prob_max = new double[2];
		final int[] start, end;
		//private Component cnvs;
	//	final int[] cnv_alias;
		int cnt_all =0;
	final	double[] probs;
		public boolean run(String id, int pos, int posl, int posr, boolean last) {
			Arrays.fill(prob_max, 0.0);
			;
		    Arrays.fill(start, 0);
		    Arrays.fill(end, 0);
			
		
			outer1:  for (int k = 0; k < indiv.size(); k++) {
				try{
			    cnt_all++;
				String indiv = ROC.this.indiv.get(k);
				Arrays.fill(prob[k], 0.0);;
				HaplotypeEmissionState baseState = ((HaplotypeEmissionState) base.dataL
						.get(indiv));
				HaplotypeEmissionState compState = ((HaplotypeEmissionState) toComp.dataL.get(indiv
				));
				
			
				if(compState==null) throw new RuntimeException("!! missing "+indiv +" in "+this.name);
				double minp = 1.0;
				for (int ik = 0; ik < 2; ik++) {

					dists[ik] = ik == 0 ? process(baseState.emissions,midB, emstsp)
							: process(compState.emissions,midC, emstsp);
					/*if(ik==1){
						compState.emissions[midC] = dists[ik];
					}
					else*/
					if(!CompareBasic.combine){
					if(ik==0 && dists[ik]!=null){
						if(dists[ik].fixedInteger()==null){
							double[] probs = dists[ik].probs();
							if(probs[Constants.getMax(probs)]<0.95){
							dists[ik]=null;
							}
						}
					}
					}
					if(ik==1){
					   compState.emissions[midC] = dists[ik];
				    }
					else{
						baseState.emissions[midC] = dists[ik];
					}
					if(dists[ik]==null ){
						if(ik==0 && !CompareBasic.combine) continue outer1;
						dists[ik] = emstsp.getHWEDist1(0.0);
						if(ik==0) baseState.emissions[midB] = dists[ik];
						else compState.emissions[midC] = dists[ik];
					//	throw new RuntimeException("is null for "+indiv+" "+id+" "+(ik==0 ? base.name : toComp.name));
					}
					/*if(dists[ik].sum()<0.99) {
						double sum = dists[ik].sum();
						throw new RuntimeException("!!");
					}*/
				   this.getProbByCN(dists[ik]);
				   this.max_cnv[ik] = Constants.getMax(this.probByCN);
				   Integer max_ind = dists[ik].fixedInteger();
				   this.max_index[ik] = max_ind==null ? Constants.getMax(dists[ik].probs()) : max_ind;
				   this.prob_max_ind[ik] = dists[ik].probs(this.max_index[ik]);
				   if(prob_max_ind[ik]<minp){
					   minp = prob_max_ind[ik];
				   }
				   this.prob_max_cnv[ik] = probByCN[max_cnv[ik]];
					for (int j1 = 0; j1 < cn.length; j1++) {
							for (int j = 0; j < geno[j1].size(); j++) {
								prob[k][ik] += dists[ik].probs(geno[j1].get(j));// score(geno[j],
							}
					}
					if(minNumberTrue>0 && ik==1 && prob[k][ik]>0 && !CompareBasic.combine){
						getProbs(probs, compState);
						double max_min=0;
						for(int ij=0; ij<midC; ij++){
							double min = probs[midC];
							for(int ij1=0; ij1<midC; ij1++){
								if(probs[ij+ij1]<0.5){
									min = 0;
								}
							}
							if(min>max_min) max_min =min;
						}
						prob[k][ik] = max_min;
						/*if(prob[k][ik]>0){
							throw new RuntimeException("!!");
						}*/
					}
					
					if (prob[k][ik] > prob_max[ik]){
						prob_max[ik] = prob[k][ik];
					}
				}
				int base_ind = prob_max_cnv[0]>0.5 ? Math.min(max_cnv[0],4)  : 5;
				int comp_ind = prob_max_cnv[1]>0.9 ?Math.min(max_cnv[1],4)  : 5;
				confMatrix[base_ind][comp_ind]++;
				/*if(base_ind<2){
					pw_conf.println("mismatch "+id+" "+pos+" "+indiv+" "+toComp.name+" "+base_ind+" "+comp_ind);
				}*/
				if(max_cnv[0]==max_cnv[1]  && this.cnv_ind[max_cnv[0]]){
				//	if(max_cnv[0]==1){
					//	System.err.println("h");
					//}
					int ind =max_cnv[0];	
					int[] conc_ = this.conc[ind].get(minp);
					if(conc_==null){
						this.conc[ind].put(minp,conc_ = new int[2]);
					}
					if(max_index[0]==max_index[1]){
						conc_[0]++;
					}
					else{
						conc_[1]++;
					}
				}
				}
			catch(Exception exc){
				System.err.println(base.name+" "+toComp.name+" "+indiv.get(k)+" "+id);
				
				exc.printStackTrace();
			}
			}
		
			l_max.put(prob_max, "all", pos, posl, posr, last, id);
		
			if(prob_max[0]>prob_upper || true ){//&& prob_max[1]>0.99){
				for (int k = 0; k < indiv.size(); k++) {
					//if(prob[k][0]-prob[k][1]>0.01){
					//	System.err.println("h");
					//}
					l_indiv[k].put(prob[k], indiv.get(k), pos, posl, posr,last, id);
					l_left[k].put(diff(prob[k], prob_prev_indiv[k]),prob[k][0], indiv.get(k), pos, posl, posr, last, id);
					l_right[k].put(diff( prob_prev_indiv[k],prob[k]), prob_prev_indiv[k][0], indiv.get(k), pos, posl, posr, last, id);
				}
				if(prob_max[0]>prob_upper ){
				
						Double res = HLALD.calcLD(base, toComp, outpDir, avgForR2, midB, midC);
					//	double res1 = res;
					
					
						l_max.corr.add(res);
						if(printDetails){
						 pw_ld.print(String.format("%5.3g", res).trim()+"\t");
						 pw_ld.flush();
						 if(last)pw_ld.println(id+"\t"+chr+":"+pos);
						}
				}
			}
			transfer(prob, prob_prev_indiv);
			return prob_max[0]>prob_upper;
		}
private PseudoDistribution process(PseudoDistribution[] emissions, int midB, CompoundEmissionStateSpace emstsp) {
  if(emissions.length==1) return emissions[midB];	
  if(emissions[midB]==null) return null;
	double[] d = new double[emstsp.genoListSize()];
		Arrays.fill(d,0.0);
		int[] cn = new int[emstsp.copyNumber.size()];
		double[] max = new double[cn.length];
		Arrays.fill(max,0.0);
		for(int i=0; i<cn.length; i++){
			int[] cni = emstsp.getGenoForCopyNo(i);
		
			for(int k=0; k<emissions.length; k++){
				if(emissions[k]==null) continue;
				Integer fixed = emissions[k].fixedInteger();
				double val =0;
				if(fixed==null){
					double[] probs = emissions[k].probs();
					for(int j=0; j<cni.length; j++){
						int[]  hp = emstsp.getHaplopairForGeno(cni[j]);
						for(int j1=0; j1<hp.length; j1++){
							if(hp[j1]<probs.length){
								val = val+probs[hp[j1]];
							}
						}
					}
				}
				else if(emstsp.getCN(fixed.intValue())==i){
					val =  1.0;
				}
				
				if(val>max[i]){
					//if(k>0) {
					//	System.err.println("h");
					//}
					max[i] = val;
				}
			}
		}
		double rem = 1.0;
		for(int i=0; i<cn.length; i++){
			int[] cni = emstsp.getHaplopairForGeno(emstsp.getGenoForCopyNo(i)[0]);
			double val = Math.min(rem, max[i]);
			d[cni[0]] = val;
			rem = rem -val;
		}
		double sum = Constants.sum(d);
		if(sum<0.99){
			double sum1 = Constants.sum(max);
			throw new RuntimeException("!!");
		}
		return new SimpleExtendedDistribution(d,  Double.POSITIVE_INFINITY);
		}

private void getProbs(double[] probs2, HaplotypeEmissionState compState) {
	Arrays.fill(probs2,0);
		for(int i=0; i<probs2.length; i++){
			
			PseudoDistribution dist = compState.emissions[i];
			if(dist==null) continue;
			for (int j1 = 0; j1 < cn.length; j1++) {
				for (int j = 0; j < geno[j1].size(); j++) {
					probs2[i] += dist.probs(geno[j1].get(j));
				}
			}
			
		}
		
			
		}

private double[] tmp = new double[2];
		
/*left is going into cnv */
private double[] diff(double[] current, double[] prev) {
	
		tmp[0] = Math.max(current[0] -prev[0],0.0 );
		tmp[1] =  Math.max(current[1] -prev[1],0.0 );
	
		return tmp;
		}

		//File out_max, out_all, out_allR,out_allL ;
		private final  Mapper[] mapper;

		// PrintWriter pw_max, pw_all;
		public XYSeries[] getSeriesLength(Mapper list) throws Exception {
			
			list.finish();
			list.makeCumulative();
			XYSeries[] series1 = new XYSeries[fpr.length];
			for (int i = 0; i < fpr.length; i++) {
				list.makeCumulativeByLength(i);
				XYSeries series = new XYSeries(toComp.name+"_"+list.fpValue[i]);
				series1[i] = series;
				for (Iterator<Entry<Integer, Integer>> it = list.cumulativeByLength
						.entrySet().iterator(); it.hasNext();) {
					Entry<Integer, Integer> nxt = it.next();
					series.add(nxt.getKey(), nxt.getValue());
				}

			}
			return series1;
		}
		
		
		
		public XYSeries getLengthHistogram(Mapper list){
			return ROC.getLengthHistogram(list.lengthCount, base.name);
		}

		/** 0 = max, 1 = indiv, then cn+2 */
		public XYSeries[][] getSeries_(int type) throws Exception {
			
			Mapper list = mapper[type];
				
			XYSeries[][] ser1 = getSeries(list, type);
			XYSeries[] ser2 =
				getSeriesLength(list);
		    XYSeries[] ser3 = new XYSeries[] {this.getLengthHistogram(list)};
			XYSeries[][] res = new XYSeries[ser2.length + ext][];
			res[0] = ser1.length==1 ? new XYSeries[] {ser1[0][0]} : new XYSeries[] { ser1[0][0], ser1[1][0] };
			res[1] = 
				new XYSeries[] { ser1[0][1], ser1[0][2] };
			res[2] = 
				ser3;
			for (int i = 0; i < ser2.length; i++) {
				res[i + ext] = new XYSeries[] { ser2[i] };
			}
			return res;
		}

		public XYSeries[][] getSeries(Mapper list, int type) throws Exception {

			// pw.close();
			
			
			//PrintWriter pw = new PrintWriter(new File(list.outF
				//	.getAbsolutePath()
				//	+ "_ROC"));
			
			final double totalNeg = makeAsRate ? list.totNeg : 1;// list.totNeg;
			total = 0;

			list.finish();
			System.err.println("length bands are " + Arrays.asList(len_bands));
			XYSeries[][] series_ = new XYSeries[len_bands.length][3];

		//	boolean show = true;// !ROC.this.toComp.name.equals(ROC.this.base.name);
		//	if (!show || max)
		//		return series_;
			for (int i = 0; i < series_.length; i++) {
				final double totalPos = makeAsRate ? list.totalTrue[i] : 1;// list.totPos;
				XYSeries[] series = series_[i];
				for (int j = 0; j < series.length; j++) {
					series[j] = new XYSeries(toComp.name+"_"+(
							(i+1==len_bands.length) ? "all" : 
							Math.round((double)len_bands[i])));
				}
				series[2] = new XYSeries(toComp.name + "_histogram");

				// resort(1); //because 1 contains the comparison

				double totalPositive = 0;
				double totalNegative = 0;

				// Collections.sort(list, comp0);
				double prevWrong = 0;
				double prevRight = 0;
				Object[] toP = new Object[4];

				// series[0].add(totalNegative, totalPositive);
				double auc = 0;
				double totalPos0 = 0;
				double totalNeg0 = 0;
				int p_prev = 10;
				double x_prev = 0;
				double y_prev = 0;
				double frac_sum = 0;
				for (Iterator<Double> it = list.values(); it.hasNext();) {
					Double  value = it.next();
					double p1 = value;

					double numRight = list.getNumRight(value,len_bands[i]);
					double numWrong = list.getNumFalse(value);
				//	System.err.println(p1 + " " + numRight + " " + numWrong);
					if (prevWrong == 0 && numWrong == 0 && numRight > 0) {

					} else if (prevRight == 0 && numRight == 0 && numWrong > 0) {

					} else {
						double x = totalNegative / totalNeg;
						double y = totalPositive / totalPos;
						series[0].add(x, y);
						auc += y_prev * (x - x_prev)*(1.0/(double)list.totPos)*(1.0/(double)list.totNeg);
						auc += 0.5 * (y - y_prev) * (x - x_prev)*(1.0/(double)list.totPos)*(1.0/(double)list.totNeg);
						x_prev = x;
						y_prev = y;
					}
					total += numWrong;

					if ((int) Math.ceil(10 * p1) < p_prev) {

						double tot = ((totalPositive - totalPos0) + (totalNegative - totalNeg0));
						double totP = (totalPositive - totalPos0);
						double frac = totP / (double) list.totalTrue[i];
						frac_sum += frac;
						// System.err.println("summ "+((((double)p_prev)/10.0)
						// -0.05)+" "+tot+" "+totN+" "+
						// (totN/tot));
					/*	System.err.println(this.cn_string + " " + p1
								+ " p_prev " + p_prev + " current "
								+ (int) Math.ceil(10 * p1) + " " + frac + " "
								+ (1.0 - totP / tot));
						System.err.println("frac is " + frac);*/
						if (frac > 0.0) {
							series[1].add((((double) p_prev) / 10.0) - 0.05,
									totP / tot);
						}
						series[2].add((((double) p_prev) / 10.0) - 0.05, frac);
						p_prev = (int) Math.ceil(10 * p1);
						totalPos0 = totalPositive;
						totalNeg0 = totalNegative;
					}
					totalPositive += numRight;
					totalNegative += numWrong;

					prevWrong = numWrong;
					prevRight = numRight;
					/*
					 * // Double p0 = l1.next(); Integer valu = e.get(p0);
					 * if(p0> prob_upper){ // this is a positive
					 * 
					 * if(prevAcross>0){ series.add(totalNegative,
					 * totalPositive); prevAcross=0; } totalPositive+=valu;
					 * if(prevUp==0){ prevUp++; } } else if(p0 < prob_lower){ //
					 * this is a negative
					 * 
					 * if(prevUp>0){ series.add(totalNegative, totalPositive);
					 * prevUp=0; } totalNegative+=valu; if(prevAcross==0){
					 * prevAcross++; } }
					 */

				}
				// if(pw!=null) pw.close();
				double tot = ((totalPositive - totalPos0) + (totalNegative - totalNeg0));
				double totP = (totalPositive - totalPos0);
				double frac = totP / (double) list.totalTrue[i];
				frac_sum += frac;
			//	if (!countRegs && Math.abs(1.0 - frac_sum) > 0.001)
			//		throw new RuntimeException("frac does not add up "+frac_sum);
				if (Double.isInfinite(totalNegative)
						|| Double.isNaN(totalNegative))
					throw new RuntimeException("problem");
				// double remainingNeg = list.totNeg - totalNegative;
				// double remainingPos = list.totPos - totalPositive;
				// series[0].add(totalNegative/totalNeg,
				// totalPositive/totalPos);
				series[1].add(-0.05, (totP / tot));
				series[2].add(-0.05, frac);
				//System.err.println("final p_prev " + p_prev + " " + totalPos0
				//		+ " " + totalPositive + " " + list.totPos + " "
				//		+ totalNeg0 + " " + totalNegative + " " + list.totNeg);
				double x = totalNegative / totalNeg;
				double y = totalPositive / totalPos;
				series[0].add(x, y);
				auc += y_prev * (x - x_prev);
				auc += 0.5 * (y - y_prev) * (x - x_prev);
				x_prev = x;
				y_prev = y;
			//	pw.println(ROC.this.probeset + "\t" + ROC.this.base.name + "\t"
			//			+ toComp.name + "\t" + cn_string + "\t" +(this.mapper[type].type) + "\t"
			//			+ auc);
			//	pw.flush();
			}
			//pw.close();
			return series_;
		}
		

		public void check(Inner inner) {
			for(int k=0; k<l_left.length; k++){
			this.l_left[k].check(inner.l_left[k]);
			}
			this.l_max.check(inner.l_max);

		}

		public void run(String rs, int pos,  boolean last) {
		run(rs, pos, -1, -1,  last);
			
		}
		public String[] getTypes() {
			String[] res = new String[mapper.length];
			for(int i=0; i<res.length; i++){
				res[i] = mapper[i].type;
			}
			return res;
		}

		
		
		public XYSeries[] getConc() {
			
		
				SortedMap<Double, int[]>[] corr1 = conc;
				 XYSeries[] ser = new XYSeries[corr1.length];
				 for(int i=0; i<ser.length; i++){
					 ser[i] = new XYSeries(name+"_CN="+i);
					 if(corr1[i].size()>0){
						/* if(iscum[i]){
							 throw new RuntimeException("");
						 }*/
						 makeCumulative(corr1[i]);
						// iscum[i] = true;
						 int[] last = corr1[i].get(corr1[i].lastKey());
					  //   double tot1  = ROC.makeAsRate ? last[0]+last[1] : 1;
					
					   
						for(Iterator<Entry<Double, int[]>> it = corr1[i].entrySet().iterator(); it.hasNext();){
							Entry<Double, int[]> nxt= it.next();
							int[] v = nxt.getValue();
							double sum = (double)(v[0]+v[1]);
							if(sum>0){
								double y =makeAsRate ?  (double)v[1]/sum : v[1];
								double x = makeAsRate ? sum/(double)cnt_all : sum;//nxt.getKey().doubleValue();
							if(cn.length==1){
								x = makeAsRate ? 1.0 - x : cnt_all -x;
							}
								if (logconcy  && y<1e-8) y = 1e-5;
								if(x>0) {
									if(!makeAsRate || y < 0.01) ser[i].add(x, y);
								}
							}
							
						}
					 }
				 }
					return ser;
			
		}

		public XYSeries[] getCum() {
			
				SortedMap<Double, Integer>[] corr1 = corr;
				
				XYSeries[] ser = new XYSeries[corr1.length];
				for(int i=0; i<ser.length; i++){
					 if(corr1[i].size()>0){
						makeCum(corr1[i]);
				 double tot1  = ROC.makeAsRate ? corr1[i].get(corr1[i].lastKey()) : 1;
				ser[i] =  new XYSeries(name);
				 new XYSeries(name);
					for(Iterator<Entry<Double, Integer> >it = corr1[i].entrySet().iterator(); it.hasNext();){
						Entry<Double , Integer> d = it.next();
					//	if(Double.isNaN(d)) continue;
						if(d.getKey().doubleValue()>0){
						ser[i].add(d.getKey().doubleValue(),d.getValue().doubleValue()/tot1);
						}
					}
					}
				}
					return ser;
			
		}

		private boolean include(Integer integer) {
			/*
			 * inner: for(int i=0; i<this.originalProbes.length; i++){
			 * 
			 * SortedSet<Integer> ts = originalProbes[i].tailSet(integer);
			 * if(ts==null || ts.size()==0) continue inner; else if( (ts.first() -
			 * integer)<thresh) continue inner; else{ SortedSet<Integer> ss =
			 * originalProbes[i].headSet(integer+1); if(ss==null || ss.size()==0)
			 * continue inner; else if( (integer - ss.last())<thresh) continue
			 * inner; } return false; }
			 */
			return true;
		}

	}

	/* all, max */
	public Inner run(String name, Integer[] cn, boolean hwe_correct,
			PrintWriter[] pw_all, PrintWriter[] pw_breakL,PrintWriter[] pw_breakR, PrintWriter pw_max,
			PrintWriter pw_ld, PrintWriter pw_conf
	) throws Exception {
		return new Inner( name, cn, hwe_correct, pw_all, pw_breakL, pw_breakR,
				
				pw_max, pw_ld, pw_conf);

		/*
		 * inner: for(int ii=0; ii<base.loc.size(); ii++){ inner.run(ii); }
		 * 
		 * return new XYSeries[] {inner.getSeries(false),
		 * inner.getSeries(true)};
		 */
	}

	
	
	
	
	
	private static double sum(Collection<Integer> values) {
		double i=0;
		for(Iterator<Integer> it = values.iterator(); it.hasNext();){
			i+=it.next();
		}
		return i;
	}
	
	private static double sum1(Collection<int[]> values) {
		double i=0;
		for(Iterator<int[]> it = values.iterator(); it.hasNext();){
			int[] nxt = it.next();
			for(int k=0; k<nxt.length; k++){
				i+=nxt[k];
			}
		}
		return i;
	}
/*	public static XYSeries[] getQQ(Inner inner1, Inner inner2){
		SortedMap<Double, Integer>[] corr1 = inner1.corr;
		SortedMap<Double, Integer>[] corr2 = inner2.corr;
		XYSeries[] res = new XYSeries[corr1.length];
		for(int k=0; k<res.length; k++){
			 double tot1  = 1.0;//getValue(corr1, -1);
			double tot2  = 1.0;///getValue(corr2, -1);corr2.get(corr2.lastKey());
			if(tot1!=tot2) throw new RuntimeException("!!");
			SortedSet<Double> vals = new TreeSet<Double>(corr1[k].keySet());
			vals.addAll(corr2[k].keySet());
			XYSeries ser = new XYSeries(inner1.name+"_"+inner2.name);
			res[k] = ser;
			for(Iterator<Double> it = vals.iterator(); it.hasNext();){
				Double d = it.next();
				ser.add(getValue(corr1[k], d, true)/tot1, getValue(corr2[k], d, true)/tot2);
			}
		}
		return res;
	}*/
	
	/** get count greater or equal to 
	static double getValue(SortedMap<Double, Integer> corr1, Double d, boolean cum){
		if(cum){
		double cnt1 =0;
		Map<Double, Integer> tail = corr1.tailMap(d);
	for(Iterator<Integer> it = tail.values().iterator(); it.hasNext();){
	
			cnt1+=it.next();
		
	}
	
		return cnt1;
		}
		else return corr1.get(d);
	}
	*/
	static double makeCum(SortedMap<Double, Integer> corr1){
	
		int cnt1 =0;
		
			for(Iterator<Entry<Double, Integer>> it = corr1.entrySet().iterator(); it.hasNext();){
				Entry<Double, Integer> nxt = it.next();
					cnt1+=nxt.getValue();
					nxt.setValue(new Integer(cnt1));
					
				
			}
	
			return cnt1;
		}
	

	
	static void makeCumulative(SortedMap<Double, int[]> corr1){
		int[] cnt1 =new int[] {0,0};
		for(Iterator<Entry<Double, int[]>> it = corr1.entrySet().iterator(); it.hasNext();){
			Entry<Double, int[]> nxt = it.next();
			for(int k=0; k<cnt1.length; k++){
				cnt1[k]+=nxt.getValue()[k];
			}
			System.arraycopy(cnt1, 0, nxt.getValue(), 0, cnt1.length);
		}
	}

	public static void plot(double[] x, double[] y, boolean logx, boolean logy ){
		XYSeriesCollection xys = new XYSeriesCollection();
		XYSeries x1 = new XYSeries("x");
		xys.addSeries(x1);
		for(int i=0; i<x.length; i++){
			x1.add(x[i], y[i]);
			
		}
		final JFreeChart chart = ChartFactory.createXYLineChart("", "", "", xys, 	PlotOrientation.VERTICAL, false, false, false);
		if (logx)
			chart.getXYPlot().setDomainAxis(new LogarithmicAxis(""));
		if (logy)
			chart.getXYPlot().setRangeAxis(new LogarithmicAxis(""));
		ChartPanel cp = new ChartPanel(chart);
		JFrame frame = new JFrame();
		frame.getContentPane().add(cp);
		frame.setVisible(true);
		frame.pack();
	}
	
	
	/*
	 * public ChartPanel plot(boolean max,int index){
	 * 
	 * ChartPanel cp = getCollection(max,index); return cp;
	 * 
	 * }
	 */
	
	
	public static Range getRange(List<Double> l, boolean log, boolean excludeTop, double percentile, 
			double botMult, double topMult){
		Collections.sort(l);
		//double max = l.size()==0 ? max1 : max1*l.get(l.size()-1);
		if(excludeTop && l.size()>0){
			l = l.subList(0, l.indexOf(l.get(l.size()-1)));
		}
		double x_min = log ? 1e-5: 0;
		double x_max =x_min;
		if(l.size()>0){
			 x_min = l.get(0);
			 x_max = l.get((int)Math.ceil(percentile*(double) l.size())-1)+0.001;
		}
		if(x_max<x_min) x_max = x_min;
		
	
			if(!log){
				double diff = x_max - x_min;
				
				x_max = x_max + diff*topMult;
				 diff = x_max - x_min;
				x_min = x_min - diff*botMult;
			}
			else{
				//double ratio = x_max/x_min;
				x_max = x_max *(1+topMult);
				x_min = x_min*(1-botMult);
			}
		if(log && x_min>0){
			double l1 =Math.log10(x_min);
			double u1 = Math.log10(x_max);
			if(Math.ceil(u1)-Math.floor(l1)<=1){
				x_min = Math.floor(l1);
				x_max = Math.ceil(u1);
			}
		}
		return new Range(x_min,x_max);
	}
	public static Range getXMax(XYSeriesCollection[] xys,  boolean log, boolean excludeTop, double perc){
	/*	boolean[] constant= new boolean[xys.length];
		for(int i=0; i<xys.length; i++){
			constant[i] = allconstant(xys[i]);
		}*/
		
		List<Double> l = new ArrayList<Double>();
	//	List<Double> l1 = new ArrayList<Double>();
		//double x_min = 0.0;
		for(int i=0; i<xys.length; i++){
			
			for(int k=  0; k<xys[i].getSeriesCount(); k++){
			//	if(i==xys.length-1 &&)
				
				XYSeries er = xys[i].getSeries(k);
				for(int j=0; j<er.getItemCount(); j++){
				//	double y = er.getY(j).doubleValue();
					double x = er.getX(j).doubleValue();
				//	if(y<1 && x<1){
						//if(!constant[i])
							l.add(x);
					//	else l1.add(x);
					//}
					
				}
			}
		}
		
		return getRange(l, log, excludeTop, perc, ROC.bottomExtrax, ROC.topExtrax );
//				Math.min(max,x_max*1.01));
	}
	private static boolean allconstant(XYSeriesCollection seriesCollection) {
		// TODO Auto-generated method stub
		return false;
	}
	// public final Stroke[] strokes ;
	public static Range getYMax(XYSeriesCollection[] xys, boolean log, boolean excludeTop, double perc){
		
		List<Double> l = new ArrayList<Double>();
		for(int i=0; i<xys.length; i++){
			for(int k=0; k<xys[i].getSeriesCount(); k++){
			//	if(i==xys.length-1 &&)
				XYSeries er = xys[i].getSeries(k);
				for(int j=0; j<er.getItemCount(); j++){
					double y = er.getY(j).doubleValue();
					l.add(y);
					
				}
			}
		}
		return getRange(l, log, excludeTop, perc,0.05,ROC.topExtra);
	}
	public static double topExtra = 0.02;
	public static boolean nameCurves = true;
	public static double overlap=0.1;
	public static ChartPanel getChartPanel(XYSeriesCollection[] xys,
			List<Color> cols, List<Stroke> strokes, String probeset,
			String result_set_names, Integer[] cn, int max, int total,
			String chr, String xaxis1, String yaxis1, String title1,
			boolean logx, boolean logy, boolean errorGraph, String level, double[] x_max1, boolean dots) {
		int ind = chr.indexOf(".zip");
		errorGraph = title1.equals("Error rate ");
		/*
		 * for(int i=cols.size()-1; i>=0; i--){
		 * if(cols.get(i).equals(Color.white)){ xys[0].removeSeries(i);
		 * cols.remove(i); } }
		 */
		if (xys.length > 1 && xys[0].getSeriesCount() > 1 &&errorGraph) {
			xys[0].removeSeries(0);
			xys[1].removeSeries(0);
		}
		for (int i = 0; i < xys.length; i++) {
			adjustForZero(xys[i], logx, logy);
		}
		if (logx || logy) {
			System.err.println("correcting for log");
			for (int i = 0; i < xys.length; i++) {
				correctForLog(xys[i], logx, logy);
			}
		}
		
		String chr1;
		if (ind >= 0) {
			chr1 = chr.substring(0, ind);
		} else {
			chr1 = chr;
		}

		String type = cn[0] == 0 || cn[0] == 1 ? "deletions" : "amplifications";
		String xaxis = xaxis1.replaceAll("aberration", type).replaceAll("ions",
				"ion");
		String yaxis = yaxis1.replaceAll("aberration", type).replaceAll("ions",
				"ion");
		String title = yaxis1.replaceAll("aberration", type).replaceAll("ions",
				"ion");
	
		final JFreeChart chart = ChartFactory.createXYLineChart(
		// title+"Measuring probes on "+probeset
				// +"\n using as benchmark results  from "+result_set_names
				// +(!chr.equals("all") ? "\n on chromosome: "+chr1 : "")
				title1
				// +"\n identification of "+type+" at "+level+" level"
				// +(deletion ? "\n prediction of deletions  " :
				// "\n prediction of amplifications")+

				, xaxis, // domain axis label
				yaxis, // range axis label
				xys[0], // data
				PlotOrientation.VERTICAL, 
				legend,
				// !type.equals("deletions") && !level.equals("population"), //
				// include legend
				false, // tooltips?
				false // URL generator? Not required...
				);
	//	if(true) return new ChartPanel(chart);
		XYPlot plot = chart.getXYPlot();
		
		// this.subtitles.add(legend);
		// legend.addChangeListener(this);
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
		plot.setBackgroundPaint(Color.WHITE);
		if (logx){
			LogarithmicAxis loge = new LogarithmicAxis(xaxis);
			loge.setExpTickLabelsFlag(true);
			chart.getXYPlot().setDomainAxis(loge);
		}
		if (logy){
			LogarithmicAxis loge = new LogarithmicAxis(yaxis);
			loge.setExpTickLabelsFlag(true);
			chart.getXYPlot().setRangeAxis(loge);
		}
			
		// XYItemRenderer rend = chart.getXYPlot().getRenderer();
		final XYLineAndShapeRenderer[] renderer = new XYLineAndShapeRenderer[xys.length];
			
		for(int i=0; i<renderer.length; i++){
			renderer[i] = new XYLineAndShapeRenderer();
			renderer[i].setShapesFilled(false);
			renderer[i].setShapesVisible(dots);
			renderer[i].setLinesVisible(Boolean.TRUE);
			chart.getXYPlot().setRenderer(i, renderer[i]);
			if(i<strokes.size()) renderer[i].setBaseStroke(strokes.get(i));
			//else renderer[i].setSeriesStroke(0,new BasicStroke(0.8f,
				//	BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,
					//CompareBasic.generate(5), 1.0f));
			
		}
		
		if (xys.length > 1 && errorGraph) {
			final XYBarRenderer renderer1 = new XYBarRenderer(0.9);
			renderer1.setDrawBarOutline(true);
		//	renderer1.set
			for(int i=0; i<renderer.length; i++){
				renderer[0].setSeriesPaint(0, Color.black);
				renderer[0].setLinesVisible(false);
				renderer[0].setSeriesOutlinePaint(0, Color.black);
				renderer[0].setSeriesFillPaint(0, Color.black);
			}
			for(int i=0; i<cols.size(); i++){
		
			renderer1.setSeriesFillPaint(i, Color.LIGHT_GRAY);
			renderer1.setSeriesPaint(i, Color.LIGHT_GRAY);
			renderer1.setSeriesOutlinePaint(i, cols.get(i));
			}
			chart.getXYPlot().setRenderer(1, renderer1);
			ValueAxis axis = chart.getXYPlot().getDomainAxis();
			axis.setAutoRange(false);
			axis.setRange(new Range(-0.1, 1.1));
		}

		// renderer1.set

		AffineTransform at = new AffineTransform();
		at.setToScale(0.5, 0.5);
		Shape[] sh = DefaultDrawingSupplier.DEFAULT_SHAPE_SEQUENCE;
		   double size = 2.0*resolution;
	    double delta = size / 2.0;
		sh[0] = new Ellipse2D.Double(-delta, -delta, size, size);
		for (int i1 = 0; i1 <= 1; i1++) {
			for (int i = 0; i < cols.size(); i++) {

				for(int k=0; k<renderer.length; k++){
					boolean white =false;// k>0 && i==0;
				{
				
				renderer[k].setSeriesOutlinePaint(i,white ? Color.WHITE : cols.get(i));
				if(!errorGraph){
					renderer[k].setSeriesFillPaint(i, white ? Color.WHITE : cols.get(i));
				    renderer[k].setSeriesPaint(i,white ? Color.WHITE : cols.get(i));
				    if(strokes.size()>i)renderer[k].setSeriesStroke(i,strokes.get(i));
				}
				else{
					renderer[k].setSeriesFillPaint(i, Color.WHITE );
					
					//renderer[k].setSeriesPaint(i, Color.WHITE : cols.get(i));
				}}
				if(title.startsWith("Concord")){
					renderer[k].setSeriesShape(i, at.createTransformedShape(renderer[k].getSeriesShape(i)));
				}
				if(setShape)renderer[k].setSeriesShape(i, sh[k % sh.length]);
				}

			}
		}
		boolean error = false;
		
		for (int i = 0; i < xys.length; i++) {
			plot.setDataset(i, xys[i]);
			if(errorGraph){
			if (i == 1 ) {
				error = true;
				NumberAxis numA = logy ? new LogarithmicAxis(
						"Fraction of all true " + type
								+ " which are in probability bin")
						: new NumberAxis("Fraction of all true " + type
								+ " which are in probability bin");
				numA.setRange(0.0, 1.0);
				plot.setRangeAxis(i, numA);
				plot.mapDatasetToRangeAxis(i, i);
			} else {
				if (logy) {
					NumberAxis numA = (NumberAxis) plot.getRangeAxis();
					Range range = numA.getRange();
					numA.setAutoRange(false);
					numA.setRange(new Range(1, 100000));
				}
			}
			}
			// plot.mapDatasetToRangeAxis(i, i);

		}
		if (error) {
			plot.getRangeAxis(0).setAutoRange(false);
			plot.getRangeAxis(0).setRange(0.0, 1.0);
		}
		boolean excludeTop = !title.startsWith("Conc") && !title.startsWith("Cumulative");
		if(errorGraph) excludeTop = false;
		double perc = ROC.yperc;
		{
			if(!logy){
		           plot.getRangeAxis(0).setAutoRange(false);
		         plot.getRangeAxis(0).setRange( getYMax(xys, logy, excludeTop, perc));
			}
		{
			//plot.getDomainAxis(0).setAutoRange(false);
			plot.getDomainAxis(0).setRange(getXMax(xys,logx, logx ? false : excludeTop, logx ? 1.0 : perc), true, true);
		}
		/*if(logx){
			LogarithmicAxis dom = ((LogarithmicAxis)chart.getXYPlot().getDomainAxis());
			Range range =((LogarithmicAxis)chart.getXYPlot().getDomainAxis()).getRange(); 
			double diff = range.getUpperBound() - range.getLowerBound();
			double scale =Math.pow(10,Math.floor(Math.log10(diff/10.0)));
			dom.setTickLabelsVisible(true);
//			NumberTickUnit nu = dom.getTickUnit();
//			dom.setTickUnit(new NumberTickUnit(scale), true, true);
			//dom.setTickUnit(new NumberT)
			
		}*/
		}
		int ik=xys.length-1;
		
		//double[] max_x =x_max;
		/*double[]  minx = new double[xys[ik].getSeriesCount()];
		for(int k=0; k<xys[ik].getSeriesCount(); k++){
			double min_x = 1e5;
			XYSeries series = xys[ik].getSeries(k);
			for(int j=0; j<series.getItemCount(); j++){
				double x = series.getX(j).doubleValue();
				if(x>max){
					max_x = x;
					
				}
				if(x<min_x){
					min_x = x;
				}
			}
			minx[k] = min_x;
			
		}
		Arrays.sort(minx);
		if(minx.length>0){
		double min_x = minx.length==1 ? minx[0] : minx[1];
		if(max_x < min_x){
			max_x = min_x+1;
		}*/
		/*if(max_x!=null && !logx && false){
		plot.getDomainAxis().setAutoRange(false);
		plot.getDomainAxis().setRange(new Range(max_x[0], max_x[1]));
		}*/
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		chart.getTitle().setFont(font10);
		
		((XYPlot) chart.getPlot()).getDomainAxis().setLabelFont(font10);
		((XYPlot) chart.getPlot()).getDomainAxis().setTickLabelFont(font10);
		((XYPlot) chart.getPlot()).getRangeAxis().setLabelFont(font10);
		((XYPlot) chart.getPlot()).getRangeAxis().setTickLabelFont(font10);
		int width = (int) Math.floor(680*resolution);
		int height = (int) Math.floor(420*resolution);
	//	Constants.scatterWidth = 
		// axis.setRange(0, 0.1*axis.getRange().getUpperBound());
		ChartPanel cp = new ChartPanel(chart,
				width, height,width, height,width, height,
			//	Constants.scatterWidth(),Constants.scatterWidth(),
			//	Constants.scatterWidth(),Constants.scatterWidth(),
			//	Constants.scatterWidth(),Constants.scatterWidth(),
   	          false,  true,      true,  true,    true,   true);
		SortedMap<Double,Double>d = new TreeMap<Double, Double>();
		Range xrange = plot.getDomainAxis().getRange();
		double x_mid = logx ? Math.sqrt(xrange.getLowerBound()*xrange.getUpperBound())  : xrange.getCentralValue();

		
			//SortedMap<Double, int[]> order = getAscending(xys[0],x_mid);
		
		if(legend) chart.getLegend().setItemFont(font10);
		//chart.getLegend().s
		//[]msi =new MaxSepInd[xys.length];
		boolean isroc =ROC.nameCurves  && CompareBasic.combine;// title1.toLowerCase().indexOf("roc")>=0;
		if(xrange.getLength()>0 && title1.toLowerCase().indexOf("error")<0 && isroc){
			System.err.println("title "+title1);
			//List<String>nme =new ArrayList<String>(0);
			
			
		for(int kk=0; kk<xys.length;kk++){
			MaxSepInd  msi = new MaxSepInd((XYSeries[])xys[kk].getSeries().toArray(new XYSeries[0]),xrange,plot.getRangeAxis().getRange(), logx,logy
					,ROC.useSlope
			);
			msi.run();
		if(msi!=null){
			for(int k=0; k<msi.ser.length; k++){
				//Entry<Double, int[]> nxt = it.next();
				XYSeries series = xys[kk].getSeries(k);
					String name = series.getKey().toString();
					boolean above = msi.above(k);
					double[] xyz = msi.getBest(k);
				
					XYTextAnnotation annot = new XYTextAnnotation(name, xyz[0], xyz[1]);
					if(k<cols.size()) annot.setPaint(cols.get(k));
				//	boolean decr = false;
				//	if(series.getItemCount()>1)
				//	 decr = series.getY(0).doubleValue()> series.getY(series.getItemCount()-1).doubleValue();
					annot.setTextAnchor(above ? TextAnchor.BOTTOM_LEFT : TextAnchor.TOP_RIGHT);
					annot.setRotationAnchor(above ? TextAnchor.BOTTOM_LEFT : TextAnchor.TOP_RIGHT);
					if(series.getItemCount()>5    ){
						
					
					annot.setRotationAngle( -1*xyz[2] );
					}
					annot.setFont(font10);
						plot.addAnnotation(annot);
				}
			//	}
			//}
		}
		}
		}
		// chart.setLegendVisible(false);
		return cp;
	}



	 public static String process(String string1) {
		 if(string1.indexOf("_")<0 && string1.indexOf('-')<0) return string1;
		String string = string1.toLowerCase().replaceAll("train", "244k+1m");
		if(string.indexOf("penn")>=0) return "PennCNV(1M)";
		else if(string.indexOf("cnvp")>=0) return "cnvPartition(1M)";
		else if(string.indexOf("quanti")>=0  && string.indexOf("2160809")>=0) return "QuantiSNP(1M)(b)";
		else if(string.indexOf("quanti")>=0) return "QuantiSNP(1M)";
		else if(string.indexOf("agilent")>=0){
			if(string.indexOf("244")>=0) 	return "ADM(244k)";
			else if(string.indexOf("185k")>=0)
			return "ADM(185k)";
		}
		else if(string1.indexOf("hap")>=0){
			if(string.indexOf("317k")>=0) return "cnvHap+snph(317k)";
			else if(string.indexOf("1m")>=0)return "cnvHap+snph(1M)";
		}
		else{
			String pref = "cnvHap(";
			String mid = 
				(string.indexOf("244k")>=0 ? "244k+" :"") +
				(string.indexOf("1m")>=0 ? "1M+" :"") +
				(string.indexOf("185k")>=0 ? "185k+" :"") +
				(string.indexOf("317k")>=0 ? "317k+" :"");
			if(mid.length()>0)
			
			return pref + mid.substring(0,mid.length()-1)+")";
			
		}
		return string1;
	}

	private static SortedMap<Double, int[]> getAscending(XYSeriesCollection x, double xmid) {
		SortedMap<Double, int[]> m = new TreeMap<Double, int[]>();
	//	x.getS
		for(int i=0; i<x.getSeriesCount(); i++){
			double[] y = getMidy(x.getSeries(i), xmid);
			
			while(m.containsKey(y[0])) y[0] +=0.00001;
			m.put(y[0], new int[] {i,(int) (y[1]*(180.0/Math.PI))});
		}
		return m;
	}
	
public	static class Line{
		XYSeries series;
		
		SortedMap<Double, Double> xy = new TreeMap<Double, Double>();
		public Line(XYSeries series, boolean logx, boolean logy){
			this.series = series;
			for(int i=0; i<series.getItemCount(); i++){
				XYDataItem item = series.getDataItem(i);
				double x = item.getX().doubleValue();
				double y = item.getY().doubleValue();
				if(logx) x = Math.log(x);
				if(logy) y = Math.log(y);
				//if(x>=start && x<=end){
					xy.put(x,y);
				//}
			}
		}
		
		public Line(SortedMap<Double, Double> xy, boolean logx, boolean logy){
			//this.series = series;
			this.xy = xy;
			
		}
		public double getStart(){
			return xy.firstKey();
		}
		public double getLast(){
			return xy.lastKey();
		}
		public double mid(){
			return (getStart()+getLast())/2.0;
		}
		public double get(double x){
			if(xy.size()==0) return 0;
			SortedMap<Double, Double> tail = xy.tailMap(x);
			if(tail.size()==0) return xy.get(xy.lastKey());
			 double x2 = tail.firstKey();
			double y2 = tail.get(x2);
			if(Math.abs(x2-x)<1e-8) return y2;
			SortedMap<Double, Double> head = xy.headMap(x);
			if(head.size()==0) return y2;
			double x1 = head.lastKey();
			double y1 = head.get(x1);
			double sl =  (y2-y1)/(x2-x1);
			return y1 + sl*(x-x1);
		}
		public double getSlopeRight(double x){
			SortedMap<Double, Double> tail = xy.tailMap(x);
			if(tail.size()<4) return 0;
			else{
				Iterator<Entry<Double, Double>> it = tail.entrySet().iterator();
				Entry<Double, Double> first = it.next();
				Entry<Double, Double> last=null;
				for(int i=0; i<3; i++){
					last = it.next();
				}
				return Math.atan(last.getValue() - first.getValue())/(last.getKey()-first.getKey());
			}
		}
		public boolean excludeSlope(double x, double thresh){
			SortedMap<Double, Double> tail = xy.tailMap(x);
			SortedMap<Double, Double> head = xy.headMap(x);
			if(head.size()==0|| tail.size()==0) return true;
		
			else{
				
				double x1 = head.lastKey();
				double x2 = tail.firstKey();
				double y2 = xy.get(x2);
				double y1 = xy.get(x1);
				double sl =  (y2-y1)/(x2-x1);
				return Math.abs(Math.atan(sl))>thresh;
			}
		}
		public boolean localSlope(double x, double y, double d, double[] xs, double[] ys, boolean right){
			if(right){
				xs[0] = x;
				ys[0] = y;
				xs[1] = d+x;
				ys[1] = this.get(xs[1]);
				
			}
			else{
				xs[1] = x;
				ys[1] = y;
				xs[0] = x-d;
				ys[0] = this.get(xs[0]);
			}
			return true;
			/*
			if(right){
			SortedMap<Double, Double> tail = xy.tailMap(x);
			if( tail.size()<2) return false;
				xs[0] = tail.firstKey();
				tail = tail.tailMap(xs[0]+d);
				if(tail.size()==0) return false;
				xs[1] = tail.firstKey();
			}
			else{
				SortedMap<Double, Double> head = xy.headMap(x+0.000001);
				if(head.size()<2) return false;
				xs[1] = head.lastKey();
				head = head.headMap(xs[1]-d);
				if(head.size()==0) return false;
				xs[0] = head.lastKey();
			}
		
				 ys[0] = xy.get(xs[0]);
			     ys[1] = xy.get(xs[1]);
				return true;*/
			
		}

		public String getKey() {
			return this.series.getKey().toString();
		}

		public double getLt(double x) {
			return xy.get(this.xy.headMap(x).lastKey());
		}
		
	}
	
	static class MaxSepInd{
		//int[] len,currInd;
	//	double[] maxSep;
	//	boolean[] top;
		//Double[] currX,currY, nxtX, nxtY;
		 double yu,yl,xu,xl;
		final double start, end; 
		XYSeries[] ser;
		Line[] line;
		double[] xs = new double[2];
		double[] ys = new double[2];
	boolean logx, logy;
	double[] x,y;
	final boolean useSlope;
		MaxSepInd(XYSeries[] ser, Range xrange,Range yrange, boolean logx, boolean logy, boolean useSlope){
			this.ser = ser;
			this.useSlope = useSlope;
			this.logx = logx;
			this.logy = logy;
			this.x = new double[ser.length];
			this.y = new double[ser.length];
			this.line = new Line[ser.length];
		 xl =  xrange.getLowerBound();
			 xu =  xrange.getUpperBound();
			 yl = yrange.getLowerBound();
			 yu = yrange.getUpperBound();
			if(logx){
				xl =Math.log(xl);
				xu = Math.log(xu);
			}
			if(logy){
				yl =Math.log(yl);
				yu = Math.log(yu);
			}
			
			start = xl;
			end = xu;
			this.xlen = end - start;
			this.ylen = yu - yl;
			for(int i=0; i<ser.length; i++){
				//len[i] = ser[i].getItemCount();
				line[i] = new Line(ser[i],logx,logy);
			
			}
		
		}
		public double[] getBest(int k) {
			double x  = this.abovex[k];
			//double slope = line[k].localSlope(x);
		
			return new double[] {logx ? Math.exp(x) : x, logy ? Math.exp(abovey[k]) : this.abovey[k],
					slope[k],
					line[k].getSlopeRight(x)};
			
		}
		public boolean above(int k) {
			return this.above[k];
		}
		
	//	public int[] ind;
	//	double x_current =0;
		
		double[][] vals;
		double[] maxAbove,  abovex, abovey;
		boolean[] above;
		//maxBelow,belowx;
		final double ylen,xlen;
		
		 double[] slope;
		public void run(){
			vals = new double[ser.length][2];
			maxAbove = new double[ser.length];
			Arrays.fill(maxAbove, -1.0);
			//maxBelow = new double[ser.length];
			abovex = new double[ser.length];
			abovey = new double[ser.length];
			above = new boolean[ser.length];
			slope = new double[ser.length];
			double mid = (xu+xl)/2.0;
			//belowx = new double[ser.length];
			Set<Integer> lo = new HashSet<Integer>();
			for(int k=0; k<ser.length; k++){
				lo.add(k);
			}
			int len = vals.length-1;
			{
				
			//	double mid = (start+end)/2.0;
				Arrays.fill(abovex,mid);
				double ycorr = ylen/xlen;
				//Arrays.fill(belowx,mid);
			//	double end1 = end - 0.1*(xlen);
				double angThresh = (30./180.0) *Math.PI;
				List<Double > xl = getList(mid, this.xl + ROC.overlap*xlen, xu - ROC.overlap*xlen,1000);
				double maxDiffInd=mid;
				double maxDiff =0;
				if(true){
				for(int ii=0; ii<xl.size(); ii++){
					double x = xl.get(ii);
					for(int j=0; j<ser.length; j++){
						vals[j][0] = this.line[j].get(x);
						vals[j][1] = j;
					}
					Arrays.sort(vals,len2);
					double totalDiff = 0;
					for(int j=0; j<vals.length; j++){
						double y = vals[j][0];
						double bel = j==0 ? Double.POSITIVE_INFINITY : y-vals[j-1][0];
						double ab = j==len ? Double.POSITIVE_INFINITY : vals[j+1][0]-y;
						double max = Math.max(bel, ab);
						if(j>1) totalDiff +=bel;
						int ind = (int)vals[j][1];
						double mind=Double.POSITIVE_INFINITY;
						if(max > maxAbove[ind]){
							//double slope = 
							boolean exclude=false;//this.line[ind].excludeSlope(x, angThresh);
							
							if( (bel>ab && x < this.xl+0.05*xlen)){
								exclude =true;
							}
							 if( (ab> bel && x > this.xu-0.05*xlen)){
								 exclude =true;
							 }
							 if(y> this.yu - 0.2*ylen) exclude = true;
							inner: for(int j1=vals.length-1; !exclude && j1>=0; j1--){
								if(j1==j) continue inner;
								int ind1 = (int) vals[j1][1];
								double x1 = abovex[ind1];
								double y1 = abovey[ind1];
								double diffy = Math.abs(y-y1);
								double diffx = Math.abs(x-x1);
								double m = Math.max(diffx/xlen, diffy/ylen);
								if(diffx < ROC.overlap*xlen && diffy < ROC.overlap*ylen){
									exclude = true; 
								
								}
								if(ROC.useEnd){
									x1 = x;
									y1 = this.line[ind1].get(x1);
									if(ab>bel && y1>y && y1-y<ROC.overlap) exclude =true;
									else if(ab<bel && y1<y && y-y1< ROC.overlap) exclude = true;
								}
								if(m<mind){
									 mind = m;
								}
								
							}
							if(!exclude ){
								boolean hasSlope =
									useSlope ? mind >= ROC.overlap && line[ind].localSlope(x, y,xlen/10.0,xs,ys, ab>bel) : false;
								slope[ind] = hasSlope ? calcSlope(xs, ys,  ycorr) : 0;
								double slopet = (180.0/Math.PI)*slope[ind];
								lo.remove(ind);
								maxAbove[ind] = ab;
								abovex[ind] = x;
								abovey[ind] = y;
								above[ind] = ab> bel;
							
							
								boolean maxx =y > this.yu  - 0.1 *this.ylen;
								boolean minxx = y < this.yl+0.1*this.ylen;
						//	System.err.println(ser[ind].getKey()+" "+abovex[ind]+" "+abovey[ind]+" "+above[ind]+" "+slopet);
								if(maxx || minxx) slope[ind] =0;
							}
						}
						
					}
				}
				}
				if(lo.size()>0){
					maxDiffInd = ROC.useEnd ? (this.xu*8+this.xl*2)/10.0 : (this.xu+this.xl)/2.0;
					if(!ROC.useEnd){
					for(int ii=0; ii<xl.size(); ii++){
						double x = xl.get(ii);
						for(int j=0; j<ser.length; j++){
							vals[j][0] = this.line[j].get(x);
							vals[j][1] = j;
						}
						Arrays.sort(vals,len2);
						double totalDiff = 0;
						for(int j=1; j<vals.length; j++){
							double y = vals[j][0];
							double bel = j==0 ? Double.POSITIVE_INFINITY : y-vals[j-1][0];
							if(j>1) totalDiff +=bel;
						}
						if(totalDiff>maxDiff && x > this.xl+0.1*xlen && x < this.xu-0.1*xlen && mindist(x,abovex)>2*ROC.overlap){
							maxDiff = totalDiff;
							maxDiffInd = x;
							
						}
					}
					}
					int kk=0;
				final	double x = maxDiffInd;
					for(int j=0; j<ser.length; j++){
						vals[j][0] = this.line[j].get(x);
						vals[j][1] = j;
					}
					Arrays.sort(vals,len2);
					for(int j=0; j<vals.length; j++){
						double y = vals[j][0];
						double bel = j==0 ? Double.POSITIVE_INFINITY : y-vals[j-1][0];
						double ab = j==len ? Double.POSITIVE_INFINITY : vals[j+1][0]-y;
						double max = Math.max(bel, ab);
						int nxt = (int)vals[j][1];
						if(lo.contains(nxt)){
							System.err.println("assigning to mid "+maxDiffInd);
							abovey[nxt] = y;
							abovex[nxt] = maxDiffInd;
							above[nxt] =  ab> bel;
							boolean hasSlope =false;
//								useSlope ? line[nxt].localSlope(maxDiffInd, abovey[nxt],xlen/10.0,xs,ys, above[nxt]) : false;
							slope[nxt] = hasSlope ? calcSlope(xs, ys,  ycorr) : 0;
						}
					}
				}
			}
			//public void 
			
		}
		private double mindist(double x2, double[] abovex2) {
			double mind = Double.POSITIVE_INFINITY;
			for(int i=0; i<abovex2.length; i++){
				if(Math.abs(abovex2[i]-x2)<mind){
					mind = Math.abs(abovex2[i]-x2);
				}
			}
			return mind;
		}
		private double calcSlope(double[] xs2, double[] ys2, double ycorr) {
			return Math.atan2((ys2[1] - ys2[0])/ycorr, xs2[1] - xs2[0]);
		}
		private static List<Double> getList(double mid, double start1, double end1, int cnt) {
			
			double incr =( end1 - start1)/(cnt);
			List<Double> l = new ArrayList<Double>();
			for(int i=0; incr*i + start1 < end1; i++){
				l.add(incr*i + start1);
			}
	/* incr =(  mid-start1)/(cnt);
			  for(int i=0; i<cnt; i++){
				l.add(mid-incr*i);
			}*/
			return l;
		}
		
	
	}
	
	
	

	private static double[] getMidy(XYSeries series, double xmin) {
		for(int i=0; i<series.getItemCount(); i++){
			if(series.getX(i).doubleValue()>xmin && i>0){
				double y1 = series.getY(i-1).doubleValue();
				double y2 = series.getY(i).doubleValue();
				double x1 = series.getX(i-1).doubleValue();
				double x2 = series.getX(i).doubleValue();
				double slope = 0;
				if(series.getItemCount()>5 && i>3){
				XYDataItem start = series.getDataItem(Math.max(3,i-3));
				XYDataItem end = series.getDataItem(Math.min(series.getItemCount()-2,i+5));
				
			slope = Math.atan((end.getYValue()-start.getYValue())/(end.getXValue()-start.getXValue()));
				
				}
				double y = y1+ ((y2-y1)/(x2-x1))*(xmin-x1);
				return new double[] {y,slope};
			}
		}	
		if(series.getItemCount()>0) return new double[] {series.getY(0).doubleValue(),0.0};
		else return new double[] {0.0,0.0};
	}

	static double getAngle(XYSeries series, int i){
		int start_ind = Math.max(3,i-3);
		int end_ind = Math.min(series.getItemCount()-2,i+5);
		if(end_ind - start_ind<2) return 0;
		XYDataItem start = series.getDataItem(start_ind);
		XYDataItem end = series.getDataItem(end_ind);
		
		double ang = Math.atan((end.getYValue()-start.getYValue())/(end.getXValue()-start.getXValue()));
		//double ang1 = (180/Math.PI)*ang;
		return ang;
	}
	private static double[] getAnnotCoords(XYSeries series, Range xrange,
			Range yrange, SortedMap<Double, Double> prevv, double y, boolean xlog, boolean ylog) {
		double xu = xrange.getUpperBound();
		double xl = xrange.getLowerBound();
		double yu = yrange.getUpperBound();
		double yl = yrange.getLowerBound();
		if(xlog){
			xu = Math.log(xu);
			xl = Math.log(xl);
		}
		if(ylog){
			yu = Math.log(yu);
			yl = Math.log(yl);
			y = Math.log(y);
		}
		
		
		double xlen = xu-xl;
		double xmin = (xl+xu)/2.0;
		double ylen = yu-yl;
		
				//double yl = yrange.getLowerBound()-ylen;
			
		
				
				SortedMap<Double, Double> head = prevv.tailMap( y - 0.1*ylen);
				if(head.size()>0){
					head = head.headMap( y+0.1*ylen);
				}
				if(head.size()>0){
					SortedSet<Double> l = new TreeSet<Double>(head.values());
					inner: for(int k=0; ;k++){
						double xless = xmin - 0.25*k*xlen;
						double xmore = xmin + 0.25*k*xlen;
						if(!l.contains(xless)){
							xmin = xless;
							break inner;
						}
						else if(!l.contains(xmore)){
							xmin = xmore;
							break inner;
						}
					}
					
				}
			
				
				prevv.put(y, xmin);
				return new double[] {xlog ? Math.exp(xmin) : xmin, ylog ? Math.exp(y) : y};
		
	}

	private static Shape resize(Shape seriesShape, double d) {
		// TODO Auto-generated method stub
		return null;
	}
	private static void correctForLog(XYSeriesCollection xys, boolean logx,
			boolean logy) {
		Logger.global.info("correcting for log");
		for (int i = 0; i < xys.getSeriesCount(); i++) {
			XYSeries series = xys.getSeries(i);
			List<double[]> toadd = new ArrayList<double[]>();
			for (int k = series.getItemCount() - 1; k >= 0; k--) {
				XYDataItem di = series.getDataItem(k);
				Number x = di.getX();
				Number y = di.getY();
				if ((logx && x.doubleValue() <= 1e-8)) {
					series.remove(k);
					if (y.doubleValue() != 0.0)
						toadd.add(new double[] { 1, y.doubleValue() });
				} else if (logy && y.doubleValue() <= 1e-8) {
					series.remove(k);
				}
			}
			for (int k = 0; k < toadd.size(); k++) {
				double[] ti = toadd.get(k);
				series.add(ti[0], ti[1]);
			}
		}

	}
	private static void adjustForZero(XYSeriesCollection xys, boolean logx,
			boolean logy) {
		Logger.global.info("correcting for log");
		for (int i = 0; i < xys.getSeriesCount(); i++) {
			XYSeries series = xys.getSeries(i);
			//List<double[]> toadd = new ArrayList<double[]>();
			double y_max = Double.NEGATIVE_INFINITY;
			boolean add = false;
			for (int k = series.getItemCount() - 1; k >= 0; k--) {
				XYDataItem di = series.getDataItem(k);
				Number x = di.getX();
				Number y = di.getY();
				if ( x.doubleValue() <= 1e-10){
					add=true;
					if(y.doubleValue()>y_max){
						y_max = y.doubleValue();
					}
					series.remove(k);
					
				} 
			}
			if(add){
				series.add(0,y_max);
			}
			
		}

	}

	public static ChartPanel getChartPanel1(XYSeriesCollection xys,
			List<Color> cols, List<Stroke> strokes, String probeset,
			String result_set_names, Integer[] cn, boolean max, int total,
			String chr) {
		int ind = chr.indexOf(".zip");
		String chr1;
		if (ind >= 0) {
			chr1 = chr.substring(0, ind);
		} else {
			chr1 = chr;
		}
		final JFreeChart chart = ChartFactory.createXYLineChart("probes on: "
				+ probeset + "\n chromosome: " + chr1
				+ "\n Using results generated from " + result_set_names
				+ "\n  cn = " + Arrays.asList(cn) + "\n  max = " + max,
				"Probability ", // domain axis label
				"Accuracy", // range axis label
				xys, // data
				PlotOrientation.VERTICAL, true, // include legend
				false, // tooltips?
				false // URL generator? Not required...
				);
		// chart.getXYPlot().setDomainAxis(new
		// LogarithmicAxis("False positive"));
		// XYItemRenderer rend = chart.getXYPlot().getRenderer();
		final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setShapesFilled(Boolean.TRUE);
		renderer.setLinesVisible(Boolean.TRUE);

		ValueAxis axis = chart.getXYPlot().getDomainAxis();

		chart.getXYPlot().setRenderer(renderer);
		for (int i = 0; i < cols.size(); i++) {

			if (strokes.get(i) != null)
				renderer.setSeriesStroke(i, strokes.get(i));
			if (cols.get(i) != null)
				renderer.setSeriesPaint(i, cols.get(i));
		}
		chart.getTitle().setFont(font10);
		((XYPlot) chart.getPlot()).getDomainAxis().setLabelFont(font8);
		((XYPlot) chart.getPlot()).getDomainAxis().setTickLabelFont(font8);
		axis.setRange(0, 0.1 * axis.getRange().getUpperBound());
		return new ChartPanel(chart);
	}

	public static Font font10 = new Font("Arial", Font.BOLD, (int)Math.round(10.0*resolution));
	public static Font font8 = new Font("Arial", Font.PLAIN, 8);
	public static Font font5 = new Font("Arial", Font.PLAIN, 5);
	public static int mapLength = 3;
	String[] header = new String[] { "pos/neg", "index_p", "k_p", "loc",
			"snpid", "individ", "length" };
	String form = "%10s %5.3g %5.3g %10s %10s";
	String form_h = "%10s %5s %5s  %10s %10s";

	static class Entr implements Comparable {
		// final Integer loc;
		String str;

		Entr(double[] probs) {
			System.arraycopy(probs, 0, prob, 0, probs.length);
			// this.loc = loc;
		}

		double[] prob = new double[2];

		public String toString() {
			return str + "_" + prob[0] + "_" + prob[1];
		}

		public int compareTo(Object o) {
			Entr ent = ((Entr) o);

			// for(int i=1;i>=0; i--)
			{
				double d2 = ent.prob[1];
				double d1 = prob[1];
				if (d1 != d2) {
					return d1 < d2 ? 1 : -1;
				}
			}
			{
				double d2 = ent.prob[0];
				double d1 = prob[0];
				if (d1 != d2) {
					return d1 > d2 ? 1 : -1;
				}
			}
			return 0;
			/*
			 * int res = ent.str.compareTo(str); if(res==0){ throw new
			 * RuntimeException("!!"); } return res;
			 */
		}
	}

	//boolean equal = false;

	public void setComp(DataCollection toComp) {
		this.toComp = toComp;
		if(toComp instanceof LightWeightDataCollection){
			midC =((LightWeightDataCollection)toComp).midpoint;
		}
		else{
			midC = ((LightWeightDataCollection)((LightWeightMergedDataCollection)toComp).ldl[0]).midpoint;
		}
		//equal = toComp.loc().equals(this.base.loc()); // throw new
		// RuntimeException(toComp.name+" vs "+this.base.name+"\n"+toComp.loc.get(0)+"-"+toComp.loc.get(toComp.size()-1)+
		// " vs "+base.loc.get(0)+"-"+base.loc.get(base.size()-1));

	}

}
