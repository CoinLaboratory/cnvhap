package lc1.dp.appl;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import lc1.dp.data.classification.ROC;
import lc1.dp.data.classification.ROC.Inner;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.Info;
import lc1.dp.data.collection.LightWeightDataCollection;
import lc1.dp.data.collection.LightWeightMergedDataCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.swing.ColorAdapter;
import lc1.dp.swing.Headless;
import lc1.ensj.GenesForRegion;
import lc1.util.Constants;
import lc1.util.Executor;

import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.PageConstants;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.block.LengthConstraintType;
import org.jfree.chart.block.RectangleConstraint;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Size2D;

public class CompareBasic {
	//private static int[] st = new int[] {17,25, 72};
//	private static int maxRSToCount = 1000;
	public static boolean combine = false;
	public static ColorAdapt cola = new ColorAdapt();
	public static class  ColorAdapt{
		Map<String, Color> m = new HashMap<String, Color>();
		public Color getColor(String col) throws Exception{
			Color col1 = m.get(col);
			if(col1==null){
				String[] str = col.split("-");
				col1 =(Color) Color.class.getField(str[str.length-1].trim()).get(null);
				if(str.length==2){
					if(str[0].equals("LIGHT")){
						col1 = ColorAdapter.merge(new Color[] {Color.WHITE, col1});
					}
					else{
						col1 = ColorAdapter.merge(new Color[] {Color.black, col1});
					}
				}
				m.put(col, col1);
			}
			return col1;
		}
	}
	
	private static final int START =0;//st[1];
	
 private static final int END =500;//START+1;//73;//18;// 26;//*1000*1000;//73000000;

	// String[] chrom = new String[] {"X"};maxcoun
 

	public static Set<String> getChroms(File res_set, String probeset) {
		String[] f = (new File(res_set, probeset)).list(new FilenameFilter() {

			public boolean accept(File pathname, String nme) {
				return nme.endsWith("zip");
			}

		});
		return new HashSet<String>(Arrays.asList(f));
	}

	public static Set<String> getCommon(File[] res_set, String[] probeset) {
		Set<String> res = getChroms(res_set[0], probeset[0]);
		for (int i = 0; i < res_set.length; i++) {
			for (int j = 0; j < probeset.length; j++) {
				res.retainAll(getChroms(res_set[i], probeset[j]));

			}

		}
		return res;

	}

	static int strlength = 1;

	public static String getExpt(String st) {
		if (true)
			return st;
		// if(!useAll) return null;
		if (st.startsWith("_") && st.indexOf("__") >= 0) {
			String[] str = st.split("\\s+");
			String[] expt = str[0].split("__");
			String res = expt[1].split("_")[0];
			return res;
		} else {
			return null;
		}
	}

	// static boolean useAll = false; //use experiment structure
	public static List<String> expt;/* = Arrays.asList(
			("1M_PennCNV 244k_agilent1 1M_cnvP 244k_agilent1 185k_agilent "+
			"1M	   317k	317k_1M 1M_hap "+
				
			"244k	185k	185k_244k "+
			"244k_1M 185_317_cnvHap 185_1M_cnvHap 244_317_cnvHap ")
					.split("\\s+") );*/

	public static List<String> col ;/*=Arrays.asList(
		(	"GRAY BLACK BLUE BLACK BLUE BLACK "+
			"RED	   	CYAN	MAGENTA YELLOW ORANGE BLUE PINK "+
			"RED	GREEN	PINK "+
			"CYAN YELLOW ORANGE MAGENTA").split("\\s+")) ;*/
public static File results1;
public static int[] centromere;
	public static void main(final String[] args) {
		// System.err.println(Format.sprintf("%7.7g",new Double[]
		// {Double.NaN}));
		// System.exit(0);
		// Color.
		// Constants.modelCNP = 10;
	System.err.println(Arrays.asList(args));
		File results = new File(new File(System.getProperty("user.dir")),"results_compare");
		results.mkdir();
		ROC.makeAsRate = Boolean.parseBoolean(args[1]);
		results1 = new File(results, Arrays.asList(args).subList(1, 5).toString().replaceAll("\\s+", "_").replace("[", "").replace("]", "").replace(':', '^'));
		if(results1.exists()){
			Constants.delete(results1);
		}
		results1.mkdir();
		Arrays.fill(Constants.var_thresh, 100);
		expt = Arrays.asList(args[offset+4].split(":"));
		col = Arrays.asList(args[offset+5].split(":"));
	//	if(args.length>7+offset )maxcount = Integer.parseInt(args[offset+7]);
		System.err.println("args are "+Arrays.asList(args));
;		String[] chr = args[offset+0].split(":");
		if(args[offset+0].startsWith("split.txt")){
			try{
			BufferedReader br = new BufferedReader(new FileReader("split.txt"));
			StringBuffer sb = new StringBuffer();
			String st = "";
			for(int i=0; (st = br.readLine())!=null; i++){
				String st1 = st.split("\\s+")[1].replace(':', ';');
			//if(st1.split(";")[0].equals("3") || st1.split(";")[0].equals("10")){
				sb.append(st1);
				 sb.append(":");
		//	}
			}
			String st1  =  sb.toString();
			chr =st1.substring(0, st1.length()-1).split(":");
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
	if(chr.length==1 && chr[0].equals("all")){
		chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22".split(":");
	}
	 GenesForRegion grf = new GenesForRegion();
	centromere = new int[chr.length];
	for(int k=0; k<centromere.length; k++){
		try{
		centromere[k] = grf.getCentromere(chr[k].split(";")[0]);
		}catch(Exception exc){
			exc.printStackTrace();
			System.exit(0);
		}
	}
		boolean all = true;
		//for (int i = 0; i < chr.length; i++) {
			if (true) {
				String []chr1 = new String[chr.length];
				int[] st1 = new int[chr.length];
				int[] end1 = new int[chr.length];
				for(int i=0; i<chr.length; i++){
					String[] str = chr[i].split(";");
					chr1[i] = str[0];
					st1[i] = str.length>1 ? Integer.parseInt(str[1]) : START*1000*1000;
					end1[i] = str.length>2 ? Integer.parseInt(str[2]) : END*1000*1000;
					if(!combine)main(new String[] {chr1[i]}, new int[] {st1[i]}, new int[] {end1[i]},args, i>0);
				}
				if(combine)main(chr, st1, end1,args,false);
						
			/*} else {
				List<String> experiments = new ArrayList<String>();
				String[] str = new File(System.getProperty("user.dir"))
						.list(new FilenameFilter() {

							public boolean accept(File dir, String name) {
								if (name.indexOf(args[offset+5]) < 0
										|| name.indexOf(".zip") >= 0)
									return false;
								File f = new File(dir, name);
								if (f.isDirectory())
									return true;
								return false;
							}

						});
				for (int j = 0; j < str.length; j++) {
					String expt = getExpt(str[j]);
					if (expt != null && !experiments.contains(expt)

					)
						experiments.add(expt);
				}
				for (int j = 0; j < experiments.size(); j++) {
					try {
						main(chr[i] + ".zip", args);
					} catch (Exception exc) {
						System.err
								.println("problem with " + experiments.get(j));
						exc.printStackTrace();
					}
				}*/
		//	}
		}
	}

	/*
	 * public static void main1(){ List<String> experiments = new
	 * ArrayList<String>(); if(useAll){ String[] str = (new File(".")).list(ff);
	 * for(int i=0; i<str.length;i++){ String expt = getExpt(str[i]);
	 * if(expt==null) common.add(str[i]); if(expt!=null &&
	 * !experiments.contains(expt)
	 * 
	 * ) experiments.add(expt); } }
	 * 
	 * if(experiments.size()==0) experiments.add("base_exp"); int[] base_ind =
	 * new int[experiments.size()]; for(int i=0; i<base_ind.length; i++){
	 * base_ind[i] = i; } }
	 */
	public static JTabbedPane[] jpane;
	public static String[] probsets;
	public static String probsets1;

	static File buildF = null;// new
	private static boolean showPlots = true;

	static int offset = 2;
								// File("/home/lcoin/Data/data1/244k/build_excl.txt");
static  int maxcount =Integer.MAX_VALUE;
static int[] minCount =new int[] {0};
		//	= new int[] {10,10,0};
static int[] minLength= new int[] {0};
//			=new int[] {0,0,0};
public static int buff = LightWeightDataCollection.midpoint;	

static String[] rsbuff = new String[2*buff+1];
static String rsmid=null;
 static boolean[] canProcess = new boolean[2*buff+1];


	private static List<String> getNme(List<String> s,
			 String toch) {
		boolean all = true;
		List<String> nme = new ArrayList<String>();
		for(int k=0; k<s.size(); k++){
			String nme1 = ROC.process(s.get(k));
			if(!nme1.startsWith(toch)) all  = false;
			if(nme.contains(nme1)){
				nme.add(nme1+"(b)");
			}
			nme.add(nme1);
		}
		if(all){
			for(int k=0; k<nme.size(); k++){
				String str = nme.get(k);
				nme.set(k,str.substring(toch.length()+1, str.length()-1));
			}
		}
		// TODO Auto-generated method stub
		return  nme;
	}
 public static void main(final String[] chr,final int[] start, final int[] end,  final String[] args1, boolean append
			) {
		try {
			
		//
			probsets = args1[1+offset].split(":");
			probsets1 = args1[1+offset].replace(':', '_');
			//Constants.build = "build36_317k.txt";
			buildF = args1.length<=(6+offset) || args1[6+offset].equals("null") ? null :  new File(args1[0]+"/"+args1[6+offset]);
			if(args1.length>7+offset){
				if(!args1[7+offset].equals("null")){
				String[] str = args1[7+offset].split(":");
				ROC.len_bands = new int[str.length];
				for(int k=0; k<str.length; k++){
					ROC.len_bands[k] = Integer.parseInt(str[k]);
				}
				}
				if(args1.length>8+offset){
					ROC.countRegs = Boolean.parseBoolean(args1[8+offset]);
					
				}
				if(args1.length>9+offset){
				//	ROC.min_len=Integer.parseInt(args1[9+offset]);
					
					String[] str1 = args1[9+offset].split(":");
			
					minCount =convert(str1[0].split(";"));
					if(str1.length>1){
					minLength = convert(str1[1].split(";"));
					}
					//LightWeightDataCollection.midpoint = minCount;
					//LightWeightDataCollection.last = 2*minCount+1;*/
				}
				
			//	ROC.minNumberTrue = ROC.countRegs ? 1 :0;
			}
		/*	minCount = 0;//Integer.parseInt(args1[9+offset]);
			LightWeightDataCollection.midpoint = minCount;
			LightWeightDataCollection.last = 2*minCount+1;*/
			if(probsets.length>1){
				ROC.combine = false;
			}
		//	System.err.println(System.getProperties());
			String headless = System.getProperty("java.awt.headless");
			Constants.loess = null;
			Constants.gc = null;
			Constants.modify0 = new char[Constants.inputDir.length][];
			Arrays.fill(Constants.modify0, new char[] { '0', '1', '2' ,'3'});
		
			
		//	boolean restrict = true;
			final String[] exp = args1[4+offset].split(":");
//			final String date = args1[3].split("-")[1];;

			
			//List<String> result_set_name = new ArrayList<String>();
			List<Color> result_set_color = new ArrayList<Color>();
		//	List<Stroke> result_set_stroke = new ArrayList<Stroke>();
			final File dir1 = new File(args1[0]);
			String[] cols = args1[5+offset].split(":");
			List<String> str_base = new ArrayList<String>(Arrays.asList(args1[3+offset].split(":")));
			List<String> str_comp = Arrays.asList(args1[4+offset].split(":"));
			if(str_base.size()<str_comp.size()){
				int len1 = str_base.size();
				for(int i=len1; i<str_comp.size(); i++){
					str_base.add(str_base.get(len1-1));
				}
			}
			List<String> str_comp_alias = getNme(str_comp, "cnvHap");
			
			//List<String> cols = Arrays.asList(args1[5].split(":"));
			List<String> all = new ArrayList<String>(str_comp);
			List<Stroke> result_set_stroke = new ArrayList<Stroke>();
			int[] alias_base = new int[str_base.size()];
			for(int i=0; i<alias_base.length; i++){
				int ind = all.indexOf(str_base.get(i));
				if(ind>=0){
					alias_base[i] = ind;
				}
				else{
					alias_base[i] = all.size();
					all.add(str_base.get(i));
				}
			}
			for(int i=0; i<str_comp.size(); i++){
				/*result_set_stroke.add( new BasicStroke(0.8f,
						BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,
						generate(5), 1.0f));*/
				if(i<cols.length) {
					try{
					Color c = cola.getColor(cols[i]);//.getColor(cols[i]);
					//Color c1 = Color.getColor(cols[i].toLowerCase());
					result_set_color.add(c);
					}catch(Exception exc){
						exc.printStackTrace();
					}
				}
				else result_set_color.add(Color.getHSBColor((float)Math.random(),(float) Math.random(), (float)Math.random()));
			}
			int szec = result_set_color.size();
			for(int i=szec; i<cols.length; i++){
				Color col = cola.getColor(cols[i]);
				result_set_color.add(col);
			}
			String[] noc_ = args1[2+offset].split(":");
			Integer[][] nocop = new Integer[noc_.length][];
			Set<Integer> cops = new HashSet<Integer>();
			for (int i = 0; i < noc_.length; i++) {
				String[] noci = noc_[i].split(";");
				nocop[i] = new Integer[noci.length];
				for (int j = 0; j < noci.length; j++) {
					nocop[i][j] = Integer.parseInt(noci[j]);
					cops.add(nocop[i][j]);
				}
			}
			DataCollection[][] baseComp = new DataCollection[str_base
			                           					.size()][chr.length];
			DataCollection[][] toComp = new DataCollection[str_comp
					.size()][chr.length];
		
			
		/*	if(toComp.length>8){
				ROC.len_bands = new int[] {ROC.len_bands[ROC.len_bands.length-1]};
			}*/
			Set<String> indiv = new HashSet<String>();
			DataCollection[][] allD = new DataCollection[all.size()][chr.length];
			int non_neg_ind = -1;
		//	for (int i = 0; i < probsets.length; i++) {
		
			List[][] indiv_ = new List[all.size()][chr.length];
		
			
		    int nonnull = -1;
			for(int km =0; km<chr.length; km++){
				 int[] mid = new int[] { start[km] ,end[km] };	
				 if(!combine){
				 List<String > l  = readBuild(buildF, "chr"+chr[km], start[km], end[km]);
					double sze = l.size();
					if(sze==0){
						System.err.println("something wrrong");
						System.exit(0);
					}
					LightWeightDataCollection.midpoint = (int) Math.floor(sze/2.0);
					buff = LightWeightDataCollection.midpoint ;
					canProcess = new boolean[2*buff+1];
					rsbuff = new String[2*buff+1];
					 rsmid = l.get(buff);
				 
					LightWeightDataCollection.last = 2*buff+1;
				 }
				List<String>[] snpid = new List[probsets.length];
				List<Integer>[] locs = new List[probsets.length];
				List<Integer> locsComb = null;
				List<String> snpsComb = null;
				Info[][]map=null;
					Constants.mid = new String[][] {new String[] {chr[km].split("\\.")[0]}};
				for (int k1 = 0; k1 < all.size(); k1++) {
					try {
						String[] res_name1 = all.get(k1).split("\\s+");
						LightWeightDataCollection[] dcin = new LightWeightDataCollection[probsets.length];
						for(int jk=0; jk<dcin.length;jk++){
							LightWeightDataCollection[] dcin1 = new LightWeightDataCollection[res_name1.length];
							for(int jk1 = 0; jk1<dcin1.length; jk1++){
							File res_set = new File(dir1, res_name1[jk1]
									+ "/res/");
							File f1 = new File(res_set, probsets[jk] + "/" + chr[km]+".zip");
							System.err.println("reading file "+f1.getAbsolutePath());
							if (!f1.exists()) {
								//System.err.println("WARNING DID NOT EXIST" + f1);
								File[] f2 = res_set.listFiles();
								f1 = new File(f2[0],chr[km]+".zip");
								if (!f1.exists()) {
								/*	File gp = f1.getParentFile().getParentFile().listFiles(new FileFilter(){

										@Override
										public boolean accept(File pathname) {
											return pathname.isDirectory();
										}
										
									})[0];
									if(!gp.exists()) throw new RuntimeException("!!");
									f1 = new File(gp, f1.getName());*/
									
									f1=null;
								}
							//	continue;
							}
							if(f1!=null){
							if(k1==0){
								dcin1[jk1] =null;/* new  LightWeightDataCollection(
									f1, (short) 0, chr[km].indexOf('X') >= 0
											|| chr[km].indexOf('Y') > 0 ? 1 : 2,
									new int[][] { mid },buildF, buff);
								snpid[jk] = dcin1[jk1].snpid;
								locs[jk] = dcin1[jk1].loc;*/
							}
							else{
								dcin1[jk1] =null;/* new LightWeightDataCollection(
										f1, (short) 0, chr[km].indexOf('X') >= 0
												|| chr[km].indexOf('Y') > 0 ? 1 : 2,
										new int[][] { mid },locs[jk], snpid[jk], buff);*/
							}
							dcin1[jk1].alleleA.clear();
							dcin1[jk1].alleleB.clear();
							dcin1[jk1].strand.clear();
							}
							else{
								System.err.println("did not exist");
							}
							}
							
							if(dcin1.length==1) {
								dcin[jk]= dcin1[0];
							}
							else{
								dcin[jk] = new LightWeightDataCollection2(dcin1);
							}
						}
						DataCollection ldc=null;
						
						if(dcin.length==1) ldc = dcin[0];
						else{
							try{
							if(k1==0){
								ldc = new LightWeightMergedDataCollection(dcin, "merged");
								map = ((MergedDataCollection)ldc).map;
							}
							else{
								ldc = new LightWeightMergedDataCollection(dcin, "merged", map);
							}
							}catch(Exception exc){
								exc.printStackTrace();
								System.err.println("prob with "+all.get(k1)+" "+chr[km]);
								System.exit(0);
							}
						}
						if(ldc!=null){
						if(locsComb==null){
							locsComb = ldc.loc;
							snpsComb = ldc.snpid;
						}
						//System.err.println("indiv are "+ldc.indiv());
						if(k1<toComp.length) toComp[k1][km] = ldc;
						if(non_neg_ind<0 && ldc!=null){
							non_neg_ind = k1;
						}
						allD[k1][km] = ldc;
						
						indiv_[k1][km] = allD[k1][km].indiv();
						if(indiv_[k1][km].size()==0){
							throw new RuntimeException("should not be zero");
						}
						if (km==0 && nonnull < 0 && indiv_[k1] != null)
							nonnull = k1;
						ldc.name = convert(all.get(k1));
						}
						
					} catch (Exception exc) {
						exc.printStackTrace();
					}
				}
				}
				
				indiv = new HashSet<String>(indiv_[nonnull][0]);
				for (int i1 = 0; i1 < all.size(); i1++) {
				
					for(int kk=0; kk<indiv_[i1].length; kk++){
					if (indiv_[i1][kk] != null && indiv_[i1][kk].size()>0)
					//	System.err.println(indiv_[i1][kk]);
						//if(indiv_[i1][kk].size()>0){// throw new RuntimeException("!! size is 0 at "+allD[i1][kk].name+" "+chr[kk]);
						indiv.retainAll(indiv_[i1][kk]);
						//}
					}
				}
				if(indiv.size()==0){
					throw new RuntimeException("size is 0");
				}
				/*if (toLimitIndiv.length > 0) {
					indiv.retainAll(Arrays.asList(toLimitIndiv));
				}*/
				for(int km =0; km<chr.length; km++){
					for (int k1 = 0; k1 < all.size(); k1++) {
						if(allD[k1][km]!=null) allD[k1][km].restricToAlias(indiv);
					}
				//}
					for(int i=0; i<str_base.size(); i++){
						baseComp[i] = allD[alias_base[i]];
					}
				}
			File compa = new File(results1, "comparison_");
			compa.mkdir();
			XYSeriesCollection[][][][] dc = new XYSeriesCollection[nocop.length][ROC.mapLength][2][];
			int[][] total = new int[nocop.length][ROC.mapLength];
			cp = new ChartPanel[nocop.length][ROC.mapLength][];
			// cp1 = new
			// ChartPanel[probsets.length][base_indices.size()][nocop.length][2];
			String[] type=null;
	
		//	for (int i = 0; i < probsets.length; i++) {
				List<Inner>[][] l = new ArrayList[nocop.length][chr.length];
				for(int k=0; k<l.length; k++){
					for(int kk=0; kk<l[k].length; kk++){
						l[k][kk] = new ArrayList<Inner>();
					}
				}
				for(int k=0; k<str_comp_alias.size(); k++){
					for(int k1=0; k1<toComp[k].length; k1++){
						toComp[k][k1].name = str_comp_alias.get(k);
					}
				}
				ROC.Inner[][][] dc1 = new ROC.Inner[nocop.length][toComp.length][chr.length];
				PrintWriter[][] out_all;// = new PrintWriter[base_indices.size()][][];
				PrintWriter[]out_conf;// = new PrintWriter[base_indices.size()][];
				PrintWriter[][]out_breakL;//= new PrintWriter[base_indices.size()][][];;
				PrintWriter[][] out_breakR;//= new PrintWriter[base_indices.size()][][];
				PrintWriter[]out_max;//= new PrintWriter[base_indices.size()][];
				PrintWriter[]out_ld;// = new PrintWriter[base_indices.size()][];
				//inner1: for (int k1 = 0; k1 < base_indices.size(); k1++)
				{
					//if (baseComp[k1].loc.size() == 0)
				//		continue inner1;
					// File[] toComp =new File[f.length];
					File outpDir = new File(compa, probsets1);
					if (!outpDir.exists())
						outpDir.mkdir();

					// toComp[i] = new DataCollection[result_set.length];
				
					 out_all = getPw1(nocop, outpDir, "roc_base_"
							+ ("all_") + "_index_"
							+  "_line"
							 + ".txt", allD[0][0].indiv(), append);
				 out_conf = getPw(nocop, outpDir, "out_conf_"+
							 "_index_"
						 + "_line_"
							 + ".txt", append);
				 out_ld = getPw(nocop, outpDir, "out_ld_"+
						 "_index_"
					+ "_line_"
						 + ".txt", append);
					 out_breakL =getPw1(nocop, outpDir, "roc_base_"
							+ ("break") + "_index_"
							+  "_line"
							 + ".txt", allD[0][0].indiv(), append);
					 out_breakR =getPw1(nocop, outpDir, "roc_base_"
								+ ("break") + "_index_"
								+  "_line"
								 + ".txt", allD[0][0].indiv(), append);
					 out_max = getPw(nocop, outpDir, "roc_base_"
							+ ( "max_" ) + "_index_"
							+  "_line"
							 + ".txt", append);	
					
					for(int kl=0; kl<chr.length; kl++){
					inner: for (int k = 0; k < toComp.length; k++) {
						
						
						// if(k == base_indices.get(k1)) continue inner;
						if (toComp[k] == null)
							continue inner;
						try {
							DataCollection base =  baseComp[k][kl];
							if(base!=null){
							ROC roc = new ROC(baseComp[k][kl],
									outpDir, chr[kl], centromere[kl], minCount(k), minLength(k));
							try {
								roc.setComp(toComp[k][kl]);

								for (int i1 = 0; i1 < nocop.length; i1++) {
									
									ROC.Inner inner = roc.run(str_comp_alias.get(k),nocop[i1],
											roc.toComp.name.toLowerCase()
													.indexOf("penn") < 0, out_all[i1], out_breakL[i1], out_breakR[i1],out_max[i1],
													out_ld[i1], out_conf[i1]);// hwe_correct[kk].get(k));
									//inner.name = str_comp_alias.get(kl);
									dc1[i1][k][kl] = inner;
									print(out_all[i1],(inner.name+"\t"));
									out_max[i1].print(inner.name+"\t");
									out_ld[i1].print(inner.name+"\t");
									print(out_breakR[i1],inner.name+"\t");
									print(out_breakL[i1],inner.name+"\t");
									l[i1][kl].add(inner);
									if(type==null) type = inner.getTypes();
								}
							} catch (Exception exc) {
								exc.printStackTrace();
							}
							// roc.aucpw.close()
							}
						} catch (Exception exc) {
							exc.printStackTrace();
						}
					}
					

				}
				}
				String toP = "posneg\tscore\tpos\tid\trsid\n";
				// for(int k1=0; k1<out_all.length; k1++){
				 for(int k=0; k<out_all.length; k++){
						print(out_all[k],toP);
						out_max[k].print(toP);
						out_max[k].flush();
						out_ld[k].print(toP);
						out_max[k].flush();
						print(out_breakL[k],toP);
						print(out_breakR[k],toP);
					}
				// }
				//System.err.println( baseComp.length);
				
				File outpDir = new File(compa, probsets1);
			
			
			
				
				PrintWriter[] exclude = new PrintWriter[allD.length];//
				for(int k=0; k<exclude.length; k++){
					exclude[k] = new PrintWriter(new File(outpDir, "exclude_"+allD[k][0].name));
				}
				String baseName = allD[non_neg_ind][0].name;
				outer: for(int kl=0; kl<chr.length; kl++){
					if(allD[non_neg_ind][kl]==null) continue outer;
				int sze = Math.min(maxcount,allD[non_neg_ind][kl].loc.size());
				
				List<String> snpid = allD[non_neg_ind][kl].snpid;
			//	PrintWriter[] include = new PrintWriter[new PrintWriter(new File(outpDir, "include_"));
			
				inner: for (int ik = -buff; ik < sze; ik++) {
					
					
					shuffle(rsbuff,ik+buff<sze ? snpid.get(ik+buff) : null);
					//List<String> rs = snpid.subList(Math.max(0, ik-buff), Math.min(ik+buff, snpid.size()));// it.next();
					int loc = ik>=0 ? allD[non_neg_ind][kl].loc.get(ik) : -1;
					// String chrom = baseComp[i][0].chrom;
					int pos = ik>=0 ? allD[non_neg_ind][kl].loc.get(ik) : -1;// baseComp[i][0].loc.get(baseComp[i][0].snpid.indexOf(rs));
					int posl =ik-1>=0 ?  allD[non_neg_ind][kl].loc.get(ik-1) : -1;
					int posr = ik+1<sze && ik+1>=0 ?   allD[non_neg_ind][kl].loc.get(ik+1) : -1;
					// boolean all = true;
					boolean excluded = false;
					String rs  = rsbuff[rsbuff.length-1];
					
					/*for (int k = 0; k < allD.length; k++) {
						if (allD[k][kl] != null
								&& !allD[k][kl].canProcess(rs)){
							
						exclude[k].println(rs+" "+allD[k][kl].name+" "+rs);
							excluded = true;
						}
						
					}
					shuffle(canProcess,excluded );*/
				//	System.err.println("updating "+rs);
					for (int k = 0; k < allD.length; k++) {
						if(!allD[k][kl].updateIndex(rs) || (rsmid!=null && ! rs.endsWith(rsmid))){
							if(rs!=null && rsmid!=null && ! rs.endsWith(rsmid)){
							exclude[k].println(rs+" "+allD[k][kl].name+" "+pos+" was empty");
							}
							excluded = true;
						}
					}
					shuffle(canProcess,excluded);
					//if(excluded) continue inner;
					
					
					int cnt=0;
					if(ik>=0 && !canProcess[buff]){
					for(int kkk=0; kkk<l.length; kkk++){
					Iterator<Inner> it1 = l[kkk][kl].iterator();
					if(it1.hasNext()){
						
						for (;it1.hasNext();cnt++) {
							Inner nxt = it1.next();
							boolean  last = !it1.hasNext();
							nxt.run(rsbuff[buff], pos, posl, posr,  last);
						}
						
					}
					}
					}
					/*for(int k=0; k<out_all.length; k++){
						for(int k1=0; k1<out_all[k].length; k1++){
						out_all[k][k1].println("posneg\t"+pos+"\t"+id);
						out_max[k][k1].println("posneg\tpos\tid");
						out_break[k][k1].println("posneg\tpos\tid");
						}
					}*/
				}
				}
				for(int ik=0; ik<exclude.length; ik++){
					exclude[ik].close();//include.close();
				}
				
				for(int kkk=0; kkk<l.length; kkk++){
					for(int kl=0;kl<l[kkk].length; kl++){
				for (Iterator<Inner> it1 = l[kkk][kl].iterator(); it1.hasNext();) {
					it1.next().finish();
				}
					}
				}
				for(int k=0; k<out_all.length; k++){
					//for(int k1=0; k1<out_all[k].length; k1++){
					close(out_all[k]);
					out_max[k].close();
					close(out_breakL[k]);
					close(out_breakR[k]);
					//
					out_ld[k].close();
					//}
				}
			//	for (int i2 = 0; i2 < dc.length; i2++)
				{
					for (int i3 = 0; i3 < dc.length; i3++) {
						ROC.Inner[][] inner = dc1[i3];
						ROC.Inner[] inner1 = new ROC.Inner[inner.length];
						for(int k=0; k<inner.length; k++){
							inner1[k] = inner[k][0];
							for(int j=1; j<inner[k].length; j++){
								if(inner[k][j]!=null) inner1[k].add(inner[k][j]);
							}
							inner1[k].printConf(inner1[k].pw_conf);
							inner1[k].pw_conf.flush();
						}
						
						// if(inner!=null) {
						
						
							///series[2][0] 
							          XYSeriesCollection[] conc=CompareBasic.getSeriesConc(inner1,baseName);
							//series[3][0] 
							          XYSeriesCollection[] cum =CompareBasic.getSeriesCum(inner1,baseName);
							        //  XYSeriesCollection[] conc1=CompareBasic.getSeriesHist(inner1,baseName);
					     for(int kkk=0; kkk< dc[i3].length; kkk++){
					    	 dc[i3][kkk] = getSeries(inner1, kkk, baseName);
					    	 dc[i3][kkk][3] = conc;
					    	 dc[i3][kkk][4] = cum;
					    	 int tot=0;
					    	 for(int kl=0; kl<inner[0].length; kl++){
					    		 if(inner[0][kl]!=null)
					    		 tot += inner[0][kl].total;
					    	 }
					    	 total[i3][kkk] = tot;
					    	
					     }
					    	
						
						
						// }

					}
				}
				for(int k=0; k<out_all.length; k++){
					out_conf[k].close();
				}
	//		}
		
			String[][] titles = ROC.getSeriesTitles();
			double[][] x_max = ROC.x_max();
			boolean[] dots = ROC.dots();
			boolean[][] log = ROC.getLogAxes();
			
			//for (int i1 = 0; i1 < dc.length; i1++) {
				//for (int i2 = 0; i2 < dc.length; i2++)
				{
					for (int i3 = 0; i3 < dc.length; i3++) {
						for (int i4 = 0; i4 < dc[i3].length; i4++) {
							cp[i3][i4] = new ChartPanel[titles.length];
							for (int kk = 0; kk < titles.length; kk++) {
								cp[i3][i4][kk] = 
									ROC.getChartPanel(
										dc[i3][i4][kk],
										result_set_color, result_set_stroke,
										probsets1, convert(str_comp.get(0)),
												
										nocop[i3], i4,
										total[i3][i4], chr.length>=22 ? "all" : Arrays.asList(chr).toString(),
										titles[kk][1], titles[kk][2],
										titles[kk][0], log[kk][0], log[kk][1],false, type[i4], x_max[kk], dots[kk]);
							}
							/*
							 * cp[i1][i2][i3][i4][1] =
							 * ROC.getChartPanel(dc[i1][i2][i3][i4][1],
							 * result_set_color, result_set_stroke,
							 * probsets[i1],
							 * convert(result_set_name.get(base_indices
							 * .get(i2))),nocop[i3], i4==1
							 * ,total[i1][i2][i3][i4],chr, "Prob", "Accuracy",
							 * "Accuracy: "); if(dc[i1][i2][i3][i4].length>2){
							 * cp[i1][i2][i3][i4][2] =
							 * ROC.getChartPanel(dc[i1][i2][i3][i4][2],
							 * result_set_color, result_set_stroke,
							 * probsets[i1],
							 * convert(result_set_name.get(base_indices
							 * .get(i2))),nocop[i3], i4==1
							 * ,total[i1][i2][i3][i4],chr, "Prob",
							 * "True Positive", "True Positive: ");
							 * cp[i1][i2][i3][i4][3] =
							 * ROC.getChartPanel(dc[i1][i2][i3][i4][2],
							 * result_set_color, result_set_stroke,
							 * probsets[i1],
							 * convert(result_set_name.get(base_indices
							 * .get(i2))),nocop[i3], i4==1
							 * ,total[i1][i2][i3][i4],chr, "Prob",
							 * "False Positive", "False Positive: "); }
							 */

						}

					}

				}
			//}
			boolean inclMax = false;
			jpane = new JTabbedPane[cp[0].length]; 
			for(int i=0; i<jpane.length; i++){
				jpane[i] = new JTabbedPane();
			}
		
			for(int l1 =0; l1<jpane.length; l1++){
			jpane[l1].setName(chr.length>=22 ? "all" : Arrays.asList(chr).toString());
			// String[][] titles = ROC.getSeriesTitles();
			for (int i1 = 0; i1 < titles.length; i1++) {
				JPanel tab = new JPanel();
				// JComponent jin1 = new JScrollPane(tab,
				// JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				// JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
				jpane[l1].add(tab, titles[i1][0]);

				tab.setLayout(new BoxLayout(tab, BoxLayout.Y_AXIS));

				//for (int i2 = 0; i2 < probsets.length; i2++) {
					JPanel[] tab1 = new JPanel[inclMax ? 2 : 1];
					// Legend leg=new Legend();
					for (int ii = 0; ii < tab1.length; ii++) {
						tab1[ii] = new JPanel();
						tab1[ii].setLayout(new BoxLayout(tab1[ii],
								BoxLayout.X_AXIS));
						tab.add(tab1[ii]);
					}

					// JPanel tab2 = new JPanel();
					// tab2.setLayout(new BoxLayout(tab2, BoxLayout.Y_AXIS));
					// tab.add(tab2);

					Dimension d = null, d1 = null, d2 = null;
					//for (int i3 = 0; i3 < base_indices.size(); i3++) 
					{
						for (int i4 = 0; i4 < nocop.length; i4++) {
							if (cp[i4][l1][i1] != null) {
								// JPanel tab12 = new JPanel();
								// tab12.setLayout(new BoxLayout(tab12,
								// BoxLayout.X_AXIS));
								// tab1[i4].add(tab12);
								ChartPanel c_i = cp[i4][l1][i1];
								d = new Dimension(2 * c_i.getSize().width, c_i
										.getSize().height);
								d1 = new Dimension(2 * c_i.getSize().width,
										2 * c_i.getSize().height);
								d2 = new Dimension(2 * c_i.getSize().width, 20);
								/*
								 * if(i3==0 && i4==0) { LegendTitle legend = new
								 * LegendTitle(c_i.getChart().getXYPlot());
								 * legend.setMargin(new RectangleInsets(1.0,
								 * 1.0, 1.0, 1.0)); // legend.setFrame(new
								 * LineBorder());
								 * legend.setBackgroundPaint(Color.GREEN);
								 * //legend.setPosition(RectangleEdge.BOTTOM);
								 * leg.set(legend);
								 * 
								 * leg.setSize(d1); leg.setMinimumSize(d1);
								 * 
								 * }
								 */

								tab1[0].add(c_i);

								tab1[0].setSize(d);
								if (tab1.length > 1) {
									tab1[1].add(cp[i4][1][i1]);
									tab1[1].setMinimumSize(d);
								}
								if(!ROC.countRegs)redrawAxes(tab1);
								tab.setSize(d1);
								tab.setMinimumSize(d1);
								// tab2.setSize(d1);
							}

						}
					}
					// tab.add(leg);
					// jpane.add(leg);

			//	}
			}
			
			if(headless==null || !headless.equals("true")){
			JFrame jf = new JFrame(chr.length>=22 ? "all" : Arrays.asList(chr).toString());
			jf.setName(chr.length>=22 ? "all" : Arrays.asList(chr).toString());
			jf.setJMenuBar(getJMenuBar());
			jf.getContentPane().add(jpane[l1]);

			jf.pack();
			jf.setVisible(true);
			}
			else{
				Headless jf;
				  jf = new Headless(jpane[l1]);
					jf.pack();
					jf.setVisible(true);
			}
			File dirOut = new File(results1,"graphs_"+args1[0+offset]+"_"+args1[3+offset]+
					"_"+(args1[2+offset].replaceAll(";", "").replaceAll(":", "")));
			dirOut.mkdir();
			writeToFile(dirOut,l1);
			}
			es.shutdown();
			// }
		} catch (Exception exc) {
			exc.printStackTrace();
		}
	}

	private static int minCount(int k) {
	if(k<minCount.length) return minCount[k];
	else return minCount[0];
}
	
	private static int minLength(int k) {
		if(k<minLength.length) return minLength[k];
		else return minLength[0];
	}

	private static int[] convert(String[] split) {
	int[] res = new int[split.length];
	for(int i=0; i<res.length; i++){
		res[i] = Integer.parseInt(split[i]);
	}
	return res;
}

	private static List<String>  readBuild(File buildF2, String chr, int st, int end) {
	try{
		List<String> l = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(buildF2));
		String st1 = "";
		while((st1 = br.readLine())!=null){
			if(st1.startsWith(chr)){
				
				String[] str = st1.split("\\s+");
				if(str[0].equals(chr)){
				int start = Integer.parseInt(str[1]);
				if(start >= st && start <=end){
					l.add(str[3]);
				}
				}
			}
		}
	return l;
	}catch(Exception exc){
		exc.printStackTrace();
	}
	return null;
}

	private static void shuffle(String[] rsbuff2, String snpid) {
	for(int i=1; i<rsbuff2.length; i++){
		rsbuff2[i-1] = rsbuff2[i];
	}
	rsbuff2[rsbuff2.length-1] = snpid;
	
}
	
	private static void shuffle(boolean[] rsbuff2, boolean snpid) {
		for(int i=1; i<rsbuff2.length; i++){
			rsbuff2[i-1] = rsbuff2[i];
		}
		rsbuff2[rsbuff2.length-1] = snpid;
		
	}

	private static DataCollection[] trans(DataCollection[][] dcin) {
		DataCollection[] dc = new DataCollection[dcin.length * dcin[0].length];
		for(int k=0; k<dcin.length; k++){
			for(int j=0; j<dcin[0].length; j++){
				dc[k*dcin[0].length+j] = dcin[k][j];
			}
		}
		return dc;
	}

	public static void redrawAxes(JPanel[] tab1) {
		List<ValueAxis> x = new  ArrayList<ValueAxis>();
		List<ValueAxis> y = new  ArrayList<ValueAxis>();
		double l_y = Double.POSITIVE_INFINITY;
		double l_x = Double.POSITIVE_INFINITY;
		double u_y = Double.NEGATIVE_INFINITY;
		double u_x = Double.NEGATIVE_INFINITY;
		for(int k=0; k<tab1.length; k++){
			for(int l=0; l<tab1[k].getComponentCount(); l++){
				
				ChartPanel cp = (ChartPanel)tab1[k].getComponent(l);
				ValueAxis x_i =cp.getChart().getXYPlot().getDomainAxis();
				ValueAxis y_i =cp.getChart().getXYPlot().getRangeAxis(); 
				x.add(x_i);
				y.add(y_i);
				if(x_i.getRange().getLowerBound() < l_x) l_x = x_i.getRange().getLowerBound();
				if(x_i.getRange().getUpperBound() >u_x) u_x = x_i.getRange().getUpperBound();
				if(y_i.getRange().getLowerBound() < l_y) l_y = y_i.getRange().getLowerBound();
				if(y_i.getRange().getUpperBound() >u_y) u_y = y_i.getRange().getUpperBound();
			}
		}
		if(x.size()>1 && l_x<u_x && l_y < u_y){
		Range xr = new Range(l_x, u_x);
		Range yr = new Range(l_y, u_y);
		for(int i=0; i<x.size(); i++){
			x.get(i).setRange(xr);
			y.get(i).setRange(yr);
		}
		}
		
	}

	private static void close(PrintWriter[] p) {
		for(int i=0; i<p.length; i++){
			p[i].close();
		}
		
	}

	private static void print(PrintWriter[] printWriters, String string) {
		for(int i=0; i<printWriters.length; i++){
			printWriters[i].print(string);
		}
		
	}

	private static PrintWriter[] getPw(Integer[][] nocop, File outpDir, String string, boolean append) throws Exception {
		PrintWriter[] pw = new PrintWriter[nocop.length];
		for(int i=0; i<pw.length; i++){
			pw[i] =  new PrintWriter((new FileWriter(new File(outpDir,ROC.getCNString(nocop[i])+"_"+ string),append)));
		}
		return pw;
	}
	
	private static PrintWriter[][] getPw1(Integer[][] nocop, File outpDir, String string, List<String>ids, boolean append) throws Exception {
		PrintWriter[][] pw = new PrintWriter[nocop.length][ids.size()];
		for(int i=0; i<pw.length; i++){
			File f = new File(outpDir,ROC.getCNString(nocop[i])+"_"+ string);
			f.mkdir();
			for(int j=0; j<ids.size(); j++){
				pw[i][j] =  new PrintWriter(new BufferedWriter((new FileWriter(new File(f,ids.get(j)),append))));
			}
			
		}
		return pw;
	}

	static class Legend extends JPanel {
		LegendTitle legend;

		Legend() {

		}

		public void set(LegendTitle legend2) {
			this.legend = legend2;

		}

		@Override
		public void paint(Graphics g) {
			double ww = this.getWidth();
			double hh = this.getHeight();
			Graphics2D g2 = (Graphics2D) g;
			RectangleConstraint constraint = new RectangleConstraint(ww,
					new Range(0.0, ww), LengthConstraintType.RANGE, hh,
					new Range(0.0, hh), LengthConstraintType.RANGE);
			Size2D d2 = legend.arrange(g2, constraint);

			legend.draw(g2, this.getBounds());
		}
	}

	static ExecutorService es = Executor.getEs(CompareBasic.class, 1);
/*
private static XYSeriesCollection getSeriesCorr(Inner[] inner, String baseName){
	XYSeriesCollection xys = new XYSeriesCollection();
	double min =Double.MAX_VALUE;
	double max = Double.MIN_VALUE;
	for(int i=0; i<inner.length; i++){
		if(inner[i].name.equals(baseName)) continue;
	
		for(int j=i+1; j<inner.length; j++){
			if(inner[i].name.indexOf("Penn")<0 &&inner[j].name.indexOf("Penn")<0 ) continue;
			//if(inner[j].name.equals(baseName)) continue;
			XYSeries[] corr = ROC.getQQ(inner[i], inner[j]);
			for(int k1=0; k1<corr.length; k1++){
				xys.addSeries(corr[k1]);
			for(int k=0; k<corr[k1].getItemCount(); k++){
				double x = corr[k1].getX(k).doubleValue();
				double y = corr[k1].getY(k).doubleValue();
				if(Math.max(x, y)>max) max = Math.max(x, y);
				if(Math.min(x, y)<min) min = Math.min(x, y);
			}
			
			}
			
		}
	}
	XYSeries line = new XYSeries("line");
	xys.addSeries(line);
	line.add(min, min);
	line.add(max,max);
	return xys;
}*/
private static XYSeriesCollection[] getSeriesCum(Inner[] inner, String baseName){
	XYSeries[][] corr1 = new XYSeries[inner.length][];
	for(int i=0; i<inner.length; i++){
	
			corr1[i] =inner[i].getCum();
	}
	return getCollection(corr1);
	
	
}

private static XYSeriesCollection[] getSeriesConc(Inner[] inner, String baseName){
	
	XYSeries[][] corr1 = new XYSeries[inner.length][];
	for(int i=0; i<inner.length; i++){
	
			corr1[i] =inner[i].getConc();
	}
	return   getCollection(corr1);
	
}

static XYSeriesCollection[] getCollection(XYSeries[][]corr1){
	XYSeriesCollection []xys = new XYSeriesCollection[corr1[0].length];
	for(int i=0; i<xys.length; i++){
		xys[i] = new XYSeriesCollection();
	}
	
		
		for(int i=0; i<corr1.length; i++){
			boolean skip = false;
			for(int k=0; k<corr1[0].length; k++){
				if(corr1[i][k]!=null && corr1[i][k].getItemCount()>0) skip = false;
			}
			if(!skip){
			for(int k=0; k<corr1[0].length; k++){
			if(corr1[i][k]!=null){
				xys[k].addSeries(corr1[i][k]);
			}
			}
		}
	}
	return xys;
}
	
	private static XYSeriesCollection[][] getSeries(Inner[] inner,
			final int type, String baseName) throws Exception {
		String[][] titles = ROC.getSeriesTitles();
		int[] cnts = ROC.getTypes();
     
		XYSeriesCollection[][] series = new XYSeriesCollection[titles.length][];
		for (int k = 0; k < cnts.length; k++) {
			series[k] = new XYSeriesCollection[cnts[k]];
		}
		for (int i = 0; i < series.length; i++) {
			for (int k = 0; k < series[i].length; k++) {
				series[i][k] = new XYSeriesCollection();
			}
		}
		List<XYSeries[][]> l = new ArrayList<XYSeries[][]>();

		for (int i = 0; i < inner.length; i++) {

			final Inner in = inner[i];
			if (inner[i] != null)
				try {
					inner[i].check(inner[0]);
				} catch (Exception exc) {
					System.err.println("prob with ");
					exc.printStackTrace();
					System.exit(0);
				}
			Callable call = new Callable() {

				public Object call() throws Exception {
					if (in != null) {
						// XYSeries[][] ser = in.getSeries(b);
						XYSeries[][] ser1 = in.getSeries_(type);
						return ser1;
						// XYSeries[] sum = new XYSeries[]{ser[0][0],ser1[0],
						// ser1[1]};
						/*
						 * for(int i=0; i<sum.length; i++){ sum[i] = ser[i][0];
						 * }
						 */
						// return ser1;
					} else
						return null;
				}

			};
			l.add((XYSeries[][])call.call());

		}
		/*for(int i=0; i<l.size(); i++){
			l.get(i).call();
		}
		List<Future> fut = es.invokeAll(l);*/
		for (int i = 0; i < inner.length; i++) {
			XYSeries[][] ser = l.get(i);
			//Object obj = fut.get(i).get();
			if (ser != null) {
				for (int i1 = 0; i1 < series.length; i1++) {
					//XYSeries[][] ser = ((XYSeries[][]) l.get(i));
					// XYDataItem di1 = ser[0][0].getDataItem(0);
					// XYDataItem di2 = ser[0][0].getDataItem(1);
					// System.err.println(di1.getXValue()+" "+di2.getXValue()+" "+di1.getYValue()+" "+di2.getYValue());
					if (ser[i1] != null) {// ser[0][0].getDataItem(0).getYValue()<0.999
											// || true){
						for (int k = 0; k < ser[i1].length; k++) {
							if (ser[i1][k] != null)
								series[i1][k].addSeries(ser[i1][k]);
						}
					}
				}

			}

		}
		 if(true){
	    	  XYSeriesCollection maxLenDists = ROC.getLenghtDists(inner, type);
	    	  series[2] =new XYSeriesCollection[] { maxLenDists};
	      }
		return series;
	}

	static Map<String, String> names = new HashMap<String, String>();

/*	private static Color convert1(String name) {

		try {
		  int i =    expt.indexOf(name);
		if (i<0) return Color.black;
			//if (nme.length == 3) {
				Color nme1 = (Color) Color.class.getField(
						col.get(i)).get(null);

				return nme1;
			//}
		} catch (Exception exc) {
		}
		;
		return Color.black;
	}*/

	private static String convert(String name) {
		/*if (names.containsKey(name))
			return names.get(name);
		if (true) {
			// if(names.)
			return name.split("-")[0];
		}
		names.put(name, name);*/
		return name;
	}

	static ChartPanel[][][] cp;

	// static ChartPanel[][][][] cp1;
	private static JMenuBar getJMenuBar() {
		JMenuBar jmb = new JMenuBar();
		JMenu opts = new JMenu("Options");
		jmb.add(opts);
		opts.setEnabled(true);
		opts.add(print_item = new MyMenuItem("Save picture"));
		return jmb;
	}

	static MyMenuItem print_item;

	/*public static float[] generate(int i) {
		float[] res = new float[i];
		for (int k = 0; k < res.length; k++) {
			res[k] = Constants.nextInt(20);
		}
		return res;
	}*/

	static JFileChooser filechooser = new JFileChooser();

	public static class MyMenuItem extends JMenuItem implements ActionListener {

		public MyMenuItem(String string) {
			super(string);
			this.addActionListener(this);
		}

		public void actionPerformed(ActionEvent e) {
			filechooser.setCurrentDirectory(new File(System
					.getProperty("user.dir")));
			filechooser.setMultiSelectionEnabled(false);
			filechooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			filechooser.showOpenDialog(null);
			File outFile = filechooser.getSelectedFile();
			try {
				writeToFile(outFile,0);
				writeToFile(outFile, 1);
			} catch (Exception exc) {
				exc.printStackTrace();
			}
		}

	}

	public static void writeToFile(File out, int l1) throws Exception {
		boolean emf = false;
		/*
		 * if(out.getName().endsWith(".emf")) emf = true; else if
		 * (!out.getName().endsWith(".pdf")){ out = new File(out.getParent(),
		 * out.getName()+".pdf"); }
		 */
		// (Iterator<String> phen = pheno.iterator(); phen.hasNext();i++){
		// String phent = phen.next();
		out.mkdir();

		for (int i3 = 0; i3 < jpane[l1].getTabCount(); i3++) { // types
			File outD1 = new File(out, probsets1);
			outD1.mkdir();
			File outD = new File(outD1, jpane[l1].getTitleAt(i3).replaceAll("\\s+",
					""));
			JPanel cp_i = (JPanel) (jpane[l1].getComponent(i3));// .getComponent(0);
			outD.mkdir();
			boolean png = true;
			File out1 = new File(outD, png ? "graph"+l1+".png" : "graph"+l1+".pdf");
			File out2 = new File(outD,"graph"+l1);
			out2.mkdir();
			Properties p = new Properties();
			p.setProperty("PageSize", "A4");
			// ExportDialog ep;
			p.setProperty(PageConstants.ORIENTATION, PageConstants.LANDSCAPE);
			ImageGraphics2D g = new ImageGraphics2D(out1, cp_i,
					ImageConstants.PNG);
			// ImageGraphics2D g = emf ? new EMFGraphics2D(out,cp_i.getSize()) :
			// new PDFGraphics2D(out1,
			// cp_i.getSize());
			g.setProperties(p);
			g.startExport();
			cp_i.print(g);
			g.endExport();
			for (int i = 0; i < cp_i.getComponentCount(); i++) {
				JPanel jp1 = (JPanel) cp_i.getComponent(i);
				for (int k = 0; k < jp1.getComponentCount(); k++) {
					JFreeChart jf = ((ChartPanel) jp1.getComponent(k))
							.getChart();
					File nf = out2;//new File(out2, jf.getTitle().getText().replace(
							//':', ' ').replaceAll("\\s+", ""));
					//nf.mkdir();
					for(int jj=0; jj<jf.getXYPlot().getDatasetCount();jj++){
					XYDataset xyp = jf.getXYPlot().getDataset(jj);
					for (int kk = 0; kk < xyp.getSeriesCount(); kk++) {
						PrintWriter pw = new PrintWriter(new FileWriter(
								new File(nf, xyp.getSeriesKey(kk).toString()+"_"+i+"_"+k+"_"+jj+"_"+kk)));
						for (int kkk = 0; kkk < xyp.getItemCount(kk); kkk++) {
							Number x = xyp.getX(kk, kkk);
							Number y = xyp.getY(kk, kkk);
							pw.println(x + "\t" + y);
						}
						pw.close();
					}
					}
				}
			}

		}
	}

	private static float[] parse(String[] strings) {
		float[] res = new float[strings.length];
		for (int i = 0; i < res.length; i++) {
			res[i] = Float.parseFloat(strings[i]);
		}
		return res;
	}

	private static int whichIndex(String[] result_set_name, String string) {
		for (int i = 0; i < result_set_name.length; i++) {
			if (result_set_name[i].startsWith(string))
				return i;
		}
		return -1;
	}
}
