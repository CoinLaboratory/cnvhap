package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.states.State;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class Constants {
	public static double samplePermute = 0;
	
	public static Boolean[] cnvP = new Boolean[]{null};
	public static Boolean cnvP(int i){
		if(i<cnvP.length) return cnvP[i];
		else return cnvP[0];
	}
	
	
//public static boolean modifyExponent =false;
//	public static double[] priorMod = null;
	static abstract class AbstractDataParams{
		//String[] phenosToInclude = null;
		String[] outputRow;
		String[] firstRow;
		int[][] inputToOutput;
		String fieldToSplit = "phenoToInclude";
		
		void initialise(String[] input) throws Exception{
			inputToOutput = new int[input.length][];
			int cnt=0;
			List<String> incl = new ArrayList(Arrays.asList(Constants.include));
			String[] inputDir = this.getRow_(0);
			List<String> firstRow = new ArrayList<String>();
			for(int k=0; k<input.length; k++){
				String[] fields = input[k].split(";");
				inputToOutput[k] = new int[fields.length];
				for(int j=0; j<fields.length; j++){
				  inputToOutput[k][j] = cnt;
				  cnt++;
				  String inputId = j==0 ? inputDir[k] : inputDir[k]+"_"+j;
				  firstRow.add(inputId);
				  if(j>0 && incl.contains(inputDir[k])){
					  incl.add(inputId);
				  }
				}
			}
			Constants.include = incl.toArray(new String[0]);
			outputRow  = new String[cnt];
			this.firstRow = firstRow.toArray(new String[0]);
		}
		
		void initialise() throws Exception{
			String[] col = getColumn(0);
			for(int i=0; i<col.length; i++){
				if(col[i].equals(fieldToSplit)){
					String[] input = getRow_(i);
					initialise(input);
					return;
				}
			}
			initialise(getRow_(0));
		//	System.err.println("not found "+fieldToSplit);
		}
		
		abstract int getRows();

		abstract String getCell(int i, int i2) throws Exception;

		abstract String[] getRow_(int inputDir_index) throws Exception;

		abstract String[] getColumn(int i) throws Exception;

		public String[] getRow(int i)  throws Exception{
			if(i==0) return firstRow;
			String[] row_ = getRow_(i);
			outputRow[0] = row_[0];
			for(int k=1; k<inputToOutput.length; k++){
				int[] inds = inputToOutput[k];
				for(int j=0; j<inds.length; j++){
					if(row_[0].equals(fieldToSplit)){
						outputRow[inds[j]] = row_[k].split(";")[j];
					}else{
						//System.err.println("ROW " + row_[k]);
						outputRow[inds[j]] =  row_[k];
					}
				}
			}
			
			// TODO Auto-generated method stub
			return outputRow;
		}
		
	}
	public static String[] rplot_names = new String[] {"all"};
	
	public static String[] rplot_names(){
		return rplot_names;
	}
	
	
	static class CSVDataParams extends AbstractDataParams{
		File ff;
		
		int ind=-1;
		final String sep = ",";
		
		final boolean transpose;
		List<String[]> res = new ArrayList<String[]>();
	public CSVDataParams(File ff) throws Exception {
		this.ff =ff;
		BufferedReader br = new BufferedReader(new FileReader(ff));
		String st =br.readLine();
		String[]str = st.replaceAll("\"", "").split(sep);
		transpose = str[1].equals("format") || str[1].equals("inputDirLoc");
		res.add(str);
		
		for(; (st=br.readLine())!=null;){
			String[] str1 = st.replaceAll("\"", "").split(sep);
			if(str1.length==0) break;
			res.add(str1);
		}
		initialise();
		
	}
	
	public int getRows() {
		return transpose ? res.get(0).length : res.size() ;
	}
	public String getCell(int i, int i2) throws Exception {
		return !transpose ? res.get(i2)[i2] : res.get(i)[i2];
	}
	public String[] getRow_(int i) throws Exception{
		if(transpose) return getColumnInt(i);
		else return getRowInt(i);
	}
	public String[] getColumn(int i) throws Exception{
		if(!transpose) return getColumnInt(i);
		else return getRowInt(i);
	}
	public String[] getRowInt(int i)  throws Exception{
		return res.get(i);
	}
	
	public String[] getColumnInt(int i) throws Exception {
		String[] res1 = new String[res.size()];
		for(int k=0; k<res1.length; k++){
			res1[k] = res.get(k)[i];
		}
		return res1;
	}
		
	}
	static class DataParams extends AbstractDataParams{
		Sheet[] sh;
		int[] cumlen;
		Cell[][] cells;
		boolean transpose;
	public DataParams(File ff, String[] sheet) throws Exception {
		Workbook wb = 
			Workbook.getWorkbook(ff);
		if(sheet==null || sheet.length==0 || wb.getSheet(sheet[0])==null){
			sheet = wb.getSheetNames();
		}
		sh = new Sheet[sheet.length];
		cumlen = new int[sheet.length];
		int cuml =0;
		cells = new Cell[sheet.length][];
		for(int i=0; i<sheet.length; i++){
			sh[i] = wb.getSheet(sheet[i]);
			if(i==0) transpose = sh[0].getCell(1, 0).getContents().replace("\"","").equals("format");
			int len = (transpose ? sh[i].getRows() : sh[i].getColumns());;
			cuml +=len;
			cumlen[i] = cuml;
		}
		{
		List<String> row = Arrays.asList(this.getColumn(0,0));
		for(int k=1; k<sheet.length; k++){
			List<String> row1 = Arrays.asList(this.getColumn(k, 0));
			if(!row.equals(row1)){
				for(int i=0; i<row.size(); i++){
					int i1 = row1.indexOf(row.get(i));
					if(i1!=i){
						System.err.println("problem "+i+" "+i1+" "+row.get(i));
					}
				}
				throw new RuntimeException("!! first column of each sheet included must be same !! "+sheet[0]+" vs "+sheet[k]);
			}
		}
		}
		initialise();
	}
	int getSheetIndex (int i1){
		for(int i=0; i<cumlen.length; i++){
			if(i1<cumlen[i]){
				return i;
			}
		}
		return -1;
	}
	public String[] getColumn(int k, int i1) throws Exception{
		return trans(transpose ? this.sh[k].getRow(i1) : this.sh[k].getColumn(i1));
	}
	public String[] getColumn(int i1) throws Exception {
		int i = getSheetIndex(i1);
		int i2 = i1 - (i==0 ? 0:  cumlen[i-1]); 
				return getColumn(i,i2);
	}
	public int getRows() {
		// TODO Auto-generated method stub
		return transpose ? sh[0].getColumns() : sh[0].getRows();
	}
	public String getCell(int i, int i2) {
		// TODO Auto-generated method stub
		return transpose ? sh[getSheetIndex(i)].getCell(i2,i).getContents().replace("\"","") : sh[getSheetIndex(i2)].getCell(i,i2).getContents().replace("\"","");
	}
	private String[]  trans(Cell[] cell){
		String[]res = new String[cell.length];
		for(int k=0; k<cell.length; k++){
			res[k] = cell[k].getContents().replace("\"","");
		}
		return res;
	}
	private String[]  trans(Cell[][] cell){
		String[]res = new String[this.cumlen[cumlen.length-1] - (cumlen.length-1)]; //correction factor not to have multiple headers
		int j=0;
		for(int k0=0; k0<cumlen.length; k0++){
			for(int k=k0==0 ? 0 : 1; k<cell[k0].length; k++){
				res[j] = cell[k0][k].getContents().replace("\"","");
				j++;
			}
		}
		return res;
	}
	public String[] getRow_(int inputDir_index) {
		// TODO Auto-generated method stub
		for(int j=0; j<this.cumlen.length; j++){
			this.cells[j] = getRow_(j,inputDir_index);
			
		}
		return this.trans(cells); 
	}
	public Cell[] getRow_(int k, int inputDir_index) {
		// TODO Auto-generated method stub
		return transpose ?sh[k].getColumn(inputDir_index) : sh[k].getRow(inputDir_index);
		
	}
		
	}

	public static String[] sheet = null;

	public static double samplePermute() {
		return samplePermute;
	}
public static boolean penaliseEnd = true;
public static boolean penaliseEnd(){
	return penaliseEnd;
}
	public static boolean savePhasedConfiguration = false;

	public static boolean savePhasedConfiguration() {
		return savePhasedConfiguration;
	}

	// public static String[] chrom = new String[] {""};
	public static String[][] mid = null;// new String[][] {new String[] {"0",
										// (Integer.MAX_VALUE-1)+""}};
	public static String[] chrom = null;

	public static double[] pseudoMod, pseudoMod1;
	
	public static double pseudoMod(int i){
		if(pseudoMod==null) {
			return 1.0;
		}
		return pseudoMod[i];
	}
	public static double pseudoMod1(int i){
		if(pseudoMod1==null) return 1.0;
		if(i>=pseudoMod1.length) return pseudoMod1[0];
		return pseudoMod1[i];
	}
	
	
	
	public static String[][][] regionsToInclude; // regions to include by data
													// set
	public static String[][][] regionsToExclude; // regions to exclude by data
													// set

	public static String[][][][] dataSetRegion(boolean include) {
		return new String[][][][] {include ? regionsToInclude : regionsToExclude};
	}

	public static String strand[] = null;

	public static String strand(int i) {
		if (strand == null || strand[i] == null || strand[i].equals("null"))
			return null;
		return strand[i];
	}

	// public static String[] mid1 = new String[] {"0",
	// Integer.MAX_VALUE-1+""};// target region for LD analysis
	public static int offset = 0;

	
	
	public static double[] expModelIntHotSpot1(int i) {
		double[][] res = expModelIntHotSpot1;
		
		
		if (i >= res.length)
			return null;
		if(i==0){
			double[]res_ = res[i];
			double[] res1 = new double[Constants.modify0.length+1];
			System.arraycopy(res_, 0, res1, 0, Math.min(res_.length, res1.length));
			return res1;
		}
		return res[i];
	}

	public static double[][] expModelIntHotSpot1 = new double[0][]; // modifies
																	// by state

	public static int[] restrict = new int[] { Integer.MAX_VALUE - 10 };
	public static int[] maxIndiv = new int[] { Integer.MAX_VALUE };
	public static String[][] restrictKb = null;
	public static int[] numIt = new int[] { 15 };// 0;//15;
	public static int numRep = 1;
	public static String outputDir = "."; // where program is run (should have
											// param files etc)
	public static String baseDir = "."; // where data is stored
										public static int maxCoordDiff=3*1000*1000;

	public static String outputDir() {
		// String out = outputDir;
		return outputDir;
	}
public static int maxCoordDiff(){
	return maxCoordDiff;
}
	public final static boolean CHECK = false	;
	// final public boolean CHECK1 =false;;
	// public static boolean realRandom =false;
	public static long seed = 23436;
	public static Random rand = new Random(seed);// );
	private static boolean fast = true;
	private static boolean cache = true;
	private static boolean trainWithGenotypes = true;

	// public static int modelCNP = 6;
	public static String[] core;

	public static int[] core() {
		if (core == null) {
			int[] res_ = new int[2];
			int[][] mid = mid();
			int[] rest =  restrictKb(0);
			res_[0] = mid[0][0] -rest[0];
			res_[1] = mid[mid.length - 1][1] + rest[1];
			return res_;
		} else {
			int[] mid_ = new int[2];
			mid_[0] = convert(core[0]);
			mid_[1] = convert(core[1]);
			return mid_;
		}

	}

	public static int[][] mid() {
		String[][] orig = mid;

		if (mid == null)
			return null;
		String[][] mid1;
		if (mid[0].length == 1) {
			mid1 = new String[1][mid.length];

			for (int i = 0; i < mid.length; i++) {
				mid1[0][i] = mid[i][0];
			}

			mid = mid1;
		}
		int index =  mid[0][0].indexOf('_');
		/////neww
		index = -1;
		Logger.global.info("here "+mid[0][0]+" "+index);
		if (index>=0){
			
			String[][] midn = new String[mid[0].length][3];
			for(int k=0; k<mid[0].length; k++){
				midn[k] = mid[0][k].split("\\_");
			}
			mid = midn;
		}
		int[][] mid_ = new int[mid.length][2];
		for (int i = 0; i < mid.length; i++) {
			System.err.println("mid is "+Arrays.asList(mid[i]));
			if(mid[i].length==2 && (mid[i][1].startsWith("p") || mid[i][1].startsWith("q"))){
				try{
					String ip = Constants.inputDir(i);
				int cent = Constants.getKaryoFile(new File(ip));
				if(mid[i][1].equals("p")){
					mid[i] = new String[] {mid[i][0], "0", (cent-10)+""};
					
				}else if(mid[i][1].equals("q")){
					mid[i] = new String[] {mid[i][0], (cent+10)+"", ""+(Integer.MAX_VALUE-1)};
				}else{
					mid[i] = new String[] {mid[i][0],0+"", ""+(Integer.MAX_VALUE-1)};
				}
				}catch(Exception exc){
					exc.printStackTrace();
				}
			}
			mid_[i][0] = convert(mid[i][1]);
			mid_[i][1] = mid[i].length == 2 ? mid_[i][0] : convert(mid[i][2]);
		}
		if (chrom == null) {
			chrom = new String[mid_.length];
			for (int i = 0; i < mid_.length; i++) {
				chrom[i] = mid[i][0];
			}
		}
		return mid_;
	}

	public static Object append(Object left, Object right, Class componentType) {
		// Class componentType = Array.get(left, 0).getClass();
		int left1 = Array.getLength(left);
		int right1 = Array.getLength(right);
		int length = right1 + left1;
		Object res = Array.newInstance(componentType, length);
		System.arraycopy(left, 0, res, 0, left1);
		System.arraycopy(right, 0, res, left1, right1);
		return res;
	}

	public static int convert(String string) {
		String st1 = string.replaceAll(",", "");
		int gb_index = st1.indexOf("gb");
		if (gb_index > 0) {
			return (int) Math.ceil(Double.parseDouble(st1
					.substring(0, gb_index)) * 1e9);
		}
		int mb_index = st1.indexOf("mb");
		if (mb_index > 0) {
			return (int) Math.ceil(Double.parseDouble(st1
					.substring(0, mb_index)) * 1e6);
		}
		int kb_index = st1.indexOf("kb");
		if (kb_index > 0) {
			return (int) Math.ceil(Double.parseDouble(st1
					.substring(0, kb_index)) * 1e3);
		} else
			return Integer.parseInt(st1);
	}

	public static String[] build = new String[] {"build37.txt"};

	public static String build(int i) {
		if(i>=build.length) return build[0];
		return build[i];
	}
public static double[][] transitionMatrix = null;
public static double[][] transitionMatrix(){
	//if(Constants.trainGlobal) return null;
	/*if(!Constants.readGlobalClusterFile() && transitionMatrix==null){
	File clust = new File(Constants.outputDir(), "/clust_in/transitionModel.txt/");
	//System.err.println("clustf "+clust);
File f1 = Constants.getClusterBase(clust);
if(f1.exists())
try{
	BufferedReader br = new BufferedReader(new FileReader(f1));
	String st = "";
	while((st =br.readLine())!=null){
		if(st.indexOf("transitionMatrix")>=0){
			transitionMatrix = Constants.parse_d(st.split("\\s+")[1].split(":"));
		}
	}
}catch(Exception exc){
	exc.printStackTrace();
}
	}*/
	return transitionMatrix;
}
private static void printRec(File dir) {
	System.err.println(dir.getAbsolutePath());
	if(dir.isDirectory()){
	File[] dir1 = dir.listFiles();
	for(int i=0; i<dir1.length; i++){
		printRec(dir1[i]);
	}
	}
	
}

public static File getClusterBase(File clust){
	File clust1=null;
	if(clust.exists()){
		final int[] mid = Constants.mid()[0];
		final String chr = Constants.chrom(0);
		File[] f = clust.listFiles(new FilenameFilter(){

			@Override
			public boolean accept(File pathname, String name) {
				String[] str = name.split("_");try{
				if(str[str.length-3].equals(chr) ){
					
					int start = Constants.convert(str[str.length-2]);
					int end = Constants.convert(str[str.length-1]);
					int overl = Math.min(end - mid[0], mid[1] - start);
					return overl>=0;
					
				}
				return false;
				}catch(Exception exc){
					return false;
				}
			}
			
		});
		if(f.length>0){
			clust1 = f[0].listFiles()[0];
		}
		else{
			clust1 = clust.listFiles()[0].listFiles()[0];
		}
	}
	if(clust1!=null){
		System.err.println("CLUSTER FILE "+clust1.getAbsolutePath());
	}
	return clust1;	
}
public static File getClusterFile(int index){
	File clust = new File(Constants.outputDir(), "/clust_in/"+Constants.inputDir[index]+"/clusters.txt/");
	System.err.println("clustf "+clust);
	return getClusterBase(clust);
}

	public static double probCrossOverBetweenBP = 1.1e-8;
	public static int index = 0; // if a field is repeated, which one to use
	public static String indexControl = null; // which field to apply index to,
												// or all if null
	// simulation params

	public static int[] noCopies = new int[] { 2 }; // no copies for x
	public static char[][] modify0 = null;											// chromosome data
	public static char[][] modify1 = null;// new double[] {0.0};
	//public static char[] modify1 = null;
	public static double[] modifyFrac0 = null;
	public static double[] modifyFrac2 = new double[] {0.1,0.8,0.1};
	public static double[] modifyFrac3 = null;
	public static double[] modifyFrac1 = null;
	public static double[] modifyFracStart = null;
	// data modification params

	private static Boolean male = null;
	// public static boolean addNullStates = false;
	public static boolean keepSamples = false;

	public static boolean keepSamples() {
		return keepSamples;
	}

	// training params

	// public static int numFounders =8;

	// public static char[] al0 = new char[] {'A', 'B','A','B', 'N', 'X', 'Y',
	// 'Z'};
	// public static char[] al1 = null;//new char[] {'A', 'B','A','B', 'N', 'X',
	// 'Y', 'Z'};
	public static int[] transMode0 = new int[] { 4 };
	public static int[] transMode1 = null;
	public static int[] transMode2 = null;// use if we expand model - i.e.
											// hierarchical

	public static double[] u_global = new double[] { 10.0, 10.0, 10.0, 100.0 }; // em,
																				// trans,
																				// exp,
																				// background
	public static double[] u_global1 = new double[] { 10.0, 10.0, 10.0, 100.0 }; // em,
																					// trans,
																					// exp,
																					// background
 static int noPseudoCountParams = 6;
	public static PseudoIterator pseudo__ = new PseudoIterator( false);
	//public static PseudoIterator pseudoTG__ = new PseudoIterator( false);
	public static PseudoIterator pseudo2__ = new PseudoIterator(false);
	public static PseudoIterator pseudo3__ = new PseudoIterator(false);
	public static PseudoIterator pseudo1__ = new PseudoIterator(false);
	public static PseudoIterator pseudo__(boolean expanded) {
		/*if(Constants.trainGlobal && ! expanded) {
			return pseudoTG__;
		}
		else{*/
			return pseudo__;
		//}*/
			}
	// initialisation
	// public static double[] stateProbOfNull = new double[] {0.5, 0.5, 0.0};
	// public static double modify1 = 0;//0.9; //negative is modifying only for
	// copy number class
	private static boolean initialise = false;

	public static boolean fillGaps = true;
	public static int noSamplesFromHMM = 1;
	private static boolean unwrapForSampling = false;// true;//
	// public static boolean viterbi = false;
	public static boolean sample = true;

	private static boolean sampleWithPedigree = false;
	private static boolean trainWithPedigree = false;

	public static String[] format = new String[] { "hapmap" };
	public static String[] emissionGroup = null;
	public static double[] fillLikelihood = null;
	public static double[] dropFracSites = null;
	public static double[] exponentB = new double[5];
	static{
		Arrays.fill(exponentB,1);
	}
	public static double[] exponentR;
	// file, misc params
	// public static boolean search = true;
	// public static boolean all =false;
	public static String out = "summary80.txt";
	// public static String type = "";
	public static boolean overwrite = true;

	public static boolean overwrite() {
		return overwrite;
	}

	// public static String dir;
	public static boolean runFastPhase = false;
	// public static int rep = 1;
	private static String bin = "";

	public static int writeHMM = 1;
	private static double[] hotspot = new double[] { 100, 1, 1 };

	// public static int[] trainData = new int[] {100,10};

	public static boolean onlyCopyNo = false;
	public static boolean resample = false;
	private static boolean readPedigree = false;
	public static double[] expModelIntHotSpot = new double[] { 100, 100, 100 };

	// public static double trainThresh = 5;
	private static boolean annotate = false;
	private static double sampleThresh = 0.5;
	public static double precision = 0.1;
	public static double[] var_thresh = new double[] { 0.1 };
	private static int indexToTrainSWHMM;
	private static double u_exp = 0.001;
	private static double initExpTrans[] = null;// new double[] {0.99, 0.99,
												// 0.9999};
	// public static String[] inputFile = new String[] {"all_x.csv.txt"};
	public static String[] inputDir = new String[] { "./" };
	public static boolean keepBest = false;

	public static int prime = 0;

	private static boolean xchrom = false;
	private static double pseudoCountWeightClumping = 0;// this controls how
														// much indicator states
														// are encouraged to be
														// close to 1 or 0

	public static double initialConcentration = 0.0;
	public static double exclude = 0.0;// false; //whether we should exclude
										// [AA, -] [AB, -] as haplolist
										// decompositions
	public static int modifyWithData = 0;

	public static int end = Integer.MAX_VALUE;
	// parameters for cgh
	public static double[] cgh_x = // new double[] {0, 0.25, 1.0/3.0, 0.5,
									// 2.0/3.0,0.75, 1, 4.0/3.0, 1.5, 2.0, 3.0,
									// 4.0};
	new double[] { 0, 0.5, 1, 1.5, 2.0, 3.0, 4.0 };
	public static double[] cgh_mean = new double[] { -2.5, -1, 0, 0.58, 1.0,
			1.58, 2 };
	// static boolean useSkew = true;
	public static double[] cgh_var = new double[] { 0.3, 0.3, 0.3, 0.3, 0.3,
			0.3, 0.3 };
	public static double[] cgh_skew = new double[] { -0.01, -0.01, 0, 0, 0.01,
			0.01, 0.01 };
	// parameters for r distribution

	public static double[][][] r_var = new double[0][0][0] ;
    public static double[][][] r_mean = new double[0][0][0] ;
	
	public static double[][][] b_mean;
	public static double[][][] b_var1;

	public static double[][] b_skew;
	public static double round = 1000;
	public static double[][] r_train = new double[0][0];
	public static double[][] r_trainVar = new double[0][0];
	public static double[][] b_train = new double[0][0];
	public static double[][] b_trainVar = new double[0][0];
	
	
	public static double[][] r_train2 = new double[0][0];
	public static double[][] r_train2Var = new double[0][0];
	public static double[][] b_train2 = new double[0][0];
	public static double[][] b_train2Var = new double[0][0];
	public static double[][] rho_train = new double[0][0];;
	public static double[] rho_train2 = new double[0];
	//public static double[][] r_train0 = new double[0][0];
	// private static double[][] b_train0 =new double[][] {new double[]
	// {1e10,1e10
	//public static double[] rho_train0 = new double[0];

	/* 0 is mean, 1 is variances */
	public static double[] r_train(int i, int type) {
		// expand(r_train,5);
		if (i == 0)
			return r_train[type];
		else if (i == 1){
//			double[][] r_Train2 = r_train2;
			return r_train2[0];
		}
		else
			throw new RuntimeException("!!");
	}

	public static double[] b_train(int i, int type) {
		if (i == 0)
			return b_train[type];
		else if (i == 1)
			return b_train2[0];
		else
			throw new RuntimeException("!!");
	}

	

	public static double[][] r_train(int i) {
		if (i == 0)
			return r_train;
		else if (i == 1)
			return r_train2;
		else
			throw new RuntimeException("!!");
	}

	public static double[][] b_train(int i) {
		if (i == 0)
			return b_train;
		else if (i == 1)
			return b_train2;
		else
			throw new RuntimeException("!!");
	}

	public static boolean trainCGH = false;
	// public static double[] meanvarskewprior = new double[] {1.0, 10.0, 0.1};
	// //which components to train
	public static double[][] r_prior = new double[][] { { 1.0 } }; // first
																	// data-type
																	// then
																	// state

	public static double[][] r2_prior = new double[][] { { 1.0 } }; // first
																	// data-type
																	// then
	
	/*this fixes r_var, r_mean, b_mean, b_var in the case that we are training individual clusters
public static void resetIndices(){
	
}*/

	
	public static double[] r_var(int i, int i1) {
		double[][] v =  r_var[i>=r_var.length ? 0 : i];
		return v[(int)Math.min(v.length-1, i1)];
	}

	
	public static double[] r_mean(int i, int i1) {
		double[][] v =  r_mean[i>=r_mean.length ? 0 : i];
		return v[(int)Math.min(v.length-1, i1)];
	}

	public static double[] b_mean(int i, int i1) {
		double[][] v =  b_mean[i>=b_mean.length ? 0 : i];
		return v[(int)Math.min(v.length-1, i1)];
	}
	public static double[] b_var1(int i, int i1) {
		double[][] v =  b_var1[i>=b_var1.length ? 0 : i];
		return v[(int)Math.min(v.length-1, i1)];
	}
	

	// public static double[] meanvarskewprior(){
	// return meanvarskewprior;
	// }

	public static boolean onlyCopyNo() {
		return onlyCopyNo;
	}

	public static boolean resample() {
		return resample;
	}

	static void setGapParams() {

	}
	/*public static boolean modExponent(int i){
		//if(format[i].startsWith("geno") || format[i].startsWith("hap") || format[i].st)
		double fi = fillLikelihood1(i);
		if(fi<1-1e-3 && fi>1e-3){
			return true;
		}
		else return false;
	}*/
	/*public static double fillLikelihood(int i) {
		if(!modExponent(i)){
			return fillLikelihood1(i);
		}
		else return 0;
	}*/
	public static double fillLikelihood(int i) {
		if(fillLikelihood!=null){
			if(i<fillLikelihood.length) {
				return fillLikelihood[i];
			}
			else {
				return fillLikelihood[0];
			}
		}
		else{
			return 0;
		}
	}

	public static boolean readPedigree() {
		return readPedigree;
	}

	public static double expModelIntHotSpot(int i) {
		return expModelIntHotSpot[i];
	}

	private static double[] parse_(String[] cop1) {
		if (cop1.length == 1 && cop1[0].equals("null"))
			return null;
		String [] cop = extractFromX(cop1);
		 cop = extractFromI(cop);
		// String[] cop = params.getValues();
		Logger.global.info("parsing " + Arrays.asList(cop));
		List<Double> res = new ArrayList<Double>();
		double sum = 0;
		int nullCount = 0;
		for (int i = 0; i < cop.length; i++) {
			if (cop[i].length() == 0) {
				res.add(null);
				nullCount++;
				continue;
			}
			int ind = cop[i].indexOf('/');
			/*
			 * int ind1 = cop[i].indexOf(','); if(ind1>=0){ cop[i] =
			 * cop[i].replace(',', '.'); }
			 */
			if (ind >= 0) {
				double d1 = Double.parseDouble(cop[i].substring(0, ind));
				double d2 = Double.parseDouble(cop[i].substring(ind + 1));
				res.add(d1 / d2);
				sum += d1 / d2;
			} else {
				if (cop[i].startsWith("^")) {
					cop[i] = cop[i].substring(cop[i].lastIndexOf('^') + 1);
				}
				if (cop[i].indexOf("^") > 0) {
					String[] spl = cop[i].split("\\^");
					int rep = Integer.parseInt(spl[1]);
					double val = Double.parseDouble(spl[0].replace('E', 'e'));
					for (int j = 0; j < rep; j++) {
						res.add(val);
						sum += val;
					}

				} else {
					double d = cop[i].equals("null") ? null : cop[i]
							.equals("Inf") ? Double.POSITIVE_INFINITY : 
								cop[i].equals("NA") ? Double.NaN:
								Double
							.parseDouble(cop[i]);
					res.add(d);
					sum += d;
				}

			}
		}
		double[] res1;
		if (res.size() == 1 && res.get(0) == null)
			return null;
		res1 = new double[res.size()];
		for (int i = 0; i < res1.length; i++) {
			Double v = res.get(i);
			if (v == null) {
				res1[i] = (1.0 - sum) / (double) nullCount;
			} else
				res1[i] = v;

		}
		return res1;
	}

	private static String[] extractFromI(String[] cop) {
		int last =0;
		for(int i=0; i<cop.length; i++){
			if(cop[i].trim().equals("i")){
				int j=i+1;
				inner: for(; j<cop.length; j++){
					if(!cop[j].equals("i")){
						break inner;
					}
				}
				double st = Double.parseDouble(cop[last]);
				double end = Double.parseDouble(cop[j]);
				double frac = (double) (i-last) / (double)(j-last);
				cop[i] = (st + (end-st) * frac)+"";
			}
			else{
				last = i;
			}
		}
		return cop;
	}

	private static String[] extractFromX(String[] cop1) {
		int cnt =0;
		for(int i=0; i<cop1.length; i++){
			int ind = cop1[i].indexOf('x');
			if(ind>=0){
				cnt+=Integer.parseInt(cop1[i].substring(ind+1));
			}
			else{
				cnt+=1;
			}
		}
		if(cnt==cop1.length) return cop1;
		else{
			String[] res = new String[cnt];
			cnt=0;
			for(int i=0; i<cop1.length; i++){
				int ind = cop1[i].indexOf('x');
				if(ind>=0){
					int cnt1 = Integer.parseInt(cop1[i].substring(ind+1));
					for(int k=0; k<cnt1; k++){
						
						res[k+cnt] = cop1[i].substring(0,ind);
					}
					cnt+=cnt1;
				}
				else{
					res[cnt] = cop1[i];
					cnt+=1;
				}
			}
			return res;
		}
	
		
	}

	private static double[][] parse_d(String[] cop, String internalSplit) {
		// String[] cop = params.getValues();
		Logger.global.info("parsing " + Arrays.asList(cop));
		double[][] res = new double[cop.length][];
		for (int i = 0; i < cop.length; i++) {
			res[i] = parse_(cop[i].split(internalSplit));
		}
		return res;
	}

	private static int[][] parse_i(String[] cop, String intSplit) {
		// String[] cop = params.getValues();
		Logger.global.info("parsing " + Arrays.asList(cop));
		int[][] res = new int[cop.length][];
		for (int i = 0; i < cop.length; i++) {
			if(cop[i].equals("null")){
				res[i] = null;
			}else{
			res[i] = cop[i].length() == 0 ? new int[0] : parse1(cop[i]
					.split(intSplit));
			}
		}
		return res;
	}

	private static int[] parse1(String[] cop) {
		// = params.getValues();
		int[] noCopies = new int[cop.length];
		for (int i = 0; i < cop.length; i++) {
			String c = cop[i].startsWith("^") ? cop[i].substring(1) : cop[i];
			noCopies[i] = // cop[i].length()==0 ? -1 :
			Integer.parseInt(c);
		}
		return noCopies;
	}

	private static char[] parse2(String[] cop) {
		// String[] cop = params.getValues();
		List<Character> res = new ArrayList<Character>();

		for (int i = 0; i < cop.length; i++) {
			if (cop[i].indexOf("^") > 0) {
				String[] spl = cop[i].split("\\^");
				int rep = Integer.parseInt(spl[1]);
				char val = spl[0].charAt(0);
				for (int j = 0; j < rep; j++) {
					res.add(val);
				}

			} else {
				res.add(cop[i].charAt(0));
			}
			// noCopies[i] =cop[i].charAt(0);
		}
		char[] noCopies = new char[res.size()];
		for (int i = 0; i < noCopies.length; i++) {
			noCopies[i] = res.get(i);
		}
		return noCopies;
	}

	private static boolean[] parse3(String[] cop) {
		// = params.getValues();
		boolean[] noCopies = new boolean[cop.length];
		for (int i = 0; i < cop.length; i++) {
			noCopies[i] = cop[i].toLowerCase().charAt(0) == 't';
		}
		return noCopies;
	}

	private static int[] parse1(String[] cop, String def) {

		// String[] cop = params.getValues();
		if (cop == null)
			cop = def.split(":");
		if (cop[0] == "null")
			return null;
		int[] noCopies = new int[cop.length];
		for (int i = 0; i < cop.length; i++) {
			noCopies[i] = Integer.parseInt(cop[i]);
		}
		return noCopies;
	}

	public static boolean sample() {
		return sample;
	}

	public static boolean initialise() {
		return initialise;
	}

	public static CommandLine parse(String[] args, Integer column, int rep,
			String repControl) throws Exception {
		CommandLine res = parse(args, OPTIONS, column, rep, repControl, null);
		Constants.column = column;
		Constants.index = rep;
		Constants.indexControl = repControl;
		return res;
	}

	public static CommandLine parse(String[] args) throws Exception {
		return parse(args, OPTIONS, 1, 1, null, null);
	}

	public static CommandLine parseSimple(String[] args) throws Exception {
		Parser parser = new PosixParser();
		return parser.parse(OPTIONS, args, false);
	}

	public static boolean writeExtractedFile = false;
	public static String paramFile = "log.txt";
	public static int column = -1;
	public static String data_params;
	public static String loc_params;
	public static String fixed_params;
	public static String[] experiment = new String[0];

	public static Object parse(String st, String[] str, Class type, String internalSplit, String internalSplit1) {
		
		
		if (type.equals(double.class)) {

			if (st.startsWith("^")) {
				st = st.substring(1);
			}
			try {
				if(st.equals("false")) return 0;
				
				return st.equals("NA") ? Double.NaN : Double.parseDouble(st);
			} catch (Exception exc) {
				return Double.parseDouble(st.replace(',', '.'));
			}
		} else if (type.equals(double[].class)) {
			System.err.println(Arrays.asList(str));
			return parse_(str);

		} else if (type.equals(double[][].class)) {
			return parse_d(str, internalSplit);
		} else if (type.equals(String[].class)) {
			return str;
		} else if (type.equals(String[][].class)) {
			String[][] res = new String[str.length][];
			for (int i = 0; i < res.length; i++) {
				res[i] = str[i].split(internalSplit);
			}
			return res;
		} else if (type.equals(char[][].class)) {
			char[][] res = new char[str.length][];
			for (int i = 0; i < res.length; i++) {
				String[] strs = str[i].split(internalSplit);
				res[i] = new char[ strs.length];
				for(int k=0; k<strs.length; k++){
					res[i][k] = strs[k].trim().charAt(0);
				}
			}
			return res;
		}  
		else if (type.equals(String[][][].class)) {
	    	String[][][] res = new String[str.length][][];
		    for (int i = 0; i < res.length; i++) {
		    	String[]res1 = str[i].split(internalSplit);
			    res[i] = new String[res1.length][];
			    if(res1.length==1 && res1[0].equals("null")) res[0][0] = new String[] {"null"};
			    else if(res1.length==1 && res1[0].equals("mid")) res[0][0] = new String[] {"mid"};
			    else if(res1.length==1 && res1[0].equals("all")) res[0][0] = new String[] {"all"};
			    else{ for(int ij =0; ij < res1.length; ij++){
			    	res[i][ij] = res1[i].split(internalSplit1);
			    }
			    }
		    }
		    return res;
	     } 
		else if (type.equals(char[].class)) {
			return parse2(str);
		} else if (type.equals(boolean[].class)) {
			return parse3(str);
		} else if (type.equals(int[].class)) {
			return parse1(str);
		} else if (type.equals(int[][].class)) {
			return parse_i(str,internalSplit);
		} else if (type.equals(int.class)) {
			return Integer.parseInt(st);
		} else if (type.equals(boolean.class)) {
			return parseBoolean(st);

		} else if (type.equals(Boolean.class)) {
			return parseBoolean(st);
		} else if (type.equals(long.class)) {
			return Long.parseLong(st);
		} else if (type.equals(String.class)) {
			return st;
		} else
			return null;

	}

	private static Boolean parseBoolean(String st) {
		if (st.toLowerCase().equals("true"))
			return true;
		else if (st.toLowerCase().equals("false"))
			return false;
		else {
			return null;
//			throw new RuntimeException("is not boolean " + st);
		}
	}

	public static boolean plotFlux = true;

	public static boolean plotFlux() {
		return plotFlux;
	}

	private static String[] originalInput;
	public static String[] include;

	public static String[] originalNames() {
		return originalInput;
	}

	public static void parseInner(Option opt, Options options, int column,
			int index, String repControl, CommandLine paramsa, CommandLine paramsb) throws Exception {
		String argName = opt.getLongOpt();
		
		int underscore_index = argName.indexOf("__");
		if (underscore_index >= 0) {
			String base = argName.substring(0, underscore_index + 2);
			String extension = argName.substring(underscore_index + 2);
			String[] st = opt.getValues();
			Field field = Constants.class.getField(base);
			PseudoIterator res = (PseudoIterator) field.get(null);
			if (res == null) {
				field.set(null, res = new PseudoIterator(st.length, true));
			}
			res.set(extension, parse_(opt.getValues()));
			return;
		}
		
		Field field = Constants.class.getField(argName);
		Class type = field.getType();
		if (type.equals(String.class)) {
			String val = opt.getValue();
			if (field.getName().startsWith("fixed_params")
					|| field.getName().startsWith("loc_params")) {
				String[] args1 = Constants.read(new File(opt.getValue()),
						column, index, repControl,paramsa);
				Constants.parse(args1, options, 1, 1, repControl, paramsa);
			}
		 
			 else if (field.getName().startsWith("data_params")) {
				 parseDataParams(val, paramsa, paramsb);
			 }
			field.set(null, val);
		} else {
			String val = opt.getValue();
			String[] vals = opt.getValues();
			if(argName.equals("modify0")){
				String[] inputDir = Constants.inputDir;
				
				if(val.indexOf(';')<0){
				   String str1 = Arrays.asList(vals).toString().replaceAll(",",";");
				   str1 = str1.substring(1,str1.length()-1);
				   vals = new String[inputDir.length];
				   Arrays.fill(vals, str1);
				}
			}
			if(opt.getDescription().equals("restrictKb") && vals.length==2){
				vals = new String[] {vals[0]+";"+vals[1]};
			}
			try{
			Object valu = Constants
					.parse(val, vals, type,";","~");
			if (valu != null) {
				field.set(null, valu);
				
			}
			}catch(Exception exc){
				System.err.println("prob with "+val+" "+field.getName());
				exc.printStackTrace();
				System.exit(0);
			}

		}
		if(field.getName().startsWith("seed")){
			
				Constants.rand = new Random(seed);
			}
		Object res = field.get(null);
		if(field.getName().equals("probeOnly")){
			System.err.println("ḧ");
		}
		System.err.println("Set " + field.getName() + " : " + (getString(res)));
		
	}

	
	private static void parseDataParams(String val, CommandLine paramsa,CommandLine paramsb)  throws Exception{
		{
			 {
					File ff = new File(val);
					if(!ff.exists()){
						ff = new File("../"+ff.getName());
						
					}
					
					// String[] exp = experiment;
					AbstractDataParams ws = getDataParams(ff,sheet);
					
					
					// BufferedReader br = new BufferedReader(new FileReader(new
					// File(val)));
					List<String> include = new ArrayList<String>(); // list of
																	// all files
																	// to
																	// include
					List<Integer> alias = new ArrayList<Integer>(); // list of
																	// alias to
																	// columns
																	// in
																	// spreadsheet
					int inputDir_index = -1;
					{
						String[] row0 = ws.getRow(0);
						for (int i = 0; i < ws.getRows(); i++) {
							String nme =row0[i].trim();
							if (nme.equals(
									"inputDir")) {
								inputDir_index = i;
								break;
							}
						}
						if (inputDir_index < 0)
							throw new RuntimeException(
									"need to have a row called inputDir in data spreadsheet");
					}
					String[] vals_ = new String[length(ws
							.getRow(inputDir_index))];
					String[] inputDir = (String[]) getVals(process(ws
							.getRow(inputDir_index), vals_), String.class, null);

					int[] toInclude = new int[inputDir.length];// (Boolean[])
																		// getVals(process(ws.getRow(0),
																		// vals_),
																		// Boolean.class,
																		// null);
					// String[] inputDir = (String[])
					// getVals(process(ws.getRow(1),vals_), String.class,
					// toInclude);
//					int[] toincl = toInclude;
					originalInput = inputDir;
					List<String> origOrder = Arrays.asList(inputDir);
					List<Integer> indices = new ArrayList<Integer>();
					
					List<String> toIncludeSet = Arrays
							.asList(Constants.include);
					toInclude = new int[toIncludeSet.size()];
					for(int i=0; i<toInclude.length; i++){
						toInclude[i] = origOrder.indexOf(toIncludeSet.get(i));
						if(toInclude[i]<0){
							throw new RuntimeException("nothing found for "+toIncludeSet.get(i));
						}
					}
					/*for (int i = 0; i < toInclude.length; i++) {
						toInclude[i] = toIncludeSet.indexOf(inputDir[i]);
						if(toInclude[i]>=0)toIncludeSet.set(toInclude[i], "done");
						if (toInclude[i]<0) {
							inputDir[i] = null;
							
						} 
					}*/
					System.err.println("input dir " + getString(inputDir));
					originalInput = inputDir;
					// expand(inputDir, null, include, alias);
					// inputDir = include.toArray(new String[0]);
					System.err
							.println("input dir after " + getString(inputDir));

					System.err.println("to include  " + getString(toInclude));
					String[] inputDirNew = (String[]) thin(inputDir, toInclude,
							String.class);
					originalInput = (String[]) thin(originalInput, toInclude,
							String.class);
					System.err.println("set inputDir" + " : "
							+ getString(inputDirNew));
					Constants.class.getField("inputDir").set(null, inputDirNew);
					inner: for (int i = 0; i < ws.getRows(); i++) {// (st
																	// =br.readLine())!=null){
						String[] cell = ws.getRow(i);
						if (cell.length == 0
								|| cell[0].trim().startsWith("#")
								|| cell[0].trim().length() == 0)
							continue;
						if (cell[0].trim().equals("inputDir")
								|| cell[0].trim().equals(
										"include"))
							continue;
						// if(st.startsWith("#")) continue;
						// String[] header = st.split("\t");
						try{
						Field f = Constants.class.getField(cell[0]
								);
						if (f.getName().equals("ldtype")) {
							System.err.println("ignoring " + ldtype);
							continue inner;
						}
						Class cl = f.getType().getComponentType();
						System.err.println(f.getName());
						Object vals = getVals(process(cell, vals_), cl,
								toInclude);
						// Object va1 = expand(vals, alias, cl);
						//Object va = thin(vals, toInclude, cl);
						if(!paramsa.hasOption(f.getName()) && !paramsb.hasOption(f.getName())){
							
						/*if(f.getName().equals("b_var1")){
							System.err.println("h");
						}*/
						System.err.println("set " + f.getName() + " : "
								+ getString(vals));
						if(!f.getName().equals("restrictKb") || restrictKb==null){
						f.set(null, vals);
						}
						}
						else{
							System.err.println("ignoring "+f.getName()+" because defined in commandline");
						}
						}catch(Exception exc){
							Logger.global.warning(exc.getMessage()+": "+cell[0]);
						}
					}
				}
			 }

		
		
	}
	private static int length(String[] row) {
		for (int i = 0; i < row.length; i++) {
			if (row[i].length() == 0
					|| row[i].trim().startsWith("#")) {
				return i;
			}
		}
		return row.length;
	}

	private static String[] getFiles(File file) {
		List<File> l = new ArrayList<File>();
		listFiles(file, l);
		System.err.println("list of files is " + l);
		int path = file.getAbsolutePath().length() + 1;
		String[] st = new String[l.size()];
		for (int i = 0; i < st.length; i++) {
			st[i] = l.get(i).getAbsolutePath().substring(path);
		}
		System.err.println("files are " + Arrays.asList(st));
		return st;
	}

	/** returns true if added files */
	private static boolean listFiles(File file, List<File> l) {
		File[] f = file.listFiles(new FileFilter() {
			public boolean accept(File pathname) {
				return pathname.isDirectory();

			}
		});
		if (f.length == 0)
			return false;
		for (int i = 0; i < f.length; i++) {

			if (!listFiles(f[i], l)) {
				l.add(f[i]);
			}
		}
		return file.isDirectory();

	}

	public static String regionFile = null;

	private static void expand(String[] inputDir2, String[] strings,
			List<String> include, List<Integer> alias) {
		for (int i = 0; i < inputDir2.length; i++) {
			include.add(inputDir2[i]);
			alias.add(i);
		}
		/*
		 * outer: for(int i=0; i<strings.length; i++){ String str_ = strings[i];
		 * 
		 * 
		 * for(int j=0; j<inputDir2.length; j++){ int ind =
		 * inputDir2[j].indexOf('*'); if((ind>=0 &&
		 * str_.startsWith(inputDir2[j].substring(0,ind)) &&
		 * str_.split("/").length==inputDir2[j].split("/").length ) ||
		 * strings[i].equals(inputDir2[j])){ include.add(strings[i]);
		 * alias.add(j); continue outer; }else if(ind>=0){ //
		 * System.err.println(str_+"\n"+inputDir2[j].substring(0,ind)); } } }
		 */
	}

	private static String getString(Object res) {
		if (res == null)
			return "null";
		if (res.getClass().isArray()) {
			int len = Array.getLength(res);
			if (len == 0)
				return "[]";
			else {
				StringBuffer sb = new StringBuffer("[");
				for (int i = 0; i < len - 1; i++) {
					sb.append(getString(Array.get(res, i)) + ",");
				}
				sb.append(getString(Array.get(res, len - 1)) + "]");
				return sb.toString();
			}

		} else
			return res.toString();
	}

	public static void copy(Object[] format_new, String name) throws Exception {
		Field field = Constants.class.getField(name);
		Object obj = field.get(null);
		if (obj == null)
			return;
		Class inner = Array.get(obj, 0).getClass();
		if (inner.equals(Double.class))
			inner = double.class;
		if (inner.equals(Integer.class))
			inner = int.class;
		Object novo = Array.newInstance(inner, format_new.length);
		for (int i = 0; i < format_new.length; i++) {
			Array.set(novo, i, format_new[i]);
		}
		field.set(null, novo);

	}

	private static Object thin(Object value, int[] toInclude2, Class cl) {

		int cnt = toInclude2.length;
		
		Object value1 = Array.newInstance(cl, cnt);
		//cnt = 0;
		for (int i = 0; i < toInclude2.length; i++) {
		
				Array.set(value1,i,  Array.get(value, toInclude2[i]));
			//	cnt++;

			
		}
		return value1;
	}

	private static Object expand(Object vals, List<Integer> alias, Class cl) {
		Object value = Array.newInstance(cl, alias.size());
		for (int i = 0; i < alias.size(); i++) {
			Array.set(value, i, Array.get(vals, alias.get(i)));
		}
		return value;
	}

	private static Object getVals(String[] header, Class cl) throws Exception {
		// Field f = Constants.class.getField(header[0]);
		// Class cl = f.getType().getComponentType();

		Object value = Array.newInstance(cl, header.length - 1);

		for (int ik = 0; ik < header.length - 1; ik++) {
			Array.set(value, ik, parse(header[ik + 1], header[ik + 1]
					.split(":"), cl,";","~"));
		}
		return value;
	}

	private static String[] process(String[] header, String[] vals) {
		for (int i = 0; i < vals.length; i++) {
			vals[i] = header[i];
		}
		return vals;
	}

	private static Object getVals(String[] header, Class cl, int[] b)
			throws Exception {
		// Field f = Constants.class.getField(header[0]);
		// Class cl = f.getType().getComponentType();
		int len =b==null? header.length-1 :  b.length;
		
		// if(len<header.length && header[len].getContents().length()==0)
		// len--;
		// if(header[0].getContents().equals("standardVar")){
		// Logger.global.info("h");
		// }
		Object value = Array.newInstance(cl, len);
		// if(cl.getName().equals("st))
		for (int ik = 0; ik < len ; ik++) {
			//if (b == null || b[ik]>=0) {
				String cnts = b==null ? header[ik+1] : header[b[ik]+1];
				//if (cnts.length() == 0)
				//	return value;
				Array.set(value, ik, parse(cnts, cnts.split(";"), cl,"~",","));
			//}
		}
		return value;
	}

	public static CommandLine parse(String[] args, Options opti,
			Integer column, int rep, String repControl, CommandLine params1) throws Exception {
		Parser parser = new PosixParser();
		final CommandLine params = parser.parse(opti, args, false);
		if (params.hasOption("paramFile")) {
			if (params.getOptionValue("paramFile").equals("many"))
				return null;

			File user;
			if (params.hasOption("dir")
					&& !params.getOptionValue("dir").equals(".")) {
				user = new File(params.getOptionValue("dir"));
			} else
				user = new File(System.getProperty("user.dir"));
			String[] args1 = read(new File(user, params
					.getOptionValue("paramFile")), column, rep, repControl, params);
			Constants.column = column;
			Constants.index = rep;
			Constants.indexControl = repControl;
			Option[] options = params.getOptions();
			
			final List<String> firstItems = Arrays.asList("include:sheet:baseDir:data_params".split(":"));
			Arrays.sort(options, new Comparator<Option>() {

				public int compare(Option o1, Option o2) {
					String st1 = o1.getLongOpt();
					String st2 = o2.getLongOpt();
					int i1 = firstItems.indexOf(st1);
					int i2 = firstItems.indexOf(st2);
					i1 = i1<0 ? firstItems.size() : i1;
					i2 = i2<0 ? firstItems.size() : i2;
					if(i1 < i2) return -1;
					else if (i2 < i1) return 1;
					else return st1.compareTo(st2);
					
				}

			});
			for (int i = 0; i < options.length; i++) {
				String nme = options[i].getDescription();
				String val = options[i].getValue();
				if (nme.equals("modify0")) {
					System.err.println(val);
				}
				parseInner(options[i], opti, column, rep, repControl,params1, params);

			}
			return parse(args1, opti, column, rep, repControl,params);
		}
		// if(params.hasOption("type")){
		// type = params.getOptionValue("type");
		// }
		Option[] options = params.getOptions();
		Arrays.sort(options, new Comparator<Option>() {

			public int compare(Option o1, Option o2) {
				String st1 = o1.getLongOpt();
				String st2 = o2.getLongOpt();
				if (st1.startsWith("include"))
					return -1;
				else if (st2.startsWith("include"))
					return +1;
				else if (st1.startsWith("sheet"))
					return -1;
				else if (st2.startsWith("sheet"))
					return +1;
				else if (st1.startsWith("baseDir"))
					return -1;
				else if (st2.startsWith("baseDir"))
					return +1;
				else
					return st1.compareTo(st2);
			}

		});
		for (int i = 0; i < options.length; i++) {
			parseInner(options[i], opti, column, rep, repControl,params1, params);

		}

		return params;
	}

	private static String[] read(File file, int col, int rep, String repControl1, CommandLine params)
			throws Exception {
		// List<String> l = new ArrayList<String>();
		String repControl = repControl1 == null ? null : "--" + repControl1;
		BufferedReader br = new BufferedReader(new FileReader(file));
		String st = "";
		Pattern pat = Pattern.compile("-[0-9]");
		
		// int cnt=0;
		// String prev = "";
		Map<String, List<String>> m = new HashMap<String, List<String>>();
		while ((st = br.readLine()) != null) {
			if (st.trim().startsWith("#"))
				continue;
			//st = st.split(" #")[0].trim();
			if(st.length()==0) continue;
			String[] str = st.trim().split("\\s+");
			
			if (st.lastIndexOf('-') >= 2) {
				System.err.println(st);
				
				for (int i = 1; i < str.length; i++) {
					Matcher mat = pat.matcher(str[i]);
					
					mat.useAnchoringBounds(false);
					boolean match =  mat.find();
					str[i] = str[i].replaceAll("E", "e");
					if(!str[0].equals("annotateSamples") && (match || str[i].matches("e-"))){
						str[i] = str[i].replaceAll("e-", "e_").replaceAll("-", "^-").replaceAll("e_", "e-");
					}
				}
			}
			// System.err.println(str[0]+" ");
			// System.err.println(str[1]);
			if (str.length < 2)
				continue;
			if(params.hasOption(str[0].substring(2))){
				System.err.println("ignoring "+str[0]+" because defined in commandline");
				continue;
			}
			List<String> l = m.get(str[0].trim());
			if (l == null)
				m.put(str[0].trim(), l = new ArrayList<String>(1));
			/*
			 * if(str[0].equals(prev)){ cnt++; } else{ cnt=0; }
			 */
			// prev = str[0];
			// if(cnt==0)l.add(str[0].trim());//st.substring(0,sp));
			String toadd = "";
			if (st.startsWith("--dir")) {
				StringBuffer toA = new StringBuffer();
				for (int i = 1; i < str.length - 1; i++) {
					toA.append(str[i]);
					toA.append("%");
				}
				toA.append(str[str.length - 1]);
				String s = toA.toString();
				s = s.replace(':', ';');
				toadd = s;
			} else {
				int ind = str.length == 2 ? 1 : col;
				toadd = str[Math.min(ind, str.length - 1)].trim();// st.substring(sp+1));
			}
			// if(cnt==0
			// )
			l.add(toadd);
			/*
			 * else if(cnt<rep ){ if (repControl==null ||
			 * repControl.equals(str[0].trim())){ l.set(l.size()-1,toadd); }
			 * else{
			 * System.err.println("no match "+str[0].trim()+" "+repControl); } }
			 */
		}
		//List l1_ = m.get("--regionsToExclude");
		List<String> l1 = new ArrayList<String>();
		for (Iterator<Entry<String, List<String>>> it = m.entrySet().iterator(); it
				.hasNext();) {
			Entry<String, List<String>> nxt = it.next();
			try{
			l1.add(nxt.getKey());
			// if(nxt.getKey().equals("--experiment")){
			// Logger.global.info("h");
			// }
			List<String> l2 = nxt.getValue();
			if (l2.size() == 1) {
				l1.add(l2.get(0));
			} else if (repControl != null && repControl.equals(nxt.getKey())) {

				l1.add(l2.get(rep - 1));
			} else if (repControl == null) {
				l1.add(l2.get(rep - 1));
			} else {
				l1.add(l2.get(0));
			}
			}catch(Exception exc){
				exc.printStackTrace();
				System.err.println(nxt);
			}
		}
		// TODO Auto-generated method stub
		br.close();
		return l1.toArray(new String[0]);
	}

	public static boolean fast() {
		// TODO Auto-generated method stub
		return fast;
	}

	public static boolean fillGaps() {
		return fillGaps;
	}

	public static int getMax(Double[] d) {
		double max = d[0];
		int max_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] > max) {
				max_id = i;
				max = d[i];
			}
		}
		return max_id;
	}

	public static int getMax(Integer[] d) {
		int max = d[0];
		int max_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] > max) {
				max_id = i;
				max = d[i];
			}
		}
		return max_id;
	}

	public static int getMax(double[] d) {
		double max = Double.NEGATIVE_INFINITY;
		int max_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (!Double.isNaN(d[i]) && d[i] > max) {
				max_id = i;
				max = d[i];
			}
		}
		return max_id;
	}

	public static int getMax(int[] d) {
		int max = d[0];
		int max_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] > max) {
				max_id = i;
				max = d[i];
			}
		}
		return max_id;
	}

	public static int getMin(double[] d) {
		double min = d[0];
		int min_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] < min) {
				min_id = i;
				min = d[i];
			}
		}
		return min_id;
	}

	public static int getMin(Double[] d) {
		double min = d[0];
		int min_id = 0;
		for (int i = 0; i < d.length; i++) {
			if (d[i] < min) {
				min_id = i;
				min = d[i];
			}
		}
		return min_id;
	}

	// public static double[] getThresh() {
	// return thresh;
	// }
	public static int maxIndiv(int i) {
		return i<maxIndiv.length ? maxIndiv[i] : maxIndiv[0];
	}

	public static String[] phenToPlot = new String[0];

	public static String[] phenToPlot() {
		return phenToPlot;
	}

	public static int segments = 4;

	public static int nextInt(int tot) {
		return rand.nextInt(tot);
		// return (int)Math.floor(rand.nextDouble()*(double)tot);
	}

	public static int[] noCopies() {
		if (noCopies.length < format.length)
			return Constants.extend(noCopies, format.length);
		else
			return noCopies;
	}

	public static int numF(int i) {
		return modify(i).length;
		// if(modify(1)!=null) res*=modify(1).length;
		// return res;

	}

	public static int[] transMode(int i) {
		if (i == 2)
			return transMode2;
		return i == 0 ? transMode0 : transMode1;
	}

	public static int[] numIt() {
		return numIt;
	}

	public static int numRep() {
		return numRep;
	}

	public static double switchU() {
		return 1e10;
	}

	public static int sample(double[] d) {
		return sample(d, 1.0);
	}

	public static int sample(double[] d, double sum) {
		double ra = rand.nextDouble() * sum;
		double cum = 0;
		for (int i = 0; i < d.length; i++) {
			if(!Double.isNaN(d[i])){
			cum += d[i];
			if (ra <= cum)
				return i;
			}
		}
		throw new RuntimeException("!!");
	}
	
	

	static public double[] u_global(int i) {
		if (i == 0)
			return u_global;
		else
			return u_global1;
	}

	/*
	 * public static double modelHemizygous() { return
	 * probHomozygousIsHemizygous; }
	 */

	public static int writeHMM() {
		// TODO Auto-generated method stub
		return writeHMM;
	}

	public static int[] restrict() {
		return restrict;
	}

	public static int noSamples() {
		return noSamplesFromHMM;
	}

	/*
	 * public static double[] stateProbOfNull() { return stateProbOfNull; }
	 */
	public static char[] modify(int i) {
		if (i <modify0.length)
			return modify0[i];
		else
			return modify0[0];

		// transMode[i]<5 ?1.0 : transMode[i]==5 ? 0 : -1;
	}
	public static char[] modify1(int i) {
		if (i <modify1.length)
			return modify1[i];
		else
			return modify1[0];

		// transMode[i]<5 ?1.0 : transMode[i]==5 ? 0 : -1;
	}

	public static double[] extend(double[] d, int len) {
		double[] res = new double[len];
		System.arraycopy(d, 0, res, 0, d.length);
		for (int i = d.length; i < res.length; i++) {
			res[i] = d[d.length - 1];
		}
		return res;
	}

	public static int[] extend(int[] d, int len) {
		int[] res = new int[len];
		System.arraycopy(d, 0, res, 0, d.length);
		for (int i = d.length; i < res.length; i++) {
			res[i] = d[d.length - 1];
		}
		return res;
	}

	public static double[] modifyFrac(int i) {
		// return i==0 ?
		// (modifyFrac0==null ? null : modifyFrac0)
		// : modifyFrac1==null ? null : modifyFrac1;
		double[] res = new double[Constants.modify(i>=2? 0 : i).length];
		double sum = 0;
		double[] mod = i == 0 ? modifyFrac0 : (i==1 ? modifyFrac1 : (i==2 ? modifyFrac2 : modifyFrac3));
		int k = 0;
		for (; k < res.length; k++) {

			res[k] = mod[k];

			sum += res[k];
		}
		if (k < res.length) {
			double rem = (1.0 - sum) / (res.length - k);
			for (; k < res.length; k++) {
				res[k] = rem;
			}
		}
		Constants.normalise(res);
		return res;

		// transMode[i]<5 ?1.0 : transMode[i]==5 ? 0 : -1;
	}

	public static double[] modifyFracStart() {
		// return i==0 ?
		// (modifyFrac0==null ? null : modifyFrac0)
		// : modifyFrac1==null ? null : modifyFrac1;
		double[] res = new double[Constants.modify(0).length];
		double sum = 0;
		double[] mod = modifyFracStart;
		int k = 0;
		for (; k < res.length; k++) {

			res[k] = mod[k];

			sum += res[k];
		}
		if (k < res.length) {
			double rem = (1.0 - sum) / (res.length - k);
			for (; k < res.length; k++) {
				res[k] = rem;
			}
		}
		Constants.normalise(res);
		return res;

		// transMode[i]<5 ?1.0 : transMode[i]==5 ? 0 : -1;
	}

	public static double[][] r_state_mean = new double[0][0];// first index is
																// by data type
																// second by
																// state
	public static double[][] r_state_var = new double[0][0];
	public static double[][] r_state_skew = new double[0][0];
	public static double[][] r2_state_mean = new double[0][0];// first index is
																// by data type
																// second by
																// state
	public static double[][] r2_state_var = new double[0][0];
	public static double[][] r2_state_skew = new double[0][0];

	/** k is state. returns [datatype, mean/var/skew] */
	public static double[][] meanvarskew(int state_index, boolean single) {
		double[][] r_sm = single ? r_state_mean : r2_state_mean; // first
																	// data_type,
																	// second
																	// state
		double[][] r_svar = single ? r_state_var : r2_state_var; // first
																	// data_type,
																	// second
																	// state
		double[][] r_sskew = single ? r_state_skew : r2_state_skew; // first
																	// data_type,
																	// second
																	// state
		double[][] r_sprior = single ? r_prior : r2_prior; // first data_type,
															// second state
		double[][] res = new double[Constants.format().length][4];
		// double[] mod = Constants.r_state_var_mod;
		// if(mod.length!=res.length) throw new RuntimeException("!!");
		// int k=0;
		if (true) {
			for (int j = 0; j < res.length; j++) {
				int j1 = Math.min(r_sm.length - 1, j);
				if (r_sm[j1] != null) {
					int state_index1 = Math.min(r_sm[j1].length - 1,
							state_index);
					res[j][0] = r_sm[j1][state_index1];
					res[j][1] = r_svar[j1][state_index1];
					res[j][2] = r_sskew[j1][state_index1];
					res[j][3] = r_sprior[j1][state_index1];
				} else {
					res[j] = null;
				}
			}
		} else {
			throw new RuntimeException("!!");
			/*
			 * int k1 = r_sm[0].length-1; for(int j=0; j<res.length; j++){
			 * res[j][0] = r_sm[j][k1]; res[j][1] = r_svar[j][k1]; res[j][2] =
			 * r_sskew[j][k1]; res[j][3] = r_sprior[j][k1]; }
			 */
		}
		return res;

		// transMode[i]<5 ?1.0 : transMode[i]==5 ? 0 : -1;
	}

	public static boolean sampleWithPedigree() {
		// TODO Auto-generated method stub
		return sampleWithPedigree;
	}

	public static boolean trainWithPedigree() {
		return trainWithPedigree;
	}

	public static boolean unwrapForSampling() {
		// TODO Auto-generated method stub
		return unwrapForSampling;
	}

	public static final Options OPTIONS = new Options() {
		{
			String[] optional_extensions = PseudoIterator.optional_extensions;
			Field[] f = Constants.class.getFields();
			for (int i = 0; i < f.length; i++) {
				if (Modifier.isStatic(f[i].getModifiers())) {
					if (f[i].getName().endsWith("__")) {
						for (int k = 0; k < optional_extensions.length; k++) {
							String nme = f[i].getName()
									+ optional_extensions[k];
							this.addOption(OptionBuilder.withLongOpt(nme)
									.withDescription(nme).withValueSeparator(
											':').hasArgs().create());

						}
					}
					// System.err.println(f[i].getName());
					else {
						this.addOption(OptionBuilder
								.withLongOpt(f[i].getName()).withDescription(
										f[i].getName()).withValueSeparator(':')
								.hasArgs().create());
					}
				} else {
					System.err.println("excluded " + f[i]);
				}
			}
		}
	};

	public static void printOptions(PrintWriter pw, String end) {
		printOptions(pw, end, null);
	}

	public static void printOptions(PrintWriter pw, String end,
			CommandLine options) {
		try {
			// Calendar cal = new GregorianCalendar();
			// URL fi = new URL(dir, "log.txt");
			// File f1 = new File();
			// File fi =
			// pw.println(cal.getTime());
			Field[] f = Constants.class.getFields();
			/*
			 * Arrays.sort(f, new Comparator<Field>(){
			 * 
			 * public int compare(Field o1, Field o2) { return
			 * o1.getName().compareTo(o2.getName()); }
			 * 
			 * });
			 */
			for (int i = 0; i < f.length; i++) {
				if (options != null && !options.hasOption(f[i].getName()))
					continue;
				if (Modifier.isFinal(f[i].getModifiers()))
					continue;
				if (Modifier.isStatic(f[i].getModifiers())) {
					if (f[i].getType().equals(Random.class))
						continue;

					else if (f[i].getType().equals(double[].class)) {
						double[] val = (double[]) f[i].get(null);
						if (val != null) {
							Double[] d = new Double[val.length];
							StringBuffer sb = new StringBuffer();
							for (int ik = 0; ik < d.length; ik++) {
								d[ik] = val[ik];
								if (d[ik] < 0)
									sb.append("^");
								sb.append(d[ik]);
								sb.append((ik < d.length - 1 ? ":" : ""));
							}
							pw.print(" --" + f[i].getName() + " "
									+ sb.toString() + end);
						}
					} else if (f[i].getType().equals(double[][].class)) {

						double[][] val = (double[][]) f[i].get(null);
						if (val != null) {
							for (int k = 0; k < val.length; k++) {
								if (val[k] != null) {
									Double[] d = new Double[val[k].length];
									StringBuffer sb = new StringBuffer();
									for (int ik = 0; ik < d.length; ik++) {
										d[ik] = val[k][ik];
										if (d[ik] < 0)
											sb.append("^");
										sb.append(d[ik]);
										sb
												.append((ik < d.length - 1 ? ":"
														: ""));
									}
									pw.print(" --" + f[i].getName() + " "
											+ sb.toString() + end);
								} else {
									pw.print(" --" + f[i].getName() + " null");
								}
							}
						}
					} else if (f[i].getType().equals(int[].class)) {
						int[] val = (int[]) f[i].get(null);
						if (val != null) {
							Integer[] d = new Integer[val.length];
							StringBuffer sb = new StringBuffer();
							for (int ik = 0; ik < d.length; ik++) {
								d[ik] = val[ik];
								sb.append(d[ik]
										+ (ik < d.length - 1 ? ":" : ""));
							}
							pw.print(" --" + f[i].getName() + " "
									+ sb.toString() + end);
						}
					} else if (f[i].getType().equals(char[].class)) {
						char[] val = (char[]) f[i].get(null);
						if (val != null) {
							StringBuffer sb = new StringBuffer();
							for (int ik = 0; ik < val.length; ik++) {
								sb.append(val[ik]
										+ (ik < val.length - 1 ? ":" : ""));
							}
							pw.print(" --" + f[i].getName() + " "
									+ sb.toString() + end);
						}
					} else if (f[i].getType().equals(boolean[].class)) {
						boolean[] val = (boolean[]) f[i].get(null);
						if (val != null) {
							StringBuffer sb = new StringBuffer();
							for (int ik = 0; ik < val.length; ik++) {
								sb.append(val[ik]
										+ (ik < val.length - 1 ? ":" : ""));
							}
							pw.print(" --" + f[i].getName() + " "
									+ sb.toString() + end);
						}
					} else if (f[i].getType().equals(String[].class)) {
						String[] val = (String[]) f[i].get(null);
						if (val != null) {
							StringBuffer sb = new StringBuffer();
							for (int ik = 0; ik < val.length; ik++) {
								sb.append(val[ik]
										+ (ik < val.length - 1 ? ":" : ""));
							}
							pw.print(" --" + f[i].getName() + " "
									+ sb.toString() + end);
						}
					} else if (f[i].getType().equals(double.class)
							|| f[i].getType().equals(Double.class)) {
						double val = (Double) f[i].get(null);
						pw.print(" --" + f[i].getName() + " ");
						if (val < 0)
							pw.print("^" + val);
						else
							pw.print(val);
						pw.print(end);
					} else {
						pw.print(" --" + f[i].getName() + " " + f[i].get(null)
								+ end);
					}
				}
			}

		} catch (Exception exc) {
			exc.printStackTrace();
		}
	}

	public static boolean runFastPhase() {
		return runFastPhase;
	}

	public static Boolean male() {
		// TODO Auto-generated method stub
		return male;
	}

	/**
	 * index is [over long evolutionary distance, long time but from gap to non
	 * gap, from parents to offspring]
	 */
	// simulation only
	public static double hotspot(int index) {
		// TODO Auto-generated method stub
		return hotspot[index];
	}

	public static String inputFile(int i, int j) {
		return inputDir(i)+"/" + Constants.chrom(j) + ".zip";
	}
	
	

	public static String[] inputDirLoc = null;
	
	public static String inputDir(int i) {
		return baseDir + "/" + 
		(inputDirLoc==null ? inputDir[i] :  inputDirLoc[i]);
	}

	public static boolean keepBest() {
		// TODO Auto-generated method stub
		return keepBest;
	}

	public static Integer maxCopies = null;

	public static int maxCopies() {
		if (maxCopies == null) {
			int val = 1;
			// int val_min =1;
			char[][] c = Constants.modify0;
			for(int k=0; k<c.length; k++){
			for (int i = 0; i < c[k].length; i++) {
				try{
					if(c[k][i]!='a'){
				int v = Integer.parseInt(c[k][i] + "");
				if (v > val) {
					val = v;
				}
					}
				}catch(Exception exc){
					
				}
				// if(v<val_min) val_min = v;
			}
			}
			maxCopies = new Integer(val);
		}
		return maxCopies;
	}

	static Integer modelCNP = null;

	

	public static String[] getFormat() {
		return format;
	}

	// public static boolean sampleFromState = false;

	public static int offset() {
		// TODO Auto-generated method stub
		return offset;
	}

	public static int prime() {
		// TODO Auto-generated method stub
		return prime;
	}

	/**
	 * index ==0 is for normal exp_transitions index ==1 is for multi-class
	 * exp_transitions (e.g. gap to non gap) index ==2 is for switching in the
	 * child of a trio staet
	 */
	public static double initExpTrans(int index) {
		// TODO Auto-generated method stub
		return initExpTrans[index];
	}

	public static boolean cache() {
		// TODO Auto-generated method stub
		return cache;
	}

	public static double sum(double[] counts) {
		double sum = 0;
		for (int i = 0; i < counts.length; i++) {
			sum += counts[i];
		}
		return sum;
	}

	public static double sum(Double[] counts) {
		double sum = 0;
		for (int i = 0; i < counts.length; i++) {
			sum += counts[i];
		}
		return sum;
	}

	public static boolean xchrom() {
		// TODO Auto-generated method stub
		return xchrom;
	}

	public static double u_exp() {
		return u_exp;
	}

	public static int indexToTrainSWHMM() {
		return indexToTrainSWHMM;
	}

	public static double precision() {
		return precision;
	}

	public static double sampleThresh() {
		return sampleThresh;
	}

	public static double pseudoCountWeightClumping() {
		return pseudoCountWeightClumping;
	}

	public static Double initialConcentration() {
		return initialConcentration;
	}

	public static boolean annotate() {
		// TODO Auto-generated method stub
		return annotate;
	}

	public static double exclude() {
		// TODO Auto-generated method stub
		return exclude;
	}

	public static boolean trainWithGenotypes() {
		return trainWithGenotypes;
	}

	public static int numItSum() {
		int sum = numIt[0];
		for (int i = 1; i < numIt.length; i++) {
			sum += numIt[i];
		}
		return sum;
	}

	/*
	 * public static double trainThresh() { return trainThresh; }
	 */

	/* minimum number */
	public static int modifyWithData() {
		// TODO Auto-generated method stub
		return modifyWithData;
	}

	public static String bin() {
		return bin;
	}

	public static int index() {
		// TODO Auto-generated method stub
		return index;
	}

	public static boolean trans1 = true; // whether to use
											// BeetweenWithinTransitionProbs1

	public static boolean trans1() {
		// TODO Auto-generated method stub
		return trans1;
	}

	/*
	 * public static double trainZeroVariance = -1; //-1 means off; 0 means on
	 * but global; anything else is local public static int trainZeroVariance()
	 * { // TODO Auto-generated method stub return (int) trainZeroVariance; }
	 */

	public static String[] format() {
		return format;
	}

	public static String chrom0() {
		// if(chrom!=null) return chrom;
		return mid[0][0];
	}

	public static String chrom(int i) {
		return chrom == null ? chrom0() : chrom[i];
	}

	/*
	 * public static String[] chrom() { // TODO Auto-generated method stub //
	 * String[] chr = chrom; return chrom; }
	 */
	public static Integer end() {
		return end;
	}

	public static int sum(int[] no_copies) {
		int sm = 0;
		for (int i = 0; i < no_copies.length; i++) {
			sm += no_copies[i];
		}
		return sm;
	}

	public static int product(int[] emStSp) {
		int res = emStSp[0];
		for (int i = 1; i < emStSp.length; i++) {
			res *= emStSp[i];
		}
		return res;
	}

	/**
	 * returns the init u for transitions out of begin state public static
	 * double initTransU() { return 100; }
	 */

	public static double bg0 = 1.0;

	public static double bg1 = 1.0;

	// new double[] {2.5e-3,2.5e-3, 0.99, 2.5e-3,2.5e-3};
	public static double[] bg(int i) {
		double bg_ = i == 0 ? bg0 : bg1;
		double[] res = new double[5];
		Arrays.fill(res, (1.0 - bg_) / (double) res.length);
		res[2] += bg_;
		return res;

	}

	public static double bg = 1.0;
	public static int bgtype = 0; // 0 normal emissions, 1 hwe only, 2 both

	public static int bgtype() {
		return bgtype;
	}

	public static double bg() {
		return bg;
	}

	private static double fixedThresh = 1.0;

	public static double fixedThresh() {
		return fixedThresh;
	}

	public static boolean useFree = true;

	public static boolean useFree() {
		// TODO Auto-generated method stub
		return useFree;
	}

	public static int[] getMax2(double[] emiss) {
		int i = Constants.getMax(emiss);
		double max = Double.NEGATIVE_INFINITY;
		int r = -1;
		for (int j = 0; j < emiss.length; j++) {
			if (j == i)
				continue;
			if (emiss[j] > max) {
				max = emiss[j];
				r = j;
			}
		}
		return new int[] { i, r };
	}
	
	public static int getMax2(double[] emiss, int i) {
		
		double max = Double.NEGATIVE_INFINITY;
		int r = -1;
		for (int j = 0; j < emiss.length; j++) {
			if (j == i)
				continue;
			if (emiss[j] > max) {
				max = emiss[j];
				r = j;
			}
		}
		return r;
	}

	public static boolean illBgAvgOfFg = true;

	public static boolean illuminaBgIsAverageOfFg() {
		// TODO Auto-generated method stub
		return illBgAvgOfFg;
	}

	private static double trainThresh = 1e-3;

	public static double trainThresh() {
		// TODO Auto-generated method stub
		return trainThresh;
	}

	public static String[][] allowSharedIds = null;

	public static String[][] allowOverlaps() {
		return allowSharedIds;
	}

	public static double countThresh =1e-3;
	public  static double countThresh1 = 1e-3;
	public  static double countThresh2 = 1e-3;

	public static double countThresh() {
		// TODO Auto-generated method stub
		return countThresh;
	}
	public static double countThresh1() {
		// TODO Auto-generated method stub
		return countThresh1;
	
	}
	public static int numToSample = 1; //should be set to one for regression
	public static int numToSample(){
		return numToSample;
	}

	public static int radix() {
		return 32;
	}

	public static int[] restrictKb(int i) {
		// TODO Auto-generated method stub
		if(restrictKb==null ) return new int[] {0,0};
		String[] restrictKb1 = i < Constants.restrictKb.length ? Constants.restrictKb[i] : 
			Constants.restrictKb[0];
		int[] res = new int[2];
	
		res[0] = convert(restrictKb1[0]) ;
		res[1] = convert(restrictKb1[1]) ;
		return res;
	}
	
	public static int[] restrictKbMax() {
		// TODO Auto-generated method stub
		int[] res = new int[] {0,0};
		if(restrictKb==null) return res;
		for(int k=0; k<restrictKb.length; k++){
			int[] r1 = restrictKb(k);
		
			if(r1[0] >res[0]) res[0] = r1[0];
			if(r1[1] > res[1]) res[1] = r1[1];
		}
		return res;
	}

	public static boolean reverse = false;

	public static boolean reverse() {
		return reverse;
	}

	public static boolean drop = false;

	public static boolean drop() {
		// TODO Auto-generated method stub
		return drop;
	}

	public static boolean trainCGH() {
		return trainCGH;
	}

	private static double bwThresh = 1e-7;

	public static double bwThresh() {
		// TODO Auto-generated method stub
		return bwThresh;
	}

	public static boolean noHMM() {
		return false;
	}

	public static boolean run = true;
	public static boolean no_hmm = false; // option to run to get max likelihood
											// values without hmm based training

	public static boolean no_hmm() {
		return no_hmm;
	}

	public static boolean run() {
		// TODO Auto-generated method stub
		return run;
	}

	public static String column() {
		return "_" + column;
	}

	public static CommandLine parseInitial(String[] args) throws Exception { 
		for (int i = 0; i < args.length; i++) {
			args[i] = check(args[i]);
		}
		Parser parser = new PosixParser();
		// Integer[] cols = para.getOptionValue("column");
		final CommandLine params = parser.parse(Constants.OPTIONS, args, false);

		return params;
	}

	private static String check(String string) {
		if (string.indexOf('-') < 0 || string.indexOf("--") >= 0)
			return string;
		StringBuffer sb = new StringBuffer();
		String[] str = string.split(":" + "" + "");
		for (int i = 0; i < str.length; i++) {
			String[] str1 = str[i].split("-");
			int st = Integer.parseInt(str1[0]);
			int end = Integer.parseInt(str1[1]);
			for (int k = st; k <= end; k++) {
				sb.append(k);
				sb.append(":");
			}
		}
		String st1 = sb.toString();
		return st1.substring(0, st1.length() - 1);
	}

	public static String[] getCols(String[] args) throws Exception {
		final CommandLine params = parseInitial(args);
		String[] cols = new String[] { 1 + "" };
		if (params.hasOption("column")) {
			cols = params.getOptionValues("column");
		}
		return cols;
	}

	public static double[] softenHapMap = new double[] {0.0};

	public static double softenHapMap(int i) {
		// TODO Auto-generated method stub
		return softenHapMap[i];
	}

	public static double cn_ratio = 1.0;

	public static double cn_ratio() {
		return cn_ratio;
	}

	// public static int[] specialTrans = new int[0];
	// public static int[] specialTrans() {
	// TODO Auto-generated method stub
	// return specialTrans;
	// }
	// public static String[] codesToRestrict = new String[0];
	/*
	 * public static Collection<Character> codesToRestrict() { String[] codes =
	 * codesToRestrict; Set<Character> res = new HashSet<Character>(); for(int
	 * i=0; i<codes.length; i++){ res.add(codes[i].charAt(0)); } return res; }
	 */
	public static void delete(File submissions) {
		if (submissions.isDirectory()) {
			File[] f = submissions.listFiles();
			for (int i = 0; i < f.length; i++) {
				delete(f[i]);
			}
			submissions.delete();
		} else {
			submissions.delete();
		}

	}

	public static int[] backgroundCount = new int[] { 2 };

	public static int backgroundCount(int i) {
		// TODO Auto-generated method stub
		return backgroundCount[i];
	}

	public static boolean trainData = false;

	public static boolean trainData() {
		// TODO Auto-generated method stub
		return trainData;
	}

	public static boolean saveStatePath = false;
	public static boolean saveStates() {
		return saveStatePath || Constants.plot() > 0 || Constants.calcAssoc;
	}
	public static boolean saveMostLikely(){
		return saveMostLikely;
	}

	// public static boolean allowCloning = false;
	// public static boolean allowCloning() {
	// TODO Auto-generated method stub
	// return allowCloning;
	// }

	public static PseudoIterator weights__ = null;

	public static PseudoIterator weights() {
		// TODO Auto-generated method stub
		return weights__;
	}

	public static double round() {
		return round;
	}

	public static int plot = 0;
	public static boolean[] plotI = new boolean[0];

	public static int plot() {
		return plot;
	}

	public static boolean plot(int i) {
		if (i < plotI.length) {
			return plotI[i];
		} else
			return true;
	}

	public static boolean pause = false;

	public static boolean pause() {
		return pause;
	}

	public static int trainEnsemble = 2; // 0 - no ensemble 1 - additive 2- exp

	public static int trainEnsemble() {
		if (true)
			throw new RuntimeException("!!");
		return -1;
	}

	/** modifies the variance by data type */
	// public static double[] r_state_var_mod = new double[] {1.0};
	public static int[] cnStatesInCommonClass = new int[] { 1 };

	public static int[] cnStatesInCommonClass() {
	//	if (true)
		//	throw new RuntimeException("!!");
		return cnStatesInCommonClass;
	}

	public static double base = 2;// Math.exp(1);

	public static double base() {
		return base;
	}

	public static int numThreads = 0;

	public static int numThreads() {
		return numThreads;
	}

	public static String[] toDel = new String[0];
	public static String[] toInclude = new String[0];
	public static String[] toReport = null;

	public static boolean toReport(int i) {
		return Boolean.parseBoolean(toReport[i]);
	}

	public static boolean dropFixed = false;

	public static boolean dropFixed() {
		return dropFixed;
	}

	public static Set<Integer> reportIds() {
		if (toReport == null)
			return null;
		String[] toR = toReport;
		Set<Integer> s = new HashSet<Integer>();
		for (int i = 0; i < toReport.length; i++) {
			if (Boolean.parseBoolean(toReport[i])) {
				s.add(i);
			}
		}
		return s;
	}

	public static String[] rOutlier = new String[0];

	public static double[] rOutlier(int i) {
		double[] res = null;
		if (i >= rOutlier.length || rOutlier[i].equals("null"))
			return null;
		String[] str = rOutlier[i].split(";");
		res = new double[str.length];
		for (int i1 = 0; i1 < res.length; i1++) {
			res[i1] = Double.parseDouble(str[i1]);
		}
		return res;
	}

	public static String toDel(int i) {
		// TODO Auto-generated method stub
		return toDel[Math.min(toDel.length-1,i)].replaceAll("\\s+", "");
	}

	public static String toInclude(int i) {
		// TODO Auto-generated method stub
		String res =  toInclude[i>=toInclude.length ? toInclude.length-1:i].replaceAll("\\s+", "");
		
		return res;
	}
	public static List<String> read(List<String> l1, int index) {
		try{
		if(l1.size()==1 && l1.get(0).endsWith(".txt")){
		
		List<String> l = new ArrayList();
		File loc = new File(Constants.inputDir(index));
		File inp = new File(loc, l1.get(0));
		if(!inp.exists()){
			inp = new File(l1.get(0));
		}
		BufferedReader br = new BufferedReader(new FileReader(inp));
		String st = "";
		while((st = br.readLine())!=null){
			if(!st.startsWith("#")) l.add(st.split("\t")[0]);
		}
		return l;
		}
		}catch (Exception exc){
			exc.printStackTrace();
		}
		return l1;
	}
	public static List<String[]> read1(List<String[]> l1) {
		try{
		if(l1.size()==1 && l1.get(0)!=null && l1.get(0)[0].endsWith(".txt")){
		
		List<String[]> l = new ArrayList();
		BufferedReader br = new BufferedReader(new FileReader(new File(l1.get(0)[0])));
		String st = "";
		while((st = br.readLine())!=null){
			l.add(st.split("\t"));
		}
		return l;
		}
		}catch (Exception exc){
			exc.printStackTrace();
		}
		return l1;
	}


	public static boolean supressR = false;

	public static boolean suppressR() {
		return supressR;
	}

	public static boolean supressB = false;

	public static boolean suppressB() {
		return supressB;
	}
	
	public static String[] homDelBdist = new String[] {"U0-1"};
	
	public static String homDelBdist(int i){
		return homDelBdist[Math.min(homDelBdist.length-1, i)];
	}

	public static boolean logplot = false;

	public static boolean logplot() {
		return logplot;
	}

	public static boolean[] indexToRestrict = null;

	public static boolean[] indexToRestrict() {
		return indexToRestrict;
	}

	public static boolean[] loess = new boolean[] { true };
	public static boolean[] median_correction = new boolean[] { true };

	public static boolean loess(int i) {
		if(loess==null) return false;
		if (i >= loess.length) {
			return loess[0];
		}
		return loess[i];
	}

	public static boolean[] gc = new boolean[] { true };

	public static boolean gc(int i) {
		if(gc==null) return false;
		if (i >= gc.length) {
			return gc[0];
		}
		return gc[i];
	}

	public static boolean median_correction(int i) {
		if (i >= median_correction.length)
			return median_correction[0];
		return median_correction[i];
	}

	// public static int[] numLevels =new int[0]; //a numLevel of zero indicates
	// continuous variable
	public static int segments() {
		return segments;
	}

	public static final boolean scoreDT = false;
	public static final boolean countDT = false;
	public static final FileFilter ZIP_FILTER = new FileFilter() {

		public boolean accept(File pathname) {
			return pathname.getName().endsWith(".zip");
		}

	};

	public static boolean writeExtractedFile() {
		// TODO Auto-generated method stub
		return writeExtractedFile;
	}

	public static double var_thresh(int i) {
		// System.err.println("var thresh "+var_thresh);
		if (i < var_thresh.length)
			return var_thresh[i];
		else {
			return var_thresh[0];
		}
	}

	public static String print(double[] distribution) {

		StringBuffer sb = new StringBuffer("");
		sb.append(String.format("%5.3g", distribution[0]));
		for (int i = 1; i < distribution.length; i++) {

			sb.append(String.format("\t%5.3g", distribution[i]));
		}
		return sb.toString();
	}

	public static String print(int[] distribution) {

		StringBuffer sb = new StringBuffer("");
		sb.append(distribution[0]);
		for (int i = 1; i < distribution.length; i++) {

			sb.append("\t" + distribution[i]);

		}
		return sb.toString();
	}

/*	public static double exponentB(int i) {
		if(modExponent(i)) return 1-fillLikelihood1(i);
		else return 1.0;
		//return exponentB[i];
	}*/

	public static double bRandom = 1.0;

	public static double bRandom() {
		return bRandom;
	}

	public static double exponentR(int i) {
		if(true) throw new RuntimeException("!!");
		return exponentR[i];
	}

	public static String[] intensityOnlyIds = new String[] { "cnv", "A_" };

	public static String[] intensityOnly() {
		return intensityOnlyIds;
	}

	//public static boolean printPlots = false;

	public static boolean printPlots() {
		return Constants.plot==1;
	}
public static boolean showScatter = true;
public static boolean showScatter(){
	return showScatter;
}
	public static String[] pdf = new String[] { "pdf", "pdf" };// plot,
																// haplotype

	public static String pdf(int i) {

		return pdf[(int) Math.min(pdf.length - 1, i)];
	}

	public static double showEnsThresh = 500000; // only show if distance from
													// end to beginning is less
													// than this, set 0 to turn
													// off

	public static double showEnsThresh() {
		// TODO Auto-generated method stub
		return showEnsThresh;
	}

	public static int indexLength = 2;

	// length of index to return for data index track in haplotype panel;
	public static int indexLength() {
		return indexLength;
	}

	public static String restrictIndivTo = "all";
	public static String restrictIndivTo1 = null;
	public static String restrictIndivTo1(){
		return restrictIndivTo1;
	}
	public static String restrictIndivTo(String[] inp) {
		// TODO Auto-generated method stub
		if(restrictIndivTo==null || restrictIndivTo.equals("null") || restrictIndivTo.equals("all")) return null;
		String st = restrictIndivTo.replaceAll("\\|\\|", "\\|");
		st = st.replaceAll("\\&\\&", "\\&").replace(':', '&');
		
		for(int i=0; i<inp.length; i++){
			st = st.replace(inp[i], i+"");
		}
		st = "("+st+");";
		int i=0;
		Comparator comp = new Comparator<Integer>(){

			@Override
			public int compare(Integer o1, Integer o2) {
				return -1*o1.compareTo(o2);
			}
			
		};
		SortedMap<Integer, String> m = new TreeMap<Integer, String>(comp);
		while(i>=0){
			int i1 = st.indexOf('|', i+1);
			int i2 = st.indexOf('&', i+1);
			i = i1<0 ? i2 : (i2< 0 ? i1 : Math.min(i1,i2));
			if(i>=0){int j = findClosingBracket(st,i+1);
			m.put(j+1, i==i1 ? ":1" : ":0");
			}
		}
		StringBuffer sb= new StringBuffer(st);
		for(Iterator<Integer> it = m.keySet().iterator(); it.hasNext();){
			Integer key = it.next();
			String v = m.get(key);
			sb.insert(key, v);
		}
		String res = sb.toString().replace('&', ',').replace('|', ',');
		return res;
	}

	private static int findClosingBracket(String st2, int i) {
		int cnt =-1;
		for(int k=i; k<st2.length(); k++){
			if(st2.charAt(k)==')') cnt++;
			else if(st2.charAt(k)=='(') cnt--;
			if(cnt==0) return k;
		}
		return -1;
	}
	public static Boolean writeRes =null;  //if null, write nothing, if true do states

	public static Boolean writeRes() {
		// TODO Auto-generated method stub
		return writeRes;
	}

	public static boolean demergeToPrint = false;

	public static boolean demergeToPrint() {

		return demergeToPrint;
	}

	public static String baseDir() {
		return baseDir;
	}
	
	public static String suffix(short i) {
		StringBuffer suff = new StringBuffer();
		if(!Constants.toInclude(i).equals("all")){
		    suff.append("_i_"+Constants.toInclude(i).length());
	    }
	   if((!Constants.platesToInclude(i)[0][0].equals("all"))){
		suff.append("_i_"+Arrays.asList(Constants.platesToInclude(i)[0]));
	   }
	   if(Constants.phenoToInclude(i)[0]!=null && (!Constants.phenoToInclude(i)[0][0].equals("all"))){
		 suff.append("_i_"+Arrays.asList(Constants.phenoToInclude(i)[0]));
	  }
	   if(!Constants.toDel(i).equals("null")){
		    suff.append("_e_"+Constants.toDel(i));
	    }
	   if((!Constants.platesToExclude(i)[0][0].equals("null"))){
		suff.append("_e_"+Arrays.asList(Constants.platesToExclude(i)[0]));
	   }
	   if((!Constants.phenoToExclude(i)[0][0].equals("null"))){
		 suff.append("_e_+"+Arrays.asList(Constants.phenoToExclude(i)[0]));
	  }
		return suff.toString();
	}

	public static String getDirName(File outp1, String loc_string) {
		
		
		
		
		File nme = new File(outp1, Constants.experiment[Math.min(
				Constants.experiment.length - 1, Constants.column - 1)]
				+ "_" + Constants.column + "_" + loc_string);
		if (!Constants.overwrite()) {
			for (int i = 0; nme.exists(); i++) {
				nme = new File(outp1, Constants.experiment[Math.min(
						Constants.experiment.length - 1, Constants.column - 1)]
						+ "_" + Constants.column + "_" + loc_string+"_" + i);
			}
		} else if (nme.exists()) {
			delete(nme);
		}
		return nme.getName();
	}

	public static Map<String, Integer> readReps(File file) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String st = "";
		Map<String, Integer> m = new HashMap<String, Integer>();
		while ((st = br.readLine()) != null) {
			if (st.startsWith("#") || st.length() < 2)
				continue;
			String[] str = st.split("\\s+");
			String str1 = str[0].substring(2);
			Integer value = m.get(str1);
			m.put(str1, value == null ? 1 : value + 1);
		}
		for (Iterator<Map.Entry<String, Integer>> it = m.entrySet().iterator(); it
				.hasNext();) {
			Map.Entry<String, Integer> st_ = it.next();
			if (st_.getValue() == 1) {
				it.remove();
			}
		}
		return m;

	}

	public static int[] expand = new int[] { 0 };
	public static double[] expand_init_prior = new double[] { 1000,10 }; // emiss,
																		// trans

	public static int[] expand() {
		return expand;
	}

	public static int[] group = null;

	public static int[] group() {
		// TODO Auto-generated method stub
		return group;
	}

	public static double expand_init_prior(int i) {
		// TODO Auto-generated method stub
		if(i>=expand_init_prior.length) return expand_init_prior[expand_init_prior.length-1];
		return expand_init_prior[i];
	}

	public static void normalise(double[] init_B) {
		double sum = Constants.sum(init_B);
		for (int i = 0; i < init_B.length; i++) {
			init_B[i] = init_B[i] / sum;
		}

	}

	public static double sum(double[] probs, int[] b) {
		double sum = 0;
		for (int i = 0; i < probs.length; i++) {
			sum += probs[i] * b[i];
		}
		return sum;
	}

	public static boolean inversion = false;

	public static boolean inversion() {
		// TODO Auto-generated method stub
		return inversion;
	}

	public static boolean showHaps = false;

	public static boolean showHaps() {
		return showHaps;
	}

	public static double[] soften = new double[] { 0 }; // same
														// direction/normal/opposite
														// direction

	public static double soften(int i) {
		// TODO Auto-generated method stub
		return soften[i];
	}

	public static String[][] restrictMarker = new String[0][0];

	public static String[] restrictMarker(int i) {
		if (i < restrictMarker.length)
			return restrictMarker[i];
		else
			return new String[] { "null" };
	}

	public static String[][] plotType = new String[0][0];

	public static String[] plotType(int i) {
		if (i < plotType.length)
			return plotType[i];
		else
			return new String[] { "XYLineAndShapeRenderer", "dot", "noline" };
	}

	public static int distanceTo = 0;

	public static int distanceTo() {
		// TODO Auto-generated method stub
		return distanceTo;
	}

	public static boolean calcLD = false;
	// type 0==no_copies; type 1 noA; type 2 noB ;type==3 all states type ==4
	// intensity type 5 = BAF; type 6 = snp (excluding cn)
	public static  String[] ldtype = new String[]{ "0", "6"} ;

	public static String[] ld() {
		return ldtype;
	}

	public static boolean calcLD() {
		// TODO Auto-generated method stub
		return calcLD;
	}

	public static String snpsToPlot = "all";

	/*public static String[] snpsToPlot() {
		if (snpsToPlot.equals("all")) {
			return null;
		}else if (snpsToPlot.equals("none")){
			return new String[1];
		}
		else
			return new String[Constants.expandScatter * 2 + 1];
	}*/

	public static boolean addAnnotationToDistGraphs() {
		// TODO Auto-generated method stub
		return true;
	}

	public static double minR(int i) {
		if(i<minR.length) return minR[i] - 0.01;
		else return minR[0] - 0.01;
	}

	public static double minB(int i) {
		if(i<minB.length) return minB[i] ;
		else return minB[0];
	}

	public static double maxB(int i) {
		if(i<maxB.length) return maxB[i];
		else return maxB[0];
	}

	public static double[] minB = new double[] {0.0};
	public static double[] maxB = new double[] {1.0};
	// public static double minR=-8;
	// public static double maxR = 8;
	public static double[] mixCoeff = new double[] { 1e-5 };

	// public static double[] mixCoeffB = new double[] {1e-3};
	public static double maxR(int i) {
		if(i<maxR.length) return maxR[i] + 0.01;
		else return maxR[0] +0.01;
	}

	public static double mixCoeff(int i) {
		// TODO Auto-generated method stub
		return mixCoeff[Math.min(i, mixCoeff.length - 1)];
	}

	public static double plotThresh = 1e-7;
	public static double state_brightness = 0.33;
	public static int[] rs = new int[0];

	public static int[] rs() {
		return rs;
	}

	public static double plotThresh() {
		// TODO Auto-generated method stub
		return plotThresh;
	}

	public static boolean writeFree() {
		// TODO Auto-generated method stub
		return false;
	}

	// public static double mixCoeffB(int i) {
	// return mixCoeffB[Math.min(i, mixCoeffB.length-1)];
	// }
	public static int getPreferredIndex() {
		return 1;
	}

	public static double[] weight = new double[] { 1 };

	public static double[] weight() {
		return weight;
	}

	public static String[] writeAverages = null; // "countA", "countB",
															// "data" is null
															// not write
	static List<String> before = Arrays.asList(new String[] { "BAF", "logR" });

	public static String[] writeAverages(boolean bef) {
		if(writeAverages==null) return null;
		String[] res = writeAverages;
		List<String> res1 = new ArrayList<String>();
		for (int i = 0; i < writeAverages.length; i++) {
			if (bef && before.contains(writeAverages[i])) {
				res1.add(writeAverages[i]);
			} else if (!bef && !before.contains(writeAverages[i])) {
				res1.add(writeAverages[i]);
			}
		}
		return res1.toArray(new String[0]);
	}

	public static String[] annotateSamples = null;

	public static String[] annotateSamples() {
		return annotateSamples;
	}

	public static double median(double[] mcs, double len) {
		Arrays.sort(mcs);
		double midp = (len - 1.0) / 2.0;
		int mid = (int) Math.floor(midp);
		if (Math.IEEEremainder(len - 1, 2) == 0) {

			return mcs[mid];
		} else {
			return (mcs[mid] + mcs[mid + 1]) / 2.0;
		}
	}

	public static boolean[] standardVar = null;

	public static boolean standardiseVariance(short index2) {
		if (standardVar == null)
			return false;
		else
			return standardVar[index2];
	}

	public static double[] imputedThresh = new double[] { 0.9, 0.9 };
	public static double[] imputedThreshGraph = new double[] { 0.0, 0.0 };

	public static double imputedThresh(int i) {
		if (i < imputedThresh.length)
			return imputedThresh[i];
		else
			return imputedThresh[0];
	}
	
	public static double imputedThreshGraph(int i) {
		if (i < imputedThreshGraph.length)
			return imputedThreshGraph[i];
		else
			return imputedThreshGraph[0];
	}

	public static double ld_dist_thresh = 200 * 1000;

	public static double ld_dist_thresh() {
		return ld_dist_thresh;
	}

	public static double ld_r2_thresh =1.01;

	public static double ld_r2_thresh() {
		return ld_r2_thresh;
	}

	public static boolean likelihoodInput = false;
	public static boolean writeMergedAverageFiles = true;
	public static String[] splitbyPheno;
	public static double[] minR = new double[] {-8};
	public static double[] maxR = new double[] {5};

	public static boolean likelihoodInput() {
		// TODO Auto-generated method stub
		return likelihoodInput;
	}

	public static String experiment() {
		return experiment[0];
	}

	public static boolean writeMergedAverageFiles() {
		return Constants.writeMergedAverageFiles;
	}

	public static String splitByPheno(int i) {
		// TODO Auto-generated method stub
		if (i >= splitbyPheno.length || splitbyPheno[i].length() == 0
				|| splitbyPheno[i].equals("null")) {
			return null;
		} else {
			String res = splitbyPheno[i];
			return res;
		}
	}

	public static boolean ascn() {
		return true;
	}

	public static int scatterWidth() {
		return scatterWidth;
	}
	public static int scatterWidth=500;

	public static int expandScatter = 5;

	public static int expandScatter() {
		return expandScatter;
	}

	public static boolean fixedBG = false;

	public static boolean fixedBG() {
		return fixedBG;
	}

	public static boolean includeLegend = false;

	public static boolean includeLegend() {
		return includeLegend;
	}

	public static boolean calcAssoc = true;

	public static boolean calcAssoc() {
		// TODO Auto-generated method stub
		return calcAssoc;
	}

	public static boolean modelbg = false;

	public static boolean modelbg() {
		// TODO Auto-generated method stub
		return modelbg;
	}

	public static double hapl_cert_thresh = 0.0;

	public static double hapl_cert_thresh() {
		// TODO Auto-generated method stub
		return hapl_cert_thresh;
	}

	public static boolean trainSkew = true;

	public static boolean trainSkew() {
		// TODO Auto-generated method stub
		return trainSkew;
	}

	public static double haplotypeHeight = 10;

	public static double haplotypeHeight() {
		return haplotypeHeight;
	}

	public static boolean showHMM = false;

	public static boolean showHMM() {
		return showHMM;
	}

	public static double annotationTail = 0.2;

	public static double[] annotationInterval() {
		// TODO Auto-generated method stub
		return new double[] { annotationTail, 0.5, 1 - annotationTail };
	}

	public static boolean sepModels = false;

	public static boolean sepModels() {
		// TODO Auto-generated method stub
		return sepModels;
	}

	public static boolean trainTransitions = true;
	public static boolean trainDists = true;
	public static boolean trainDists(){
		return trainDists;
	}
	public static boolean trainTransitions() {
		// TODO Auto-generated method stub
		return trainTransitions;
	}

	public static boolean[] showAll = new boolean[] { false, false};

	public static boolean showAll(int i) {
		// TODO Auto-generated method stub
		return i>=showAll.length ? showAll[0] : showAll[i];
	}

	public static boolean halfNormal = false;

	public static boolean halfNormal() {
		// TODO Auto-generated method stub
		return halfNormal;
	}

	public static double shapeSize = 10;

	public static double shapeSize() {
		// TODO Auto-generated method stub
		return shapeSize;
	}
	public static double shapeSize1 = 5; //scaling the font

	public static double shapeSize1() {
		// TODO Auto-generated method stub
		return shapeSize1;
	}

	public static boolean export = false;

	public static boolean export() {
		// TODO Auto-generated method stub
		return export;
	}

	public static final boolean joint = true;
	public static String[] toInclude1 = new String[0];
	
	// public static boolean joint() {
	// // TODO Auto-generated method stub
	// return true;/
	// }

	public static int[] maxG = new int[] { 3 };

	public static int maxG(int i) {
		// TODO Auto-generated method stub
		if (i >= maxG.length)
			return maxG[maxG.length - 1];
		return maxG[i];
	}

	public static int[][] transitions2;
	public static int[][] transitions1;

	public static int[][] transitions(int i) {
		if (i == 0)
			return transitions1;
		return transitions2;
	}

	public static boolean orthogonal = true;

	public static boolean orthogonal() {
		// TODO Auto-generated method stub
		return orthogonal;
	}

	

	public static boolean convertIds = false;

	public static boolean convertIds() {
		return convertIds;
	}

	/*public static double[][] b_mean1(int i) {
		return b_mean[i>=b_mean.length ? 0 : i];
	}*/

	public static boolean b_panel = true;
	public static boolean r_panel = true;

	public static boolean b_panel() {
		return b_panel;
	}

	public static boolean r_panel() {
		return r_panel;
	}

/*	public static double[] b_var(int i, boolean first) {
		double[][] b_var = first? b_var1 : b_var4;
		if (b_var[i].length < b_mean[i].length) {
			double[] newr = new double[b_mean[i].length];
			System.arraycopy(b_var[i], 0, newr, 0, b_var[i].length);
			for (int k = b_var[i].length; k < newr.length; k++) {
				newr[k] = b_var[i][b_var[i].length - 1];
			}
			b_var[i] = newr;
		}
		return b_var[i];
	}*/

	public static double[] b_skew(int i) {
		if (b_skew[i].length < b_mean[i].length) {
			double[] newr = new double[b_mean[i].length];
			System.arraycopy(b_skew[i], 0, newr, 0, b_skew[i].length);
			for (int k = b_skew[i].length; k < newr.length; k++) {
				newr[k] = b_skew[i][b_skew[i].length - 1];
			}
			b_skew[i] = newr;
		}
		return b_skew[i];
	}

	public static int[] regress = new int[0];

	public static int[] regress() {
		// if(regress.length==1 && regress[0]==0) return new int[0];
		return regress;
	}

	public static int[] mltrain = new int[0];

	public static int[] mltrain() {
		// if(mltrain.length==1 && mltrain[0]==0) return new int[0];
		return mltrain;
	}



	public static double rhoMax = 0.5;

	public static double rhoMax() {
		// TODO Auto-generated method stub
		return rhoMax;
	}

	public static double getMinValue(double[][] r_train4) {
		double d = Double.POSITIVE_INFINITY;
		for (int i = 0; i < r_train4.length; i++) {
			double min = r_train4[i][Constants.getMin(r_train4[i])];
			if (min < d) {
				d = min;
			}
		}
		return d;
	}

	public static int[] stop = new int[0];
	private static Set<Integer> st = null;

	public static boolean stop(int i1) {
		if (st == null) {
			st = new HashSet<Integer>();
			for (int i = 0; i < stop.length; i++) {
				st.add(stop[i]);
			}
		}
		return st.contains(i1);
	}

	public static double getMinValue(double[] rho_train3) {
		return rho_train3[getMin(rho_train3)];
	}

	public static String experimentPhenoFile = "pheno_all.txt";

	public static String experimentPhenoFile() {
		return experimentPhenoFile;
	}

	public static String[][] phenoToAssoc = new String[][] {{"null"}};

	public static String[] phenoToAssoc(int i) {
		String[] res = i<0 || phenoToAssoc.length<=i ? phenoToAssoc[0] : phenoToAssoc[i];
		if (res.length == 1
				&& res[0].equals("null"))
			return null;
		return res;
	}

	public static String[][] quantiles = null;

	public static PhenoGroup[] quantiles(int i) {
		if(quantiles==null || i>=quantiles.length || quantiles[i]==null || quantiles[i][0].equals("null")) return null;
		PhenoGroup[] res = new PhenoGroup[quantiles[i].length];
		for(int k=0; k<res.length; k++){
			res[k] = new PhenoGroup(quantiles[i][k]);
		}
		return res;
	}
public static boolean saveModel = true;;

	public static boolean saveModel() {
		// TODO Auto-generated method stub
		return saveModel;
	}

	public static String hmmFile[] = null;

	public static String useHMMFile(int i) {
		// TODO Auto-generated method stub
		if (hmmFile == null || hmmFile.length<=i || hmmFile[i].equals("null"))
			return null;
		return hmmFile[i];
	}

	public static int[] transferHMMStates = null;

	public static int[] transferHMM() {
		// TODO Auto-generated method stub
		return transferHMMStates;
	}

	public static int maxPloidy1 = 1;

	public static int maxPloidy() {
		int maxp =  Math.max(maxPloidy1(), Constants.noCopies[Constants.getMax(Constants.noCopies)]);
		
		return maxp;
	}

	public static double product(double[] sc) {
		double res = sc[0];
		for (int i = 1; i < sc.length; i++) {
			res *= sc[i];
		}
		return res;
	}

	// public static double mixPrior = 1.0;

	public static double[] skewTransitions = new double[] { 0.1, 0.1 };

	public static double skewTransitions(int i) {
		// TODO Auto-generated method stub
		return skewTransitions[i];
	}

	public static String[] toDropSnpId = new String[0];

	public static String[] toDropPrefix() {
		// TODO Auto-generated method stub
		return toDropSnpId;
	}

	public static String[][] toJoin = null;
	public static String rsid;

	public static String rsId() {
		return rsid;
	}

	public static String[] toJoin(int i) {
		// TODO Auto-generated method stub
		if (toJoin == null || toJoin[i] == null || toJoin[i].length == 0
				|| toJoin[i][0].equals("null"))
			return null;

		return toJoin[i];
	}

public static int assocTest = 0;// linear armitage chisq

	public static int assocTest() {
		return assocTest;
	}

	public static double qualityThresh = 0.05;

	public static double qualityThresh() {
		return qualityThresh;
	}

	public static String rsid() {
		// TODO Auto-generated method stub
		return null;
	}

	public static boolean haploidMale = false;
	public static String splString = "\t";
	
	public static String splString(){
		return splString;
	}

	public static boolean haploidMale() {
		// TODO Auto-generated method stub
		return haploidMale;
	}

	public static double priorOdds;

	public static double priorOdds() {
		// TODO Auto-generated method stub
		return priorOdds;
	}

	

	

	public static int getFinalModelLength() {
		int[] numIt = Constants.numIt;
		if (numIt.length == 1)
			return Constants.modify(0).length;
		else {
			int num = 0;
			for (int i = 0; i < Constants.modify(0).length; i++) {
				int nocop = Integer.parseInt(""+Constants.modify(0)[i]);
				num += getNoCop(Constants.expand[i],nocop);
			}
			return num;
		}

	}

	private static int getNoCop(int i, int nocop) {
		if(nocop<=1) return i;
		else if(nocop==2){
			double i1 = i;
			return (int) Math.round(i1*(i1+1)/2.0);
		}
		else if(nocop==3){
			if(i==1) return 1;
			else if(i==2) return 4;
			else if(i>3) throw new RuntimeException("still need to work out formula"); 
		}
		else if(i==1) return 1;
		else if(i==2) return nocop+1;
		else 	throw new RuntimeException("still need to work out formula"); 
		return -1;
	}

	public static String stratify = null;
	public static String[] restrictPheno = null;

	public static int scatterPlotsPerPage = 20;

	public static String stratify() {
		return stratify;
	}

	public static int scatterPlotsPerPage() {
		// TODO Auto-generated method stub
		return scatterPlotsPerPage;
	}

	public static int maxSNPS = 20;

	public static int maxSNPS() {
		return maxSNPS;
	}

	public static boolean replaceAB = false;

	public static boolean replaceAB() {
		return replaceAB;
	}

	public static String restrictPheno(short index2) {
		if (restrictPheno == null || restrictPheno.length == 0
				|| restrictPheno[index2] == null
				|| restrictPheno[index2].equals("null"))
			return null;
		else
			return restrictPheno[index2];
	}

	public static int[] cnToAnnotate = null;
	
	public static int[] cnToAnnotate() {
		return cnToAnnotate;
	}
public static int[] bafToAnnotate = null;
	
	public static int[] bafToAnnotate() {
		return bafToAnnotate;
	}

	public static boolean singleBClust = false;

	public static boolean singleBClust() {
		return singleBClust;
	}

	public static double minBafVar = 0.005;
	public static double minLRRVar = 0.03;

	public static double minBafVar() {
		// TODO Auto-generated method stub
		return minBafVar;
	}

	public static double minLRRVar() {
		// TODO Auto-generated method stub
		return minBafVar;
	}

	public static int changePriorCount = 2;

	public static int changePriorCount() {
		return changePriorCount;
	}
public static boolean phaseInner = false;
	public static boolean phaseInner() {
		// TODO Auto-generated method stub
		return phaseInner;
	}
public static double maxLogP = 30;
	public static double maxLogP() {
		return maxLogP;
	}

	

	public static String[] emissionGroup() {
		if(emissionGroup==null)
		{
			 String[] nme = new String[inputDir.length];
		        for(int i=0; i<inputDir.length; i++){
		        	
		        	 nme[i] =inputDir(i).split("_")[0];
		        }
		        emissionGroup = nme;
		}
		 return emissionGroup;
	}
	public static String[] equaliseGroup = null;
	public static int[][] equaliseGroup(){
		if(equaliseGroup==null) return null; 
		Map<String,List< Integer>> m = new HashMap<String, List<Integer>>();
		for(int i=0; i<equaliseGroup.length; i++){
			List<Integer> l = m.get(equaliseGroup[i]);
			if(l==null){
				m.put(equaliseGroup[i], l = new ArrayList<Integer>());
			}
			l.add(i);
		}
		int[][]res   = new int[m.size()][];
		int k=0;
		for(Iterator<List<Integer>> it = m.values().iterator();it.hasNext();k++){
			List<Integer> l = it.next();
			res[k] = new int[l.size()];
			for(int j=0; j<l.size(); j++){
				res[k][j] = l.get(j);
			}
		}
		
		return res;
	}
public static int plotHeight = 600;
	public static int plotHeight() {
		return plotHeight;
	}

	public static String print(Double[] doubles) {
		return Arrays.asList(doubles).toString();
	}
public static double[] expSd = new double[] {1e-5,1e-5};
	public static double expSd(int i) {
		return expSd[i];
	}
public static double printThresh = 0.01;
	public static double printThresh() {
		// TODO Auto-generated method stub
		return printThresh;
	}

	
	public static boolean useLocInHap = false;
	public static boolean useLocInHap() {
		// TODO Auto-generated method stub
		return useLocInHap;
	}
	public static String haploImageType = "emf";
	public static String haploImageType() {
		return haploImageType;
	}
public static double haplotypeWidth = 5;

public static boolean showHapAllele = false;

public static boolean[] deletionIsNa;

public static String hapAllele="x";

public static boolean includeNullPhenotypeValues = false;

private static boolean newTrans = true;

public static boolean newTrans(){
	return newTrans;
}

//public static double[] pseudoMod;
	public static double haplotypeWidth() {
		return haplotypeWidth;
	}

	public static boolean showHapAllele() {
		// TODO Auto-generated method stub
		return showHapAllele ;
	}

	public static boolean useDeletion(int i) {
		boolean res = false;
		if(deletionIsNa==null){
			return false;
		}
		else if(i>=deletionIsNa.length) {
			res =  deletionIsNa[deletionIsNa.length-1];
		}
		res =  deletionIsNa[i];
		if(res) throw new RuntimeException("!!");
		return res;
	}

	public static char hapAllele() {
		return hapAllele.charAt(0);
	}

	public static boolean includeNullValues() {
		// TODO Auto-generated method stub
		return includeNullPhenotypeValues ;
	}

	
public static double[] modifyDirectCounts = new double[] {1.0,1.0};

	public static double[] modifyDirectCounts() {
		return modifyDirectCounts;
	}

	public static boolean ignoreHWEInData() {
		// TODO Auto-generated method stub
		return false;
	}
public static boolean trainGlobal =false;
	public static boolean measureGlobal() {
		return trainGlobal;
	}

	public static double[] range() {
		return new double[] {-1,1,-0.1,1.1};  //xmin, xmax, ymin, ymax
//		miny = -0.1; maxy = 1.1;
	}

	public static boolean regressVariance(int i) {
		// TODO Auto-generated method stub
		return regressVariance[i];
	}
	public static boolean regressMean(int i) {
		// TODO Auto-generated method stub
		return regressMean[i];
	}
	public static boolean[] regressVariance = new boolean[] {false, false}; //[r b]
	public static boolean[] regressMean = new boolean[] {true, true}; //[r b]

	public static int expandCN(char i) {
		int k=0; 
		for(;k<Constants.modify(0).length; k++){
			if(modify(0)[k]==i){
			 return Constants.expand[k];	
			}
		}
		return -1;
	}
public static double bmid=0.5;
public  static double svnTransM=0;
	public static double bmid() {
		return bmid;
	}

	public static double svnTransM() {
		// TODO Auto-generated method stub
		return svnTransM;
	}
	public static Boolean readClusters = false;

	public static boolean readClusters() {
		// TODO Auto-generated method stub
	//	if(readClusters==null) readClusters = !trainGlobal;
		return readClusters ;
	}
private static double[][][]priors;

	public static double[][]priors(int index1) {
		//if(priors==null) 
		if(priors==null){
			priors = new double[Constants.r_train.length][][];
			for(int index=0; index<priors.length; index++){
			priors[index] = new double[][]{Constants.r_train[index],Constants.r_trainVar[index],
					Constants.b_train[index], Constants.b_trainVar[index],Constants.rho_train[index]};
//					Constants.r_train2[index], Constants.r_train2Var[index], Constants.b_train2, Constants.b_train2Var,
	//				new double[][] {Constants.rho_train2}};
			
			}
			if(expPriors){
				 pow(priors,10);
			}
		}
		return priors[index1 >= priors.length ? 0 : index1];
	}
	
	private static void pow(double[][][] priors2, int pow) {
		for(int i=0; i<priors2.length; i++){
			for(int j=0; j<priors2[i].length; j++ ){
				for(int k=0; k<priors2[i][j].length; k++){
					priors2[i][j][k] = Math.pow(pow, priors2[i][j][k]);
				}
			}
		}
		
	}
	public static boolean expPriors = false; //whether to raise the r_train etc to power 10
	
	

	public static boolean maximiseGlobalRegress() {
		// TODO Auto-generated method stub
		return maximiseGlobalRegress;
	}
	public static boolean maximiseGlobalRegress = false;

	public static String print(double[][] countE) {
		StringBuffer sb = new StringBuffer(print(countE[0]));
		for(int i=1; i<countE.length; i++){
			sb.append("\n"+print(countE[i]));
		}
		return sb.toString();
	}
	public static boolean transferGlobalLocal =true;
	public static boolean transferGlobalLocal() {
		// TODO Auto-generated method stub
		return transferGlobalLocal;
	}
public static double[] gammaRate;//, gammaRate1;
public static String dcFile = null;
public  static boolean readGlobalClusterFile = true;


public  static boolean useSameModelForAllCN = true;;

	public static double gammaRate(int i) {
		// TODO Auto-generated method stub
		if(gammaRate==null ) return 1e-5;
		return Constants.gammaRate[i];
	}

	public static String dcFile() {
		// TODO Auto-generated method stub
		return Constants.dcFile ;
	}

	public static boolean readGlobalClusterFile() {
		// TODO Auto-generated method stub
		return readGlobalClusterFile ;
	}

	public static double[] transform(Double[] dist) {
		double[] res = new double[dist.length];
		for(int i=0; i<res.length; i++){
			res[i] = dist[i];
		}
		return res;
	}

	static class IntDouble implements Comparable{
		int i;
		double d;
		public IntDouble(int k, double e) {
			i = k; d = e;
		}
		@Override
		public int compareTo(Object o) {
			IntDouble o1 = (IntDouble) o;
			if(d > o1.d) return -1;
			else if (d < 01.d) return 1;
			else return 0;
		}
	}
	
	public static int[] getOrder(double[] probs) {
		IntDouble[] p = new IntDouble[probs.length];
		for(int k=0; k<p.length; k++){
			p[k] = new IntDouble(k,probs[k]);
		}
		Arrays.sort(p);
		
		
		int[] res = new int[probs.length];
		for(int k=0; k<p.length; k++){
			res[k] = p[k].i;
		}
		return res;
		
	}
	
	

	public static boolean useSameModelForAllCN() {
		// TODO Auto-generated method stub
		return useSameModelForAllCN ;
	}

	/*public static int whichInd(Integer[][] nocop, Integer[] cn) {
		for(int k=0; k<nocop.length; k++){
			if(nocop[k].equals(cn)){
				return k;
			}
		}
		return -1;
		
		
	}*/
/*public static boolean orthogMix = fa
	public static boolean orthogMix() {
		// TODO Auto-generated method stub
		return orthogMix;
	}
*/
public static boolean allowComponent = false;

	public static boolean allowComponent() {
		// TODO Auto-generated method stub
		return allowComponent;
	}
public static boolean orthogMix = false;
	public static boolean orthogMix() {
		// TODO Auto-generated method stub
		return orthogMix;
	}
public static boolean transferEquilToStart = true;
public static boolean fbHapPanel = false;
	public static boolean transferEquilToStart(int numsite) {
		
		return numsite > 200 && transferEquilToStart;
	}

	public static boolean fbHapPanel() {
		// TODO Auto-generated method stub
		return fbHapPanel ;
	}
public static boolean collapseScatterInd = false;
	public static boolean collapseScatterInd() {
		return collapseScatterInd;
	}
	public static boolean diffRatesPerState = false;
	public  static String[] rsToAnnotate = null;
	public static boolean diffRatesPerState() {
		// TODO Auto-generated method stub
		return diffRatesPerState;
	}

	public static String[] rsToAnnotate() {
		// TODO Auto-generated method stub
		return rsToAnnotate ;
	}
	
public static boolean run2 = false;
	public static boolean run2() {
		if(!run2) return false;
		else{
			if(Constants.fillLikelihood[0] <1){
				for(int k=1; k<Constants.fillLikelihood.length; k++){
					if(Constants.fillLikelihood[k]<1) return false;
				}
				return true;
			}
			else return false;
		}
	}

	public static void normaliseL(double[] probs) {
		double sum=0;
		for(int i=0; i<probs.length; i++){
			probs[i] = Math.exp(probs[i]);
			sum+=probs[i];
		}
		for(int i=0; i<probs.length; i++){
			probs[i] =probs[i]/sum;
		}
	}

	public static void normalise1(double[] probs) {
		double sum = 0;
		int zeroCount=0;
		for(int i=0; i<probs.length; i++){
			if(probs[i]>0){
				sum+=probs[i];
				
			}
			else{
				zeroCount++;
			}
		}
		if(sum<1){
			double diff = (1.0-sum)/(double) zeroCount;
			for(int i=0; i<probs.length; i++){
				if(probs[i]==0){
					probs[i]=  diff;
				}
			}
		}
		
		
		
	}
	public static double dropFracSites(int i) {
		// TODO Auto-generated method stub
		if(dropFracSites==null || i>=dropFracSites.length) {
			return 0.0;
		}
		return dropFracSites[i];
	}
	public static boolean calcHWE = false;
	public static boolean hwe() {
		return calcHWE;
	}
	public static boolean bufferCompress = false;
	public static boolean bufferCompress() {
		// TODO Auto-generated method stub
		return bufferCompress;
	}
	public static boolean writePhased() {
		// TODO Auto-generated method stub
		return writePhased;
	}
	public static boolean writePhased = true;
	public static boolean printAb = true;
	public static boolean printAb() {
		// TODO Auto-generated method stub
		return printAb;
	}
	public static boolean impute() {
		// TODO Auto-generated method stub
		return false;
	}
	public static int hmmFontSize = 8;
	
	public static int hmmFontSize(int i){
		if(i<2){
			return hmmFontSize;
		}
		return 6;
	}
	public static int hmmFontSize() {
		return hmmFontSize;
	}
	
	public static double hmmBubblePow = 0.25;


	public static int[] alleles ;
	public static double hmmBubblePow() {
		// TODO Auto-generated method stub
		return hmmBubblePow;
	}
	/*static String[] dataSpecificFields = new String[] {
		"exponentB", "includeSNPS","pseudoMod1","regionsToInclude", "regionsToExclude"
	};
	public static Map<String, List> getSpreadsheetVals() {
		Map<String, List> m = new HashMap<String, List>();
		try{
		AbstractDataParams sheet = getDataParams(new File(Constants.data_params), Constants.sheet);
		
		String[] cell = sheet.getColumn(0);
		for(int i=0; i<cell.length; i++){
			String nme = cell[i];
			if(nme.length()>0){
				m.put(nme, new ArrayList());
			}
		}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		for(int i=0; i<dataSpecificFields.length; i++){
			m.put(dataSpecificFields[i], new ArrayList());
		}
		
		return m;
	}*/
	private static AbstractDataParams getDataParams(File ff, String[] sheet)  throws Exception{
		return ff.getName().endsWith(".xls") ? new DataParams(ff, sheet) : new CSVDataParams(ff);
	}
	public static int alleles(int i) {
		return alleles[i];
	}
	
	public static boolean globalTrans = true; //whether to use a global transition matrix as base for transition matrices


	public static double oneeighty = 1.0; // value to use for theta max


	public static boolean transformRToLogR = true;
	/*public static boolean transformRToLog(int i){
		if(i<transformRToLogR.length) return transformRToLogR[i];
		else return transformRToLogR[0];
	}*/
	public static boolean globalTrans() {
		// TODO Auto-generated method stub
		return globalTrans;
	}
	public static double exponentB(short data_index) {
		
		return exponentB==null ? 1: exponentB[data_index];
	}
	public static boolean globalRange() {
		// TODO Auto-generated method stub
		return true;
	}
	public static int lengthMod() {
		// TODO Auto-generated method stub
		//assumed gap if none given
		return 1000;
	}
	public static double pi  = 0.01;
	public static double pi() {
		return pi;
	}
	public static boolean sameSizeWithin = true;
	public static boolean sameSizeWithin() {
		// TODO Auto-generated method stub
		return sameSizeWithin;
	}
	public static boolean rotate() {
		// TODO Auto-generated method stub
		return rotate;
	}
	public static boolean rotate = false;
	public static boolean[] includeSNPS;

	public static boolean[] includeSNPS() {
		// TODO Auto-generated method stub
		return includeSNPS;
	}
	public static boolean allowFlips =true;
	public static double alleleDiffThresh=0.5;
	public static boolean allowFlips() {
		// TODO Auto-generated method stub
		return allowFlips;
	}
	public static double alleleDiffThresh() {
		// TODO Auto-generated method stub
		return alleleDiffThresh;
	}
	public static int updateEmissionsIt = -1;
	public static int updateEmissionsIt() {
		return updateEmissionsIt;
	}
	public static String[][] regionsToInclude(int i1) {
		if(regionsToInclude==null || i1>= regionsToInclude.length) return new String[][]{ {"all"}};
		else{
			return parseRegions(regionsToInclude,i1);
		}
		}
	private static String[][] parseRegions(String[][][] regions, int i1){
			File f = null;
			if(regions[i1].length==1 && ((f=new File(regions[i1][0][0])).exists())){
				List<String[]> l = new ArrayList<String[]>();
				try{
				BufferedReader br = new BufferedReader(new FileReader(f));
				String st = "";
				while((st = br.readLine())!=null){
					if(st.startsWith("chr")){
					l.add(st.replaceAll("chr","").split("\\s+"));
					}
					
				}
				}catch(Exception exc){
					exc.printStackTrace();
				}
				regions[i1] = l.toArray(new String[0][]);
			}
			return regions[i1];
		
	}
	public static String[][] regionsToExclude(int i1) {
		if(regionsToExclude==null || i1>= regionsToExclude.length) return new String[][] {{"null"}};
		else return parseRegions(regionsToExclude,i1);
		}
	public static boolean annotateName = false;
	public static boolean annotateName() {
		// TODO Auto-generated method stub
		return annotateName;
	}
	public static float[] randomFloat(int len, float sc) {
		float[] res = new float[len];
		for(int k=0; k<res.length; k++){
			res[k] = Constants.rand.nextFloat()*sc;
		}
		return res;
	}
	public static boolean annotateClusterPosition = true;
	public static boolean annotateClusterPosition() {
		return annotateClusterPosition;
	}
	public static boolean drawAnnotationLines = true;
	public static boolean drawAnnotationLines() {
	return drawAnnotationLines;
	}
	public static String[][] platesToExclude(int ind) {
		// TODO Auto-generated method stub
		return (String[][]) get(platesToExclude,ind);
	}
	public static String[][] phenoToExclude(short index2) {
		// TODO Auto-generated method stub
		return (String[][]) get(phenoToExclude,index2);
	}
	public static String[][] phenoToInclude(short index2) {
		// TODO Auto-generated method stub
		return (String[][]) get(phenoToInclude,index2);
	}
	
	
	private static Object get(Object array, int ind) {
		try{
		if(ind < Array.getLength(array)) return Array.get(array, ind);
		else  return Array.get(array, Array.getLength(array)-1);
		}catch(Exception exc){
			exc.printStackTrace();
			return null;
		}
	}
	public static String[][][] platesToExclude = new String[][][] {{{"null"}}};
	public static String[][][] phenoToExclude = new String[][][] {{{"null"}}};
	public static String[][][] phenoToInclude = new String[][][] {{{"null"}}};
	public static String[][][] platesToInclude = new String[][][] {{{"all"}}};
	public static boolean flip = true;
	public static String[][] platesToInclude(int ind) {
		
		return (String[][]) get(platesToInclude,ind);
	}
	public static int r_panel_width = 800;
//	public static int r_panel_height = 400;
	public static int r_panel_width() {
		return r_panel_width;
	}
	public static int r_panel_height() {
		return Constants.plotHeight;
	}
	
	public static boolean calcAssocAvg = true;
	public static boolean calcAssoc1() {
		// TODO Auto-generated method stub
		return calcAssocAvg;
	}
public static String[] duplicates;
public static String duplSep="###"; //
public static double joinStrokeWidth=0.2;
public static double decayStrokeWidth=1.0;

	
	public static String duplicates(int i){
		if(duplicates==null) return "firstOnly";
		if(i<duplicates.length) return duplicates[i];
		else return duplicates[0];
	}
	
	//public static int extra =0;// (1000*1000)/50;
	/*public static int extra() {
		// TODO Auto-generated method stub
		return extra;
	}*/
	public static String[] thin(short index2) {
		if(index2<thin.length)
		return thin[(int)index2];
		else return thin[0];
	}
	public static String[][] thin = new String[][] {new String[]{"1"}};
public static boolean transformTheta =false;
	public static boolean transformTheta() {
		return transformTheta;
	}
	public static double NAthresh=1.0;
	public static double NAthresh(){
		return NAthresh;
	}
	public static Iterator<Integer> iterator(final int size) {
		return  new Iterator<Integer>(){
           int i=0;
			@Override
			public boolean hasNext() {
				// TODO Auto-generated method stub
				return i<size;
			}

			@Override
			public Integer next() {
				Integer val = i;
				i++;
				return val;
			}

			@Override
			public void remove() {
				// TODO Auto-generated method stub
				
			}
			
		};
	}

	public static boolean showNaN = false;
	public static boolean showNaN() {
		// TODO Auto-generated method stub
		return showNaN;
	}
	public static double rweight = 1.0;
	public static double rweight() {
	return rweight;
	}

	public static Object[] drop(
			Object[] r, List<Integer> toDrop) {
		if(r.length==0 || toDrop.size()==0) return r;
		// TODO Auto-generated method stub
		List l = new ArrayList();
		for(int i=0; i<r.length; i++){
			if(!toDrop.contains(i)){
				l.add(r[i]);
			}
		}
		return l.toArray((Object[]) Array.newInstance(r[0].getClass(), 0));
	}
	public static double printRoundThresh = 0.01;
	public static double printRoundThresh() {
		return printRoundThresh;
	}
	
	public static boolean allowLocalDist = true;
	public static boolean allowLocalDist() {
		// TODO Auto-generated method stub
		return allowLocalDist;
	}
	public static int dataPanelToPlot = 0;
	public static int dataPanelToPlot() {
		// TODO Auto-generated method stub
		return dataPanelToPlot;
	}
	public static boolean trainB() {
		// TODO Auto-generated method stub
		return trainB;
	}
	public static boolean trainR() {
		// TODO Auto-generated method stub
		return trainR;
	}
	public static boolean trainR = true;
	public static boolean trainB = true;
public static boolean trainEmissions = true;
	public static boolean trainEmissions() {
		// TODO Auto-generated method stub
		return trainEmissions;
	}
	
	public static int[][]stateCNOffset;
	public static int[] stateCNOffset(int i) {
		int[] res = new int[inputDir.length];
		for(int k=0; k<res.length; k++){
			res[k] = (stateCNOffset == null || k>stateCNOffset.length) ? 0 : stateCNOffset[k][i];

		}
		return res;
	}
	public static int[][][] mixStates = null;
	private static Boolean mixStatesAllNull=null;
	public static int[][][] mixStates(CompoundMarkovModel hmm) {
		if(mixStates!=null && mixStatesAllNull==null){
			;
			List<String>nme = new ArrayList<String>();
			for(Iterator<State> it = hmm.states(); it.hasNext();){
				nme.add(it.next().getName());
			}
			List<String> nme1 = new ArrayList<String>();
			for(Iterator<State> it = hmm.getMarkovModel(0).states(); it.hasNext();){
				nme1.add(it.next().getName());
			}
			boolean allnull = true;
			for(int k=0; k<mixStates.length; k++){
				if(!(mixStates[k]==null || mixStates[k].length==0 || mixStates[k].length==1 && mixStates[k][0]==null)){
					allnull = false;
					List<int[]> indices = new ArrayList<int[]>();
					for(int j=0; j<mixStates[k].length; j++){
						String from = nme1.get(mixStates[k][j][0]+1);
						String to = nme1.get(mixStates[k][j][1]+1);
						for(int kk=0; kk<nme.size(); kk++){
							int from1 = kk;
							if(nme.get(kk).indexOf(from)>=0){
								
								int to1 = nme.indexOf(nme.get(kk).replaceAll(from, to));
								indices.add(new int[] {from1,to1});
							}
						}
					}
					mixStates[k] = new int[indices.size()][2];
					for(int j=0; j<mixStates[k].length; j++){
						mixStates[k][j] = indices.get(j);
					}
				}
				if(mixStates[k].length==1 && mixStates[k][0]==null) mixStates[k]=null;
			}
			mixStatesAllNull  = allnull;
			if(allnull) {
				mixStates = null;
			}
		}
		// TODO Auto-generated method stub
		return mixStates;
	}
	/*public static double scoreThresh = 1e-5;
	public static double scoreThresh() {
		// TODO Auto-generated method stub
		return scoreThresh;
	}*/
	public static boolean convertAvgToZip = true;
	public static boolean convertAvgToZip() {
		// TODO Auto-generated method stub
		return convertAvgToZip;
	}
	public static double adjustR(double d, int i) {
	//	if(d< minR[i]) minR[i] = d;
//		else if ( d>maxR[i]) maxR[i] = d;
		return Math.min(Math.max(Constants.minR[i], d),Constants.maxR[i]);

	}
	public static String specialCode(short index2) {
		if(Constants.inputDir[index2].equals("317k")) return "#";
		return null;
	}
	
	public static int depthThresh = 1;
	public static int depthThresh() {
		// TODO Auto-generated method stub
		return depthThresh;
	}
	public static double backgroundPower = 0.0; //this is the power which quality score is raised to (and colours out points
	//if set to 0.0 should have no effect.s
	public static double backgroundPower() {
		return backgroundPower;
	}
	public static Boolean[] probeOnly = new Boolean[] {true};
	public static Boolean probeOnly(int index2) {
		// TODO Auto-generated method stub
		return probeOnly[index2>= probeOnly.length ? 0 : index2];
	}
	public static double transformVariance(double d) {
		return Math.pow(d,2);
	}
	/*public static boolean[] drawMeanForm = new boolean[] {true,true};
	public static boolean drawMeanForm(int i) {
		// TODO Auto-generated method stub
		return drawMeanForm[i];
	}*/
	public static String[][] basisNme=new String[0][0];
	public static String[] basisNme(int index2) {
		// TODO Auto-generated method stub
		return basisNme[index2>=basisNme.length ? 0 : index2];
	}
	public static double countThresh2() {
		return countThresh2;
	}
	public static boolean trainIndivClusters = true;
	
	public static boolean trainIndiv() {
		return trainIndivClusters;
	}
	
	public static int[][][] cnToSplit = new int[][][]{ {{0,1,2,3,4,5,6,7,8}}};
	
	
	
	public static int[][] cnToSplit(int index) {
		// TODO Auto-generated method stub
		
		return cnToSplit[index];
	}
	public static Boolean[] beta;
	public static Boolean beta(int index2) {
		return beta[index2];
	}
	public static Double strokeWidth = 10.0;
	public static float strokeWidth() {
		return strokeWidth.floatValue();
	}
	public static Double bafPseudoCount = 500.0;
	public static Double bafPseudo() {
		// TODO Auto-generated method stub
		return bafPseudoCount;
	}
	public static Double countThresh3 = 0.0;
	public static double countThresh3() {
		return countThresh3;
	}
	
	public static String[] useDataAsModel = null; 
	
	public static String[] useDataAsModel() {
		
		if(Constants.useDataAsModel!=null && Constants.modify0[0].length<=Constants.useDataAsModel.length){
			for(int k=0; k<Constants.modify0.length; k++){
				char[] mod0 = new char[Constants.useDataAsModel.length];
				Arrays.fill(mod0, Constants.modify0[k][0]);
				Constants.modify0[k] = mod0;
				Constants.modifyFrac0 = new double[mod0.length];
				Constants.modifyFrac1 = new double[mod0.length];
				Constants.modifyFrac2 = new double[mod0.length];
			//	Constants.modifyFrac3 = new double[mod0.length];
				Constants.modifyFracStart = new double[mod0.length];
				Arrays.fill(modifyFrac0, 1.0/(double)mod0.length);
				Arrays.fill(modifyFrac1, 1.0/(double)mod0.length);
				Arrays.fill(modifyFrac2, 1.0/(double)mod0.length);
				//Arrays.fill(modifyFrac3, 1.0/(double)mod0.length);
				Arrays.fill(modifyFracStart, 1.0/(double)mod0.length);
			}
		}
		// TODO Auto-generated method stub
		return useDataAsModel;
	}
//<<<<<<< .mine
	public static boolean saveMostLikely = false;
	public static boolean mostLikely() {
		// TODO Auto-generated method stub
		return saveMostLikely;
	}
//=======
	
	public static boolean indelFilter() {
		// TODO Auto-generated method stub
		return false;
	}
	//public static boolean removeAmbiguousStrand = false;
	/*public static boolean removeAmbiguousStrand() {
		// TODO Auto-generated method stub
		return removeAmbiguousStrand;
	}*/
	public static List<Integer> indexOf(List strand2, Object object) {
		List<Integer> res = new ArrayList<Integer>();
		for(int i=0; i<strand2.size(); i++){
			if(object==null){
			   if(strand2.get(i)==null) res.add(i);	
			}
			else{
				if(strand2.get(i)!=null && strand2.get(i).equals(object)) res.add(i);
			}
		}
		return res;
	}
	
	public static List<Integer> indexOf(List strand2, double object, double thresh) {
		List<Integer> res = new ArrayList<Integer>();
		for(int i=0; i<strand2.size(); i++){
			
				if(Math.abs(((Number)strand2.get(i)).doubleValue()-object)<thresh) res.add(i);
			
		}
		return res;
	}
	public static boolean swapAllelesDuringFlip = false; // instead of using complementarity, should it just flip alleles
	//useful for ancestral allele vs ref allele correction
	public static boolean switchAllelesDuringFlip() {
		// TODO Auto-generated method stub
		return swapAllelesDuringFlip;
	}
	public static int equaliseGroupMode = 0; //0 means remove things which are in one but not other;  1 means set the missing one to have zero 
	//minor allele frequency (i.e. all A or AA)
	public static int equaliseGroupMode() {
		// TODO Auto-generated method stub
		return equaliseGroupMode;
	}
//<<<<<<< .mine
	public static boolean standardiseLRR = false;
	public static boolean standardiseLRR() {
		return standardiseLRR;
	}
	public static double calcMean(double[] lrrVals) {
		double m = 0;
	 double cnt=0;
		for(int i=0; i<lrrVals.length; i++){
			if(!Double.isNaN(lrrVals[i])) {
				m+=lrrVals[i];
				cnt++;
			}
		}
		return m/cnt;
	}
	public static double calcSd(double[] lrrVals, double mu) {
		double m = 0;
	 double cnt=0;
		for(int i=0; i<lrrVals.length; i++){
			if(!Double.isNaN(lrrVals[i])) {
				m+=Math.pow(lrrVals[i]-mu,2);
				cnt++;
			}
		}
		return Math.sqrt(m/cnt);
	}
//=======
//>>>>>>> .r224
	public static double longThresh = 0.5;
	public static double longThresh() {
		return longThresh;
	}
	public static int getKaryoFile(File parentFile) throws Exception {
		File karyoDir = parentFile;
		final String buildf1 = Constants.build(0).split("\\.")[0];
	    File karyo = null;
		for(int i=0; i<3 && (karyo==null || karyo.length()==0) && !karyoDir.getAbsolutePath().equals("/"); i++){
//		File karyoDir = f.getParentFile().getParentFile();
		File[] karyo1 = karyoDir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.getName().indexOf("karyotype")>=0 && pathname.getName().indexOf(buildf1)>=0;
			}

			
		});
		if(karyo1!=null){
		 karyo = karyo1.length ==0 ? null : karyo1[0];
		}
		 karyoDir = karyoDir.getParentFile();
		}
		
		if(karyo!=null ){
			
			BufferedReader br = new BufferedReader(new FileReader(karyo));
			String st = "";
			while((st = br.readLine())!=null){
				String[] str = st.split("\t");
				if(str[0].equals(Constants.chrom0())){
					return Integer.parseInt(str[1]);
					
					}
				}
		}
		return -1;
	}
	public static boolean[] dropMonomorphic = new boolean[]{false};
//>>>>>>> .r241
	public static boolean dropMonomorphic(short index2) {
		return dropMonomorphic[(int) Math.min(index2, dropMonomorphic.length-1)];
	}
	
	public static double[] coverage = new double[]{Double.NaN};
	public static double coverage(int dataIndex) {
		// TODO Auto-generated method stub
		return dataIndex < coverage.length ? coverage[dataIndex] : coverage[0];
	}
	public static int sum(boolean[] cnts) {
		int cnt=0;
		for(int i =0; i<cnts.length; i++){
			if(cnts[i])cnt++;
		}
		return cnt;
	}
	
	public static boolean useOriginalLikelihoodAsUncertainty = false;
	public static boolean useOriginalLikelihoodAsUncertainty() {
		return useOriginalLikelihoodAsUncertainty;
	}
	public static boolean removeIndels = false;

	public static double excludeBafThresh=0.1;
	public static boolean removeIndels() {
		// TODO Auto-generated method stub
		return removeIndels;
	}
	public static double excludeBafThresh() {
		// TODO Auto-generated method stub
		return excludeBafThresh;
	}
	public static String[] color() {
	
		return color;
	}
	public static int numPcs = 0;
	
	public static String[] color="pink:red:gray~0.5:green:cyan:blue:yellow:orange:MAGENTA:green~0.5:green~0.3:yellow~0.5:yellow~0.3".split(":");
	public static int numPcs() {
		// TODO Auto-generated method stub
		
		return numPcs;
	}
	public static double fgPower = 1.0 ; //0 no effect, 1 max effect
	public static double foregroundPower() {
		// TODO Auto-generated method stub
		return fgPower;
	}
	public static boolean usePhenotypeVals = false;
	public static boolean usePhenoTypeVals() {
		// TODO Auto-generated method stub
		return usePhenotypeVals;
	}
	public static boolean plotDistributions(){
		return usePhenotypeVals;
	}
	public static boolean showNA = false;
	public static boolean showNA() {
		// TODO Auto-generated method stub
		return showNA;
	}
	public static boolean equalRank = true;
	public static boolean equalRank() {
		// TODO Auto-generated method stub
		return equalRank;
	}
	public static String[][] split(String[] sb, String c) {
		// TODO Auto-generated method stub
		String[][] res = new String[sb.length][];
		for(int i=0; i<res.length; i++){
			res[i] = sb[i].split(c);
		}
		return res;
	}
	public static double min(double[] pseudo) {
		double min =pseudo[0];
		for(int i=1; i<pseudo.length; i++){
			if(pseudo[i] < min) min = pseudo[i];
		}
		return min;
	}
	public static int cumulativeR(int i) {
		// TODO Auto-generated method stub
		return i<cumulativeR.length ? cumulativeR[i] : 1;
	}
public static int[] cumulativeR = new int[] {1};
public static boolean plotMerged=false;
public static boolean plotMerged() {
	// TODO Auto-generated method stub
	return plotMerged;
}
public static double lrr0 = -5.0;
public static double lrr0() {
	// TODO Auto-generated method stub
	return lrr0;
}
public static double seqQualityThresh = 1e-3;  //below this, the GT will be taken to be 100% correct
public static double seqQualityThresh() {
	// TODO Auto-generated method stub
	return seqQualityThresh;
}
public static boolean quantileNormalise = false;
public static boolean quantileNormalise() {
	// TODO Auto-generated method stub
	return quantileNormalise;
}
public static boolean excludeMultiAllelicSites=true;
public static boolean excludeMultiAllelicSites() {
	// TODO Auto-generated method stub
	return excludeMultiAllelicSites;
}
public static boolean isLogProbs = false;
public static boolean isLogProbs() {
	// TODO Auto-generated method stub
	return isLogProbs;
}
public static double priorWeight = 1e5;

public static double[] initialCellularity=new double[] {0.99,1.0};
public static double priorWeight() {
	// TODO Auto-generated method stub
	return priorWeight;
}
public static double[] initialCellularity() {
	// TODO Auto-generated method stub
	return Constants.initialCellularity;
}
public static double continueThresh = 0.5; //put to zero to keep singletons, and to not read through singletons
public static double continueThresh() {
	// TODO Auto-generated method stub
	return continueThresh;
}
public static double[] scaleLoc = null;
public static double[] scaleLoc() {
	// TODO Auto-generated method stub
	return scaleLoc;
}
public static int noSampleFromBeta = 1;
public static int noSampleFromBeta() {
	// TODO Auto-generated method stub
	return noSampleFromBeta;
}
public static boolean updateCellularity = true;
public static boolean updateCellularity() {
	
	// TODO Auto-generated method stub
	return updateCellularity ;
}
public static double[] downsample = new double[] {1,1}; //tumour then ref
public static double downsample(int i) {
	return Constants.downsample[i];
}
public static double[] dilute = new double[] {1.0,0.0}; //tumour then ref

public  static double proportionOfDistanceToNextProbeToExtend = 0.0;


public static double[] dilute() {
	// TODO Auto-generated method stub
	return dilute;
}
public static double proportionOfDistanceToNextProbeToExtend() {
	// TODO Auto-generated method stub
	return proportionOfDistanceToNextProbeToExtend;
}
//public static boolean removeSingletonDropouts = false;
//public static boolean removeSingletonDropouts() {
	// TODO Auto-generated method stub
	//return removeSingletonDropouts;
//}
public static boolean useAvgDepth = false;
public static boolean useAvgDepth() {
	// TODO Auto-generated method stub
	return useAvgDepth;
}
public static boolean makeNewRef = false;
public static boolean makeNewRef() {
	// TODO Auto-generated method stub
	return makeNewRef;
}

public static int mergeSampleLen = 0;
public static int mergeSamplesLen() {
	// TODO Auto-generated method stub
	return mergeSampleLen;
}
public static String reference = "REFALL";
public static String reference() {
	// TODO Auto-generated method stub
	return reference;
}
public static Boolean ratioAsLevels = null;
public static Boolean ratioAsLevels() {
	// TODO Auto-generated method stub
	return ratioAsLevels;
}
public static int maxPloidy1() {
	// TODO Auto-generated method stub
	return ratioAsLevels==null ? 1 : maxPloidy1;
}
public static int trainCellularity = 2;
public static int trainCellularity() {
	// TODO Auto-generated method stub
	return trainCellularity;
}
public static int betaDownWeight=1;
public static int betaDownWeight() {
	return betaDownWeight;
}


	
	
	

}
