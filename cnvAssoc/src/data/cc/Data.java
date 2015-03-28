package data.cc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import Jama.Matrix;

import data.util.Utils;




public class Data extends AbstractData{
static double var_thresh = 1e-90;

//static boolean adjustVarianceThresh = true;
	static final String build = "build36";
	public static double[] perc = new double[] {0,1};//0.05;  //variance bound
	//int[] inds;
	String type;
	
public void reorder(List<Integer> s) {
		RealMatrix matr = this.data;
		RealMatrix matr_ = matr.copy();
		for(int k=0; k<s.size(); k++){
			int ind = s.get(k);
			matr.setRow(k, matr_.getRow(ind));
		}
	}
	
	
	
	public void logData(){
		for(int k=0; k<this.probes.size(); k++){
			if(probes.get(k).equals("ones")) continue;
			else{
				for(int j=0; j<this.data.getRowDimension(); j++){
					this.data.setEntry(j, k, Math.log(data.getEntry(j, k)));
				}
			}
		}
	}
	
	public void makebinary(double log) {
		for(int k=0; k<this.probes.size(); k++){
			if(probes.get(k).equals("ones")) continue;
			else{
				for(int j=0; j<this.data.getRowDimension(); j++){
					this.data.setEntry(j, k, data.getEntry(j, k)<log ? 0 : 1);
				}
			}
		}
		
	}
//	List<String> bins;
	
	/*public static void main(String[] args){
		try{
		Data d = new Data(new File(args[0]),new String[]{"Original"});
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}*/
	public void getPCs(){
		  Jama.SingularValueDecomposition svd =     new Jama.SingularValueDecomposition(
				  new Jama.Matrix(data.transpose().getData()));//alg.transpose(X));
			 double[] d =  svd.getSingularValues();
			 Matrix v = svd.getV();
			Matrix ev = (new  Jama.EigenvalueDecomposition(v)).getV();
//				  new SingularValueDecomposition(X.transpose());//alg.transpose(X));
			  System.err.println("h "+svd.rank());
	}
	

	public int[] nonNaCols(boolean all,boolean remZero){
		boolean[] include = new boolean[data.getColumnDimension()];
		int cntInc = 0;
	//	boolean meanCentre = false;
	//	boolean setNaAsAvg = false;
		for(int i=0; i<this.data.getColumnDimension(); i++){
			if(this.probes.get(i).equals("ones")){
				include[i] = true;
				cntInc++;
				continue;
			}
			double[] d = data.getColumn(i);
			
			// double[] weights = new double[genotypes[0].length];
			//
		//	for (int i = 0; i < include.length; i++) {
				double sum = 0;
				double cnt = 0;
				boolean hasNA = false;
				boolean allNA = true;
				for (int j = 0; j < d.length; j++) {
					if (!Double.isNaN(d[j])) {
						sum += d[j];
						cnt++;
						allNA = false;
					} else {
						hasNA = true;
					}

				}
			double mean = sum / cnt;
			cnt = 0;
				sum = 0;
				for (int j = 0; j < d.length; j++) {
					if (!Double.isNaN(d[j])) {
					sum += Math.pow(d[j] - mean, 2);
						cnt++;
					//	if(meanCentre){
					//		data.setEntry( j, i, d[j]-mean);
					//	}
				//	} else if (setNaAsAvg) {
						
					//	data.setEntry( j, i,meanCentre ? 0 : mean);
					}

				}
				double var = sum / cnt;
				include[i] = (var<=var_thresh && remZero) ||  (all && allNA) || (!all && hasNA) ? false
						: true;
				if (include[i])
					cntInc++;
				else{
					System.err.println(this.type+" removed column "+" "+allNA+" "+hasNA+" "+var+" "+this.probes.get(i));
				}
		}
			int[] alias = new int[cntInc];
			cntInc = 0;
			for (int k = 0; k < include.length; k++) {
				if (include[k]) {
					alias[cntInc] = k;
					cntInc++;
				}
			}
			return alias;
		}
	
	public int[] nonNaRows(boolean all){
		boolean[] include = new boolean[data.getRowDimension()];
		int cntInc = 0;
		int k_start =0;
		if(probes.get(0).equals("ones")) k_start = 1;
		for(int i=0; i<this.data.getRowDimension(); i++){
			double[] d = data.getRow(i);
			if(d.length<=k_start) throw new RuntimeException("no cols");
			// double[] weights = new double[genotypes[0].length];
			//
		//	for (int i = 0; i < include.length; i++) {
				double sum = 0;
				double cnt = 0;
				boolean hasNA = false;
				boolean allNA = true;
				for (int j = k_start;j < d.length; j++) {
					if (!Double.isNaN(d[j])) {
						sum += d[j];
						cnt++;
						allNA = false;
					} else {
						hasNA = true;
					}

				}
				double mean = sum / cnt;
				cnt = 0;
				sum = 0;
				for (int j = k_start; j < d.length; j++) {
					if (!Double.isNaN(d[j])) {
						sum += Math.pow(d[j] - mean, 2);
						cnt++;
						
					} 
					

				}
			//	double var = sum / cnt;
				include[i] =  (all && allNA || !all && hasNA) ? false
						: true;
				if (include[i])
					cntInc++;
				else{
					System.err.println("removed indiv "+this.name+" "+" "+allNA+" "+hasNA+" "+this.indiv.get(i));
				}
		}
			int[] alias = new int[cntInc];
			cntInc = 0;
			for (int k = 0; k < include.length; k++) {
				if (include[k]) {
					alias[cntInc] = k;
					cntInc++;
				}
			}
			return alias;
		}
	
	/** new probes */
	public void extend(Data d){
		int len = probes.size();
		int len1 = d.probes.size();
		boolean ones = d.probes.get(0).startsWith("ones");
		int start = ones ? 1 : 0;
		this.probes.addAll(d.probes.subList(start, len1));
		if(d.loc.size()>0) this.loc.addAll(d.loc.subList(start, d.loc.size()));
	
	
		RealMatrix data1 = new Array2DRowRealMatrix(indiv.size(), this.probes.size());
		
		for(int j=0; j<indiv.size(); j++){ 
			for(int k=0; k<len; k++){
				data1.setEntry(j, k, data.getEntry(j, k));
			}
			int j1 = d.indiv.indexOf(this.indiv.get(j));
			for(int k=start; k<len1; k++){
				data1.setEntry(j, k+len - start, j1<0 ? Double.NaN : d.data.getEntry(j1, k));
			}
		}
		this.data = data1;
		if(probes.size()!=data.getColumnDimension()) {
			int d1 = data.getColumnDimension();
			int d2 = d.probes.size();
			throw new RuntimeException("!!");
		}
		
		
	}
	
	//final public String chrom;
	/** allows overlaps in d */
	public void append(Data d){
		if(d==null) return;
		List<Integer> alias = new ArrayList<Integer>();
		List<String> indiv1  = new ArrayList<String>();
		for(int i=0; i<d.indiv.size(); i++){
			if(!this.indiv.contains(d.indiv.get(i))){
				alias.add(i);
				indiv1.add(d.indiv.get(i));
			}
		}
		
		int len1 = indiv.size();
		int len2 = indiv1.size();
		this.indiv.addAll(indiv1);
		List<String> probes1 = new ArrayList<String>(probes);
		probes1.retainAll(d.probes);
		int[] inds = this.getAlias(d.probes, probes1);
		int[] inds1 = this.getAlias(probes,probes1);
		
		if(loc.size()>=inds.length) this.loc =
			new ArrayList<Integer>(Arrays.asList((Integer[])sublist(this.loc.toArray(new Integer[0]), inds1)));
		
		
		RealMatrix data1 = new Array2DRowRealMatrix(indiv.size(), probes1.size());
		int len =0;
		for(int k=0; k<probes1.size(); k++){
			len += (inds[k]<0) ? 0 :1;
			CC.checkAndFix(data, d.data, inds1[k], inds[k], probes1.get(k));
			for(int j=0; j<len1; j++){
				data1.setEntry(j, k, data.getEntry(j, inds1[k]));
			}
			for(int j=0; j<len2; j++){
				data1.setEntry(j+len1, k,d.data.getEntry(alias.get(j), inds[k]));
			}
		}
		//int[] alias_c = new int[len];
		this.data = data1;
		this.probes = probes1;
	//	int[] alias_r = nonNaRows(false);
	}
	
	public String toString(){
		return this.type_name;
	}
	
private static int[] getAlias(List<String> subset, List<String> superset ) {
	int[] inds = new int[superset.size()];
		for(int k=0; k<superset.size(); k++){
			int index = subset.indexOf(superset.get(k));
			inds[k] = index;
			if(index<0){
				throw new RuntimeException("!!");
			}
		}
		
		return inds;
	}
	
public Data(Data d){
	this.probes = new ArrayList<String>(d.probes);
	//this.chrom = d.chrom;
	this.indiv = new ArrayList<String>(d.indiv);
	this.data = d.data.copy();
	this.name = d.name;
	this.type_name = d.type_name;
	this.type = d.type;
	this.doubleBin = d.doubleBin;
}

public Data(File file, String[] type, String chr, int[] start, boolean b,
		boolean c, String[] strings2, boolean d, String string, boolean e,
		boolean f) throws Exception {
	this(file.getParentFile(),type, Arrays.asList(new String[] {chr}), Arrays.asList(new Integer[] {start[0]})
			, Arrays.asList(new Integer[] {start[1]}),b,c,strings2,d,string,e,f);
}
public Data(File file, String[] type, List<String> chr, int[] start, boolean b,
		boolean c, String[] strings2, boolean d, String string, boolean e,
		boolean f) throws Exception {
	this(file.getParentFile(),type, chr, Arrays.asList(new Integer[] {start[0]})
			, Arrays.asList(new Integer[] {start[1]}),b,c,strings2,d,string,e,f);
}

public Data(File dir, String[] type_, List<String> chr, List<Integer> start, List<Integer> end, boolean modify, boolean addConstant,
			String[] rest_, boolean only, String type_name, boolean removeZeroVar, boolean correct) throws Exception{
		//System.err.println("opening "+f.getAbsolutePath());
		this.type_name = type_name;
		List<String> rest = rest_==null ? null : Arrays.asList(rest_);
		this.name = dir.getName()+(type==null ?"":Arrays.asList(type));
	//	ZipFile zf = new ZipFile(f);
		//this.chrom = chr;
		List<Integer> chromPos = new ArrayList<Integer>();
		List<String> chroms = new ArrayList<String>();
		
		probes = this.getProbes(dir, chr, start, end,rest, only,this.loc,  chroms,chromPos);
		if(probes==null || probes.size()==0 || probes.size()==1 && probes.get(0).equals("ones")){
			probes = this.getProbes1(dir, chr, start, end, rest, only, this.loc, chroms,chromPos);
		}
		System.err.println("HERE "+dir+" "+chr+" "+start+" "+end+" "+rest+"\n"+probes);
		//File[] f = this.getFiles(dir, chr);
		File tof = new File(dir,chroms.get(0)+".zip");
		System.err.println("opening "+tof.getAbsolutePath());
		ZipFile zf = new ZipFile(tof);
		this.indiv = read(zf,"Samples",0);
		List<String[]> names = readSplit(zf,"Name");
		int inds =0;
		int length =0;
		if(names!=null){
			inds = -1;
			List<String> l =Arrays.asList( names.get(0));
			length = l.size();
			for(int k=0; k<l.size(); k++){
				l.set(k, l.get(k).trim());
			}
		//	for(int k=0; k<inds.length; k++){
			for(int kk=0; kk<type_.length && inds <0; kk++){
				String type1 = type_[kk].split("__")[0];
				inds = l.indexOf(type1);
				
				if(inds>=0) {
					
					this.type = type_[kk];
					System.err.println(type);
					break;
				}
			}
			if(inds<0) {
				throw new RuntimeException("!! "+Arrays.asList(type_)+" "+Arrays.asList(l));
			}
		}
		/*if(names==null){
			inds = new int[] {0};
		}
		else{
			inds = new int[type.length];
			
		}
		if(inds.length>1) throw new RuntimeException("!!");
		*/
	
		/*if(type_name.equals("baf")){
			for(int k=probes.size()-1; k>=0;k--){
				String pr = probes.get(k);
				if(pr.toLowerCase().startsWith("cnv")){
					probes.remove(k);
					if(k<loc.size()) loc.remove(k);
				}
			}
		}*/
		File correction = new File(dir,"correction.txt");
		data = new Array2DRowRealMatrix(indiv.size(), probes.size());
		//Correction corr =  MergedData.correct && correct  && correction.exists() && correction.length()>0 ?  new Correction(correction, this.indiv,length, Integer.MAX_VALUE) : null;
//		ZipFile zf = null;
		int currPos = 0;
		int nxtPos = currPos+1 < chromPos.size() ? chromPos.get(currPos+1 ) : -1;
		for(int i=0; i<probes.size(); i++){
			if(nxtPos==i){
				currPos++;
				nxtPos = currPos+1 < chromPos.size() ? chromPos.get(currPos+1 ) : -1;
				zf.close();
				zf = new ZipFile(new File(dir, chroms.get(currPos)+".zip" ));
			}
			
			if(probes.get(i).equals("ones")){
				for(int j=0; j<indiv.size(); j++){
					data.setEntry(j, i, 1.0);
				}
			}
			else{
				
			String[] probe = probes.get(i).split("__");
			List<String[]> l  = readSplit(zf,probe[0]);
			/*if(corr!=null){
				corr.correct1(l, inds, false);
			}*/
			for(int j=0; j<l.size(); j++){
				String[] str = l.get(j);
			//	for(int k=0; k<inds.length; k++){
					int i1 =i;
					String val = str.length==1 && str[0].equals("null") ? "NaN": str[inds];
					double v;
					if(val.equals("NC") || val.equals("NA")){
						v = Double.NaN;
					}
					else if(type!=null && this.type.endsWith("ype__A")){
						 v = count(val.toCharArray(),'A');
					}
					else if(type!=null && this.type.endsWith("ype__B")){
						 v = count(val.toCharArray(),'B');
					}
					else if(probe.length>1){
						
						v = getBin(probe[1], val);
						if(v>2){
							System.err.println("h");
						}
						//System.err.println(probe[1]+" "+val+" "+v)
					}
					else{
						if(val.equals("NC") || val.equals("NA")){
							v = Double.NaN;
						}
						
						else{
							v = 	Double.parseDouble(val);
						}
					}
				
					data.setEntry(j, i1,v );
				//}
			}
			}
		}
/*		if(inds.length>1){
		List<String> pr = new ArrayList<String>();
			for(int k=0; k<inds.length;k++){
				for(int j=0; j<this.probes.size(); j++){
					pr.add(probes.get(j)+"##"+type[k]);
				}
			}
			this.probes = pr;
		}*/
		if(modify){
			int[] rows = this.nonNaRows(true);

			int[] alias = this.nonNaCols(true, removeZeroVar);
//		for(int k=0; k<rows.length; k++){
//			rows[k] = k;
//		}
	if(rows.length==0) throw new RuntimeException("removed all rows "+this.type);
    if(alias.length==0) throw new RuntimeException("removed all cols "+this.type);	
			this.data = data.getSubMatrix(rows, alias);
		this.probes =
			new ArrayList<String>(Arrays.asList((String[])sublist(this.probes.toArray(new String[0]), alias)));
		if(loc.size()>=alias.length) this.loc =
			new ArrayList<Integer>(Arrays.asList((Integer[])sublist(this.loc.toArray(new Integer[0]), alias)));
		this.indiv =
			new ArrayList<String>(Arrays.asList((String[])sublist(this.indiv.toArray(new String[0]), rows)));
		
		}
		//this.getPCs();
		
	}
	
	

	

	private Object sublist(Object genotypes2, int[] alias) {
		Object obj = Array.newInstance(Array.get(genotypes2, 0).getClass(),
				alias.length);
		for (int k = 0; k < alias.length; k++) {
			Array.set(obj, k, Array.get(genotypes2, alias[k]));
		}
		return obj;
	}
	
	public static int period =1;
	
	public static List<String> getProbes(File zf1, String string, int i, int j,
			Object loc1, boolean only, List<Integer> loc) throws Exception{
		// TODO Auto-generated method stub
	return Data.getProbes(zf1, Arrays.asList(string.split(":")), Arrays.asList(new Integer[] {i}),
			Arrays.asList(new Integer[] {j}),null,only,
			loc, new ArrayList<String>(), new ArrayList<Integer>());
	}
	public static List<String> getProbes(File dir,
			List<String> chr, List<Integer> start, List<Integer> end, List<String> restrict, boolean only, 
			List<Integer> loc, List<String> chrom_new, List<Integer> chromStart)  throws Exception{
		System.err.println("getting probes");
		System.err.println(dir.getAbsolutePath());
		System.err.println(chr);
		File[] f = getFiles(dir, chr);
		if(dir.getName().indexOf("vntr")>=0 || dir.getName().indexOf("watson")>=0 ) return null;
		if(f.length==0)
			return null;
		else if(f.length==1){
			ZipFile zf1 = new ZipFile(f[0]);
			boolean nu = zf1.getEntry("SNPS")==null;
			zf1.close();
			//if(nu){
			//	return null;
			//}
			
		}
		File in = new File(dir, build+".gc.var.txt");
		
		BufferedReader br =Utils.getBufferedReader(in);
		if(br==null){
			in = new File(dir, build+".var.txt");
			br = Utils.getBufferedReader(in);
		}
		ZipFile zf = null;
		if(br==null && f.length==1){
		zf = new ZipFile(f[0]);
			br = Utils.getBufferedReader(zf, "SNPS");
		}
		if(br==null) return null;
		double var_low = 0;
		double var_high = Double.POSITIVE_INFINITY;
		int var_index = -1;
		boolean hasVar = false;
		if(in.exists()){
			File histFile = new File(dir, in.getName().replace(".var.txt", ".hist.txt"));
			
			double percLow = perc[0];
			double percHigh = perc[1];
			if(histFile.exists() && histFile.length()>0 && perc[1]<1 && perc[0]>0){
				BufferedReader br1 = Utils.getBufferedReader(histFile);
				String st = "";
				hasVar = true;
				while((st = br1.readLine())!=null){
					String[] str = st.split("\\s+");
					double p = Double.parseDouble(str[0]);
					if( p <=percHigh){
						var_high = Double.parseDouble(str[1]);
					}
					if(p<=percLow){
						 var_low = Double.parseDouble(str[1]);
						//if(var_low)
					}
				}
			}
		}
		
		
		//ZipEntry snps =zf.getEntry("SNPS");
		int cnt=0;
		List<String> probes = new ArrayList<String>();
		if(CC.addOne){
			loc.add(0);
			probes.add("ones");
		}
		
		//	List<String > l1 =  Arrays.asList(Compressor.read(zf, "Name").get(1).split("\\s+"));
			int snp_index =3;//l1.indexOf("id");
		
			//BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(snps)));
			String st = "";
			String chrom = "";
			int chrom_ind = -1;
			
			//String chr1 = "chr"+chr;
			for(int kk=0;(st = br.readLine())!=null;kk++){
				String[] str = st.split("\\s+");
				if(kk==0 && hasVar){
					var_index = in.getName().indexOf("var")>=0 ? str.length-1 : -1;
				}
				String chr1 = str[0].substring(3);
				if(!chr1.equals(chrom)){
					chrom = chr1;
					chrom_ind = chr.indexOf(chrom);
					if(chrom_ind>=0){
						chrom_new.add(chrom);
						chromStart.add(probes.size());
					}
				}
				if(chrom_ind<0) continue;
				if(var_index>=0){
					double v = Double.parseDouble(str[var_index]);
					if(v<0) throw new RuntimeException("variance cannot be negative");
					if(v>var_high || v < var_low) {
						continue;
					}
				}
//				chromStart.indexOf();
				{
					int pos = Integer.parseInt(str[1]);
					int chrom_ind1 = Math.min(chrom_ind, start.size()-1);
					if(pos>=start.get(chrom_ind1) && pos<=end.get(chrom_ind1)){

						if(restrict==null || restrict.size()==0 || (restrict.contains(str[snp_index]) && only ) 
								||(!restrict.contains(str[snp_index]) && !only) ){
						//	if(cnt == period) cnt=0;
							//if(cnt==0){
								probes.add(str[snp_index]);
								loc.add(Integer.parseInt(str[1]));
						//	}
						//	cnt++;
						}

					}
				}
			}
			br.close();
		if(zf!=null) zf.close();
		return probes;
	}
	
	private static File[] getFiles(File dir, final List<String> chr) {
		File[] f = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File arg0) {
				String name = arg0.getName();
				if(name.endsWith(".zip")){
					if(chr.contains(name.substring(0,name.length()-4))){
						return true;
					}
				}
				return false;
			}
			
		});
		return f;
	}
	public static List<String> getProbes1(File dir, List<String> chr, List<Integer> start, List<Integer> end, List<String> restrict, boolean only, 
			List<Integer> loc, List<String> chroms, List<Integer> chromStart)  throws Exception{

	
		List<String> probes = new ArrayList<String>();
		if(CC.addOne){
			loc.add(0);
			probes.add("ones");
		}
		for(int kk=0; kk<chr.size(); kk++){
			File f1 = new File(dir, chr.get(kk)+".zip");
			chromStart.add(probes.size());
			chroms.add(chr.get(kk));
			if(f1.exists()){
				ZipFile zf = new ZipFile(f1);
				Enumeration it = zf.entries();
			
				while(it.hasMoreElements()){
					String name = ((ZipEntry)it.nextElement()).getName();
					if(restrict!=null && restrict.size()>0 && (restrict.indexOf(name))<0) continue;
					//if(name.startsWith("VNTRA")) continue;
					if(!name.startsWith("bins_") && !name.startsWith("Name") && ! name.equals("SNPS") && ! name.startsWith("Sample")){
						String binname = "bins_"+name;
						File bn1 = new File(binname);
						List<String[]> bins;
						if(  bn1.exists() && bn1.length()>0){
							bins = readSplit(bn1, binname);
						}
						else bins = readSplit(zf, binname);
						if(bins!=null){
							for(int k=0; k<bins.size(); k++){
								probes.add(name+"__"+getName(bins.get(k)));
							}
						}
						else probes.add(name);
					}
				}
			}
		}
		return probes;
	}

	private static String getName(String[] strings1) {
		String[] strings;
		if(strings1.length>1){
			strings = strings1;
		}
		else{
			strings = strings1[0].split("<");
			strings[0] = strings[0]+("<=");
			strings[1] = strings[1].replaceAll("=", "");
			strings[2] = "<"+strings[2];
		}
	return strings[0]+"_"+strings[1]+"_"+
	 (strings.length<=2 ? "<Infinity" : strings[2]);
	}

	 Boolean doubleBin=null;
	private double getBin(String bin, String val) {
		if(val.indexOf("FAIL")>=0 || val.trim().length()==0 || !this.type_name.equals("plate") && val.indexOf(';')<0){
			return Double.NaN;
		}
		
		String[] str1 = val.split(";");
		int cnt=0;
		//Arrays.fill(res, -1);
		outer: for(int k=0; k<str1.length; k++){
			String[] str = bin.split("_");
			if(doubleBin==null){
				try{
					Double v = Double.parseDouble(str1[k]);
				}catch(Exception exc){
							doubleBin= false;
						 
				}
				if(doubleBin==null)doubleBin=true;
			}
			
			if(doubleBin){
		str1[k] = str1[k].replaceAll(",", "");
		Double v = Double.parseDouble(str1[k]);
		//for(int i=0; i<bins.size(); i++){
			
		/*if(MergedData.lowerB!=null){
		   if(str.length <3 || evalLt(str[2], MergedData.lowerB)){
			   if(v>=MergedData.lowerB && evalLt(str[2],v)){
				   cnt++;
			   }
		   }
		   else{
			   if(v<=MergedData.lowerB && evalGt(str[0],v)){
				   cnt++;
			   }
		   }
		}
		
		else{*/
			if(evalGt(str[0],v) &&  evalLt(str[2],v)){
				cnt++;
			//	continue outer;
			}
		//}
		
			}
			if(!doubleBin){
				if(str1[k].equals(str[1])){
					 cnt++;
				 }
			}
			
		//}
		}
		return cnt;
	}

	private boolean evalGt(String string, Double v) {
		int ind = string.indexOf('<');
		double v1 = Double.parseDouble(string.substring(0,ind));
		boolean ge = string.endsWith("<=");
		if(ge && v1<=v) return true;
		else if(!ge && v1<v) return true;
		else return false;
	}
	
	private boolean evalLt(String string, Double v) {
	
		boolean le = string.startsWith("<=");
		if(le){
			String r = string.substring(2);
			if(r.equals("Inf")) {
				return true;
			}
			double v1 = Double.parseDouble(r);
			return (v<=v1);
		}
		else{
			double v1 = Double.parseDouble(string.substring(1));
			return v<v1;
		}
	}

	private int count(char[] charArray, char ch) {
		// TODO Auto-generated method stub
		int cnt=0;
		for(int k=0; k<charArray.length; k++){
			if(charArray[k]==ch) cnt++;
		}
		return cnt;
	}

	public static List<String> read(ZipFile zf,  String probe, int i) throws Exception {
		ZipEntry ent = zf.getEntry(probe);
		if(ent==null) return null;
		List<String> res = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
		String st = "";
		while((st =br.readLine() )!=null){
			
			res.add(st.split("\t")[0].split("#")[0]);
		}
		br.close();
		return res;
		
	}
	public static List<String[]> readSplit(ZipFile zf,  String probe) throws Exception {
		ZipEntry ent = zf.getEntry(probe);
		if(ent==null) return null;
		List<String[]> res = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
		String st = "";
		for(int i=0; (st =br.readLine() )!=null; i++){
			String[] str = st.split("\t");
			 res.add(str);
		}
		br.close();
		return res;
		
	}
	public static List<String[]> readSplit(File zf,  String probe) throws Exception {
		
		List<String[]> res = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new FileReader(zf));
		String st = "";
		for(int i=0; (st =br.readLine() )!=null; i++){
			String[] str = st.split("\t");
			 res.add(str);
		}
		br.close();
		return res;
		
	}
	private static Collection<? extends String[]> getAll(String[] split, int num) {
		double st = Double.parseDouble(split[0].substring(0,split[0].length()-1));
		double end = Double.parseDouble(split[2].substring(2));
		
		List<String[]> res = new ArrayList<String[]>();
		if(num>20){
		for(double k=st; k<end; k++){
			res.add(new String[] {k+"<",split[1]+(""+(int)k),"<="+(k+1)});
		}
		}
		else{
			double incr = (end-st)/num;
			for(double k=0; k<num; k++){
				double k1 = st+incr*k;
				double start = k1;
				double end1 = k1+incr;
				String[] str =new String[] {start+"<",split[1]+(""+k1),"<="+end1}; 
				res.add(str);
			}
		}
		return res;
	}

	public void applyRecoding(int[][] coding1, boolean keepOrig) {
		int len = this.data.getColumnDimension()-1;
		int offset = (keepOrig ? 1 : 0);
		double[][] genotypes = this.data.transpose().getData();
		String[] pi = new String[1+len * (offset + coding1.length)];
		double[][] gen1 = new double[1+len * (offset + coding1.length)][data.getRowDimension()];
		if(keepOrig){
			System.arraycopy(this.probes.toArray(new String[0]), 0, pi, 0, len+1);
			System.arraycopy(genotypes, 0, gen1, 0, len+1);
		}
		for (int k1 = 0; k1 < coding1.length; k1++) {
			int[] coding = coding1[k1];
			StringBuffer cod = new StringBuffer();
			for (int kj = 0; kj < coding.length; kj++) {
				cod.append("." + coding[kj]);
			}
			for (int j = 1; j < genotypes.length; j++) {
				int j1 = (k1 + offset) * len + j;
				pi[j1] = probes.get(j) + "." + cod;
				//loc[j1] = this.loc.get(j);
				//maf[j1] = this.maf[j];
				for (int k = 0; k < genotypes[j].length; k++) {
					double gen = genotypes[j][k];
					if (!Double.isNaN(gen)) {
						int lb = (int) Math.floor(gen);
						int ub = (int) Math.ceil(gen);
						if (ub == lb) {
							gen1[j1][k] = coding[ub];
						} else {
							double w1 = (gen - lb);
							double w2 = ub - gen;
							gen1[j1][k] = (coding[lb] * w1 + coding[ub] * w2);
						}
					}
					else{
						gen1[j1][k] = Double.NaN;
					}
				}
			}
		}

		this.data = (new Array2DRowRealMatrix(gen1)).transpose();
		;
		this.probes = new ArrayList<String>(Arrays.asList(pi));
		//this.loc = new ArrayList<String>(Arrays.asList(loc));
		//this.maf = maf;
	}


	public void replace(String string, String string2) {
		int ind = probes.indexOf(string);
		if(ind>=0)probes.set(ind, string2);
		
	}

	

	
	
	
}
