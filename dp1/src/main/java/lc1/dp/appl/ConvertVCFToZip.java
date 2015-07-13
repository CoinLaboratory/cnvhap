package lc1.dp.appl;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import conversion.Utils;

import lc1.util.CompressDir;

public class ConvertVCFToZip {
	static boolean replaceDotWithZero = false; 
	static boolean replaceNumWithAB = false;
static String split="\\s+";
static int buffer =1;
static boolean includeRef = false; //SET TO TRUE TO ADD UP DEPTHS!!!
 static boolean useID = false;
static  boolean noskip = true;
//static boolean mergeNeighbourin = true;	

	public static void main(String[] args){
	 try{
		// if(true)System.exit(0);
		 File f = new File( System.getProperty("user.dir"));
		// File f1 = new File(f,args[0]);
		// if(!f1.exists()) f1.mkdir();
		 String[] args1 = args;
		 if(args.length>0){
		     args1 = args;
		 }
		 else{
			 args1 = new String[] {"","build36","Y"};
		 }
		 String headerSt = args1[2]; //"#chr";
		 Integer buffer = 1;
		 if(args1.length>3) buffer = Integer.parseInt(args1[3]);
				 //f[0].getName().indexOf("vcf")>=0 ?  "#chr" :"rs#";
		 main(f, f,false, args1[0], args1[1], args.length>4 && !args1[4].equals("null") ? args1[4]:null, args1.length>5 ? args1[5] : null, headerSt, buffer);
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
 }
	final String buildF;
	private int[] format_string_index;

private static void  getKary(File karyfile, String substring, List<Integer> res, List<String> app, String chrom) throws Exception{
	String st = "";
	BufferedReader br = getBR(karyfile);
	while((st = br.readLine())!=null){
		String[] str = process(st, chrom);//.split(split);
		if(str[0].equals(substring)){
			//Integer[] res = new Integer[str.length-1];
			for(int k=1; k<str.length; k++){
				res.add( Integer.parseInt(str[k]));
				if(k==1) app.add("p");
				else app.add(k+"");
			}
			res.add(Integer.MAX_VALUE-1);
			app.add("q");
			return;
		}
	}
	
	
}
public void close() throws Exception{
		 snps.close();
		 compress.run();
		 compress.close();
}
//String currSubDir="";
public void getNew() throws Exception{
	String prevChromT = currentChromT;
	if(chr_index>=0){
	currentChrom = this.currString[chr_index];
	}else{
		currentChrom="all";
	}
	currentChromT = currentChrom;
	
	/*if(prevChrom==null || !prevChrom.equals(currentChrom)){
		currSub=0;
	}else{
		currSub = currSub+1;
	}*/
	String nme = this.hasChrPrefix ? currentChrom.substring(3) : currentChrom;
	 File dir = new File(this.dir1,nme);//f[k].getName().split("\\.")[0]);
	 
	 dir.mkdir();
	// if(currSub==0){
		 if(snps_global!=null) snps_global.close();
		 File glob = new File(dir1,buildF+"_"+nme+".txt");
	this.snps_global = new PrintWriter(new FileWriter(glob));
			
	 //}
	 //currSubDir = currentChrom+"_"+currSub;
	 //if(this.hasChrPrefix)currSubDir = currSubDir.substring(3);
	 //File f = new File(dir,currSubDir);
     
	 compress = (new CompressDir(dir));
	
		snps = compress.getWriter("SNPS", false);
		String[] str = this.currString;
		if(format_index>=0){
			format_string = str[format_index];
		}
}


void writeSamples() throws Exception{
	String sample_line = "sample";
    {
   	String[] str = samples.get(0).split("=");
		for(int k=1; k<str.length; k++){
			sample_line = sample_line+"\tinfo_"+k;
		}
    }
	 OutputStreamWriter pw = compress.getWriter("Samples",true);
		for(int i=0; i<samples.size(); i++){
		
				
			
			pw.write(samples.get(i).replace('=', '\t')+"\t"+this.sums[i]+"\n");
		}
		compress.closeWriter(pw);
	 pw = 	compress.getWriter("Name",true);
	if(format_string!=null){
		String[] format = format_string.split(":");
		this.format_stringL = format_string.split(":");
		this.format_string_index = new int[format_stringL.length];
		
	
	//int[] alias = new int[header.size()];
		if(replaceNumWithAB) format[0] =  format[0].replaceAll("GT","genotype");
		for(int k=0; k<format.length; k++){
			pw.write(format[k]+(k<format.length-1 ? "\t":"\n"));
		}
	}else{
		pw.write("Genotype\n");
	}
	pw.write("chr\tstart\tend\tsnpid");
//		pw.write("countAll\tstate.0\tstate.1\tstate.2\n");
	/*	for(int k=0; k<offset; k++){
			if(k!=chr_index && k!=pos_index && k!=format_index){
				pw.write("\t"+header.get(k));
			}
		
		}*/
		if(this.alleles_index>=0){
			pw.write("\talleleA\talleleB");
		}
		pw.write("\n");
		//
		pw.write(sample_line+"\n");
		compress.closeWriter(pw);
}
 

ConvertVCFToZip(BufferedReader br, BufferedReader br1, File dir1, List<String> header, String buildF, String chrom, String maxcnt, String filename, int buffer) throws Exception{
	this.br = new BufferedReader1(br, buffer);
	this.br1 = br1;
	this.chrom = chrom;
	target = this.chrom;//==null ? "" : this.chrom;
	if(target==null) noskip=true;
	else noskip = false;
	this.dir1 = dir1;
	if(maxcnt!=null){
		String[] str = maxcnt.split(":");
		this.max_cnt = Integer.parseInt(str[0]);
		this.max_rep = Integer.parseInt(str[1]);
	}else{
		max_cnt = Integer.MAX_VALUE-1;
		max_rep = Integer.MAX_VALUE-1;
	}
	
	 id_index = header.indexOf("ID");
	 if(id_index<0) 	 id_index =header.indexOf("rs#");
	 if(id_index<0) 	 id_index =header.indexOf("Name");
	 if(id_index<0) 	 id_index =header.indexOf("#ID");
	 
	 chr_index = find(header,"chrom", true);
	 pos_index = find(header,"pos",true);
	 end_index = find(header,"end",true);
	

	 if(pos_index<0) pos_index = find(header,"start",true);
	 if(pos_index<0) pos_index = find(header,"MapInfo",true);
	 if(chr_index<0) chr_index = find(header,"Chr",true);
	 format_index = find(header,"format",true);
	 ref_index = header.indexOf("REF");
	 alleles_index = header.indexOf("alleles");
	 alt_index = header.indexOf("ALT");
	 if(format_index>=0){
		 offset = format_index+1;
	 }else if( header.indexOf("QCcode")>=0){
		offset =  header.indexOf("QCcode")+1;
	 }else if(find(header,"EndPosition",true)>=0){
			offset = find(header,"EndPosition",true)+1;
		
      }else if(find(header,"probe",true)>=0){
	offset = find(header,"probe",true)+1;
	 if(id_index<0) 	 id_index =offset-1;
 }
	 if(offset<0 && end_index>=0) offset = this.end_index+1;
	 if(offset<0) offset = this.pos_index+1;
	 this.header = header;
	 firstLine();
	 this.format_string = format_index>=0 ? this.currString[format_index] : "DP";
	
	this.buildF = buildF;
	if(offset<0) offset =  1+find(header,"endPos",true);
	 samples = new ArrayList<String>();
		for(int k=offset; k<header.size(); k++){
			samples.add( header.get(k));
		}
		if(includeRef){
			samples.set(samples.size()-1, "REFALL");
		}
		sums = new double[samples.size()];
}


private int find(List<String> header2, String string, boolean lower) {
	String mtc = lower? string.toLowerCase() : string;
	for(int i=0; i<header2.size(); i++){
		String st = header2.get(i);
		if(lower) st = st.toLowerCase();
		if(st.indexOf(mtc)>=0) return i;
	}
	return -1;
}
String[] dp_format = new String[] {"DP"};
String chrom;
final File dir1;
static int id_index, chr_index,  format_index, ref_index, alt_index, pos_index, alleles_index, end_index;
int offset = -1;
//int max_cnt = 1000;
class BufferedReader1 {
	int buffer;
	double[][] sums;
	double[]sumsall;
	String nxtLine;
	BufferedReader1(BufferedReader fir, int buffer) throws IOException{
		this.buffer= buffer;
		this.br =fir;
		nxtLine = br.readLine();
	}
	
	BufferedReader br;
	public String readLine() throws IOException{
		String st = nxtLine;
		if(st!=null) nxtLine = br.readLine();
		return st;
	}
	public String join(double[] d){
		StringBuffer sb = new StringBuffer(d[0]+"");
		for(int k=1; k<d.length; k++){
			sb.append(':'+d[k]);
		}
		return sb.toString();
	}
	public String readLine1() throws IOException{
	  
		String st = nxtLine;
		  if(sums==null){
			  String[] str = nxtLine.split(split);
			  String[] format = format_index>=0 ?str[format_index].split(":") : dp_format;
			  sums = new double[str.length-offset+1][];
			  sumsall = new double[format.length];
			  for(int j=0; j<sums.length;j++){
				  sums[j] = new double[format.length];
			  }
		  }
		if(st==null) return null;
		else if(!includeRef && buffer==1){
			nxtLine = br.readLine();
			return st;//.split(split);
		}
		String chrom=""; 
		String[] str=null;
		int i=0;
		Arrays.fill(this.sumsall, 0);
		for(int k=0;k<sums.length; k++){
			Arrays.fill(sums[k],0);
		}
		String[] res=null;
		for(; i<buffer && st!=null && ((str =st.split(split))[chr_index].equals(chrom) || i==0 ); i++){
			if(i==0) {
				res = str;
				chrom = str[chr_index];
			}else{
				if(end_index>0) res[end_index] = str[end_index];
						
			}
              for(int k=offset; k<str.length; k++ ){
            	 if(!str[k].equals(".") && str[k].indexOf("/")<0){
            	 String[] vals = str[k].split(":");
            	 for(int j=0; j<vals.length; j++){
            		 double d = Double.parseDouble(vals[j]);;
            		 sums[k-offset][j] += d;
            		 sumsall[j]+=d;
            	 }
            	 }
      		   }
          	  st = br.readLine();
          	  System.err.println(st);
		}
		StringBuffer sb = new StringBuffer(res[0]);
		for(int k=1; k<offset; k++){	
			sb.append("\t");sb.append(res[k]);
		}
		for(int k=offset; k<str.length; k++ ){
 			 sb.append("\t"); sb.append(join(sums[k-offset]));
 		}
		if(includeRef) sb.append("\t"+join(sumsall));
		nxtLine = st;
		return sb.toString();
	}
	public void close() throws IOException {
		this.br.close();
		
	}
}

final BufferedReader1 br;
final BufferedReader br1;
String currentChrom;
String format_string="DP";

String[] format_stringL;

String currentChromT="NA";
List<String> header,samples;
double[] sums; //keeps a running depth sum across all sites
String[] currString;
//int currSub=-1;
boolean first = true;
OutputStreamWriter snps;
PrintWriter snps_global;
/** returns if chrom is same */

boolean hasChrPrefix = false;

//String[]res;
String[] readLine(){
	try{
   String	currentLine = br.readLine1();

  if(currentLine==null) {
	 return null;
 }
else if(br1!=null){//for broad omni intensities
	
	String[] currentLine0 = process(currentLine, chrom);//.split(split);
	String[] currentLine1 = process(br.readLine1(), chrom);//.split(split);
	String id1 = currentLine0[0].split("-")[0];
	String id2 = currentLine0[0].split("-")[0];
	if(!id1.equals(id2)) throw new RuntimeException(id1+" "+id2);
	/* String[] currentSNP = br1.readLine().split(split);
	 while(!currentSNP[0].equals(id1)){
		 System.err.println(currentSNP[0]);
		 currentSNP = br1.readLine().split(split);
		
	 }*/
	
	//int snpl = currentSNP.length;
	//int snpl =0;
	int cl= currentLine0.length;
	currentLine0[0] = id1;
	//String[] res = new String[snpl+cl];
//	System.arraycopy(currentSNP, 0, res, 0, snpl);
	for(int j = 1;j<currentLine0.length; j++){
		currentLine0[j] = currentLine0[j]+";"+currentLine1[j];
	}
	return currentLine0;
}
else{
	//System.err.println(currentLine.substring(0, 100));
	return process(currentLine, chrom);//.split(split);
}
	}catch(Exception exc){
		exc.printStackTrace(); return null;
	}
	
}

Boolean firstLine(){
	try{
	target = this.chrom;//==null ? "all" : this.chrom;
	currString = readLine();
	int h = header.size();
	for(int k=h; k<currString.length; k++){
		header.add("S"+(k-offset));
	}
	if(chr_index>=0){
    if(currString[this.chr_index].startsWith("chr")){
		hasChrPrefix = true;
		if(chrom!=null){
			chrom = "chr"+chrom;
		target = this.chrom;
		}else{
			chrom = currString[chr_index];
			 target = this.chrom;
		}
	}else if(currString[this.chr_index].startsWith("Chr")){
		hasChrPrefix = true;
		if(chrom!=null){
			chrom = "Chr"+chrom;
		    target = this.chrom;
		}else{
			chrom = currString[chr_index];
			 target = this.chrom;
		}
	}
	
	}else{
		target = "all";
	}
	//if(target==null) target = "all";
	if(!target.equals("all") &&!target.equals("chrall") && !currString[chr_index].equals(target)){
		while((currString = readLine())!=null && !currString[chr_index].equals(target)){}
	}
	if(currString==null) return null;
	return target.equals("all") || currString[chr_index].equals(currentChromT);
	}catch(Exception exc){
		exc.printStackTrace();
		return null;
	}
	
}
 String target;
Boolean nextLine(){
	try{
	while((currString = readLine())!=null){
		if(noskip || target.equals("all") || target.equals("chrall") || currString[chr_index].equals(target)  ){
			break;
		}else{
			System.err.println(Arrays.asList(currString));
		}
	//	System.err.println(Arrays.asList(currString));
	}
	if(currString==null){
		return null;
	}
	return target.equals("all") || target.equals("chrall") ||currString[chr_index].equals(currentChromT);
	}catch(Exception exc){
		exc.printStackTrace();
		return null;
	}
	
}

public void run() throws Exception{
	Boolean nxtLine = true;
	int i=0;
	for(; nxtLine!=null && i<max_rep;i++){
	 this.getNew();
	// this.currentChrom.equals
	 nxtLine = this.mainTranspose();
	 writeSamples();
	 close();
	 if(currString==null){
		 break;
	 }else if(chr_index>=0){
		 currentChromT = currString[chr_index];
		 target = currentChromT;
	 }
	 if(currentChromT.equals(target) && target.length()>0 &&(chr_index<0 || !currString[chr_index].equals(target))){
		 break;
	 }
	}
	this.snps_global.close();
	this.br.close();
}
final int max_cnt,max_rep;
private Collection<String> ids = new HashSet<String>();;

Boolean mainTranspose() throws Exception{
		String st = "";
		Boolean nxtLine = true;
		boolean reformat = false;;
		Arrays.fill(sums,0.0);
		outer: for(int cnt=0;  nxtLine!=null && nxtLine && cnt<max_cnt; cnt++){
			if(cnt % 1000  ==0) System.err.println(cnt/1000);
			try{
			String[] str =this.currString;
			
				if(format_index>=0 && !str[format_index].equals(format_string)){
					nxtLine = this.nextLine();
					continue outer;
					/*this.format_string = str[format_index];
					reformat = true;
					List<String> str_ = Arrays.asList(str[format_index].split(":"));
					for(int k=0; k<this.format_string_index.length; k++){
						format_string_index[k] = str_.indexOf(this.format_stringL[k]);
					}
//					throw new RuntimeException("format changed:\n"+str[format_index]+"\n"+format_string);*/
				}else{
					reformat = false;
				}
				String rsid;
				if(pos_index>=0){
			       int sta = (int) Math.round(Double.parseDouble(str[pos_index]));
			       int end = end_index>=0 ? (int) Math.round(Double.parseDouble(str[end_index])): sta+20;
		
			    
				if(!this.hasChrPrefix){
					rsid = str[chr_index]+"_"+str[pos_index];
					str[chr_index] = str[chr_index];
				}else{
					rsid = str[chr_index].substring(3)+"_"+str[pos_index];
				}
				if(this.id_index>=0 && useID) rsid = str[id_index];
				if(ids.contains(rsid)){
					String rsid1 = rsid;
					for(int i=0; ids.contains(rsid); i++){
						rsid = rsid1+"."+i;
					}
				}
				
			snps.write(str[chr_index]+"\t"+sta+"\t"+end+"\t"+rsid);
			//System.err.println((str[chr_index]+"\t"+sta+"\t"+end+"\t"+rsid));
			this.snps_global.write(str[chr_index]+"\t"+sta+"\t"+end+"\t"+rsid);
			for(int k=0; k<offset; k++){
				if(k!=chr_index && k!=pos_index && k!=format_index){
					if(k<ref_index) snps.write("\t"+str[k]);
					if(k<ref_index) snps_global.write("\t"+str[k]);
				}
				
			}
				}else{
					rsid = str[id_index];
					if(ids.contains(rsid)){
						System.err.println("warning, rsid already present "+rsid);
						String rsid1 = rsid;
						for(int i=0; ids.contains(rsid); i++){
							rsid = rsid1+"."+i;
						}
					}
					snps.write("NA"+"\t"+"NA"+"\t"+"NA"+"\t"+rsid);
				}
				ids.add(rsid);
			String[] alleles = null;
			if(alleles_index>=0){
				alleles = str[alleles_index].split("/");
				snps.write("\t"+alleles[0]+"\t"+alleles[1]);
			}
			snps.write("\n");
			snps_global.write("\n");
			String comment =  null;
			if(ref_index>=0 && alt_index>=0){
				comment = str[ref_index]+str[alt_index]+","+str[id_index];
			}
			OutputStreamWriter pw = compress.getWriter(rsid,true, comment);
			/*if(mergeNeighbouring){
				String[] str1 = new String[offset+(str.length-offset)/2];
				for(int k=offset; k<str.length; k=k+2){
					str1[offset+(k-offset)/2] = str[k]+":"+str[k+1];
				}
				str = str1;
			}*/
			for(int k=offset; k<str.length; k++){
				String[] str2 = str[k].split(":");
				if(format_string.startsWith("DP")){
					sums[k-offset]+=Double.parseDouble(str2[0]);
				}
				if(alleles!=null){
					str2[0] =str2[0].replace(alleles[0].charAt(0), 'A').replace(alleles[1].charAt(0), 'B');
					if(this.chrom.equals("chrY")){
						if(str2[0].charAt(0)!=str2[0].charAt(1)){
							System.err.println("warning at "+str[this.pos_index]+" "+samples.get(k-offset));
							str2[0] = "N";
						}
						else str2[0] = str2[0].charAt(0)+"";
					}
				}
				if(reformat){
					for(int j=0; j<this.format_string_index.length; j++){
						int j1 = format_string_index[j];
						pw.write((j1< 0 ? "." : str2[j1])+(j<format_string_index.length-1 ? "\t":"\n"));
						
					}
				}else{
					if(replaceNumWithAB){
						str2[0]   = str2[0].replaceAll("0", "A").replaceAll("1", "B").replaceAll("/","");
					}
					for(int j=0; j<str2.length; j++){
						if(replaceDotWithZero  && str2[j].equals(".")) str2[j] = "0";
						pw.write(str2[j]+(j<str2.length-1 ? "\t":"\n"));
					}
				}
			}
			compress.closeWriter(pw);
			nxtLine = this.nextLine();
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		
		return nxtLine;
}
static String chr;	
public static CompressDir compress=null; 

private static BufferedReader getBR(File f) throws Exception{
	return f.getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f)))) : 
		 new BufferedReader(new FileReader(f));
}

static boolean countQCMG = false;
static Map<String, Integer> chrToNumeric = new HashMap<String, Integer>();
private static  String[] process(String str, String chrom){
	if(countQCMG){
		if(str.startsWith("#")){
			 str = str.replace("Format", "Format").replace("Counts", "P100\tREFALL");
			 return str.split(split);
		}
		else{
			str = str.replace("Ref", "DP").replace(':', '\t');
			String[] st = str.split(split);
			 String chr_ = st[chr_index].replaceAll("chr", "");
			st[id_index] = chr_+"_"+st[pos_index];
          
         Integer chr = chrToNumeric.get(chr_);//chr_.indexOf('X')>=0 ? 23.0: chr_.indexOf('Y')>=0 ? 24 : chr_.indexOf('M')>=0 ? 25 : Double.parseDouble(chr_);
         if(chr==null){
        	 chr  = chrToNumeric.size()+1;
        	 chrToNumeric.put(chr_,chr);
        	 System.err.println(chr_+" "+chr);
         }
         if(chrom!=null && chrom.equals("all")){
          double pos1 = Double.parseDouble(st[pos_index]);
          double pos2 =Double.parseDouble(st[end_index]);
 
          int pos11 = (int)Math.round(((pos1/3e8)+chr.doubleValue())*5e7);
          int pos21 =(int)  Math.round(((pos2/3e8)+chr.doubleValue())*5e7);
          if(pos11 < 0 ){
        	  System.err.println("h");
          }
        st[chr_index] = "all";
          st[pos_index] = pos11+"";
          st[end_index] = pos21+"";
		}
			 return st;
		}
	}
    return str.split(split);
}


private static void main(File dir, File dirout, boolean remove, final String prefix,final String build, final String chrom, final String maxcnt, final String headerSt, int buffer) {
	 try{
		// if(true) return;
		//	File dir = new File(System.getProperty("user.dir")+"/"+args[0]);
	File  avg = dir;
	File avg1 = dirout;
	 File[] f = avg.listFiles(new FileFilter(){

		@Override
		public boolean accept(File pathname) {
		 return (pathname.getName().indexOf("_fwd")>=0  || pathname.getName().indexOf("vcf")>=0 || pathname.getName().indexOf("counts")>=0|| pathname.getName().indexOf("depth")>=0   )
		 && pathname.getName().startsWith(prefix) && !pathname.isDirectory();// && pathname.getName().endsWith(".gz");
		}
		 
	 });
	
	
//	 chr = "chr"+dir.getName().split("_")[2];
	 for(int k=0; k<f.length; k++){
		
		 if(f[k].length()==0) continue;
		
		 BufferedReader br = 
			 f[k].getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f[k])))) : 
			 new BufferedReader(new FileReader(f[k]));
			 
			 File snpsF = new File(f[k].getParent(), f[k].getName().substring(0,f[k].getName().length()-4)+".snps");
			   BufferedReader br1 = null;
			   List<String>headerSNP = new ArrayList<String>();
			 if(snpsF.exists()){
				    br1 = 	 new BufferedReader(new FileReader(snpsF));
						//headerSNP.addAll(Arrays.asList(br1.readLine().split("\\s+")));		
			 }
		 String str =  "";
		 if(f[k].getName().endsWith(".counts")){
			 countQCMG = true;
			 useID=true;
			 for(int i=1; i<=22; i++){
				 chrToNumeric.put(i+"", i);
			 }
			 chrToNumeric.put("X", 23);
			 chrToNumeric.put("Y", 24);
			 chrToNumeric.put("MT", 25);
		 }
		 if(f[k].getName().endsWith(".depth")){
			 headerSNP.addAll(Arrays.asList("#Chrom:POS".split(":")));
			 includeRef = true;
			 File list = new File(avg,"list");
			 if(list.exists());
			 headerSNP.addAll(Utils.readStringInfo(list, 0, false));
			System.err.println(headerSNP);
		 }
		 else{
		 while(!(str=br.readLine()).toLowerCase().startsWith(headerSt.toLowerCase())){
			 System.err.println(str);
		 }
		 String[] rsid =process(str, chrom);
		 headerSNP.addAll(Arrays.asList(rsid));
		 }
		 int ind = f[k].getName().lastIndexOf('.');
		// File avg1 = new File(avg.getParentFile(),avg.getName()+"1");
		// avg1.mkdir();
		// System.err.println(avg1.getAbsolutePath());
		 File dir1 = new File(avg1,f[k].getName().substring(0,ind));
		 if(buffer>1){
			 dir1 = new File(avg1,f[k].getName().substring(0,ind)+"."+buffer);
		
		 }
		 dir1.mkdir();
		 System.err.println(dir1.getAbsolutePath());
		 ConvertVCFToZip cvz = new ConvertVCFToZip(br,br1, dir1, headerSNP, build, chrom, maxcnt, f[k].getName(), buffer);
		 cvz.run();
			if(remove)f[k].delete();
	 }
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
	
 }
}
