package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import conversion.CompressDir.OutputStreamWriter1;

public class ConvertLong {
	
	public static void main(String[] args){
		try{
		//	if(true)System.exit(0);
			File inputDir = new File(args[0]);
			File dir = new File(System.getProperty("user.dir"));
		
			String[] fields = "SNP Name:Sample ID:Chr:Sample Name:Position".split(":");
		System.err.println("fields "+Arrays.asList(fields));
			//	"ProbeName:null:SystematicName:null:SystematicName".split(":");
			String[] toKeepFields = null;//"LogRatio".split(":");
			File	in1 = new File(inputDir, args[1]);
			
			File extFile = args.length >3 ? new File(args[3]) : null;
			ConvertLong cl = new ConvertLong(dir, in1,extFile, 
					args[2],fields, toKeepFields);
			cl.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	public static String[] getKaryoBand(File f, String band) throws IOException{
		BufferedReader br =new BufferedReader(new FileReader(f));
		String st = br.readLine();
		int ind = band.indexOf('p');
		if(ind<0) ind = band.indexOf('q');
		String chr = band.substring(0,ind);
		String nme = band.substring(ind);
		String[] str = st.split("\t");
//		#chrom  chromStart      chromEnd        name    gieStain
		int start=-1;
		int end = -1;
		while((st = br.readLine())!=null){
			str = st.substring(3).split("\t");
			if(str[0].equals(chr)){
				if(str[3].startsWith(nme)){ 
					if(start<0){
						start = Integer.parseInt(str[1]);
						end = Integer.parseInt(str[2]);
					}else{
						end = Integer.parseInt(str[2]);
					}
				}
			}
		}
		return new String[] {chr,start+"",end+""};

	}
	
	File in;
	List<String> snps;
	File outDir;
	int offset=-2;
	final PrintWriter err;

	Map<String, OutputStreamWriter1> pw ;
	
	OutputStreamWriter1 samples;
	
	
	String sep = "\t";
	
	static String defaul = "SNP Name	Sample ID	Allele1 - Top	Allele2 - Top	GC Score	Theta	R	X	Y	X Raw	Y Raw	" +
	"B Allele Freq	Log R Ratio	CNV Value	CNV Confidence";  //this should normally not be used, just here for ep data
	
	
	static String[] empty_line = ("Sample ID	Array Type	SNP rsID	" +
			"Chromosome	Position	NaN	NaN	" +
			"NaN	NaN	NaN	NaN").split("\t");
	public static String filter ="";//"genotype";
	public static int header = -1; //number of lines to discard at top of zip entries

	
	//int start_col_index=6;  //6
	int[] cols;
	int id_index = 0;
	int snp_index = 0;
	int snp_index1 = 1;
	
	
	BufferedReader in1, in_snps;
	String st = "";
	String st_snp = "";
	OutputStreamWriter1 SNPS;
	int start,end;
	public ConvertLong(File dir, File in, File extFile,  String chrom
			, String[] fields, String[] toKeepFields) throws Exception{
		//offset = new int[in.length];
		
		outDir = new File(dir,chrom);
		err = new PrintWriter(new FileWriter(new File(dir,"log.txt")));
		outDir.mkdir();
		File karyo = new File(dir,"karyo_b36.txt");
		if(karyo.exists()){
			String[] band = getKaryoBand(karyo,chrom);
			this.chrom = band[0];
			start = Integer.parseInt(band[1]);
			end = Integer.parseInt(band[2]);
		}else{
			start = 0;
			end = Integer.MAX_VALUE-1;
			this.chrom = chrom;
		}
		//this.start_col_index = start_col_index;
		this.id_index = id_index;
		//this.snp_index = snp_index-1;
		this.snp_index1 = snp_index;
		this.in = in;
	
		this.compress = (new CompressDir(this.outDir));
	
				
	
		
		String snpName = fields[0];
		String sampleId = fields[1];
		String chrSt = fields[2];
		String sampleName = fields[3];
         String pos = fields[4];
	//	for(int k=0; k<offset.length; k++){
			 in1 = Utils.getBufferedReader(in);
			 in_snps = Utils.getBufferedReader(extFile);//Utils.getBufferedReader(in[k]);
	
			inner: for(int j=0; j<20 && (st=in1.readLine())!=null; j++){
				if(st.indexOf(snpName)>=0 && (st.indexOf(sampleName)>=0 ||st.indexOf(sampleId)>=0 ) ){
						offset = j;
						//String st1 = in1.readLine();
						//System.err.println(st1);
					System.err.println("offsets "+offset);
					List<String> l = Arrays.asList(st.split(sep));
					NAstring = new String[l.size()];
					Arrays.fill(NAstring,"NA");
					this.snp_index1 = l.indexOf(snpName);
					this.id_index = l.indexOf(sampleId);
					int chr_ind = l.indexOf(chrSt);
					List<Integer> a = new ArrayList<Integer>();
					if(snp_index1>=0) a.add(snp_index1);
					int tmp = Math.min(l.indexOf(sampleName), this.id_index);
					if(tmp<0) tmp = Math.max(l.indexOf(sampleName), this.id_index);
					if(tmp>=0) a.add(tmp);
					if(chr_ind>=0) a.add(chr_ind);
    				Collections.sort(a);
					if(snp_index<0 && extFile!=null){
						List< String> l2 = new ArrayList<String>();
						for(int kk=0; kk<a.size(); kk++){
							l2.add(l.get(a.get(kk)));
						}
					//	List< String> l2 = Arrays.asList(new String[] {a[0]>=0 ? l.get(a.get()) : null,a[1]>=0 ?  l.get(a[1]) : null, a[2]>=0 ? l.get(a[2]) : null});
						this.snp_index = l2.indexOf(snpName);
						System.err.println("snp_index is "+snp_index);
					//	snp_index1 = snp_index;
						
					}
					List< String> l2;
					if(toKeepFields!=null){
						l2 = Arrays.asList(toKeepFields);
					}
					else{
					l2 = new ArrayList<String>(l);
						l2.remove(sampleName);
						l2.remove(sampleId);
						l2.remove(snpName);
						l2.remove(pos);
						l2.remove(chrSt);
					}
					cols = getCols(l2,l);
					defaul=st;
					break inner;
				}
			}
			
			in1.close();
		//}
		OutputStreamWriter1 name = compress.getWriter("Name", true);
		name.printLine(defaul.split(sep), cols);
		name.println("chr	start	end	id	A	B");
		name.println("id");
		name.close();
		
		samples = compress.getWriter("Samples", false);	
		SNPS = compress.getWriter("SNPS", false);	
	}
String[] NAstring;
	
	private int[] getCols(List<String> l2, List<String> l) {
		int[] res = new int[l2.size()];
		for(int k=0; k<res.length; k++){
			res[k] = l.indexOf(l2.get(k));
		}
		return res;
	}


	//Sample ID       SNP Name        Chr     Position        Allele1 - Forward       Allele2 - Forward       X Raw   Y Raw   X       Y 
    //Theta   B Allele Freq   Log R Ratio


	public void close() throws Exception{
		SNPS.close();
		samples.close();
		compress.close();
		err.close();
	}
	
	CompressDir compress;
	final String chrom;
	
	
	
	public void run() throws Exception{
		st_snp = in_snps.readLine();
	   for(int k=0;st_snp!=null ;k++){
			try{
				run(k);
				//st_snp = in_snps.readLine();
			}catch(Exception exc){
				err.println(exc.getMessage());
				err.flush();
				exc.printStackTrace();
			}
	   }
		this.close();

		

	}
	
	
	
	
	
  List<String> samplesL = new ArrayList<String>();
	
	public void run(int k) throws Exception{
		int[] cols1 = new int[] {1,2,-1,0};
		pw =  new HashMap<String, OutputStreamWriter1>();
		
		//String firstSNP = null;
		while(st_snp!=null && pw.size()<maxsize){
			String[] strs = st_snp.split(sep);
			if(strs[1].equals(chrom) ){
				if(Integer.parseInt(strs[2])>=start){
					if(Integer.parseInt(strs[2])<=end){
				    //  if(firstSNP==null) firstSNP = strs[0];
				
						 SNPS.printLineSNP(strs, Integer.parseInt(strs[2])+40,cols1);
			
				      this.pw.put(strs[0], compress.getWriter(strs[0], false));
				//System.err.println(strs[0]+" "+pw.size());
			}
				}
			}
			st_snp = in_snps.readLine();
		}
		
		in1 = Utils.getBufferedReader(in);
			for(int j=0; j<=offset; j++){
				this.st = in1.readLine();
			}
			//System.err.println("h");
			System.err.println(st);
			int kk;
		//boolean writeSamples = k==0;
		for( kk=0; (st=in1.readLine())!=null; kk++){
			String[] str = st.split(sep);
			if(k==0 && !samplesL.contains(str[id_index])){
				System.err.println("added "+str[id_index]);
				System.err.println(st);
	    		samplesL.add(str[id_index]);
	    		samples.println(str[id_index]);
	    		samples.flush();
	    	}
			//boolean writeSamples  = str[0].equals(firstSNP);
			OutputStreamWriter1 osw = this.pw.get(str[0]);
		    if(osw!=null){
		    		String sampleCheck = samplesL.get(osw.curr_index());
		    		int ind1 = samplesL.indexOf(str[id_index]); //target index
		    		while(osw.curr_index()<ind1){
		    			err.println("adding NA "+str[0]+" "+osw.curr_index()+" "+samplesL.get(osw.curr_index()));
		    			osw.printLine(this.NAstring,cols);
		    		}
		    		if(osw.curr_index()>ind1){
		    			err.println(samplesL);
		    			err.println(samplesL.size());
		    			err.println(st);
		    			
		    			err.println("not matched "+str[id_index]+" "+sampleCheck+" "+osw.curr_index()+" "+ind1);
		    			err.println("excluding "+st);
		    			err.flush();
		    			//System.exit(0);
		    		}else{
		    			osw.printLine(str,cols);
		    		}
		    }
		}
		for(Iterator<OutputStreamWriter1> it = this.pw.values().iterator(); it.hasNext();){
			OutputStreamWriter1  osw1 = it.next();
			while(osw1.curr_index()<samplesL.size()){
			//	System.err.println("adding NA "+osw1.toString()+" "+osw1.curr_index()+" "+samplesL.get(osw1.curr_index()));
				osw1.printLine(this.NAstring,cols);
			}
			osw1.close();
		}
		System.err.println("done with "+pw.keySet());
		
		System.err.println("read "+k+" "+kk);
	}
	static int maxsize = 80;
	
	
}
