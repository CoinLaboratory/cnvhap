package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Executors;
import java.util.logging.Logger;

public class CompressIlluminaTable {
//    final String chromToDo;
   // final PropertyChangeSupport pcs;
	
	static String ref_string = null;//"0"; need to change this for chrM conversion!!!
	
	final boolean append;
    public static void main(String[] args){
        try{
        	//if(ConvertIllumina.es==null){
        		ConvertIllumina.es = OptionBuild.numThreads== 0 ? null:   Executors.newFixedThreadPool(OptionBuild.numThreads);;
    		//	
        	
        	 File dir1 = new File(System.getProperties().getProperty("user.dir"));
        	    File[] dir = dir1.listFiles(new FileFilter(){
        	        public boolean accept(File pathname) {
        	           return pathname.getName().endsWith("txt");
        	        }
        	    });
        	   File op = null;
        	  //  for(int i=0; i<dir.length; i++){
        	    	CompressIlluminaTable chmp = new CompressIlluminaTable(dir, new File("output"), 
        	    			new File(OptionBuild.SNPS),
        	    			"build36", false);
        	    	 chmp.run();
        	    	 op = chmp.dir;
        	  //  }
               
                CompressDir1.compress(op);
                ConvertIllumina.es.shutdown();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    public static File run(File[] dir, File outputDir, final String build, final String phenFile){
    	try{
    	
           CompressIlluminaTable chmp = new CompressIlluminaTable(dir, outputDir,	new File(OptionBuild.SNPS), build, false);
         chmp.run();
         CompressDir1.compress(chmp.dir);
        // Utils.delete(chmp.dir);
           return chmp.dir;
    	}catch(Exception exc){
    		exc.printStackTrace();
    		return null;
    	}
    }
    
    File dir;
    List<String>[] ids;
    List<Integer>[] idOffset;
    
    PrintWriter osw;
   final Map<String, PrintWriter> snps = new HashMap<String, PrintWriter>();
    
   // Map<String, String[]> samples = new HashMap<String, String[]>();
   
   // List<String> samplesL = null;
    
  
    	public void readSNP(){
    		
    	}
    	
	    
	//    String prefix;
    	 public File writeHeader(String chr) throws Exception{
             File out = new File(dir, chr);//
             if(!out.exists()) {
	             out.mkdir();
	             
	             {
		             File nme = new File(out, "Name");
		             PrintWriter osw = new PrintWriter(new BufferedWriter(new FileWriter(nme)));
		            
		             for(int i=0; i<header.length; i++){
		                 osw.write(header[i]);
		                 osw.write("\n");
		             }
		            osw.flush();
		           osw.close();
	             }
	             {
		             File nme = new File(out, "Samples");
		             PrintWriter osw = new PrintWriter(new BufferedWriter(new FileWriter(nme)));
		            for(int j=0; j<ids.length; j++){
		             for(int i=0; i<ids[j].size(); i++){
		            	 
		                 osw.write(ids[j].get(i));
		                 osw.write("\n");
		             }
		            }
		            osw.flush();
		           osw.close();
	             }
             }
             if(snps.get(chr)==null){
            	 PrintWriter snps = new PrintWriter(new BufferedWriter(new FileWriter(new File(out, "SNPS1"))));
	             this.snps.put(chr, snps);
             }
             return out;
         //}
     }
	   
	 //   boolean finished = false;
	    public void finish() throws Exception{
	    	for(Iterator<PrintWriter> it = this.snps.values().iterator(); it.hasNext();){
	    	 it.next().close();
	    	
	    	}
	    }
	  
	      
	     
	    
	     
	     
	     
		
    
    ///String chr;
    List<String>[] order; //ordder of info on snps
   // String currentString;
   // String[][] currentSpl;
 //  final  String currentChrom;
   List<String>[] inputHeader;
    BufferedReader[] br;
    int chrom_index;
    final String[] header;
     File[] f;
     List<String> order_snp = new ArrayList<String>();
// int b_index;
 //int r_index;
 //int geno_index;
 //int top_index;
 int snp_index;
 int ref_index, info_index;
 int loc_index ; 
 
 List<Integer> rem_indices = new ArrayList<Integer>();
 
 String currentChom;
 String[] currentString;
 final String[][] stri;
private String overallChrom=null;
BufferedReader SNPS = null;
String snpLine = null;
 boolean readLine() throws Exception{
	 for(int i=0; i<f.length; i++){
		 currentString[i] = br[i].readLine();
		
	 }
	 if(SNPS!=null) snpLine = SNPS.readLine();
	return (currentString[0]!=null);
 }
 List<String> SNPHeader;   
 
 CompressIlluminaTable(File[] f, File output, File SNPs, String build, boolean append) throws Exception{
    	header =   new String[] {
    	            "Genotype\tB allele\tLog R\tTop allelles",
    	            "chr\t"+build+"_start\t"+build+"_end\t id",
    	            "id"};
    	if(SNPs.exists()) this.SNPS = new BufferedReader(new FileReader(SNPs));
    	this.append = append;
   // 	this.chromToDo= todo;
    	this.dir = output;
    	Utils.delete(dir);
        dir.mkdir();
        
        this.f =f;
     //   this.currentChrom = todo;
        br =new BufferedReader[f.length];//
        currentString = new String[f.length];
        this.stri = new String[f.length][];
        for(int i=0; i<br.length; i++){
        	System.err.println("opengin "+f[i].getAbsolutePath());
        	br[i] =Utils.getBufferedReader(f[i]);
        }
      //  this.currentSpl = new String[br.length][];
      
       readLine();
    //   readLine();
       // for(int k=0; k<br.length; k++){
        
        	inputHeader = new List[f.length];
        	for(int i=0; i<inputHeader.length; i++){
        		inputHeader[i] = Arrays.asList(currentString[i].split(OptionBuild.splStr));
        	}
        		if(SNPS!=null) SNPHeader = Arrays.asList(this.snpLine.split(OptionBuild.splStr));
        	if(SNPHeader==null) SNPHeader = inputHeader[0];
        		chrom_index = SNPHeader.indexOf("chr");
        		snp_index = SNPHeader.indexOf("Name");
        		loc_index = SNPHeader.indexOf("Position");
        		ref_index = SNPHeader.indexOf("ref");
        		info_index = SNPHeader.indexOf("info");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("rsID");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("CNV");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("rs#");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("X");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("ProbeName");
        		if(snp_index<0) snp_index = SNPHeader.indexOf("PATIENT");
        		if(loc_index< 0) loc_index = SNPHeader.indexOf("Coord");
        		if(loc_index< 0) loc_index = SNPHeader.indexOf("position_b36");
        		if(loc_index< 0) loc_index = SNPHeader.indexOf("pos");
        		if(loc_index< 0) loc_index = SNPHeader.indexOf("Start");
        		if(loc_index< 0) loc_index = SNPHeader.indexOf("start");
        		if(chrom_index< 0) chrom_index = SNPHeader.indexOf("ChrName");
        		if(chrom_index< 0) chrom_index = SNPHeader.indexOf("chr");
        		if(chrom_index< 0) chrom_index = SNPHeader.indexOf("Chr");
        		if(chrom_index<0){
        			String prefix = f[0].getName();
        			   int indofChr = prefix.indexOf("chr");
        			   int end = Math.min(prefix.indexOf('.', indofChr),prefix.indexOf('_', indofChr));
        			   if(end<0) end =  prefix.length();
        			   if(indofChr<0){
        				   overallChrom = prefix.split("\\.")[0];
        			   }
        			   else{
        		        this.overallChrom = prefix.substring(indofChr+3,
        		        		end);
        			   }
        		}
       // }
        	//	if(snp_index<0) snp_index = loc_index;
           StringBuffer  rem = new StringBuffer();
        	 order = new List[f.length];//new ArrayList<String>(); //has order of b_freq log r, etc
        
        	ids = new List[f.length];//new ArrayList<String>();
        	idOffset = new List[f.length];//new ArrayList<Integer>();
        	String[] dataIndicators = new String[] {".","+","-"};
        	int i=0;
        	
        	List<Integer> captured_inds = new ArrayList<Integer>();
        	captured_inds.add(loc_index); captured_inds.add(chrom_index); captured_inds.add(snp_index); captured_inds.add(ref_index); 
        	for(;matches(inputHeader[0].get(i), dataIndicators)<0; i++){
        	//&& (inputHeader[0].get(i).indexOf("NA")<0) && inputHeader[0].get(i).indexOf('+')<0;i++){
        			this.order_snp.add(inputHeader[0].get(i));
        			if(!captured_inds.contains(i)){
        				rem_indices.add(i);
        				rem.append(inputHeader[0].get(i)+"\t");
        			}
        		
        	}
        	
        	
        	this.alleles = new boolean[inputHeader.length][];
        	boolean hasAllele = false;
        	for(int j=0; j<inputHeader.length; j++){
        	
        		order[j] = new ArrayList<String>();
        		ids[j] = new ArrayList<String>();
        		idOffset[j] = new ArrayList<Integer>();
		        for(int k=i; k<inputHeader[j].size(); k++){
		        	String stri = inputHeader[j].get(k);//;
		        	int index = stri.lastIndexOf('.');
		        	if(index<0) index = stri.lastIndexOf('-');
		        	String[] st;
		        	if(index>=0) st = new String[] {stri.substring(0,index).replace(' ','_'), stri.substring(index+1, stri.length())};
		        	else st = new String[] {stri, OptionBuild.default_type};
		        	
		        	if(!this.ids[j].contains(st[0])){
		        		idOffset[j].add( k);
		        		ids[j].add(st[0]);
		        	}
		        	String type =  st[1];
		        	if(!order[j].contains(type)){
		        		order[j].add(type);
		        	}
		        }
		        alleles[j] = new boolean[order[j].size()];
		        for(int kk=0; kk<alleles[j].length; kk++){
		        	alleles[j][kk] = order[j].get(kk).indexOf("Alleles")>=0;
		        	hasAllele = alleles[j][kk] || hasAllele;
		        }
		        if(hasAllele) this.header[1] =  "chr\t"+build+"_start\t"+build+"_end\tid\tAlleleA\tAlleleB";
		        if(j>0 && !order[j].equals(order[0])) throw new RuntimeException("cols not in same order between input files");
        	}
        
      
       header[0] = order[0].toString().replace(',', '\t');
       header[0] = header[0].substring(1,header[0].length()-1);
     String remst = rem.toString();
      if(remst.length()>0) header[1] = header[1]+("\t"+remst);
      
        
    }
    
    
    
   



	
  
  
  
  
    private int matches(String st, String[] dataIndicators) {
	for(int k=0; k<dataIndicators.length; k++){
		if(st.indexOf(dataIndicators[k])>=0) return k;
	}
	return -1;
}

	public void split(){
	  for(int i=0; i<currentString.length; i++){
		  stri[i] = currentString[i].split(OptionBuild.splStr);
		/*  if(i>0 && !stri[i][loc_index].equals(stri[0][loc_index])) {
			  throw new RuntimeException("files must be in same order");
		  }*/
	  }
	  if(this.snpLine!=null){
		  strisnp = snpLine.split(OptionBuild.splStr);
	  }
	  else strisnp = this.stri[0];
  }
    String[] strisnp;
   int overall_index=0;
 //  String alleleA, alleleB;
  SortedSet<Character> alleleSet = new TreeSet<Character>();
final boolean[][] alleles;
   public void run() throws Exception{
	   //	String str = "";
	    outer: while(this.readLine()){
		overall_index++;
		alleleSet.clear();
		//if(strisnp==null || this.strisnp[0].indexOf("reserved")>=0) continue outer;
	  //   String[][] stri = split();
	    	this.split();
	     //System.err.println(Arrays.asList(this.currentString));
	     //System.err.println(chrom_index);
	     String chrom = chrom_index>=0 ? this.strisnp[this.chrom_index] : this.overallChrom;
	     File out = writeHeader(chrom);
	     PrintWriter snps = this.snps.get(chrom);
	     String ref = ref_index>=0 ? strisnp[ref_index] : ref_string;
	     String info = info_index >=0 ? strisnp[info_index] : null;
	     
	     String snpid = snp_index>=0  ? strisnp[snp_index] :
	    	( loc_index>=0  && chrom_index >=0? strisnp[chrom_index]+"_"+strisnp[loc_index] :overall_index+"");
	     File outF = new File(out, snpid);
	 	if(!append){
	 		
	 		for(int i=0; outF.exists(); i++){
	 			System.err.println("input file must be sorted according to rs number "+snpid);
	 			outF = new File(out, snpid+"_"+i);
	 		}
	 	}
	     PrintWriter osw = new PrintWriter(new FileWriter(outF,append));
	   String st =  loc_index >= 0 ? strisnp[loc_index] : overall_index+"";
	   st = st.split("\\.")[0].split("-")[0];
      String end=st;
      try{
    	  end = "" +(Integer.parseInt(st)+20);
      }
      catch(Exception exc){
    	  
      }
        for(int j=0; j<f.length; j++){
	     for(int i=0; i<idOffset[j].size(); i++){
	    		  int offset = idOffset[j].get(i);
	    			  for(int kk=0; kk<order[j].size(); kk++){
	    				if(alleles[j][kk]){
	    					char[] ch = stri[j][offset+kk].toCharArray();
	    					for(int ikk =0; ikk<ch.length; ikk++){
	    					alleleSet.add(ch[ikk]);
	    					}
	    				}
	    				String top = stri[j][offset+kk];
	    				if(ref!=null){
	    					top = top.equals("-") ? top : top.equals(ref) ? "A" : "B";
	    				}
	    				  osw.print(top);
	    				 if( kk<order[j].size()-1 ){osw.print("\t");}
	    				 else osw.print("\n");
	    				//  osw.print(? "\t":"\n");
	    		
	    			  }
	     }
		 
      }
        //System.err.println(snpid+" "+alleleSet);
       StringBuffer rems = new StringBuffer();
       if(this.rem_indices.size()>0){
    	   for(Iterator<Integer> it = rem_indices.iterator(); it.hasNext();){
    		   rems.append("\t"+strisnp[it.next()]);
    	   }
       }
        if(info!=null){
        	//if(info.indexOf('G')>=0){
        	//	System.err.println("h");
        	//}
        	String[] allele = info.replaceAll("!","").replaceAll("\\(","").replaceAll("\\)", "").split("[0-9]");
        	info = info.indexOf('!')>=0 ? "FALSE" : "TRUE"; //for chrM data only
        	alleleSet.add('-');
        	alleleSet.add(allele.length>0 && allele[allele.length-1].length()>0 ? allele[allele.length-1].charAt(0):'-');
        }
        if(alleleSet.size()>0 ){
        	 snps.println("chr"+chrom+"\t"+st+"\t"+end+"\t"+snpid+"\t"+alleleSet.first()+"\t"+alleleSet.last()+(info!=null ? ("\t"+info):"")+rems.toString());
        }
        else{
         snps.println("chr"+chrom+"\t"+st+"\t"+end+"\t"+snpid+(ref!=null ? "\t"+ref :"")+(info!=null ? ("\t"+info):"")+rems.toString());
        }
       
	     osw.close();
	     
   		
}
   this.finish(); 
   this.sortSNPs();
}

public void sortSNPs() throws Exception{
	File[] f = this.dir.listFiles(new FileFilter(){

		public boolean accept(File pathname) {
			return pathname.isDirectory();
		}
		
	});
	for(int i=0; i<f.length; i++){
		File snpsinF = new File(f[i], "SNPS1");
		BufferedReader snps_in = new BufferedReader(new FileReader(snpsinF));
		String st = "";
		List<Comparable[]> l = new ArrayList<Comparable[]>();
		while((st = snps_in.readLine())!=null){
			String[] str = st.split("\t");
			Comparable[] res = new Comparable[str.length];
			for(int i1=0; i1<str.length; i1++){
				res[i1] = i1==1 ? (str[i1]==null || str[i1].equals("NA") || str[i1].equals("null")? -1 :Integer.parseInt(str[i1])) : str[i1];
			}
			l.add(res);
		}
		Collections.sort(l, compa);
		PrintWriter snps_out = new PrintWriter(new BufferedWriter(new FileWriter(new File(f[i], "SNPS"))));
		for(Iterator<Comparable[]> it = l.iterator(); it.hasNext();){
			Comparable[] res = it.next();
			snps_out.print(res[0]);
			for(int i1=1; i1<res.length; i1++){
				snps_out.print("\t");
				snps_out.print(res[i1]);
			}
			snps_out.print("\n");
		}
		snps_out.close();
		snpsinF.delete();
	}
}
final static Comparator<Comparable[]> compa = new Comparator<Comparable[]>(){

	public int compare(Comparable[] o1, Comparable[] o2) {
		return o1[1].compareTo(o2[1]);
	}
	
};
	


}