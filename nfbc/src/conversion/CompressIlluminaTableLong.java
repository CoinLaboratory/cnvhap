package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Executors;

import lc1.util.CompressDir;

public class CompressIlluminaTableLong {
//    final String chromToDo;
   // final PropertyChangeSupport pcs;
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
        	    for(int i=0; i<dir.length; i++){
        	    	CompressIlluminaTableLong chmp = new CompressIlluminaTableLong(dir[i], new File("output"), "build36", i>0);
        	    	 chmp.run();
        	    	 op = chmp.dir;
        	    }
               
                CompressDir.compress(op);
                ConvertIllumina.es.shutdown();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    public static void run(File dir, File outputDir, final String build, final String phenFile){
    	try{
    	
           CompressIlluminaTableLong chmp = new CompressIlluminaTableLong(dir, outputDir, build, false);
         chmp.run();
        //   CompressDir.compress(chmp.dir);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    }
    
    File dir;
    List<String> ids;
   // List<Integer> idOffset;
    
    //PrintWriter osw;
    OutputStreamWriter snps;
   //final Map<String, PrintWriter> snps = new HashMap<String, PrintWriter>();
    
   // Map<String, String[]> samples = new HashMap<String, String[]>();
   
   // List<String> samplesL = null;
    
  
    	public void readSNP(){
    		
    	}
    	
	    
	//    String prefix;
	  
    	public void finish() throws Exception{
    		if(compress !=null){
       		 compress.closeWriter(snps);
       		 compress.run();
       	 }
    	}
	    public void writeHeader(String chr) throws Exception{
	             File out = new File(dir, chr);//
	             
	             if(compress==null || !compress.inDir.equals(out)){
	            	 this.finish();
	            	 if(out.exists()) throw new RuntimeException("chroms must be in order");
	            	 compress = (new CompressDir(out));
	            	 OutputStreamWriter osw = compress.getWriter("Name", true);
			             for(int i=0; i<header.length; i++){
			                 osw.write(header[i]);
			                 osw.write("\n");
			             }
			             compress.closeWriter(osw);
			          // File nme = new File(out, "Samples");
			             osw = compress.getWriter("Samples", true);
			            
			               for(int i=0; i<this.ids.size(); i++){
			               	osw.write(ids.get(i)+"\n");
			               }
			               compress.closeWriter(osw);
			           snps = compress.getWriter("SNPS", false);  
	             }
	         //}
	     }
	 //   boolean finished = false;
	    
	  
	      
	     
	    
	     
	     
	     
		
    
    ///String chr;
  //  List<String> order; //ordder of info on snps
   // String currentString;
   // String[][] currentSpl;
 //  final  String currentChrom;
   List<String> inputHeader;
    BufferedReader br;
    int chrom_index;
    int sample_index;
    int[] toprint_ind;
    final String[] header;
     File f;
     List<String> order_snp = new ArrayList<String>();
// int b_index;
 //int r_index;
 //int geno_index;
 //int top_index;
 int snp_index;
 int loc_index ; 
 String currentChom;
 String currentString;
 
 Map<String, String[]>output = new HashMap<String, String[]>();
 String rs = "";
 
 CompressDir compress;
 
 public void parseNext() throws Exception{
	 this.output.clear();
 	String[] str = currentString.split(OptionBuild.splStr);
 	String rs = str[snp_index];
 	  String chrom = chrom_index<0 ? "null" : str[this.chrom_index];
	     writeHeader(chrom);
	    // PrintWriter snps = this.snps.get(chrom);
	   
	     int st = loc_index < 0 ? 0 :Integer.parseInt(str[loc_index]);
     snps.write("chr"+chrom+"\t"+st+"\t"+(st+20)+"\t"+rs+"\n");
	    
 	while(str[snp_index].equals(rs)){
 		//if(makeIdList)this.ids.add(str[sample_index]);
 		
 		this.output.put(str[sample_index], getPrint(str));
 		currentString = br.readLine();
 		if(currentString==null) break;
 		 str = currentString.split(OptionBuild.splStr);
 	}
 	{
// 		File outF = new File(out, rs);
 		//if(!append && outF.exists()) throw new RuntimeException("input file must be sorted according to rs number "+outF.getName());
 	OutputStreamWriter osw = compress.getWriter(rs, true);
// 		PrintWriter osw = new PrintWriter(new FileWriter(outF,append));
 	  for(int i=0; i<this.ids.size(); i++){
 		  String[] toP = this.output.get(ids.get(i));
 		  if(toP==null) {
 			  System.err.println("was null for "+ids.get(i)+" "+rs);
 			  osw.write("null");
 		  }else{
 		  osw.write(toP[0]);
 		  for(int k=1; k<toP.length; k++){
 			  osw.write("\t");
 			  osw.write(toP[k]);
 		  }
 		  }
 		  osw.write("\n");
 		
 	  }
 	  compress.closeWriter(osw);
 	}
	
 }
 CompressIlluminaTableLong(File f, File output, String build, boolean append) throws Exception{
    	header =   new String[] {
    	            "Genotype\tB allele\tLog R\tTop allelles",
    	            "chr\t"+build+"_start\t"+build+"_end\t id",
    	            "id"};
    	this.append = append;
   // 	this.chromToDo= todo;
    	this.dir = output;
        dir.mkdir();
        
        this.f =f;
     //   this.currentChrom = todo;
        br =Utils.getBufferedReader(f);
    
      //  this.currentSpl = new String[br.length][];
      
         currentString = br.readLine();
         System.err.println(currentString);
       // for(int k=0; k<br.length; k++){
        List<String> order = new ArrayList<String>();
        StringBuffer header1 = new StringBuffer();
        	inputHeader = Arrays.asList(currentString.split(OptionBuild.splStr));
        	System.err.println(Arrays.asList(inputHeader));
        		chrom_index = inputHeader.indexOf("Chr");
        		snp_index = inputHeader.indexOf("SNP Name");
        		loc_index = inputHeader.indexOf("Position");
        		this.sample_index = inputHeader.indexOf("Sample ID");
        		if(sample_index<0) this.sample_index = inputHeader.indexOf("Sample Id");
        		//Assume everything after pos is to print
        		int topminindex = Math.max(Math.max(sample_index,chrom_index), Math.max(loc_index, snp_index));
        		this.toprint_ind = new int[inputHeader.size()-(topminindex+1)]; 
        		for(int i=0; i<toprint_ind.length; i++){
        			toprint_ind[i] = topminindex+i+1;
        			order.add(map(inputHeader.get(toprint_ind[i])));
        			if(i>0) header1.append("\t");
        			header1.append(map(inputHeader.get(toprint_ind[i])));
        		}
       // }
       
        //	 order = new ArrayList<String>(); //has order of b_freq log r, etc
        
        	ids = new ArrayList<String>();
       	 currentString = br.readLine();
      
       
       
     
       header[0] = header1.toString();
       
      
        
    }
    
    
    
   



	
  
  
 static Map<String, String> map = new HashMap<String, String>();
 

 static{
	 String from1 =  
		 ",       ,A      A,A     A,X     X,X     X,Y     A,Y     X,Z     A,B     A,Z     Y,Z     " +
		 ",B      B,B     B,Z     Z,Z     ,X      ,Y      ,Z      B,X     B,Y     Y,Y";
	 String[] from = from1.split("\\s+");
	 for(int i=0; i<from.length; i++){
		 if(i==0) map.put("-",from[i]);
		 else{
			 
			 if(from[i].charAt(0)==',' && from[i].length()>=2){ 
				 char char2 = from[i].charAt(1);
				 if(char2=='Y' || char2=='Z' || char2=='X') continue;
			 }
					 
			char[] st1 = from[i].replaceAll(",", "").replaceAll("X", "AA").replaceAll("Y", "AB").replaceAll("Z", "BB").toCharArray();
			Arrays.sort(st1);
			map.put(new String(st1), from[i]);
		 }
	 }
 }
 
 private String map(String string) {
	if(map.containsKey(string)) {
		String str =  map.get(string);
		return str;
	}
	else return string;
}

private String[] getPrint(String[] str) {
		String[] res = new String[this.toprint_ind.length];
		for(int i=0; i<res.length; i++){
			res[i] = str[toprint_ind[i]];
		}
		return res;
	}


public void getIds() throws Exception{
	String[] str = currentString.split(OptionBuild.splStr);
	String chrom = chrom_index<0 ? "null" : str[this.chrom_index];
  //  File out = writeHeader(chrom);
	File samples = new File(this.dir,"Samples");
	
	if(samples.exists() && samples.length()>0){
	   BufferedReader br = new BufferedReader(new FileReader(samples));
	   String st = "";
	   while((st = br.readLine())!=null){
		   this.ids.add(st);
	   }
	}
	else{//if(makeIdList)this.ids.add(str[sample_index]);
	
    Set<String> ids1 = new HashSet<String>();
    while(this.currentString!=null){
    	str = currentString.split(OptionBuild.splStr);
    	ids1.add(str[this.sample_index]);
    	currentString = br.readLine();
    }
    this.br.close();
    ids.addAll(ids1);
   
   	         // }
   	          br =Utils.getBufferedReader(f);
   	       
   	       //  this.currentSpl = new String[br.length][];
   	       
   	          currentString = br.readLine();
   	       currentString = br.readLine();
	}    
	
	 
}
public void run() throws Exception{
	this.getIds();
	   	String str = "";
	 	 this.parseNext();
	    outer: while(currentString!=null){
		
	    this.parseNext();
	   
	     
   		
}
   this.finish(); 
 //  this.sortSNPs();
}

/*public void sortSNPs() throws Exception{
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
				res[i1] = i1==1 ? Integer.parseInt(str[i1]) : str[i1];
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
		snpsinF.delete();
		snps_out.close();
	}
}*/
final static Comparator<Comparable[]> compa = new Comparator<Comparable[]>(){

	public int compare(Comparable[] o1, Comparable[] o2) {
		return o1[1].compareTo(o2[1]);
	}
	
};
	


}