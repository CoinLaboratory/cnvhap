package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Executors;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class CompressIlluminaTableNoFiles {
//    final String chromToDo;
   // final PropertyChangeSupport pcs;
	final boolean append;
	
	 FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter osw;
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
        	    	CompressIlluminaTableNoFiles chmp = new CompressIlluminaTableNoFiles(dir, new File("output"), "build36", false);
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
    	
           CompressIlluminaTableNoFiles chmp = new CompressIlluminaTableNoFiles(dir, outputDir, build, false);
         chmp.run();
      //   CompressDir.compress(chmp.dir);
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
    
  // final Map<String, PrintWriter> snps = new HashMap<String, PrintWriter>();
    
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
          /*  if(snps.get(chr)==null){
            	 PrintWriter snps = new PrintWriter(new BufferedWriter(new FileWriter(new File(out, "SNPS1"))));
	             this.snps.put(chr, snps);
             }*/
             return out;
         //}
     }
    	 public void writeHeader(File f) throws Exception{
 	    	if(f.isDirectory()){
 	    		File[] f1 = f.listFiles();
 	    		for(int i=0; i<f1.length; i++){
 	    			writeHeader(f1[i]);
 	    		}
 	    	}
 	    	else{
 	        ZipEntry headings = new ZipEntry(f.getAbsolutePath().substring(len));
 	        BufferedReader br = new BufferedReader(new FileReader(f));
 	        outS.putNextEntry(headings);
 	          String str = "";
 	          while((str = br.readLine())!=null){
 	        	  osw.write(str);osw.write("\n");
 	          }
 	            osw.flush();
 	           outS.closeEntry();
 	           br.close();
 	    	}
 	    	
 	           
 	    }
	   
	 //   boolean finished = false;
	  /*  public void finish() throws Exception{
	    	for(Iterator<PrintWriter> it = this.snps.values().iterator(); it.hasNext();){
	    	 it.next().close();
	    	
	    	}
	    }*/
	  
	      
	     
	    
	     
	     
	     
		
    
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
 int snp_index, strand_index;
 int loc_index ; 
 String currentChom="";
 String[] currentString;
 final String[][] stri;
private String overallChrom=null;
 boolean readLine() throws Exception{
	 for(int i=0; i<f.length; i++){
		 currentString[i] = br[i].readLine();
		
	 }
	return (currentString[0]!=null);
 }
 int len;   
 CompressIlluminaTableNoFiles(File[] f, File output, String build, boolean append) throws Exception{
    	header =   new String[] {
    	            "Genotype\tB allele\tLog R\tTop allelles\tstrand",
    	            "chr\t"+build+"_start\t"+build+"_end\t id",
    	            "id"};
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
       // for(int k=0; k<br.length; k++){
        
        	inputHeader = new List[f.length];
        	for(int i=0; i<inputHeader.length; i++){
        		inputHeader[i] = Arrays.asList(currentString[i].split(OptionBuild.splStr));
        	}
        		
        	
        		chrom_index = inputHeader[0].indexOf("Chr");
        		snp_index = inputHeader[0].indexOf("Name");
        		loc_index = inputHeader[0].indexOf("Position");
        		strand_index = inputHeader[0].indexOf("strand");
        		if(snp_index<0) snp_index = inputHeader[0].indexOf("rsID");
        		if(snp_index<0) snp_index = inputHeader[0].indexOf("rs#");
        		if(loc_index< 0) loc_index = inputHeader[0].indexOf("Coord");
        		if(loc_index< 0) loc_index = inputHeader[0].indexOf("position_b36");
        		if(loc_index< 0) loc_index = inputHeader[0].indexOf("pos");
        		if(chrom_index<0){
        			String prefix = f[0].getName();
        			   int indofChr = prefix.indexOf("chr");
        			   int end = Math.min(prefix.indexOf('.', indofChr),prefix.indexOf('_', indofChr));
        			   if(end<0) end =  prefix.length();
        		        this.overallChrom = prefix.substring(indofChr+3,
        		        		end);
        		}
       // }
       
        	 order = new List[f.length];//new ArrayList<String>(); //has order of b_freq log r, etc
        
        	ids = new List[f.length];//new ArrayList<String>();
        	idOffset = new List[f.length];//new ArrayList<Integer>();
        	int i=0;
        	for(;inputHeader[0].get(i).indexOf('.')<0 && (inputHeader[0].get(i).indexOf("NA")<0);i++){
        			this.order_snp.add(inputHeader[0].get(i));
        		
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
		        if(hasAllele) this.header[1] =  "chr\t"+build+"_start\t"+build+"_end\tid\tAlleleA\tAlleleB"+"\tstrand";
		        if(j>0 && !order[j].equals(order[0])) throw new RuntimeException("cols not in same order between input files");
        	}
        
      
       header[0] = order[0].toString().replace(',', '\t');
       header[0] = header[0].substring(1,header[0].length()-1);
       
      
        
    }
    
    
    
   



	
  
  
  
  
    public void split(){
	  for(int i=0; i<currentString.length; i++){
		  stri[i] = currentString[i].split(OptionBuild.splStr);
		  if(i>0 && !stri[i][loc_index].equals(stri[0][loc_index])) {
			  throw new RuntimeException("files must be in same order");
		  }
	  }
  }
   int overall_index=0;
 //  String alleleA, alleleB;
  SortedSet<Character> alleleSet = new TreeSet<Character>();
final boolean[][] alleles;
PrintWriter snps;
File chromDir;
public void newChrom(String chrom) throws Exception{
	if(!this.currentChom.equals("")){
		 this.sortSNPs();
		   this.writeHeader(new File(chromDir, "Samples"));
		   this.writeHeader(new File(chromDir, "SNPS"));
		   this.writeHeader(new File(chromDir, "Name"));
		   this.outS.close();
		   this.dest.close();
		  // Utils.delete(chromDir);
	}
	this.currentChom = chrom;
	if(chrom!=null){
		  chromDir =  writeHeader(chrom);
			if(chromDir!=null && chromDir.exists()) throw new RuntimeException("should not exist "+chromDir );
		   
			 snps_f = new File(chromDir, "SNPS1");
		     snps = new PrintWriter(new BufferedWriter(new FileWriter(snps_f)));
		 
		     dest = Compressor.getOS(new File(dir, chrom+".zip"));
		 	len = chromDir.getAbsolutePath().length()+1;
		     checksum = new   CheckedOutputStream(dest, new Adler32());
		     outS = new 
		     ZipOutputStream(new 
		       BufferedOutputStream(checksum));
		     osw = new OutputStreamWriter(outS);
		     outS.setMethod(ZipOutputStream.DEFLATED);
		
	}
}

public void run() throws Exception{
	   //	String str = "";
	    outer: while(this.readLine()){
		overall_index++;
		alleleSet.clear();
	  //   String[][] stri = split();
	    	this.split();
	     //System.err.println(Arrays.asList(this.currentString));
	     //System.err.println(chrom_index);
	     String chrom = chrom_index>=0 ? stri[0][this.chrom_index] : this.overallChrom;
	     if( !chrom.equals(this.currentChom)){
	    	 this.newChrom(chrom);
	     }
	  
	       String snpid = snp_index>=0  ? stri[0][snp_index] : overall_index+"";
	     ZipEntry entry = new ZipEntry(snpid);
	     outS.putNextEntry(entry);
	  //   File outF = new File(out, snpid);
	 	//if(!append && outF.exists()) throw new RuntimeException("input file must be sorted according to rs number");
	  //   PrintWriter osw = new PrintWriter(new FileWriter(outF,append));
	     int st = Integer.parseInt(stri[0][loc_index]);
        String strand = strand_index>=0 ? stri[0][strand_index] : "null";
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
	    				  osw.write(stri[j][offset+kk]);
	    				 if( kk<order[j].size()-1 ){osw.write("\t");}
	    				 else osw.write("\n");
	    				//  osw.print(? "\t":"\n");
	    		
	    			  }
	     }
		 
      }
        //System.err.println(snpid+" "+alleleSet);
        if(alleleSet.size()>0){
        	
        	 snps.println("chr"+chrom+"\t"+st+"\t"+(st+20)+"\t"+snpid+"\t"+alleleSet.first()+"\t"+alleleSet.last()+"\t"+strand);
        }
        else{
              snps.println("chr"+chrom+"\t"+st+"\t"+(st+20)+"\t"+snpid+"\t"+strand);
        }
        osw.flush();
        outS.closeEntry();
    
//	     osw.close();
	     
   		
}
newChrom(null);
  // this.finish(); 
  
}

File snps_f;
   
   public void sortSNPs() throws Exception{
	  this.snps.close();
		File snpsinF = this.snps_f;
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
		PrintWriter snps_out = new PrintWriter(new BufferedWriter(new FileWriter(new File(snps_f.getParentFile(), "SNPS"))));
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
final static Comparator<Comparable[]> compa = new Comparator<Comparable[]>(){

	public int compare(Comparable[] o1, Comparable[] o2) {
		return o1[1].compareTo(o2[1]);
	}
	
};
	


}