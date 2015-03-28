package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public  class CompressDepth {

	public static void main(String[] args) {
	    try{
	    //	if(true) return;
	    	
	    	  File dir1 = new File(System.getProperties().getProperty("user.dir"));
	  
	   // 	File snps = new File("SNPS");
	      String[] chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22"
	    	 // "X_M"
	    	  .split(":");
	      for(int i=0; i<chr.length; i++){
	    	  try{
	    	  System.err.println("doing "+chr[i]);
	    	CompressDepth chmp = new CompressDepth(new File(dir1, "cnv_yh"), chr[i]);
	          
	            
	            	
	            		/*		: new CompressPennCNV(new File("Pennraw.txt"), chr[i], raw)) : 
	            			(a185 ? new CompressedAgilent185CNV(new File("185k_French.txt"), chr[i], new File("../185k")):
	            			new CompressedAgilentCNV(new File("CombinedProbeIntervals.txt"), chr[i], raw));*/
	            chmp.run();
	    	  }catch(Exception exc){
	    		  exc.printStackTrace();
	    	  }
	      }
	    }catch(Exception exc){
	        exc.printStackTrace();
	    }
	}

	ZipFile rawData;
	
	protected File dir;
	protected String prefix;
	protected FileOutputStream dest;
	protected CheckedOutputStream checksum;
	protected ZipOutputStream outS;
	protected OutputStreamWriter osw;
	String[] header = new String[] {
	            "Depth",
	            "chr\tstart\tend\tid",
	            "id"};
	protected File f;
	List<String> snps = new ArrayList<String>();
	//Elemen probe = new Elemen(0,0,0,"");

	
	//protected File snps;

	public void writeHeader() throws Exception {
	    ZipEntry headings = new ZipEntry("Name");
	    outS.putNextEntry(headings);
	        for(int i=0; i<header.length; i++){
	            osw.write(header[i]);
	            osw.write("\n");
	        }
	        osw.write(String.format("%5.3g",this.meanD).trim()+"\n");
	        osw.flush();
	       outS.closeEntry();
	}

	public void writeSnps(List<String> l, String name) throws Exception {
		
	    ZipEntry headings = new ZipEntry(name);
	    outS.putNextEntry(headings);
	  //  String st = "";;
	   for(int i=0; i<l.size();i++ ){
		  
	            osw.write(l.get(i));
	            osw.write("\n");
	        }
	   
	        osw.flush();
	       outS.closeEntry();
	     //  br.close();
	}

	protected List<String> samplesM = new ArrayList<String>();
	BufferedReader br;


	

	double meanD;
	public CompressDepth(File f, String chr) throws Exception {
		 this.dir = f.getParentFile();
	        this.f =f;
	        prefix = ">chr"+chr.split("_")[0];
	        br = new BufferedReader(new FileReader(f));
	        String st = "";
	       this.chrom  = chr.split("_")[0];
	        while((st = br.readLine())!=null){
	        	if(st.equals(prefix)) {
	        		break;
	        	}
	        }
	      
	      meanD= Double.parseDouble(br.readLine().split(":")[1].trim());
	
	        dest = Compressor.getOS(new File(dir, chr+".zip"));
	        checksum = new   CheckedOutputStream(dest, new Adler32());
	        outS = new 
	        ZipOutputStream(new 
	          BufferedOutputStream(checksum));
	        osw = new OutputStreamWriter(outS);
	        outS.setMethod(ZipOutputStream.DEFLATED);
	        samplesM = Arrays.asList(new String[] {"4"});
	        writeSnps( 
	        	samplesM,"Samples");
	    
	}
	

	/*protected BufferedReader getBR(String name) throws Exception{
		return new BufferedReader(new InputStreamReader(rawData.getInputStream(rawData.getEntry(name))));
	}*/
	
 //int geno_index;
	  
	final String chrom;
	public void run() throws Exception {
	   writeHeader();
	 
	   String st = "";;
	   int prev=1;
	   int k1=0;
	   for(int k=0;(st = br.readLine())!=null; k++){
		   if(st.startsWith(">")) break;
		 String[] str = st.split("\\s+");
		 int pos = Integer.parseInt(str[0]);
		 if(pos>prev){
			 Logger.global.info("h");
		 }
		 int pos1 = Integer.parseInt(str[2])+pos;
		 prev = pos1;
		 pos1 = Math.max(pos+1, pos1-5);
		 int[] snps1 = new int[] {pos,pos1};
	
		int steps =(int) Math.floor((double) (pos1-pos)/5000.0);
		
		for(int j=0; j<steps; j++){
			int pos_ = pos + j*5000;
			snps.add("chr"+chrom+"\t"+pos_+"\t"+pos_+"\t"+(k1));
			  ZipEntry headings = new ZipEntry(k1+"");
		
		      outS.putNextEntry(headings);
		      for(int kk=0; kk<this.samplesM.size(); kk++){
		    	  osw.write(str[1]+ (k<samplesM.size()-1 ? "\t" : "\n"));
		      }
			
			    osw.flush();
		      outS.closeEntry();
		      k1++;
		}
	
		 
	
	{
		int pos_ = pos1;
		if(pos_==0) continue;
		snps.add("chr"+chrom+"\t"+pos_+"\t"+pos_+"\t"+(k1));
		  ZipEntry headings = new ZipEntry(k1+"");
	
	      outS.putNextEntry(headings);
	      for(int kk=0; kk<this.samplesM.size(); kk++){
	    	  osw.write(str[1]+ (k<samplesM.size()-1 ? "\t" : "\n"));
	      }
		
		    osw.flush();
	      outS.closeEntry();
	      k1++;
	}
	
		 
		
		  //writeSNP
	  }
	br.close();
	    writeSnps( snps,"SNPS");
	    outS.close();
	}

	private String getString(Map<String, Character> m) {
		StringBuffer sb = new StringBuffer();
		for(Iterator<Entry<String, Character>> it = m.entrySet().iterator(); it.hasNext();){
			Entry<String, Character> ent = it.next();
			sb.append("\t"+ent.getKey()+"="+ent.getValue());
		}
		return sb.toString();
	}

}