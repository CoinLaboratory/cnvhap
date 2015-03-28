package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;


public abstract class CompressCNVResults {
File rawDir;
File buildFile;
static String build = "build35.txt";	
public static void main(String[] args) {
	    try{
	    	//if(true) return;
	    	String clazz = "conversion.CompressCooperSeq";//args[0];f
	    	String inpFile = args[0];
	    	build = args[1];
	    	Class claz = Class.forName(clazz);
	    	File raw = new File("/home/lcoin/Data/data1/1M");
	    	  File dir1 = new File(System.getProperties().getProperty("user.dir"));
	    //	boolean penn = dir1.getName().toLowerCase().indexOf("penn")>=0;
	    //	boolean a185 = dir1.getName().toLowerCase().indexOf("185")>=0;
	    //	boolean ds = dir1.getName().toLowerCase().indexOf("ds")>=0;
	    	File samples = new File("Samples");
	    	File snps = new File("SNPS");
	      String[] chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22"
	    	 // "X_M"
	    	  .split(":");
	      for(int i=0; i<chr.length; i++){
	    	  try{
	    	  System.err.println("doing "+chr[i]);
	    	
	            CompressCNVResults chmp = new CompressCooperSeq(new File(inpFile),chr[i], raw);
//	            	(CompressCNVResults) claz.getConstructor(new Class[] {File.class, String.class, File.class}).newInstance(
	//            		new Object[] {new File(inpFile),chr[i], raw});
	            
	            	
	            		/*		: new CompressPennCNV(new File("Pennraw.txt"), chr[i], raw)) : 
	            			(a185 ? new CompressedAgilent185CNV(new File("185k_French.txt"), chr[i], new File("../185k")):
	            			new CompressedAgilentCNV(new File("CombinedProbeIntervals.txt"), chr[i], raw));*/
	            chmp.initialise();
	            	chmp.run();
	    	  }catch(Exception exc){
	    		  exc.printStackTrace();
	    	  }
	      }
	    }catch(Exception exc){
	        exc.printStackTrace();
	    }
	}

	//ZipFile rawData;
	
	protected File dir;
	protected String prefix;
	protected FileOutputStream dest;
	protected CheckedOutputStream checksum;
	protected ZipOutputStream outS;
	protected OutputStreamWriter osw;
	String[] header = new String[] {
	            "Genotype",
	            "chr\tstart\tend\tid",
	            "id"};
	protected File f;
	Elemen probe = new Elemen(0,0,0,"");

	protected SortedSet<Elemen> byStart = new TreeSet<Elemen>(new Comparator<Elemen>(){
		public int compare(Elemen o1, Elemen o2) {
			if(o1.start!=o2.start) return o1.start < o2.start ?  -1 : 1;
			if(o1.end!=o2.end) return o1.end < o2.end  ? -1 : 1;
			else return o1.id.compareTo(o2.id);
		}
		  
	  });
	protected SortedSet<Elemen> byEnd = new TreeSet<Elemen>(new Comparator<Elemen>(){
			public int compare(Elemen o1, Elemen o2) {
				if(o1.end!=o2.end) return o1.end < o2.end  ? -1 : 1;
				if(o1.start!=o2.start) return o1.start < o2.start ?  -1 : 1;
				else return o1.id.compareTo(o2.id);
			}
			  
		  });
	//protected File snps;

	public void writeHeader() throws Exception {
	    ZipEntry headings = new ZipEntry("Name");
	    outS.putNextEntry(headings);
	        for(int i=0; i<header.length; i++){
	            osw.write(header[i]);
	            osw.write("\n");
	        }
	        osw.flush();
	       outS.closeEntry();
	}

	public void writeSnps(String start, BufferedReader br, String name, List<String> m) throws Exception {
		
	    ZipEntry headings = new ZipEntry(name);
	    outS.putNextEntry(headings);
	    String st = "";;
	   for(; (st = br.readLine())!=null; ){
		   if(start==null || st.startsWith(start)){
		   if(m!=null){
			   m.add(st.split("\t")[0].split("#")[0]);
		   }
	            osw.write(st);
	            osw.write("\n");
	        }
	   }
	        osw.flush();
	       outS.closeEntry();
	       br.close();
	}
public void writeSnps(String start,  String name, List<String> m) throws Exception {
		
	    ZipEntry headings = new ZipEntry(name);
	    outS.putNextEntry(headings);
	   // String st = "";;
	   for(int i=0; i<m.size(); i++ ){
		   String st = m.get(i);
		   if(start==null || st.startsWith(start)){
		  
	            osw.write(st);
	            osw.write("\n");
	        }
	   }
	        osw.flush();
	       outS.closeEntry();
	   //    br.close();
	}

	protected List<String> samplesM = new ArrayList<String>();
	BufferedReader br1;


	
public abstract String getChrom(String[] str);
public void initialise() throws Exception{
	
     String st = "";
   
    
  //   String start = "chr"+chr;
    
     while((st = br1.readLine())!=null){
     	String[] str = st.replaceAll("NA", "NaN").replaceAll("amb", "NaN").split("\t");
     	if(getChrom(str).equals(prefix) ){
     		
     		boolean exclude = exclude(str);
     		if(!exclude){
	     		Elemen ele = new Elemen(getStart(str),
	     					getEnd(str), 
	     					getCN(str), 
	     					getId(str));
	     		byEnd.add(ele);
	     		byStart.add(ele);
     		}
     	}
     	
     }
}

BufferedReader br;
public CompressCNVResults(File f, String chr, boolean header, File rawDir, List samples) throws Exception {
		 this.dir = f.getParentFile();
	        this.f =f;
	      br1 = new BufferedReader(new FileReader(f));
	        prefix = "chr"+chr.split("_")[0];
	        if(header){
		    	   String[] str = br1.readLine().split("\t");
		    	 setHeader(str);
		       }
	       
	        File f1 = new File(rawDir, chr+".zip");
	        if(!f1.exists()) {
	        	throw new RuntimeException("!!");
	        }
	        System.err.println("openging "+f1.getAbsolutePath());
	        this.rawDir = rawDir;
	        this.buildFile = new File(rawDir, build);
	  
	        dest = Compressor.getOS(new File(dir, chr+".zip"));
	        checksum = new   CheckedOutputStream(dest, new Adler32());
	        outS = new 
	        ZipOutputStream(new 
	          BufferedOutputStream(checksum));
	        osw = new OutputStreamWriter(outS);
	        outS.setMethod(ZipOutputStream.DEFLATED);
	        samplesM.addAll(samples);
	        writeSnps( null, 
	        	"Samples", samplesM);
	   
	        writeSnps(prefix,     getBuildBR(),"SNPS", null);
	}
void setHeader(String[] str) {
	// TODO Auto-generated method stub
	
}

	boolean exclude(String[] str) {
		// TODO Auto-generated method stub
		return false;
	}

	public CompressCNVResults() {
		// TODO Auto-generated constructor stub
	}

/*	protected BufferedReader getBR(String name) throws Exception{
		return new BufferedReader(new InputStreamReader(rawData.getInputStream(rawData.getEntry(name))));
	}*/
	 public abstract int getStart(String[]str);
	  public abstract int getEnd(String[]st1);
	  public abstract int getCN(String[] str);
	  public abstract String getId(String[] str);
 //int geno_index;
	public void run() throws Exception {
	   writeHeader();
	   BufferedReader br  = getBuildBR();
	   String st = "";;
	   while((st = br.readLine())!=null){
		   if(!st.startsWith(prefix)) continue;
		 String[] str = st.split("\\s+");
		 int pos = Integer.parseInt(str[1]);
		 int pos1 = Integer.parseInt(str[2]);
		 probe.set(pos);
		 probe.start = pos-40;
		 probe.end = pos+40;
		  ZipEntry headings = new ZipEntry(str[3]);
		//  BufferedReader br1 = getBR(str[3]);
	      outS.putNextEntry(headings);
		 SortedSet<Elemen> start = byStart.headSet(probe);
		 boolean X = str[0].indexOf('X')>=0;
		 //SortedSet<Elemen> end = byEnd.tailSet(probe);
		 Map<String, Elemen> eles = new HashMap<String,Elemen>();
		 for(Iterator<Elemen> it = start.iterator(); it.hasNext();){
			 Elemen next = it.next();
			 if(probe.overlap(next)>10){
				 eles.put(next.id, next);
			 }
		 }
	//	 if(eles.size()>0){
		//	 System.err.println("could detect at "+st+" "+eles);
	//	 }
		// String st1 = "";
		// System.err.println(pos+" "+eles.keySet());
		  if(eles.size()==0){
			  for(int i=0; i<samplesM.size(); i++){
			//	  String genot = st1.split("\t")[geno_index];
				  if(X){
					  osw.write("A\n");
				  }
				  else{
					  osw.write("AA\n");
				  }
			  }
		  }
		  else{
			  for(int i=0; i<samplesM.size(); i++){
				//  String[] str1 = st1.split("\t");
				//  String genot = str1[geno_index];
				  Elemen ele = eles.get(samplesM.get(i));
				  String toP;
				  if(X){
					  if(ele==null || ele.cn==1){
						  toP = "A";
					  
					  }
					  else if(ele.cn==0){
						  toP = "_";
					  }
					  else if(ele.cn==2){
						  toP = "X";
					  }
					  else {
						  throw new RuntimeException("!!");
					  }
					  
				  }
				  else{
					//  if(true) throw new RuntimeException("!!");
				  if(ele==null || ele.cn==2){
					  toP =  "AA" ;
				  }
				  else if(ele.cn==0){
					  toP = "__";
				  }
				  else if(ele.cn==1){
					  toP = "A_";
				  }
				  else if(ele.cn==3) {
					  toP = "AX";
				  }
				  else if(ele.cn==4) 
					  {
					  toP = "XX";
					  }
				  else throw new RuntimeException("");
				  }
			          osw.write(toP+"\n");
			    
			  
				  
			  }
		  }
		  
		    osw.flush();
	      outS.closeEntry();
		  
		 
		
		  //writeSNP
	  }
	
	
	
	    outS.close();
	}

	 BufferedReader getBuildBR() throws Exception {
	return Utils.getBufferedReader(buildFile);
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