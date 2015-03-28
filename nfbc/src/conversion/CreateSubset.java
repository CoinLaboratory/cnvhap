package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public  class CreateSubset {

	public static void main(String[] args) {
	    try{
	    	  File dir1 = new File(System.getProperties().getProperty("user.dir"));
	    	boolean penn = dir1.getName().toLowerCase().indexOf("penn")>=0;
	    	boolean a185 = dir1.getName().toLowerCase().indexOf("185")>=0;
	    	boolean ds = dir1.getName().toLowerCase().indexOf("ds")>=0;
	    	File samples = new File("Samples");
	    	File snps = new File("SNPS");
	      String[] chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22".split(":");
	      for(int i=chr.length-1; i>=0; i--){
	    	  System.err.println("doing "+chr[i]);
	    	
	            CreateSubset chmp = new CreateSubset(new File("1M"), new File("610k"), chr[i]);
	            chmp.copy();
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
	            "Genotype",
	            "chr\tstart\tend\tid",
	            "id"};
	protected File f;
	
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

	protected List<String> samplesM = new ArrayList<String>();
	BufferedReader br;


	

	public CreateSubset(File dirIn, File dirOut, String chrom)  throws Exception{
		this.in = dirIn;
		this.out = dirOut;
		this.chr = chrom;
		File f = new File(dirOut, "610.txt");
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dirOut, chrom+"_miss.txt"))));
		BufferedReader br = new BufferedReader(new FileReader(f));
		String st = "";
		//Map<String, Set<Integer>> m = new HashMap<String, Set<Integer>>();
		SortedMap<Integer, String> s = new TreeMap<Integer, String>();
		BufferedReader br1 = new BufferedReader(new FileReader(new File(dirIn, "build36.txt")));
		while((st = br1.readLine())!=null){
			String[]str = st.split("\\s+");
			String chr = str[0].substring(3);
			if(!chr.equals(chrom)) continue;
			Integer pos = Integer.parseInt(str[1]);
			s.put(pos, str[3]);
		/*	if(s.contains(pos)){
				subset.add(str[3]);
				s.remove(pos);
			}*/
			/*else{
				System.err.println("did not contain "+str[3]+" "+s.headSet(pos));
			}*/
			
		}
	
		br1.close();
		
		while((st = br.readLine())!=null){
			String[]str = st.split("\\s+");
			String chr = str[1];
			if(!chr.equals(chrom)) continue;
			Integer pos = Integer.parseInt(str[2]);
			if(s.containsKey(pos)){
				subset.add(s.get(pos));
			}
			else{
				SortedMap<Integer, String> tail = s.tailMap(pos);
				SortedMap <Integer, String> head = s.headMap(pos);
				pw.print("did not contain "+chr+" "+pos);
				if(tail.size()>0){
					pw.print("\t tail "+(tail.firstKey() - pos));
				}
				if(head.size()>0){
					pw.print("\t head "+(head.lastKey()- pos));
				}
				pw.println();
			}
			
		}
		br.close();
		pw.close();
		
	}

	
	
	final File in;
	final File out;
	final String chr;
	final Set<String> subset = new HashSet<String>();
	public  void copy() throws Exception{
	    //List<String> excs = Arrays.asList(exception);
	    File f1 = new File(in,chr+".zip" );
	    ZipFile zf = new ZipFile(f1);
	   
	//    String chr = f1.getName().substring(0,f1.getName().indexOf(".zip"));
	   // File tmp = new File(f1.getParentFile(), chr+System.currentTimeMillis()+".zip");
	    
	    FileOutputStream dest = Compressor.getOS(new File(out, chr+".zip"));
	    checksum = new   CheckedOutputStream(dest, new Adler32());
	    outS = new 
	     ZipOutputStream(new 
	       BufferedOutputStream(checksum));
	   osw = new OutputStreamWriter(outS);
	   outS.setMethod(ZipOutputStream.DEFLATED);
	    for(Enumeration en = zf.entries(); en.hasMoreElements();){
	        ZipEntry ent = (ZipEntry) en.nextElement();
	       String name = ent.getName();
	       if(name.startsWith("Name") || name.startsWith("Sample") || subset.contains(name)){
	        List<String> res = Compressor.getIndiv(zf, ent.getName());
	           ZipEntry headings = new ZipEntry(ent.getName());
	           outS.putNextEntry(headings);
	           
	               for(int i=0; i<res.size(); i++){
	                   osw.write(res.get(i));
	                   osw.write("\n");
	               }
	          
	           osw.flush();
	           outS.closeEntry();
	       }
	       if(name.startsWith("SNPS")){
	    	   List<String> res = Compressor.getIndiv(zf, ent.getName());
	           ZipEntry headings = new ZipEntry(ent.getName());
	           outS.putNextEntry(headings);
	           
	               for(int i=0; i<res.size(); i++){
	            	   String str = res.get(i);
	            	   String snpid = str.split("\\s+")[3];
	                   if(subset.contains(snpid)){
	                	   osw.write(str);
	                	   osw.write("\n");
	                   }
	                  
	               }
	          
	           osw.flush();
	           outS.closeEntry();
	       }
	    }
	    outS.close();
	    zf.close();
	  /* boolean deleted = f1.delete();
	   if(deleted) tmp.renameTo(f1);
	   else{
	       throw new RuntimeException("could not move");
	   }*/
	}
	

}