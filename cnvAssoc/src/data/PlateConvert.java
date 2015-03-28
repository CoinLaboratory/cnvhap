package data;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

public class PlateConvert {
	public static void main(String args[]){
		try{
			File base = new File(args[0]);
			
			String[] type = args[1].split(":");
			File[] corrections = new File[type.length];
			String chr = args[3];
			String logr_ = args[4];
		
			File[] dir = new File[type.length];
			File[] outp = new File[type.length];
		
			File dirF = new File(System.getProperty("user.dir"));
			
			for(int k=0; k<dir.length; k++){
				File out = new File(type[k]);
				out.mkdir();
				dir[k] = new File(base, type[k]+"/"+chr+".zip");
				File outdir = new File(out,"output");
				outdir.mkdir();
				outp[k] =  new File(outdir,chr+".zip");
				corrections[k] = new File(base,type[k]+"/"+args[2]);
			}
		//	File corrections = new File(args[5])
			int max = args.length>5 ? Integer.parseInt(args[5]) : Integer.MAX_VALUE;
			PlateAdjustment pa = new PlateAdjustment(corrections[0], dir[0],logr_.replace('_', ' '),max) ;
			PlateConvert pc = new PlateConvert(pa, outp[0]);
			pc.run();
			pa.zf.close();
				
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	PlateAdjustment pa;
	FileOutputStream dest;
    CheckedOutputStream checksum;
	ZipOutputStream outS;
	OutputStreamWriter osw;
	
	PlateConvert(PlateAdjustment pa,File output) throws Exception{
		  dest = new FileOutputStream(output);
		  this.pa = pa;
          checksum = new   CheckedOutputStream(dest, new Adler32());
          outS = new 
          ZipOutputStream(new 
            BufferedOutputStream(checksum));
         osw = new OutputStreamWriter(outS);
         outS.setMethod(ZipOutputStream.DEFLATED);
     	this.transfer("Name",pa.zf, new int[][] {null,  new int[] {0,1,2,3},new int[] {0}});
     	this.transfer("Samples", pa.zf,new int[] {0});
        this.transfer("SNPS", pa.zf, new int[] {0,1,2,3});
	}
	List<String> reserved = Arrays.asList("SNPS:Samples:Name".split(":"));
	
	public void run() throws Exception{
		for(Enumeration en = pa.zf.entries();en.hasMoreElements();){
			String name =((ZipEntry)en.nextElement()).getName();
			if(reserved.indexOf(name)<0){
				this.process(name);
			}
		}
		this.osw.close();
		this.outS.close();
		
		this.dest.close();
		this.checksum.close();
		
		pa.zf.close();
	}
	
	private void write(String string, List<String> l2) throws Exception {
		 ZipEntry headings = new ZipEntry(string);
		 outS.putNextEntry(headings);
		 for(int k=0; k<l2.size(); k++){
			 osw.write(l2.get(k));
			 osw.write("\n");
		 }
        osw.flush();
        outS.closeEntry();
        
		
	}
	
	public void process(String snpid) throws Exception{
		String[][]input = pa.readAdjusted(snpid);
		 ZipEntry headings = new ZipEntry(snpid);
         outS.putNextEntry(headings);
		for(int k=0; k<input.length; k++){
	        writeLine1(osw, input[k]);
		}
		   osw.flush();
           outS.closeEntry();
	}
	private void transfer(String string, ZipFile zf) throws Exception {
		 ZipEntry headings = new ZipEntry(string);
		 ZipEntry ent = zf.getEntry(string);
		 outS.putNextEntry(headings);
			//List<String> res = new ArrayList<String>();
			BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
			String st = "";
		for(int k=0;(st =br.readLine() )!=null;k++){
			 osw.write(st);
			 osw.write("\n");
		 }
       osw.flush();
       outS.closeEntry();
	}
	
	final String splStr = "\t";
	private void transfer(String string, ZipFile zf, int[] k1) throws Exception {
		 ZipEntry headings = new ZipEntry(string);
		 ZipEntry ent = zf.getEntry(string);
		 outS.putNextEntry(headings);
			//List<String> res = new ArrayList<String>();
			BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
			String st = "";
		for(int k=0;(st =br.readLine() )!=null;k++){
			String[] str = st.split(splStr);
			 osw.write(str[k1[0]]);
			 for(int j=1; j<k1.length; j++){
				 osw.write("\t"+str[k1[j]]);
			 }
			 osw.write("\n");
		 }
      osw.flush();
      outS.closeEntry();
	}
	
	private void transfer(String string, ZipFile zf, int[][] k_1) throws Exception {
		 ZipEntry headings = new ZipEntry(string);
		 ZipEntry ent = zf.getEntry(string);
		 outS.putNextEntry(headings);
			//List<String> res = new ArrayList<String>();
			BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
			String st = "";
		for(int k=0;(st =br.readLine() )!=null;k++){
			int[] k1 = k_1[k];
			if(k1==null){
				osw.write(st);
			}
			else{
			String[] str = st.split(splStr);
			 osw.write(str[k1[0]]);
			 for(int j=1; j<k1.length; j++){
				 osw.write("\t"+str[k1[j]]);
			 }
			}
			 osw.write("\n");
		 }
     osw.flush();
     outS.closeEntry();
	}
	
	private void writeLine1(OutputStreamWriter osw2, String[] inp2) throws Exception{
		osw2.write(inp2[0]);
		for(int k=1; k<inp2.length; k++){
			osw2.write("\t");
			osw2.write(inp2[k]);
		}
		osw2.write("\n");
		
	}
}
