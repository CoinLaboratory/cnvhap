package data;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public class PlateAdjustment {

	
	
	Correction correction;
	ZipFile zf;
	//String[][] input;
	
//	int[] alias;
	
	
	
	
	
	
	List<String> indv = new ArrayList<String>();
	int lrr_id;
	int len;
	
	
	
	public PlateAdjustment(File corrections, File zf, String logr, int max) throws Exception{
		
		indv = new ArrayList<String>();
		this.zf = new ZipFile(zf);
		BufferedReader br1 = new BufferedReader(new InputStreamReader(this.zf.getInputStream(this.zf.getEntry("Samples"))));
		String st = "";
		while((st = br1.readLine())!=null){
			indv.add(st);
		}
		
		List<String > l2 =  read1(this.zf, "Name");
		List<String> l20 =Arrays.asList(l2.get(0).split("\t")); 
		fix(l20);
		lrr_id = l20.indexOf(logr);
		if(lrr_id < 0) throw new RuntimeException("did not find "+logr+" in "+l20);
		len = l20.size();
		correction = new Correction(corrections,indv,len,max);
	
	}
	public static List<String> read1(ZipFile zf,  String probe) throws Exception {
		List<String> res = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st);
		}
		br.close();
		return res;
		
		
	}
	
	


	public static void fix(List<String> l21) {
		// TODO Auto-generated method stub
		for(int k=0; k<l21.size(); k++){
			l21.set(k, l21.get(k).trim());
		}
	}
	
	
	
	

	
	
	private String[][] read(ZipFile zf,  String probe) throws Exception {
		ZipEntry ent = zf.getEntry(probe);
		//List<String> res = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
		
	
		
		this.correction.correct(br, this.lrr_id);
		return this.correction.input;
	}
	//String[] input;	
	public String[][] readAdjusted(String snpid) throws Exception {
		return this.read(zf, snpid);
	}
	
}
