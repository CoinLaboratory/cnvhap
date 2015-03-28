package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipOutputStream;

public class CompressedAgilent185CNV extends CompressCNVResults {
	private String currentId;
	public CompressedAgilent185CNV(File f, String chr,  File rawDir) throws Exception {
		super();
		boolean header = true;
		 this.dir = f.getParentFile();
	        this.f =f;
	        prefix = "chr"+chr;
	        
	        BufferedReader br = new BufferedReader(new FileReader(f));
	        String st = "";
	        for(int i=0; i<4; i++){
	        	br.readLine();
	        }
	        if(header){
		    	   String[] str = br.readLine().split("\t");
		    	 System.err.println(str);
		       }
	     //   String start = "chr"+chr;
	        while((st = br.readLine())!=null){
	        	String[] str = st.split("\t");
	        	if(str.length==1 || str[1].length()<=1){
	        		currentId = str[0].split("\\+")[0];
	        	}
	        	else if(getChrom(str).equals(prefix)){
	        		
	        		
	        		
	        		Elemen ele = new Elemen(getStart(str),
	        					getEnd(str), 
	        					getCN(str), 
	        					currentId);
	        		byEnd.add(ele);
	        		byStart.add(ele);
	        	}
	        	
	        }
	       
	       // chr = prefix.split("_")[1].substring(3);#
	      //  this.rawData = new ZipFile(new File(rawDir, chr+".zip"));;
	   /*     BufferedReader br1 = this.getBR("Name");
	        String[] str = br1.readLine().split("\t");
	        br1.close();
	     /*   geno_index = -1;
	        for(int i=0; i<str.length; i++){
	        	if(str[i].toLowerCase().indexOf("geno")>=0){
	        		geno_index = i;
	        		break;
	        	}
	        }*/
	        dest = Compressor.getOS(new File(dir, chr+".zip"));
	        checksum = new   CheckedOutputStream(dest, new Adler32());
	        outS = new 
	        ZipOutputStream(new 
	          BufferedOutputStream(checksum));
	        osw = new OutputStreamWriter(outS);
	        outS.setMethod(ZipOutputStream.DEFLATED);
	        writeSnps( null, 
	        	getBuildBR(),"Samples", samplesM);
	        writeSnps(prefix,getBuildBR(),"SNPS", null);
	}
	

	@Override
	public int getCN(String[] str) {
		double d1 = Double.parseDouble(str[6]);
		double d2 = Double.parseDouble(str[7]);
		double d = Math.abs(d1)>Math.abs(d2) ? d1  : d2;
		if(d<0.3 && d>-0.3) return 2;
		 if(d<-1.0) return 0;
		else if(d<0) return 1;
		else if (d>0.6) return 4;
		else return 3;
	}

	@Override
	public int getEnd(String[] st1) {
		return Integer.parseInt(st1[4]);
	}

	@Override
	public String getId(String[] str) {
		return currentId;
	}

	@Override
	public int getStart(String[] str) {
		return Integer.parseInt(str[3]);
	}
	@Override
	public String getChrom(String[] str){
		//if(str.length<=1) return "";
		return str[1];
	}

}
