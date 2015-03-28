package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.util.CompressDir;

public class CompressCnvPipeOutput  {
	
	public static void main(String[] args){
		try{
		 File f = new File( System.getProperty("user.dir"));
		 CompressCnvPipeOutput ccpo = new  CompressCnvPipeOutput(f);
		 ccpo.compress();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
    CompressCnvPipeOutput(File dir) throws Exception{
    	this.dir = dir;
    	String[] outnmes = dir.getName().split("_");
    	String outnme = outnmes[2]+"_"+outnmes[3]+"_"+outnmes[4];
    	this.outdir =new File(dir,outnme);
    	outdir.mkdir();
    	File samps = (new File(dir, "genos_0_1_3_4.txt"));
    	 br = new BufferedReader(new FileReader(samps));
    	 String st = br.readLine();;
    	 header = Arrays.asList(st.split("\\s+"));
    	 snpid = header.indexOf("regionId");
    	 chrid = header.indexOf("Chr");
    	 startid = header.indexOf("FirstProbe");
    	 endid = header.indexOf("LastProbe");
    	samples = header.subList(snpid+1, header.size());
    }
    int snpid,chrid,startid,endid;
	BufferedReader br;
	File dir,outdir;
	CompressDir compress;
	List<String> header,samples;
	void compress() throws Exception{
		 compress = (new CompressDir(outdir));
		 OutputStreamWriter pw = compress.getWriter("Samples",true);
		 for(int i=0; i<samples.size(); i++){
				pw.write(samples.get(i)+"\n");
			}
			compress.closeWriter(pw);
			
			pw = 	compress.getWriter("Name",true);
			pw.write("CN\n");
			pw.write("chr\tstart\tend\tid\n");
			pw.write("id\n");
			compress.closeWriter(pw);
			OutputStreamWriter snps = compress.getWriter("SNPS", false);
			
			String st = "";
	    	while((st = br.readLine())!=null){
	    		String[] str = st.split("\\s+");
	    		int start = Integer.parseInt(str[startid])-1;
	    		int end = Integer.parseInt(str[endid])+1;
	    		snps.write(str[chrid]+"\t"+start+"\t"+(start+20)+"\t"+str[snpid]+".s"+"\n");
	    		snps.write(str[chrid]+"\t"+end+"\t"+(end+20)+"\t"+str[snpid]+".e"+"\n");
	    		pw = compress.getWriter(str[snpid]+".s",true);
	    		for(int k=snpid+1; k<str.length; k++){
	    			pw.write(str[k]+"\n");
	    		}
	    		compress.closeWriter(pw);
	    		pw = compress.getWriter(str[snpid]+".e",true);
	    		for(int k=snpid+1; k<str.length; k++){
	    			pw.write(str[k]+"\n");
	    		}
	    		compress.closeWriter(pw);
	    	}
	    	br.close();
	    	compress.closeWriter(snps);
	    	compress.run();
	    	compress.close();
	}
}
