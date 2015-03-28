package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ConvertSeq {
	
	 File all;
	public static void main(String[] args){
		try{
			ConvertSeq cs = new ConvertSeq(new File("."));
			cs.run();
			CompressDir1.compress(cs.dirout);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	List<String[]> snps = new ArrayList<String[]>();
	List<String> samples = new ArrayList<String>();
	
	BufferedReader data;
	
	File dirout;
	
	ConvertSeq(File dir) throws Exception{
		{
		BufferedReader snp = new BufferedReader(new FileReader(new File(dir, "snps.csv")));
		snp.readLine();
		String st = "";
		dirout = new File(dir, "out");
		Utils.delete(dirout);
		dirout.mkdir();
		while((st = snp.readLine())!=null){
			String[] str = st.split("\\t");
			String id = str[0];
			String chr = str[1];
			String start = str[2];
			String end = str[3];
			snps.add(new String[]{"chr"+chr, start, end,id});
		}
		snp.close();
		}
		{
			BufferedReader samplesr = new BufferedReader(new FileReader(new File(dir, "Samples")));
			String st = "";
			while((st = samplesr.readLine())!=null){
				String[] str = st.split("\\t");
				samples.add(str[0]);
			}
			samplesr.close();
		}
		data = new BufferedReader(new FileReader(new File(dir, "data.csv")));
		 all = new File(dirout, "all");
		mknewDir(all);
		//data.readLine();
	}
	
	public void mknewDir(File dir) throws Exception{
		if(!dir.exists()){
			dir.mkdir();
			PrintWriter sampleo =new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Samples"))));
			for(int i=0; i<samples.size(); i++){
				sampleo.println(samples.get(i));
			}
			sampleo.close();
			
			PrintWriter nameo =new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Name"))));
			nameo.println("distr");//Genotype
			nameo.println("chr\tstart\tend\tid");
			nameo.println("id");
			nameo.close();
		}
	}
	public void run() throws Exception{
		String st = "";
		for(int i=0; (st= data.readLine())!=null; i++){
			processLine(st,i);
		}
		data.close();
	}
	
	
	
	public void processLine(String line, int i) throws Exception{
		String[] str = snps.get(i);
		String chr = str[0];
		String[] str1 = line.split("\\s+");
		List<String> ids = Arrays.asList(str1).subList(1, str1.length);
		File outChr = new File(dirout, chr.substring(3));
		mknewDir(outChr);
		PrintWriter snps = new PrintWriter(new BufferedWriter(new FileWriter(new File(outChr, "SNPS"),true)));
		snps.print(str[0]);
		PrintWriter snpsall = new PrintWriter(new BufferedWriter(new FileWriter(new File(all, "SNPS"),true)));
		snpsall.print(str[0]);
		for(int k=1; k<str.length; k++){
			snps.print("\t"+str[k].replaceAll(",", ""));
			snpsall.print("\t"+str[k].replaceAll(",", ""));
		}
		snps.println();
		snps.close();
		snpsall.println();
		snpsall.close();
		PrintWriter data = new PrintWriter(new BufferedWriter(new FileWriter(new File(outChr, str[3]),true)));
		PrintWriter allout = new PrintWriter(new BufferedWriter(new FileWriter(new File(all, str[3]),true)));
		for(int k=0; k<samples.size(); k++){
			String samplk = samples.get(k);
			if(ids.contains(samplk)){
				data.println(",=0.5;,A=0.5");
				allout.println(",=0.5;,A=0.5");
			}
			else{
				data.println("AA");
				allout.println("AA");
			}
		}
		data.close();
		allout.close();
	}
	
	
	
}
