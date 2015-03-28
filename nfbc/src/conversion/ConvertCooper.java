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

public class ConvertCooper {
	
	 File all;
	public static void main(String[] args){
		try{
			ConvertCooper cs = new ConvertCooper(new File("."));
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
	
	ConvertCooper(File dir) throws Exception{
		{
		BufferedReader snp = new BufferedReader(new FileReader(new File(dir, "snps.csv")));
		//snp.readLine();
		String st = "";
		dirout = new File(dir, "out");
		Utils.delete(dirout);
		dirout.mkdir();
		PrintWriter build = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "build.txt"))));
		for(int i=0;(st = snp.readLine())!=null; i++){
			String[] str = st.split("\\t");
			String id = "id"+i;
			String chr = str[0];
			String start = str[1];
			String end = str[2];
			snps.add(new String[]{"chr"+chr, start, end,id});
			if(rep<=1){
			build.println("chr"+st+"\t"+id);
			}
			else{
				double sta_ = Double.parseDouble(start);
				double  end_ = Double.parseDouble(end);
				for(int k=0; k<rep; k++){
					int stk = (int) (sta_ + ((end_-sta_)/(double)rep) * (double)k);
					int endk = stk+20;
					build.println("chr"+str[0]+"\t"+stk+"\t"+endk+"\t"+id+"_"+k);
				}
			}
		}
		snp.close();
		build.close();
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
		data = new BufferedReader(new FileReader(new File(dir, "ref.csv")));
		 all = new File(dirout, "all");
		mknewDir(all);
		//data.readLine();
	}
	static int rep = 101;
	
	public void mknewDir(File dir) throws Exception{
		if(!dir.exists()){
			dir.mkdir();
			PrintWriter sampleo =new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Samples"))));
			for(int i=0; i<samples.size(); i++){
				sampleo.println(samples.get(i));
			}
			sampleo.close();
			
			PrintWriter nameo =new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Name"))));
			nameo.println("Genotype");//Genotype
			nameo.println("chr\tstart\tend\tid");
			nameo.println("id");
			nameo.close();
		}
	}
	public void run() throws Exception{
		String st = "";
		for(int i=0; (st= data.readLine())!=null; i++){
			if(rep<=1){
			processLine(st,i,"");
			}
			else{
				for(int k=0; k<rep; k++){
					processLine(st,i,"_"+k);
				}
			}
		}
		data.close();
	}
	
	
	
	public void processLine(String line, int i, String suff) throws Exception{
		String[] str = snps.get(i);
		String str_k = str[3]+suff;
		String chr = str[0];
		String[] str1 = line.split("\\s+");
		List<String> ids = Arrays.asList(str1);//.subList(1, str1.length);
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
		PrintWriter data = new PrintWriter(new BufferedWriter(new FileWriter(new File(outChr, str_k),true)));
		PrintWriter allout = new PrintWriter(new BufferedWriter(new FileWriter(new File(all, str_k),true)));
		for(int k=0; k<samples.size(); k++){
			int v = Integer.parseInt(ids.get(k));
			String samplk = samples.get(k);
			if(v<0){
				data.println("NC");//(",=0.5;,A=0.5");
				allout.println("NC");//,=0.5;,A=0.5");
			}
			else{
				String st = v==0 ? "__" : 
					v==1 ? "A_" : "AA"; 
				data.println(st);
				allout.println(st);
			}
		}
		data.close();
		allout.close();
	}
	
	
	
}
