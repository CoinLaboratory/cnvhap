package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class CompressPheno {
	static String spl = "\t";
	public static void main(String[] args){
		try{
			//if(true) return;
			boolean factor = Boolean.parseBoolean(args[1]);
			File dir = new File(args[0].split("\\.")[0]);
			if(dir.exists() && false) {
				throw new RuntimeException("existed");
			}
			else{
				File inp = new File(args[0]);
				BufferedReader br = new BufferedReader(new FileReader(inp));
			dir.mkdir();
			String[] head = br.readLine().split(spl);
			Printer[] pw = new Printer[head.length];
			for(int k=0; k<pw.length; k++){
				String nme = head[k].replace('_', '_').replace('/', '_');
				if(nme.equals("PATIENT")){
					nme = "Samples";
				}
			
				pw[k] =new PrintW(dir,nme,factor && !nme.equals("Samples"));
			}
			run(br,pw);
			}
			if(dir.isDirectory()){
				(new CompressDir1(dir)).run();
			}
		}catch(Exception exc){
		exc.printStackTrace();	
		}
	}
	
	static interface Printer{
		public void println(String str);
		public void close();
	}
	static class PrintW implements Printer{
		PrintWriter pw;
		boolean fact = false;
		List<String> factors = new ArrayList<String>();
		File dir; 
		String name;
	//	List<Integer> vals = new ArrayList<Integer>();
		PrintW(File dir, String name, boolean fact) throws Exception{
			this.pw =   new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name)))); 
			this.fact = fact;
			this.dir = dir;
			this.name = name;
		}
		@Override
		public void close() {
		pw.close();
		try{
			if(factors.size()>0){
				pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir,"bins_"+ name)))); 
				for(int i=0; i<factors.size(); i++){
					pw.println(i+"<="+factors.get(i)+"<"+(i+1));
				}
				pw.close();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		}

		@Override
		public void println(String str) {
			String str1 = str;
			if(fact){
				int ind = factors.indexOf(str);
				if(ind <0){
					ind = factors.size();
					factors.add(str);
				}
				//str1 = ind+"";
			}
			
			pw.println(str1);
			
		}
		
	}
	
	public static void run(BufferedReader br, Printer[] pw) throws Exception{
		
		
		
		String st = "";
		
		while((st =br.readLine() )!=null){
			String[] str = st.split(spl);
			for(int k=0; k<pw.length; k++){
				if(k<str.length)
				pw[k].println(str[k]);
				else
					pw[k].println("");
			}
		}
		for(int k=0; k<pw.length; k++){
			pw[k].close();
		}
	}
	
}
