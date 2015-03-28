import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import conversion.Utils;


public class MergeCoords {
	
	public static void main(String[] args){
		try{
			MergeCoords mc1 = new MergeCoords(new File(args[0]),new File(args[1]));
			mc1.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	PrintWriter out, out1;
	BufferedReader b1; //informative
	BufferedReader b2; //primary

	String spl1=",";
	String spl2="\t";
	
	int ind1=0;
	int ind2=0;
	
	//assumes f1 f2 sorted on ind1 ind2
	public MergeCoords(File f1, File f2) throws Exception{
		b1 = Utils.getBufferedReader(f1);
		b2 = Utils.getBufferedReader(f2);
		out = new PrintWriter(new BufferedWriter(new FileWriter("out")));
		out1 = new PrintWriter(new BufferedWriter(new FileWriter("snps.txt")));
	}
	
	public void run() throws Exception{
		String st1=b1.readLine();
		String st2="";
		String[] str2=null;
		String[] str1=st1.split(spl1);
		String id1=str1[ind1];
		String id2 = "";
		
		while((st2=b2.readLine())!=null){
		    str2 = st2.split(spl2);
			id2 = str2[ind2];
			int compa = id1.compareTo(id2);
			while(st1!=null && compa<0){
				st1=b1.readLine();
				str1 = st1.split(spl1);
				id1=str1[ind1];
				compa =  id1.compareTo(id2);
			}
			
		
			if(compa==0){
				
				out.println("chr"+str1[2]+"\t"+str1[3]+"\t"+(Integer.parseInt(str1[3])+20)+
						"\t"+st2);
				out1.println(str2[0]+"\t"+str1[2]+"\t"+str1[3]);
				
			}
			else{
				out.println("chrNA"+"\tNA"+"\tNA"+
						"\t"+st2);
				out1.println(str2[0]+"\tNA"+"\tNA"
						);
			}
		}
		b1.close();
		b2.close();
		out.close();
		out1.close();
		
	}
}
