package cnvtrans;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import data.util.CompressDir;

public class CompressGenos {
	public static void main(String[] args){
		//if(true) System.exit(0);
		try{
			//genoFile, inputfiles, outputfile
			
			String genoF = args[0];
			String[] files = args[1].split(":");
			File dirF = new File(System.getProperty("user.dir"));
			File[] f = new File[files.length];
			for(int i=0; i<f.length; i++){
				f[i] = new File(dirF, files[i]+"/"+genoF );
			}
		  	CompressGenos cg = new CompressGenos(f, new File(dirF,args[2]));
		  	cg.run();	
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	/*public static void main(File f){
		try{
		  	File dirF = new File(System.getProperty("user.dir"));
		  	CompressGenos cg = new CompressGenos(f);
		  	cg.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}*/
	
	
	
	
 CompressDir f;
 BufferedReader[] br;
 List<String> header;
 List<String> samples;
 int st_id, end_id, rs_id;
 String[][] current;
 
 public boolean readLine() throws Exception{
	 boolean hasMore = true;
	 for(int i=0; i<br.length; i++){
		 String st =   br[i].readLine();
		 hasMore = st!=null;
		 if(hasMore){
			 current[i] = st.split("\t");
		 }
	 }
	 return hasMore;
 }
 CompressGenos(File[] f1, File outp) throws Exception{
	br = new BufferedReader[f1.length];
	current = new String[f1.length][];
	samples = new ArrayList<String>();
	for(int i=0; i<br.length; i++){ br[i] = 
		new BufferedReader(new FileReader(f1[i])); 
		header = Arrays.asList(br[0].readLine().split("\t"));
		samples.addAll(header.subList(9, header.size()));	
	}
	
	
	st_id = 2;
	end_id = 3;
	rs_id = 8;
	f = new CompressDir(new File(outp, f1[0].getName().split("\\.")[0]));
 }

 
 public void run() throws Exception {
	OutputStreamWriter osw =  f.getWriter("Samples", true);
	for(int i=0; i<samples.size(); i++){
		osw.write(samples.get(i)+"\n");
	}
	f.closeWriter(osw);
	osw =  f.getWriter("Name", true);
	osw.write("countAll\n");
	osw.write("chr\tstart\tend\tsnpid\n");
	osw.write("sample\n");

	f.closeWriter(osw);
	osw = f.getWriter("SNPS", false);
	while(readLine()){
	String[] str = current[0];
		String snpid = str[1]+"_"+str[st_id];
		osw.write("chr"+str[1]+"\t"+str[st_id]+"\t"+str[end_id]+"\t"+snpid+"\n");
		OutputStreamWriter osw1 = null;
		for(int kk=0; osw1==null; kk++){
			try{
			osw1 = f.getWriter(snpid+(kk==0 ? "" : ""+kk), true);
			}catch(Exception exc){
				osw1=null;
			}
		}
		for(int ik=0; ik<current.length; ik++){
			str = current[ik];
		for(int i= 9; i<str.length; i++){
			String genot = str[i];
			if(genot.equals(".")) genot="2";
			osw1.write(genot+"\n");
		}
		}
		f.closeWriter(osw1);
		
	}
	f.closeWriter(osw);
	f.run();
	 
 }
 
 
 
 
}
