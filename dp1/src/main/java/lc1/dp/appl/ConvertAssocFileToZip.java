package lc1.dp.appl;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import lc1.util.CompressDir;
import lc1.util.Compressor;
import lc1.util.Constants;

public class ConvertAssocFileToZip {
 public static void main(String[] args){
	 try{
		 File f = new File( System.getProperty("user.dir"));
		 File f1 = new File(f,args[0]);
		 if(!f1.exists()) f1.mkdir();
		if((new File(f,"avg")).exists()){
		 main(f, false);
		}else{
			mainImpute(f1,args[1],args[2], Integer.parseInt(args[3]));
			
		}
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
 }
	
 public static void main(BufferedReader br, File dir1, String[] rsid) throws Exception{
		dir1.mkdir();
		String[] loc = br.readLine().split("\t");
		OutputStreamWriter pw = 	compress.getWriter("Name",true);
		String[] tags = Constants.writeAverages;
		pw.write(tags[0]);
		for(int k=1; k<tags.length; k++){
			pw.write("\t"+tags[k]);
		}
		pw.write("\n");
		//pw.write("countAll\tstate.0\tstate.1\tstate.2\n");
		pw.write("chr\tstart\tend\tsnpid\n");
		pw.write("sample\n");
		compress.closeWriter(pw);
		List<Double>[][] l = new List[loc.length-1][4];
	
		 pw = 	compress.getWriter("SNPS",true);
		for(int kk=1; kk<rsid.length; kk++){
			for(int j=0; j<4; j++){
				l[kk-1][j] = new ArrayList<Double>();
			}
			pw.write(chr+"\t"+loc[kk]+"\t"+(Integer.parseInt(loc[kk])+20)+"\t"+rsid[kk]+"\n");
		}
		compress.closeWriter(pw);
		String st = br.readLine();
		List<String> samples = new ArrayList<String>();
		while((st)!=null){
			String[] str = st.split("\t");
			 String indiv = str[0];
			int ind = indiv.lastIndexOf('_');
			
			indiv = indiv.substring(0,ind);
			samples.add(indiv);
			
			inner1: for(int j=0; st!=null && st.startsWith(indiv); j++){
				str = st.split("\t");
				 ind = str[0].lastIndexOf('_');
				 String indiv1 = str[0].substring(0,ind);
				if(!indiv1.equals(indiv)) {
					System.err.println(indiv);
					System.err.println(indiv1);
					break inner1;
				}
				String type = str[0].substring(ind+1);
				//System.err.writeln(type+" "+indiv);
				for(int kk=1; kk<loc.length; kk++){
					if(kk<str.length)
					l[kk-1][j].add(Double.parseDouble(str[kk]));
					else
						l[kk-1][j].add(Double.NaN);
				}
				st = br.readLine();
			}
		}
		pw = compress.getWriter("Samples",true);
		for(int i=0; i<samples.size(); i++){
			pw.write(samples.get(i)+"\n");
		}
		compress.closeWriter(pw);
		for(int kk=0; kk< l.length; kk++){
			pw = 	compress.getWriter(rsid[kk+1],true);
			for(int i=0; i<samples.size(); i++){
			inner: for(int j=0; j<4; j++){
					pw.write(l[kk][j].get(i)+"");
				    if(j<3 && l[kk][j+1].size()>0){
					   pw.write("\t");
				    }
				    else{
				      break inner;	
				    }
			}
			pw.write("\n");
			}
			compress.closeWriter(pw);
		}
		
 }
 
 public static void compressImpute(BufferedReader br, File dir,  String nme,  File karyfile, File samplesFile, int max) throws Exception{
	 
	 ZipFile zipSamples = new ZipFile(samplesFile);
	 List<Integer> kary = new ArrayList<Integer>();
	 List<String> append = new ArrayList<String>();
		 getKary(karyfile,nme.substring(3), kary, append);
	int[] top = new int[3];
	 String st = br.readLine();
	outer: for(int ind =0; st!=null; ind++){
		int start =0;
	 outer1: for(int app1=0; st!=null && start < kary.get(ind); app1++){
	 // app1 = (int)Math.floor((double)kk/ (double)max);
	  File dir1 = new File(dir,nme+append.get(ind)+( app1));
	  int start_next = kary.get(ind);
	  compress = (new CompressDir(dir1));
	
	  OutputStreamWriter snps = compress.getWriter("SNPS", false);
		 String chr = dir1.getName();
		inner: for(int kk1=0; (st )!=null ; kk1++){
			String[] str = st.split(" ");
			//System.err.println(str[2]);
			 start = Integer.parseInt(str[2]);
			if(start > start_next || kk1>=max) break inner;
			//snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"\t"+str[0]+"\n");
			OutputStreamWriter pw = null;
			try{
			  pw= compress.getWriter(str[1], true);
			  snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"\t"+str[0]+"\n");
			}catch(ZipException exc){
				exc.printStackTrace();
				for(int i=1; pw==null; i++){
				 try{
				 pw= compress.getWriter(str[1]+"_"+i, true);
				 snps.write(chr+"\t"+start+"\t"+(start+20)+"\t"+str[1]+"_"+i+"\t"+str[0]+"\n");
				 }catch(ZipException exc1){
					 exc1.printStackTrace();
				 }
				}
			}
			for(int k=5;k<str.length; k+=3){
				StringBuffer sb = new StringBuffer();
				int sum=0;
				int max_ind = 0;
				for(int j=0; j<3; j++){
					top[j] = (int) Math.round(Double.parseDouble(str[k+j])*1000.0);
					sum+=top[j];
					if(top[j]>top[max_ind]) max_ind = j;
				}
				if(sum!=1000) top[max_ind] = top[max_ind] + (1000-sum);
				for(int j=0; j<3; j++){
					if(j==max_ind){
						 sb.append("-");
					}else if(top[j]>0){
					   sb.append(top[j]);
					}
					if(j<2) sb.append(",");
				}
				sb.append("\n");
				pw.write(sb.toString());
			}
			compress.closeWriter(pw);
			st=br.readLine();
		}
		snps.close();
		String sampsLine = writeSamples(compress,zipSamples,nme+".sample");
		OutputStreamWriter pw = 	compress.getWriter("Name",true);
		pw.write("geno\n");
		pw.write("chr\tstart\tend\tsnpid\timputed\n");
		pw.write(sampsLine+"\n");
		pw.write("1000,0,0\n");
		compress.closeWriter(pw);
		compress.run();
	}
	}
		br.close();
		zipSamples.close();
		
		//compress.run();
 }
 
 private static String writeSamples(CompressDir compress2, ZipFile samplesFile, String ent) throws Exception {
	 String entry = "";
	 for(Enumeration en =  samplesFile.entries();en.hasMoreElements();){
		 String nme = ((ZipEntry) en.nextElement()).getName();
		if(nme.indexOf(ent)>=0){
			entry = nme;
			break;
		}
	 }
	System.err.println("matching entry is "+entry+" for "+ent);
	BufferedReader br =Compressor.getBufferedReader(samplesFile, entry);
	String str = br.readLine().replace(' ','\t');
	OutputStreamWriter ow = compress2.getWriter("Samples", true);
	String st = br.readLine();
	while((st=br.readLine())!=null){
		ow.write(st.replace(' ','\t')+"\n");
	}
	compress2.closeWriter(ow);
	return str;
}

private static void  getKary(File karyfile, String substring, List<Integer> res, List<String> app) throws Exception{
	String st = "";
	BufferedReader br = getBR(karyfile);
	while((st = br.readLine())!=null){
		String[] str = st.split("\\s+");
		if(str[0].equals(substring)){
			//Integer[] res = new Integer[str.length-1];
			for(int k=1; k<str.length; k++){
				res.add( Integer.parseInt(str[k]));
				if(k==1) app.add("p");
				else app.add(k+"");
			}
			res.add(Integer.MAX_VALUE-1);
			app.add("q");
			return;
		}
	}
	
	
}

public static void mainTranspose(BufferedReader br, File dir1, String[] rsid) throws Exception{
		
		OutputStreamWriter pw = 	compress.getWriter("Name",true);
		pw.write("countAll\tstate.0\tstate.1\tstate.2\n");
		pw.write("chr\tstart\tend\tsnpid\n");
		pw.write("sample\n");
		compress.closeWriter(pw);
		List<String>samples = new ArrayList<String>();
		for(int k=2; k<rsid.length; k+=4){
			samples.add(rsid[k].substring(0,rsid[k].lastIndexOf('_')));
		}
		pw = compress.getWriter("Samples",true);
			
			//new PrintWriter(new BufferedWriter(new FileWriter(new File(dir1,"Samples"))));
		for(int i=0; i<samples.size(); i++){
			pw.write(samples.get(i)+"\n");
		}
		compress.closeWriter(pw);
		List<String> snps = new ArrayList<String>();
		String st = "";
		while((st=br.readLine())!=null){
		//	System.err.write((st);
			if(st.length()==0) {
				continue;
			}
			String[] str = st.split("\\s+");
			snps.add(chr+"\t"+str[1]+"\t"+(Integer.parseInt(str[1])+20)+"\t"+str[0]);
			pw = compress.getWriter(str[0],true);
				
			for(int k=2; k<str.length; k+=4){
				for(int j=k; j<k+4; j++){
					pw.write(j<str.length ? str[j]:"NaN");
					if(j<k+3) pw.write("\t");
				}
				pw.write("\n");
			}
			compress.closeWriter(pw);
		}
		pw = 	compress.getWriter("SNPS",true);
		for(int i=0; i<snps.size(); i++){
			pw.write(snps.get(i)+"\n");
		}
		compress.closeWriter(pw);
		
}
static String chr;	
public static CompressDir compress=null; 

public static void  mainImpute(File dir, String nme, String sampleName, Integer max) throws Exception{
	
	 File kary = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
			 return  pathname.getName().startsWith("kary");
			}
			 
		 })[0];
	 File samples = new File(sampleName);
	 //for(int k=0; k<f.length; k++){
		 BufferedReader br = getBR(new File(dir,nme));
		  compressImpute(br, dir,nme.split("\\.")[0]	,kary, samples,max);
			
		 
	 //}
}

private static BufferedReader getBR(File f) throws Exception{
	return f.getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f)))) : 
		 new BufferedReader(new FileReader(f));
}

public static void main(File dir, boolean remove) {
	 try{
		// if(true) return;
		//	File dir = new File(System.getProperty("user.dir")+"/"+args[0]);
	 File avg = new File(dir,"avg");
	 if(!avg.exists()) return;
	 File[] f = avg.listFiles(new FileFilter(){

		@Override
		public boolean accept(File pathname) {
		 return pathname.getName().endsWith(".txt") && ! pathname.getName().startsWith("res")&& ! pathname.getName().startsWith("hwe");
		}
		 
	 });
	 chr = "chr"+dir.getName().split("_")[2];
	 for(int k=0; k<f.length; k++){
		 if(f[k].length()==0) continue;
		
		 BufferedReader br = 
			 f[k].getName().endsWith(".gz") ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f[k])))) : 
			 new BufferedReader(new FileReader(f[k]));
		 String str =  br.readLine();
		
		 String[] rsid =str.split("\t");
		
		 File dir1 = new File(avg,f[k].getName().split("\\.")[0]);
		  compress = (new CompressDir(dir1));
		  if(rsid.length>1 && rsid[1].equals("loc")) 
			 mainTranspose(br, dir1, rsid);
		 else
		 main(br, dir1, rsid);
		
			if(remove)f[k].delete();
		compress.run();
	 }
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
	
 }
}
