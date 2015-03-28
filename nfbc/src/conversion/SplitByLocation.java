package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public class SplitByLocation {
public static void main(final String[] args){
	try{
	//	if(true) System.exit(0);
	File user = new File(args[0]);
	File output = new File(args[1]);
	if(!output.exists()){
		output.mkdir();
	}
	File buildF = new File(args[2]);
	File karotypes = new File(buildF.getParentFile(), "karyotypes_"+buildF.getName().replace(".gz", "") );
	String length = args.length<=3 ? null : args[3];
	if(length!=null){
	if(length.endsWith("kb")){
		mult = 1000;
		length = length.substring(0,length.length()-2);
	}
	else if(length.endsWith("mb")){
		mult = 1000*1000;
		length = length.substring(0,length.length()-2);
	}
	}
	
	SplitByLocation spl = new SplitByLocation(user, buildF, output,karotypes,length==null ? null: Integer.parseInt(length));
	spl.run();
	}catch(Exception exc){
		exc.printStackTrace();
	}
	/* File excel = user.listFiles(new FileFilter(){

	public boolean accept(File pathname) {
	return pathname.getName().endsWith(".xls");
	}
	 
 })[0];*/
	
	
}

public static void split(String tod, File di, File buildF, File karyotypes, boolean deleteOriginal) throws Exception {
	// TODO Auto-generated method stub
	File user = new File(di, tod+".zip");
	File output =di;
	SplitByLocation spl = new SplitByLocation(user, buildF,  output,karyotypes, null);
	spl.run();
	if(deleteOriginal){
		File f1 = new File(di, tod+"q.zip");
		File f2 = new File(di, tod+"p.zip");
		if(f1.exists() &&f1.length()>0 && (new ZipFile(f1)).size() >4){
			//long len = f1.length();
			
			user.delete();
		}
		else if(f2.exists() &&f2.length()>0 && (new ZipFile(f2)).size() >4){
			user.delete();
		}
	}
}

public  String prepareBuild(BufferedReader br) throws Exception{
	String buildString = br.readLine();
	String pref = "chr"+chr;
	while(!buildString.startsWith(pref)){
		buildString = br.readLine();
	}
	return buildString;
}

BufferedReader br;
File outDir;
Integer mb;
String buildString = "";
String[] build;
FileOutputStream dest;
CheckedOutputStream checksum ;
ZipOutputStream outS ;
OutputStreamWriter osw ;
ZipFile zf;
String chr;
PrintWriter snps;
File snpsF;
List<String> indiv, names;
int prev;
static String sep = "\t+";
static int pos = 1;
static int id = 3;
int max, max_bp;
static int chrid =0;
static int mult = 1000;
String name;
String chrname;




boolean leftRestOnly = false;

SplitByLocation(File in, File buildF, File outDir, File karyotypes, Integer mb) throws Exception{
	this.outDir = outDir;
	  chr = in.getName().toString().split("\\.")[0];
	  chrname = "chr"+chr;
	
	
	  zf = new ZipFile(in);
	  br = Utils.getBufferedReader(zf, "SNPS");
	  if(br==null){
	  br = Utils.getBufferedReader(buildF);
		 if(br==null){
			br = Utils.getBufferedReader(new File(buildF.getParentFile(),buildF.getName().replaceAll(".gc", "")));
		 }
	  }
	  buildString = prepareBuild(br);
		 build = buildString.split(sep);
	  indiv = Compressor.getIndiv(zf, "Samples");
	   names =Compressor.getIndiv(zf, "Name");
	   this.mb = mb;
	   if(mb==null){
		   Integer centromere=null;
		  /* if(!karyotypes.exists() || karyotypes.length()==0){
			   GenesForRegion.writeKaryotypeFile(karyotypes);
			   
		   }*/
		   
		  centromere = readCentromere(karyotypes).get(chr);
		   if(centromere==null){
			   throw new RuntimeException("no centromere "+chr+" "+karyotypes);
		   }
		 		  
		   this.mb = (int) Math.floor(centromere.doubleValue()/mult);
		   leftRestOnly=true;
	   }
}

public static Map<String, Integer> readCentromere(File karyotypes) throws Exception{
	if(!karyotypes.exists()) return null;
	 BufferedReader br = Utils.getBufferedReader(karyotypes);
	   String st = "";
	   Map<String, Integer> m = new HashMap<String, Integer>();
	   inner: while((st = br.readLine())!=null){
		   if(st.startsWith("M\t")) st = st.replace("M", "Mt");
		   String[] str = st.split("\t");
		   m.put(str[0], Integer.parseInt(str[1]));
		 
		   
	   }
	   br.close();
	   return m;
	
}

public void transferSNPS() throws Exception{
	snps.close();
	String oldname = "chr"+chr;
	String newname =  "chr"+name;
	List<String> resu = Compressor.getIndiv(snpsF);
	for(int i=0; i<resu.size(); i++){
		resu.set(i, resu.get(i).replace(oldname,newname));
	}
	Compressor.writeEntry("SNPS", resu, outS, osw);
}

int cnt=0;

public void mknewOutput() throws Exception{
	if(outS!=null){
		transferSNPS();
		outS.close();
		osw.close();
	}
	prev = (int) Math.floor(Double.parseDouble(build[pos])/((double)mult));
	max =this.leftRestOnly  ? ( cnt==0  && mb >=prev ? mb : Integer.MAX_VALUE/mult) :  prev + mb;
	name = chr+(this.leftRestOnly  ? (cnt==0 && mb >=prev? "p" :"q" ): +prev+"-"+max);
	cnt++;
	max_bp = mult * max;
	
	snpsF = new File(outDir, "snps");
	snps = new PrintWriter(new BufferedWriter(new FileWriter(snpsF)));

	File outF = new File(outDir, name+".zip");
	  dest = Compressor.getOS(outF);
      checksum = new CheckedOutputStream(dest, new Adler32());
      outS = new 
       ZipOutputStream(new 
         BufferedOutputStream(checksum));
     osw = new OutputStreamWriter(outS);
     outS.setMethod(ZipOutputStream.DEFLATED);
     if(names!=null) Compressor.writeEntry("Name", names, outS, osw);
     if(indiv!=null)  Compressor.writeEntry("Samples", indiv, outS, osw);
     
}

public void run() throws Exception{
	while(this.build!=null  && build[chrid].equals(chrname)){
		mknewOutput();
		while(build!=null&& Integer.parseInt(build[pos]) < max_bp && build[chrid].equals(chrname)){
		try{
		  ZipEntry ent = (ZipEntry) zf.getEntry(build[id]);
		  if(ent!=null){
				snps.println(buildString);
	        List<String> res = Compressor.getIndiv(zf, ent.getName());
	            Compressor.writeEntry(ent.getName(), res,outS, osw);
	      
		  }
		  else{
			  System.err.println("zip was missing "+build[id]);
		  }
		}catch(Exception exc){
			exc.printStackTrace();
		}
		  this.buildString = br.readLine();
	       this.build = buildString==null ? null : buildString.split(sep);
		}
	}
	transferSNPS();
		outS.close();
		osw.close();
	zf.close();
	br.close();
	snps.close();
	snpsF.deleteOnExit();
	
}




}
