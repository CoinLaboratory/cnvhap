package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class AffyBuildPrep {

	ZipFile zf;
	PrintWriter pw;
	Map<String, Map<Integer, String>> m = new HashMap<String, Map<Integer, String>>();
	public void read() throws Exception{
		for(Enumeration en = zf.entries(); en.hasMoreElements();){
      	  ZipEntry ent = (ZipEntry) en.nextElement();
      	//  String[] str = st.split("\\s+");
      	 // String rsid = str[3];
      	 List<String> l =  Compressor.read(zf, ent);
      	 
      	 if(l.size()>6){
      		 try{
      	 int pos = Integer.parseInt(l.get(5));
      	 String chr = l.get(2);
      	 if(!this.chr.contains(chr)) continue;
      	;
      	 String rsid = ent.getName();
      	 String res = "chr"+l.get(2)+"\t"+pos+"\t"+(pos+40)+"\t"+rsid+"\t"+l.get(3)+"\t"+l.get(4)+"\t"+l.get(6)+"\t"+l.get(1);
      
      	pw.println(res);
      		 }catch(Exception exc){
      			 exc.printStackTrace();
      		 }
      /*	 Map<Integer, String> m1 = m.get(chr);
      	if(m1==null){
      		m.put(chr, m1 = new TreeMap<Integer, String>());
      	}
      	m1.put(pos, res);*/
      	 }
      	 }
		pw.close();
	}
	
	public void print(){
		System.err.println("printing");
		for(Iterator<Map<Integer, String>> it = m.values().iterator(); it.hasNext();){
			for(Iterator<String> it1 = it.next().values().iterator(); it1.hasNext();){
		 pw.println(it1.next());
			}
		}
	        pw.close();
	}
	final Set<String> chr;
	public AffyBuildPrep(File dir1, String[] chr) throws Exception{
		  zf = new ZipFile(new File(dir1, "SnpToGene-Affy.zip"));
	         // String st = "";
	        pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir1, "build35-out.txt"))));
	        this.chr = new HashSet(Arrays.asList(chr));
	}
	public static void main(String[] args){
		try{
			 File dir1 = new File(System.getProperties().getProperty("user.dir"));
	      //    File in = new File(dir1, "build35.txt");
	  //  if(!in.exists()) return;
	        //  BufferedReader br = new BufferedReader(new FileReader(in));
	        AffyBuildPrep afbp = new AffyBuildPrep(dir1, "8".split(":"));
	          afbp.read();//afbp.print();
	        //  br.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
}
