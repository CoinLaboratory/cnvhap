package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class ConvertToTable {
	// public static String dir = System.getProperty("user.dir");
//	 public static String prefix = "";
	//  public static String[] chr = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22".split(":");
/* first is location, second is file */
	  public static void main(String[] args){
	try{
		if(args.length==1 || true){
		run(new File(args[0]), null, null, Integer.parseInt(args[1]), Integer.parseInt(args[2]));
		}
		else if(args.length==2){
			run(new File(args[0]), new File(args[1]), null,null, null);
		}
		else if(args.length==3){
			run(new File(args[0]), new File(args[1]), new File(args[2]), null, null);
		}
	/*	File [] f = (new File(dir)).listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.getName().startsWith(prefix) && pathname.isDirectory();
			}
			
		});
		for(int i=0; i<f.length; i++){
			run(f[i]);
		}*/
	
		
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
	  
	  public static void printline(PrintWriter pw, String[] l, Set<Integer> indiv_exclude) throws Exception{
		  for(int i=0; i<l.length; i++){
			  if(!indiv_exclude.contains(i)){
				  pw.print(" "+l[i]);
			  }
		  }
		  pw.println();
	  }
	  public static int count(char[] res, char ch){
		  int cnt=0;
		  for(int i=0; i<res.length; i++){
			  if(res[i]==ch) cnt++;
		  }
		  return cnt;
	  }
	  public static void transform(String[] res){
		  Set<Character> chars = new HashSet<Character>();
		  for(int i=0 ;i < res.length; i++){
			  char[] ch = res[i].toCharArray();
			  for(int j=0; j<ch.length; j++){
				  chars.add(ch[j]);
			  }
		  }
		  char ch  = (new ArrayList<Character>(chars)).get(0);
		  
		  for(int i=0 ;i < res.length; i++){
			  res[i] = count(res[i].toCharArray(), ch)+"";
		  }
		  
	  }
static int col_ind = 3;
static boolean geno = false;
	  
	  public static void run(File name, File samples, File excludedSamples,Integer start, Integer end) throws Exception{
	String nme = name.getName().split("\\.")[0];
	PrintWriter res = new PrintWriter(new BufferedWriter(new FileWriter(new File(name.getParentFile(),nme+".txt"))));
	
	String[] toIgnore = "SNPS:Samples:Name".split(":");
	Set<String> exclude = new HashSet<String>();
	exclude.addAll(Arrays.asList(toIgnore));
		System.err.println("opening "+name.getAbsolutePath());
		ZipFile zf = new ZipFile(name);
		List<String>indiv = new ArrayList<String>();
		Set<Integer> indiv_exclusion = new HashSet<Integer>();
		if(samples!=null && samples.exists()){
			Utils.readPosInfo(samples, new int[] {0}, false, new List[]{indiv},new Class[] {String.class});
		}
		else{
			Compressor.read(zf, "Samples",indiv,  0);
		}
		if(excludedSamples!=null && excludedSamples.exists()){
			List<String> excluded_indiv = new ArrayList<String>();
			Utils.readPosInfo(excludedSamples, new int[] {0}, false, new List[]{excluded_indiv},new Class[] {String.class});
			for(int i=0; i<excluded_indiv.size(); i++){
				indiv_exclusion.add(indiv.indexOf(excluded_indiv.get(i)));
			}
			System.err.println("excluding indiv :"+excluded_indiv+" "+indiv_exclusion);
		}
		res.print("id");
		
		printline(res, indiv.toArray(new String[0]), indiv_exclusion);
		String[]ele = new String[indiv.size()];
		if(start!=null){
	BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
	String st = "";
	String chr = "chr"+name.getName().split("\\.")[0].split("p")[0].split("q")[0];
	while((st = br.readLine())!=null){
		String[] str = st.split("\\s+");
		if(str[0].split("p")[0].split("q")[0].equals(chr)){
			int pos = Integer.parseInt(str[1]);
			if(pos>= start && pos<=end){
				ZipEntry nxt = (ZipEntry) zf.getEntry(str[3]);
				if(!exclude.contains(nxt.getName())){
					Compressor.readZip(zf, nxt.getName(),ele,  col_ind, "\t");
					if(geno)transform(ele);
					res.print(nxt.getName());
					printline(res, ele, indiv_exclusion);
				}
			}
		}
	}
		}
		else{
	for(Enumeration en = zf.entries(); en.hasMoreElements();){
		ZipEntry nxt = (ZipEntry) en.nextElement();
		if(!exclude.contains(nxt.getName())){
			Compressor.readZip(zf, nxt.getName(),ele,  col_ind, "\t");
			if(geno)transform(ele);
			res.print(nxt.getName());
			printline(res, ele, indiv_exclusion);
		}
	}
		}
	
	res.close();
	zf.close();
}

}
