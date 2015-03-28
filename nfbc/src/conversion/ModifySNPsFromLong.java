package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

public class ModifySNPsFromLong {
	public static void main(String[] args){
		try{
			ModifySNPsFromLong msfl = new ModifySNPsFromLong(new File("snps.txt"), new File("diff.txt"));
			msfl.run(new File("build36.txt"), new File("out.txt"));
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	Map<String, Double> gc = new HashMap<String, Double>();
	Map<String, SortedMap<Integer, String>> missed = new HashMap<String, SortedMap<Integer,String>>();
	
	ModifySNPsFromLong(File snps, File diff) throws Exception{
		
		
			Set<String> missl = new HashSet<String>(Compressor.getIndiv(diff, 0));
			
		
		
		BufferedReader br = Utils.getBufferedReader(snps);
		String st = "";
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			if(missl.contains(str[0])){
				String chr = str[1];
				SortedMap< Integer, String> m = missed.get(chr);
				if(m==null){
					missed.put(chr, m = new TreeMap<Integer,String>());
				}
				m.put(Integer.parseInt(str[2]), str[0]);
			}
		}
		br.close();
		
	}
	
	public void run(File in, File outf) throws Exception{
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outf)));
		BufferedReader br = Utils.getBufferedReader(in);
		String st = "";
		SortedMap<Integer, String> currMap=null;
		String currentChr = "";
		Integer first = null;
		String prev = null;
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			String chr = str[0];
			if(!chr.equals(currentChr)){
				currentChr = chr.substring(3);
				currMap = missed.get(currentChr);
				if(currMap!=null && currMap.size()>0){
					first = currMap.firstKey();
				}
			}
			if(first!=null){
				Integer pos = Integer.parseInt(str[1]);
				
				if(pos>=first){
					//print
					String idn = currMap.get(first);
					if(pos==first){
						System.err.println(idn+" "+str[3]);
						
					}
					out.print("chr"+currentChr+"\t"+pos+"\t"+(pos+20)+"\t"+idn+"-\t-\tNaN");

					currMap = currMap.tailMap(first+1);
					if(currMap.size()>0) first = currMap.firstKey();
					else first = null;
				}
			}
			if(prev!=null){
				out.println(prev);
			}
			prev = st;
		}
		out.println(prev);
		br.close();
	}
	//Buffered
	
}
