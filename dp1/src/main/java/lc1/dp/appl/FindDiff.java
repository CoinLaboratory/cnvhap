package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class FindDiff {

	
	public static void main(String[] args){
		try{
			boolean[] tf = new boolean[] {false, true};
			for(int i=0; i<tf.length; i++){
				for(int j=0;j<tf.length; j++){
		FindDiff fd = new FindDiff(tf[i], tf[j]);
		fd.run();
				}
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	BufferedReader[] br= new BufferedReader[2];
	Map<String, String>[]map = new Map[2];
	PrintWriter[] out = new PrintWriter[2];
	
	FindDiff(boolean all, boolean del) throws Exception{
		String base = "roc_base_";
		base = base +(all ? "all_" : "max_");
		base = base +(del ? "0.1." :"3.4.");
		base = base +"_index_244k_ADM2_line_";
		String[] base1 = new String[] { base+"1M_PennCNV.txt",base+"1M_cnvHap.txt"};
	 
	
		
		for(int i=0; i<br.length; i++){
			File f = new File(base1[i]);
			System.err.println(f);
			this.br[i] = new BufferedReader(new FileReader(base1[i]));
			map[i] = new HashMap<String, String>();
			this.out[i] = new PrintWriter(new BufferedWriter(new FileWriter(new File(base1[i]+".out"))));
		}
		
		
		
		
	}
	public void run() throws Exception{
		read(0); read(1);
		findDiff(0); findDiff(1);
	}
	private void findDiff(int comp) {
		for(Iterator<String> it = this.map[comp].keySet().iterator(); it.hasNext();){
			String nxt = it.next();
			if(!map[1-comp].containsKey(nxt)){
				this.out[comp].println(map[comp].get(nxt));
			}
		}
		out[comp].close();
	}
	public void read(int i) throws Exception{
		String st = br[i].readLine();
		System.err.println(st);
		while((st = br[i].readLine())!=null ){
		//	System.err.println(st);
			st = st.trim();
			if(st.startsWith("positive")){
				String[] str = st.split("\\s+");
				double pp = Double.parseDouble(str[2]);
				if(pp>0.5){
					map[i].put(str[3],st);
				}
			}
		}
		br[i].close();
	}
	
}
