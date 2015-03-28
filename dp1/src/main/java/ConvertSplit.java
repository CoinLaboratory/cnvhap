import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import lc1.dp.appl.CNVHap;


public class ConvertSplit {
public static void main(String[] args){
	try{
		String[][] res1 = CNVHap.getMid(new File("split.txt"), null);
		String[][] res = mergeSplit(1000*1000,res1);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("split_out.txt"))));
		for(int i=0; i<res.length; i++){
			pw.println("--mid  "+res[i][0]+":"+res[i][1]+":"+res[i][2]);
		}
		pw.close();
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

public static String[][] mergeSplit(int dist, String[][]res) throws Exception{
	Map<String, List<String[]>> m = new HashMap<String, List<String[]>>();
	List<String[]> l_new = new ArrayList<String[]>();
	for(int k=0; k<res.length; k++){
		String chr = res[k][1].split(":")[0];
		List<String[]> l = m.get(chr);
		if(l==null){
			m.put(chr, l = new ArrayList<String[]>());
		}
		l.add(res[k][1].split(":"));
		Collections.sort(l, new Comparator<String[]>(){

			@Override
			public int compare(String[] arg0, String[] arg1) {
				int i1 = Integer.parseInt(arg0[1]);
				int i2 = Integer.parseInt(arg1[1]);
				if(i1!=i2) return i1<i2 ? -1 : 1;
				return 0;
			}
			
		});
	}
	for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
		String key = it.next();
		List<String[]> l = m.get(key);
		
		int pos =0;
		while(pos<l.size()){
			String[] str = l.get(pos);
			int p1 = Integer.parseInt(str[1]);
			int p2 = p1;
			//int loc = new int[]{p1,p1};
			boolean br=false;
			for(; !br && pos<l.size(); pos++){
				int j=pos;
				 str = l.get(j);
				int locn = Integer.parseInt(str[1]);
				if(locn < p2+dist){
					p2 = locn;
				}
				else{
					pos--;
					//pos = j+1;
					br=true;
				}
			}
			l_new.add(new String[] {key, p1+"", p2+""});
			
		}
		
		
	}
	return l_new.toArray(new String[0][]);
}
}
