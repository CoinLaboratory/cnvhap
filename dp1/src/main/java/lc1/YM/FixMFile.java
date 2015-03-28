package lc1.YM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;


public class FixMFile {
	
	public static void main(String[] args){
		try{
			FixMFile fmf  = new FixMFile(new File("chrM.tree.matrix"));
			fmf.print();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
  BufferedReader br;
  String[] basic = new String[] {"chr","pos"};
  PrintWriter out;
  SortedMap<Integer, List<String>> m = new TreeMap<Integer, List<String>>();
  public void print(){
	   basic = new String[] {"MT","0"};
outer:	  for(Iterator<Integer> it = m.keySet().iterator(); it.hasNext();){
		 Integer key = it.next();
		  List<String> st = m.get(key);//
		  basic[1] = key+"";//.split("[A-Z]")[0].replace("!", "").split("\\.")[0];
		  if(st.size()==1){
			  String[] str_  = st.get(0).split("\t");
			
			if(!str_[0].startsWith("(") && str_[0].indexOf('d')<0 && str_[0].indexOf('.')<0){
				out.println(getStr(basic,str_,""));
			}
		  }else if(st.size()==2){
			  String[][] str = new String[st.size()][];
			  int[][] vals = new int[st.size()][];
			  int[] sum = new int[st.size()];
			  int rev_ind = -1;
			  int anc_ind = 0;
			  for(int k=0; k<str.length; k++){
				  str[k] = st.get(k).split("\t");
				  vals[k] = new int[str[k].length-1];
				  for(int j=1; j<str[k].length; j++){
					  vals[k][j-1] = Integer.parseInt(str[k][j]);
					  sum[k]+=vals[k][j-1];
				  }
				  if(str[k][0].indexOf('!')>=0){
					  if(rev_ind>=0) continue;
					  rev_ind = k;
				  }
				  if(sum[k] > sum[anc_ind]) anc_ind = k;
			  }
			  if(rev_ind>=0){
				  String[] res = new String[str[0].length];
				  String[] der = str[1-anc_ind];
				  String[] anc = str[anc_ind];
				  System.arraycopy(anc, 0, res, 0, res.length);
				  for(int i=1; i<res.length; i++){
					  if(der[i].equals("1")){
					     if(anc[i].equals("1") ){
						     res[i] = "0";
					     }else if(anc[i].equals("0") ){
					    	 //res[i] = "0";
					    	 System.err.println("cannot explain "+key);
						  //  throw new RuntimeException("!!");
					    	 continue outer;
					     }
					  }
				  }
				  out.println(getStr(basic,res,""));
			  }
			  
		  }
	  }
	  out.close();
  }
  
  
  private String getStr(String[] basic, String[] res, String suff) {
	StringBuffer sb = new StringBuffer(basic[0]);
	for(int i = 1; i<basic.length; i++){
		sb.append("\t"+basic[i]);
	}
	for(int i = 0; i<res.length; i++){
		if(i==0) sb.append("\t"+res[i]);
		else sb.append("\t"+res[i]+suff); 
	}
	return sb.toString();
}


public FixMFile(File f) throws Exception{
	  br = new BufferedReader(new FileReader(f));
	  out = new PrintWriter(new FileWriter(new File(f.getParentFile(), f.getName()+".out")));
	  String st = br.readLine().replaceAll("path_point", "info");
	  //st = st.replaceAll("haplo", "geno");
	  out.println(getStr(basic,st.split("\t"),".geno"));
	  List<String> h = Arrays.asList(st.split("\t"));
	  while((st=br.readLine())!=null){
		  String[] str = st.split("\t");
	  if(str[0].indexOf('-')>=0 || str[0].startsWith("res") 
			  ) continue;
		  else{
			  
			  Integer rep = Integer.parseInt(str[0].split("[A-Z]")[0].replaceAll("!", "").
					  replaceAll("d", "").replaceAll("\\(","").replaceAll("\\)","").replaceAll("\\.", ""));
			  //String rep = getString(str1,str1.length-1);
			  List<String> l = m.get(rep);
			  if(l==null) m.put(rep, l=new ArrayList<String>(2));
			  l.add(st);
		  //}
	  }
	  }
  }


private String getString(String[] str1, int i) {
	StringBuffer sb  = new StringBuffer();
	for(int k=0; k<i; i++){
		sb.append(str1[k]);
	}
	return sb.toString();
}
}
