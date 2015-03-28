package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

public class SubtractBuild {
	public static void main(String[] args){
		try{
			SubtractBuild sb = new SubtractBuild(new File("317k/build36.txt"), new File("1M/build36.txt"));
			sb.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	Set<String> rsToSubtract = new HashSet<String>();
	BufferedReader toP;
	PrintWriter out;
	SubtractBuild(File toSubtract, File toP) throws Exception{
		if(!toSubtract.exists() || !toP.exists()) throw new RuntimeException("!!");
		BufferedReader br = new BufferedReader(new FileReader(toSubtract));
		String st = "";
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			this.rsToSubtract.add(str[3]);
		}
		br.close();
		this.toP = new BufferedReader(new FileReader(toP));
		out = new PrintWriter(new BufferedWriter(new FileWriter(new File(toP.getParentFile(),toP.getParentFile().getName()+"_minus_"+toSubtract.getParentFile().getName()))));
	}

	public void run() throws Exception{
		String st = "";
		while((st = toP.readLine())!=null){
			String[] str = st.split("\\s+");
			if(!this.rsToSubtract.contains(str[3])){
				out.println(st);
			}
		}
		toP.close();
		out.close();
	}
}
