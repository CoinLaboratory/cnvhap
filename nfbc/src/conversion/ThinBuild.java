package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class ThinBuild {
	public static void main(String[] args){
		try{
			File f = new File("build36.txt");
			ThinBuild tb = new ThinBuild(f, 0.9);
			tb.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	BufferedReader br;
	PrintWriter pw;
	PrintWriter pw_excl;
	double acc_prob;
	public ThinBuild(File f, double acc_prob) throws Exception{
		br = new BufferedReader(new FileReader(f));
		 pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(f.getParentFile(), "build_out.txt"))));
		pw_excl = new PrintWriter(new BufferedWriter(new FileWriter(new File(f.getParentFile(), "build_excl.txt"))));
		this.acc_prob = acc_prob;
	}
	public void run() throws Exception{
		String st = "";
		while((st = br.readLine())!=null){
			if(Math.random()<acc_prob){
				pw.println(st);
			}
			else{
				pw_excl.println(st);
			}
		}
		pw.close();
		pw_excl.close();
	}
	
}
