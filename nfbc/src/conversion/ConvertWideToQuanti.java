package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class ConvertWideToQuanti {
BufferedReader br;

PrintWriter[] out;
public static void main(String[] args){
	try{
		File in = new File(args[0]);
		new ConvertWideToQuanti(in, new File("out"));
		
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

ConvertWideToQuanti(File in, File outDir) throws Exception{
	this.br = Utils.getBufferedReader(in);
	outDir.mkdir();
	String[] head = br.readLine().split("\t");
	List<String> samples = new ArrayList<String>();
	List<String> type = new ArrayList<String>();
	 getSamples(head, 3, samples, type);
	 this.out = new PrintWriter[samples.size()];
	 indices = new int[out.length][5];
	 for(int k=0; k<out.length; k++){
		 out[k] = new PrintWriter((new FileWriter(new File(outDir, samples.get(k)))));
		  indices[k][0] = 0;
		  indices[k][1] = 1;
		  indices[k][2] = 2;
		  indices[k][3] = 3 +1 + type.size()*k;
		  indices[k][4] = 3 + 2 + type.size()*k;
		//  print(out[k], head, indices[k]);
		  out[k].print(head[indices[k][0]]);
		  out[k].print("\t"+head[indices[k][1]]);
		  out[k].print("\t"+head[indices[k][2]]);
		  out[k].print("\t"+head[indices[k][3]].split("\\.")[1]);
		  out[k].print("\t"+head[indices[k][4]].split("\\.")[1]);
		  out[k].println();
	 }
	 String st = "";
	 while((st = br.readLine())!=null){
		 String[] str = st.split("\t");
		 for(int k=0; k<out.length; k++){
			 print(out[k], str, indices[k]);
		 }
	 }
	 for(int k=0; k<out.length; k++){
		out[k].close();
	 }
	 br.close();
	 
}

private void print(PrintWriter pw, String[] head, int[] is) {
	pw.print(head[is[0]]);
	for(int k=1; k<is.length; k++){
		pw.print("\t"+head[is[k]]);
	}
	pw.println();
	
}

int[][] indices;

int[] indToPrint = new int[] {1,2};
private void getSamples(String[] head, int st, List<String> samples,
		List<String> type) {
	int i =st;
	String firstName = head[i].split("#")[0];
	while(head[i].split("#")[0].equals(firstName)){
		type.add(head[i].split("\\.")[1]);
		i++;
	}
	for(int j=st; j<head.length; j+=i){
		samples.add(head[j].split("#")[0]);
	}
	
	
}

}
