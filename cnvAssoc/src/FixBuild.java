import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;


public class FixBuild {
public static void main(String[] args){
	try{
		BufferedReader br = new BufferedReader(new FileReader(new File(args[0])));
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[0]+".out"))));
		String st = "";;
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			int p1 = Integer.parseInt(str[1]);
			int p2 = Integer.parseInt(str[2]);
			if(p2 < p1){
				str[1] = p2+"";
				str[2] = p1+"";
			}
			pw.print(str[0]);
			for(int k=1; k<str.length; k++){
				pw.print("\t"+str[k]);
			}
			pw.println();
		}
		pw.close();
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
