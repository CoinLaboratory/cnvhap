import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;


public class ConvertPhenoNames {
public static void main(String[] args){
	try{
		//if(true) System.exit(0);
		String str = "1\t\t\t1";
		String[] st = str.split("\t");
	//	if(true) return;
	ConvertPhenoNames cpn =
		new ConvertPhenoNames(new File("plate.txt"));//b35-NFBC-release3.original.fam"));
	cpn.convert(new File(args[0]), new File("variance.txt"));
	}catch(Exception exc){
		exc.printStackTrace();
	}
	
}
	Map<String, String> conversion = new HashMap<String, String>();
	Map<String, String> conversionRev = new HashMap<String, String>();
	Map<String, String> variance = new HashMap<String, String>();
	
	public static Map<String, String> readVar(File in) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(in));
		String st = br.readLine();
		String[] str = st.split("\\s+");
		Map<String, String> variance = new HashMap<String, String>();
		
		while((st = br.readLine())!=null){
			str = st.split("\\s+");
			variance.put(str[0], str[1]);
		}
		br.close();
		return variance;
	}
	ConvertPhenoNames(File in) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(in));
		String st = br.readLine();
		String[] str = st.split("\\s+");
		int ind1 = Arrays.asList(str).indexOf("PATIENT");
		int ind2 = Arrays.asList(str).indexOf("ALTERNATIVE");
		while((st = br.readLine())!=null){
			str = st.split("\\s+");
			this.conversion.put(str[ind2], str[ind1]);
			this.conversionRev.put(str[ind1], str[ind2]);
		}
		br.close();
	}
	public void convert(File in, File var) throws Exception{
		variance = readVar(var);
		BufferedReader br = new BufferedReader(new FileReader(in));
		File out = new File(in.getParentFile(), in.getName()+"_1");
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
		String st = br.readLine();
		pw.print(st);
		pw.println("\tVariance");
		int len = st.split("\t").length;
		inner: for(int i=0; (st = br.readLine())!=null; i++){
			String[] str = st.split("\t");
			//int len = str.length;
			
			System.err.println(len);
			String conv = conversion.get(str[0]);
			if(conv==null){
			//	System.err.println("failed for "+str[0]);
				continue inner;
			}
			pw.print(str[0]);
			for(int k=1; k<str.length; k++){
				pw.print("\t");
				pw.print(str[k].trim().length()==0 ? "NA" : str[k]);
			}
			for(int k=str.length; k<len; k++){
				pw.print("\t");
				pw.print("NA");
			}
			pw.println("\t"+variance.get(conv));
		}
		pw.close();
		br.close();
	}
	
	
	
}
