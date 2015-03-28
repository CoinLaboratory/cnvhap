import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;


public class MakeNewBuild {
	
	public static void main(String[] args){
		File dir = new File(".");
		MakeNewBuild mnb = new MakeNewBuild(dir, args[0], args[1]);
		mnb.run();
	}

	File build1;
	File build2;
	
	PrintWriter b1minusb2;
	
	Set<String> l2 = new HashSet<String>();
	
	/** build2 smaller */
	MakeNewBuild(File dir, String build1, String build2){
		 this.build1 = new File(dir, build1+".txt");
		 this.build2 = new File(dir,build2+".txt");
		 try{
			 this.b1minusb2 = 
				 new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, build1+"_"+build2))));
			 BufferedReader br1 = new BufferedReader(new FileReader(this.build2));
			 String st = "";
			 while((st = br1.readLine())!=null){
				 l2.add(st.split("\\s+")[3]);
			 }
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
	}
	
	public void run(){
		try{
			 BufferedReader br1 = new BufferedReader(new FileReader(this.build1));
			 String st = "";
			 while((st = br1.readLine())!=null){
				 if(!l2.contains(st.split("\\s+")[3])){
					b1minusb2.println(st); 
				 }
			 }
			 b1minusb2.close();
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
	}
	
}
