package lc1.dp.appl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;

import lc1.util.Constants;

public class MakeReps {
	public static void main(String[] args){
		try{
		 Map<String, Integer> m = Constants.readReps(new File(args[0]));
		 PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(args[1]))));
		 pw.println("--index 1");
		 for(Iterator<Map.Entry<String, Integer>> it = m.entrySet().iterator(); it.hasNext();){
				Map.Entry<String,  Integer> st_ = it.next();
				for (int i=1; i<st_.getValue(); i++){
					pw.println("--index "+(i+1)+" --indexControl "+st_.getKey());
				}
			}
		 pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
}
