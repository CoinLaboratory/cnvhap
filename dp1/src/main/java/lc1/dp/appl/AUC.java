package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AUC {
	
	public static void main(String[] args){
		try{
		File f =  new File(System.getProperty("user.dir"));
		runAll(f);
	/*	AUC auc = new AUC(f);
		auc.run();*/
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	File dir;
	PrintWriter pw;
	
	
	static void runAll(File dir) throws Exception{
		if(dir.getName().startsWith("ROC")) {
			System.err.println("running "+dir.getAbsolutePath());
			AUC auc = new AUC(dir);
			auc.run();
		}
		else if(dir.isDirectory()){
			File[] f = dir.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.isDirectory();
				}
				
			});
			for(int i=0; i<f.length; i++){
				runAll(f[i]);
			}
		}
	}
	AUC(File dir) throws Exception{
		this.dir = dir;
		pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "auc.txts"))));
	}
	
	
	void run() throws Exception{
		File[] f = dir.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.isDirectory();
			}
			
		});
		if(f.length>0){
			List<String> head;
			{
		Map<String, Double> m = runFile(f[0]);
	 head = new ArrayList<String>(m.keySet());
			}
		pw.println(String.format("%-7s", "type")+String.format(getFormat("%7s", head.size()), head.toArray()));
		//pw.println(String.format("%-7s", f[0].getName())+String.format(getFormat("%5.3g", head.size()), getArray(head, m)));
		for(int i=0; i<f.length; i++){
			Map<String, Double> m = runFile(f[i]);
			String st_ =  f[i].getName().substring(9);
			int ind = st_.indexOf("in");
			String[] str =new String[] { st_.substring(0, ind), st_.substring(ind, st_.length())};
			pw.println(String.format("%-7s",
					str[0]+"/"+str[1])+String.format(getFormat("%5.3g", head.size()), getArray(head, m)));
		}
		pw.close();
		}
	}
	
	private Object[] getArray(List<String> head, Map<String, Double> m) {
		Object[] res = new Object[head.size()];
		for(int i=0; i<head.size(); i++){
			res[i] = m.get(head.get(i));
		}
		return res;
	}


	private String getFormat(String string, int i) {
		StringBuffer sb = new StringBuffer(string);
		for(int k=1; k<i; k++){
			sb.append(" "+string);
		}
		return sb.toString();
	}


	Map<String, Double> runFile(File dir) throws Exception{
		File[] f = dir.listFiles();
		Map<String, Double> m = new HashMap<String, Double>();
		for(int i=0; i<f.length; i++){
			
			m.put(f[i].getName(),run(f[i]) );
		}
		return m;
		
	}
	
	double run(File f) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(f));
		String st = "";
		double x_prev =0, y_prev =0, auc =0;
		while((st = br.readLine())!=null){
			String[] str = st.split("\\s+");
			double x=  Double.parseDouble(str[0]);
			double y=  Double.parseDouble(str[1]);
			auc += y_prev *(x - x_prev);
			auc+= 0.5* (y - y_prev)*(x-x_prev);
				x_prev = x;
			y_prev = y;
		}
		auc = auc /(x_prev*y_prev);
		br.close();
		return auc;
	}
	
}
