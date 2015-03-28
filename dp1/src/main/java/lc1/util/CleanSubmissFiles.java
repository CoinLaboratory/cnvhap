package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;

public class CleanSubmissFiles {
public static void main(String[] args){
	try{
	       File user = new File(System.getProperty("user.dir"));
		File[] f = user.listFiles(new FileFilter(){

			public boolean accept(File pathname) {
				return pathname.getName().indexOf(".e")>=0;
			}
			
		});
		File error = new File(user, "error");
		error.mkdir();
		File done = new File(user, "done");
		done.mkdir();
		for(int i=0; i<f.length; i++){
			String name = f[i].getName();
			final String name1 = name.substring(0, name.indexOf(".e"));
			File[] f1 = user.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.getName().startsWith(name1);
				}
				
			});
		
			boolean haserror  = hasError(f[i]);
			if(haserror){
				move(f1, error);
			}
			else{
				move(f1, done);
			}
		}
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

private static void move(File[] f1, File done) {
	for(int i=0; i<f1.length; i++){
		File newF = new File(done, f1[i].getName());
		if(newF.exists()) newF.delete();
		f1[i].renameTo(newF);
	}
}

private static boolean hasError(File file) {
	try{
		BufferedReader br = new BufferedReader(new FileReader(file));
		String st = "";
		while((st = br.readLine())!=null){
			if(st.indexOf("Exc")>=0){
				br.close();
				return true;
			}
		}
	}catch(Exception exc){
		exc.printStackTrace();
	}
	return false;
}
}
