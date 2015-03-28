import java.io.File;
import java.io.FileFilter;


public class FileConvert {
	public static void main(String[] args){
		try{
			File user = new File(System.getProperty("user.dir"));
			File[] f = user.listFiles(new FileFilter(){

				public boolean accept(File pathname) {
					return pathname.getName().endsWith("sdf");
				}
				
			});
			for(int i=0; i<f.length; i++){
				String name = f[i].getName();
				name = name.substring(0,name.indexOf(".sdf"));
				final String name1 = name;
				File f_new = new File(name);
				f_new.mkdir();
				File[] f1 = user.listFiles(new FileFilter(){

					public boolean accept(File pathname) {
						return pathname.getName().startsWith(name1);
					}
					
				});
				for(int j=0; j<f1.length; j++){
					File dest = new File(f_new, f1[j].getName());
					f1[j].renameTo(dest);
				}
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
}
