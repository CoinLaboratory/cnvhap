package conversion;

import java.io.File;
import java.io.FilenameFilter;

public class ConvertToByIndiv {
    File out;
    public ConvertToByIndiv(File f) {
        out = new File(f.getParentFile(),f.getName().substring(0, f.getName().indexOf("zip")));
       out.mkdir();
    }

    public static void main(String[] args){
         File user = new File(System.getProperty("user.dir"));
         File[] f =  user.listFiles(new FilenameFilter(){

             public boolean accept(File arg0, String arg1) {
               return arg1.endsWith("zip");
             }
             
         });
         for(int i=0; i<f.length; i++){
             ConvertToByIndiv ctbi = new ConvertToByIndiv(f[i]);
             //ctbi.run();
             
         }
    }
}
