package conversion;
import java.io.File;
import java.io.FileFilter;
import java.io.StringWriter;
import java.util.Arrays;

import exec.Exec;


public class Sorter {
    
    
    public static void main(String[] args){
       try{
        File user = new File(System.getProperty("user.dir"));
        File[] f = user.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
               return pathname.getName().endsWith("1M.txt");
            }
            
        });
        for(int i=0; i<f.length; i++){
            Sorter.sort(f[i]);
            if(true) System.exit(0);
        }
       }catch(Exception exc){
           exc.printStackTrace();
       }
    }
    File f;
    final int index =4;
   static String bin = "";//"c:/cygwin/bin/";
  static  int[] sort = new int[] {1,4};
    public static void sort(File f) throws Exception{
        File f_out = new File(f.getParentFile(),f.getName().split("_")[0]+"_data1M.txt");
        if(!f_out.exists() || f_out.length()==0){
            String[] command = new String[] {bin+"sort",
                     "--k",
                     sort[0]+","+sort[1],
//                     "4,1",
                     "-n", 
                     "-o", f_out.getName(),
                     f.getName()}; 
             System.err.println(Arrays.asList(command));
             StringWriter writer = new StringWriter();
             Exec.exec(command, null, null, writer);
             System.err.println(writer.getBuffer().toString());
    }
    }
   
   
}
