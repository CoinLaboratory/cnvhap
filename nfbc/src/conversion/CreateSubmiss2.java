package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;

public class CreateSubmiss2 {
    
    public static final Options OPTIONS  = new Options(){
        {
                    this.addOption( OptionBuilder.withLongOpt( "todo" ).withDescription( "todo").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "build" ).withDescription( "build").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "rlib" ).withDescription( "rlib").withValueSeparator( ':' ).hasArgs().create());
                    this.addOption( OptionBuilder.withLongOpt( "prefix" ).withDescription( "prefix").withValueSeparator( ':' ).hasArgs().create());
                  
                    
        }
    };    
    public static int mb = 900;
    public static String time = "0:59:00" ;
    public static void main(String[] args){
        try{
            Parser parser= new PosixParser();
           
            final CommandLine params = parser.parse(OPTIONS, args, false);
            File dir1 =  new File(System.getProperty("user.dir"));
           
            final String[] prefix = params.getOptionValues("prefix");
            for(int j=0; j<prefix.length; j++){
                final int j1 = j;
            submiss = new File(dir1, "submissions/"+Arrays.asList(prefix).toString());
            submiss.getParentFile().mkdir();
            final String build = params.getOptionValue("build", "build35");
            String[] todo = params.getOptionValues("todo");
            File[] dirs = dir1.listFiles(new FileFilter(){

                public boolean accept(File pathname) {
                  return pathname.getName().startsWith(prefix[j1]) && pathname.isDirectory();
                }
                  
              });
            for(int i=0; i<dirs.length; i++){
                File dir = dirs[i];
            for(int k=0; k<todo.length; k++){
                print(build,  todo[k], dirs[i]);
            }
            }
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    static File submiss;
    public static void  print( String build, String todo, File dir) throws Exception{
        String absPath = dir.getAbsolutePath();
        String postData =absPath.substring(absPath.lastIndexOf("data")+5, absPath.indexOf(dir.getName())-1);
        String code = dir.getName()+"_"+todo;
      submiss.mkdir();
     // postData.replaceAll("\", "\\/");
        File out = new File(submiss,code +".sub.txt"); 
        String prefix = dir.getName();
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
       pw.println("#!/bin/sh");
       pw.println("#PBS -l walltime="+time);
       pw.println("#PBS -l mem="+mb+"mb");
       pw.println("#PBS -l ncpus=1");
       pw.println("module load R/2.4.1");
       pw.println("module load java");
       pw.println("mkdir "+dir.getName());
       pw.println("cp $WORK/data/"+postData+"/"+dir.getName()+"/"+todo+".zip "+dir.getName());
       pw.println("cp $WORK/data/"+postData+"/"+dir.getName()+"/"+build+"* "+dir.getName());
       pw.println("cp $WORK/data/"+postData+"/"+build+"* "+dir.getName());
       pw.println("java -Xmx"+(mb-100)+"m -server -cp /home/lcoin/nfbc.jar conversion.ConvertToWide --build "+build+" --todo "+todo+
                       " --prefix "+prefix+" --rlib /home/lcoin/R --mode 0");
       pw.println("java -Xmx"+(mb-100)+"m -server -cp /home/lcoin/nfbc.jar conversion.ConvertToWide --build "+build+" --todo "+todo+
               " --prefix "+prefix+" --rlib /home/lcoin/R --mode 1 --runR false");
       pw.println("R --vanilla < gcplot_tmp.R");
       pw.println("java -Xmx"+(mb-100)+"m -server -cp /home/lcoin/nfbc.jar conversion.ConvertToWide --build "+build+" --todo "+todo+
               " --prefix "+prefix+" --rlib /home/lcoin/R --mode 2");
       if(todo.equals("1")){
           pw.println("cp $WORK/data/"+postData+"/"+dir.getName()+"/Samples.txt modified/"+dir.getName());
       }
       pw.println("zip -r "+code+".zip modified/ pdf/");
       pw.println("ls -R");
       pw.println("cp "+code+".zip $WORK/data/"+postData);
       pw.close();
       
   }
}
