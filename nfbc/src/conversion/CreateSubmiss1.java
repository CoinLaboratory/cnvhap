package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

public class CreateSubmiss1 {
    public static int mb = 1800;
    public static String time = "3:59:00" ;
    public static void main(String[] args){
        try{
        for(int i=0; i<100; i++){
       print(i);
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    public static void  print( int i) throws Exception{
        File out = new File(i+".sub.txt"); 
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
       pw.println("#!/bin/sh");
       pw.println("#PBS -l walltime="+time);
       pw.println("#PBS -l mem="+mb+"mb");
       pw.println("#PBS -l ncpus=1");
       pw.println("module load R/2.4.1");
       pw.println("mkdir perm"+i);
       
       pw.println("echo 'perms = "+i+"' > analysis.txt");
       pw.println("cat $WORK/NFBC/analysis.txt >> analysis.txt");
       pw.println("cp -rf $WORK/NFBC/toCopy/* $TMPDIR");
       pw.println("R --vanilla < analysis.txt");
       pw.println("tar -czf perm"+i);
       pw.println("mv -f perm* $WORK/NFBC/");
       pw.close();
       
   }
}
