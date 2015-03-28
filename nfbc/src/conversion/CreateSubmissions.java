package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;

public class CreateSubmissions {
    
    public static void main(String[] args){
        try{
           CreateSubmissions cs = new CreateSubmissions(new File(args[0]));
           cs.run();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
//    File outdir; 
    File dir;
    File[] f;
    CreateSubmissions(File dir){
      /*  outdir =  new File(dir, "submiss");
        if(!outdir.exists()){
            outdir.mkdir();
        }
        else{
            File[] f= outdir.listFiles();
            for(int i=0; i<f.length; i++){
                f[i].delete();
            }
        }*/
       // String d = System.getProperty("user.dir");
        this.dir = dir;
        f = dir.listFiles(new FileFilter(){

            public boolean accept(File arg0) {
             return  (arg0.getName().endsWith("data1.txt.gz")) ;
            }
            
        });
    }
    
    public void run() throws Exception{
        for(int i=0; i<f.length; i++){
            printInfinium(f[i]);
        }
    }
    public static int mb = 1800;
    public static String time = "0:59:00" ;
    public void  print( int i) throws Exception{
        File out = new File(i+".sub.txt"); 
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
       pw.println("#!/bin/sh");
       pw.println("#PBS -l walltime="+time);
       pw.println("#PBS -l mem="+mb+"mb");
       pw.println("#PBS -l ncpus=1");
       pw.println("module load R");
       pw.println("mkdir perm"+i);
       
       pw.println("echo perms = "+i+" > analysis.txt");
       pw.println("cat $WORK/NFBC/analysis.txt >> analysis.txt");
       pw.println("R --vanilla < analysis.txt");
       pw.println("tar -czf perm"+i);
       

       pw.println("mv -f perm* $WORK/NFBC/");
       pw.close();
       
   }
    
    public void  print1(File f) throws Exception{
        String name = f.getName();
        name = name.substring(0, name.indexOf(".gz"));
        File out = new File(f.getAbsolutePath()+".sub.txt");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
       pw.println("#!/bin/sh");
       pw.println("#PBS -l walltime="+time);
       pw.println("#PBS -l mem="+mb+"mb");
       pw.println("#PBS -l ncpus=1");
       pw.println("module load java");
       pw.println("cp "+f.getAbsolutePath()+" $TMPDIR");
       pw.println("gunzip "+f.getName());
       pw.println("sort --k 1,2 -n -o tmp "+name);
      
       pw.println("mv -f tmp "+name);
       pw.println("java -Xmx"+(mb-100)+"m -server -cp $HOME/nfbc.jar conversion.Compressor "+name);
       pw.println("mv *zip $WORK");
       pw.close();
       
   }
    
    public void  printInfinium(File f) throws Exception{
        String name = f.getName();
        name = name.substring(0, name.indexOf(".gz"));
        File out = new File(f.getAbsolutePath()+".sub.txt");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
       pw.println("#!/bin/sh");
       pw.println("#PBS -l walltime="+time);
       pw.println("#PBS -l mem="+mb+"mb");
       pw.println("#PBS -l ncpus=1");
       pw.println("module load java");
       pw.println("cp $WORK/French/data/"+f.getName()+" $TMPDIR");
    
       pw.println("java -Xmx"+(mb-100)+"m -server -cp $HOME/nfbc.jar conversion.Compressor "+f.getName());
       pw.println("mv *zip $WORK");
       pw.close();
       
   }
}
