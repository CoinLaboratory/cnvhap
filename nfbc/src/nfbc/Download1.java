package nfbc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JOptionPane;
import javax.swing.JPasswordField;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;


public class Download1 {
    
    public static final Options OPTIONS  = new Options(){
        {
            Field[] f = Download1.class.getFields();
            for(int i=0; i<f.length; i++){
                if(Modifier.isStatic(f[i].getModifiers())){
                    this.addOption( OptionBuilder.withLongOpt( f[i].getName() ).withDescription( f[i].getName()).withValueSeparator( ':' ).hasArgs().create());
                }
                else{
                    System.err.println("excluded "+f[i]);
                }
            }
        }
    };
    
    public static String user;
    public static String password;
    public static String  experiment;
    
    public static void setOptions(String[] args){
            try{
            Parser parser= new PosixParser();
            final CommandLine params = parser.parse(OPTIONS, args, false);
            Field[] f = Download1.class.getFields();
            for(int i=0; i<f.length; i++){
                if(Modifier.isStatic(f[i].getModifiers()) && !f[i].getName().equals("OPTIONS")){
                        try{    
                            //String[] val = params.getOptionValues(f[i].getName());
                            f[i].set(null, params.getOptionValue(f[i].getName()));
                        }catch(Exception exc){
                            exc.printStackTrace();
                        }
                }
               
            }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    };
    
    File[] f1;
    File dir;
    
    static FileFilter assoc = new FileFilter(){

        public boolean accept(File arg0) {
          return  arg0.getName().endsWith(".dat");
        }
        
    };
    static FileFilter cover = new FileFilter(){

        public boolean accept(File arg0) {
          return  arg0.getName().startsWith("cover");
        }
        
    };
    
    public static final String base =  "/work/lcoin/diab_results/";
   
    public static void ftp() {
        boolean error = false;
        FTPDownload ftp=null;
        try{
           
            ftp = new FTPDownload("login.cx1.hpc.imperial.ac.uk", user,password);
               ftp.getFiles2(base,
            		   experiment);
            ftp.disconnect();
        }catch(Exception exc){
            exc.printStackTrace();
            if(ftp!=null)  ftp.disconnect();
        }
        
    }
    public static void main(String[] args){
        try{
            if(true){
                setOptions(args);
                ftp();
            }
           
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   File outputDir;
   
 
    Download1(File f) throws Exception{
        dir = f;
    
    }
  
}
