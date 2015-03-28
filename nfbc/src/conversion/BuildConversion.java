package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipFile;

public class BuildConversion {
    
    public static void main(String[] args){
        try{
         makeBuildFile(new File(System.getProperty("user.dir")),"build36");
           // postProcess("build.txt", "build35.txt");
          //  postProcess("liftover.txt", "build36.txt");
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }

    public static void postProcess(String inS, String outS) throws Exception{
        File user = new File(System.getProperty("user.dir"));
        File in = new File(user, inS);
        File out = new File(user, outS);
        BufferedReader br = new BufferedReader(new FileReader(in));
        String st = "";
        PrintWriter res35 = new PrintWriter(new BufferedWriter(new FileWriter(out)));
        while((st = br.readLine())!=null){
            String[] str = st.split("\t");
          
              
               int mid = getMid(str, 1, 2);//str[1],str[2]);
               
                res35.println(str[0]+"\t"+mid+"\t"+str[3]+"\t"+str[4]+"\t"+str[5]);
            
         }
        res35.close();
        br.close();
    }
    
    private static int getMid(String[] str,int i1, int i2) {
        int start = Integer.parseInt(str[i1]);
        int end = Integer.parseInt(str[i2]);
        return (int) Math.round(((double)start+(double)end)/2.0);
    }
static String[] nonnumeric = new String[] {"X", "Y","M"};
static boolean matches(String st){
	for(int i=0; i<nonnumeric.length; i++){
		if(st.indexOf(nonnumeric[i])>=0) return true;
	}
	return false;
}
    public static File makeBuildFile(File user, String build) throws Exception{
    	File buildF = new File(user, build+".txt.gz");
        PrintWriter res = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(buildF)))));
       // File user = new File(System.getProperty("user.dir"));
        File[] f = user.listFiles(new FilenameFilter(){

            public boolean accept(File arg0, String arg1) {
              return arg1.endsWith("zip");
            }
            
        });
        Arrays.sort(f, new Comparator<File>(){

            public int compare(File arg0, File arg1) {
               String name1 = arg0.getName();
               String name2 = arg1.getName();
               name1 = name1.substring(0,name1.indexOf(".zip"));
               name2 = name2.substring(0,name2.indexOf(".zip"));
               if(name1.equals(name2)) return 0;
               else if(matches(name1) || matches(name2)){
                   return name1.compareTo(name2);
               }
               else{
            	   try{
                   Integer i1 = Integer.parseInt(name1);
                   Integer i2 = Integer.parseInt(name2);
                   return i1.compareTo(i2);
            	   }catch(Exception exc){
            		   //exc.printStackTrace();
            		   return name1.compareTo(name2);
            	   }
               }
            }
            
        });
        for(int i=0; i<f.length; i++){
        	String nme = "chr"+getChrom(f[i].getName());
        	nme = nme.substring(0,nme.indexOf(".zip"));
        	System.err.println(f[i].getName());
        	String chrom = getChrom(f[i].getName());
        	List<String> l = Compressor.getIndiv(new ZipFile(f[i]), "SNPS");
        	if(l==null){
        		l = Compressor.getIndiv(new ZipFile(f[i]), chrom.split("\\.")[0]+"/SNPS");
        	}
           for(Iterator<String> it =l.iterator(); it.hasNext();){
              // res.println(it.next());
               String[] str = it.next().split("\\s+");
              // if(str[1].equals("NA")) continue;
               String start = str[1].trim();
            	   //Integer.parseInt(str[1].trim());
               String end = str[2];
//               int end = str[2].equals("NA") ? start+20 : Integer.parseInt(str[2].trim());
               /*if(end<start){
                   System.err.println("switched  at "+Arrays.asList(str));
                   int tmp = start;
                   start = end;
                   end = tmp;
               }*/
               res.print(nme+"\t"+start+"\t"+end);
               for(int j=3; j<str.length; j++){
                   res.print("\t"+str[j]);
               }
               res.println();
               
           }
        }
        res.close();
        return buildF;
    }

	private static String getChrom(String name) {
		String st = name;
		if(name.startsWith("chr")) st =  name.substring(3);
		if(name.endsWith("b")) return st.substring(0,st.length()-2);
		return st;
	}
}
