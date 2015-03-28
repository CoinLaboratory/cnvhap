package lc1.assoc;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LikelihoodDataCollection;
import lc1.util.Constants;
import lc1.util.Executor;

public class CalcAssociation {
	
	
	
	   static  Map<String, List<String>> sex = new HashMap<String,List< String>>();
    public static void main(final String[] args){
        try{
            final int start = 0;
            int[] mid = new int[] {0,Integer.MAX_VALUE};
           
           final  String[] tag_pheno = new String[] {"chisq_state", "armitage_state","odds_state","regress", "regrP",
                    "cases_state",
                    "controls_state"};
            final String[] tag_report =   new String[] {"snpid","loc", "index",   "maf"};
        Constants.plot = 2;
        Constants.numThreads = 4;
       // Constants.topBottom = false;
        String col = "1";
        final String sta = args.length>0 ? args[0] : "";
        File user = new File(System.getProperty("user.dir"));
        File sexF = new File(user, "sex.txt");
        while(!sexF.exists()){
     	   sexF = new File(sexF.getParentFile().getParentFile(), "sex.txt");
        }
        BufferedReader br = new BufferedReader(new FileReader(sexF));
        String st = br.readLine();
       String[] sexP = new String[] {"both"};//,"male", "female" };
       for(int l=0; l<sexP.length; l++){
     	 sex.put(sexP[l], new ArrayList<String>());
       }
        while((st = br.readLine())!=null){
     	   String[] str1 = st.split("\\s+");
     	   List<String> l = sex.get(str1[1]);
     	   if(l!=null)l.add( str1[0]);
     	   sex.get("both").add(str1[0]);
        }
      
        
        File[] fs =  user.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
             boolean res =  pathname.getName().endsWith("zip");
             boolean res1 =  pathname.getName().startsWith(sta);
            //System.err.println(pathname.getName()+" "+res+" "+res1+" "+args[0]);
             return res && res1;
            }
            
        });
        for(int i=0; i<fs.length; i++){
        	try{
               Constants.mid = 
            	   new String[][] {
            	   new String[] {fs[i].getName().split("_")[0], "0", ""+(Integer.MAX_VALUE -1)}
               };
                ZipFile zf = new ZipFile(fs[i]);
                Enumeration en = zf.entries();
                en.nextElement();
                System.err.println(((ZipEntry)en.nextElement()).getName());
                String nme1 = fs[i].getName();
                nme1 = nme1.substring(0,nme1.indexOf(".zip"));
                String nme = nme1+"/phased2_"+col+".txt";
                System.err.println(nme);
                String resEntry = findEntry(zf, fs[i].getName().split("_")[0]+".zip");
                File result=null;
                File dirOut = new File(user, nme1);
                if(resEntry!=null){
                   
                    dirOut.mkdir();
                    result = expand(dirOut, zf,  resEntry);
                }
                if(result==null){
                	throw new RuntimeException("null result for "+fs[i].getName());
                }
                final DataCollection dc  = new LikelihoodDataCollection(result,(short) 0, 2, new int[][]{mid},  null, null);
                final  int end  = dc.loc.size();
                ZipEntry ent1 = zf.getEntry(nme1+"/phased1_"+col+".txt") ;
                BufferedReader br1 =new BufferedReader(new InputStreamReader( zf.getInputStream(ent1)));
                String str = br1.readLine();
             //   DataCollection.scoreRegression = true;
                for(int l=0; l<sexP.length; l++){
                {
                	String out_name =  "assoc.txt";
                	if(!sexP[l].equals("both")){
                		out_name="assoc_"+sexP[l]+".txt";
                	}
                	DataCollection dc1 = l==0 ? dc : (DataCollection) dc.clone();
                	if(l>0) dc1.dropIndiv(sex.get(sexP[l]).toArray(new String[0]));
                    PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(dirOut, out_name))));
               //     dc1.printHapMapFormat(pw_hap1, null, true, 
                  //       tag_report,
                  //          tag_pheno,
                           //  : new String[] {"snpid", "loc"}, 
                 //           new String[] {}, "%7s", start, end, new int[] {0,1,2,3}, new int[] {0}, true);
                    pw_hap1.close();
                		
                }
                }
                int numRand =0;
                for(int ik=0; ik<numRand; ik++){
                    PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(dirOut, "assoc"+ik+".txt"))));
                    dc.randomizePhenotypes();
                   //     dc.printHapMapFormat(pw_hap1, null, true, 
                   //     tag_report,
                     //      tag_pheno,
                          //  : new String[] {"snpid", "loc"}, 
                      //     new String[] {}, "%7s", start, end, new int[] {0,1,2,3}, new int[] {0}, true);
                    pw_hap1.close();
                }
        	}catch(Exception exc){
        		exc.printStackTrace();
        	}
             //   if(true) System.exit(0);
               
        }
        Executor.shutdown();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    public static File expand(File dir, ZipFile zf,  String resEntry) throws Exception {
        // TODO Auto-generated method stub
        String[] f = resEntry.split("/");
        File par = dir;
        for(int i=1; i<f.length-1; i++){
            par = new File(par, f[i]);
            par.mkdir();
        }
        String name = f[f.length-1];
        File outF = new File(par, name);
        if(outF.exists() && outF.length()>0) return outF;
        BufferedInputStream in = new BufferedInputStream(zf.getInputStream(zf.getEntry(resEntry)));
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outF));
        byte[] buf = new byte[1024];
        int len = 0;
        int reallen = 0;

        //System.out.println(file+":"+getLocalPath()+outfile);
        while(true)
        {
            len = in.read(buf);

            //System.out.print(".");
            if(len == StreamTokenizer.TT_EOF)
            {
                break;
            }

            out.write(buf, 0, len);
            reallen += len;

        }
        in.close();
        out.close();
       return  outF;
    }
    public static String findEntry(ZipFile zf, String string) {
    	String nme1 = "";
        for(Enumeration en = zf.entries(); en.hasMoreElements();){
           ZipEntry ent =  (ZipEntry) en.nextElement();
           String nme = ent.getName();
           if(nme.indexOf(string)>=0 && nme.length() > nme1.length() ){
        	   String[] spl = nme.split("/");
        	   if( spl.length>2)
                nme1 = nme;
           }
        }
        return nme1.length()==0 ? null: nme1;
      }
}
