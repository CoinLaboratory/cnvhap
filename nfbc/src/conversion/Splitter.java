package conversion;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


public class Splitter {
    Map<String, PrintWriter> m = new HashMap<String, PrintWriter>();
    BufferedReader[] br;
    static final Comparator<String[]> comp = new Comparator<String[]>(){
        final int i =3;
        public int compare(String[] o1, String[] o2) {
            Integer i1 = Integer.parseInt(o1[i]);
            Integer i2 = Integer.parseInt(o2[i]);
           return i1.compareTo(i2);
        }

     
        
    };
    
    public static void main(String[] args){
        try{
            File user = new File(System.getProperty("user.dir"));
            File[] files = user.listFiles(new FileFilter(){

				public boolean accept(File arg0) {
					return arg0.getName().startsWith("All");
				}
            	
            });
            Splitter sp = new Splitter(files, true);
            sp.run();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    
    static int readFirst = 9;
    static String identifier = "1M.txt";
    Splitter(File[] f, boolean hasHeader) throws Exception{
    	br = new BufferedReader[f.length];	
    	this.inclList =null;// this.getInclList(f[0].getParentFile());
    	PrintWriter excl = new PrintWriter(new BufferedWriter(new FileWriter(new File(f[0].getParentFile(), "excl.txt"))));
    	PrintWriter inclp = new PrintWriter(new BufferedWriter(new FileWriter(new File(f[0].getParentFile(), "incl.txt"))));
    	for(int i=0; i<f.length; i++){
    		if(f[i].getName().endsWith("zip")){
    			ZipFile zf = new ZipFile(f[i]);
    			ZipEntry ent = zf.entries().nextElement();
    			br [i]= new BufferedReader(new InputStreamReader(zf.getInputStream(ent)));
    		}
    		else{
	        br[i] = 
	        	new BufferedReader(new FileReader(	f[i]   ));
	        		//"Affx_20060707fs1_SNP.txt" 
	           
    		}
    		 for(int ik=0; ik<readFirst; ik++){
                 System.err.println( br[i].readLine());
              }
    		if(hasHeader){
    		    header = new String[f.length][];
    		    incl = new boolean[f.length][];
    	        this.header[i] = br[i].readLine().split("\t");
    	        this.incl[i] = new boolean[header[i].length];
    	        Arrays.fill(incl[i], false);
    	        for(int j=0; j<header[i].length; j++){
    		        if(i==0){
    	    			for(int k=0; k<initTags.length; k++){
    	    				if(header[i][j].startsWith(initTags[k])){
    	    					incl[i][j] = true;
    	    				}
    	    			}
    	    		}
    	    		for(int k=0; k<tags.length; k++){
    	    			if(header[i][j].endsWith(tags[k])){
    	    				String id = header[i][j].split("#")[0];
    	    				if(inclList==null || this.inclList.contains(id)){
    	    					incl[i][j] = true;
    	    					inclp.println(id+"\t"+header[i][j]);
    	    				}
    	    				else{
    	    					excl.println("excluded "+header[i][j]);
    	    				}
    	    			}
    	    		}
    	        }
    		}
    		
    	}
    	excl.close();
    	inclp.close();
    	
    }
  
    public PrintWriter getPw(String st) throws Exception{
        PrintWriter pw = m.get(st);
        if(pw==null){
            m.put(st, pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(st+"_"+identifier)))));
            if(header!=null){
                List<String> l = new ArrayList<String>();
                for(int i=0; i<header.length; i++){
                	getString(header[i], i, l);
                }
                print(pw, l);
            }
        //    pw.print(getString(header, true));
        }
        return pw;
        
    }
    private Set<String> getInclList(File dir) throws Exception{
    	Set s  = new HashSet<String>();
    	{
	    	BufferedReader br = new BufferedReader(new FileReader(new File(dir, "list300k_31Aug2007.txt")));
	    	String st = "";
	    	while((st = br.readLine())!=null){
	    		s.add(st.split("\\t")[1]);
	    	}
    	}
    	{
	    	BufferedReader br = new BufferedReader(new FileReader(new File(dir, "samples_to_remove_for_300k_reanalysis.txt")));
	    	String st = "";
	    	while((st = br.readLine())!=null){
	    		if(st.length()==0) continue;
	    		s.remove(st.split("\\s+")[1]);
	    	}
    	}
    	return s;
    }
    private void print(PrintWriter pw, List<String> l) {
    	int len = l.size();
    //	System.err.println("size is "+len);
		for(int i=0; i<len-1; i++){
			pw.print(l.get(i)+"\t");
		}
		pw.println(l.get(l.size()-1));
	}
	int index = 2;
    int matchIndices = 1;
    String[] initTags = new String[] {"Sample ID", "SNP Name", "Index", "Name", "Address", "Chr", "Position"};
    String[] tags = new String[] {"GType", "B Allele Freq","Log R Ratio"};
    String[][] header;
    boolean[][]incl;
    Collection<String> inclList ;
    public void run() throws Exception{
        String st = "";
        while(true){
        	String st0 = br[0].readLine();
        	if(st0==null) break;
        	String[] str0 = st0.split("\\t");
	        PrintWriter pw = getPw(str0[index]);
	        List<String> res = new ArrayList<String>();
	        getString(str0, 0, res);
        	for(int i=1;i<br.length; i++){
        		st = br[i].readLine();
	            String[] str = st.split("\\s+");
	            if(!str[matchIndices].equals(str0[matchIndices])) throw new RuntimeException(" !! "+str[matchIndices]+" "+str0[matchIndices]);
	        	getString(str, i, res);
	         //   PrintWriter pw = getPw(str[index]);
	          
	        	}
        	 print(pw, res);
        	 pw.flush();
        }
        for(int i=0; i<br.length; i++){
        	br[i].close();
        }
    }
    public void  getString(String[] str, int i, List<String> res){
    	for(int j=0; j<str.length; j++){
    		if(incl==null || incl[i][j]){
    			res.add(str[j]);
    		}
    		
    	}
    }
}
