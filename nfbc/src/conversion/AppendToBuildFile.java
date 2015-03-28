package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

public class AppendToBuildFile {
public static void main(String[] args){
    try{
        File dir =  new File(System.getProperty("user.dir"));
        AppendToBuildFile atbf = new AppendToBuildFile(dir);
        atbf.run();
    }catch(Exception exc){
        exc.printStackTrace();
    }
}
    BufferedReader in;
    PrintWriter out;
    Map<String, ZipFile> m = new HashMap<String, ZipFile>();
    
    AppendToBuildFile(File dir) throws Exception{
        File[] f = dir.listFiles(new FilenameFilter(){
            public boolean accept(File arg0, String arg1) {
                return arg1.endsWith(".zip");
            }
            
        });
        m = new HashMap<String, ZipFile>();
       for(int i=0; i<f.length; i++){
           String nme = f[i].getName();
           m.put("chr"+nme.substring(0,nme.indexOf(".zip")), new ZipFile(f[i]));
       }
       in = new BufferedReader(new FileReader("build36.txt"));
    out = new PrintWriter(new FileWriter("build36_mod.txt"));
    }
    
    
    public void run() throws Exception{
        String st = "";
        while((st = in.readLine())!=null){
            String[] str = st.split("\t");
            Character[] ch = getAB(str[0],str[3]);
            if(ch[0]==null) ch[0] = '-';
            if(ch[1]==null) ch[1] = '-';
            out.println(st+"\t"+ch[0]+" "+ch[1]);
        }
        out.close();
        in.close();
    }
    
    Character[] tmp = new Character[2];
    public  Character[] getAB(String chr_st, String rs) throws Exception{
        ZipFile zf = m.get(chr_st);
        List<String> l = Compressor.getIndiv(zf, rs);
        Arrays.fill(tmp, null);
        if(rs.startsWith("cnv")) {
            Arrays.fill(tmp, '-');
            return tmp;
        }
        for(int i=0; i<l.size(); i++){
            String[] str = l.get(i).split("\t");
            char[] ab = str[0].toCharArray();
            char[] bases = str[3].toCharArray();
            for(int k=0; k<bases.length; k++){
                if(ab[k]=='A') tmp[0] = bases[k];
                else if(ab[k]=='B') tmp[1] = bases[k]; 
            }
            if(tmp[0]!=null && tmp[1]!=null) return tmp;
        }
        
        return tmp;
    }

    
}
