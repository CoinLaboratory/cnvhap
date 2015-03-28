package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;
import java.util.zip.ZipFile;

public class CreatePhenoTables {

    
    public static void main(String[] args){
        try{
            hapmapTable("Case");
        }catch(Exception exc){
            exc.printStackTrace();
        }
      
    }
    
    public static void hapmapTable(String hap) throws Exception{
        ZipFile zf = new ZipFile(new File("1.zip"));
        ZipFile zf1 = new ZipFile(new File("2.zip"));
        List<String> samples = Compressor.getIndiv(zf, "Samples");
        List<String> samples1 = Compressor.getIndiv(zf, "Samples");
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("Samples.txt")));
        pw.println("id\tT2Diabetes");
        for(int i=0; i<samples.size(); i++){
            if(!samples.get(i).equals(samples1.get(i))) throw new RuntimeException("!!");
            pw.println(samples.get(i).split("\t")[0]+"\t"+hap);
        }
        pw.close();
    }
    
    
     //Map<String, String[]> phens = new HashMap<String, String[]>();
     
     
}
