package lc1.dp.illumina;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class SplitIllumina {
    public static void main(String[] args){
        try{
            int ind = 3;
        BufferedReader br = new BufferedReader(new FileReader(new File("French_raw.txt")));
        String header = br.readLine();
        String st = br.readLine();
         while(st!=null){
            String chrom = st.split("\\s+")[ind];
            File dir = new File(chrom);
            dir.mkdir();
            File f = new File(dir, "data.txt");
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
          pw.println(header);
            while(st!=null && chrom.equals(st.split("\\s+")[ind])){
                pw.println(st);
                st = br.readLine();
            }
            pw.close();
           // break;
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
}
