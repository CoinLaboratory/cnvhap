package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class FixStartEndBuild {

    public static void main(String[] args){
        try{
            File in = new File("build35.txt");
            File out = new File("build35.out");
            BufferedReader br = new BufferedReader(new FileReader(in));
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
            String st = "";
            while((st = br.readLine())!=null){
                String[] str = st.split("\\s+");
               int start = Integer.parseInt(str[1]);
               int end = Integer.parseInt(str[2]);
               if(end<start){
                   String tmp = str[1];
                   str[1] = str[2];
                   str[2] = tmp;
               }
               pw.print(str[0]);
               for(int i=1; i<str.length; i++){
                   pw.print("\t"+str[i]);
               }
               pw.println();
            }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
}
