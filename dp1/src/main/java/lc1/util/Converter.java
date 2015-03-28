package lc1.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class Converter {
    
    
    public static void main(String[] args){
        try{
        Converter conv = new Converter(new File("lhoods.chr22.merged.txt"));
        conv.write();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    final Double[][][] store; 
    final String outpString;
    final String headerLine;
    final int noIndiv = 62;
    final int noSnps = 759;
    final String[] names = new String[] {"A", "B", "AA", "AB", "BB"};//, "AAA", "AAB","ABB", "BBB"};
    public  Converter(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        String st = "";
        StringBuffer sb = new StringBuffer();
        StringBuffer header = new StringBuffer();
        store = new Double[names.length][noIndiv][noSnps];
        for(int j=0; j<noSnps; j++){
            header.append("\""+j+ "\" ");
            if(j<noSnps-1) sb.append("%5.3f ");
        }
        headerLine = header.toString();
        sb.append("%5.3f");
        outpString = sb.toString();
        for(int i=0; i<noSnps; i++){
            for(int j=0; j<noIndiv; j++){
                for(int k=0; k<names.length; k++){
                  //  System.err.println(i+" "+j+" "+k);
                    store[k][j][i] = Double.parseDouble(br.readLine());
                }
            }
        }
        if((st =br.readLine())!=null) throw new RuntimeException(st);
        br.close();
    }
    
    public void write() throws Exception{
        for(int k=0; k<names.length; k++){
            File out = new File(names[k]+".lhood");
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
            pw.println(headerLine);
            for(int j=0; j<noIndiv; j++){
                pw.print(String.format("%-6s", new String[] {"\""+j+"\""}));
                pw.println(String.format(outpString, store[k][j]));
            }
            pw.close();
        }
    }
}
