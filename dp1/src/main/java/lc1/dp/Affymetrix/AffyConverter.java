package lc1.dp.Affymetrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;

public class AffyConverter {
    
    
    public static void main(String[] args){
        try{
        AffyConverter conv = new AffyConverter(new File("calls.chr22.merged.txt"));
        conv.write();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    final DataCollection data;
    final int noIndiv = 62;
    final int noSnps = 761;
    final String[] names = new String[] {"A", "B", "N"};//, "AAA", "AAB","ABB", "BBB"};
    public  AffyConverter(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
        String st = "";
        final PIGData[] store; 
        store = new PIGData[noIndiv];
        for(int j=0; j<noIndiv; j++){
            store[j] =SimpleScorableObject.make(j+"", noSnps,
                    Emiss.getEmissionStateSpace(1),(short)-1);
        }
            for(int j=0; j<noIndiv; j++){
                for(int i=0; i<noSnps; i++){
                    String line = br.readLine().trim();
                    System.err.println(line);
                    Comparable ma, ma2; 
                    if(line.equals("A")){
                        ma = Emiss.a(); 
                        ma2 = Emiss.a();
                    }
                    else if(line.equals("B")){
                        ma = Emiss.b();
                        ma2 = Emiss.b();
                    }
                    else if (line.equals("AB")){
                        ma = Emiss.a();
                        ma2 = Emiss.b();
                    }
                    else if(line.startsWith("N")){
                        ma = Emiss.N();
                        ma2 = Emiss.N();
                    }
                    else throw new RuntimeException("!!");
                    ComparableArray comp = ComparableArray.make(ma, ma2);
                  //  System.err.println(i+" "+j+" "+k);
                    store[j].addPoint(i,comp);
            }
        }
            if((st =br.readLine())!=null) throw new RuntimeException(st);
        data    = new SimpleDataCollection(Arrays.asList(store));
        br.close();
    }
    
    public void write() throws Exception{
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("affy_calls.txt"))));
        
     //  data.writeFastphase( pw,  false);
       pw.close();
    }
}
