package lc1.CGH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

public class LocCompare {
    public String toString(){
        return name1+" "+name2+" "+Arrays.asList(oneVtwo)+" "+Arrays.asList(ratios(oneVtwo))+" "+Arrays.asList(twoVone)+" "+Arrays.asList(ratios(twoVone));
    }
    private Double[] ratios(Integer[] r) {
        return new Double[] {(double)r[0]/(double)r[1], (double)r[2]/(double)r[3]};
    }
    Integer[] oneVtwo;
    Integer[] twoVone;
    int[] noProbes;
    String name1, name2;
    PrintWriter out1;
    PrintWriter out2;
    static int overlp = 2;
    LocCompare(Locreader loc1, Locreader loc2, File logDir) throws Exception{
        this.name1 = loc1.name;
        this.name2 = loc2.name;
        PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(logDir, name1+"v"+name2))));
        PrintWriter out2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(logDir, name2+"v"+name1))));
        oneVtwo = loc1.detected(loc2,overlp, out1); //of loc2, how many detected by loc1
        twoVone = loc2.detected(loc1,overlp, out2); //of loc1 how many detected by loc2
        noProbes = new int[] {loc1.probes.size(), loc2.probes.size()};
        
        System.err.println(this);
        System.err.println(noProbes[0]+" "+noProbes[1]);
        out1.flush();
        out2.flush();
        out1.close();
        out2.close();
    }
}
