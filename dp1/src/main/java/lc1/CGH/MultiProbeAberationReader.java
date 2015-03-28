package lc1.CGH;


/** designed for the file MultiProbeByIndividual.txt*/

public class MultiProbeAberationReader extends AbstractAberatiionReader {
 //MultiprobeAberationReader
    //Array Barcode AbsPol  Cy3 Cy5 DLRSpread   HybDate Name    Polarity    Popn    chr Start (bp)  End (bp)    size (bp)   size +60bp  >1kb        start identical end identical               #probes Mean LR isNested?   Height above parent Score   pVal    probe start (bp)    probe end (bp)  probe before start (bp) probe after end (bp)    IsAluIRsaISNP

    //SingleProbeIntervals
    //Array Barcode AbsPol  Cy3 Cy5 DLRSpread   HybDate Name    Polarity    Popn    chr Start (bp)  End (bp)    size (bp)   #probes Mean LR isNested?   Height above parent Score   pVal    probe start (bp)    probe end (bp)  probe before start (bp) probe after end (bp)    IsAluIRsaISNP
    
    
 static String[] cols1 = new String[] {"Name", "chr", "Start", "End", "Mean", "#probes"};
public String[] cols;
 MultiProbeAberationReader(String[] cols, long length, String name){
  super(length, name);
     this.cols = cols;
 }
 public  MultiProbeAberationReader(long length, String name){
     this(length, cols1, name);
 }
public  MultiProbeAberationReader(long length, String[] cols, String name){
   super(length, name);
   this.cols = cols;
 }
 
@Override
public String[] getCols() {
   return cols;
}
   
public String getName(String[] str){
    return col[0]>=0 ? str[col[0]] : "";
 }
 public String getChr(String[] str){
     return str[col[1]];
 }
 public String getStart(String[] str){
     return str[col[2]];
 }
 public String getEnd(String[] str){
     return str[col[3]];
 }
 public int getNoProbes(String[] str){
     if(col[5]>=0)
     return Integer.parseInt(str[col[5]]);
     else return 0;
 }  
 public double getNoCopy(String[] str){
     if(col[4]<0) return 0;
     else return Double.parseDouble(str[col[4]]);
             
  }
 
 
}
