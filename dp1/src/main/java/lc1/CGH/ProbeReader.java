package lc1.CGH;


/** designed for the file MultiProbeByIndividual.txt*/

public class ProbeReader extends AbstractAberatiionReader {
 //ProbeReader
   // AMADID  014068  hg17 ProbeName   ChrName Start   Stop
   // 
 String[] cols =  new String[] {"aaa", "ChrName", "Start", "Stop", "aaaa"};

 ProbeReader(String[] cols){
  super(1000000,"");
     this.cols = cols;
 }
public  ProbeReader(){
   super(1000000,"");
 }
 
@Override
public String[] getCols() {
   return cols;
}
   
public String getName(String[] str){
    return  "";
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
     return -1;
 }  
 @Override
 public boolean exclude(String[] str){
     return false;
 }
 @Override
 public double getNoCopy(String[] str){
     return 0;
             
  }
 
 
}
