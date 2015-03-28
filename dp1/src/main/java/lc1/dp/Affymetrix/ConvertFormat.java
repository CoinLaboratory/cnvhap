package lc1.dp.Affymetrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;

import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.util.Constants;

/** takes file with one column in format
 * genotypes/snps/individuals
 * 
 * and returns file with many columns  with each row a snp and [genotypes individials]
 * @author lcoin
 *
 */
public class ConvertFormat {
    
    double[][][] l;
    int numGen = 15;
    
  final   int numIndiv;
  final  int numSnp;
    Integer[] pos;
    
 public static void main(String[] args){
   //  Constants.modelCNP = 6;
     Constants.offset = 0;
   String f = "L.unphased.txt";
       //"L.SU1.genos.haps.1.txt";
   String f1 = "sim.unphased.txt";
  try{ 
   ConvertFormat conv = new ConvertFormat(f, f1);
 }catch(Exception exc){
     exc.printStackTrace();
 }
 }
 
 ConvertFormat(String file, String posFSt) throws Exception{
         BufferedReader br = new BufferedReader(new FileReader(new File(file)));
         PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(new File("L.50K_Merged.CEPH.Affy.txt"))));
         String st = "";
         File posF = new File(posFSt);
         SimpleDataCollection data = new SimpleDataCollection();
         String[] names = SimpleDataCollection.readAffy(posF.getParentFile(), posF.getName(), false, data,2);
         PrintWriter delPos = new PrintWriter(new BufferedWriter(new FileWriter(new File("cnv_orig.txt"))));
        // data.printDeletedPositions(delPos);
         delPos.close();
         List<Integer> pos1=  data.loc;
         pos = pos1.toArray(new Integer[0]);
         numSnp = pos.length;
       //  List<String>  names = data.getNames();
         numIndiv = names.length;
         System.err.println("num indiv "+numIndiv+" "+numSnp+" "+numGen+" "+(numIndiv*numSnp*numGen));
         EmissionStateSpace emiss =Emiss.getEmissionStateSpace(1);
         l = new double[numIndiv][numSnp][numGen];
        // pos = new int[numSnp];
         for(int i=0; i<numIndiv; i++){
            PIGData dat =  data.get(names[i]);
           
             for(int j=0; j<numSnp; j++){
                 int rightIndex = emiss.getGenotype(dat.getElement(j));
                 for(int k=0; k<numGen; k++){
                     l[i][j][k] = Double.parseDouble(br.readLine());
                 }
                 int maxIndex = Constants.getMin(l[i][j]);
                 //if ( l[i][j][rightIndex]-l[i][j][maxIndex] >3){
                 //    throw new RuntimeException("!! "+l[i][j][rightIndex]+" "+l[i][j][maxIndex]);
                     
               //  }
                 //else System.err.println("ok");
             }
         }
         String str1 = br.readLine();
         if(str1!=null) throw new RuntimeException("not null "+str1);
         br.close();
     
       
         /*BufferedReader br1 = new BufferedReader(new FileReader(posF));
         for(int j=0; j<numSnp; j++){
             String[] st1 = br1.readLine().trim().split("\\s+");
             pos[j] = Integer.parseInt(st1[2]);
         }*/
         out.print("chrom  pos  ");
         for(int i=0; i<numIndiv; i++){
             String nme = names[i];
             for(int k=0; k<numGen; k++){
                 out.print(nme+"_"+k+" ");
             }
         }
         out.println();
         for(int j=0; j<numSnp; j++){
             out.print("T\t"+pos[j]+"\t");
             for(int i=0; i<numIndiv; i++){
                 for(int k=0; k<numGen; k++){
                     out.print(l[i][j][k]+" ");
                 }
             }
             out.println();
         }
         
       out.close();
         
         PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(new File("calls.50K_Merged.CEPH.Affy.txt"))));
      PrintWriter indiv = new PrintWriter(new BufferedWriter(new FileWriter(new File("indiv.txt"))));
       for(int i=0; i<numIndiv; i++){
           indiv.println(names[i]);
        }
         indiv.close();
       /*  SimpleDataCollection data = DataCollection.readMarchiniFormat(new File[] {new File("genos.haps.1")}, posF);*/
         out1.print("chrom  pos  ");
         for(int i=0; i<numIndiv; i++){
             String nme = names[i];
                 out1.print(nme+" ");
         }
         out1.println();
         for(int j=0; j<numSnp; j++){
             out1.print("T\t"+pos[j]+"\t");
             for(int i=0; i<numIndiv; i++){
                 String key = names[i];//"INDIV"+(i+1);
                 PIGData d = data.get(key);
                ComparableArray comp = (ComparableArray) d.getElement(j);
               out1.print("\"");
                out1.print(((Emiss)comp.get(0)).toStringShort()+((Emiss)comp.get(1)).toStringShort());
                out1.print("\" ");
             }
             out1.println();
         }
             out1.close();
    
 }
}
