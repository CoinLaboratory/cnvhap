package lc1.dp.Affymetrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

import lc1.dp.data.collection.LikelihoodDataCollection;
import lc1.ensj.PedigreeDataCollection;

public class MakePedigreeFile {
   public static void main(String[] args){
       try{
       File fi = new File("L.50K_Merged.CEPH.Affy.txt");
     
                BufferedReader br = new BufferedReader(new FileReader(fi));
                String[] st = br.readLine().trim().split("\\s+");
                String[] cat = LikelihoodDataCollection.cat;
                System.err.println(st.length-2+" "+cat.length);
              //  System.err.println
                if(Math.IEEEremainder(st.length-2, cat.length)!=0) throw new RuntimeException("!!");
                int noIndiv = (int) Math.floor((double) (st.length-2)/(double)cat.length);
                if(noIndiv*cat.length!=st.length-2) throw new RuntimeException("!!"+noIndiv);
                PrintWriter indiv = new PrintWriter(new BufferedWriter(new FileWriter("indiv.txt")));
                PrintWriter ped_pw = new PrintWriter(new BufferedWriter(new FileWriter("ped.txt")));
              String[] ldata = new String[noIndiv];
              for(int i=0; i<noIndiv; i++){
                  ldata[i] = "INDIV"+i;
                  indiv.println(ldata[i]);
              }
                PedigreeDataCollection ped = new PedigreeDataCollection();
                for(int i=0; i<ldata.length; i+=3){
                    String father = ldata[i];
                    String mother = ldata[i+1];
                    String child = ldata[i+2];
                    
                   ped.mother.put( child, mother);
                   ped.father.put( child,father);
              
                }
                ped.print(ped_pw);
                ped_pw.close();
                indiv.close();
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
}
