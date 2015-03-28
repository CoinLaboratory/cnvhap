package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import lc1.dp.data.representation.Emiss;
import lc1.dp.states.EmissionState;

public class FindAncestral {
    
   
   static List<String> data = new ArrayList<String>();
    static File ill = new File("illumina_snps_v1.txt");
    static File ancestral = new File("ancestral_state.txt");
   // static File allele = new File("Allele.bcp");
   
    
   {
        File user = new File(System.getProperty("user.dir"));
        
    }
   public static void main(String [] args){
       try{
           for(int i=9; i>=1; i--){
               data.add(i+"");
           }
           data.add("X");
           for(int i=0; i<data.size(); i++){
           FindAncestral fa = new FindAncestral(data.get(i));
           fa.readIll(data.get(i));
           }
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
   // Map<String, String> alleleM;
    Map<String, String> ancestralM;
   List<Integer>loc;
   List<String>snpid;
    PrintWriter out_pw;
    EmissionState maf;
   // SimpleDataCollection sdt;
    Set<String> rs_id = new HashSet<String>();
    FindAncestral(String chr) throws Exception{
    
        File dir = new File(System.getProperty("user.dir"));
        File inp = 
            new File(dir,   "data");
        File out = new File(inp, chr+"_anc.txt");
        if(false){
                 /*   IlluminaRDataCollection sdt = new IlluminaRDataCollection(new File(inp, chr+"_data.txt"), (short) 0);
            //        SimpleDataCollection.readAffy(inp, chr+"_data.txt", false, sdt, 2);
                 sdt.dataL.clear();
                maf =  sdt.calculateMaf1();
                  loc = sdt.loc;
                  snpid = sdt.snpid;*/
        }
     //   alleleM = map(allele, 10, null);
        ancestralM = map(ancestral, Integer.MAX_VALUE, chr);
        out_pw = new PrintWriter(new BufferedWriter(new FileWriter(out)));
    }
    
    public String map(String st){
        if(st.equals("T")) return "A";
        else if(st.equals("G")) return "C";
        else return st;
    }
    public void readIll(String chr) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(ill)) ;
        String st = "";
       st =  br.readLine();
       int a_ind=0; int b_ind=1;
       if(snpid!=null){
         a_ind = maf.getEmissionStateSpace().get(Emiss.a());
         b_ind = maf.getEmissionStateSpace().get(Emiss.b());
       }
       Object[] toPrint = new Object[5];
       String fs = "%7s %7s %7s %7s %5.3g";
        while((st = br.readLine())!=null){
            String[] str = st.split("\\s+");
            if(str[2].equals(chr)){
                double mf=-1;;
                if(snpid!=null){
                    int i = snpid.indexOf(str[1]);
                    int loc1 = loc.get(i);
                    int loc2 = 
                    Integer.parseInt(str[3]);
                  
                    Integer best_index = maf.getFixedInteger(i);
                    if(best_index!=null){
                        if(best_index==a_ind) mf =1.0;
                        else if(best_index==b_ind) mf = 0.0;
                        else throw new RuntimeException("!!");
                    }
                    else mf = maf.getEmiss(i)[a_ind];
                }
                toPrint[0] = str[1]; toPrint[1] = str[2];toPrint[2] = str[3];
                
              //  out_pw.print(str[1]+"\t"+str[2]+"\t"+str[3]+"\t");
                if(this.ancestralM.containsKey(str[1])){
                    String anc = ancestralM.get(str[1]);
                    String a = map(str[4]);
                    String b = map(str[5]);
                    if(anc.equals(a)){
                        toPrint[3] = "A"; toPrint[4] = mf;
                      //  out_pw.print("A\t"+mf);
                       
                    }
                    else if(anc.equals(b)){
                        toPrint[3] = "B"; toPrint[4] = 1-mf;
                      //  out_pw.print("B\t"+(1-mf));
                    }
                    else{
                        toPrint[3] = "mism"; toPrint[4] = -1;
                        out_pw.print("mismatch\t-");
                    }
                    out_pw.println(String.format(fs, toPrint));
                    out_pw.flush();  
                }
                else{
                    toPrint[3] = "null"; toPrint[4] = -1;
                }
               // out_pw.println("\t"+loc1+"\t"+loc2);
               
            }
        }
        out_pw.close();
    }
    public Map<String, String> map(File f, int lim, String chr) throws Exception{
        Map<String, String> m = new HashMap<String, String>();
        BufferedReader br = new BufferedReader(new FileReader(f)) ;
        br.readLine();
        String st = "";
        for(int i=0; (st = br.readLine())!=null && i<lim; i++){
            String[] str = st.trim().split("\\s+");
            if(str[1].equals(chr)){
             //   System.err.println("mapping "+str[0]+" "+str[1]);
                m.put(str[0], map(str[3]));
            }
         // else{
          //      System.err.println("excluded "+Arrays.asList(str));
       //  }
        }
        return m;
    }
   static class RS{
       String rs_id;
       int loc;
       String chr;
       MAF MAF_ceu;
       public String toString(){
           return chr+" "+rs_id+" "+loc+" "+MAF_ceu.toString();
       }
       public boolean get(String anc){
           if(MAF_ceu.A.equals(anc)) return true;
           else if (MAF_ceu.B.equals(anc)) return false;
           else throw new RuntimeException( "!! "+anc+" "+MAF_ceu.A+" "+MAF_ceu.B+" "+this.toString());
       }
       RS(String[] str){
          
           this.rs_id = str[0];
           this.loc = Integer.parseInt(str[2]);
           this.chr = str[1].substring(3);
           String[] st2 = str[str.length-1].split(";");
           int index =0;
           for(;index < st2.length; index++){
               if(st2[index].indexOf("CEU")>=0){
                   break;
               }
           }
           try{
           MAF_ceu = new MAF(st2[index], "CEU");
           }catch(Exception exc){
               inner: for(int i=0; i<st2.length; i++){
                   try{
                       MAF_ceu = new MAF(st2[i], "CEU");
                   }catch(Exception exc1){
                       
                   }
                   break inner;
               }
               if(MAF_ceu==null){
               System.err.println("prob with "+Arrays.asList(st2));
               System.err.println(st2[0]);
               exc.printStackTrace();
               }
           }
       }
       
   }
   
   static class MAF{
       String A,B;
       double mafA;
       public String toString(){
           return "MAF "+A+" "+B+" "+mafA;
       }
       MAF(String st, String pop){
          // if(!st.startsWith("CSHL-HAPMAP")) throw new RuntimeException("!! "+st);
           String[] st1 = st.split(":");
        //   String pop1 = st1[1].split("-")[1].split("\\s+")[0];
         //  if(!pop.equals(pop1)) throw new RuntimeException("!!");
           String[] st2 = st1[2].split(",");
          String[] stA =  st2[0].split("\\s+");
          String[] stB =  st2[1].split("\\s+");
          A = stA[0];
          B = stB[0];
          mafA = Double.parseDouble(stA[1]);
       }
   }
}
