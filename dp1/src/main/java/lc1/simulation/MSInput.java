package lc1.simulation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.PushbackReader;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.util.Constants;

import org.biojava.utils.ProcessTools;

import pal.alignment.Alignment;
import pal.alignment.SimpleAlignment;
import pal.alignment.StrippedAlignment;
import pal.datatype.DataType;
import pal.datatype.Nucleotides;
import pal.datatype.TwoStates;
import pal.misc.IdGroup;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeUtils;
import pal.tree.Tree;

public class MSInput {
    static PrintWriter log;
    public static void main(String[] args){
    
       
        try{
        	//BigInteger bi = new BigInteger("100",2);
        	//BigInteger bi1 = new BigInteger("010",2);
        	//BigInteger bi2 = new BigInteger("001",2);
            Constants.parse(args);
         //   for(int i=0; i<1; i++){
                String filename = Constants.outputDir();
                File f = new File(filename);
                if(!f.exists()) f.mkdir();
                log =  new PrintWriter(new BufferedWriter(new FileWriter(new File(f, filename+".log_ms .txt"))));
                Field[] fields = MSInput.class.getFields();
                for(int ik=0; ik<fields.length; ik++){
                    if(Modifier.isStatic(fields[ik].getModifiers())){
                        log.println(fields[ik].getName()+"\t"+fields[ik].get(null));
                    }
                }
                log.flush();
                PrintWriter enm =  new PrintWriter(new BufferedWriter(new FileWriter(new File(f, "summaryFile.txt"))));
                MSInput ms = new MSInput(enm);
              ExtNodeMap enM = new ExtNodeMap(ms);
              ExtNodeMap enM1 = new ExtNodeMap(ms.align, ms.position, enM);
              enM.check(enM1);
              for(Iterator<BigInteger> it = enM1.keySet().iterator(); it.hasNext();){
            	  enM1.inferred.remove(it.next());
              }
            
              enM.check(enM1.inferred);
      //            Node[] ext = NodeUtils.getExternalNodes(ms.trees.get(0).getRoot());
                System.err.println("h");
                  // enm_out.println("Node map");
                 //enm.print(enm_out);
                 //enm_out.close();
           //     DataCollection d = ms.executeMS(enm);
              /*  Alignment align = ms.align;
                System.err.println("before");
                System.err.println(align.getSiteCount());
                align = trim(align);
                System.err.println("after ");
                System.err.println(align.getSiteCount());*/
              //  log.println(ms.mutationTrees.size()+" marginal trees and "+align.getSiteCount()+" seg sites");
               
             
              /*  if(!f.exists()) f.mkdir();
               
                File fi = new File(f, "sim.phased.txt");
                File fi1 = new File(f, "sim.unphased.txt");
                    (new File(f, "sim.phased.txt.loc")).delete();
                    (new File(f, "sim.unphased.txt.loc")).delete();*/
               // DataCollection d = ms.print();
               
             //   d.printHapMapFormat(fi, "sim");
//                ((SimpleDataCollection)d).transform(new double[] {0.0, 1.0}, new ArrayList(d.data.keySet()), ".");
               // d.printHapMapFormat(fi1, "sim");
           //     pw.close();
               // double[] gaps = ms.getFracGaps();
              //  ((SimpleDataCollection)d).split();
              //  double[] ld_av = ((SimpleDataCollection)d).calcLDAverage();
             //   log.println("ld average is "+ld_av[0]+" "+ld_av[1]);
             //   System.err.println("ld average is "+ld_av[0]+" "+ld_av[1]);
               // log.println("frac gaps (overall and seg sites"+gaps[0]+" "+gaps[1]);
                log.close();
           // }
           // log.println(align);
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   private static Alignment  trim(Alignment align2) {
       DataType dt = align2.getDataType();
        StrippedAlignment al = new StrippedAlignment(align2);
        for(int i=0; i<align2.getSiteCount(); i++){
            double count0=0;
            double count1 =0;
            for(int j=0; j<align2.getSequenceCount(); j++){
                int ch = dt.getState(align2.getData(j, i));
                if(ch==0 || ch==2) count0++;
                else if(ch==1 || ch==3 ) count1++;
            }
            double total = count0+count1;
            if(count0/total < maf_min || count1/total < maf_min){
                al.dropSite(i);
            }
        }
       return al;
        
    }
public static  int numSamples = 1,
                      numAlleles = 11;
static {
    Constants.probCrossOverBetweenBP = 1.1e-8;
}
   //public static int numAllelesToPrint =200;
   //0.66, 0.14
   public static double  N0 = 1e4, //10e3,  //equivalent population
                        //prob cross over between base pairs per generation 
                         probMutationPerBp = 2.2e-8, //prob of mutation per site per generation
                         length = 68012*5,//1,//9,//19,//19;//222022;//2220222;//500000;//1000000; //bp
                         probMutationPerBp_del = 2.2e-8, //prob of mutation per site per generation
                         probMutationPerBp_copy = 2.2e-8, //prob of mutation per site per generation
                         rho = 4 *N0 *Constants.probCrossOverBetweenBP*(length-1),
                         theta = 4 *N0 *probMutationPerBp*(length-1),
                         theta_del = theta /4000,
                         theta_copy =theta/4000 ,
                         deletion_length_geomean = 7000,
                         deletion_length_geo_dev = 4,
                         copy_length_geomean = 7000,
                         copy_length_geo_dev = 4,
                         density = 1000,// 6187,//1271,//838;//724.58;
                         maf_min =0.01;// 0.0227;
    

    List<Integer> end;// = new ArrayList<Integer>();
    List<Integer> startL;// = new ArrayList<Integer>();
    Double[] maf;
    Alignment align;
    Integer[] position;
    List<Tree> trees;
 
   
   public MSInput(PrintWriter pw) throws Exception{
      
       String[] st = new String[] {"ms"
       , numAlleles+" ",
                                   numSamples+" ", "-t ", theta+"",
                                   "-r ", rho+" ", length+" ", "-T"
    		   };
//       PipedInputStream pi = new PipedInputStream();
 //      PipedOutputStream po = new PipedOutputStream();
  //     po.connect(pi);
  //	 final    OutputStreamWriter osw = new OutputStreamWriter(po);
   //osw.write(str)
  //	 osw.write("");
     //  PrintWriter pw = new PrintWriter(osw);
    //   BufferedReader in = new BufferedReader(new InputStreamReader(pi));
	 // Writer err = new FileWriter(new File("tmp.out"));//new OutputStreamWriter(System.err);
		  
   Writer std_out = new StringWriter();
       StringBuffer sb = new StringBuffer();
      log.println("rho, theta = "+0+" "+theta);
      
       ProcessTools.exec(st, null, std_out,std_out);
       st[4] = theta+" ";
       for(int i=0; i<st.length; i++){
           sb.append(st[i]);
       }
       log.println("executing\n "+sb.toString());
       System.err.println("executing\n "+sb.toString());
   //  System.out.println(std_out.getBuffer().toString());
    //  System.err.println(std_err.getBuffer().toString());
       BufferedReader br = new BufferedReader(new StringReader(std_out.toString()));
       startL = new ArrayList<Integer>();
       end = new ArrayList<Integer>();
       trees = read(br, (int) length, pw, startL, end);
       br = new BufferedReader(new StringReader(std_out.toString()));
       List<Integer> sites = new ArrayList<Integer>();
      align = readAlignment(br,length,sites);
      this.position = sites.toArray(new Integer[0]);
    //   br.close();
     /*  align = conc();
       log.println("average r2 "+calculateR2());
       position = getPosition();
       log.println(" number of sites is "+align.getSiteCount());
      // thin(maf_min);
       log.println("average r2 "+calculateR2());
       int no_sites = align.getSiteCount();
       log.println(" number of sites is "+align.getSiteCount());
       int target_no_sites = (int) Math.round((double)length/density);
       if(target_no_sites < no_sites){
           List<Integer> keptSites = this.thinRandom(target_no_sites);
           this.thin(keptSites);
       }
      log.println("target no sites "+target_no_sites);
       log.println("average r2 "+calculateR2());
       log.println("new number of sites is "+align.getSiteCount()+" "+target_no_sites);*/
  
     //  return  null;
   }
   
  

private Alignment readAlignment(BufferedReader br, double length, List<Integer> posi) throws Exception{
	String st = "";
	while(!(st = br.readLine()).startsWith("positions")){}
	String[] pos = st.split(":")[1].trim().split("\\s+");
	for(int k=0; k<pos.length; k++){
		posi.add((int)Math.round(length*Double.parseDouble(pos[k])));
	}
	List<String> l = new ArrayList<String>();
	while((st = br.readLine())!=null){l.add(st);}
	Identifier[] id = new Identifier[l.size()];
	for(int k=0; k<id.length; k++){
		id[k] = new Identifier(k+"");
	}
	br.close();
	DataType dt = TwoStates.DEFAULT_INSTANCE;
	return new SimpleAlignment(id,l.toArray(new String[0]),dt);
	
}


//Map<Long, BigInteger> m = new HashMap<Long, BigInteger>();
  int segsites=-1;
   public List<Tree> read(BufferedReader br,int totlen, PrintWriter enm_out, List<Integer> startL, List<Integer> end) throws Exception{
       String st = "";
       List<Tree> tree = new ArrayList<Tree>();
      
       while(!(st= br.readLine()).startsWith("//")){}
       int start =0;
       String prevTreeSt = "";
       enm_out.println("trees");
       outer: while((st= br.readLine())!=null){
           if(st.startsWith("segsites")){
               segsites = Integer.parseInt(st.split("\\s+")[1]);
               System.err.println("seg sites is "+segsites);
               break;
           }
           String[] str = st.split("]");
           Integer pos = null;
           String treeSt = null;
           if(str.length>1){
               pos = Integer.parseInt(str[0].substring(1));
               treeSt = str[1];
           }
           else{
               pos = totlen;
               treeSt = str[0];
           }
           if(!prevTreeSt.equals(treeSt)){
               startL.add(start);
               enm_out.print(start+" "+(start+pos)+" "+treeSt);
               PushbackReader pr = new PushbackReader(new StringReader(treeSt));
               Tree tr = new pal.tree.ReadTree(pr);
               for(int i=0; i<tr.getIdCount(); i++){
               //  ****  tr.getIdentifier(i).setName(Long.parseLong(tr.getIdentifier(i).getName()));
               }
               for(int i=0; i<tr.getInternalNodeCount(); i++){
                   //****tr.getInternalNode(i).getIdentifier().setID((long)(i+tr.getExternalNodeCount()+1));
                   tr.getInternalNode(i).getIdentifier().setName(""+(i+tr.getExternalNodeCount()+1));
               }
               tree.add(tr);
               end.add(start+pos);
               prevTreeSt = treeSt;
           }
        //   else{
           //    System.err.println("identical trees");
         //  }
         
           start+=pos;
         
         
        /*  MutationTree mutations = new MutationTree(tr,  theta * (pos / length), 1,1.00001,start, start+ pos, "mut");
          // MutationTree mutations_del  = new MutationTree(tr,theta_del* (pos / length), deletion_length_geomean, deletion_length_geo_dev, start, start+pos, "del");
          // MutationTree mutations_copy  = new MutationTree(tr,theta_copy* (pos / length), copy_length_geomean, copy_length_geo_dev, start, start+pos, "copy");
           IdGroup idg = mutations.getIDG();
           mutations.simulate(mutations.mutations, idg);
           mutations_del.simulate(mutations.mutations, idg);
           mutations_copy.simulate(mutations.mutations, idg);
           mutationTrees.add(mutations);
           deletionTrees.add(mutations_del);
           copyTrees.add(mutations_copy);
          
           */
       }
     
       br.close();
     //  LocationMap lm = new LocationMap(ext, enm);
       return  tree;//lm.getDataCollection();
       /*for(int i=0; i<tree.size(); i++){
           Tree tr = tree.get(i);
           MutationTree1 mutT = new MutationTree1(tr, startL.get(i), end.get(i), "mut");
           mutT.fillMutations(enm, mutT.getIDG());
           mutT.convertDoublePosToOrder();
           mutT.simulate();
           this.mutationTrees.add(mutT);
           
         
       }*/
   }
   
 
   
}
/*
public Alignment getMergedAlignment(Alignment mut, Alignment del, Alignment copy){
    DataType dt = new Nucleotides();   
    char[][] c = new char[mut.getSequenceCount()][mut.getSiteCount()];
    for(int j=0; j<c.length; j++){
        for(int i=0; i<c[j].length; i++){
            char c1 = mut.getData(j, i);
            if(del.getData(j, i)=='1'){
                c[j][i] ='-';
            }
            else if(copy.getData(j, i)=='1'){
                if(c1=='0'){
                    c[j][i] = dt.getChar(2);
                }
                else{
                    c[j][i] = dt.getChar(3);   
                }
            }
            else{
                if(c1=='0'){
                    c[j][i] = dt.getChar(0);
                }
                else{
                    c[j][i] = dt.getChar(1);   
                }
            }
        }
    }
    return new SimpleAlignment((IdGroup)mut,c,"-",dt );
}
*/

/*
public double[] getFracGaps(){
    double cnt = 0;
    Set<Integer> gapPos = new HashSet<Integer>();
    for(int j=0; j<align.getIdCount(); j++){
        for(int i=0; i<align.getSiteCount(); i++){
            if(align.GAPS.indexOf(align.getData(j,i))>=0){
                cnt++;
                gapPos.add(i);
            }
        }
    }
    return 
    new double[] {
    cnt / (double)(align.getIdCount()*align.getSiteCount()),
    (double)gapPos.size() / (double) align.getSiteCount()
    };
}

public void addRandomGaps(double perc){
    char[][] c = new char[align.getIdCount()][align.getSiteCount()];
    for(int j=0; j<align.getIdCount(); j++){
        for(int i=0; i<align.getSiteCount(); i++){
            if(Constants.rand.nextDouble()<perc){
                c[j][i] = '?';
            }
            else{
                c[j][i] = align.getData(j,i);
            }
        }
    }
    this.align = new SimpleAlignment((IdGroup)align,c, align.GAPS, align.getDataType());
    log.println("h");
}
public void thin(List<Integer> keptSites){
    List<Double> maf_new = new ArrayList<Double>();
    List<Integer> pos_new = new ArrayList<Integer>();
    StringBuffer[] sb = new StringBuffer[align.getIdCount()];
    for(int i=0; i<sb.length; i++){
        sb[i] = new StringBuffer();
    }
    log.println("kept sites "+keptSites);
    for(int i1=0; i1<keptSites.size(); i1++){
        int i = keptSites.get(i1);
           maf_new.add(maf[i]);
           pos_new.add(position[i]);
           for(int j=0; j<sb.length; j++){
               sb[j].append(align.getData(j,i));
           }
      
   }
   this.maf = maf_new.toArray(new Double[0]);
   this.position = pos_new.toArray(new Integer[0]);
   String[] st = new String[sb.length];
   for(int i=0; i<sb.length; i++){
       st[i] = sb[i].toString();
   }
   this.align = new SimpleAlignment((IdGroup)align,st,  align.GAPS, align.getDataType() );
}
public List<Integer> thinRandom(int target){
    List<Integer> allSites = new LinkedList<Integer>();
    for(int i=0; i<align.getSiteCount(); i++){
        allSites.add(i);
    }
    List<Integer> keptSites = new ArrayList<Integer>();
    while(keptSites.size()<target){
        keptSites.add(allSites.remove(Constants.nextInt(allSites.size())));
    }
    Collections.sort(keptSites);
    return keptSites;
   
}

public void thin(double maf_thresh){
    List<Double> maf_new = new ArrayList<Double>();
    List<Integer> pos_new = new ArrayList<Integer>();
    StringBuffer[] sb = new StringBuffer[align.getIdCount()];
    for(int i=0; i<sb.length; i++){
        sb[i] = new StringBuffer();
    }
    for(int i=0; i<align.getSiteCount(); i++){
        if(this.maf[i]>maf_thresh){
            maf_new.add(maf[i]);
            pos_new.add(position[i]);
            for(int j=0; j<sb.length; j++){
                sb[j].append(align.getData(j,i));
            }
        }
       
    }
    this.maf = maf_new.toArray(new Double[0]);
    this.position = pos_new.toArray(new Integer[0]);
    String[] st = new String[sb.length];
    for(int i=0; i<sb.length; i++){
        st[i] = sb[i].toString();
    }
    this.align = new SimpleAlignment((IdGroup)align,st,  align.GAPS, align.getDataType() );
    //   for(int i=0; i<maf.l)
}*/

/*
public double calculateR2(){
       Double[] res = new Double[align.getSiteCount()];
       Double[] res1 = new Double[align.getSiteCount()];
       this.maf = this.getMaf(res, res1);
       return this.calculateAverageR(res, res1);
   }
   
   String getPrintString(int no, String st){
       StringBuffer sb = new StringBuffer();
       for(int i=0; i<no; i++){
           sb.append(st);
       }
       return sb.toString();
   }
  */ 
 /*  private static Comparable translate(char c, DataType dt) {
       int num =dt.getState(c); 
       if(num==0) return Emiss.b();
       else if(num==1) return Emiss.a();
       else if(num==2) return Emiss.bb();
       else if(num==3) return Emiss.aa();
     //  else if(num==4)
      // else if(c=='9'){
     //      if(first) return Emiss.A();
     //      else return Emiss.B();
     //  }
      else return Emiss.N();
   }
   public DataCollection print(){
       DataCollection d = new SimpleDataCollection(align.getSiteCount());
       d.loc = Arrays.asList(this.position);
       int noSites = align.getSiteCount();
           for(int i=0; i<align.getSequenceCount(); i++){
               PIGData data = SimpleScorableObject.make(i+"", align.getSiteCount(),    
            		   Emiss.getEmissionStateSpace(1),(short)-1);
             
               for(int j=0; j<noSites; j++){
                  char c = align.getData(i, j);
                 
           //        Comparable c1 = translate(c,align.getDataType());
                   
                   //     data.addPoint(j,ComparableArray.make(c1));
                 
                }
               d.data.put(data.getName(), data);
           }
            //  char c = align.getData(k,i);
            // 
                   
             //  if(i<align.getSiteCount()-1)pw.print("\t");
             //  else pw.println();
         //  }
       //}
    return d;
   }*/
  
   /*private Integer[] getPosition() {
       List<Integer> res = new ArrayList<Integer>();
       for(int i=0; i<this.mutationTrees.size(); i++){
        //   double frac = (double)(this.end.get(i) - startPos)/(double)length;
           res.addAll(Arrays.asList(this.mutationTrees.get(i).getPositions()));
       }
      
       return res.toArray(new Integer[0]);
   }*/
   /*
   public double calculateAverageR(Double[] maf, Double[] pair_maf){
       double sum=0;
       int cnt =0;
       for(int i=1; i<maf.length; i++){
           if(maf[i]!=null && maf[i-1]!=null && pair_maf[i]!=null && maf[i]!=0 && maf[i-1]!=0){
               double r = Math.pow(pair_maf[i] - maf[i]*maf[i-1],2) / 
               ((maf[i]*maf[i-1])*(1-maf[i])*(1-maf[i-1]));
               sum+=r;
               cnt++;
           }
           
       }
       return sum/(double)cnt;
   }

   public Double[] getMaf(Double[] res, Double[] res1){
       Character[] minor = new Character[res.length];
       for(int i=0; i<res.length; i++){
           int count0 =0;
           int count1 =0;
           for(int j=0; j<align.getIdCount(); j++){
               if(align.getData(j,i)=='0')count0++;
               else if(align.getData(j,i)=='1')count1++;
           }
           double total = count0+count1;
           if(count0<count1){
               minor[i] = '0';
               res[i] =(double) count0 / (double)total;
           }
           else {
               minor[i]  = '1';
               res[i] = (double)count1/(double)total;
           }
           int countMaf =0;
           int totalMaf =0;
           if(i>0){
               for(int j=0; j<align.getIdCount(); j++){
                   char ch = align.getData(j,i);
                   char ch1 = align.getData(j,i-1);
                   if(align.GAPS.indexOf(ch)>=0 || align.GAPS.indexOf(ch1)>=0) continue;
                   if(ch==minor[i] && ch1 == minor[i-1])countMaf++;
                   totalMaf++;
               }
               res1[i] = (double)countMaf/(double)totalMaf;
           }
           
           
       }
       return res;
   }
   
/*public Alignment conc(List<MutationTree1> trees, int cat){
       Alignment[] list1 = new Alignment[trees.size()];
       for(int i=0; i<list1.length; i++){
           list1[i] = trees.get(i).align[cat];
       }
       return new ConcatenatedAlignment(list1);
   }
List<MutationTree1> mutationTrees = new ArrayList<MutationTree1>();
private Alignment conc() {
   Alignment align_mut = conc(this.mutationTrees, 0);
   Alignment align_del = conc(this.mutationTrees, 1);
   Alignment align_copy = conc(this.mutationTrees, 2);
   return getMergedAlignment(align_mut, align_del, align_copy);
}*/
