package lc1.assoc;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
//import lc1.dp.data.collection.LightWeightDataCollection;
//import lc1.dp.data.collection.LightWeightIlluminaDataCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.Regression;
import lc1.util.Constants;
import lc1.util.Executor;

public class CalcAssociation1 {
	
   public static void main(String[] args){
	   try{
	   OptionBuild.getOptions(args);
	   main1();
	   }catch(Exception exc){
		   exc.printStackTrace();
	   }
   }
    public static final boolean illumina = false;
	  // static  Map<String, List<String>> sex = new HashMap<String,List< String>>();
    public static void main1(){
        try{
        //	OptionBuild.getOptions(args);
         
            final String[] tag_report =   new String[] {"snpid","loc", "index",   "maf"};
        Constants.plot = 2;
        Constants.var_thresh= new double[] {1e3};
        Constants.numThreads = OptionBuild.numThreads;
      //  Constants.topBottom = false;
        Constants.loess =new boolean[] { false};
        Constants.median_correction = new boolean[] {false};
     //   String col = "1";
      // final String sta = "XY.zip";
      //  DataCollection.scoreRegression = 1.0;
        for(int kk=0; kk<OptionBuild.experiment().length; kk++){
        	 List<String>toexcl = new ArrayList<String>();
             
           	  if(OptionBuild.exclList.length>0 && OptionBuild.exclList[kk]!=null) {
	           	  BufferedReader br = DataCollection.getBufferedReader(new File(OptionBuild.user(), "assoc/"+OptionBuild.experiment[kk]+"/"+OptionBuild.exclList[kk]));
	           	  DataCollection.readPosInfo(br, new int[] {0}, false, new List[] {toexcl}, new Class[] {String.class});
	           	  br.close();
             }
             final Set<String> toexclS = new HashSet<String>(toexcl);
        	final int kk1 =kk;
        	
		      final  Map<String, String> exclude = new HashMap<String, String>();
		        if(OptionBuild.excl.length>0 && OptionBuild.excl[kk]!=null){
		        	String[] args1 = OptionBuild.excl[kk].split(";");
		        for(int i=0; i<args1.length; i++){
		        		String[] str = args1[i].split(",");
		        		exclude.put(str[0], str[1]);
		        }
		        }
        //final String sta = args.length>0 ? args[0] : "";
       List l = new ArrayList();
      
       final String[] inputFiles = OptionBuild.input();
       final File[] fs =  (new File(OptionBuild.user(), inputFiles[0])).listFiles(new FileFilter(){

            public boolean accept(File pathname) {
             boolean res =  pathname.getName().endsWith("zip");// && VirtualDataCollection.check(pathname);
             boolean res1 = false;
             for(int i=0; i<OptionBuild.chr.length && !res1; i++){
            	 res1 = res1 ||(
            			 pathname.getName().equals(OptionBuild.chr[i]+".zip")
            			 && VirtualDataCollection.check(pathname)
            	 );
             }
            //System.err.println(pathname.getName()+" "+res+" "+res1+" "+args[0]);
             return res && res1;
            }
            
        });
        for(int i1=0; i1<fs.length; i1++){
        	final int i = i1;
        	
        		Callable call = new Callable(){
        			Set<String>nullcnt = new HashSet<String>();
        		public Object call(){
        			try{
        		String nme = fs[i].getName();
        		nme = nme.substring(0, nme.indexOf(".zip"));
               Constants.mid = 
            	   new String[][]{
            	   new String[] {nme, "0", ""+(Integer.MAX_VALUE -1)}};
           
               
                File dirOut1 = new File(OptionBuild.user(), "assoc");
                dirOut1.mkdir();
             
                final  int[] mid = new int[] {0,Integer.MAX_VALUE};
                final DataCollection[] dc  = new DataCollection[inputFiles.length];
                for(int kkk =0; kkk<dc.length; kkk++){
                	dc[kkk] = null;//
//                		illumina ? 
 //               				new LightWeightIlluminaDataCollection(new File(OptionBuild.user(), inputFiles[kkk]+"/"+fs[i].getName()),(short) 0, 2, new int[][] {mid},  null): 
  //              		new LightWeightDataCollection(new File(OptionBuild.user(), inputFiles[kkk]+"/"+fs[i].getName()),(short) 0, 2, new int[][] {mid},  null,0) ;
                	//		new LikelihoodDataCollection(new File(OptionBuild.user(), inputFiles[kkk]+"/"+fs[i].getName()),(short) 0, 2, mid,  null);
                }
                
               DataCollection dc1 = dc.length==1 ? dc[0] :  new MergedDataCollection(dc, "merged", null, true);
                //
                final  int end  = dc1.loc.size();
              //  String str = br1.1readLine();
               List<Integer> phenIndex = new ArrayList<Integer>();
                for(int ii=0; ii<dc1.pheno.phen.size(); ii++){
                	if(!exclude.containsKey(dc1.pheno.phen.get(ii)))
                		phenIndex.add( ii);
                }
               // DataCollection.scoreRegression = OptionBuild.regr;
              //  for(int l=0; l<sexP.length; l++){
               // {
                String[] expt = OptionBuild.experiment();
                	//String out_name =  "assoc-"+OptionBuild.experiment()[kk1]+".txt";
                	
                //	DataCollection dc1 = dc;//l==0 ? dc : (DataCollection) dc.clone();
                	List<String> toD = new ArrayList<String>();
                	for(Iterator<EmissionState> it = dc1.dataL.values().iterator(); it.hasNext();){
            			HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
            			if(toexclS.contains(nxt.getName())){
            				Arrays.fill(nxt.phenValue(), null);
            				System.err.println("excluding "+nxt.getName()+" in excl list");
            				//pw.flush();
            				nullcnt.add(nxt.getName());
            				//toD.add(nxt.getName());
            			}
            		}
                	for(int kk=0;kk<dc1.pheno.phen.size(); kk++)
                	{                	
                		String excl = exclude.get(dc1.pheno.phen.get(kk));
                		if(excl!=null){
                			if(dc1.pheno.type()[kk]==0){
                				double exclval = Double.parseDouble(excl);
                				for(Iterator<EmissionState> it = dc1.dataL.values().iterator(); it.hasNext();){
    	                			HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
    	                			Double phenv = nxt.phenValue()[kk];
    	                			if(phenv!=null && phenv>exclval){
    	                				Arrays.fill(nxt.phenValue(), null);
    	                				nullcnt.add(nxt.getName());
    	                			System.err.println("excluding "+nxt.getName()+" matches "+excl+" "+dc1.pheno.phen.get(kk));// in excl list");
    	                				//pw.flush();
    	                				//toD.add(nxt.getName());
    	                			}
    	                			else{
    	                				System.err.println("included "+dc1.pheno.phen.get(kk)+" "+phenv);
    	                				//pw.println("included "+dc.pheno.phen.get(kk)+" "+phenv);
    	                			}
    	                		}
                			}
                			else{
	                		int exclval = dc1.pheno.phenVals[kk].get(excl);
	                		for(Iterator<EmissionState> it = dc1.dataL.values().iterator(); it.hasNext();){
	                			HaplotypeEmissionState nxt = (HaplotypeEmissionState)it.next();
	                			Double phenv = nxt.phenValue()[kk];
	                			if(phenv!=null && Math.abs(phenv-exclval)<0.0001){
	                				Arrays.fill(nxt.phenValue(), null);
	                				nullcnt.add(nxt.getName());
	                				//pw.println("excluding "+nxt.getName()+" matches "+exclval+" "+dc.pheno.phen.get(kk));// in excl list");
	                				System.err.println("excluding "+nxt.getName()+" matches "+exclval+" "+dc1.pheno.phen.get(kk));// in excl list");
	                				//pw.flush();
	                				//toD.add(nxt.getName());
	                			}
	                		}}
                		}
                	}
                	
                	System.err.println(nullcnt.size()+" nullcnt is "+nullcnt);
                	// dc1.dropIndiv(toD.toArray(new String[0]));
                	File dirOut = new File(dirOut1, expt[kk1]);
                	dirOut.mkdir();
                	PrintWriter exclF = new PrintWriter(new FileWriter(new File(dirOut, "excl.txt")));
                	for(Iterator it = nullcnt.iterator(); it.hasNext();){
                		exclF.println(it.next());
                	}
                	exclF.close();
                	 File outF = new File(dirOut, nme);
                	 int start = -1;
                	 if(outF.exists() && OptionBuild.append){
                		 BufferedReader br = new BufferedReader(new FileReader(outF));
                		 String st = br.readLine();
                		//start++;
                		 while((st = br.readLine())!=null){
                			 start++;
                		 }
                	 }
                	boolean printh = start<0;
                	if(start<0)start =0;
                    PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter(outF, OptionBuild.append)));
                    dc1.printHapMapFormat(pw_hap1,
                    		dc1.indiv(),
                    		null, true, 
                         tag_report,
                         OptionBuild.tag_pheno,
                           //  : new String[] {"snpid", "loc"}, 
                            new String[] {}, "%7s", dc1.loc.get(start), dc1.loc.get(end-1)+1, new int[] {0,1,2}, phenIndex.toArray(new Integer[0]), printh);
                    pw_hap1.close();
                		if(dc1.log!=null){
                			dc1.log.close();
                		}
                //}
               /* int numRand =0;
                for(int ik=0; ik<numRand; ik++){
                
                    PrintWriter pw_hap1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(dirOut, "assoc"+ik+".txt"))));
                    dc.randomizePhenotypes();
                        dc.printHapMapFormat(pw_hap1, null, true, 
                        tag_report,
                           tag_pheno,
                          //  : new String[] {"snpid", "loc"}, 
                           new String[] {}, "%7s", start, end, new int[] {0,1,2,3}, new int[] {0}, true);
                    pw_hap1.close();
                }*/
        			}catch(Exception exc){
                		exc.printStackTrace();
                		Executor.shutdown();
                	}
                    return null;
        		}
        		};
        	l.add(call);
             //   if(true) System.exit(0);
               
        }
      //  pw.close();
        
        Regression.es.invokeAll(l);
        }
       Executor.shutdown();
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
   
	public static File expand(File dir, ZipFile zf,  String resEntry) throws Exception {
        // TODO Auto-generated method stub
        String[] f = resEntry.split("/");
        File par = dir;
        for(int i=1; i<f.length-1; i++){
            par = new File(par, f[i]);
            par.mkdir();
        }
        String name = f[f.length-1];
        File outF = new File(par, name);
        if(outF.exists() && outF.length()>0) return outF;
        BufferedInputStream in = new BufferedInputStream(zf.getInputStream(zf.getEntry(resEntry)));
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outF));
        byte[] buf = new byte[1024];
        int len = 0;
        int reallen = 0;

        //System.out.println(file+":"+getLocalPath()+outfile);
        while(true)
        {
            len = in.read(buf);

            //System.out.print(".");
            if(len == StreamTokenizer.TT_EOF)
            {
                break;
            }

            out.write(buf, 0, len);
            reallen += len;

        }
        in.close();
        out.close();
       return  outF;
    }
    static String findEntry(ZipFile zf, String string) {
        for(Enumeration en = zf.entries(); en.hasMoreElements();){
           ZipEntry ent =  (ZipEntry) en.nextElement();
           String nme = ent.getName();
           if(nme.indexOf(string)>=0){
               return nme;
           }
        }
        return null;
      }
    
    /*public synchronized String getPhenInfo(String string , int pos_index, int phenIndex, int type){
        int numCl = numClasses(phenIndex);
       if(pos_index!=currentPosScIndex || currentPhenScIndex!=phenIndex ){
       	Arrays.fill(regression, 0.0);
       	Arrays.fill(regressionSig, 1.0);
          // System.err.println(scoreRegression);
           if(OptionBuild.scoreRegression())
           		this.scoreRegression(pos_index, phenIndex);
           if(OptionBuild.scoreChi())   this.scoreChi1(pos_index, true, phenIndex);
           
           currentPhenScIndex = phenIndex;
           currentPosScIndex = pos_index;
           currentType = type;
       }
        if(string.startsWith("chisq")){
            
            if(numCl==2) return Format.sprintf("%5.3g", new Double[] {getSignificance(false, type)});
            else return "";
        }
        else if(string.startsWith("armitage")){
        //    this.scoreChi1(pos_index, true, string.endsWith("state"), phenIndex);
            if(numCl==2)  return Format.sprintf("%5.3g",new Double[] { getSignificance(true, type)});
            else return "";
        }
        else if(string.startsWith("regress")){
            if(regression==null) return "null";
           return Format.sprintf("%5.3g", new Double[] {this.regression[type]});
           
        }
        else if(string.startsWith("regrP")){
            if(regressionSig==null) return "null";
            return Format.sprintf("%5.3g", new Double[] {this.regressionSig[type]});
            
         }
        else{
           
            if(string.startsWith("odds_")){
            ///  System.err.println(Arrays.asList(odds[i]));
                if(numCl==2)  return Format.sprintf(formatOdds[type], this.odds[type]);
                else return "";
            }
            else if(string.startsWith("cases_")){
              
                if(numCl==2)  return Format.sprintf(formatOdds[type], this.cases[type]);
                else return "";
            }
            else if(string.startsWith("controls_")){
                if(numCl==2)  return Format.sprintf(formatOdds[type], this.controls[type]);
                else return "";
            }
            else throw new RuntimeException( "!! " +string);
        }
       
    }*/
}
