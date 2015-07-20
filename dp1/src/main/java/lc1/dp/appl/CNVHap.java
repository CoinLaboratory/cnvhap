package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.io.PushbackReader;
import java.io.StringReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import lc1.CGH.Location;
import lc1.CGH.Locreader;
import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.data.collection.ASCNDataCollection;
import lc1.dp.data.collection.ArmitageCalculator;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.HLADataCollection;
import lc1.dp.data.collection.HWECalculator;
import lc1.dp.data.collection.HapDataCollection;
import lc1.dp.data.collection.IlluminaRDataCollection;
import lc1.dp.data.collection.IlluminaXYDataCollection;
import lc1.dp.data.collection.LikelihoodDataCollection;
import lc1.dp.data.collection.MatchedSequenceDataCollection;
import lc1.dp.data.collection.MergedDataCollection;
import lc1.dp.data.collection.MergedDistributionCollection;
import lc1.dp.data.collection.SequenceDataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.collection.SimpleDataCollection1;
import lc1.dp.data.collection.VNTRDataCollection;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.external.Fastphase;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.states.EmissionState;
import lc1.possel.HLALD;
import lc1.sequenced.Convert;
import lc1.util.Constants;
import lc1.util.Executor;

import org.apache.commons.cli.CommandLine;

import pal.tree.Node;
import pal.tree.ReadTree;


public class CNVHap {
    public static CNVFrame cnvf =null;
 
    public static  int rep ;
    
    static String  summaryFile;
    public static void main(String[] args){
    	System.err.println("RUNNING");
        try{
        	
            if(args.length==0){
                args = new String[] {"--paramFile", "param.txt", "--column", "1"};
            }
         
        	   CommandLine par1 = Constants.parseInitial(args);
        	   String[] cols = par1.getOptionValues("column");
        	   if(cols==null || cols.length==0) cols = new String[] {"1"};
        	  
        	   String[] index = par1.getOptionValues("index");
        	   if(index==null || index.length==0) index = new String[] {"1"};
        	   
        	   String[] repControl = par1.hasOption("indexControl") ? par1.getOptionValues("indexControl") : null;
        	   String[][] mid = getMid(new File("split.txt"), par1);
        	 //  for(int k=0; k<mid.length; k++){
        		   
        	  
            for(int ij=0; ij<cols.length; ij++){
            	System.err.println("COLUMN IS "+cols[ij]);
            	 for(int ij1=0; ij1<index.length; ij1++){
            		 System.err.println("INDEX IS "+index[ij1]);
            CommandLine para = Constants.parse(args, Integer.parseInt(cols[ij]),Integer.parseInt(index[ij1]), repControl==null ? null: repControl[ij1]);
        if(Constants.mid!=null && Constants.mid.length>1 && Constants.mid[0].length==1){
        	int len = Constants.mid.length;
        	String [][] orig = Constants.mid;
        	Constants.mid = new String[1][len];
        	for(int k=0; k<len; k++){
        		Constants.mid[0][k] = orig[k][0];
        	}
        }
            if(mid!=null && mid[0]!=null)
            {
            	for(int i=0; i<mid[0].length/2; i++){
	            String[][]sb = new String[mid.length][];
	            for(int k=0; k<mid.length; k++){
	            	sb[k] = mid[k][i*2+1].split(":");
//	            	sb.append(mid[k][1].replaceAll(":", ";"));
	//            	if(k<mid.length-1) sb.append(":");
	            }
	            System.err.println(mid[0][i*2].substring(2));
	            Field f = Constants.class.getField(mid[0][i*2].substring(2));
	        	f.set(null, f.getName().equals("mid") || f.getName().equals("restrictKb") ? sb : f.getName().equals("rsid") ? sb[0][0] : 
	        		f.getName().equals("allowLocalDist") ? Boolean.parseBoolean(sb[0][0]):
	        		f.getName().equals("phenoToAssoc") ? Constants.split(sb[0],";"):
	        		sb[0]);
            	}
            }
      int[][] mid1 = Constants.mid();
            String column  = Constants.column();
            int index1 = Constants.index;
            summaryFile = Constants.out;
          
            String dirF = Constants.baseDir;
            if(dirF.equals(".")){
            	 dirF = System.getProperty("user.dir");
            }
           File dir1 = new File(dirF);
               for(int i=0; i<mid1.length; i++){
            //	   Constants.mid = new String[][] {mid1[i]};
            //	   Constants.chrom = new String[]{mid1[i][0]};
                   if(false){//Constants.type.equals("simulate")){
                       simulate(dir1, null,(short) i);
                   }
                   else if(para ==null){
                	   String[] str = dir1.list(new FilenameFilter(){

						public boolean accept(File dir, String name) {
							return name.endsWith("_log.txt");
						}
                		   
                	   });
                	   Arrays.sort(str, new Comparator<String>(){

						
						public int compare(String o1, String o2) {
							String[] st1 =o1.split("_");
							String[] st2 = o2.split("_");
							for(int i=0; i<2; i++){
								Integer i1 = Integer.parseInt(st1[i]);
								int comp = i1.compareTo(Integer.parseInt(st2[i]));
								if(comp!=0) return comp;
							}
							return 0;
						}
                		   
                	   });
                	   inner: for(int ik=0; ik<str.length; ik++){
                		   String st = str[ik];
                		   File f = new File(dir1, st.substring(0,st.indexOf("_log.txt")));
                		   if(f.exists()){
                			   File p1 = new File(f, "raw.txt");
                			   if(p1.exists() && p1.length()>0) continue inner;
                		   }
                		   f.mkdir();
                		   File newLog = new File(f, "log.txt");
                		   
                		   copy(new File(dir1, st), newLog);
                		   String[] args1 = new String[] {"--paramFile", "log.txt", "--dir", f.getName()};
                		   Constants.parse(args1);
                		  // Constants.dir = f.getAbsolutePath();
                		  // Constants.paramFile = "log.txt";
                		   try{
                		   run(f, null);
                		 
                		 
                		   }catch(Exception exc){
                			   exc.printStackTrace();
                		   }
                	   }
                   }
                   else{
                            File outputDir = run(dir1, null);
                            if(Constants.convertAvgToZip()){
                         	  // File avg = new File(outputDir,"avg");
                         	   ConvertAssocFileToZip.main(outputDir,true);
                 		   }
                   }
                   
               }
            }
            }
        	  // }
            if(Constants.plot()==1) System.exit(0);
            Executor.shutdown();
            if(Constants.plot()==1) System.exit(0);
            //if(BaumWelchTrainer.es!=null)    BaumWelchTrainer.es.shutdown();
          //  if(ProbMultivariate.es!=null)    ProbMultivariate.es.shutdown();
        }catch(Exception exc){
            exc.printStackTrace();
        }
        Executor.shutdown();
    }

  
    
  /*  private static boolean summarise(File dir) throws Exception {
        File inp = new File(dir, dir.getName()+".dick.txt");
        if(!inp.exists()){
            Logger.global.warning("input did not exist" +inp);
            return false;
        }
        DataCollection obj = DataCollection.readDickFormat(inp);
        DataCollection obj1 = DataCollection.readFastPhaseOutput(new File(dir, "phased.txt"));
       if(obj1.size()!=obj.size()) throw new RuntimeException("!!");
       if(obj1.length()!=obj.length()){
           throw new RuntimeException("!!");
       }
       Map<Integer, List<String>> nullPos = new TreeMap<Integer, List<String>>();
       for(int j=0; j<obj.size(); j++){
           PIGData pig1 = (PIGData) obj1.get(j);
           PIGData pig2 =(PIGData)  obj.get(j);
           for(int i=0; i<obj1.get(0).length(); i++){
               Comparable[] o1 =( Comparable[])pig1.getElement(i);
               Comparable[]  o2 = ( Comparable[])pig2.getElement(i);
               if(o1[0]!=o2[0]){
                  if(o1[0]!=null && o2[0]!=null) throw new RuntimeException("!!");
               }
               if(o1[0]==null){ //ie predicted a gap
                   if(o2[0]!=null) throw new RuntimeException("!!");
                   List<String> ids = nullPos.get(i);
                   if(ids==null){
                       nullPos.put(i, ids = new ArrayList<String>());
                   }
                   ids.add(pig2.getName());
               }
           }
      }
       File output = new File(dir, "gap_poly_loc.txt");
       PrintWriter pw = new PrintWriter(new FileWriter(output));
       for(Iterator<Entry<Integer, List<String>>> it = nullPos.entrySet().iterator(); it.hasNext();){
           Entry<Integer, List<String>> nxt = it.next();
           pw.print(names.get(nxt.getKey()));
           pw.print("\t");
           pw.print(pos.get(nxt.getKey()));
           pw.print("\t");
           pw.println(nxt.getValue());
       }
       pw.close();
       return true;
        
    }*/



    public  static String[][] getMid(File file, CommandLine params) throws Exception {
		if(!file.exists() || (params!=null && params.hasOption("mid"))) {
	
		return new String[][] {null};
		}
		else{
			BufferedReader br = new BufferedReader(new FileReader(file));
			//List<String> header = Arrays.asList(br.readLine().split("\t"));
			//int loc_id = -1;
			
			//int chr_id = -1;
			//for(int i=0; i<header.size(); i++){
			//	if(header.get(i).indexOf("Loc")>=0) loc_id = i;
			//	if(header.get(i).indexOf("Chr")>=0) chr_id = i;
			//}
			String st = "";
			List<String[]> res = new ArrayList<String[]>();
			while((st = br.readLine())!=null){
				if(st.indexOf("#FINISH#")>=0) break;
				if(st.startsWith("#") || st.length()==0) continue;
				String[] str = st.split(" #")[0].split("\\s+");
				res.add(str);
				//return res.toArray(new String[0][]);
			}
			return res.toArray(new String[0][]);
		}
	}



	private static void copy(File inputFile, File outputFile) throws Exception{
    	
    	    FileReader in = new FileReader(inputFile);
    	    FileWriter out = new FileWriter(outputFile);
    	    int c;

    	    while ((c = in.read()) != -1)
    	      out.write(c);

    	    in.close();
    	    out.close();
    	
	}



	


public static void simulate(File dir, PrintWriter log,short index) throws Exception{
    File inp = new File(dir, Constants.chrom0()+".zip");
       System.err.println(inp.getAbsolutePath());
    IlluminaRDataCollection coll   = 
        new IlluminaRDataCollection(inp, index,Constants.noCopies()[0], new int[][] {new int[] {0, Integer.MAX_VALUE}}, null,null);
    coll.simulate();
//        IlluminaRDataCollection coll = new IlluminaRDataCollection(obj);
        File out = new File("data.txt");
        (new File("data.txt.loc")).delete();
        coll.writeCompressed(new File( "simulation"),false);
}
public static boolean simulateMissing = false;
public static DataCollection read(final File dir) throws Exception{
	//   Map<String, List> newvals = Constants.getSpreadsheetVals();
 //   double cn_ratio = Constants.cn_ratio();
   final  String[] format = Constants.format();
    final int[] no_copies = Constants.noCopies();
   final File compFiles = new File("extracted");
	//final  List format_new = new ArrayList(); 
	
   //	 format_new = exponentB_new =exponentR_new = snpsToPlot_new = new ArrayList();
   // 
    int[][] mid1_ = Constants.mid();
   
    if(mid1_==null){
    	String rs = Constants.rsId();
    	  outer: for(int i=0; i<format.length && mid1_==null; i++){
    		  File inp_i =  new File(Constants.inputDir(i));
    	        //DataCollection res[i];
    	        File buildF = new File(inp_i, Constants.build(i));
    	        BufferedReader br = DataCollection.getBuildReader(inp_i, buildF, null);
    	       
    	        if(br==null){
    	        	buildF = new File(inp_i, Constants.build(i).split("\\.")[0]+".gc.txt");
    	        	 br = DataCollection.getBuildReader(inp_i, buildF, null);
    	    	      
    	        }
    	        if(br!=null){
    	        	String st = "";
    	        	while((st = br.readLine())!=null){
    	        		if(st.indexOf(rs+"\t")>=0){
    	        			String[] str = st.split("\t");
    	        			int pos = Integer.parseInt(str[1]);
    	        			String chr = str[0];
    	        			if(chr.startsWith("chr")) chr = chr.substring(3);
    	        			mid1_ = new int[][] {new int[] {pos-1,pos+1}};
    	        			Constants.mid = new String[][] {new String[]{chr, (pos)+"",(pos)+""
    	        			}};
    	        			Constants.chrom = new String[] {chr};
    	        			br.close();
    	        			break outer;
    	        			}
    	        		}
    	        	 br.close();
    	        	}
    	       
    	        }
    }
   // final  String[] inpF =
    final File[] inp = new File[format.length];
    
   final double[] weight = Constants.weight();
   // Locreader[]locR = new Locreader[res.length];
   
    List<Callable> tasks = new ArrayList();
    int i3=0;
    String[] chrom = Constants.chrom;
    {
    	List l = Arrays.asList(chrom);
    	 if(l.indexOf("X")>=0 || l.indexOf("X_F")>=0 || l.indexOf("X_M")>=0){
    		 Constants.haploidMale  = true;
    	 }
    }
    /**** need to look into this for reading across boundary */
   /*if(chrom.length==1 && chrom[0].equals("X")){
	   Constants.chrom = new String[] {"X","XY"};
	   chrom = Constants.chrom;
   }*/
//   Constants.alleles = getNumbAlleles(format);
  int max = 2;
  final DataCollection[][] res = new DataCollection[chrom.length][format.length];
 
  String[] groups = Constants.equaliseGroup; //groups for which we have the same snps
  Map<String, Collection<String>> groupSNPs = new HashMap<String, Collection<String>>();
  
  inner: for(int kk=0; kk<res.length; kk++){
    for(int i1=0; i1<format.length; i1++){
    	try{
       final int i = i1;
       final int[] kb = Constants.restrictKb(i1);
       final int[][] mid_ = new int[mid1_.length][];
       for(int i2=0; i2<mid1_.length;  i2++){
       	mid_[i2]  = new int[]{ Math.max(0, mid1_[i2][0] - kb[0]),
       Math.min(Integer.MAX_VALUE, mid1_[i2][1]+kb[1])};
     }
    //   final int i2 = i3;
     // tasks.add(new Callable(){
      //  	public Object call(){
     //   		 try{
      //note need to revisit for joining x and XY
    	   try{
    		   
    		   
    		   int[][] mid = calculateRegion(new int[][] {mid_[kk]}, chrom[kk], Constants.regionsToInclude(i1), Constants.regionsToExclude(i1));
        inp[i] =  new File( Constants.inputFile(i,kk));
    /*    if(!inp[i].exists()) {
        	inp[i] = new File(inp[i].getParentFile(), "chr"+inp[i].getName());
        	 if(!inp[i].exists())
        	continue inner;
        }*/
        DataCollection resu;
        Collection<String> snpidrest = null;//groupSNPs.get(groups[i]);
        
        
        File buildF = chrom[kk].equals("pheno") ? null : new File(inp[i].getParentFile(), Constants.build(i));
        if(format[i].startsWith("geno1")) {
            resu = new SimpleDataCollection1(inp[i],(short)i, no_copies[i], mid,buildF, snpidrest);
       //     EmissionState emst = res[0].dataLvalues().next();
        }
        else if(format[i].startsWith("geno")) {
        	
            
            if(Constants.extraChrom!=null){
            	resu = new SimpleDataCollection(new File( Constants.inputDir(i)+"/" + Constants.extraChrom(0) + ".zip"),(short)i, no_copies[i], mid,buildF, snpidrest);
            	Fastphase.marks = new HashSet<Integer>();
            	if(Constants.reverse[0]) resu.reverse();
            	for(int jj=1; jj<Constants.extraChrom.length; jj++){
            		Fastphase.marks.add(resu.length);
            		DataCollection resu1 = new SimpleDataCollection(new File( Constants.inputDir(i)+"/" + Constants.extraChrom(jj) + ".zip"),(short)i, no_copies[i], mid,buildF, snpidrest);
            		if(Constants.reverse[jj]) resu1.reverse();
            		resu.addCollection(resu1,(int)Constants.initalSeparation());
            	}
            	}else{
            		resu = new SimpleDataCollection(inp[i],(short)i, no_copies[i], mid,buildF, snpidrest);
            	}
            String[] join = Constants.toJoin(i);
            if(join!=null &&  join.length>0){
            	((SimpleDataCollection)res[kk][i]).join(join);
            }
            //     EmissionState emst = res[0].dataLvalues().next();
        }
        else if(format[i].startsWith("alleleCount")) {
            resu = new lc1.dp.data.collection.SequenceAlleleDataCollection(inp[i],(short)i, no_copies[i], mid,buildF, snpidrest);
            String[] join = Constants.toJoin(i);
            if(join!=null &&  join.length>0){
            	((SimpleDataCollection)res[kk][i]).join(join);
            }
            //     EmissionState emst = res[0].dataLvalues().next();
        }
        else if(format[i].startsWith("vntr")) {
            resu = new VNTRDataCollection(inp[i],(short)i, no_copies[i], mid,buildF, snpidrest);
           // HaplotypeEmissionState ems =(HaplotypeEmissionState)  resu.dataL.get("1410");
          //  PseudoDistribution[] emissions = ems.emissions;
            String[] join = Constants.toJoin(i);
            if(join!=null &&  join.length>0){
            	((SimpleDataCollection)res[kk][i]).join(join);
            }
            //     EmissionState emst = res[0].dataLvalues().next();
        }
        else if(format[i].startsWith("hap")) {
            resu = new HapDataCollection(inp[i],(short)i, no_copies[i], mid,buildF, snpidrest);
            String[] join = Constants.toJoin(i);
            if(join!=null &&  join.length>0){
            	((SimpleDataCollection)res[kk][i]).join(join);
            }
            //     EmissionState emst = res[0].dataLvalues().next();
        }
        else if(format[i].startsWith("illuminaxy")){
            
            resu =  new IlluminaXYDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
        
       
    }
      
        else if(format[i].startsWith("illumina")){
            
                resu =  new IlluminaRDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
            
           
        }
        else if(format[i].startsWith("ascn")){
            resu =  new ASCNDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
    }
        else if(format[i].startsWith("depth")){
        	resu =  new SequenceDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
//            resu =  new DepthDataCollection(inp[i],(short)i, no_copies[i], mid, buildF);
    }
        else if(format[i].startsWith("matcheddepth")){
        	resu =  new MatchedSequenceDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
//            resu =  new DepthDataCollection(inp[i],(short)i, no_copies[i], mid, buildF);
    }
        else if(format[i].startsWith("HLA")){
            resu =  new HLADataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
        }
        else if(format[i].startsWith("likelihood")){ //used for simulation
        	resu = new LikelihoodDataCollection(inp[i],(short) i, no_copies[i], mid,  buildF, snpidrest);
        	//((SimpleDataCollection)res[i]).transform(new double[] {0.0, 1.0}, new ArrayList<String>(res[i].getKeys()), ",");
        //	((SimpleDataCollection)res[i]).dropSingle();
        //	resu = new IlluminaRDataCollection(res[i]);
        //	int size = res[i].dataL.size();
        //	int size1 = res[i].data.size();
        //	if(size!=size1) throw new RuntimeException("!!");
        	
        }
       // else if(format[i].startsWith("chiamo")) resu =  new ChiamoDataCollection(inp[i]);
        else if(format[i].startsWith("seq")){
         //   locR[i] = new Locreader(Integer.MAX_VALUE, "seq");
           Location[] ni =  Convert.getDeletions(inp[i], 
                    Integer.parseInt(Constants.chrom0()),mid[0]);
            resu = Convert.getLikelihoodDataCollection(ni);
            
           // for(int j=0; j<ni.length; j++){
           //     locR[i].add(new Location(ni[j].chr+"", ni[j].min, ni[j].max));
          //  }
            System.err.println("format seq");
        }
        else if(format[i].startsWith("cghall")){
          //  int[] mid = Constants.mid();
          //  int[] kb = Constants.restrictKb();
           // if(format[i].startsWith("seq")){
          //      kb = new int[] {0,0};
           // }
        	
            resu = new CGHDataCollection(inp[i],(short)i, no_copies[i], mid, buildF, snpidrest);
            //EmissionState ems = resu.dataL.values().iterator().next();
            //System.err.println("h "+ems);
            //   locR[i] = ((CGHDataCollection)res[i]).locR;
        }
//        else if(format[i].startsWith("depth")){
//        	resu =  new SequenceDataCollection(inp[i],(short)i, no_copies[i], mid, buildF);
//        }
        else throw new RuntimeException("!!"+format[i]);
       // if(kk==0){
        	res[kk][i] = resu;
        	if(snpidrest==null){ //update the groupSNPs for this group
        		groupSNPs.put(groups[i], new HashSet<String>(resu.snpid));
        	}else{
        		snpidrest.retainAll( new HashSet<String>(resu.snpid));
        	}
       /* }
        else{
        
        	DataCollection dc = DataCollection.append(res[i], resu);
        	if(dc instanceof MergedDataCollection){
        	res[i].transfer(dc);
        	((DistributionCollection)res[i].dc).append((DistributionCollection)resu.dc,
        			((MergedDataCollection)dc).map);
        	}
        
        }*/
    	   }
        catch(Exception exc){
        	exc.printStackTrace();
        }
   // }
       if(res[kk][i]==null){
    	   throw new RuntimeException("could not find "+inp[i]);
       }
        res[kk][i].weight = weight[i];
        res[kk][i].dropIndiv(Constants.toDel(i).split(";"));
        
        res[kk][i].dropROutliers(Constants.rOutlier(i));
        if( Constants.fillLikelihood(i)>0){
        	res[kk][i].fillLikelihoodData(
        		//false,
       			Constants.fillLikelihood(i)>=1 ? true : false, 
        					new double[] {1-Constants.fillLikelihood(i),Constants.fillLikelihood(i)});
        	/*if(Constants.fillLikelihood(i)>=1.0){
        		res[i].dc = null;
        	}*/
        	if(Constants.fillLikelihood(i)>=1 ){
        		res[kk][i].dc=null;
        	}
        }
        if(Constants.dropFracSites(i)>0 && false){
        	res[kk][i].fillLikelihoodData(
            		//false,
           			 true ,
            					new double[] {0,1}, getRandomSites(res[kk][i].loc.size(), Constants.dropFracSites(i)));
        }
/*<<<<<<< .mine
=======
       // res[kk][i].keepIndiv(Constants.toInclude(i).split(";"));
>>>>>>> .r224*/
        boolean writeExtracted  = Constants.writeExtractedFile();
        if(writeExtracted && res[kk][i].size()>0){
        	compFiles.mkdir();
        
        	res[kk][i].writeCompressed(compFiles,false);
           
        }
       /* if(Constants.fillLikelihood(i).length>0 && !Constants.fillLikelihood(i)[0].equals("0")){
        	//if(true) throw new RuntimeException("!!");
        	res[i].fillLikelihoodData(Constants.fillLikelihood(i), null);
        }*/
       
      if(Constants.restrictMarker(i)!=null)
         res[kk][i].restrictToMarkers(Constants.restrictMarker(i))	;
        //}
       // if(res[kk][i]!=null && Constants.splitByPheno(i)==null) {
        	//res[kk][i].trim(Constants.maxIndiv(i));
        	
        //}
     
        
        
      
       
        /*if(res[i].size()==0){
        	res[i]=null;
        }*/
       // if(res[kk][i].dataL.size() ==0|| res[kk][i].loc.size()==0){
       // 	res[kk][i] = null;
        //}
        //else{
        	i3++;
        //}
        Logger.global.info("free memory "+Runtime.getRuntime().maxMemory());
     //   }catch(Exception exc){
     //       exc.printStackTrace();
         //   res[i] = null;
     //   }
      //  return null;
       // 	}
      //  });
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    }
   
    BaumWelchTrainer.involeTasks(tasks, true);
   
    String str = Constants.restrictIndivTo(Constants.inputDir);
    String str1 = Constants.restrictIndivTo1();
    if(str!=null || str1!=null ){
    	 HashSet<String> indiv = new HashSet<String>();
    	    for(int i1=0; i1<format.length; i1++){
    	    	if(res[kk][i1]!=null){
    	    	indiv.addAll(res[kk][i1].dataL.keySet());
    	    	}
    	    }
    	    if(str!=null){
    	    	System.err.println("reading "+str);
    	 
    Node node = (new ReadTree(
    		new PushbackReader(new StringReader(str)))).getRoot();
   
    for(Iterator<String> it = indiv.iterator(); it.hasNext();){
    	try{
    	if(!contains(node, it.next(), res[kk])) {
    		it.remove();
    	}
    	}catch(Exception exc){
    		System.err.println("--include  did not contain "+node.getIdentifier().getName());
    		exc.printStackTrace();
    	}
    }
    	    }
    	    if(str1!=null){
    	    	List<String> str1_ = Arrays.asList(str1.split(":"));
    	    	indiv.retainAll(str1_);
    	    }
    	    //for(int k=0; k<res[kk].length; k++){
    	   // boolean cont = indiv.contains("NA19239");
   //List<String> restrictIndivTo = Arrays.asList();
    for(int i1=0; i1<format.length; i1++){
    	//for(int j=0; j<format.length; j++){
    		//if(restrictIndivTo.indexOf(Constants.inputDir[j])>=0 && res[j]!=null && res[i1]!=null){
    			if(res[kk][i1].restricToAlias(indiv)==0){
    				res[kk][i1]=null;
    			}
    		//}
    	//}
    }
    }
    
     boolean[] indexToRestrict = Constants.indexToRestrict();
     if(indexToRestrict!=null){
         Locreader[] locs = new Locreader[res.length];
         for(int i=0; i<res[kk].length; i++){
             if(indexToRestrict[i]){
            	 if(true) throw new RuntimeException("!!");
                 locs[i] = res[kk][i].getMergedDeletions(true, true);
             }
         }
         for(int i=0; i<res.length; i++){
             for(int j=0; j<locs.length; j++){
                 if(locs[j]!=null){
                	 if(true) throw new RuntimeException("!!");
                	// res[i].fix(locs[j], 0);
                 }
             }
         }
     }
  }
     DataCollection result;
  
     List<DataCollection> res_ = Arrays.asList(res[0]);//new ArrayList<DataCollection>();
     
   
     
    
     if(res_.size()==0){
    	 throw new RuntimeException("has no size");
     }
     else if(res_.size()==1){
    	//
       //  if(Constants.splitByPheno(i)==null || Constants.splitByPheno().equals("")){
        	 result = res_.get(0);
        	 result.trim(Constants.maxIndiv(0));
       //  }
        
        
     }
    else{
    	 
    //	EmissionState em = res_.get(0).dataL.get("NA19257");
    	for(int k=0; k<res_.size(); k++){
    		int maxInd = Constants.maxIndiv(k);
    		if(res_.get(k)!=null) res_.get(k).trim(maxInd);
    	}
        result =  MergedDataCollection.getMergedDataCollection(res_.toArray(new DataCollection[0]), Constants.experiment(), Constants.allowOverlaps(), true);
    }
   
    DistributionCollection.dc = result.dc;
    DataCollection.datC = result;
 //   EmissionState ems = result.dataL.values().iterator().next();
   // EmissionStateSpace.emStSp = result.getEmStSpace();
//result.trim(Constants.maxIndiv(0));
    return result;
}
    private static boolean contains(Node node, String next, DataCollection[] res) {
	if(node.isLeaf()){
		return res[Integer.parseInt(node.getIdentifier().getName())].dataL.containsKey(next);
	}
	else{
		int cnt = node.getChildCount();
		boolean or = node.getBranchLength()>0;
		boolean resu = or ?  false : true;
		for(int k=0; k<cnt; k++){
			boolean a = contains(node.getChild(k),next, res);
			resu = or ? (a || resu) : (a && resu);
		}
		return resu;
	}
}



	private static List<Integer> getRandomSites(int size, double dropFracSites) {
    	List<Integer> l = new ArrayList<Integer>();
    	for(int k=0; k<size; k++){
    		if(Constants.rand.nextDouble()<dropFracSites){
    			l.add(k);
    		}
    	}
    	return l;
}



	private static int[][] calculateRegion(int[][] mid_, String chrom,
		String[][] toinclude, String[][] toexclude) {
		String chrom1 = chrom;
		if(chrom.startsWith("chr"))chrom1 = chrom.substring(3);
    	List<int[]> mids = new ArrayList<int[]>();
    	for(int i=0; i<mid_.length; i++){
    		List<int[]> mid_new = new ArrayList<int[]>();
    		int[] mid1 = mid_[i];
    	
    		for(int k=0; k<toinclude.length; k++){
    			if(toinclude[k][0].equals("mid")){
    				toinclude[k]  = Constants.mid[0];
    			}
    			if(toinclude[k][0].equals("all")){
    				mid_new.add(mid1);
    			}
    			else if(toinclude[k][0].equals(chrom) || toinclude[k][0].equals(chrom1)){
    				int start = Constants.convert(toinclude[k][1]);
    				int end = Constants.convert(toinclude[k].length>2 ? toinclude[k][2]: toinclude[k][1]);
    				if(end<start) throw new RuntimeException("end must be greater than start "+start+" "+end);
    				int start1 = Math.max(start, mid1[0]);
    				int end1 = Math.min(end, mid1.length==1 ? mid1[0] : mid1[1]);
    				if(start1<=end1){
    					mid_new.add(new int[] {start1,end1});
    				}
    			}
    		}
    	
    		if(mid_new.size()>0){
    			
    			for(int k=0; k<toexclude.length; k++){
    				if(toexclude[k][0].equals("mid")){
        				toexclude[k]  = Constants.mid[0];
        			}
        			if(toexclude[k][0].equals("null")){
        				continue;
        			}
        			else if(toexclude[k][0].equals(chrom)|| toexclude[k][0].equals(chrom1)){
        				int start =Constants.convert(toexclude[k][1]);
        				int end = Constants.convert(toexclude[k][2]);
        				if(end<start) throw new RuntimeException("end must be greater than start "+start+" "+end);
            			
        				mid_new = remove(mid_new, start, end);
        			}
        		}
    		}
    		mids.addAll(mid_new);
    	}
	return mids.toArray(new int[0][]);
}
    
    private static List<int[]> remove(List<int[]>reg, int start, int end){
    	List<int[]> newL = new ArrayList<int[]>();
    	for(int i=0; i<reg.size(); i++){
    		int[] mid1 = reg.get(i);
    		int start0 = mid1[0];
    		int end0 = mid1.length==1 ? mid1[0] : mid1[1];
    		int start1 = Math.max(start, start0);
			int end1 = Math.min(end, end0);
			if(start1<=end1){
				if(start0<start){
					newL.add( new int[] {start0, start});
				}
				if(end0>end){
					newL.add( new int[] {end, end0});
				}
			}
			else{
				newL.add(mid1);
			}
    	}
    	return newL;
    }



	public static File run(File dir, PrintWriter log) throws Exception{
        File summ = new File(dir, summaryFile);
      //  SignificancePlot.reset();
        DataCollection obj = read(dir);
       
     //   EmissionState st = obj.dataL.get("21603");
   if(Constants.indelFilter()){    IndelFilter indelFilter = new IndelFilter(obj);
       indelFilter.filter();
   }
        if(Constants.plot()>=1){
        	//SwingUtilities.in
        	Runnable run = new Runnable(){
        		
        	public void run(){
        				cnvf = new CNVFrame();
        	}
        	};
        	Thread th = new Thread(run);
        	th.run();
        	
        		
        }
        File pedigree = new File(dir,  "ped.txt");
        if(!pedigree.exists()){
            pedigree = new File(dir.getParent(),  "ceu_ped.txt");
        }

       File outp1 = new File(Constants.outputDir);
       outp1.mkdir();
       StringBuffer lo_st = new StringBuffer();
       for(int i=0; i<Constants.mid.length; i++){
    	   if(i>0) lo_st.append("__");
    	   lo_st.append(Arrays.asList(Constants.mid[i]).toString().replaceAll("\\s+","").replaceAll("\\[", "").replaceAll("\\]", "").replace(',','_'));
     
       }
       String loc_string = lo_st.toString();
     File outputDir =  new File(outp1,Constants.getDirName(outp1, loc_string));
     System.err.println("OUTPUT DIR IS "+outputDir.getAbsolutePath());
     outputDir.mkdir();
     File paramFile = new File(Constants.outputDir(), Constants.paramFile);
   //  File dataFile = new File(Constants.outputDir(), Constants.data_params);
     File locFile = Constants.loc_params==null ? null  : new File(Constants.outputDir, Constants.loc_params);
     copy(paramFile, new File(outputDir,Constants.paramFile));
    // copy(dataFile, new File(outputDir,Constants.data_params));
     if(locFile!=null) copy(locFile, new File(outputDir,Constants.loc_params));
  //   PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputDir, "log_out.txt"))));
   //  Constants.printOptions(pw, "\n");
    // System.err.println("OUTPUTFILE IS "+outputDir);
    // pw.close();
    if(Constants.export()){
    	File export = new File(outputDir, "export");
    	export.mkdir();
    	obj.printWide(export);
    	return outputDir;
    }
     if(Constants.run()){
    	 if(obj.loc.size()==0) {
    		 throw new RuntimeException("region is too small - no probes found");
    	 }
    	 if(Constants.run2() && obj instanceof MergedDataCollection){
    		 MergedDistributionCollection mdc = (MergedDistributionCollection) obj.dc;
    		 DataCollection obj1 = ((MergedDataCollection)obj).ldl[0]; 
    		 DataCollection.datC = obj1;
    		 DistributionCollection.dc =mdc.dc[0];
    		 int plot = Constants.plot;
    		 if(plot==1) Constants.plot = 0;
    		 
    		 RunFastPhase.run(obj1, summ, new File(dir, "clusters.txt"+Constants.column()),outputDir); ;
    		EmissionState dat1 = obj1.dataL.get("22394");
    		DataCollection.datC = obj;
    		DistributionCollection.dc = mdc;
    		Constants.plot = plot;
    		if(Constants.numIt.length==2){
    			Constants.numIt = new int[] {0,5};
    		}else{
    			Constants.numIt = new int[] {5};
    		}
    		 RunFastPhase.run(obj, summ, new File(dir, "clusters.txt"+Constants.column()),outputDir);
    	 }
    	 else RunFastPhase.run(obj, summ, new File(dir, "clusters.txt"+Constants.column()),outputDir);
     }
     
   
     if(Constants.calcAssoc() && Constants.sample() && Constants.run()){
    	File  pw = new File(outputDir,"assoc");
    //	int assocTest = Constants.assocTest();
    	ArmitageCalculator.printResults(pw, obj);
     }
     if(Constants.sample() && Constants.run()){
    	 File  pw1 = new File(outputDir,"hwe");	
    	 HWECalculator.printResults(pw1, obj);
     }
     if(Constants.calcAssoc1() ){
    	 
     }
        if(log!=null){
            log.println("directory \t"+dir.getAbsolutePath());
            log.flush();
        }
       return outputDir;
    }
    
    public static String getString(String loc, List dat, List<String> ids){
        StringBuffer sb = new StringBuffer();
        int i=0;
        while(i<ids.size() && !ids.get(i).equals(loc)){
            i++;
        }
        if(i==ids.size()) throw new RuntimeException("!!");
        for(int j=0; j<dat.size(); j++){
            Object obj =((SimpleScorableObject) dat.get(j)).getElement(i);
            if(obj instanceof Comparable[]){
            Comparable[] comp = (Comparable[])obj;
            for(int k=0; k<comp.length; k++){
                sb.append(comp[k]);
            }
            }else sb.append(obj==null ? "null" : obj.toString());
            sb.append(" ");
        }
        return sb.toString();
    }
    
    public static Integer[][] calculateAverage(File summary, PrintWriter log1) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(summary));
        String st = "";
       Integer[] ff = new Integer[] {0,0,0,0};
        Integer[] sampl = new Integer[] {0,0,0,0};
        Integer[] count = new Integer[] {0, 0};
        List<String> st_sample = new ArrayList<String>();
        List<String> st_ff = new ArrayList<String>();
         while((st = br.readLine())!=null){
            if(st.startsWith("comparing with fastphase")){
                br.readLine();
                  st = br.readLine();
                  if(st!=null) {
                      st_ff.add(st);
                  count[0]++;
                  }
            }
            if(st.startsWith("comparing with sampling avg")){
                br.readLine();
                 st = br.readLine();
                 if(st!=null){
                     st_sample.add(st);
                     count[1]++;
                 }
                 
            }
        }
        for(int i=0; i<st_sample.size(); i++){
                 add(ff, st_ff.get(i));
                 add(sampl,st_sample.get(i));
        }
        log1.println(summary.getName());
        log1.println("totals ff are"+Arrays.asList(ff));
        log1.println("totals sampling are"+Arrays.asList(sampl));
        log1.flush();
        br.close();
        return new Integer[][] {ff, sampl, count};
    }
private static void add(Integer[] ff, String st) {
       if(st==null) return;
       else{
           String[] str = st.split("\\s+");
           int[] ind = new int[] {1,3,5,7};
           for(int i=0; i<ind.length; i++){
              ff[i]+=Integer.parseInt(str[ind[i]]);
           }
       }
        
    }

public static void main1(String[] args){
    try{
       
    }catch(Exception exc){
        exc.printStackTrace();
    }
}

static void searchAll() throws Exception{
    File dir1 = new File(System.getProperties().getProperty("user.dir"));
    File[] dir = dir1.listFiles(new FileFilter(){
        public boolean accept(File pathname) {
           return pathname.isDirectory();
        }
    });
    File logfile = new File(dir1, "log.txt");
   
    PrintWriter log =new PrintWriter(new BufferedWriter(new FileWriter(logfile, false)));
    for(int j=0; j<1; j++){
        Logger.getAnonymousLogger().info("doing round "+j);
        log.println("round "+j);
        log.flush();
        for(int i1=0; i1<dir.length; i1++){
            File f = dir[i1];
                if(run(f, log)==null){
                    System.err.println("already done "+f.getName());
                }
         }
    }
    log.close();
}

static void score() throws Exception{
    File dir1 = new File(System.getProperties().getProperty("user.dir"));
    File logfile1 = new File(dir1, "log_summary.txt");
    PrintWriter log1 =  new PrintWriter(new BufferedWriter(new FileWriter(logfile1, true)));
    File[] dir = dir1.listFiles(new FileFilter(){
        public boolean accept(File pathname) {
           return pathname.isDirectory();
        }
    });
    log1.println(RunFastPhase.cal.getTime().toString().toUpperCase()+"-----------------");
    Integer[] avg_ff = new Integer[] {0,0,0,0};
    Integer[] avg_sampl = new Integer[] {0,0,0,0};
    Integer[][] counts = new Integer[2][dir.length];
    for(int i1=0; i1<dir.length; i1++){
        counts[0][i1] = 0;
        counts[1][i1] = 0;
        // int i = dir.length-i1-1;
         File f = dir[i1];
             File summ = new File(f, summaryFile);
             if(!summ.exists()) continue;
             Integer[][] res = calculateAverage(summ, log1);
             //if(res[2][0]!=res[2][1]) throw new RuntimeException(res[2][0]+" "+res[2][1]);
             counts[0][i1]+=res[2][0];
             counts[1][i1]+=res[2][1];
            add(res, new Integer[][]{avg_ff, avg_sampl});
        
      }
   
         log1.println("total ff is "+Arrays.asList(avg_ff) );
         log1.println("total sampl is "+Arrays.asList(avg_sampl) );
         log1.println("counts are "+Arrays.asList(counts[0]));
         log1.println("counts are "+Arrays.asList(counts[1]));
         log1.flush();
         log1.close();
}

    private static void add(Integer[][] is, Integer[][] is2) {
    for(int i=0; i<is2.length; i++){
        for(int j=0; j<is[i].length; j++){
            is2[i][j]+=is[i][j];
        }
    }
    
}

    
    
   
    
    /*
     *  
     */
}
