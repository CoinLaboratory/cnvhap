package lc1.dp.data.collection;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.logging.Logger;

import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.states.BackgroundEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;

//todo : we are keeping the  map to be able to map back, but probably can get rid of it if it takes too much memory
public class MergedDataCollection extends LikelihoodDataCollection {
  // public DataC[] ldl;
//private  final double[] weight;
	//EmissionState[] bgs;
	/** need to fix this for DataCollections with different underlying snp list */
	
	/*public DataCollection[] split(){
		this.get
		Map<Integer, DataCollection> m = new HashMap<Integer>
		for(Iterator<EmissionState> it = this.dataL.values().iterator(); it.hasNext();){
			HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
		}
	}*/
	@Override
	 protected void makeBackgroundState() {
		 
	 }
	@Override
	public void setDC(AbstractDistributionCollection res) {
		super.setDC(res);
		for(int i=0; i<ldl.length; i++){
			ldl[i].setDC(((MergedDistributionCollection)res).dc[i]);
		}
		
	}

	private List<AssociationCalculator> getMergedCalc(
			List<AssociationCalculator>[] lists) {
		List<AssociationCalculator> l = new ArrayList<AssociationCalculator>();
		LinearRegressionCalculator[] li = new LinearRegressionCalculator[lists.length-1];
		for(int i=0; i<lists[0].size(); i++){
			for(int j=0; j<li.length; j++){
				li[j] = (LinearRegressionCalculator) lists[j].get(i);
			}
			l.add(new LinearRegressionCalculator(li, this));
		}
		return l;
	}
	@Override
	public int getIndex(String name){
	       
		for(int i=0; i<ldl.length; i++){
			if(ldl[i].dataL.containsKey(name)) return i;
		}
		return -1;
    }
	@Override
	public void initialisationStep() {
		// TODO Auto-generated method stub
		for(int i=0; i<ldl.length; i++){
			ldl[i].initialisationStep();
		}
		
	}
	@Override
	public PseudoDistribution getBGDist(int data_index, int pos){
		return this.ldl[data_index].getBGDist(0, this.map[data_index][pos].relative_position);
	}
	@Override
	public Integer getBGCount(int data_index, int pos){
		if(data_index<0) return Constants.backgroundCount(data_index);
		return ldl[data_index].getBGCount(0, pos);
		//return this.ldl[data_index].getBGCount(0, this.map[data_index][pos].relative_position);
	}
	//public double[] getBGCountDist(int data_index, int pos){
		//return this.ldl[data_index].getBGCountDist(0, this.map[data_index][pos].relative_position);
	//}
	
	
	public void writeCompressed(File dir, boolean writeDistr){
		
		if(Constants.demergeToPrint()){
		this.updateOriginal();
	    for(int i =0; i<this.ldl.length; i++){
	    	ldl[i].writeCompressed(dir,writeDistr);
	    }
		}
		else{
			  dir.mkdir();
			  	String[] origI = Constants.originalNames();
			  	
			for(int i =0; i<origI.length; i++){
			      // String inp = Constants.inputDir[i].replace('/', '_');
				   String[] inp1 = origI[i].split("/");;;
				   String inp = inp1[inp1.length-1];
				   int ind = inp.indexOf('*');
				   if(ind>=0) inp = inp.substring(0,ind);
				   int cnt =0;
				   List<List<String>> indiv = new ArrayList<List<String>>();
				   for(int k=0; k<ldl.length; k++){
					   if(ldl[k].loc.size()>0){
						   if(ldl[k].name.indexOf(inp)>=0){
							   indiv.add(ldl[k].indiv());
							   cnt+=ldl[k].indiv.size();
						   }
					   }
				   }
				  if(cnt==0) continue;
				    File out = new File(dir, inp);
				    try{
				    	int[] core = Constants.core();
				        CompressDC cc = new CompressDC( out, this);
				        cc.run(indiv, core[0], core[1]);
				    }catch(Exception exc){
				        exc.printStackTrace();
				    }
				
		    }
			
			
		}
	   
	}
	@Override
	public boolean hasIntensity(int i) {
		Info[] ma = this.map[i];
		for(int k=0; k<ma.length; k++){
			if(ma[k]!=null && this.ldl[k].hasIntensity(ma[k].relative_position)) return true;
		}
		return false;
	}
	
	
	boolean[] contains;
	public   void setData(String key, PIGData dat, HaplotypeEmissionState datL, PIGData datvit, HaplotypeEmissionState datvitL, double[] certainty){
	  if(contains==null)contains = new boolean[ldl.length];
	   for(int i=0; i<ldl.length; i++){
		   contains[i] = ldl[i].containsKey(key);
	   }
		if(Constants.writeMergedAverageFiles()) super.setData(key, dat, datL, datvit, datvitL, certainty);
		else{
			for(int i=0; i<ldl.length; i++){
		
		//	System.err.println("key");
		
			if(contains[i]){
				/*if(groupToPrint!=null){
					String nme = ldl[i].name.split("_")[0];
					i = groupToPrint.get(nme).get(0);
				}*/
				ldl[i].setData(key, dat, datL, datvit, datvitL, certainty);
			
			}
		}
		}
	   
	}
	
	//Map<String, List<Integer>> groupToPrint = null;
	
	public void printcompressed(List<String> keys) throws Exception {
		if(Constants.writeMergedAverageFiles()) super.printcompressed(keys);
		else{
			for(int i=0; i<ldl.length; i++){
		
			ldl[i].printcompressed(keys);
			}
		}
	}
	@Override
	public void initialisePrinting(File parent) throws Exception{
		//if(Constants.prin)
		if(Constants.writeMergedAverageFiles()) super.initialisePrinting(parent);
		else{
		//	this.groupToPrint = this.getGroups();
		//	if(groups)
			for(int i=0; i<this.ldl.length; i++){
				ldl[i].initialisePrinting(parent);
				ldl[i].loc = this.loc;
				ldl[i].snpid = this.snpid;
			}
		}
	}
	 public void printWide(File file) throws Exception{
		 for(int i=0; i<this.ldl.length; i++){
			 ldl[i].printWide(file);
		 }
	 }
	public void writeResults(String key) throws Exception{
		if(Constants.writeMergedAverageFiles()) {
			super.writeResults(key, this.alleleA, this.alleleB);
		}
		this.updateOriginal(key);
		for(int i=0; i<ldl.length; i++){
			ldl[i].writeResults(key, this.alleleA, this.alleleB);
		}
	}
	public void finishedPrinting(){
		super.finishedPrinting();;
		for(int i=0; i<ldl.length; i++){
			ldl[i].finishedPrinting();
		}
	}
	@Override
	public void writeAverages(File file, String[] strings) throws Exception {
		if(Constants.writeMergedAverageFiles()){
			super.writeAverages(file, strings);
		}
		else{
			//Map<String,List<String>> groups = this.getGroups();
			file.mkdir();
			
			System.err.println("made dir "+file.getAbsolutePath());
			Object[] obj = getIndiv();
			Map<String, PIGData>[] data1 = (Map<String, PIGData>[])obj[0];
			Map<String, double[]>[] uncPh1 =(Map<String, double[]>[]) obj[1];
			String[] inp_name = (String[]) obj[2];
			int[] core = Constants.core();
	//		file.mkdir();
			for(int i=0; i<ldl.length; i++){
				//String key = it.next();
				List<String> indiv = new ArrayList<String>(ldl[i].dataL.keySet());
				
				if(indiv.size()>0){
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(file,ldl[i].name+".txt"))));
				printHapMapFormat( 
		                 pw, indiv, null, true, 
		                 new String[] {"snpid", "loc"}, 
		                 new String[0],
		                  strings, "%7s", core[0], core[1]);
			    pw.close();
				}
			}
			PrintWriter pw_hap2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(file,this.name+".txt"))));
			//PrintWriter pw_hap2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(outdir, chrom+".txt"))));
			// dc.writeFastphase(pw_hap2, false);
		
	          //pw_hap1.close();
	         pw_hap2.close();
		}
		
	}
	
	
	private Map<String,List<Integer>> getGroups() {
		Map<String, List<Integer>> m = new HashMap<String, List<Integer>>();
		for(int k=0; k<this.ldl.length; k++){
			String nme1 = ldl[k].name.split("_")[0];
			List<Integer> l = m.get(nme1);
			if(l==null){
				m.put(nme1, l = new ArrayList<Integer>());
			}
			l.add(k);
		}
		return m;
	}
	public Object[] getIndiv(){
		String[] origI = Constants.originalNames();
		List<Integer>[] inp = new List[origI.length];
		String[] inp_name = new String[origI.length];
		for(int i=0; i<inp.length; i++){
			inp[i] = new ArrayList<Integer>();
			if(origI[i]==null) continue;
			String st = i+"";
			try{
			 String[] inp1 = origI[i].split("/");;;
			  st = inp1[inp1.length-1];
			}catch(Exception exc){
				exc.printStackTrace();
				System.err.println("problem with "+Arrays.asList(origI));
			}
			   int ind = st.indexOf('*');
			   if(ind>=0) st = st.substring(0,ind);
			   inp_name[i] = st;
			   for(int j=0; j<ldl.length; j++){
				   if(ldl[j].name.indexOf(st)>=0){
					   inp[i].add(j);
				   }
			   }
		}
		Map<String, PIGData>[] data1 = new Map[origI.length];
		Map<String, double[]>[] uncPh1 = new Map[origI.length];
		for(int i=0; i<data1.length; i++){
			data1[i] = new HashMap<String, PIGData>();
			uncPh1[i] = new HashMap<String, double[]>();
		}
		for(Iterator<String> it = data.keySet().iterator(); it.hasNext();){
			
			String key = it.next();
			PIGData val = data.get(key);
			for(int ind=0; ind<origI.length; ind++){
				 for(Iterator<Integer> it1 = inp[ind].iterator(); it1.hasNext();){
					 int j = it1.next();
				     if(ldl[j].indiv().contains(key)){
				    	 data1[ind].put(key, val);
				    	 uncPh1[ind].put(key, this.uncertaintyPhase.get(key));
					
				}
				 }
				
			}
//			else System.err.println("WARNING negative index for "+val.getName());
				
		}
		return new Object[] {data1, uncPh1, inp_name};
	}
	
	
	@Override 
	public PrintWriter[] getPrintWriter(File f) throws Exception{
		f.mkdir();
		PrintWriter[] pw = new PrintWriter[this.ldl.length];
		for(int k=0; k<pw.length; k++){
			pw[k] = new PrintWriter(new BufferedWriter(new FileWriter(new File(f, ldl[k].name+".txt"))));
		}
		return pw;
	}
	
	public  void writeFastphase(File pw_hap2, boolean states) throws Exception{
		//pw_hap2.mkdir();
		//super.writeFastphase(pw_hap2, states)
		if(this.loc.size()==0) return;
		Object[] obj = getIndiv();
		Map<String, PIGData>[] data1 = (Map<String, PIGData>[])obj[0];
		Map<String, double[]>[] uncPh1 =(Map<String, double[]>[]) obj[1];
		String[] inp_name = (String[]) obj[2];
		for(int i=0; i<data1.length; i++){
			if(data1[i].size()>0){
			PrintWriter[] pw = this.getPrintWriter(pw_hap2);
				//new PrintWriter(new BufferedWriter(new FileWriter(new File(pw_hap2,inp_name[i]+".txt"))));
		    if(states){
		        
		    }
//		this.writeFastphase(this.,uncertainty,uncertaintyPhase, pw,  false, true, false,calculateMaf1().getConstantPos());
		    else{
		        this.writeFastphase(data1[i],uncPh1[i], pw,  false, true, false,null);
		    }
		   for(int k=0; k<pw.length; k++){
			   pw[k].close();
		   }
			}
		}
		
	}
	
	
	
@Override
public String[] getUnderlyingDataSets() {
    return dataTypes;
 }
    public final String[] dataTypes;
  //  Map<String, Double> weights = new HashMap<String, Double>();
   
    @Override
    protected int getNumberDataTypes() {
        return dataTypes.length;
     }
   
   
    public static DataCollection getMergedDataCollection(DataCollection[] ldl1, String name, String[][] overl, boolean b){
    	  if(Constants.offsetSamples()){
      		Set<Integer> positions = new HashSet<Integer>();
      		for(int k=ldl1.length-1; k>=0 ;k--){
      			ldl1[k].offset(positions);
      		//	positions.addAll(l.get(k).loc);
      		}
      	//  System.err.println(ldl1[0].loc.get(217));  System.err.println(ldl1[1].loc.get(500));
      	}
        List<DataCollection> l = new ArrayList<DataCollection>();
        int dropped=0;
        for(int i=0; i<ldl1.length; i++){
     	   if(ldl1[i]!=null ){ //&& ldl1[i].dataL.size()>0 && ldl1[i].loc().size()>0){
     		   l.add(ldl1[i]);
     		 if(ldl1[i].name.equals(name)){
     			   name = name+"_expt";
     		   }
     		 if(dropped>0){
     			 for(Iterator<EmissionState> it = ldl1[i].dataLvalues(); it.hasNext();){
     				 EmissionState st = it.next();
     				 st.setIndex(st.getIndex()-dropped);
     			 }
     		 }
     	   }
     	   else{
     		   System.err.println("no samples for "+(ldl1[i]!=null ? ldl1[i].name : i+""));
     		   dropped++;
     	   }
        }
        if(l.size()==1) return l.get(0);
         /*for(int i=1; i<ldl1.length; i++){
             Set<String> alias = new HashSet<String>(ldl1[0].getKeys());
            Set<String>keys1 = new HashSet<String>(ldl1[i].getKeys());
            alias.removeAll(keys1);
            keys1.removeAll(ldl1[0].getKeys());
            if(keys1.size()>0 || alias.size()>0) throw new RuntimeException("\n"+keys1+"\n"+alias);
//             ldl1[i].restricToAlias(alias);
         }*/
      /*  for(int i=0; i<l.size(); i++){
     	   for(Iterator<EmissionState> it = l.get(i).dataLvalues(); it.hasNext();){
     		  ((HaplotypeEmissionState) it.next()).setDataIndex(i);
     	   }
        }*/
    //	EmissionState emst = l.get(0).dataL.values().iterator().next();
      
      //  System.err.println(l.get(0).loc.indexOf(152458882));
       // System.err.println(l.get(1).loc.indexOf(152458882));
         return new MergedDataCollection(l.toArray(new DataCollection[0]), name,overl, b);
     }
  
    
    
   // Map<String, int[]> phenAlias;
   
  //  public List<int[]> transfer;
  //  List< Map<Integer, Info>> list;
    // in case of conflicts takes earliest one first (although it will report conflicts)
   
    public double[] weights() {
		double[] res = new double[ldl.length];
		for(int i=0; i<res.length; i++){
			res[i] = ldl[i].weight;
		}
		return res;
	}

   /* public void dropIndiv(String[] toDel) {
    	throw new RuntimeException("!!");
    }*/
    
    public File[] getPhenoFiles(){
    	File[] res = new File[ldl.length];
    	for(int i=0; i <res.length; i++){
    		File[] f = ldl[i].getPhenoFiles();
    		res[i] = f==null ? null:f[0];
    			//new File(Constants.baseDir,Constants.inputDir[ldl[i].index]+"/pheno.txt");
    	}
    	return res;
    }
    public Map<String, String>[] getRenaming() {
    	Map[] res = new Map[ldl.length];
    	for(int i=0; i<res.length; i++){
    		res[i] = ldl[i].rename;
    	}
    	return res;
    }
    
    public MergedDataCollection(DataCollection[] ldl, String name, String[][] allowedOverlaps, boolean makeDc)
    {
    
    this(ldl,name, allowedOverlaps, makeDc, null);
    	//EmissionState a = ldl[0].dataL.get("T5078");
    	//EmissionState b = ldl[1].dataL.get("T5078");
    	//System.err.println("h");
    }
    
    public double baf(int i) {
    	Info[] ma = this.map[i];
    	double d = 0;
    	double cnt=0;
		for(int k=0; k<ma.length; k++){
			if(ma[k]!=null){
				d+=this.ldl[k].baf(ma[k].relative_position) ;
				cnt++;
			}
		}
    	return d/cnt;
    }
    public MergedDataCollection(DataCollection[] ldl, String name, String[][] allowedOverlaps,  boolean makeDc, Info[][]map1){
    //	EmissionState st = ldl[2].dataL.get(ldl[0].indiv.get(0));
    //	System.err.println("hj");
      //  EmissionState st = ldl[0].dataL.get("NA06985");
       // HaplotypeEmissionState st1 = (HaplotypeEmissionState) ldl[0].dataL.get("NA06984");
    	if(Constants.intersectSamples()){
    		Set<String> intersection = new HashSet<String>(ldl[0].indiv());
    		for(int k=1; k<ldl.length; k++){
    			intersection.retainAll(ldl[k].indiv());
    		}
    		for(int k=0; k<ldl.length; k++){
    			Set<String> torem = new HashSet<String>(ldl[k].indiv());
    			torem.removeAll(intersection);
    			ldl[k].dropIndiv(torem.toArray(new String[0]));
    		}
    
    	}
    	
    	
    	if(allowedOverlaps!=null){
    	for(int i=0; i<ldl.length; i++){
    		if(allowedOverlaps[i]!=null && 
    				 allowedOverlaps[i].length==1
    				&& allowedOverlaps[i][0].equals("all"))
    			
    				continue;
    		List<String> l = allowedOverlaps[i]==null ?
    				new ArrayList<String>() : Arrays.asList(allowedOverlaps[i]);
    		List<String> indiv =  ldl[i].indiv();
    		inner: for(int j=i+1; j<ldl.length; j++){
    			if(allowedOverlaps[j]!=null && 
       				 allowedOverlaps[j].length==1
       				&& allowedOverlaps[j][0].equals("all"))
       			
       				continue inner;
       		List<String> lj = allowedOverlaps[i]==null ?
       				new ArrayList<String>() : Arrays.asList(allowedOverlaps[j]);
    			
    			
    			if(l.indexOf(ldl[j].name)<0 && lj.indexOf(ldl[i].name)<0){
    				Set<String> indiv_j = new HashSet<String>( ldl[j].indiv());
    				indiv_j.retainAll(indiv);
    				if(indiv_j.size()>0){
    					ldl[j].rename(indiv_j);
    					ldl[i].rename(indiv_j);
    					indiv = ldl[i].indiv();
    				}
    			}
    		}
    	}
    	}
    	
        this.pheno = new MergedPhenotypes(ldl);
        this.name = name;
        this.chrom = ldl[0].chrom;
        this.ldl = ldl;
      
        this.dataTypes = new String[ldl.length];
    
        AbstractDistributionCollection[] dc = new AbstractDistributionCollection[ldl.length];
     
        String[] nme = new String[ldl.length];
        for(int i=0; i<dc.length; i++){
        //	 dc[i] = ldl[i].dc;
        	 nme[i] =ldl[i].name.split("_")[0];
        }
       SortedMap<String,Integer> keys = new TreeMap<String,Integer>();
       List<Integer> posToDelete = new ArrayList<Integer>();
       List[] posToDel = new List[ldl.length];
       for(int k=0; k<posToDel.length; k++){
    	   posToDel[k] = new ArrayList<Integer>();
       }
      if(map==null){
       Map<String, Integer> idToPos= getPosLoc(ldl, Constants.includeSNPS());
        compareAlleles(idToPos, posToDelete, posToDel);
     //   for(int i=0; i<posToDelete; i++)
        this.drop(posToDelete, false);
        List<Info[]> map_n = new ArrayList<Info[]>();
/*        if(posToDelete.size()>0){ //this may be a bit of a hack!!!
        	for(int k=0; k<posToDel.length; k++){
        		ldl[k].drop(posToDel[k], false);
        		posToDel[k].clear();
        	}
        	posToDelete.clear();
        	 idToPos= getPosLoc(ldl, Constants.includeSNPS());
        	 compareAlleles(idToPos, posToDelete, posToDel);
        	 EmissionState hes = ldl[1].dataL.get("NA18561");
        	 System.err.println("h");
	    //	 EmissionState hes1 = ((MergedDataCollection)DataCollection.datC).ldl[1].dataL.get("NA18561");
	    	
        	// System.err.println("h");
        }*/
        for(int i=0; i<map.length; i++){
        	//if(!posToDelete.contains(i))
        		map_n.add(map[i]);
        }
        map = map_n.toArray(new Info[0][]);
       /* for(int k=0; k<ldl.length; k++){
        	ldl[k].drop(posToDel[k], false);
        }*/
      }
      else{
        this.map = map1;
      }
      for(int i=0; i<dc.length; i++){
     	 dc[i] = ldl[i].dc;
     	
     }
        if(!allnull(dc) && makeDc)
        this.dc = new MergedDistributionCollection(dc, this.map,Constants.emissionGroup());
        
        for(int i=0; i<ldl.length; i++){
            if(ldl[i].loc.size()>0){
            	for(Iterator<String> it = ldl[i].getKeys().iterator(); it.hasNext();){
            		String key = it.next();
            		keys.put(key,ldl[i].dataL.get(key).noCop());
            	}
              
            }
           
            dataTypes[i] = ldl[i].name;
        }
        List<String> keysL = new ArrayList<String>(keys.keySet());
        List<Integer> ploidy = new ArrayList<Integer>(keys.values());
       Logger.global.info("creating merged data set ");
        this.length = loc.size();
        
     
        this.probeOnly = new Boolean[map.length];
        inner: for(int i=0; i<map.length; i++){
           Info[] m_i = map[i];
            Info first = getFirst(m_i);
            if(first==null){
            	throw new NullPointerException("problem with "+i);
            }
            int pos = first.loc;
            this.alleleA.add(first.alleleB);
            this.alleleB.add(first.alleleA);
            probeOnly[i] = this.probeOnly(i);//ldl[first.data_index].probeOnly[first.relative_position];
          
        }
        this.length = loc.size();
        this.createDataStructure(keysL, ploidy, null);
        CompoundEmissionStateSpace[] emstsp = new CompoundEmissionStateSpace[Constants.maxPloidy()];
        for(int k=0; k<emstsp.length; k++){
        	emstsp[k] = Emiss.getSpaceForNoCopies(k+1);
    		int[] noBIndices = emstsp[k].noBIndex();
    		
    		EmissionStateSpace memb = emstsp[k].getMembers()[0];
    		
    		
    	
    			//new SimpleExtendedDistribution(emstsp.defaultList.size(), Double.POSITIVE_INFINITY);
    		//probeOnly.emstsp = emstsp;
    		//unknownDist.emstsp = emstsp;
    	}
        Set<String> problems = new HashSet<String>();
   
        prepareStates(keysL, problems, emstsp);
     
      // this.calculateMaf(true);
        Logger.global.info("h");
       if(Constants.modelbg){
    	   this.bg = (BackgroundEmissionState) this.dataL.get(ldl[0].bg.name);
    	   this.bg1 = (PhasedDataState) this.data.get(ldl[0].bg.name);
       }
     
    }
   
   
  

  
protected void prepareStates(List< String > keysL, Set<String> problems, CompoundEmissionStateSpace[] emstsp) {
	
	for(int i=0; i<map.length; i++){
     	  Info[] m_i = map[i];
           Info first = getFirst(m_i);
           if(first==null) continue;
           Map<Integer, Integer> dataIndex = new HashMap<Integer, Integer>();
           List<String> probs = new ArrayList<String>();
         outer: for(int k=0; k<keysL.size(); k++){
             String key = keysL.get(k);
         //    if(perfectMatch.contains(key)) continue;
             HaplotypeEmissionState st = (HaplotypeEmissionState) this.dataL(key);
             
             st.setDataIndex(this.getIndex(st.name));
            // HaplotypeEmissionState st1 = (HaplotypeEmissionState) this.data.get(key);
             if(Constants.modelbg() && st.name.equals(ldl[0].bg.name) && allequalBG(ldl)){
             	dataL.put(st.name, ldl[0].bg);
             	data.put(st.name, ldl[0].bg1);
             	continue outer;
             }
           //  boolean contained = false;
             for(int dat_ind =0; dat_ind < ldl.length; dat_ind++){
             	//if() continue;
                // int dat_ind = it.next();
                 if(m_i[dat_ind]!=null && ldl[dat_ind].containsKey(key)){
                    
                     Integer count = dataIndex.get(dat_ind);
                     dataIndex.put(dat_ind, count==null ? 1 : count+1);
                     Info inf = m_i[dat_ind];
                     int position_within = inf.relative_position;
                     if(dat_ind!=inf.data_index) throw new RuntimeException("!!");
                     HaplotypeEmissionState within  = (HaplotypeEmissionState)ldl[dat_ind].dataL(key);
                   PseudoDistribution psd  = (within).emissions[position_within];
                   if(psd!=null)psd.setDataIndex((short)dat_ind);
                   try{
                   boolean success =   st.updateEmissions(i, psd);
                    if(!success){
                        System.err.println("problem with "+key+" at "+snpid.get(i)+" "+i+" "+loc.get(i));
                        problems.add(key);
                    }
                   }catch(Exception exc){
                	
                	   exc.printStackTrace();
                	      System.err.println("problem with "+key+" at "+snpid.get(i)+" "+i+" "+loc.get(i));
                	   probs.add(key);
                	  
                       st.emissions[i] = null;
                	
                   }
                 
                 }
            
                 
             }
            
        
             if(st.emissions[i]==null){
            //	 System.err.println(first);
            	  if(first.isProbeOnly() !=null && first.isProbeOnly()){
             		
             		st.emissions[i] = emstsp[st.noCop()-1].getHWEDist1(0.0);
             	}
             	else{
             		st.emissions[i] = emstsp[st.noCop()-1].getHWEDist1(null);
             			
             	}
                 st.emissions[i].setDataIndex((short)-2);
              //  if(st1!=null){ 
             //	   st1.emissions[i] = new IntegerDistribution(0, emStSp);
              //   st1.emissions[i].setDataIndex((short)-2);
             //   }
                // if(st.getName().startsWith("NA")){
                //     System.err.println("hap map not contained "+first);
               //  }
             }
            /*PseudoDistribution dist =  ((HaplotypeEmissionState)maf).emissions[i];
            for(Iterator<Integer> it = dataIndex.keySet().iterator(); it.hasNext();){
         	   Integer ind = it.next();
         	   Integer cnt = dataIndex.get(ind);
         	   int relpos = m_i[ind].relative_position;
         	 //  double[]  prob =  ((HaplotypeEmissionState)ldl[ind].maf).emissions[relpos].probs();
         	   //for(int kk=0; kk<prob.length; kk++){
         		//   dist.addCount(kk, prob[kk]*cnt.doubleValue());
         	   //}
            }
           // dist.transfer(0.0);
            //dist.initialise();*/
             //so there is no data for this individual at this marker
         }
           if(probs.size()>0){
          	 try{
          		 throw new RuntimeException("PROB with imputation at "+snpid.get(i)+" "+probs);
          	 }catch(Exception exc){
          		 exc.printStackTrace();
          	 }
           }
     }
	   
     {
         outer: for(int k=0; k<keysL.size(); k++){
             String key = keysL.get(k);
             HaplotypeEmissionState st = (HaplotypeEmissionState) this.dataL(key);
             if(st.dataIndex()<0) {
                 System.err.println("nothing contained "+key+" ");
             }
         }
     }
		
	}
private boolean allnull(AbstractDistributionCollection[] dc) {
	  for(int i=0; i<dc.length; i++){
		  if(dc[i]!=null) return false;
	  }
		return true;
	}
private boolean allequalBG(DataCollection[] ldl2) {
	  for(int i=1; i<ldl2.length; i++){
		  if(ldl2[i].loc.equals(ldl2[0].loc)) return false;
	  }
		return true;
	}
private lc1.dp.data.collection.Info getFirst(
		lc1.dp.data.collection.Info[] m_i) {
	for(int i=0; i<m_i.length; i++){
		if(m_i[i]!=null) return m_i[i];
	}
	return null;
}


/* returns a map from string id to relative overall position
   * also initialises loc and position
   *  */  
  private Map<String, Integer> getPosLoc(DataCollection[] ldl2, boolean[] includeSNPS) {
	  Map<String, Integer> m1= new HashMap<String, Integer>();
	SortedMap<Integer, Set<String>> m= new TreeMap<Integer, Set<String>>();
	for(int i=0; i<ldl2.length; i++){
		if(includeSNPS==null || includeSNPS[i]){
		for(int j=0; j<ldl2[i].loc.size(); j++){
			Integer pos_ij = ldl2[i].loc.get(j);
			String string_ij = ldl2[i].snpid.get(j);
		//	if(string_ij.equals("A_16_P11960872")){
		//		Logger.global.info("h");
		//	}
			Set<String> id1 = m.get(pos_ij);
			
			
			
			
			
			Integer pos1 = m1.get(string_ij);
			
			if(pos1!=null){
				int diff = Math.abs(pos1-pos_ij);
				if(diff>Constants.maxCoordDiff())throw new RuntimeException(ldl[i].snpid.get(j)+" !! coordinates do not match "+pos1+" vs "+pos_ij);
				else if(diff>0){
					 pos_ij = pos1;
					 ldl2[i].loc.set(j, pos1);
					 id1 = m.get(pos_ij);
					 
				}
			}
			if(id1==null){
				m.put(pos_ij, id1 = new HashSet<String>(1));
			}
			
			id1.add(string_ij);
			m1.put(string_ij, pos_ij);
			
			//if(id1.size()>1){
			////	System.err.println("gh");
			//}
		
		}
		}
	}
	this.loc = new ArrayList<Integer>();
	this.snpid = new ArrayList<String>();
	for(Iterator<Entry<Integer, Set<String>>> it = m.entrySet().iterator(); it.hasNext();){
		Entry<Integer, Set<String>> nxt = it.next();
		for(Iterator<String> it1 = nxt.getValue().iterator(); it1.hasNext();){
			//System.err.println(loc.size()+" "+nxt);
			//if(loc.size()>=1161){
			//	System.err.println("here");
		//	}
			if(loc.contains(nxt.getKey())){
				System.err.println("WARNING different  snps at the same location "+nxt.getKey());
				loc.add(nxt.getKey()+1);
				snpid.add(it1.next());
			}else{
				loc.add(nxt.getKey());
				snpid.add(it1.next());
			}
		}
	}
	Map<String, Integer> m2 = new HashMap<String, Integer>();
	for(int i=0; i<snpid.size(); i++){
		//if(i==1162){
		//	Logger.global.info("h");
		//}
		m2.put(snpid.get(i), i);
	}
	//int i1 = m2.get("A_16_P11960872");
	//int i2 = m2.get("A_16_P01721604");
	return m2;
}


public void updateOriginal(){
	 List<String> keysL = new ArrayList<String>(this.dataL.keySet());
	 for(int i=0; i<map.length; i++){
      Info[] m_i = map[i];
         Info first = getFirst(m_i);
       outer: for(int k=0; k<keysL.size(); k++){
           String key = keysL.get(k);
           HaplotypeEmissionState st = (HaplotypeEmissionState) this.dataL(key);
           HaplotypeEmissionState st1 = (HaplotypeEmissionState) this.data.get(key);
           for(int dat_ind=0; dat_ind<ldl.length; dat_ind++){
              // int dat_ind = it.next();
               if(ldl[dat_ind].containsKey(key)){
                   //st.setDataIndex(dat_ind);
                   Info inf = m_i[dat_ind];
                   int position_within = inf.relative_position;
                   if(dat_ind!=inf.data_index) throw new RuntimeException("!!");
                   HaplotypeEmissionState within  = (HaplotypeEmissionState)ldl[dat_ind].dataL(key);
                 
                  within.emissions[position_within] = st.emissions[i];
                 
                   HaplotypeEmissionState emst =(HaplotypeEmissionState)ldl[dat_ind].data.get(key); 
                   emst.emissions[position_within] =  (st1).emissions[i];
//                       emst.updateEmissions(position_within, );
                       
                  
                 
                //   continue outer;
               }
           }
         
          
           //so there is no data for this individual at this marker
       }
}
}

public void updateOriginal(String key){
	// List<String> keysL = new ArrayList<String>(this.dataL.keySet());
	HaplotypeEmissionState st = (HaplotypeEmissionState) this.dataL(key);
    HaplotypeEmissionState st1 = (HaplotypeEmissionState) this.data.get(key);
   PIGData vit = this.viterbi.get(key);
   HaplotypeEmissionState vit1 = this.viterbiL.get(key);
    for(int dat_ind=0; dat_ind<ldl.length; dat_ind++){
    	HaplotypeEmissionState within  = (HaplotypeEmissionState)ldl[dat_ind].dataL(key);
    	if(within!=null){
    		ldl[dat_ind].data.put(key,(PIGData)st1);
    		ldl[dat_ind].dataL.put(key, st);
    		ldl[dat_ind].viterbi.put(key, vit);
    		ldl[dat_ind].viterbiL.put(key, vit1);
	    /*PhasedDataState emst =(PhasedDataState)ldl[dat_ind].data.get(key);
	        if(emst==null){
	        	emst =  SimpleScorableObject.make(key,loc.size(), st1.getEmissionStateSpace(),this.index);
	        	ldl[dat_ind].data.put(key, emst);
	        }
		 for(int i=0; i<map.length; i++){
	     Info[] m_i = map[i];
	        Info first = getFirst(m_i);
	 
	              if(ldl[dat_ind].containsKey(key)){
	              
	                  Info inf = m_i[dat_ind];
	                  int position_within = inf.relative_position;
	                  if(dat_ind!=inf.data_index) throw new RuntimeException("!!");
	               
	                 within.emissions[position_within] = st.emissions[i];
	                
	                  
	                  emst.emissions[position_within] =  (st1).emissions[i];
	              }
	          }*/
        
    	}
}
}
    class Comparison{
        int index1;
        int index2;
        String rsid;
        Set<String> keys;
        
        
        
    }
    
    private Set<String>[][] getCommonKeys(DataCollection[] ldl2){
        Set<String>[][] res = new HashSet[ldl2.length][ldl2.length];
        for(int i=0; i<ldl2.length; i++){
            Set<String> keys1 = ldl2[i].getKeys();
            for(int j=i+1; j<ldl2.length; j++){
                res[i][j] = new HashSet<String>(ldl2[j].getKeys());
                res[i][j].retainAll(keys1);
                
            }
        }
        return res;
    }
    
    public  void switchAlleles(){
    	//EmissionState hap = ldl2[3].dataL.get("NA18515");
    //	EmissionState hap1 = ldl2[1].dataL.get("NA18515");
        for(int i=0; i<        this.switchedAlleles.size(); 
        i++){
           Info inf=  switchedAlleles.get(i);
           System.err.println("SWITCHING "+inf);
           ldl[inf.data_index].switchAlleles(inf.relative_position);
       //    PseudoDistribution di = ((HaplotypeEmissionState)ldl[1].dataL.get("NA18561")).emissions[0];
        //    System.err.println("h "+di);
        }
       // PseudoDistribution di = ((HaplotypeEmissionState)ldl[1].dataL.get("NA18561")).emissions[0];
        //System.err.println("h");
    }
    
    public static void compareDataPoints(Set<String> ck, Info inf1, Info inf2,  DataCollection ldl1, DataCollection ldl2){
    //int ind1 = (Integer) res[2];
    //int j1 = (Integer) res[3];
      
    for(Iterator<String> it = ck.iterator(); it.hasNext();){
        String st = it.next();
           int st1 = ((PhasedDataState) ldl1.data.get(st)).emissions[inf1.relative_position].fixedInteger();
           int st2 = ((PhasedDataState) ldl2.data.get(st)).emissions[inf2.relative_position].fixedInteger();
           if(st1!=st2)  try{
               throw new RuntimeException("inconsistent "+st);
           }catch(Exception exc){
               exc.printStackTrace();
           }
      
           
        }
        
    }
    
    public List<Info> switchedAlleles = new ArrayList<Info>();
   public Info[][] map;// = new HashMap<String, Map<Integer,Info>>(); 
   private boolean compare(Info[] m){
	   int i=0;
	   Info first = null;
	  
	   for(; i<m.length; i++){
		   if(m[i]!=null){
			   first = m[i];
			   if(first.strand!=null){
				   //so if anything has a strand it will be first
				   break;
			   }
		   }
	   }
	  
//	   List<Integer> toDelete = new ArrayList<Integer>();
	   
      outer: for(i=0; i<m.length; i++){
    	   try{
          if(m[i]!=null && m[i]!=first){
        	  Info curr = m[i];
        	  if(countAlleles(curr,first)>2){
        		//  toDelete.add(i);
    			  System.err.println(">2 alleles "+ldl[i].name+curr);
    			  return true;
        	  }
        	  if(m[i].strand==null && first.strand!=null){
        		  if(iscompl(curr.alleleA, curr.alleleB) || iscompl(first.alleleA,first.alleleB)
        				  ){
        			if(!this.ldl[0].strand_represents_ref) {
        			//  toDelete.add(i);
        			  System.err.println("ambiguous alleles "+ldl[i].name+curr);
        			  return true;
        				  //throw new RuntimeException("need strand for "+ldl[i].name+" "+curr);
        			  }
        		  }
        		  if(isoppstrand(first, curr) ){
        			
        			  this.ldl[curr.data_index].flipStrand(curr.relative_position);
        			  curr.flipStrand();
        			
        		  }
        	   }
        	/* 
        	 * note: removed this
        	 * if(isoppstrand(first,curr)){
        		 String rsid = this.snpid.get(i);
        		 try{
        		 throw new RuntimeException("should have flipped "+rsid);
        		 }catch(Exception exc){
        			 curr.flipStrand();
        		 }
        	 }*/
        	 // if(m[i].alleleA!=null && first.alleleA!=null && m[i].alleleB!=null && first.alleleB!=null){
        		  first.compare(m[i], this.switchedAlleles);
        	  //}
          }
    	   }catch(Exception exc){
    		 exc.printStackTrace();
    		  // exc.printS
    	   }
       }
       return false;
   }
private static List<Character> bases= Arrays.asList(new Character[] {'A','C','T','G','I','D',
																						//'.',
																						'N'});  
private static boolean[] cnts = new boolean[bases.size()];
   private int countAlleles(Info curr, Info first) {
	   Arrays.fill(cnts,false);
	   setTrue(cnts, bases.indexOf(curr.alleleA));
	   setTrue(cnts, bases.indexOf(curr.alleleB));
	   setTrue(cnts,bases.indexOf(first.alleleA));
	   setTrue(cnts,bases.indexOf(first.alleleB));
	int res = Constants.sum(cnts);
	return res;
   }

private void setTrue(boolean[] cnts2, int indexOf) {
if(indexOf>=0) cnts2[indexOf]=true;
	
}
private boolean isoppstrand(Info first, Info curr) {
	   {
		   if(first.alleleA==null && first.alleleB==null) return false;
	   Boolean matchA = first.alleleA==null || curr.alleleA==null  ? null :   first.alleleA.charValue()==curr.alleleA.charValue();
	   Boolean matchB = first.alleleB==null || curr.alleleB==null ? null :  first.alleleB.charValue()==curr.alleleB.charValue();
	   Boolean matchAB = first.alleleB==null || curr.alleleA==null ? null : first.alleleB.charValue()==curr.alleleA.charValue();
	   Boolean matchBA = first.alleleA==null || curr.alleleB==null  ? null : first.alleleA.charValue()==curr.alleleB.charValue();
	   if((matchB !=null && matchA !=null && matchA!=matchB )|| (matchAB!=null && matchBA!=null && matchBA!=matchAB)){
		   throw new RuntimeException("!!");
	   }
	   if(matchA!=null && matchA || matchAB!=null && matchAB || matchB!=null && matchB || matchBA!=null && matchBA) {
		   return false;
	   }
	   
	   }
	   {
		   Boolean matchA = first.alleleA==null || curr.alleleA ==null ? null : first.alleleA.charValue()==DataCollection.compl(curr.alleleA.charValue());
		   Boolean matchB = first.alleleB==null || curr.alleleB==null  ? null : first.alleleB.charValue()==DataCollection.compl(curr.alleleB.charValue());
		   Boolean matchAB = first.alleleB==null || curr.alleleA==null  ? null : first.alleleB.charValue()==DataCollection.compl(curr.alleleA.charValue());
		   Boolean matchBA = first.alleleA==null || curr.alleleB==null  ? null : first.alleleA.charValue()==DataCollection.compl(curr.alleleB.charValue());
		   if((matchB !=null && matchA !=null && matchA!=matchB )|| (matchAB!=null && matchBA!=null && matchBA!=matchAB)){
			   throw new RuntimeException("!! "+first+"\n"+curr.alleleB);
		   }
		   if(matchA!=null && matchA || matchAB!=null && matchAB || matchB!=null && matchB || matchBA!=null && matchBA) {
			   return true;
		   }
		   
	      
	   }
	   if(Constants.CHECK)
	   throw new RuntimeException("Insufficient information "+first+"\t"+curr+"\n"+
			   this.ldl[first.data_index].name+"\t"+this.ldl[curr.data_index].name);
	   else{
		   if(first.alleleA!=null && curr.alleleA!=null ) return !first.alleleA.equals(curr.alleleA);
		   else if(first.alleleB!=null && curr.alleleB!=null ) return !first.alleleB.equals(curr.alleleB);
		   else{
			System.err.println("Insufficient information "+first+"\t"+curr+"\n"+
					   this.ldl[first.data_index].name+"\t"+this.ldl[curr.data_index].name);   
			return false;
		   }
		   
	   }
	   
}
private boolean iscompl(Character alleleA, Character alleleB) {
	if(alleleA!=null && alleleB!=null && DataCollection.compl(alleleA).charValue() == alleleB.charValue()) return true;
	else return false;
}
private static void compareCommonStates(Map<Integer, Info> m, Set<String>[][]commonKeys, DataCollection[] ldl2){
       Iterator<Integer> it = m.keySet().iterator();
       Integer firstKey = it.next();
       Info first = m.get(firstKey);
       while(it.hasNext()){
           Info info = m.get(it.next());
           compareDataPoints(commonKeys[first.data_index][info.data_index],first, info, ldl2[first.data_index], ldl2[info.data_index] );
       }
      
   }
   public Boolean probeOnly(int i) {
	   Boolean probeOnly = true;
		for(int k=0; k<this.map[i].length; k++){
			if(map[i][k]!=null){
				Boolean po = this.ldl[k].probeOnly(map[i][k].relative_position);
			    if(po!=null && !po) return false;
			}
		}
		return true;

	}
   
  
   
  protected Info[][] getMap(DataCollection[] ldl2,Map<String, Integer> idToPos ){
	  this.map = new Info[this.loc.size()][ldl.length];
      
      for(int i=0; i<ldl2.length; i++){
    	//  if(Constants.offsetSamples() && i>1 && overlaps(ldl2[i].loc, ldl2[])
          for(int j=0; j<ldl2[i].loc.size(); j++){
              Character A_all =ldl2[i].alleleA.size()>0 ?  ldl2[i].alleleA.get(j) : null;
              Character B_all = ldl2[i].alleleB.size()>0 ?  ldl2[i].alleleB.get(j) : null;
              if(A_all!=null && (A_all.charValue()=='N' || A_all.charValue()=='-' || A_all.charValue()=='.')) A_all = null;
              if(B_all!=null && (B_all.charValue()=='N' || B_all.charValue()=='-' || B_all.charValue()=='.')) B_all = null;
              Boolean strand = ldl2[i].strand.size()>0 ? ldl2[i].strand.get(j) : null;
              String id = ldl2[i].snpid.get(j);
              Integer loc = ldl2[i].loc.get(j); 
              double maf = 0;//ldl2[i].getMaf(j)
              Info inf = new Info(id, B_all, A_all,i, j, loc, maf, strand,ldl2[i].probeOnly(j));
           //   inf.
              Integer x =  idToPos.get(id);
              if(x!=null){
              	map[x][i]= inf;
              	if(Constants.offsetSamples() && i>=1 && map[x][0]!=null) {
              		throw new RuntimeException("!!");
              	}
              }
          }
      }
//      Info[] inf = map[717];
      return map;
  }
   
  public void drop(List<Integer> toDrop, boolean setAsMissing){
	  super.drop(toDrop, setAsMissing);
	  List<Info[]> mapl= new ArrayList<Info[]>(Arrays.asList(map));
	  for(int i=toDrop.size()-1; i>=0 ; i--){
		  mapl.remove(toDrop.get(i).intValue());
	  }
	  map = mapl.toArray(new Info[0][]);
  }
  
   /*[minor major index pos_in_index, loc]*/
    protected void compareAlleles(Map<String, Integer> idToPos, List<Integer> toDelete, List<Integer>[]toDel ) {
    	  map = getMap(ldl, idToPos);
    	  
    	  
       int[][] eq = Constants.equaliseGroup();
     
       if(eq!=null){
    	   boolean modified = false;
    	   int[]eqlens = new int[eq.length];
    	   for(int i=0; i<eqlens.length; i++){
    		   eqlens[i] = eq[i].length;
    	   }
          List<Integer>[] todrop = new List[ldl.length];
          List<Info>[] toinsert = new List[ldl.length];
          List<Integer> todropall = new ArrayList<Integer>();
       for(int k=0; k<todrop.length; k++){
    	   todrop[k] = new ArrayList<Integer>();
    	   toinsert[k] = new ArrayList<Info>();
       }
       int len = ldl.length;
       int eqlen = eq.length;
       for(int i=0; i<map.length; i++){
    	   //System.err.println(i);
    	   Info[] mapi = map[i];
    	inner: for(int j=0; j<eqlen ;j++){
    		if(eqlens[j]<=1) continue inner;
    		boolean hasOne = false;
    		boolean hasAll = true;
    		boolean hasFirst = false;
    		Info nonNull = null;
    		for(int k=0; k<eqlens[j]; k++){
    			boolean notnull = mapi[eq[j][k]]!=null;
    			if(notnull) nonNull = mapi[eq[j][k]];
    			if(k==0 && notnull) hasFirst = true;
    			hasOne = hasOne || notnull;
    			hasAll = hasAll && notnull;
    		}
    		if(Constants.equaliseGroupMode()==2){
    			if(!hasFirst){
    				double baf = 0;
    				int ind = -1;
    				for(int k=1; k<eq[j].length; k++){
		    			if(mapi[eq[j][k]]!=null){
		    				baf = Math.max(baf,ldl[eq[j][k]].baf(mapi[eq[j][k]].relative_position));
		    				ind = eq[j][k];
		    			}
    				}
    				if(baf <= Constants.excludeBafThresh(ind) || baf > (1-Constants.excludeBafThresh(ind))){
	    			for(int k=1; k<eq[j].length; k++){
		    			if(mapi[eq[j][k]]!=null){
		    				
		  				  todrop[eq[j][k]].add(mapi[eq[j][k]].relative_position);
		  				  mapi[eq[j][k]]=null;
		  				  modified = true;
		  				}
	    			}
    				}
    			}else{
    				for(int k=1; k<eq[j].length; k++){
    					if(mapi[eq[j][k]]==null){
	    				  toinsert[eq[j][k]].add(nonNull);
	    				  modified = true;
	    				}
    				}
    			}
    		}
    		else{
    			if(hasOne!=hasAll){
    		
    			for(int k=0; k<eq[j].length; k++){
    				 if(Constants.equaliseGroupMode()==0){
	    				if(mapi[eq[j][k]]!=null){
	    				  todrop[eq[j][k]].add(mapi[eq[j][k]].relative_position);
	    				  mapi[eq[j][k]]=null;
	    				  modified = true;
	    				}
    				}else{
    					if(mapi[eq[j][k]]==null){
  	    				  toinsert[eq[j][k]].add(nonNull);
  	    				  modified = true;
  	    				}
    				}
    			}
    			boolean dropinAll = true;
    			for(int k=0; k<len; k++){
    				if(mapi[k]!=null ) dropinAll = false;
    			}
    			if(dropinAll){
    				modified = true;
    				todropall.add(i);
    			}
    		}
    		}
    	}
       }
       for(int j=0; j<len; j++){
    	  if(todrop[j].size()>0) ldl[j].drop(todrop[j], false);
    	  if(toinsert[j].size()>0) if(toinsert[j].size()>0) ldl[j].insert(toinsert[j]);
       }
       if(modified){
       this.map = getMap(ldl,getPosLoc(ldl, Constants.includeSNPS()));
       }
       }
      // List<Integer> toDelete = new ArrayList<Integer>();
        for(int i=0; i<map.length; i++){
        	
           if( compare(map[i])){
        	   toDelete.add(i);
        	   for(int k=0; k<map[i].length; k++){
        		   if(map[i][k]!=null){
        			   toDel[k].add(map[i][k].relative_position);
        		   }
        	   }
           }
        }
        if(Constants.allowFlips()){
        this.switchAlleles();
        }
       
     }
       
    private void print(double[] probs1) {
        Double[] d = new Double[probs1.length];
        StringBuffer sb = new StringBuffer();
       for(int i=0; i<probs1.length; i++){
           d[i] = probs1[i];
           sb.append("%5.3f \t");
       }
       System.err.println(String.format(sb.toString(), d));
        
    }
    private static char transform(char c){
        char a = Character.toLowerCase(c);
        if(a=='a' || a=='t') return 'a';
        else if(a=='g' || a=='c') return 'g';
        else {
            throw new RuntimeException("!!");
        }
    }
    
    static boolean  check(char c, char charValue) {
       return c==charValue;
        
    }
   
   
    
   public final DataCollection[] ldl;
    HaplotypeEmissionState[] tmp;
   
   // public Set<String> perfectMatch = new HashSet<String>();
    @Override
    public HaplotypeEmissionState createEmissionState(String key, int no_copies){
        if(tmp==null){
            tmp = new HaplotypeEmissionState[ldl.length];
          
        }
        for(int i=0; i<tmp.length; i++){
            tmp[i] = (HaplotypeEmissionState)ldl[i].dataL.get(key);
        }
      Double[] phenD =  ((MergedPhenotypes) this.pheno).getRecodedValues(tmp);
          int no_copy = -1;
        Class clazz=null;
        HaplotypeEmissionState match = null;
        int matchCount=0;
          for(int i=0; i<ldl.length; i++){
              if(ldl[i].containsKey(key)){
            	  EmissionState val =ldl[i].dataL.get(key);
            	  matchCount++;
            	  if(ldl[i].loc.size()==this.loc.size()){
            		 match = (HaplotypeEmissionState)val;
            	  }
                //  Double[] phenD_i = ldl[i].dataL.get(key).phenValue();
                  clazz = val.getClass();
                 if(no_copy <0) no_copy = no_copies;
                 else if(no_copy!=no_copies) throw new RuntimeException("!!");
                 
                 
                 
                 
              }
          }
          try{
        	  if(match!=null && matchCount==1){
        		//  perfectMatch.add(key);
        		  return match;
        	  }
        	  EmissionStateSpace stSp = Emiss.getSpaceForNoCopies(no_copy);
           HaplotypeEmissionState res =  
        	   //clazz.equals(ASCNEmissionState.class) ? 
        	 //(HaplotypeEmissionState)  clazz.getConstructor(new Class[] {String.class, int.class, EmissionStateSpace.class, short.class}).newInstance(
        		//	   new Object[] {key, length, stSp[no_copy-1], (short)-1}) : 
            (clazz.equals(BackgroundEmissionState.class) ?
            		  new BackgroundEmissionState(key, length, stSp.size(), stSp,(short)-1):
            		new HaplotypeEmissionState(key, this.length,  stSp, (short) -1));
           res.setPhenotype(phenD);
           return res;
          }catch(Exception exc){
        	  exc.printStackTrace();
          }
          return null;
       
    }
    
   
	@Override
	public List<AssociationCalculator>[][] getArmitageCalculator() {
		if(Constants.plotMerged()) return super.getArmitageCalculator();
		if(ac==null){
			SortedSet<Integer> ploidy = getPloidy();
			ac = new List[ploidy.last()][ldl.length];
		
		for(int i=0; i<this.ldl.length; i++){
			List<AssociationCalculator>[][] assoc2 = this.ldl[i].getArmitageCalculator();
			if(assoc2!=null){
			for(int j=0; j<assoc2.length; j++){
				
					ac[j][i] = assoc2[j][0];
				
			}
			}
			
		 }
			/*for(int j=0;j<ac.length; j++){
			    ac[j][ldl.length] = getMergedCalc(ac[j]); 
			}*/
		}
		return ac;
	}
   /* ArmitageCalculator ac;
	@Override
    public List<AssociationCalculator>[] getArmitageCalculator() {
		// TODO Auto-generated method stubt
		
		if(this.ac == null){
			ac = super.getArmitageCalculator();
			if(allnull(ac)){
				SortedSet<Integer> ploidy = getPloidy();
				ac = new ArmitageCalculator[ploidy.last()];
				
				 //String dirF = Constants.baseDir;
		        // if(dirF.equals(".")){
		         	String dirF = System.getProperty("user.dir");
		        // }
				try{
						for(Iterator<Integer> it = ploidy.iterator(); it.hasNext();){
							int ploi = it.next();
							ac[ploi-1] = new ArmitageCalculator(this, ploi);
						}
				}catch(Exception exc){
					exc.printStackTrace();
				}
				
				
			}
			
		}
		return ac;
	}*/
	private boolean allnull(AssociationCalculator[] ac) {
	for(int i=0; i<ac.length; i++){
		if(ac[i]!=null) return false;
	}
	 return true;
}
	public void removeKey(String key,boolean states){
		if(Constants.writeMergedAverageFiles()) super.removeKey(key, states);
		else {
			for(int i=0; i<ldl.length; i++){
			if(ldl[i].containsKey(key)){
			ldl[i].removeKey(key, false);
			}
		}
		}
	}
	@Override
	protected String getTypes(int i) {
		StringBuffer sb = new StringBuffer();
		for(int j=0; j<this.map[i].length; j++){
			if(map[i][j]!=null){
				sb.append(ldl[j].name);
			}
		}
		return sb.toString();
	}
  
	
	
    
}
