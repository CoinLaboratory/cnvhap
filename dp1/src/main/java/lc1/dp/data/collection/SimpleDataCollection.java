package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.logging.Logger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;
import lc1.dp.states.CachedEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.ensj.PedigreeDataCollection;
import lc1.stats.IntegerDistribution;
import lc1.stats.MixtureDistribution;
import lc1.stats.PseudoDistribution;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Compressor;
import lc1.util.Constants;

public class SimpleDataCollection extends DataCollection {
    public SimpleDataCollection(){
        
    }
    protected SimpleDataCollection(DataCollection dat){
      super(dat);
     
       
    }
    
    public void calcIndiv() {
    	this.indiv = new ArrayList<String>(this.data.keySet());
    	
    }
    public boolean hasIntensity(int i) {
    	return false;
    }
    
    public  static void readIllumina (File f, SimpleDataCollection sdt) throws Exception{
        List<Integer >loc1 = readPosInfo(f,4, true);
        int read; int end;
       if(Constants.mid()[0][0]>0){
           int[] mid = Constants.mid()[0];
           int mid_index = IlluminaRDataCollection.firstGreaterThan(loc1, mid[0]);
           int mid_index1 = IlluminaRDataCollection.firstGreaterThan(loc1, mid[1]);
           read = mid_index -(int) Math.round( (double)Constants.restrict()[0]/2.0);
           end  = mid_index1 +(int) Math.round( (double)Constants.restrict()[0]/2.0);
       }
       else{ 
           read = IlluminaRDataCollection.firstGreaterThan(loc1,  Constants.offset());
           end  = Math.min(IlluminaRDataCollection.getEndIndex(loc1), read+Constants.restrict()[0]);
         
       }
       read = Math.max(0,read);
       end =  Math.min(end,loc1.size());
        sdt.loc = new ArrayList<Integer>(loc1.subList(read, end)) ;
        sdt.snpid = new ArrayList<String>();
       sdt.length = sdt.loc.size();
       //this.make();
//       loc = new ArrayList<Integer>(loc.subList(0, length));
      //  File dir1 = new File(Constants.getDirFile());
        BufferedReader br = new BufferedReader(new FileReader(f));
        CompoundEmissionStateSpace stSp = Emiss.getEmissionStateSpace(1);
        CompoundEmissionStateSpace stSp1 = Emiss.getEmissionStateSpace(1)  ;
        EmissionStateSpaceTranslation trans = new EmissionStateSpaceTranslation(stSp, stSp1, true);
          
        String[] str = br.readLine().split("\\t");
       // System.err.println(str.length);
       // if(true)System.exit(0);
        int len = str.length-5;
        String[] indv = new String[len / 3];
        PIGData[] data = new PIGData[indv.length];
        int i1 =0;
       
       StringBuffer sb = new StringBuffer();
       for(int i=0; i<sdt.loc.size(); i++){
           sb.append("%5.3f ");
       }
    
        for(int i=5; i<str.length; i+=3){
            indv[i1] = str[i].substring(0, str[i].indexOf(".GType"));
          
            data[i1] = SimpleScorableObject.make(indv[i1], sdt.loc.size(), stSp, (short)-1);
            i1++;
        }
        for(int i=0; i<read; i++){
            br.readLine();
        }
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for(int i=0; i<sdt.length; i++){
            String[] st = br.readLine().split("\\t");
            int loc = Integer.parseInt(st[4]);
            sdt.snpid.add(st[1]);
            if(loc!=sdt.loc.get(i)) throw new RuntimeException("!!");
            for(int j=0; j<indv.length; j++){
               
                String geno = st[5+j*3];
                if(geno.equals("NC")) geno = "NN";
                double b = Double.parseDouble(st[6+j*3]);
                double r = Double.parseDouble(st[7+j*3]);
                if(r < min) min = r;
                if(r>max) max = r;
                geno = geno.replaceAll("N", "").replaceAll("_", "");
                int comp = trans.bigToSmall(stSp.getFromString(geno));
                data[j].addPoint(i, stSp1.get(comp));
            }
        }
     
       
        sdt.length = sdt.size();
        for(int i=0; i<data.length; i++){
              sdt.data.put(data[i].getName(), data[i]);
        }
      //  sdt.calculateMaf(false);
        br.close();
    }
 
    
    public static String[]  readAffy(File dir1, String name,boolean quotes, SimpleDataCollection res, int noCopies) throws Exception{
       // SimpleDataCollection res = new SimpleDataCollection();
     EmissionStateSpace emStSp =  Emiss.getEmissionStateSpace(noCopies-1); //source
      EmissionStateSpace emStSp1 = Emiss.getEmissionStateSpace(noCopies-1); //target
      EmissionStateSpaceTranslation trans = new EmissionStateSpaceTranslation(emStSp, emStSp1, true);
        File f = new File(dir1, name);
        res.loc = DataCollection.readPosInfo(f, 4, true);
        int read = IlluminaRDataCollection.firstGreaterThan(res.loc, Constants.offset());
        res.loc= res.loc.subList(read, Math.min(read+Constants.restrict()[0], res.loc.size()));
        BufferedReader br = DataCollection.getBufferedReader(f);
        String st = br.readLine().trim();
        boolean header = st.toUpperCase().indexOf("CHR")>=0;
        PIGData[] dat=null;
        String[] indiv;
        if(header){
           String[] str = st.trim().split("\\s+");
           indiv = new String[str.length-2];
           System.arraycopy(str, 2, indiv, 0, indiv.length);
            st = br.readLine();
        }
        else indiv = readIndiv(new File(dir1, "indiv.txt"));
        for(int i=0; i<read; i++){
            st = br.readLine();
        }
        
        for(int ik1=0; st!=null && ik1<res.loc.size(); ik1++){
           String[] str =  st.trim().replaceAll("No Call", "NN").split("\\s+");
           if(dat ==null){
               dat = new PIGData[str.length-2];
               for(int i=0; i<dat.length; i++){
                   dat[i] = SimpleScorableObject.make(indiv[i], 10, emStSp,(short)-1);
               }
           }
           for(int i=2; i<str.length; i++){
               String sti = str[i];
           //    sti = sti.replace('_', 'N');
               if(quotes) sti = sti.substring(1, sti.length()-1);
           char[] c = sti.toCharArray();
            Arrays.sort(c);
            //  sti = new String(c);
             StringBuffer sb = new StringBuffer();
               for(int ik=0; ik<c.length; ik++){
                   if(c[ik]!='_' && c[ik]!='-'){
                       sb.append(c[ik]);
                   }
               }
               sti = sb.toString();
               
               int comp = trans.bigToSmall(emStSp.getFromString(sti));
               dat[i-2].addPoint( ik1,emStSp1.get(comp));
           }
            st = br.readLine();
        }
        for(int i=0; i<dat.length; i++){
            res.data.put(dat[i].getName(), dat[i]);
        }
        res.length = dat[0].length();
        return indiv;
        
    }
    private static ComparableArray getArray(String sti, int noC)throws Exception {
        if(noC==2){
        if(sti.length()==0) sti = "NN";
        if(sti.length()==1) sti = new String(sti+"N");
        int half =(int)  Math.floor((double)sti.length()/2.0);
             return    ComparableArray.make(Emiss.get(sti.substring(0, half)),Emiss.get(sti.substring(half)) );
        }
        else if(noC==1){
            if(sti.length()==0) sti = "N";
        
                  return  ComparableArray.make(Emiss.get(sti) );
                  
                  
        }
        else throw new RuntimeException("!!");
    }
    public void selectMale(boolean male) {
        
        Set<String> alias = new HashSet<String>();
         for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
             PIGData nxt = it.next();
             if(nxt.noCopies()==1 || nxt.noCopies()==3 ){
                 if(male) alias.add(nxt.getName());
             }
             else{
                 if(!male) alias.add(nxt.getName());
             }
         }
           this.restricToAlias(alias);
       }
    


    
    public void selectMaleTrios() {
        Set<String> alias = new HashSet<String>();
         for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
             PIGData nxt = it.next();
             ComparableArray arr = (ComparableArray) ((ComparableArray)nxt.getElement(0)).get(2);
             if(arr.size()==1){
                 alias.add(nxt.getName());
             }
         }
           this.restricToAlias(alias);
       }
   
    public int size() {
        return this.data.size();
     }
  
    public List<String[]> getPairings( double[] fraction, List<String>keys){
        int no_copies = Constants.sample(fraction)+1;
        List<String[]>res = new ArrayList<String[]>();
        while(keys.size()>=no_copies){
            String[] keys1 = new String[no_copies];
            res.add(keys1);
            for(int i=0; i<no_copies; i++){
                int i1 = //obj.size()-1;
                    no_copies==1 ?  0 : 
                   Constants.nextInt(keys.size());
               keys1[i] = keys.remove(i1);
            }
            no_copies = Constants.sample(fraction)+1;

        }
        return res;
    }
   
   /* public   void transform ( double[] fraction, List<String>keys,  String join){
    	this.dataL.clear();
        transform(getPairings(fraction, keys).iterator(), join);
    }
    public void transform(List<String> names){
        List<String[]> res = new ArrayList<String[]>();
        for(Iterator<String> it = names.iterator(); it.hasNext();){
            res.add(it.next().split("\\."));
        }
        transform(res.iterator(), ".");
    }
    
    public void dropSingle(){
    	List<String> toD = new ArrayList<String>();
    	for(Iterator<String> it = this.data.keySet().iterator(); it.hasNext();){
    		String st = it.next();
    		PhasedDataState state = (PhasedDataState) this.data.get(st);
    		if(state.getEmissionStateSpace().noCopies()==1) toD.add(st);
    	}
    	this.dropIndiv(toD.toArray(new String[0]));
    }*/
   /* public   void transform (Iterator<String[]> pairs, String join){
      //  Map<String, PIGData> newL = new HashMap<String, PIGData>();
      
        while(pairs.hasNext()){
           
            String[] keys1 = pairs.next();
            PhasedDataState[] unit = new PhasedDataState[keys1.length];
            SortedMap<Integer, Integer>[] recSites= new SortedMap[keys1.length];
            String[] parents = new String[keys1.length];
            for(int i=0; i<unit.length; i++){
              
                String str = keys1[i];
                parents[i] = ped!=null && ped.mother!=null ? ped.mother.remove(str) : null;
                unit[i] = (PhasedDataState) data.remove(str);
                if(recSites!=null && this.recSites!=null){
                    SortedMap<Integer, Integer>[] tmp =  this.recSites.remove(str);
                    if(tmp!=null) recSites[i] =tmp[0];
                }
            }
            PIGData pid = SimpleScorableObject.make(unit, unit[0].noCopies()==1, join, this.stSp[1]);
            if(parents[0]!=null){
                ped.mother.put( pid.getName(), parents[0]);
            }
            if(parents.length>1 && parents[1]!=null){
                ped.father.put( pid.getName(), parents[1]);
            }
            if(recSites[0]!=null){
                this.recSites.put(pid.getName(), recSites);
            }
            data.put(pid.getName(), pid);
        }
      
       // this.data = newL;
    }*/
    
    public void rename(){
        int ik=0;
        Map<String, PIGData> newL = new HashMap<String, PIGData>();
        for(Iterator<String> it =data.keySet().iterator(); it.hasNext();ik++){
            String key = it.next();
            PIGData nxt = this.data.get(key);
            nxt.setName("pair"+ik);
            newL.put(nxt.getName(), nxt);
           
        }
        data = newL;
    }
    
    public void makeTrioData(){
        ped = new PedigreeDataCollection();
        List<String> l = new ArrayList<String>(data.keySet());
        int ik=0;
        for(Iterator<String> it =l.iterator(); it.hasNext();ik++){
            String key = it.next();
            PIGData dat_i = this.data.get(key);
           SortedMap<Integer, Integer> rec_sites = getRecSites(Constants.hotspot(2), dat_i.getName()+"r");
           int startingIndex = Constants.nextInt(2);
            System.err.println("sw "+ik+" " +rec_sites);
            PIGData child = dat_i.recombine(rec_sites, startingIndex);
            ped.mother.put(child.getName(), dat_i.getName());
            data.put(child.getName(), child);
            this.recSites.put(child.getName(), new SortedMap[] {rec_sites});
        }
    }
    
    
    /*
        Map<String, PIGData> newL = new HashMap<String, PIGData>();
        Map<String, PIGData> newRec = new HashMap<String, PIGData>();
        int no_copies = 2;
        List<String> keys = new ArrayList<String>(data.keySet());
        while(data.size()>=no_copies){
            CSOData[] unit = new CSOData[3];
            CSOData[] unit1 = new CSOData[2];
            CSOData[] recSites1 = new CSOData[2];
            int[] i1 = new int[2];
            String[] keys1 = new String[2];
            for(int i=0; i<2; i++){
                i1[i] =  Constants.nextInt(data.size());
                keys1[i] = keys.remove(i1[i]);
                unit[i] =  data.remove(keys1[i]);
                unit1[i] =l.remove(keys1[i]);
                recSites1[i] = recSites.remove(keys1[i]);
            }
            unit[2] = new PIGData(unit1, true, ".");
            PIGData pid = new PIGData(unit, false, ";");
            PIGData recs = new PIGData(recSites1, true,".");
            newL.put(pid.getName(), pid);
            newRec.put(recs.getName(), recs);
          //  no_copies = Constants.sample(fraction)+1;
        }
        this.data = newL;
        this.recSites = newRec;
       // Set<Integer>[] set =  newRec.get(0).getSwitches();
       // this.makeDataIndex();
        
    }*/
    
   
  
   
    public Iterator<PIGData> iterator() {
        return this.data.values().iterator();
        }
    public SimpleDataCollection(int length) {
        super(length);
        // TODO Auto-generated constructor stub
    }
  
   
   
   
    /*  private List<PIGData> data_bu;
    private List<PIGData> trio_bu;*/
    
    public SimpleDataCollection(List<PIGData> name) {
        this(name.get(0).length());
        for(Iterator<PIGData> it = name.iterator();  it.hasNext();){
            PIGData nxt = it.next();
            this.data.put(nxt.getName(), nxt);
        }
    }
    
   
   
    
  /*  public void resetToTrio(){
        if(trio_bu!=null){
            this.data = trio_bu;
        }
    }*/
    public void applyMedianCorrection(ZipFile zf) throws Exception{
    	
    }
    public  void applyVarianceThreshold(ZipFile zf, boolean standardise) throws Exception{
    	
    }
   public SimpleDataCollection(File file,short index, int no_copies, 
		   int[][] mid,  File bf,Collection<String> snpidrest) throws Exception{
       super(file, index, no_copies, mid, bf,snpidrest);
      
    }
   
    public  SimpleDataCollection clone(){
        return new SimpleDataCollection(this);
    }
    
    @Override
    protected double calcBaf(List<String> l) {
    	//int ploidy = Constants.backgroundCount(index);
    	//   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	int gl_ind = this.headsLowerCase.indexOf("gl");
    	boolean pl = false;
    	if(gl_ind < 0){
    		pl = true;
    		gl_ind = this.headsLowerCase.indexOf("pl");
    	}
    	double cnt=0;
    	double bcnt =0;
    	if(gl_ind<0){
    		int gt_ind = this.headsLowerCase.indexOf("gt");
    		if(gt_ind <0) return super.calcBaf(l);
    		else{
    			for(int k=0; k<Math.min(l.size(), Constants.includeFirstInBAF()); k++){
    				
    	    		String str = l.get(k).trim();//.split("\\s+")
    	    		if(str.equals(".")) continue;
    	    		else{
    	    			
    	    			char[]  st = str.split("\\s+")[gt_ind].replaceAll("/", ",").replaceAll("0","A").replaceAll("1", "B").replaceAll("2", "BB").toCharArray();;
    	    			cnt+=st.length;
    	    			for(int kk=0; kk<st.length ;kk++){
    	    				if(st[kk]=='B') bcnt++;
    	    			}
    	    			
    	    		 
    	    		
    	    		
    	    		}
    	    	}
    		}
    	}
    	else{
    	//double[] d = new double[3];
    	for(int k=0; k<Math.min(l.size(), Constants.includeFirstInBAF()); k++){
    		String str = l.get(k).trim();//.split("\\s+")
    		if(str.equals("./.") || str.equals(".")) continue;
    		else{
    		String[] st = str.split("\\s+")[gl_ind].split(",");
    		if(st.length==1 && (st[0].equals("./.") || st[0].equals("."))) continue;
    		for(int i=0; i<st.length; i++){
    			double val =Double.parseDouble(st[i]); 
    			double d = pl ? Math.pow(10, val/-10.0): Math.pow(10, val);
    		 cnt+=d;
    		    bcnt+= i==0 ? 0 :( i==1 ? d/2.0 :d);
    		}
    		
    		}
    	}
    	
    	}
    	return bcnt/cnt;
    }
  
   // double[] dist;
   
    public Set<String> applyLoess(ZipFile zf, Map<String, String>platem, boolean loess_, boolean gc_) throws Exception{
    	return new HashSet<String>();
    }
   
   
    public void join(String[] str){
    	if(str[0].equals("null")) return;
    //	if(true) throw new RuntimeException("!!");
    	List<String> toRemove = new ArrayList<String>();
    	Map<String, List<String>> m = new HashMap<String, List<String>>();
    	if(str.length==1){
    		String sti = str[0];
    		for(Iterator<String> it = this.dataL.keySet().iterator();it.hasNext();){
    			String st = it.next();
    			String[] stri =  st.split(sti);
    			if(stri.length!=2){
    				toRemove.add(st);
    				continue;
    			}
    			String key =stri[0];
    			List<String> res = m.get(key);
    			if(res ==null) m.put(key, res = new ArrayList<String>());
    			res.add(st);
    			if(res.size()>2){
    				throw new RuntimeException("!!");
    			}
    		}
    	}
    	else{
    		for(int i=0; i<str.length; i+=2){
    			String key = str[i]+"_"+str[i+1];
    			m.put(key, Arrays.asList(new String[] {str[i], str[i+1]}));
    		}
    	}
    	for(int i=0; i<toRemove.size(); i++){
    		dataL.remove(toRemove.get(i));
    	}
    	join(m);
    }
  
    
    public void join(Map<String, List<String>> tojoin){
    	for(Iterator<String> keys = tojoin.keySet().iterator(); keys.hasNext();){
    		String key = keys.next();
    		
    		List<String> vals = tojoin.get(key);
    		EmissionState[] state = new EmissionState[vals.size()];
    		int sum=0;
    		for(int i=0; i<state.length; i++){
    			state[i] = (HaplotypeEmissionState) dataL.remove(vals.get(i));
    			sum+=state[i].noCop();
    		}
    		CompoundEmissionStateSpace ces = Emiss.getSpaceForNoCopies(sum);
    		PairEmissionState pes = new PairEmissionState(Arrays.asList(state), ces, true);
    		CachedEmissionState ces_ = new CachedEmissionState(pes, ces.size());
    		  PhasedDataState emst = SimpleScorableObject.make(key,loc.size(), null, this.index);
             ces_.refreshSiteEmissions();
    		emst.setNoCop(sum);
    		emst.emissions = ces_.emissions;
    		for(int i=0; i<emst.emissions.length; i++){
    			if(emst.emissions[i] instanceof SimpleExtendedDistribution){
    				Logger.global.info("h");
    			}
    		}
             dataL.put(key,emst );
    	}
    		this.indiv = new ArrayList<String>(dataL.keySet());
    }
    
    
   
    PseudoDistribution[] softened;
   @Override
    public Boolean  process(String indiv,  String[] header,  String[] geno, int ii, int ploidy, double[] missing){
    	
    	if(header[0].equals("CNV Genotype")) {
    		if(geno[4].equals("NaN")){
    			 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	    		dataSt.emissions[ii] =stsp.getHWEDist1(null);
    	    	   return false;
    		
    		}
    		return this.processQuanti(indiv, header, geno, ii, ploidy);
    	}
    	else if(header[0].equals("cn")){
 		 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
 		double d = Double.parseDouble(geno[0]);
 		int v = (int) Math.round(d);
 		double prob = 1-Math.abs(v-d);
		dataSt.emissions[ii] =process(v,  prob,  1.0, ploidy, "A");
		return true;
 		}
    	else if(header[0].equals("call")) {
    		if(geno[0].indexOf("N")>=0){
   			 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
   	    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
   	    		dataSt.emissions[ii] =stsp.getHWEDist1(null);
   	    	   return false;
   		
   		}
    		
   		 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
   		 CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
   		 String st1 =
			 geno[0].equals("H") ? "AB" : geno[0].equals("F") ? "AA" : geno[0].equals("T") ? "BB" : null;
		
					 int ind = trans(st1, stsp, null, null,index) ;
		 //if(dataSt.getName().equals("91169-6")){
		//	 System.err.println("h");
		 //}
		 dataSt.emissions[ii] = stsp.getIntDist(ind);
		 int nob = stsp.getBCount(ind);
		 int noa = stsp.getCN(ind)- nob;
		 missing[2]+=noa;
		 missing[3]+=nob;
		  return false;
    	}
    	else if(header.length>1 && header[1].equals("call")) {
    		if(geno[0].equals("null") || geno[1].equals("Probability")){
    			 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	    		dataSt.emissions[ii] =stsp.getHWEDist1(null);
    	    	   return false;
    		
    		}
    		 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    		 CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    		/* if(geno[1].replaceAll("DEL", "").replaceAll("[ACGT]","").replaceAll("\\.","").length()>0){
    			 //System.err.println("h");
    			 throw new RuntimeException("!!");
    		 }*/
    		 
    		 geno[1] = 
    			
    			 geno[1].equals("DEL") ? "BB" : geno[1].indexOf("DEL")>=0 ? "AB" : "AA";
    	//	 System.err.println(geno[1]+" "+str);
    		int ind = trans(geno[1], stsp, null, null,index) ;
    		 dataSt.emissions[ii] = stsp.getIntDist(ind);
    		 int nob = stsp.getBCount(ind);
    		 int noa = stsp.getCN(ind)- nob;
    		 missing[2]+=noa;
    		 missing[3]+=nob;
    		 return false;
    		//return true;
    	}
    	int noB=-1;
    //	if(header[0].equals("Allele1 - Forward")){
    //		return this.processForward(indiv,header,geno,ii,ploidy);
    //	}
    	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
    	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
    	   int i = ii;
    	try{
        	
        	
          //  boolean doneGeno = false;
        //    EmissionStateSpace stsp = this.stSp[this.no_copies-1];
        	   boolean deletionIsNa = Constants.useDeletion(i);
        	if(deletionIsNa) throw new RuntimeException("!!");
           
           // PhasedDataState data =(PhasedDataState)  this.data.get(indiv);
         
            for(int k=0; k<header.length; k++){
            	String hk =header[k].toLowerCase(); 
            if(hk.equals("gt")&& dataSt.emissions[i]==null){
            	if(k>= geno.length || geno[k].equals("./.")|| geno[k].equals(".") || geno[k].indexOf('.')>=0) {
            		//dataSt.emissions[i] = stsp.getZeroDistProb((int)Constants.coverage(index), ploidy);
            		dataSt.emissions[i] = stsp.getHWEDist1(null);
            		return false;
            	}
            	String st =  geno[k].replaceAll("/", ",").replaceAll("0","A").replaceAll("1", "B").replaceAll("2", "BB");;
            	if(ploidy==1 && st.length()>1  && st.indexOf("A")>=0 && st.indexOf("B")>=0){
            		dataSt.emissions[i] = stsp.getHWEDist1(null);
            		return false;
            	}else{
            	int i1 = 	  stsp.getHapl(ploidy==1 ? st.substring(0,1): st);
            	dataSt.emissions[i] = new IntegerDistribution( i1,stsp.getCN(i1),stsp.getBCount(i1));
            	dataSt.emissions[i].setEmStSp(stsp);
            	
            	}
            }
            else if(hk.equals("gq") && dataSt.emissions[i]!=null && dataSt.emissions[i] instanceof IntegerDistribution){
            	 double v = Math.pow(10,Double.parseDouble(geno[k])/-10);
            	 if(v>Constants.seqQualityThresh()){
            	dataSt.emissions[i] = new MixtureDistribution(new double[] {1-v,v},new PseudoDistribution[] {dataSt.emissions[i],stsp.getHWEDist1(null)} );
            	 }
            	 
            	 }
             else if((hk.equals("gl") || hk.equals("pl")) && dataSt.emissions[i]!=null)
            		{
            	boolean pl = hk.equals("pl");
            	if(k>= geno.length || geno[k].equals("./.")|| geno[k].equals(".")) {
            		dataSt.emissions[i] = stsp.getHWEDist1(null);
            			//stsp.getZeroDistProb((int)Constants.coverage(index), ploidy);
            		return false;
            	}
            	else{
            		String st = geno[k];
            		String[] str = st.split(",") ;
            		String[] genos = null;
            		if(ploidy==1) genos = new String[] {"A",null,"B"};
            		else if (ploidy==2) genos =  new String[] {"A,A","A,B","B,B"} ;
            		else throw new Exception("need to modify for diff ploidy");
            		//PseudoDistribution hwe = stsp.getHWEDist1(null);
            		double[] probs = new double[stsp.haplopairListSize()];
            		SimpleExtendedDistribution dist =new SimpleExtendedDistribution(probs, Double.POSITIVE_INFINITY);
            		//	hwe.clone();
            		dist.setDataIndex(this.index);
            		for(int kk=0; kk<genos.length;kk++){
            			if(genos[kk]!=null){
            			int i1 = stsp.getHapl(genos[kk]);
            			double st_ =Double.parseDouble(str[kk]); 
            			double v = pl ?Math.pow(10,st_/-10.0) :  Math.pow(10,st_);
            			dist.addCount(i1, v);
            			}
            		}
            		dist.transfer(0);
            		dist.initialise();
            		dist.setEmStSp(stsp);
            		   double soften = Constants.softenHapMap(this.index);
            		   if(soften>0){
            			  dist.mix(stsp.getHWEDist1(null),Constants.softenHapMap(this.index));
            		   }
            		dataSt.emissions[i] = dist;
            		return false;
            	}
            }
           
            else if(header[k].toLowerCase().indexOf("geno")>=0 || header[k].toLowerCase().indexOf("plus_allele")>=0 ||  header[k].toLowerCase().indexOf("gtype")>=0){
            	if(geno[k].indexOf("null")>=0 ||  geno[k].indexOf('N')>=0 || geno[k].indexOf('-')>=0){
            	
                		dataSt.emissions[i] =stsp.getHWEDist1(null);
                	
            		return false;
            	}
                if(dataSt.emissions[i]==null){
                	int ind;
                	
                	boolean isint = geno[k].matches("[0-9]");
                	if(isint){
                		int cn = Integer.parseInt(geno[k]);
                		
                		ind =  stsp.getGenoForCopyNo(cn)[0];
                	}
                	else
                     ind = 
                     this.abGenos ||
                     alleleA.size()==0 && alleleB.size()==0 ? 
                    	trans(geno[k], stsp, null, null,index) :
                    		trans(geno[k], stsp,alleleA.get(i), alleleB.get(i),index)
                    		;
                    	
                    if(geno[k].equals("null")){
                    	Logger.global.info("h");
                    }
                   double soften = Constants.softenHapMap(this.index);
                  noB = stsp.getBCount(ind);
                   // dataSt.emissions[i] = new IntegerDistribution(ind,stsp);
                    if(deletionIsNa & stsp.getCN(stsp.getHaploPairFromHaplo(ind))==0){
                    	if(Constants.format.length==1){
                    		dataSt.emissions[i] =stsp.getHWEDist1(null);
                    	}
                    	else{
                    		dataSt.emissions[i] = null;
                    	}
                    }
                    else{
                        if(soften>0){
                        	
                        	 dataSt.emissions[i] = stsp.getSoftenedDist(ind, 1-soften);
                        	// dataSt.emissions[i]
                        }
                        else{
                        	
                        	
                        	 dataSt.emissions[i] = stsp.getIntDist(ind);
                        }
                        
                    }
                   if(isint) return true;
                }
            }
            
        }
           
           // data.emissions[i].setDataIndex(this.index);
        }catch(Exception exc){
            exc.printStackTrace();
            dataSt.emissions[i] =stsp.getHWEDist1(null);
        }
       return false;
    }
    
    static double pi = Constants.pi();
    
    public PseudoDistribution process(int cn, double prob_cn, double genotp, int ploidy, String genot){
      	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
      	  if(prob_cn>0.9999 && genotp>0.99 && !genot.equals("NC")){
        		 int ind = stsp.getGenoForCopyNo(cn)[0];
        	         
        	     	//trans(genot, stsp, null, null,index) ;
        	 	return stsp.getIntDist(ind);
        		 
        	   }else{
        		 
        		   int[] genos =  stsp.getGenoForCopyNo(cn);
        		   double[] d = new double[genos.length];
        		   
        		   if(genot.equals("NC")){
        			   Arrays.fill(d, 1.0/(double)d.length);
        		   }
        		   else{
        			 int ind = 
        		         
        		     	trans(genot, stsp, null, null,index) ;
        			   double rem = (1.0-genotp)/(double)(d.length-1);
        			   for(int k=0; k<genos.length; k++){
        				   if(ind==genos[k]){
        					   d[k] = genotp;
        				   }
        				   else{
        					   d[k] = rem;
        				   }
        			   }
        			   if(Constants.sum(d)<0.9999) throw new RuntimeException("!!");
        		   }
        		  double[] hwep = stsp.getHWEDist1(null).probs();
        		   double[] d1 = new double[hwep.length];
        		   double rem1 = (1-prob_cn)/((double)d1.length-genos.length);
        		   Arrays.fill(d1, rem1);
        		   for(int k=0; k<genos.length;k++){
        			   d1[genos[k]] = d[k]*(prob_cn);
        		   }
        		  if(Constants.sum(d1)<0.9999) throw new RuntimeException("!!");
        		return new  SimpleExtendedDistribution(d1, Double.POSITIVE_INFINITY,stsp);
        	   }
    }
    
    public Boolean  processQuanti(String indiv,  String[] header,  String[] geno, int ii, int ploidy){
    	int cn = Integer.parseInt(geno[4]);
    	double prob_cn = 1.0;
    	if(cn!=2){
    		double logbf =Double.parseDouble(geno[5]);
    		double bf =Math.exp(logbf);// Math.pow(2,logbf);
    		double po = bf*(pi/(1-pi));
    		prob_cn = po / (1+po);
    	}
    	String genot = geno[0];
    	double genotp = Double.parseDouble(geno[1]);
    	
    	
   	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
   	 dataSt.emissions[ii] = process(cn,  prob_cn,  genotp, ploidy, genot);
   	  
      return false;
   }
    public Boolean  processForward(String indiv,  String[] header,  String[] genots, int ii, int ploidy){
    	String genot = genots[0]+genots[1];
    	  Character alleleA_ = alleleA.get(ii);
    	  Character alleleB_ = alleleB.get(ii);
    	     
    	genot = genot.replace(alleleA_, 'A').replace(alleleB.get(ii), 'B');
    	
    	
   	 HaplotypeEmissionState dataSt =(HaplotypeEmissionState)  this.dataL.get(indiv);
   	   CompoundEmissionStateSpace stsp = Emiss.getSpaceForNoCopies(ploidy);
     	 if(genot.equals("NC")){
     	
		   }
     	 else{
   		 int ind = 
   	         
   	     	trans(genot, stsp, null, null,index) ;
   	 	 dataSt.emissions[ii] = stsp.getIntDist(ind);
   		 
   	   }
   	 
       
      return false;
   }
    
   
//static int[] alleleCounts;
    public static int trans(String st, CompoundEmissionStateSpace stSpa, 
			Character alleleA, Character alleleB, int index){
    	//alleleCounts = stSpa.tmp();
    	   // int cn_index = st.length()-1;
    	
    	st = st.replace('0', '_');
    	    //st = st.replace('N', '_').replace('C', '_').replaceAll("0", "__").replaceAll("1", "A_").replaceAll("2", "AA")
    	   // .replaceAll("3", "AX").replaceAll("4","XX");
    	   if(alleleA!=null){
    		   st = st.replace(alleleB,'B').replace(alleleA,'A' );
    	   }
    	//  if(Constants.alleles[index]<=2) {
    		
    	 // }
       // Arrays.fill(alleleCounts, 0);
    	  //}
    	  String mod = modify(st);
    	  mod = mod.replaceAll("X", "AA").replaceAll("Y", "AB").replace("Z", "BB");
    	  Integer hapl_ = stSpa.getHapl(mod);
    	  Integer hapl;
    	  if(hapl_ ==null) {
    		  hapl = stSpa.getFromString(mod);
    	  }
    	  else {
    		   hapl = stSpa.getHaploPairFromHaplo(hapl_);
    	  }
    	  if(hapl==null) {
    		  char[] ch = st.toCharArray();
    		  Arrays.sort(ch);
    		 hapl =  stSpa.getFromString(new String(ch));
    		//  throw new RuntimeException("could not convert "+mod);
    	  }
    	  if(hapl==null){
    		  hapl = stSpa.getGenoForCopyNo(Integer.parseInt(st))[0];
    	  }
    	// stSpa.getGenoForHaplopair(hapl);
        //char[] ch = st.toCharArray();
       // for(int i=0; i<ch.length; i++){
       // 	alleleCounts[Emiss.conv.get(ch[i])]++;
      //  }
/*	   int cntA = 0;
    	   
    	   int cntB = 0;
    	   int cntNull=0;
    	 
    	   for(int k=0; k<ch.length; k++){
    		   if(ch[k]=='A') cntA++;
    		   else if(ch[k]=='B') cntB++;
    		   else if(ch[k]=='_')cntNull++;
    		 
    	   }*/
    	 // if(cntA+cntB+cntNull!=ch.length) throw new RuntimeException("!!! problem with genotype "+st);
    //	  System.err.println(st+" "+hapl);
    	 return hapl;
    	}
    
    private static String modify(String st) {
    	if(st.indexOf(',')<0){
    		StringBuffer sb = new StringBuffer();
	for(int i=0; i<st.length(); i++){
		sb.append(i==0 ? st.charAt(i) :  ","+st.charAt(i));
		
	}
	return sb.toString().replaceAll("_", "");
    	}
    	return st.replaceAll("_", "");
}
	public boolean hasABGenos(ZipFile zf){
    	if(Constants.impute() || true){
    	try{
    		
    		
    	for(Enumeration en = zf.entries(); en.hasMoreElements();){
    		ZipEntry ent = (ZipEntry) en.nextElement();
    		if(ent.getName().equals("SNPS") || ent.getName().equals("Name") || ent.getName().equals("Samples")) continue;
    		List<String> l = Compressor.getIndiv(zf, ent.getName(), geno_id);
    		for(int k=0; k<l.size(); k++){
    			if(l.get(k).indexOf("B")>=0) {
    				return true;
    			}
    			else if(l.get(k).indexOf("C")>=0) {
    				return false;
    			}
    			else if(l.get(k).indexOf("T")>=0) {
    				return false;
    			}
    			else if(l.get(k).indexOf("G")>=0) {
    				return false;
    			}
    		}
    	}
    	}catch(Exception exc){
    		exc.printStackTrace();
    		return true;
    	}
    	return true;
    	}
    	return false;
    	
    	
    }
    
    public void removeKeyIfStartsWith(String st) {
        for(Iterator<String> it = data.keySet().iterator(); it.hasNext();){
            String key = it.next();
            if(key.startsWith(st) && !st.equals("NA12239")){
             recSites.remove(key);
               viterbi.remove(key);
             it.remove();
            }
        }
        
    }
    @Override
    public void createDataStructure(List<String> indiv, List<Integer>ploidy,List<Integer> samps){
        for(int i1=0; i1<samps.size(); i1++){
        	int i = samps.get(i1);
            PhasedDataState emst = SimpleScorableObject.make(indiv.get(i),loc.size(), null, this.index);
            emst.setNoCop(ploidy.get(i));
            dataL.put(indiv.get(i),emst );
          //  data.put(indiv.get(i), SimpleScorableObject.make(indiv.get(i),loc.size(), stSp[no_copies-1], this.index));
        }
    }
    
    public void  restrictDataSites(int i){
        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
          it.next().restrictSites(i);
        }
    }
    
   
    
  
   
    public void add(List<PIGData> toAdd) {
        for(Iterator<PIGData> it = toAdd.iterator(); it.hasNext();){
            PIGData nxt = it.next();
            this.data.put(nxt.getName(),nxt);
        }
          
      }
    
  
   
    
    
    public void remove(List<String> toRemove) {
       
        for(int i=toRemove.size()-1; i>=0; i--){
            this.data.remove(toRemove.get(i));
        }
    }
    /** splits the data into haplotype data */
   public void split(){
       Set<PIGData >set =new HashSet<PIGData>();
       for(Iterator<PIGData> it = this.iterator(); it.hasNext();){
           PhasedDataState da = (PhasedDataState) it.next();
           PIGData[] dat = da.split();
           for(int j=0; j<dat.length; j++){
               dat[j].setName(da.getName()+"_"+j);
               set.add(dat[j]);
           }
       }
       this.data.clear();
       for(Iterator <PIGData> it = set.iterator(); it.hasNext();){
           PIGData da = it.next();
           this.data.put(da.getName(), da);
       }
   }
   
    public void drop(List<Integer> toDrop) {
        Collections.sort(toDrop);
        this.dataL.clear();
       for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
           PIGData data = it.next();
           data.removeAll(toDrop);
       }
       for(int i=toDrop.size()-1; i>=0 ; i--){
           this.loc.remove(toDrop.get(i).intValue());
           if(snpid!=null && snpid.size()>0) this.snpid.remove(toDrop.get(i).intValue());
       }
    //  this.maf = null;
       this.length=loc.size();
        
    }
    public int restrict(int st1, int end1, int lkb, int rkb, boolean kb) {
        int il = st1;
        int ir = end1;
        if(kb){
        for(; il>=0; il--){
            if(this.loc.get(st1) - this.loc.get(il) > lkb*1000){
                il++;
                break;
            }
        }
    
        for(; ir<loc.size(); ir++){
            if(this.loc.get(ir) - this.loc.get(end1) > rkb*1000){
                //ir--;
                break;
            }
        }
        }
        else{
            il= st1 - lkb;
            ir = end1+rkb;
        }
        il = Math.max(0, il);
         ir = Math.min(loc.size()-1, ir);
         
        this.loc = this.loc.subList(il, ir);
        this.snpid = this.snpid.subList(il, ir);
        for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
            it.next().restrictSites(il, ir);
        }
        this.length = loc.size();
        return il;
    }
    public boolean cnvBiggerThan(int i, int lenThresh, EmissionStateSpace stSp) {
       for(Iterator<PIGData> it = this.data.values().iterator(); it.hasNext();){
           PIGData dat = it.next();
           Comparable comp = dat.getElement(i);
           int cnv = stSp.getCN(stSp.get(comp));
           if(cnv==1) continue;
           int cnt =1;
           
           for(int i1=i+1; i1-i<lenThresh && i1<dat.length(); i1++){
               int cnv1 = stSp.getCN(stSp.get(dat.getElement(i1)));
               if(cnv1==cnv) cnt++;
               else break;
           }
           for(int i1=i-1; i-i1<lenThresh && i1>=0; i1--){
               int cnv1 = stSp.getCN(stSp.get(dat.getElement(i1)));
               if(cnv1==cnv) cnt++;
               else break;
           }
           if(cnt>=lenThresh) return true;
       }
       return false;
    }
	
    
  
  
   
}
