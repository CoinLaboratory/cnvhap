package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.PIGData;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpaceTranslation;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.ensj.PedigreeDataCollection;
import lc1.util.Constants;

public class HLADataCollection extends DataCollection {
    public HLADataCollection(){
        
    }
    protected HLADataCollection(DataCollection dat){
      super(dat);
      
      
       
    }
    public boolean hasIntensity(int i) {
    	return false;
    }
    public  static void readIllumina (File f, HLADataCollection sdt) throws Exception{
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
   //     File dir1 = new File(Constants.getDirFile());
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
          
            data[i1] = SimpleScorableObject.make(indv[i1], sdt.loc.size(), stSp,(short)-1);
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
    //    sdt.calculateMaf(false);
        br.close();
    }
 
    
    public static String[]  readAffy(File dir1, String name,boolean quotes, HLADataCollection res, int noCopies) throws Exception{
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
   
    public   void transform ( double[] fraction, List<String>keys,  String join){
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
    
  /*  public void dropSingle(){
    	List<String> toD = new ArrayList<String>();
    	for(Iterator<String> it = this.data.keySet().iterator(); it.hasNext();){
    		String st = it.next();
    		PhasedDataState state = (PhasedDataState) this.data.get(st);
    		if(state.getEmissionStateSpace().noCopies()==1) toD.add(st);
    	}
    	this.dropIndiv(toD.toArray(new String[0]));
    }*/
    public   void transform (Iterator<String[]> pairs, String join){
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
            PIGData pid = SimpleScorableObject.make(unit, unit[0].noCop()==1, join, null);
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
    }
    
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
    public HLADataCollection(int length) {
        super(length);
        // TODO Auto-generated constructor stub
    }
  
   
   
   
    /*  private List<PIGData> data_bu;
    private List<PIGData> trio_bu;*/
    
    public HLADataCollection(List<PIGData> name) {
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
    
 
   public HLADataCollection(File file,short index, int no_copies, int[][] mid,  File bf
		   ,Collection<String> snp_ids_to_rest) throws Exception{
       super(file, index, no_copies, mid, bf, snp_ids_to_rest);
    }
   @Override 
   public int restricToAlias(Collection<String> alias){
	   super.restricToAlias(alias);
	   List<String> keys = new ArrayList<String>(this.getKeys());
       for(Iterator<String> it = keys.iterator(); it.hasNext();){
           String key = it.next();
           if(!alias.contains(key)){
           
            data2.remove(key);
           }
       }
       int sze = dataL.size();
      return sze;
   }
    public  HLADataCollection clone(){
        return new HLADataCollection(this);
    }
    SortedMap<String, List<Comparable>> data2;
    @Override
    public boolean  process(HaplotypeEmissionState state, HaplotypeEmissionState stateL, String header, String geno, int i, double[] missing){
        try{
         
         //  List<Comparable> l = new ArrayList<>
        	List<Comparable> st = data2.get(indiv);
        	st.add(ComparableArray.make(geno.charAt(0), geno.charAt(1)));
            
           
        }catch(Exception exc){
            exc.printStackTrace();
        }
        return true;
       
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
    public void createDataStructure(List<String> indiv, List<Integer> ploidy, List<Integer> sampIdToIncl){
    	data2 = new TreeMap<String, List<Comparable>>();
        for(int i=0; i<sampIdToIncl.size(); i++){
         	this.data2.put(indiv.get(sampIdToIncl.get(i)), new ArrayList<Comparable>());
        }
    }
   static public Map<Character, String>[] conversion = 
	   new Map[] {
	   new HashMap<Character, String>(),
	   new HashMap<Character, String>()
   };
   @Override
    public void process(String[] str, int k,int no,  int loc_index, int[] maf_index,int chr_index, int strand_index,int snp_index,
    		List l, List chr, List majorAllele, List alleleB, List forward, int bin){
	   super.process(str, k, no, loc_index, maf_index, chr_index, strand_index, snp_index, l, chr, majorAllele, alleleB, forward, bin);
    	for(int i=0; i<str.length; i++){
    		if(str[i].indexOf('=')>0){
    			String[] str1 = str[i].split("=");
    			conversion[k].put(str1[1].charAt(0), str1[0]);
    		}
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
     // this.maf = null;
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
	
	public List<Comparable> getGenotypes(int pos){
		   List<Comparable> res = new ArrayList<Comparable>();
		 //  Set<String> keys = data.keySet();
		   
		    for(Iterator<String> it = this.data2.keySet().iterator(); it.hasNext();){
				String st = it.next();
			     List<Comparable> str = data2.get(st);
				//for(int i=0; i<str.size(); i++){
			     ComparableArray compa = (ComparableArray) str.get(pos);
			     
					res.add(ComparableArray.make(conversion[pos].get(compa.get(0)),
							conversion[pos].get(compa.get(1))));
				//}
			}
		    return res;
		}  
  
  
   
}
