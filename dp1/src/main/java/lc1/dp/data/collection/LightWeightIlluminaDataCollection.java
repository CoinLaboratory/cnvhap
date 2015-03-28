package lc1.dp.data.collection;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.PIGData;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.IlluminaNoBg;

import org.apache.commons.compress.archivers.zip.ZipFile;
/** like likelihood datacollection, but only reads a single snp ata time for mem efficiency */
public class LightWeightIlluminaDataCollection extends IlluminaRDataCollection {
	double[] missing = new double[2];
	public  LightWeightIlluminaDataCollection (File f, short index, int no_copies, int[][] mid,
			File bf, Collection<String> snpidrest) throws Exception{
		super(f, index, no_copies, mid, bf, snpidrest);
		 zf = new ZipFile(f);
		 this.samps = this.getSamps();
		 try{
			
		super.process(snpid.get(0), 0, zf, null, samps,missing);
	//	super.calculateMaf(true, false);
		 }catch(Exception exc){
	    		System.err.println("could not process "+snpid.get(0));
	    		exc.printStackTrace();
	    	}
	}
	public String[][] getHaplotypes(int st, int end){
		
	    String[][] str = new String[data.values().size()][];
	    for(int i=0; i<str.length; i++){
	    	str[i] = new String[2];
	    }
	   
	    {
	    	 this.updateIndex(st);
			    Iterator<PIGData> it = this.data.values().iterator();
			   
			    outer: for(int i=0; it.hasNext(); i++){
			        PIGData nxt = it.next();
			        str[i][0] = nxt.getStringRep(0);//, nxt.getStringRep(end)};
			        
			    }
	    }
	    {
		    this.updateIndex(st);
		    Iterator<PIGData> it = this.data.values().iterator();
		   
		    outer: for(int i=0; it.hasNext(); i++){
		        PIGData nxt = it.next();
		        str[i][1] = nxt.getStringRep(0);//, nxt.getStringRep(end)};
		        
		    }
	    }
	    return str;
	}
	@Override
	public void append(int pos, StringBuffer[] sb){
		this.updateIndex(pos);
		for(int i=0; i<indiv.size(); i++){
    		ComparableArray comp  = (ComparableArray)(data.get(indiv.get(i)).getElement(0));
    		int no_copies = comp.size();
    		for(int k=0; k<no_copies; k++){
    			sb[i*no_copies+k].append(comp.get(k).toString());
    		}
    	}
    }
	public void getCompa(int pos, Comparable[] genotypes){
		this.updateIndex(pos);
		for(int i=0; i<indiv.size(); i++){
			genotypes[i] = (data.get(indiv.get(i)).getElement(0));
			
		}
	}
	 
	@Override
	protected Boolean process(String snpd, int i,ZipFile zf,
			List<Integer> ploidy, List<Integer> toinc, double[] missing) {
		return null;
	}

	@Override
	public HaplotypeEmissionState createEmissionState(String key, int no_copies){
        //if(stSp[1].size()==stSp1[1].size()) 
          //  return   SimpleScorableObject.make(key, 1, stSp[no_copies-1]);
		  return   new IlluminaNoBg(key, null,  this.snpid, index);
       
    }
	int currPosIndex =-1;
	public synchronized String getPhenInfo(String string , int pos_index, int phenIndex, int type){
		updateIndex(pos_index);
		if(pos_index!=currPosIndex){
		this.currentPosScIndex = -1;
		currPosIndex = pos_index;
		}
		return super.getPhenInfo(string, 0, phenIndex, type);
	}
	public String getInfo(String string, int pos_index) {
			updateIndex(pos_index);
		if(string.indexOf("maf")>=0) {
			return super.getInfo(string, 0);
		}	
		else{
			return super.getInfo(string, pos_index);
		}
	}
	final List<Integer> samps;
	public void updateIndex(int i) {
		 if(i!=currentIndex){
			 currentIndex=i;
			 try{
				 super.process(this.snpid.get(i), 0, zf, null, this.samps,missing);
				// super.calculateMaf(true,false);
			 }catch(Exception exc){
		    		System.err.println("could not process "+snpid.get(i));
		    		exc.printStackTrace();
		    	}
			 }
		
	}
	int currentIndex = 0;
	 ZipFile zf;
	 public String getCompressedString(String key, int i, boolean b,boolean b2){
		 updateIndex(i);
	    	return super.getCompressedString(key, 0, b,b2);
	 }
	/* public EmissionState makeMafState(EmissionStateSpace emStSp1){
	        if(emStSp1 instanceof CompoundEmissionStateSpace){
	            CompoundEmissionStateSpace emStSp = (CompoundEmissionStateSpace) emStSp1;
	        EmissionStateSpace[] ems = emStSp.getMembers();
	        EmissionState[] st = new EmissionState[ems.length];
	        for(int i=0; i<st.length; i++){
	            st[i] = new HaplotypeEmissionState("maf_"+i, 1, ems[i].size(), ems[i], null, null);
	        }
	        
	        this.maf =  new AlleleCopyPairEmissionState(Arrays.asList(st), emStSp, false, null);// emStSp.size(), emStSp);
	        }
	        else{
	            return  new HaplotypeEmissionState("maf_", 1, emStSp1.size(), emStSp1,  null, null);
	        }
	        maf.initialiseCounts();
	        return maf;
	    }*/
	   
	   
	
}
