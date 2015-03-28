package lc1.dp.data.collection;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.model.MarkovModel;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

public abstract class AssociationCalculator {

	final String name;
	
	public AssociationCalculator(String name){
		this.name = name;
	}
	
	public String name(){
		return name;
	}
	
	abstract public void scoreChi1(double[] prob, int i, boolean b, String name);

	abstract public void scoreChi1(StateDistribution emissionC, int i, boolean b, String name);

	 public void initialise(){}
	 MarkovModel hmm;
	  public final  void setModel(MarkovModel hmm){
		  this.hmm =hmm;
		  int modelLength = hmm.modelLength();
		  EmissionStateSpace emstsp = Emiss.getSpaceForNoCopies(this.ploidy);
		  int modelLength1 =  ((CompoundMarkovModel)hmm).getMemberModels()[0].modelLength();
		
         int new_type_len =modelLength1+this.state_in-1;
         if(Constants.saveStates() && new_type_len!=type_len){
        	 reset();
       		this.type_len=new_type_len;
         }
	  }
	protected abstract void reset() ;
	 public  static int armitage = 0;
	 public    static int beta =1;
	 public  static int chisq = 2;
	 public static int linear = 3;
	 
	 public int width =0;
	 
	public DataCollection dc1;
	   static final double log2 = Math.log(2);
	    static double maxOdds = 20;
public int type_len;
	 public int ploidy;
	  Map<String, Integer> keysToIndex = new HashMap<String, Integer>();
int[][] cn_alias;
public  int len1;
public  int state_in=1;
public  List<String>type_assoc;

public List<String> types_all = new ArrayList<String>();
public boolean[] rank; //can we use type to rank


public AssociationCalculator(AssociationCalculator[] li, DataCollection dc) {
	this.name = li[0].name;
	this.width = li[0].width;
	this.ploidy = li[0].ploidy;
	this.cn_alias = li[0].cn_alias;
	this.type_assoc = li[0].type_assoc;
	this.types_all = li[0].types_all;
	this.len1 = li[0].len1;
	this.dc1 = dc;
	int prev_count=0;
	for(int i=0; i<li.length; i++){
		Iterator<Entry<String, Integer>> it = li[i].keysToIndex.entrySet().iterator();
		while(it.hasNext()){
			Entry<String, Integer> nxt = it.next();
			this.keysToIndex.put(nxt.getKey(), nxt.getValue()+prev_count);
		}
		prev_count = this.keysToIndex.size();
	}
}

public static void printResults(File dir, DataCollection dc){
	List<AssociationCalculator>[][] ac_ = (dc).getArmitageCalculator();
	 if(ac_==null) return; 
	 for(int ii=0; ii<ac_.length; ii++){
		 if(ac_[ii]!=null){
			 for(int jj=0; jj<ac_[ii].length; jj++){
			 for(int j=0; j<ac_[ii][jj].size(); j++){
				 	ac_[ii][jj].get(j).printResults1(dir);
			 }
			 }
		 }
	 }
}
public abstract Double getSignificance(int i, int overall_ind, int j2);

public abstract Double[][] oddsRatio(int i, int overall_ind, int j2) ;

public abstract void printResults1(File dir);
public String getPheno(int kj) {
	return this.types_all.get(kj);
}
	

}