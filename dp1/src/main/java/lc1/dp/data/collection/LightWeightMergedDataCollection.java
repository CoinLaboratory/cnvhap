package lc1.dp.data.collection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import lc1.dp.emissionspace.CompoundEmissionStateSpace;

public class LightWeightMergedDataCollection extends MergedDataCollection {

	final int midpoint,last;
	  public LightWeightMergedDataCollection(LightWeightDataCollection[] ldl, String name){
		  super(ldl, name, null, false);
		 this.which = new int[ldl.length];
		 Arrays.fill(which,-1);
		 last = ldl[0].last;
		 midpoint = ldl[0].midpoint;
		  this.dataL = ldl[0].dataL;
	  }
	  public LightWeightMergedDataCollection(LightWeightDataCollection[] ldl,
			String name, Info[][] map) {
		  super(ldl, name, null, false,map);
		  last = ldl[0].last;
		  this.which = new int[ldl.length];
			 midpoint = ldl[0].midpoint;
		/*  if(!ldl[0].indiv().equals(ldl[1].indiv())){
			  List<String> ind1 = ldl[0].indiv;
			  List<String> ind2 = ldl[1].indiv;
			  throw new RuntimeException("prob with "+ldl[0].chrom+" "+ldl[0].name+" "+ldl[1].name+" "+ind1.size()+" "+ind2.size());
		  }*/
		  this.dataL = ldl[0].dataL;
		  Arrays.fill(which,-1);
		// TODO Auto-generated constructor stub
	}
	public int restricToAlias(Collection<String> alias) {
		boolean alleq = true;
		int sze =0;  
		for(int i=0; i<ldl.length; i++){
			  sze+=ldl[i].restricToAlias(alias);
			  if(i>1 && !ldl[i].indiv().equals(ldl[0].indiv)){
				  alleq = false;
			  }
		  }
		  this.sameIndiv = alleq;
		return sze;
	  }
	  
	 boolean sameIndiv = false;
	public List<String> indiv(){
		if(sameIndiv)
		  return ldl[0].indiv();
		else{
			List<String> indiv = new ArrayList<String>(ldl[0].indiv);
			for(int k=1; k<ldl.length; k++){
				indiv.retainAll(ldl[k].indiv());
			}
			return indiv;
		}
	  }
	
	
final int[] which;
	
	@Override
	  public boolean updateIndex(String rs) {
		for(int k=1; k<which.length; k++){
			which[k-1] = which[k];
		}
		which[last] = -1;
			for(int i=0; i<this.ldl.length; i++){
				if(ldl[i].updateIndex(rs)){
					which[last] = i;
				
				}
				if(ldl[i].canProcess(rs)){
					ldl[i].updateIndex(rs);
					this.dataL = ldl[i].dataL;
					return true;
				}
			}
			this.dataL  =which[midpoint]>=0 ? ldl[which[midpoint]].dataL : null;
			
			return which[last]>=0;
			
		}
	 public boolean canProcess(String rs){
		 for(int i=0; i<this.ldl.length; i++){
				if(ldl[i].canProcess(rs)){
					return true;
				}
			}
			
			return false;
		}
	  @Override
	  protected void prepareStates(List< String > keysL, Set<String> problems, CompoundEmissionStateSpace[] emstsp){
		  
	  }
	  @Override
	    public  void createDataStructure(List<String> indiv,List<Integer> ploidy, List<Integer>sampid){
		  
	  }
	  @Override
	  public  void switchAlleles(DataCollection[] ldl2){
	       
	    }
	  
	  
}
