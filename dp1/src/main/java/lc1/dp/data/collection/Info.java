/**
 * 
 */
package lc1.dp.data.collection;

import java.io.Serializable;
import java.util.List;
import java.util.logging.Logger;

import lc1.util.Constants;

public class Info implements Serializable{
	 public Boolean strand, probeOnly;
    public Info(String id, Character b_all, Character a_all, int data_index, int relative_position,
    		Integer loc2, double afreq, Boolean boolean1, Boolean probeOnly) {
        this.id = id;
       this.alleleB = b_all;
       this.alleleA = a_all;
       this.data_index = data_index;
       this.relative_position = relative_position;
       this.loc = loc2;
       this.Afreq = afreq;
       this.strand = boolean1;
       this.probeOnly = probeOnly;
     //  if(strand==null && iscompl(alleleA, alleleB)) throw new RuntimeException("need to specify strand for study "+ldl[data_index].name);
    }
    public Boolean isProbeOnly() {
    	return probeOnly;
	}
	public String toString(){
        return id+"_"+loc+"_"+data_index+"_"+relative_position+"_"+alleleA+"_"+alleleB;
    }
    String id;
    Character alleleA;
    Character alleleB;
    Integer loc;
    double Afreq;
    public int data_index;
	public int relative_position;
    public void compare(Info res, List switchedAlleles) {
        if(Math.abs(res.loc.intValue()-loc.intValue())>100){
            
            throw new RuntimeException("same id at different location " +id+" "+res.loc+" "+loc);
        }
        else {
           // System.err.println("same id loc "+id+" "+res.loc+" "+loc);
        }
        
     
        Boolean swtch = swtch(res.alleleA,alleleA, true, null);
       swtch =swtch(res.alleleA,alleleB, false, swtch);
       swtch= swtch(res.alleleB,alleleA, false, swtch);
        swtch =swtch(res.alleleB,alleleB, true, swtch);
      
            if(swtch!=null && swtch){
            	System.err.println("switch "+res+" "+this);
                switchedAlleles.add(res);
               // compare(this.Afreq, 1.0- res.Afreq, res.data_index);
            }/*else {
            	System.err.println("no switch "+res+" "+this);
            }*/
           
       // }
    }
    private Boolean swtch(Character alleleA2, Character alleleA3, boolean same, Boolean resP) {
		// TODO Auto-generated method stub
		if(alleleA2!=null && alleleA3!=null){
			boolean res = alleleA2.charValue()==alleleA3.charValue();
			boolean result = same ? !res : res;
			if(resP!=null && resP.booleanValue()!=result){
				throw new RuntimeException("!!");
			}
			return result;
		}
		return resP;
	}
	private void compare(double afreq2, double afreq3, int data_ind1) {
        Logger.global.info("freq difference "+Math.abs(afreq2 - afreq3)+" at "+this.id+" between "+data_index+" "+data_ind1);
        
    }
	
	public void swapAlleles(){
		Character tmp = alleleA;
		this.alleleA = alleleB;
		this.alleleB = tmp;
	}
	public void flipStrand() {
		this.alleleA = alleleA==null ? null : DataCollection.compl(alleleA);
		this.alleleB = alleleB==null ? null : DataCollection.compl(alleleB);
	//	if(strand!=null && strand) throw new RuntimeException("!!");
		strand = true;
	}
 }