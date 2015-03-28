package lc1.stats;

import lc1.dp.emissionspace.EmissionStateSpace;


public class ACSNDistribution extends CompoundDistribution{
 public ACSNDistribution(IlluminaRDistribution dist1, IlluminaRDistribution dist2, EmissionStateSpace emstsp) {
		super(dist1,dist2,emstsp);
		this.dist1 = dist1;
		this.dist2 = dist2;
	}
public ACSNDistribution(short index, EmissionStateSpace emstsp) {
	this(new IlluminaRDistribution(index),new IlluminaRDistribution(index), emstsp);
	this.data_index = index;
}


public void setDataIndex(short data_index2) {
    this.data_index = data_index2;
     for(int i=0; i<l.size(); i++){
    	 l.get(i).setDataIndex(data_index2);
     }
 }
@Override
public double[] probs() {
  return null;
}
public IlluminaRDistribution dist1,dist2;
public double scoreR(int bg, int k, int i){
	return 1.0;
}
public double scoreB(int j, int i){
	return 1.0;
}

public double scoreBR(int j, int i){
	throw new RuntimeException("!!");
	//return 1.0;
}
//this.emissions[i].scoreR(emStSp.getCN(j),i)
//*this.emissions[i].scoreB(j,i)
 



}
