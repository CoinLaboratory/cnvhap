package calc;

import java.util.List;

import lc1.dp.data.representation.Emiss;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class LDCalc {
	DoubleMatrix2D matrix;
	List<Emiss>[] emm;
	LDCalc(List<Emiss>[] emm){
		this.emm = emm;
	}
	
	public double calculate(int type){
		matrix = new DenseDoubleMatrix2D(emm[0].size(), emm.length);
		for(int i=0; i<emm.length; i++){
			for(int j=0; j<emm[i].size(); j++){
			//	matrix.setQuick(j, i, type ==0 ? emm[i].get(j).noCopies() :
			//		(type==1 ? emm[i].get(j).));
			}
		}
		return 0;
	}
	
}
