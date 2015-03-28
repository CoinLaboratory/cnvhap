package lc1.dp.emissionspace;

import java.util.ArrayList;
import java.util.List;

public class PosEmissionStateSpace<T extends EmissionStateSpace> {

	List<T> l = new ArrayList<T>();
	
	public T get(int i){
		return l.get(i);
	}
	
}
