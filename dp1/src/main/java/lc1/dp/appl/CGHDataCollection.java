package lc1.dp.appl;

import java.io.File;
import java.util.Collection;
import java.util.List;

import lc1.dp.data.collection.IlluminaRDataCollection;

public class CGHDataCollection extends IlluminaRDataCollection {

	public CGHDataCollection(File f, short index, int no_copies, int[][] mid, File bf, Collection<String> snpidrest) throws Exception {
		super(f, index, no_copies, mid, bf,snpidrest);
		// TODO Auto-generated constructor stub
	}
	protected double calcBaf(List<String> l) {
	return 0;
	}

	

}
