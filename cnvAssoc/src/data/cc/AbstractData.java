package data.cc;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.linear.RealMatrix;

public abstract class AbstractData {
	List<String> indiv;
	RealMatrix data;
	List<String> probes;
	String name;
	public String type_name;
	public List<Integer> loc = new ArrayList<Integer>();
	
	
	
}
