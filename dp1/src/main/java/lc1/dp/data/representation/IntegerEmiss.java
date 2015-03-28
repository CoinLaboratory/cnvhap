package lc1.dp.data.representation;

import lc1.util.Constants;

public  class IntegerEmiss implements AbstractEmiss{
public IntegerEmiss(int index) {
	this.v = index;
	}

public int noCopies(){
	return 1;
}

public int noB(){
	return 1;
}

public int noB(int i){
	return 1;
}
public Integer v;
public String toStringShort(){
	return Integer.toString(v, Constants.radix());
}
public String toString(){

	return  Integer.toString(v, Constants.radix());
}

@Override
public int compareTo(Object o) {
	if(o instanceof ComparableArray) return -1;
	else{
		return v.compareTo(((IntegerEmiss)o).v);
	}
}
}
