package lc1.dp.data.representation;

import java.io.Serializable;

public interface AbstractEmiss extends Comparable,  Serializable{
public int noCopies();

public int noB();

public String toStringShort();

public int noB(int i);
}
