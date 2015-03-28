package lc1.possel;

import java.io.File;
import java.io.FilenameFilter;

public class FilenameFilter1 implements FilenameFilter {
	FilenameFilter1(String pre){
		this.pre = pre;
	}
	String pre;
	public boolean accept(File arg0, String arg1) {
	return arg1.startsWith(pre);
	}

}
