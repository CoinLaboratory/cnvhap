package lc1.util;

import java.io.BufferedReader;
import java.io.File;
import java.util.List;

public interface ZipFileAccess {

	public abstract List<String> getIndiv(String entryName, Integer column)
			throws Exception;

	public abstract List<String> getIndiv(String entryName, Integer column,
			String spl) throws Exception;

	public abstract boolean getIndiv(String entryName, Integer column,
			String[] indiv) throws Exception;

	public abstract BufferedReader getBufferedReader(String string)
			throws Exception;

	public abstract void getAvgDepth(String pref, int avgDepthCol, List<Integer> dToInc,
			File samplesFile,List<Integer> ploidy,List header_sample,  List avgDepth);
}