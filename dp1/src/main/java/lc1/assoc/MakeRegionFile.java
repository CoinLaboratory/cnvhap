package lc1.assoc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;

public class MakeRegionFile {

	
	public static void main(String[] args){
		try{
	    File dir1 =  new File(System.getProperty("user.dir"));
		Map<String, List<Location>>m = getChrom(dir1);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("regions.txt"))));
		for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String key = it.next();
			List<Location> loc = m.get(key);
			for(int i=0; i<loc.size(); i++){
				Location lo = loc.get(i);
				pw.println("chr"+key+":"+lo.getMin()+":"+lo.getMax());
			}
		}
		pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	static Comparator compa = new Comparator< Location>(){

		public int compare(Location o1, Location o2) {
			return (new Integer(o1.getMin())).compareTo(o2.getMin());
		}
		
	};

	private static Map<String, List<Location>> getChrom(File dir1) {
		String[] str = dir1.list();
		Map<String, List<Location>> res = new HashMap<String, List<Location>>();
		outer: for(int i=0; i<str.length; i++){
			String[] st = str[i].split("_");
			List<Location> l = res.get(st[0]);
			if(l==null) res.put(st[0], l = new ArrayList<Location>());
			int arg0 = Integer.parseInt(st[1]);
			Location loc = LocationTools.makeLocation(arg0, arg0);
			for(int k=0; k<l.size(); k++){
				if(mindist(l.get(k), loc) < 100*1000){
					l.set(k, loc.union(l.get(k)));
					continue outer;
				}
			}
			l.add(loc);
		}
	
		
		return res;
	}

	private static int mindist(Location location, Location loc) {
		int d1 = Math.abs(location.getMax() - loc.getMin());
		int d2 = Math.abs(loc.getMax() - location.getMin());
		return Math.min(d1, d2);
	}
}
