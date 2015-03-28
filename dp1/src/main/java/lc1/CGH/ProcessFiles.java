package lc1.CGH;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import lc1.dp.data.collection.DataCollection;


public class ProcessFiles {
	
	public static void main(String[] args){
		processCGH();
	}
	static String[] getChr()
	{
		
		List<String> c = new ArrayList<String>();
		for(int i=22; i>=1; i--){
			c.add(""+i);
		}
		c.add("X");
		return c.toArray(new String[0]);
	}
	

	public static void processSNP(){
		try{
		   
			String abs = "MultiProbeByIndividual.txt";
			String probes = "IC_Custom244K_238459_ProbesTable.txt";
			String outFile = "MultiP";
			
			 File user = new File(System.getProperty("user.dir"));
			String[] chr = getChr();
		File outDir = new File(user, outFile);
         if(!outDir.exists()) outDir.mkdir();
		 for(int i=0; i<chr.length; i++){
			 AbstractAberatiionReader agilent =    new MultiProbeAberationReader(Long.MAX_VALUE,"");
			 AbstractAberatiionReader p244k = new AgilentProbeReader(Long.MAX_VALUE);
			 p244k.initialise(user, chr[i], null, 0, 0, probes, null);
			 agilent.initialise(user, chr[i], null, 0,0,abs, null);
			 agilent.sort();
			 Locreader sp = extractSingleProbes(agilent, p244k,false, 0,new int[] {0,Integer.MAX_VALUE});
			 print(sp, outDir, "probes_"+chr[i]+".txt");
			 print(agilent, outDir, "regions_"+chr[i]+".txt");
		 }
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	
	public static void processCGH(){
		try{
			String abs = "MultiProbeByIndividual.txt";
			String probes = "IC_Custom244K_238459_ProbesTable.txt";
			String outFile = "MultiP";
			 Set<String> individuals = new HashSet<String>();
             BufferedReader indv = new BufferedReader(new FileReader(new File("indiv_shared.txt")));
               String st1 = "";
               while((st1 = indv.readLine())!=null){
                   individuals.add(st1.trim());
               }
               indv.close();
			 File user = new File(System.getProperty("user.dir"));
			String[] chr = getChr();
		File outDir = new File(user, outFile);
         if(!outDir.exists()) outDir.mkdir();
		 for(int i=0; i<chr.length; i++){
			 AbstractAberatiionReader agilent =    new MultiProbeAberationReader(Long.MAX_VALUE,"");
			 AbstractAberatiionReader p244k = new AgilentProbeReader(Long.MAX_VALUE);
			
		
			 p244k.initialise(user, chr[i], null, 0, 0, probes, null);
			 p244k.sort();
			 agilent.initialise(user, chr[i], null, 0,0,abs, individuals);
			 agilent.sort();
			   agilent.restrictEnds(p244k);
			   agilent.sort();
			 Locreader sp = extractSingleProbes(agilent, p244k,false, 0,new int[] {0,Integer.MAX_VALUE});
			 File data = new File("../data/"+chr[i]+" ");
			 List<Integer> locs = DataCollection.readPosInfo(new File("../data/"+chr[i]+"_data.txt"), 4, true);
			 Locreader locsL = new Locreader(10000000,"");
			 for(int ik=0; ik<locs.size();ik++){
			     locsL.add(new Location(chr[i], locs.get(ik), locs.get(ik)));
			 }
			 agilent.getNoProbes(locsL);
			// print(sp, outDir, "probes_"+chr[i]+".txt");
			 print(agilent, outDir, "regions_"+chr[i]+".txt");
			
		 }
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	
	  private static void print(Locreader agilent, File outDir, String string)  throws Exception{
	   
	     List<Location> l = new ArrayList<Location>();
		 for(java.util.Iterator<Location> it = agilent.iterator(); it.hasNext();){
			 Location loc = it.next();
			 l.add(loc);
			// 
		 }
		 Collections.sort(l);
		  PrintWriter regions = new PrintWriter(new BufferedWriter(new FileWriter(new File(outDir, string))));
		 for(int i=0; i<l.size(); i++){
		     Location loc = l.get(i);
		     loc.print(regions);
		 }
	     regions.close();
		// TODO Auto-generated method stub
		
	}
	public static Locreader extractSingleProbes(Locreader ag1, Locreader agPr, boolean noExclusions, int exact, int[] maxMatchPerFeature){
	        Locreader locr = new Locreader(Integer.MAX_VALUE,"");
	       
	        if(agPr.number()==0){
	            return locr;
	        }
	        Collection<Location> probesinMultiRegion=  new ArrayList<Location>();
	        Collection<Location> multiRegionWithProbes=  new ArrayList<Location>();
	        Collection<Location> exclude = null;//AbFinder.overlap(ag1.iterator(), agPr, null,null, 0, 0, multiRegionWithProbes, probesinMultiRegion, exact,null, 
	                //maxMatchPerFeature);
	       if(noExclusions && exclude.size()>0){
	           throw new RuntimeException("!! "+exclude);
	       }
	       locr.addAll(probesinMultiRegion.iterator()); 
	       return locr;
	    }
}
