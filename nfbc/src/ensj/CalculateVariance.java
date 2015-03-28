package ensj;

/*
Copyright (C) 2001 EBI, GRL

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipFile;

import lc1.util.Compressor;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.KaryotypeBand;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Sequence;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.SequenceAdaptor;
import org.ensembl.registry.Registry;
import org.ensembl.variation.driver.VariationDriver;

import conversion.OptionBuild;


/**
* Example code demonstrating how to retrieve data from ensembl databases
* using ensj.
* 
* The examples below retrieve data from the latest human core and variation databases on 
* ensembldb.ensembl.org
* 
* @author Craig Melsopp
* @see <a href="Example.java.html">Example.java</a> source <a href="Example.java">(txt)</a>
* 
*/

public class CalculateVariance {

 public CoreDriver coreDriver;

 private VariationDriver variationDriver;

;
public static void main(String[] args){
	main(new File("."), args[0].replace('_',' '));
}
 
public static void main(File dir, String lrrst) {
	List<Double> l1_auto = new ArrayList<Double>();
	List<Double> l1_sex = new ArrayList<Double>();
	 try{
	 File in = new File(dir, OptionBuild.build+".gc.txt");
	
	 
	 BufferedReader br =  conversion.Utils.getBufferedReader(in);
	 if(br==null){
		 in =  new File(dir, OptionBuild.build+".txt");
		 br =  conversion.Utils.getBufferedReader(in);
	 }
	 String nme = in.getName();
	 nme = nme.substring(0,nme.length()-4);
  File outf = new File(dir, nme+".var.txt");
  PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outf,false)));
 String st = "";
 BufferedReader br1 =  conversion.Utils.getBufferedReader(outf);
 String st1 = br1==null ?  null: br1.readLine();
 //GenesForRegion example = new GenesForRegion();
 //example.initialise();
 //CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
 ZipFile zf = null;
 String chr1 = "";
 String[] res=null;
 double[] d = null;
 int lrr_id=0;
 List<Double> l1 = null;
 boolean skip = false;
 while((st = br.readLine())!=null){
	  if(st1!=null){
		  
		  st1 = br1.readLine();
		  continue;
		 
	  }
	  try{
	  String[] str = st.split("\\s+");
	 
	  String chr = str[0].substring(3);
	  if(!chr.equals(chr1)){
		chr1 = chr;
		if(chr1.indexOf('X')>=0 || chr1.indexOf('Y')>=0 || chr1.indexOf('M')>=0){
			l1 = l1_sex;
		}
		else{
			l1 = l1_auto;
		}
		skip = false;
		if(zf!=null) zf.close();
		File f1 = new File(dir,chr+".zip");
		if(!f1.exists() || f1.length()==0) skip =true;
		else{
		System.err.println("opening "+f1);
		zf = new ZipFile(f1);
		res = new String[Compressor.readZipFrom(zf, "Samples").size()];
		d = new double[res.length];
		List<String > l = Arrays.asList(Compressor.readZipFrom(zf, "Name").get(0).split("\\t+"));
		lrr_id = l.indexOf(lrrst);
		if(lrr_id <0) {
			throw new RuntimeException(l+"\n"+lrrst);
		}
		}
	  }
	  if(skip) continue;
//	  int pos = Integer.parseInt(str[1]);
	  String id = str[3];
	  Compressor.readZip(zf, id, res, lrr_id);
	  double[] meanvar = getVar(res, d);
	l1.add(meanvar[1]);
	  out.println(st+"\t"+String.format("%5.3g",meanvar[0]).trim()
			  +"\t"+String.format("%5.3g",meanvar[1]).trim());
	  }catch(Exception exc){
	  exc.printStackTrace();
		  out.println(st+"\t"+Double.NaN);
	  }
 }
 out.close();
 File outf1 = new File(dir, nme+".hist.txt");
 PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(outf1,false)));
 getHistogram(l1_auto, out1);
 File outf1_sex = new File(dir, nme+".sex.hist.txt");
 PrintWriter out1_sex = new PrintWriter(new BufferedWriter(new FileWriter(outf1_sex,false)));
 getHistogram(l1_sex, out1_sex);

System.err.println();//.getString());
//List l = example.coreDriver.getExternalDatabaseAdaptor().fetch();
//Map<String, List<int[]>> res = example.fetchGenesByLocation(chromosomeLoc, null);
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
}

 
 private static void getHistogram(List<Double> l1, PrintWriter pw) {
	// TODO Auto-generated method stub
	 Collections.sort(l1);
	 double step = 0.01;
	 double nextp=0;
	 nextp+=step;
	for(int k=0; k<l1.size(); k++){
		double perc = (double) k / (double) l1.size();
		if(perc>=nextp){
			pw.println(String.format("%5.3g",perc )+"\t"+String.format("%5.3g",l1.get(k)));
			nextp+=step;
		}
	}
	pw.close();
}

private static double[] getVar(String[] res, double[] d) {
	double mean =0;
	double cnt=0;
	for(int k=0; k<res.length; k++){
		String st = res[k];
		if(st.indexOf('N')>=0) st = "NaN";
		double v = Double.parseDouble(st);
		d[k] = v;
		if(!Double.isNaN(v)){
			mean+=v;
			cnt++;
		}
	}
	mean = mean/cnt;
	double var =0;
	for(int k=0; k<res.length; k++){
		double v = d[k];
		if(!Double.isNaN(v)){
			var+=Math.pow(v-mean,2);
		}
	}
	return new double[] {mean,var/cnt};
}

public void initialise(){
	 Callable run  = new Callable(){
			

			public Object call() {
				
			
// We need core and variation drivers that point to the latest human databases
// on ensembldb.ensmbl.org. The easiest way to get these is from the default registry.
// See org.ensembl.registry.Registry for more information.
  try{
Registry dr = Registry.createDefaultRegistry();

coreDriver = dr.getGroup("human").getCoreDriver();
String[] database  = coreDriver.getConfiguration().getProperty("database").split("_");
String build = database[database.length-1];
String build1 = conversion.OptionBuild.build;
//  if(!build.startsWith(build1.split("ild")[1].substring(0,2))) {
	//   throw new RuntimeException(build+" !! "+build1);
// }
variationDriver = dr.getGroup("human").getVariationDriver();
// Another approach is to use custom configuration files which can point to
// any ensembl databases.
//coreDriver = CoreDriverFactory.createCoreDriver("resources/data/example_core_database.properties");
//variationDriver = VariationDriverFactory.createVariationDriver("resources/data/example_variation_database.properties");
// the variation driver needs a sister core driver to access the core database
//variationDriver.setCoreDriver(coreDriver);


// A third way is to use the user default registry. See org.ensembl.registry.Registry
// for more information.
  }catch(Exception exc){
      exc.printStackTrace();
  }
  return null;
			}
		};
		//return run;
		
		List l = new ArrayList();
		l.add(run);
		  ExecutorService es = Executors.newFixedThreadPool(OptionBuild.numThreads);;
		        try{
		        es.invokeAll(l, 1, TimeUnit.SECONDS);
		        }catch(Exception exc){
		        	exc.printStackTrace();
		        }
		        es.shutdown();
		// System.err.println("");     
 }
 
 /**
  * Creates an instance of Example with core and varition drivers configured.
  * 
  * @throws AdaptorException
  */
 public CalculateVariance() {
	
 }


 
 public static ExecutorService es;
 public void getGC(){
	
 }
 /**
  * Fetches information about the sequence regions.
  * 
  * @throws AdaptorException
  *           if problem occurs during retrieval
  */
 public void fetchSequenceRegionsSuchAsChromosomeOrContig(double[] chrom
    ) throws AdaptorException {
	 CoordinateSystem coordinateSystem = new CoordinateSystem("chromosome");
   // What sequence regions are are available?
   // In the case of the "chromosome" coordinate system these are chromosomes
   SequenceRegion[] seqRegions = coreDriver.getSequenceRegionAdaptor()
       .fetchAllByCoordinateSystem(coordinateSystem);
   System.out.println("There are " + seqRegions.length
       + " sequence regions in the " + coordinateSystem.getName() + "."
       + coordinateSystem.getVersion() + " coordinate system.");

   SequenceRegion sr = seqRegions[0];
   System.out.println(coordinateSystem.getName() + " " + sr.getName()
       + " has length " + sr.getLength());
   StringBuffer sb1 = new StringBuffer();
   StringBuffer sb2 = new StringBuffer();
   for(int i=0; i<seqRegions.length; i++){
	   String name = seqRegions[i].getName();
	   if(name.equals("X")){
		   chrom[22] = seqRegions[i].getLength();
	   }
	   else if(name.equals("Y")){
		   chrom[23] = seqRegions[i].getLength();
	   }
	   else{
		   try{
			   int pos = Integer.parseInt(name);
			   sb1.append(name+":");
			   sb2.append(seqRegions[i].getLength()+":");
			   chrom[pos-1] = seqRegions[i].getLength();
		   }catch(Exception exc){
			   System.err.println(name);
		   }
	   }
   }
   System.err.println(sb1.toString());
   System.err.println(sb2.toString());
  // System.exit(0);
 }
 

 /**
  * Fetches information about karyotypes (chromosome bands) for the specified
  * chromosome.
  * 
  * @param coordinateSystem
  *          coordinate system containing the chromosomeName
  * @param chromosomeName
  *          name of the chromosome of interest
  * @throws AdaptorException
  *           if problem occurs during retrieval
  */
 public List  fetchKaryotypes(
     String chromosomeName) throws AdaptorException {
	 CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
	 System.err.println("fetching "+chromosomeName);
   List l = coreDriver.getKaryotypeBandAdaptor().fetch(chromosomeCS,
       chromosomeName);
   System.out.println("Chromosome " + chromosomeName + " has " + l.size()
       + " karyotypes.");
   
   KaryotypeBand kb = (KaryotypeBand) l.get(0);
   
   Location loc = kb.getLocation();
   int start = loc.getStart();
   int end = loc.getEnd();

   System.out.println("The first karyotype on chromosome " + chromosomeName
       + " is from " + start + "bp to " + end + "bp.");
   return l;
 }
 
public int getCentromere( String chromosomeName, List<String> karya) throws AdaptorException {

	List kary = this.fetchKaryotypes(chromosomeName);
	int qmin = Integer.MAX_VALUE;
	int pmax = 0;
	String curr = "";
	int start =0;
	int end=0;
	for(int i=0; i<kary.size(); i++){
		 KaryotypeBand kb = (KaryotypeBand)kary.get(i);
		 String band = kb.getBand().split("\\.")[0];
		 Location loc = kb.getLocation();
		 end = loc.getEnd();
		 if(!band.equals(curr)){
			 if(curr.length()>0){
				 karya.add(band+"="+start+"-"+end);
			 }
			 curr = band;
			 start = loc.getStart();
		 }
		
		 boolean centLeft = i>=1 && kb.getBand().startsWith("q") && ((KaryotypeBand)kary.get(i-1)).getBand().startsWith("p");
		 boolean centRight = i<kary.size()-1 && kb.getBand().startsWith("p") && ((KaryotypeBand)kary.get(i+1)).getBand().startsWith("q");
		
		if(centLeft && loc.getEnd()>pmax){
			pmax = loc.getEnd();
		}
		if(centRight && loc.getStart() < qmin){
			qmin = loc.getStart();
		}
	}
	 karya.add(curr+"="+start+"-"+end);
	if(pmax>0) return pmax +1;
	else if(qmin <Integer.MAX_VALUE){
		return qmin -1;
	}
	else return pmax+1;
}

public static void writeKaryotypeFile(File karyotypes) throws Exception{
	CalculateVariance g = new CalculateVariance();
	   g.initialise();
	   PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(karyotypes)));
	   String[] st = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y".split(":");
	   for(int k=0; k<st.length; k++){
		   List<String> l = new ArrayList<String>();
		   int centromere= g.getCentromere(st[k],l);
		   pw.println(st[k]+"\t"+centromere+"\t"+l);
	   }
	   pw.close();
	
}
 
 
 /* returns start/end location */
 

 
 
} // Example