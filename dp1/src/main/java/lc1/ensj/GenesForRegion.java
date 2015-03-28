package lc1.ensj;

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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;

import lc1.dp.data.collection.DataCollection;
import lc1.util.Constants;
import lc1.util.Executor;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Exon;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.KaryotypeBand;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.ExternalDatabaseAdaptor;
import org.ensembl.registry.Registry;
import org.ensembl.variation.driver.VariationDriver;


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

public class GenesForRegion {

 private CoreDriver coreDriver;

 private VariationDriver variationDriver;

 
 public static void main(String[] args) {
	 try{
	 File in = new File("build36.txt");
	 BufferedReader br =  DataCollection.getBufferedReader(in);
   GenesForRegion example = new GenesForRegion();
   PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("out.txt")));
  String st = "";
  
 
  CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
  while((st = br.readLine())!=null){
	  String[] str = st.split("\t");
	  String chr = str[0].substring(3);
	  int pos = Integer.parseInt(str[1]);
	  Location chromosomeLoc = new Location("chromosome:"+chr+":"+(pos-150)+"-"+(pos+150));
	  String str1 = example.coreDriver.getSequenceAdaptor().fetch(chromosomeLoc).getString();
	  int gc =0;
	  int at=0;
	  for(int i=0; i<str1.length(); i++){
		  char ch = str1.charAt(i);
		  if(ch=='G' || ch=='C') gc++;
		  else if(ch=='A' || ch=='T') at++;
	  }
	  double gcc = ((double)gc)/ ((double)gc + (double)at);
	  out.println(st+"\t"+gcc);
  }
  out.close();
  

 System.err.println();//.getString());
//List l = example.coreDriver.getExternalDatabaseAdaptor().fetch();
 //Map<String, List<int[]>> res = example.fetchGenesByLocation(chromosomeLoc, null);
	 }catch(Exception exc){
		 exc.printStackTrace();
	 }
 }

 /**
  * Creates an instance of Example with core and varition drivers configured.
  * 
  * @throws AdaptorException
  */
 public GenesForRegion() {
	//Callable run  = new Callable(){
			

		//	public Object call() {
				
			
   // We need core and variation drivers that point to the latest human databases
   // on ensembldb.ensmbl.org. The easiest way to get these is from the default registry.
   // See org.ensembl.registry.Registry for more information.
     try{
   Registry dr = Registry.createDefaultRegistry();
//    coreDriver =
	//   CoreDriverFactory.createCoreDriver("ensembldb.ensembl.org",
		//		  3308, "homo_sapiens_core_54_36p", "anonymous", null);

//   variationDriver =   VariationDriverFactory.createVariationDriver("ensembldb.ensembl.org",
	//		  3306, "homo_sapiens_variation_50_36l", "anonymous", null);
   //.getGroup("human").getVariationDriver();
  coreDriver = dr.getGroup("human").getCoreDriver();
   String[] database  = coreDriver.getConfiguration().getProperty("database").split("_");
   String build = database[database.length-1];
   String build1 = Constants.build(0);
   if(!build.startsWith(build1.split("ild")[1].substring(0,2))) throw new RuntimeException(build+" !! "+build1);
 
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
     //return null;
		//	}
	//	};
	/*	List l = new ArrayList();
		l.add(run);
		  ExecutorService es = 
		      Executor.getEs(GenesForRegion.class, 1);
		        try{
		        es.invokeAll(l, 1, TimeUnit.SECONDS);
		        }catch(Exception exc){
		        	exc.printStackTrace();
		        }*/
		 System.err.println("");     
 }

 

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
 
 public List[] fetchKaryotypes() throws AdaptorException{
	List[] str = new List[24];
	 for(int i=0; i<str.length; i++){
		 String st  = getString(i+1);
		 str[i] = this.fetchKaryotypes(st);
	 }
	 return str;
 }

public static String getString(int o1){
	if(o1==23) return "X";
	else if(o1==24) return "Y";
//	else if(o1==25) return "";
//	else if(o1==26) return "M";
	else return o1+"";
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
   if(l.size()>0){
   KaryotypeBand kb = (KaryotypeBand) l.get(0);
   
   Location loc = kb.getLocation();
   int start = loc.getStart();
   int end = loc.getEnd();

   System.out.println("The first karyotype on chromosome " + chromosomeName
       + " is from " + start + "bp to " + end + "bp.");
   }
   return l;
 }
 
public int getCentromere( String chromosomeName) throws AdaptorException {
	List kary = this.fetchKaryotypes(chromosomeName);
	int qmin = Integer.MAX_VALUE;
	int pmax = 0;
	for(int i=0; i<kary.size(); i++){
		 KaryotypeBand kb = (KaryotypeBand)kary.get(i);
		 Location loc = kb.getLocation();
		 boolean centLeft = i>=1 && kb.getBand().startsWith("q") && ((KaryotypeBand)kary.get(i-1)).getBand().startsWith("p");
		 boolean centRight = i<kary.size()-1 && kb.getBand().startsWith("p") && ((KaryotypeBand)kary.get(i+1)).getBand().startsWith("q");
		if(centLeft && loc.getEnd()>pmax){
			pmax = loc.getEnd();
		}
		if(centRight && loc.getStart() < qmin){
			qmin = loc.getStart();
		}
	}
	if(pmax>0) return pmax +1;
	else if(qmin <Integer.MAX_VALUE){
		return qmin -1;
	}
	else return pmax+1;
}
 
 
 /* returns start/end location */
 public SortedMap<String, List<int[]>> fetchGenesByLocation(final Location location, final Map<String, String>nmes) throws AdaptorException {

final SortedMap<String,List< int[]>> l = new TreeMap<String, List<int[]>>();
			if(coreDriver==null) return l;
			Callable call = new Callable(){
				
			public Object call(){
				try{
			     Iterator geneIterator = coreDriver.getGeneAdaptor().fetchIterator(location,
			         true);
			     ExternalDatabaseAdaptor edb = coreDriver.getExternalDatabaseAdaptor();
			   
			     while (geneIterator.hasNext()) {
			         Gene gene  = (Gene) geneIterator.next();
			         String nme = gene.getDisplayName();
			         String id =  gene.getAccessionID();
			       /* List l1 =  gene.getExternalRefs();
			        for(int i=0; i<l1.size(); i++){
			        	ExternalRefImpl eri = (ExternalRefImpl) l1.get(i);
			        	String name = eri.getExternalDatabase().getName();
			        	System.err.println(name);
			        }
			   
			      List l2 = edb.fetch();*/
			      if(nmes!=null) nmes.put(nme, id);
			         List<int[]> res = new ArrayList<int[]>();
			         if(nme==null){
			            nme =  gene.getAccessionID();
			             
			         }
			         l.put(nme, res);
			         List exons =(gene).getExons();
			         for(int i=0; i<exons.size(); i++){
			             Exon ex = (Exon) exons.get(i);
			            Location loc = ex.getLocation();
			            res.add( new int[] {loc.getStart(),loc.getEnd()});
			           
			         }
			     }
				}catch(Exception exc){
					exc.printStackTrace();
				}
			     return null;
			}
			
			
			};
			List l1 = new ArrayList();
			l1.add(call);
			  ExecutorService es = 
			      Executor.getEs(GenesForRegion.class, 1);
			        try{
			        es.invokeAll(l1, 2, TimeUnit.SECONDS);
			        }catch(Exception exc){
			        	exc.printStackTrace();
			        }
			     return l;
   }
 
 
} // Example