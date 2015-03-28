package lc1.dp.model;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.CachedEmissionState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.WrappedEmissionState1;
import lc1.stats.PseudoDistribution;

public class CompressHMM {

    
 
    File dir;
    FileOutputStream dest;
    CheckedOutputStream checksum;
    ZipOutputStream outS;
    OutputStreamWriter osw;
  //  PrintWriter buildPw;
    PrintWriter sampw;
    ///String chr;
    HaplotypeHMM hmm ;
    public void delete(File dir){
    	if(dir!=null && dir.exists()){
    	File[] f = dir.listFiles();
    	
    	for(int i=0; i<f.length; i++){
    		f[i].delete();
    	}
    	dir.delete();
    	}
    	
    }
  /*String[] header =   new String[] {
            "Genotype\tScore",
            "id",
            "id"};*/
    public CompressHMM(File dir, HaplotypeHMM hmm, DataCollection dc) throws Exception{
        this.dir = dir;
        this.dc = dc;
     //   this.delete(dir);
        dir.mkdir();
        System.err.println("dir is "+dir.getAbsolutePath());
        this.hmm = hmm;
    //    buildPw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(dir, Constants.build()+".gz"))))));
        
        //indiv = dc.indiv();
        dest = new FileOutputStream(new File(dir, hmm.name+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
    }
   
    public void writeHeader() throws Exception{
      /*  ZipEntry headings = new ZipEntry("Name");
        outS.putNextEntry(headings);
      String headSNP = dc.headSNP();
        osw.write(headSNP+"\n");
        
        osw.write(dc.head_snp()+"\n");
        osw.write("id\tdata_index");
    
        StringBuffer fs= new StringBuffer();
     
       osw.write("\n");
        osw.flush();
        outS.closeEntry();*/
    }
    public void write(Character o) throws Exception{
        if(o==null) osw.write("null");
        else osw.write(o);
    }
   final DataCollection dc;
    public void writeSnps() throws Exception{
        ZipEntry headings = new ZipEntry("SNPS");
        outS.putNextEntry(headings);
            for(int i=0; i<dc.snpid().size(); i++){
            //	System.err.println(i+" "+dc.snpid().size());
                osw.write("chr"+dc.chrom()+"\t");
                osw.write(dc.loc().get(i)+"\t");
                osw.write((dc.loc().get(i)+40)+"\t");
                osw.write(dc.snpid().get(i)+"\t");
                
            //    buildPw.write("chr"+dc.chrom()+"\t");
             //   buildPw.write(dc.loc().get(i)+"\t");
              //  buildPw.write((dc.loc().get(i)+40)+"\t");
             //   buildPw.write(dc.snpid().get(i)+"\t");
                
                
                if(dc.alleleA()!=null && dc.alleleA().size()>0 && dc.alleleB()!=null && dc.alleleB().size()>0){
                    Character maj = dc.alleleA().get(i);
                    Character min = dc.alleleB().get(i);
                    if(maj==null) maj = 'N';
                    if(min==null) min = 'N';
                    osw.write(maj); osw.write("\t");
                    osw.write(min);
                    
                 //   buildPw.write(maj+"\t");
                  //  buildPw.write(min);
                }//
                osw.write("\n");
               // buildPw.write("\n");
            }
            osw.flush();
           outS.closeEntry();
         //  buildPw.close();
    }
  /* public String getDataType(String key){
	   if(dc instanceof MergedDataCollection){
		   Set<Integer> s = new HashSet<Integer>();
		   int data_index = ((HaplotypeEmissionState)((MergedDataCollection)this.dc).dataL.get(key)).dataIndex();
		   if(data_index==-1)
			   return "-1";
		   return ((MergedDataCollection)dc).ldl[data_index].name;
	   }
	   else return dc.name();
   }*/
    public void writeSamples() throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
       // Phenotypes pheno = this.dc.pheno;
        StringBuffer fs= new StringBuffer();
       
      //  String fstring = fs.toString();
        outS.putNextEntry(headings);
        for(int ik=1; ik<hmm.states.size(); ik++){
      	  if(hmm.states.get(ik) instanceof EmissionState){
      	EmissionState state = (EmissionState)hmm.states.get(ik);
      	  EmissionStateSpace stsp = state.getEmissionStateSpace();
      	  Integer nocop = state.noCop();
      	  int[] indices = stsp.getGenoForCopyNo(nocop);
      	  osw.write(nocop.toString());
      	  for(int k=0; k<indices.length; k++){
      		  osw.write("\t"+stsp.get(indices[k]).toString());
      	  }
      	  osw.write("\n");
      	  }
        }
        
            osw.flush();
           outS.closeEntry();
    
    }
    
    
    
    
  
   // List<Integer> loc1 = new ArrayList<Integer>();
  //  List<String> major = new ArrayList<String>();
 //   List<String> minor =  new ArrayList<String>();
  // List<String> list1 = new ArrayList<String>();
   //List<String> list2 = new ArrayList<String>();
   
   // List<String> indiv;
    public  void run() throws Exception{
    	System.err.println("writing compressed hmm "+this.dir.getName());
       writeHeader();
       //BufferedReader br =DataCollection.getBufferedReader(f);
      // String st = "";
      for(int i=0; i<this.dc.loc().size(); i++){
          String snp_id = dc.snpid().get(i);
          ZipEntry headings = new ZipEntry(snp_id);
          outS.putNextEntry(headings);
      // 	System.err.println(i+" "+dc.snpid().size());
          for(int ik=1; ik<hmm.states.size(); ik++){
        	  if(hmm.states.get(ik) instanceof EmissionState){
        	 EmissionState state = (EmissionState)hmm.states.get(ik);
        	  EmissionStateSpace stsp = state.getEmissionStateSpace();
        	  int nocop = state.noCop();
        	  int[] indices = stsp.getGenoForCopyNo(nocop);
        	  PseudoDistribution dist;
        	  if( state instanceof HaplotypeEmissionState) 
        		  dist = ((HaplotypeEmissionState)state).emissions[i];
        	  else if(state instanceof WrappedEmissionState1){
        		  dist =   ((WrappedEmissionState1)state).emissions(i);//.emissions[i] ;
        	  }
        	  else{
        		dist =   ((CachedEmissionState)state).emissions[i] ;
        	  }
        	  osw.write(String.format("%5.3g", dist.probs(indices[0])).trim());
        	  for(int k=1; k<indices.length; k++){
        		  osw.write("\t"+String.format("%5.3g", dist.probs(indices[k])).trim());
        	  }
        	osw.write("\n");
        	  }
        	
          }
          osw.flush();
           outS.closeEntry();
        
          
       
      }
       osw.flush();
       outS.closeEntry();
        writeSnps();
       writeSamples();
        outS.close();
     	System.err.println("done writing compressed file "+this.dir.getName());
        
    }
}
