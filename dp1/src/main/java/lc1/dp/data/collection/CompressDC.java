package lc1.dp.data.collection;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.Compressor;
import lc1.util.Constants;

public class CompressDC {

    
 ZipFile originalZip = null;
	
    File dir;
    FileOutputStream dest;
    CheckedOutputStream checksum;
    ZipOutputStream outS;
    OutputStreamWriter osw;
  //  PrintWriter buildPw;
    PrintWriter sampw;
    ///String chr;
    DataCInterface dc;
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
    final PrintWriter plate_pw,pheno_pw;
    public CompressDC(File dir, DataCInterface dc) throws Exception{
      File origDir = new File(Constants.baseDir+"/"+dc.name());
      
        File zf = new File(Constants.baseDir+"/"+dc.name()+"/"+Constants.chrom0()+".zip");
        if(zf==null || ! zf.exists()){
        	zf = DataCollection.getKaryoFile(zf, dc.loc().get(0));
        }
        this.originalZip = new ZipFile(zf);
        System.err.println("dir is "+dir.getAbsolutePath());
        this.dc = dc;
        if(false && (Constants.writeExtractedFile && ! dc.name().endsWith("M")
        		&& ! dc.name().endsWith("k") && dc.name().indexOf("sanger")<0
        		&& dc.name().toLowerCase().indexOf("hapmap")<0)){
        	((DataCollection)dc).name = "study"+((DataCollection)dc).index+"";
        	//  dir = new File(dir.getParentFile(), dc.name().split("_")[0]);
        	  this.anonymise = true;
        }
        else{
        	this.anonymise = false;
        }
        this.dir = dir;
        //   this.delete(dir);
           dir.mkdir();
         
           
    //    buildPw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(dir, Constants.build()+".gz"))))));
        
        //indiv = dc.indiv();
           plate_pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "plate.txt"))));
           pheno_pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "pheno.txt"))));
           
           dest = new FileOutputStream(new File(dir, dc.chrom()+".zip"));
        checksum = new   CheckedOutputStream(dest, new Adler32());
        outS = new 
        ZipOutputStream(new 
          BufferedOutputStream(checksum));
        osw = new OutputStreamWriter(outS);
        outS.setMethod(ZipOutputStream.DEFLATED);
        plat = new File(origDir, "plate.txt");
        phen = new File(origDir, "pheno.txt");
    }
    File plat,phen;
   
    public void writeHeader() throws Exception{
    	List<String> l = null;
    	if(this.originalZip!=null){
    		l = Compressor.getEntries(originalZip, "Name");
    	}
    	
        ZipEntry headings = new ZipEntry("Name");
        outS.putNextEntry(headings);
      String headSNP =!this.writeDistribution ?  "Log R\tB allele frequency"  : dc.headSNP();
        osw.write(headSNP+"\n");
        if(l!=null){
        	osw.write(l.get(1)+"\n");
        	osw.write(l.get(2)+"\n");
        }
        else{
        osw.write(dc.head_snp()+"\n");
        osw.write("id\tdata_index");
        }
       // Phenotypes pheno = this.dc.pheno;
        StringBuffer fs= new StringBuffer();
      /*  if(pheno!=null && pheno.phen.size()>0){
            fs.append("%10s");
            for(int i=0; i<pheno.phen.size()-1; i++){
                fs.append("\t%10s");
            }
            String fstring = fs.toString();
            osw.write("\t"+Format.sprintf(fstring, pheno.phen.toArray()));
        }*/
       osw.write("\n");
        osw.flush();
        outS.closeEntry();
    }
    public void write(Character o) throws Exception{
        if(o==null) osw.write("null");
        else osw.write(o);
    }
    public void writeSnps() throws Exception{
    	   ZipEntry headings = new ZipEntry("SNPS");
    	   outS.putNextEntry(headings);
    	if(this.originalZip!=null){
    		List<String> l = Compressor.getIndiv(originalZip, "SNPS",3,"\\s+");
    		int first = l.indexOf(dc.snpid().get(0));
    		int last = l.indexOf(dc.snpid().get(dc.snpid().size()-1));
    		String st = "";
    		BufferedReader br = Compressor.getBufferedReader(originalZip, "SNPS");
    		for(int i=0; i<first; i++){
    			br.readLine();
    		}
    		for(int i= first;i<last; i++){
    			String[] str = br.readLine().split("\t");
    			if(false){str[3] = this.dc.name()+i;
    			
    			str[2] = ""+(Integer.parseInt(str[2])+10101010);
    			str[1] = ""+(Integer.parseInt(str[1])+10101010);
    			str[0] = "chr8";
    			}
    			for(int k=0; k<str.length; k++){
    				   	osw.write(str[k]+(k==str.length-1 ? "\n":"\t" ));
    			}
    			//osw.write(str+"\n");
    		}
    		br.close();
    		
    	}
    	else{
        
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
    	}
            osw.flush();
           outS.closeEntry();
         //  buildPw.close();
    }
   public String getDataType(String key){
	   if(dc instanceof MergedDataCollection){
		   Set<Integer> s = new HashSet<Integer>();
		   int data_index = ((HaplotypeEmissionState)((MergedDataCollection)this.dc).dataL.get(key)).dataIndex();
		   if(data_index==-1)
			   return "-1";
		   return ((MergedDataCollection)dc).ldl[data_index].name;
	   }
	   else return dc.name();
   }
   public final boolean anonymise; 
   public void writeSamples(List<List<String>> stl, Map<String, String> rename) throws Exception{
        ZipEntry headings = new ZipEntry("Samples");
       // Phenotypes pheno = this.dc.pheno;
        StringBuffer fs= new StringBuffer();
       
      //  String fstring = fs.toString();
        outS.putNextEntry(headings);
        List<String> l =null;
        if(originalZip!=null){
        	l = Compressor.getEntries(originalZip, "Samples");
        }
       
        
            for(int i=0; i<l.size(); i++){
                String key = l.get(i);
                if(l!=null){
                	String str_ = l.get(i);
                	if(!str_.startsWith(key)){
                		throw new RuntimeException("!!");
                	}
                	String[] str = str_.split("\t");
                	if(this.anonymise){	
                		String value = dc.name()+"-"+i;
                		rename.put(str[0], value);
                		str[0] = value;
                		
                	}
                	
                	for(int k=0; k<str.length; k++){
                		osw.write(str[k]+(k==str.length-1 ? "\n":"\t"));
                	}
                	//osw.write("\n");
                }
                else{
                    osw.write(key);
                    osw.write("\t"+getDataType(key));
                    
                
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
   boolean writeDistribution =true;
   // List<String> indiv;
    public  void run(List<List<String>> indiv1, int start, int end) throws Exception{
    	System.err.println("writing compressed file "+this.dir.getName());
       writeHeader();
       //BufferedReader br =DataCollection.getBufferedReader(f);
      // String st = "";
       List<String> indiv11 = Compressor.getIndiv(this.originalZip, "Samples",0 );
       List<String> indiv_ = indiv1.get(0);
      for(int i=0; i<this.dc.loc().size(); i++){
          String snp_id = dc.snpid().get(i);
          ZipEntry headings = new ZipEntry(snp_id);
          outS.putNextEntry(headings);
      // 	System.err.println(i+" "+dc.snpid().size());
        //  for(int ik=0; ik<indiv1.size(); ik++){
        	//  List<String> indiv11 = indiv1.get(ik);
          for(int k=0; k<indiv11.size(); k++){
              String key = indiv11.get(k);
           //  PseudoDistribution hes_ =  ((HaplotypeEmissionState)((DataCollection)dc).data.get(key)).emissions[i];
            
              String hes = indiv_.contains(key) ? dc.getCompressedString(key, i,true,false) : "NaN";
              String hes1 =indiv_.contains(key) ? dc.getCompressedString(key, i,false,false) : "NaN";
           /*  EmissionState st1 = ((DataCollection)dc).dataL.get(key);
             if(key.equals("21805")){
            	 Logger.global.info("h");
             }*/
              if(writeDistribution){
              osw.write( hes1);
              osw.write(hes==null  ? "" :"\t"+hes);
              }
              else{
            	  osw.write(hes);
              }
              osw.write("\n");
            
          }
        //  }
          osw.flush();
           outS.closeEntry();
        
          
       
      }
       osw.flush();
       outS.closeEntry();
        writeSnps();
        Map<String, String > m = new HashMap<String,String>();
       writeSamples(indiv1,m);
       PrintWriter tr = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "trans.txt"))));
       for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
    	   String key = it.next();
    	   String value = m.get(key);
    	   tr.print(key+"\t"+value);
       }
       tr.close();
       write(plat, plate_pw,  "PLATE",new HashMap(m));
       write(phen, pheno_pw,  "pheno",m);
        outS.close();
     	System.err.println("done writing compressed file "+this.dir.getName());
        
    }
private void write(File plat, PrintWriter plate_pw, String val, Map<String, String> m) throws Exception{
	 if(plat.exists() ){
	       	BufferedReader br = new BufferedReader(new FileReader(plat));
	       
	       	List<String> l = Arrays.asList(br.readLine().split("\t"));
	       	plate_pw.println("PATIENT\t"+val);
	       	int i1 = l.indexOf("PATIENT");
	       	int i2 = l.indexOf(val);
	       	if(i2<0){
	       		i2 = l.size()-1;
	       	}
	       	String st = "";
	       	while((st = br.readLine())!=null){
	       		String[] str = st.split("\t");
	       		if(str[i2].toLowerCase().equals("na")){
	       			str[i2] = "1";
	       		}
	       		if(anonymise){
	       			String st1 = m.remove(str[i1]);
	       			if(this.dir.getName().equals("SLEGEN")){
	       			if(st1==null || st1.equals("null")){
	       			
	       				st1 = getClosest(m.keySet(),str[i1].split("-Inf")[0]);
	       				if(st1!=null) st1 = m.remove(str[i1]);
	       				System.err.println("warning null "+str[i1]);
	       			}
	       			}
	       			if(st1!=null){
	       		plate_pw.println(st1+"\t"+str[i2]);
	       			}
	       		}
	       		else{
	       			plate_pw.println(str[i1]+"\t"+str[i2]);
	       		}
	       	}
	       }
	      
	       for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
	    	   String key = it.next();
	    	   if(key.toLowerCase().indexOf("cont")>=0){
	    		   plate_pw.println(m.get(key)+"\t1");
	    		   it.remove();
	    	   }
	       }
	       plate_pw.close();
	       System.err.println("left "+m.keySet());

	
}
private String getClosest(Set<String> keySet, String string) {
	for(Iterator<String> it = keySet.iterator(); it.hasNext();){
		String nxt = it.next();
		if(nxt.startsWith(string)){
			return nxt;
		}
	}
	return null;
}
}
