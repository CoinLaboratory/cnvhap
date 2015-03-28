package conversion;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

public class SplitByPheno {
public static void main(final String[] args){
	File user = new File(System.getProperty("user.dir"));
	/* File excel = user.listFiles(new FileFilter(){

	public boolean accept(File pathname) {
	return pathname.getName().endsWith(".xls");
	}
	 
 })[0];*/
	String[] ar = args[0].split(":");
	for(int i=0; i<ar.length; i++){
		split(user,new File(user,"plate.txt"), ar[i], "PLATE", "PATIENT", OptionBuild.sampleId1);
	}
	
}
public static File[] split(File user, File excel, final String chrom, final String col_header,
		final String id_header, String[] id_header1){
	try{ 
	//	System.err.println("running sep!");
    File[] f = user.listFiles(new FileFilter(){

        public boolean accept(File pathname) {
           return pathname.getName().equals(chrom+".zip") 
         //  &&   !pathname.getName().endsWith("zip")
         ;
        }
        
    });
    List<String> types = new ArrayList<String>();
	
	  Map<String, Integer> m =  readCaseControl(excel, types,col_header, id_header, id_header1);
		types.add("null");
	//	System.err.println("plate info "+m);
	  File[] outF = getOutFiles(user, types);
	  File data_f = new File(user, "data.xls");
	 /* if(!data_f.exists() && chrom.equals("1") ){
		  try{
			WritableWorkbook workbook = Workbook.createWorkbook(data_f);
			WritableSheet sh1 = workbook.createSheet("0", 0);
			String[] rows = new String[] {
					"include", "inputDir", "format", "noCopies", "indexToRestrict", 
					"maxIndiv","toDel", "exponentB", "exponentR", "toInclude", 
					"r_state_mean", "r_state_var", "r_state_skew","r_prior", "fillLikelihood", 
					"restrictIndivTo"};
			
			String[] rows1 = new String[] {
					"TRUE", "1M", "illumina", "2", "FALSE",
					"100", "null","1","1","all",
					"-3.0;1.0;0.0", "0.12;0.12;0.12", "-1.0;1.0;0.0", "1;1;1", "0",
					"null"};
			for(int i=0; i<rows.length; i++){
				WritableCell c  = new Label(0, i, rows[i]);
				sh1.addCell(c);
				for(int j=0; j<outF.length; j++){
					WritableCell c1  = new Label(j+1, i, i==1 ? outF[i].getName() : rows1[i]);
					sh1.addCell(c1);
				}
			}
			workbook.write();
			workbook.close();
		  }catch(Exception exc){
			  System.err.println("failed to create "+data_f.getName());
			  exc.printStackTrace();
		  }
		}*/
    for(int i=0; i<f.length; i++){
    	
         split(f[i], outF,m, types, i==0 );
    }
   return outF;
	}catch(Exception exc){
		exc.printStackTrace();
		return null;
	}
}
/*
public static int getHeaderIndex(Sheet ws, String val){
	
	 for(int k=0;k<ws.getColumns(); k++){
		 String cnt =ws.getCell(k, 0).getContents(); 
		 if(cnt.equals(val)) return k;
	 }
	 return -1;
}*/
/* (non-Javadoc)
 * @see lc1.dp.data.collection.DataC#readCaseControl(java.io.File)
 * 
 *  
    
  desir are the controls!
 */
public static Map<String, Integer> readCaseControl(File cca, 
      List<String> type,
        String val, String id_header, String[] id_header1) throws Exception{
	BufferedReader br = Utils.getBufferedReader(cca);
	List<String> header = Arrays.asList(br.readLine().split("\t"));
	
	
	 int k=header.indexOf(val);
	 int k1=header.indexOf(id_header);
	 if(k<0) throw new RuntimeException("no index of "+val+"\n"+header);
	 if(k1<0) throw new RuntimeException("no index of "+id_header+"\n"+header);
	int[] k2 = new int[id_header1.length];
	for(int jk=0; jk<k2.length; jk++){
		k2[jk] = 	header.indexOf(id_header1[jk]);
		 if(k2[jk]<0) throw new RuntimeException("no index of "+id_header1[jk]+"\n"+header);
	}

    HashMap<String, Integer> cc = new HashMap<String, Integer>();
    String st;
    Set<String> toappendSentrix = new HashSet<String>();
    for(int j=1; (st = br.readLine())!=null; j++){
    	String[]  str = st.split("\t");
    	System.err.println(st);
        String val_j = str[k];
        String id_j =str[k1];
        if(toappendSentrix.contains(id_j)){
        	for(int jk=0;jk<id_header1.length;jk++){
        		id_j = id_j+"_"+str[k2[jk]];
        	}
        }
        int index = type.indexOf(val_j);
        if(index<0){
        	index = type.size();
        	type.add(val_j);
        	
        }
         if(cc.containsKey(id_j) && !cc.get(id_j).equals(index)){
        	 type.remove(index);
        	 toappendSentrix.add(id_j);
        	 Integer value = cc.remove(id_j); 
        	 for(int jk=0;jk<id_header1.length;jk++){
         		id_j = id_j+"_"+str[k2[jk]];
         	}
        	 cc.put(id_j, value);
        	 //System.err.println("WARNING dupl "+id_j);
        	 //Integer prev = cc.remove(id_j);
        	//throw new RuntimeException("!!");
        	 
         }
         else{
             cc.put(id_j, index);
         }
       //System.err.println(id_j);
        //if(cc.size()!=this.dataL.size()) throw new RuntimeException("!!!");
        
    }
   // int sze = cc.size();
  // Integer v =  cc.get("EPG796224");
    return cc;

}



public static Integer get(Map<String, Integer>m, String st, List<String> types, String plate_inf){
	Integer res = m.get(st);
	if(res==null){
		int ind = st.lastIndexOf('_');
		while(res==null && ind>=0){
			String substr = st.substring(0, ind);
			res = m.get(substr);
			ind = substr.lastIndexOf('_');
		}
	}
	if(res==null ){
		int best_ind =-1;
		int bestMatch = -1;
		for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
			String id = it.next();
			
			String st1 = id.split("_")[0];
			if(st1.equals(st)){
				int ind = m.get(id);
				String plate = types.get(ind);
				int totMatch = plate_inf==null ? 0 : getMatch(plate.toLowerCase().toCharArray(), plate_inf.toLowerCase().toCharArray());
				if(totMatch>bestMatch){
					best_ind = ind;
					bestMatch = totMatch;
				}
			}
		}
		if(best_ind>=0){
			res = best_ind;
		}
	}
	return res;
}
private static int getMatch(char[] plate, char[] plate_inf) {
	int len1 = plate.length;
	int len2 = plate_inf.length;
	
	int cnt=0;
	int k1 = 0;
	int k2 = 0;
	while(k1<len1 && k2 <len2){
		if(plate[k1]==plate_inf[k2]) cnt++;
		k1++; k2++;
	}
	 k1 = len1-1;
	 k2 = len2 -1;
	while(k1>0 && k2>0){
		if(plate[k1]==plate_inf[k2]) cnt++;
		k1--; k2--;
	}
	return cnt;
	
}
public static void split(File in, File[] out, final Map<String, Integer> indices, List<String> types, boolean dosample) throws Exception{
    FileOutputStream[] dest = new FileOutputStream[out.length];
    CheckedOutputStream[] checksum = new CheckedOutputStream[out.length];
    ZipOutputStream[] outS = new ZipOutputStream[out.length];
    OutputStreamWriter[] osw = new OutputStreamWriter[out.length];
    List[] resu = new List[out.length];
    for(int i=0; i<out.length; i++){
    	resu[i] = new ArrayList();
    }
    ZipFile zf = new ZipFile(in);
    if(in.length()==0) return;
    String chr = in.getName().toString().split("\\.")[0];
    List<String> snps = Compressor.getIndiv(zf, "SNPS");
    List<String> indiv = Compressor.getIndiv(zf, "Samples",0);
   List<String>  plate1 = Compressor.getIndiv(zf, "Samples",1);
    List<String> names =Compressor.getIndiv(zf, "Name");
    if(snps==null){
    	 snps = Compressor.getIndiv(zf, chr+"/SNPS");
    }
    if(indiv==null){
    	  indiv = Compressor.getIndiv(zf, chr+"/Samples",0);
    }
    if(names==null){
    	   names =Compressor.getIndiv(zf, chr+"/Name");
    }
    if(indiv==null){
    	 if(indiv==null) Logger.global.info("could not find entry Samples "+in.getAbsolutePath()+" "+in.exists());
    	indiv = Compressor.getIndiv(new File(in.getParentFile(), "Samples"), 0);
    }
    String head = in.getName();
    Set<Integer> poss_vals = new HashSet<Integer>();
        for(int j=0; j<indiv.size(); j++){
            poss_vals.add(get(indices, indiv.get(j).split("#")[0], types, plate1.get(j)));
        }
    head = head.substring(0, head.indexOf(".zip"));
    for(int i=0; i<out.length; i++){
        if(poss_vals.contains(i)){
            if(!out[i].exists()) out[i].mkdir();
            File outF = new File(out[i], head+".zip");
            dest[i] = Compressor.getOS(outF);
            checksum[i] = new CheckedOutputStream(dest[i], new Adler32());
            outS[i] = new 
             ZipOutputStream(new 
               BufferedOutputStream(checksum[i]));
           osw[i] = new OutputStreamWriter(outS[i]);
           outS[i].setMethod(ZipOutputStream.DEFLATED);
          if(names!=null) Compressor.writeEntry("Name", names, outS[i], osw[i]);
          if(snps!=null)  Compressor.writeEntry("SNPS", snps, outS[i], osw[i]);
           resu[i] = new ArrayList();
        }
    }
    List<Integer> alias = new ArrayList<Integer>(indiv.size());
    for(int i=0; i<indiv.size(); i++){
        String ind = indiv.get(i);
        Integer index = get(indices,ind.split("#")[0], types, plate1.get(i));
        if(index==null){
        	
        	System.err.println("null index for "+ind);
        	indices.put(ind.split("#")[0], index = types.size()-1);
        	
        }
     //   System.err.println(ind+" "+index+" "+types.get(index));
        alias.add(index);
    }
    
    
    boolean doneSamp = false;
    for(Enumeration en = zf.entries(); en.hasMoreElements();){
        ZipEntry ent = (ZipEntry) en.nextElement();
        if(ent.getName().startsWith("Name") || ent.getName().startsWith("SNPS") || ent.getName().endsWith("SNPS")) continue;
        List<String> res = Compressor.getIndiv(zf, ent.getName());
        Compressor.subList(res, alias, resu);
        if(ent.getName().indexOf("Samples")>=0) doneSamp= true;
        for(int i=0; i<out.length; i++){
            if(outS[i]==null || resu[i].size()==0) continue;
         //   System.err.println(out[i].getName()+" "+types.get(i)+" "+resu[i].size());
            Compressor.writeEntry(ent.getName(), resu[i],outS[i], osw[i]);
            resu[i].clear();
        }
    }
  if(!doneSamp){
	  File sampleFile = new File(in.getParentFile(), "Samples");
	  List<String> res = Compressor.getIndiv(sampleFile, 0);
	  Compressor.subList(res, alias, resu);
    
      for(int i=0; i<out.length; i++){
          if(outS[i]==null || resu[i].size()==0) continue;
       //   System.err.println(out[i].getName()+" "+types.get(i)+" "+resu[i].size());
          Compressor.writeEntry("Samples", resu[i],outS[i], osw[i]);
          resu[i].clear();
      }
  }
    String sample_header = null;
    File sampleFile = new File(in.getParentFile(), "Samples.txt");
 //   if(!sampleFile.exists()) sampleFile = new File(in.getParentFile(), "Samples");
    if(dosample && sampleFile.exists() && sampleFile.length()>0){
        List<String> samples = new ArrayList<String>();
    
    BufferedReader br =Utils.getBufferedReader(sampleFile);
    
   sample_header =br.readLine(); 
    Utils.readPosInfo(br,new int[] { 0}, false, new List[] {samples}, new Class[] {String.class}, "________________");
    Compressor.subList(samples, alias, resu);
    br.close();
    }
    for(int i=0; i<out.length; i++){
        if(outS[i]!=null) outS[i].close();
        zf.close();
        if(sample_header!=null && outS[i]!=null){
            PrintWriter sample_pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(out[i], "Samples.txt"))));
            sample_pw.println(sample_header);
            for(int j=0; j<resu[i].size(); j++){
                sample_pw.println(resu[i].get(j));
            }
            sample_pw.close();
        }
    }
    System.err.println("done "+in);
    
}

private static File[] getOutFiles( File user, List<String> types ){
    File [] out = new File[types.size()];
    for(int i=0; i<out.length; i++){
        out[i] = new File(user, types.get(i).replaceAll("\\s+", "_"));
    }
    return out;
}
}
