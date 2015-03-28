package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.List;

import lc1.util.ApacheCompressor;
import lc1.util.Constants;

import org.apache.commons.compress.archivers.zip.ZipFile;

public class Repackage {
	
	public static void main(String[] args){
		try{
			String[] chrom = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:X:Y".split(":");
			//String[] chrom = "18".split(":");
			File  dirF = new File(System.getProperty("user.dir"));
			for(int i=0; i<chrom.length; i++){
				
				Repackage rp = new Repackage(dirF, chrom[i]);
				rp.run();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
   ZipFile zf1,zf2;
   int cent1, cent2;
   CompressDir dir1, dir2;
   BufferedReader snps1, snps2;
   OutputStreamWriter snpso1, snpso2;
   
   static String[] toExtract = new String[] {"Allele1", "Allele2",  "Allele Freq","Log R"};
  static String[] sep = new String[] {"", "\t", "\t", ""};
   
 int[] toext = null;
   String p,q;
   
   Repackage(File dir, String chr) throws Exception{
	   Constants.build = new String[] {"build36"};
	   File dirN = new File(dir, "new");
	   dirN.mkdir();
	   Constants.mid = new String[][] {{chr}};
	   cent1 = lc1.util.Constants.getKaryoFile(dir);
	   cent2 = lc1.util.Constants.getKaryoFile(dir.getParentFile());
	   p = chr+p;
	   q = chr+q;
	   dir1 = new CompressDir(new File(dirN, chr+"p"));
	   dir2 = new CompressDir(new File(dirN, chr+"q"));
	   zf1 = new ZipFile(new File(chr+"p.zip"));
	   zf2 = new ZipFile(new File(chr+"q.zip"));
	   snps1 = ApacheCompressor.readZip(zf1,"SNPS");
	   snps2 = ApacheCompressor.readZip(zf2,"SNPS");
	   snpso1 = dir1.getWriter("SNPS",false);
	   snpso2 = dir2.getWriter("SNPS",false);
	  
	   
	  if(toExtract!=null){
	   toext = new int[toExtract.length];
	  List<String> l =  Arrays.asList(ApacheCompressor.getIndiv(zf1, "Name").get(0).split("\t"));
	   for(int i=0; i<toExtract.length; i++){
		   toext[i] = find(l,toExtract[i]);
	   }
	  }
	  transfer("Name", 0, 0, toext, 0); transfer("Samples", 0,0);
	   transfer("Name", 1, 1, toext, 0); transfer("Samples", 1,1);
   }
   
   private int find(List<String> l, String string) {
	for(int i=0  ; i<l.size(); i++){
		if(l.get(i).indexOf(string)>=0) return i;
	}
	 return -1;
}

public void run() throws Exception{
	   run(0); run(1);
	   dir1.closeWriter(snpso1);
	   dir2.closeWriter(snpso2);
	   dir1.run();
	   dir2.run();
	   zf1.close();
	   zf2.close();
   }
   
   public void run(int from) throws Exception{
	   String str = "";
	   BufferedReader snps = from==0 ? snps1 : snps2;
	  for(int i=0; (str = snps.readLine())!=null  ; i++){
		   String[]st = str.split("\\s+");
		   int pos = Integer.parseInt(st[1]);
		   String id = (st[3]);
		   if(pos < cent2){
			   transfer(id, from,0, toext, null);
			   snpso1.write(str.replaceAll(q, p)); snpso1.write("\n");
		   }else{
			   transfer(id,from,1, toext, null);
			   snpso2.write(str.replaceAll(p,q)); snpso2.write("\n");
		   }
	   }
	   snps.close();
   }
   
   void transfer(String st, int from, int to) throws Exception{
		   ZipFile zf = from==0 ? zf1 : zf2;
		   CompressDir dir = to==0 ? dir1 : dir2;
		   OutputStreamWriter ow = dir.getWriter(st, true);
		   BufferedReader br = ApacheCompressor.readZip(zf,st);
		   String str = "";
		   for(int i=0; (str = br.readLine())!=null; i++){
			   ow.write(str); ow.write("\n");
		   }
		   dir.closeWriter(ow);
		   br.close();
	
   }
   
   void transfer(String st, int from, int to, int[] toExtract, Integer conditional) throws Exception{
	  
	   ZipFile zf = from==0 ? zf1 : zf2;
	   CompressDir dir = to==0 ? dir1 : dir2;
	   OutputStreamWriter ow = dir.getWriter(st, true);
	   BufferedReader br = ApacheCompressor.readZip(zf,st);
	   String str = "";
	   for(int i=0; (str = br.readLine())!=null; i++){
		   if(toExtract==null || (conditional!=null && conditional.intValue()!=i)){
			   ow.write(str); 
		   }else{
			   String[] st1 = str.split("\t");
			   int len = st1.length;
			   for(int i1=0; i1<toExtract.length; i1++){
				   int i2 = toExtract[i1];
				   ow.write(i2 < len ? st1[i2] :"NaN");  
				   ow.write(sep[i1]);  
			   }
		   }
		 ow.write("\n");
	   }
	   dir.closeWriter(ow);
	   br.close();

}
   
}
