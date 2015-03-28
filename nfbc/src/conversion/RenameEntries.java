package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class RenameEntries {
static CompressDir cd;
static ZipFile zf;
	
	public static void main(String[] args){
		//if(true) System.exit(0);
		try{
			String build = args[1];
			
			boolean changeNames = true;
			File f = new File(args[0]);
			String chr = f.getName().split("\\.")[0];
			f.mkdir();
			 cd = new CompressDir(new File(build+"/"+chr));
			 zf = new ZipFile(f);
			 transfer("Samples", "Samples", null); transfer("Name", "Name", null);
			BufferedReader br =//getBR("SNPS");
				new BufferedReader(new FileReader(new File(build+"_"+chr+".txt")));
			OutputStreamWriter osw = cd.getWriter("SNPS", false);
			String st = "";
			while((st = br.readLine())!=null){
				String[] str = st.split("\t");
				String nme = str[3];
				String nme2 = str[0].substring(3)+"_"+str[1];
			//	boolean sw = str[4].equals(".");
			//	String nme2 = sw ? str[3] : str[4];
				
				if(changeNames){
					transfer(nme,nme2,null);
					str[3] = nme2;
				}
				else transfer(str[3],str[3], str[5]+str[6]);
				
				osw.write(str[0]+"\t");
				for(int k=1; k<str.length; k++){
					int k1 = k;
//					if(sw && (k1==3 || k1==4)) k1 = k==3 ? 4 : 3;
					osw.write(str[k1]+"\t");
				}
				osw.write("\n");
			}
			cd.closeWriter(osw);
			cd.run();
			cd.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	public static BufferedReader getBR(String nme) throws Exception{
		return getBR(zf.getEntry(nme));	
	}
	public static BufferedReader getBR(ZipEntry nme) throws Exception{
		return new BufferedReader(
				new InputStreamReader(zf.getInputStream((nme))));	
	}
	public static void transfer(String nme, String nme2, String comment) throws Exception{
		ZipEntry z1 = zf.getEntry(nme);
		if(comment==null) comment = z1.getComment();
		OutputStreamWriter osw = cd.getWriter(nme2, true, comment);
		BufferedReader br = getBR(z1);
		String st = "";
		while((st = br.readLine())!=null){
			osw.write(st);osw.write("\n");
		}
		cd.closeWriter(osw);
	}
}
