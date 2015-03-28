package data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.rosuda.JRI.Rengine;

import cnvtrans.LoopCall;


public class RunR {

	BufferedReader r;
	static Rengine re = new Rengine(new String[] { "--vanilla" }, true,
			new LoopCall());
			
	int cnt =0;
	
	public int getSize(String first, File f, int k) throws Exception{
		ZipFile zf = new ZipFile(f);
		BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
		String st = "";
		int len=0;
		while((st=snps.readLine())!=null){
			len++;
		}
		snps.close();
/*		int len = zf.size();*/
		if(k==0) indiv =  read(zf, "Samples",0);
		
		
		//len=len-2;
		/*if(zf.getEntry("SNPS")!=null) len --; */
		zf.close();
		return len;
	}
	
	private int[] getAlias(List<String> indiv2, List<String> target) {
		int[] res = new int[indiv2.size()];
		for(int k=0; k<indiv2.size(); k++){
			res[k] = target.indexOf(indiv2.get(k));
		}
		return res;
	}

	//int[] alias;
	public void run() throws Exception{
		String st="";
		while((st=this.r.readLine())!=null){
			if(st.length()>0 && ! st.startsWith("#"));
			re.eval(st);
		}
		r.close();
		re.end();
		
	}
	
	public static void main(String[] args){
		try{
			File dir = new File(args[0]);
			File rscript = new File(args[1]);
			String[] chroms = args[2].split(":");
			String logr_id =args[3].replace("_", " ");
			RunR rr = new RunR(dir, chroms, logr_id, rscript,Integer.parseInt(args[4]));
			rr.run();
		}catch(Exception exc){
			exc.printStackTrace();
			System.exit(0);
		}
		System.exit(0);
	}
	public RunR(File dir, String[] chr, String logr_id, File rscript, int thin) throws Exception{
		List<File> l = new ArrayList<File>();
		for(int k=0; k<chr.length; k++){
			File f1 = new File(dir, chr[k]+".zip");
			if(f1.exists()){
				l.add(f1);
			}
		}
		this.r = new BufferedReader(new FileReader(rscript));
		
		this.st1 = readNext();
		this.st2 = readNext();
		this.st3 = readNext();
		this.st4 = readNext();
		System.err.println("first 4 lines are ");
		System.err.println(st1);
		System.err.println(st2);
		System.err.println(st3);
		System.err.println(st4);
		//this.st1 = "mat = matrix(0,ncol=sze,nrow = indvs)";
		//this.st2 = "mat[,cnt] = tmp";
		this.thin = thin;
		File[] fi = l.toArray(new File[0]);
		//this.alias = new int[fi.length][];
		this.setup(fi);
		for(int k=0; k<fi.length; k++){
			this.readData(fi[k], logr_id,k);
		}
		
		
	}
	private String readNext() throws Exception{
		String st = r.readLine().trim();;
		while(st.startsWith("#") && st.length()==0){
			st = r.readLine().trim();;
		}
		return st;
	}
	final int thin;
	public void setup(File[] f) throws Exception{
		int sze = 0;
		for(int k=0; k<f.length; k++){
			sze+=Math.ceil(this.getSize(f[0].getName(), f[k],k)/(double) thin);
		}
		re.assign("indiv", this.indiv.toArray(new String[0]));
		re.eval(st1.replaceFirst("sze", sze+"").replace("indvs",indiv.size()+""));
		re.eval(st3.replaceFirst("sze", sze+"").replace("indvs",3+""));
		//re.eval("print(mat)");
	}
	final String st1,st2,st3,st4;
	List<String> indiv;
	public void readData(File f, String id, int index) throws Exception{
		ZipFile zf = new ZipFile(f);
		int[] alias = null;
		if(index>0){
			List<String> indiv1 =  read(zf, "Samples",0);
			 if(!indiv.equals(indiv1)) {
				if(indiv1.size()!=indiv.size()) throw new RuntimeException("size has to be same");
				System.err.println("warning sample order diff betwwen "+f.getName());
				alias = getAlias(indiv1, indiv);
			//	throw new RuntimeException();
			}
		}
		List<String > l2 =  read(zf, "Name");
		List<String> l20 =Arrays.asList(l2.get(0).split("\t")); 
		List<String> l21 =Arrays.asList(l2.get(1).split("\t")); 
		fix(l20);fix(l21);
		int ind = l20.indexOf(id);
		int snp_ind = l21.indexOf("id");
	//	List<String>
		double[] res1 = new double[indiv.size()];
		int k=0;
		BufferedReader snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry("SNPS"))));
		String st = "";
		String[] names = new String[3];
		while((st = snps.readLine())!=null){
			if(k==thin) k=0;
			String[] snps_ =st.split("\\s+"); 
			//ZipEntry ent = (ZipEntry)en.nextElement();
			String name = snps_[snp_ind];
			names[0] = name;
			names[1] = snps_[0];
			names[2] = snps_[1];
			//if(name.equals("SNPS") || name.equals("Samples") || name.equals("Name")) continue;
			if(k==0){
				this.read(zf, name, ind, res1, alias);
				re.assign("tmp", res1);
				re.assign("tmp1", names);
				String st_new = st2.replaceAll("cnt", (cnt+1)+"");
			//	System.out.println(st_new);
				re.eval(st_new);
				if(st4.length()>0){
					st_new = st4.replaceAll("cnt", (cnt+1)+"");
				//	System.out.println(st_new);
					re.eval(st_new);	
				}
				cnt++;
				
			}
			k++;
		}
		snps.close();
		zf.close();
		
	}
	
	public static void fix(List<String> l21) {
		// TODO Auto-generated method stub
		for(int k=0; k<l21.size(); k++){
			l21.set(k, l21.get(k).trim());
		}
	}

	private List<String> read(ZipFile zf,  String probe, int i) throws Exception {
		List<String> res = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st.split("\t")[0]);
		}
		br.close();
		return res;
		
		
	}
	public static List<String> read(ZipFile zf,  String probe) throws Exception {
		List<String> res = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st);
		}
		br.close();
		return res;
		
		
	}
	private double[] read(ZipFile zf,  String probe, int i, double[] res, int[] alias) throws Exception {
		//List<String> res = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		for(int k=0; (st =br.readLine() )!=null;k++){
			res[alias==null ? k : alias[k]]=Double.parseDouble(st.split("\t")[i]);
		}
		br.close();
		return res;
		
	}
	
}
