package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.ZipFile;

public class MergeLongZips {

	//File[] zips;
	ZipFile current;
	BufferedReader snps;
	Map<String, File>  dirs = new HashMap<String, File>();
	
	static int chr_ind = 1;
	static int snp_ind = 0;
//	int step;
//	int noSnps;
	File dir;
	public static void main(String[] args){
		try{
		//	if(true)System.exit(0);
			File dir = new File(System.getProperty("user.dir"));
			snp_ind = Integer.parseInt(args[2]);
			chr_ind = Integer.parseInt(args[3]);
			System.err.println("set chr_ind "+chr_ind);
			System.err.println("set snp_ind "+snp_ind);
			MergeLongZips mlz = new MergeLongZips(args[0].split(":"),dir, Integer.parseInt(args[1]));
			mlz.run();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	int step;
	MergeLongZips(String[] chr, File dir, int step) throws Exception{
		this.dir = dir;
		compress = new CompressDir(dir);
		this.snps = Utils.getBufferedReader(new File(dir, "snps.txt"));
	/*	zips = dir.listFiles(new FileFilter(){
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith("zip");
			}
			
		});
		Arrays.sort(zips,new Comparator<File>(){
			@Override
			public int compare(File o1, File o2) {
				Integer i1 = Integer.parseInt(o1.getName().split("_")[1].split("\\.")[0]);
				Integer i2 = Integer.parseInt(o2.getName().split("_")[1].split("\\.")[0]);
				return i1.compareTo(i2);
			}
			
		});*/
	//	this.current = new ZipFile(zips[0]);
		//this.step=Integer.parseInt(zips[0].getName().split("_")[1]);
		this.step = step;
		for(int i=0; i<chr.length; i++){
			File out = new File(dir,chr[i]);
			out.mkdir();
			this.dirs.put(chr[i], out);
			
		}
		
	}
	
	
	
File currentF = null;
final CompressDir compress;

public void run() throws Exception{
		String st = "";
		boolean first = true;
		int rem=0;
		for(int i=0; (st = snps.readLine())!=null; i++){
			int j = (int) Math.floor((double) i/(double)step);
			if( rem==0){
				if(current!=null) this.current.close();
				int k = j*1000;
				File f = new File(dir, "out_"+k+".zip");
				if(f.exists()){
				current = new ZipFile(f);
				currentF = f;
				if(first){
					for(Iterator<File> it = this.dirs.values().iterator(); it.hasNext();){
						File out = it.next();
					try{
						write(out,"Name",current);
					}catch(Exception exc){
						System.err.println("prob with Name "+currentF);
						exc.printStackTrace();
					}
					try{
						write(out,"Samples",current);
					}catch(Exception exc){
						System.err.println("prob with Samples "+currentF);
						exc.printStackTrace();
					}
						
					}
					first = false;
				}
				}
				else{
					System.err.println("did not exist "+f.getAbsolutePath());
					current=null;
				}
				System.err.println("current "+f.getAbsolutePath());
			}
			if(current!=null){
				try{
				String[] str = st.split("\\s+");
				File dir = dirs.get(str[chr_ind]);
				//System.err.println("snp "+str[snp_ind]);
				if(dir!=null){
					
					write(dir, str[snp_ind], current);
					
				}
				}catch(Exception exc){
					
					System.err.println("prob with Name "+currentF+" "+st);
					exc.printStackTrace();
				}
			}
			rem++;
			if(rem==step) rem=0;
			
		}
		if(current!=null) this.current.close();

		compress.run();

		
	}
	
	
	void write(File dir, String nme,  ZipFile zf) throws Exception{
		//try{
		BufferedReader br = Utils.getBufferedReader(zf, nme);
		OutputStreamWriter osw = compress.getWriter(nme, true);
		//PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir,nme))));
		
		String st1 = "";
		while((st1=br.readLine())!=null){
			osw.write(st1);
			osw.write("\n");
		}
		compress.closeWriter(osw);
	/*	}catch(Exception exc){
			exc.printStackTrace();
			System.err.println("prob with "+nme);
			System.exit(0);
		}*/
	}
}
