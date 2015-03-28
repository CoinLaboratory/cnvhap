package lc1.util;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;



public class ConvertFstDirectory {

	public static void main(String[] args){
		try{
			String[] str = args[0].split(":");
			for(int i=0; i<str.length; i++){
			 File dir = new File(System.getProperties().getProperty("user.dir")+"/chr"+str[i]);
			
			ConvertFstDirectory cfd = new ConvertFstDirectory(dir);
			cfd.run();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	

	static class TxtFilter  implements FileFilter{
		final String txt;
		TxtFilter(String txt) {
			this.txt = txt;
		}
		
		@Override
		public boolean accept(File arg0) {
			// TODO Auto-generated method stub
		return arg0.getName().endsWith(txt);
		}

		
	};
	
	CompressDir compress;
	OutputStreamWriter snps;
	int index=-1;
	BufferedReader[] br;
	//BufferedReader pr;
	double[] st;
	
	List<int[]> indicesToComp = new ArrayList<int[]>();
	List<String> nme = new ArrayList<String>();
	double[] vals;
	String chr;
	public ConvertFstDirectory(File dir) throws Exception{
		File[] f = dir.listFiles(new TxtFilter(".fst"));
		//File probeloc = dir.listFiles(new TxtFilter(".probepos"))[0];
		//pr = new BufferedReader(new FileReader(probeloc));
		//pr.readLine();
		chr = dir.getName().substring(3);
		compress = new CompressDir(new File(dir, chr));
		this.br = new BufferedReader[f.length];
		st = new double[f.length];
		for(int k=0; k<f.length; k++){
			br[k] = new BufferedReader(new FileReader(f[k]));
			br[k].readLine();
		}
		OutputStreamWriter osw = compress.getWriter("Samples", true);
		for(int j=0;j<f.length; j++){
			String n1 = f[j].getName().split("\\.")[1];
			Set<String> nme1 = new HashSet<String>(Arrays.asList(n1.split("_")));
			for(int k=0; k<j; k++){
				String n2 = f[k].getName().split("\\.")[1];
				Set<String> nme2 = new HashSet<String>(Arrays.asList(n2.split("_")));
				nme2.retainAll(nme1);
				if(nme2.size()>0){
					boolean swp = false;
				
					if(n2.indexOf("ceu")>=0){
						swp = true;
					}
					String nam = n1+"v"+n2;
					if(swp){
						indicesToComp.add(new int[] {k,j});
						nam = n2+"v"+n1;
					}
					else{
						indicesToComp.add(new int[] {j,k});
					}
					
					this.nme.add(nam);
					osw.write(nam+"\n");
				}
			}
		}
		compress.closeWriter(osw);
		OutputStreamWriter osw1 = compress.getWriter("Name", true);
		osw1.write("left\tright\tLog R Ratio\n");
		osw1.write("chr\tstart\tend\tid\n");
		osw1.write("id\n");
		compress.closeWriter(osw1);
		snps = compress.getWriter("SNPS", false);
		vals = new double[nme.size()];
	}
	private void run() throws Exception {
		while(readLine()){}
	    compress.closeWriter(this.snps);
	    compress.run();
	}
	String block_id = "";
	public boolean readLine() throws Exception{
	
		for(int k=0; k<st.length; k++){
			String sti = br[k].readLine();
			if(sti==null) return false;
			String[] str =sti.split("\\s+"); 
			st[k] = Double.parseDouble(str[1]);
			if(k==0) {
				block_id = str[0];
			}
			else if(!str[0].equals(block_id)){
				throw new RuntimeException("block ids dont match "+str[3]+" "+block_id);
			}
		}
	//	String[] prob =this.pr.readLine().split("\t");
		index++;
		//if(!prob[3].equals(block_id)){
		//	throw new RuntimeException("line mismatch with probe file "+block_id+" "+prob[3]);
		//}
		String rsid = "A_"+block_id;
		OutputStreamWriter osw1 = compress.getWriter(rsid, true);
		snps.write("chr"+chr+"\t"+block_id+"\t"+block_id+"\tA_"+block_id+"\n");//prob[1]+"\t"+prob[2]+"\t"+rsid+"\n");
		for(int k=0; k<this.vals.length; k++){
			int[] pos = indicesToComp.get(k);
			//double sum = st[pos[0]]+st[pos[1]];
			//double diff = st[pos[0]]-st[pos[1]];
			double left = st[pos[0]];
			double right = st[pos[1]];
			osw1.write(String.format("%5.3g",left).trim() +"\t");
			osw1.write(String.format("%5.3g",right).trim() +"\t");
			osw1.write(String.format("%5.3g",(left-right)/(left+right)).trim() +"\n");
		}
		compress.closeWriter(osw1);
		return true;
	}
	
}
