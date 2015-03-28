package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipFile;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.util.Compressor;
import lc1.util.Constants;

public class ExtractHaplotypes {
	
	public static void main(String[]args){
		try{
			
			//String chrom = "11";
		
			File dirF =new File( System.getProperty("user.dir"));
			File in_file = new File(dirF, args[0]);
			BufferedReader br = new BufferedReader(new FileReader(in_file));
			String st = br.readLine();
			Map<String, List<Object[]>> map  = new HashMap<String, List<Object[]>>();
			while((st = br.readLine())!=null){
				String[] str = st.trim().split("\\s+");
				String chrom = str[2];
				Integer loc = Integer.parseInt(str[1]);
				String rs = str[0];
				//if(rs.length()==0  || chrom.indexOf('.')>=0) continue;
				List<Object[]> l = map.get(chrom);
				if(l==null){
					map.put(chrom, l = new ArrayList<Object[]>());
				}
				System.err.println(chrom+" "+rs);
				l.add(new Object[] {rs, loc});
			}
			//File bf = new File(dirF, "build36.txt");
			//	BufferedReader buildF =    DataCollection.getBufferedReader(bf);
        
		//	Constants.modelCNP = 10;
			PrintWriter failed = new PrintWriter(new BufferedWriter(new FileWriter("failed.txt")));
			for(Iterator<String> it = map.keySet().iterator(); it.hasNext();){
				String chrom = it.next();
				List<Object[]> l = map.get(chrom);
				try{
				ExtractHaplotypes ec = new ExtractHaplotypes(new File(dirF, chrom+".zip"), chrom);
				for(int i=0; i<l.size(); i++){
					try{
					System.err.println("DOING "+chrom+" "+l.get(i)[0]);
					ec.run((String)l.get(i)[0], (Integer) l.get(i)[1]);
					
					}catch(Exception exc){
						exc.printStackTrace();
						failed.println(chrom+" "+l.get(i)[0]+" "+l.get(i)[1]+" "+exc.getMessage());
					//	exc.printStackTrace();
					}
				}
				ec.zf.close();
				}catch(Exception exc){
					exc.printStackTrace();
					for(int i=0; i<l.size(); i++){
						failed.println(chrom+" "+l.get(i)[0]+" "+l.get(i)[1]+" "+exc.getMessage());
					}
					//ec.zf.close();
				}
				
				
			}
			failed.close();
			
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
		DataCollection dc;
		ZipFile zf;
		
		List<String>  header_geno, header_sample, indiv;
	//	File bf;
		File inputF;
		ExtractHaplotypes(File f,  String chrom) throws Exception{
			System.err.println("opening "+f);
			this.dir = f.getParentFile();
			   zf = new ZipFile(f);
			//   this.bf = bf;
			   this.inputF = f;
			  List<String> headers=  Compressor.getIndiv(zf, "Name", null);
			//  String[]   header = headers.get(0).split("\t");
		         header_geno = Arrays.asList(headers.get(0).split("\t"));
		         header_sample = Arrays.asList(headers.get(2).split("\t"));
		         for(int i=0; i<header_geno.size(); i++){
		        	 header_geno.set(i, header_geno.get(i).trim());
		         }
		         indiv = Compressor.getIndiv(zf, "Samples", header_sample.indexOf("id"));
		         this.readBuildFile(Compressor.getBufferedReader(zf, "SNPS"), chrom);
		
			
			
		       //  
		}
		
		final File dir;
		public void run(String rs, Integer pos) throws Exception{
			 dc = new SimpleDataCollection(inputF,(short)0, 2, new int[][] {this.getLocs(pos)},null,null);
			 dc.name = "hap";
			File pw_hap2 = new File(dir, rs);
			pw_hap2.mkdir();
			 dc.writeFastphase(pw_hap2, false);
              //pw_hap1.close();
          //   pw_hap2.close();
             dc.writeSNPFile(new File(pw_hap2,"snp.txt"), Constants.chrom0(), false, null);
			
		}
		
		static int chr_index = 0;
		static int pos_index = 1;
		static int rs_index = 3;
		static class Posi{
			int start;
			String chrom;
			String rs;
			Posi(String[] pos){
				this.chrom = pos[chr_index];
				this.start = Integer.parseInt(pos[pos_index]);
				this.rs = pos[rs_index];
			}
		}
		List<String> locs = new ArrayList<String>();
		List<Integer> pos = new ArrayList<Integer>();
		public  void readBuildFile(BufferedReader br,  String chrom1) throws Exception{
		    locs = new ArrayList<String>();
		    String st = "";
		  String chrom = "chr"+chrom1;
		   for(int i=0;(st = br.readLine())!=null; i++){
			  
		     String[] str = st.split("\\s+");  
		       if( str[chr_index].equals(chrom)){
		    	   locs.add(str[rs_index]);
		    	   pos.add(Integer.parseInt(str[pos_index]));
		        
		        }
		    }
		   br.close();
		   // return readZip(zf, l);
		    
		}
		
		
		public int[] getLocs(Integer location) throws Exception{
			int posi = pos.indexOf(location);
			if(posi<0){
				
				inner: for( posi=0; posi<pos.size(); posi++){
					
					if(pos.get(posi)>location){
						int dist = pos.get(posi) - location;
						if(posi>1 && location - pos.get(posi-1) < dist){
							posi = posi - 1;
						
						}
						break inner;
					}
				}
				if(posi==pos.size()) posi--;
			}
			Set<Integer> set1 = new HashSet<Integer>();
			int geno_ind = header_geno.indexOf("Genotype");
			List<String> l = Compressor.getIndiv(zf, locs.get(posi), header_geno.indexOf("Genotype"));
			for(int i=0; i<l.size(); i++){
				if(iscn(l.get(i))) set1.add(i);
			}
			List<Integer> positions = new ArrayList<Integer>();
			positions.add(pos.get(posi));
			//to right
			Set<Integer> set = new HashSet<Integer>(set1);
			for(int pos1 = posi+1; set.size()>0 && pos1 < locs.size(); pos1++){
				String rs1 = locs.get(pos1);
				positions.add(pos.get(pos1));
				l =  Compressor.getIndiv(zf, rs1, geno_ind);
				for(int i=0; i<l.size(); i++){
					if(!iscn(l.get(i))) set.remove(i);
				}
			}
			//to left
			set = new HashSet<Integer>(set1);
			for(int pos1 = posi-1; set.size()>0 && pos1 >=0; pos1--){
				String rs1 = locs.get(pos1);
				positions.add(0,pos.get(pos1));
				l =  Compressor.getIndiv(zf, rs1, geno_ind);
				for(int i=0; i<l.size(); i++){
					if(!iscn(l.get(i))) set.remove(i);
				}
			}
			return new int[]{positions.get(0)-10, positions.get(positions.size()-1)+10};
		}


		private boolean iscn(String string) {
			return string.length()!=2 || string.indexOf('_')>=0 || string.indexOf('X')>=0 || string.indexOf('Y')>=0|| string.indexOf('Z')>=0; 
		}
}
