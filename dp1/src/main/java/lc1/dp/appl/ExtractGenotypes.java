package lc1.dp.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import lc1.dp.data.collection.LightWeightDataCollection;
import lc1.util.Constants;
import lc1.util.Executor;

public class ExtractGenotypes {
	static class FileFilter1 implements FileFilter {
		final Pattern st;
		FileFilter1(Pattern p){
			this.st =p;
		}
		public boolean accept(File pathname) {
			Matcher mat = st.matcher(pathname.getName());
		
			mat.useAnchoringBounds(false);
			return mat.find();
			
			}
	}

	public static void main(String[] args){
		//Constants.topBottom = false;
		Constants.loess = new boolean[]{false};
		Constants.median_correction = new boolean[] {false};
		Constants.writeAverages = new String[] {"countAll", "countA"};
		Constants.var_thresh = new double[] {Double.POSITIVE_INFINITY};
		
		try{
		
		File f = new File(args[0]);
		/*		List<File> f1 = listFilesRecursively(f, Arrays.asList(args[1].split("/")));
			String[] arg1 = new String[f1.size()];
			for(int i=0; i<arg1.length; i++){
				arg1[i] = f1.get(i).getAbsolutePath();
			}*/
			run(f,args);
			es.shutdown();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	static Map<String, Integer> readInfile() throws Exception{
		Map<String, Integer> res = new HashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(new File("infile.txt")));
		String st = br.readLine();
		while((st=br.readLine())!=null){
			String[] str = st.split("\\s+");
			res.put(str[0], Integer.parseInt(str[1]));
		}
		return res;
	}
	
	 static ExecutorService es = Executor.getEs(ExtractGenotypes.class, 4);
	private static List<File> listFilesRecursively(File f, List<String> string) {
			if(!f.isDirectory()) {
				throw new RuntimeException("!!");
			}
			final String st = string.get(0);
			Pattern p = Pattern.compile(st);
		//	Pattern.compile(st).;
			List<File> res = new ArrayList<File>();
			boolean last = string.size()==1;
		    File[] fs1 = f.listFiles(new FileFilter1(p));
			for(int i=0; i<fs1.length; i++){
				if(last)res.add(fs1[i]);
				
				else if(fs1[i].isDirectory()) res.addAll(listFilesRecursively(fs1[i], string.subList(1, string.size())));
			}
			return res;
		
	}
	public static void run(final File dirF, String[] args) throws Exception{
	//	Map<String, Integer> len = (new File("infile.txt")).exists() ? ExtractGenotypes.readInfile() : null;
		//for(int ii=0; ii<args.length; ii++){
			List l = new ArrayList();
		//	final String arg1 = args[ii];
			
		try{
		//	final File dirF =new File(arg1);
			final File bf = new File(dirF, "build36.txt");
		//		Constants.modelCNP = 10;
			File outdir = new File(dirF, "genotypes");
			outdir.mkdir();
			//PrintWriter failed = new PrintWriter(new BufferedWriter(new FileWriter(new File(outdir, "failed.txt"))));
			File[] zip = dirF.listFiles(Constants.ZIP_FILTER);
			String chr = null;
			final int[] startend ;
			final String[] type;
			
			if(args.length>1){
				String[] args1 = args[1].split(":");
				chr = args1[0];
				startend = new int[] {Constants.convert(args1[1]), Constants.convert(args1[2])};
				if(args.length>2){
					type = args[2].split(":");
				}
				else{
					type =   new String[] { "countAll", "countA"};
				}
			}
			else{
				type =    new String[] { "countAll", "countA"};
				startend=null;
			}
			
			inner: for(int i=0; i<zip.length; i++){
				final String chrom = zip[i].getName();
				if(chr!=null && !chrom.equals(chr+".zip")){
					System.err.println("skip "+chrom+" "+chr);
					continue inner;
				}
				//if(!chrom.equals("chr11")) continue;
				Callable call = new Callable(){
					public Object call(){
				
				try{
				ExtractGenotypes ec = new ExtractGenotypes(dirF, chrom,  bf);
					try{
					//	if()
						
					ec.run(startend, type);
					}catch(Exception exc){
				
						exc.printStackTrace();
					}
				
				}catch(Exception exc){
					exc.printStackTrace();
				
				}
				return null;
					}
					};
					l.add(call);
				}
			try{
				es.invokeAll(l);
				}catch(Exception exc){
					exc.printStackTrace();
				}	
					
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
			
			
	//	}
		
	}
		LightWeightDataCollection dc;
		//ZipFile zf;
		
		List<String>  header_geno, header_sample, indiv;
	//	File bf;
		File parentDir;
		final String chrom;
		final int start,end;
		
		int[] mid = new int[] {0,Integer.MAX_VALUE-1};
		ExtractGenotypes(File f,   String chrom,  File bf) throws Exception{
			System.err.println("opening "+f);
			this.chrom = chrom.split("\\.")[0];
			this.start = mid[0];
			this.end = mid[1];
			 //  zf = new ZipFile(f);
			//   this.bf = bf;
			 // DataCollection[] dc1 = new DataCollection[dirs.size()];
			
				 dc = new LightWeightDataCollection(new File(f, chrom),(short)0, 2, new int[][] {mid},bf,0,null);
			
			   this.parentDir = f;
			
		       //  
		}
		
		public void run(int[] startend, String[] outp) throws Exception{
			File outdir = new File(parentDir, "genotypes");
			outdir.mkdir();
				PrintWriter pw_hap2 = new PrintWriter(new BufferedWriter(new FileWriter(new File(outdir, chrom+".txt"))));
			// dc.writeFastphase(pw_hap2, false);
				if(startend!=null){
			dc.printHapMapFormat( 
                     pw_hap2, dc.indiv(), null, true, 
                     new String[] {"snpid", "loc"}, 
                     new String[0],
                    outp, "%7s",startend[0], startend[1]);
				}
				else{
					dc.printHapMapFormat( 
		                     pw_hap2, dc.indiv(), null, true, 
		                     new String[] {"snpid", "loc"}, 
		                     new String[0],
		                      outp, "%7s");
				}
              //pw_hap1.close();
             pw_hap2.close();
          //   dc.writeSNPFile(new File("snp_"+rs+".txt"), Constants.chrom0(), false, null);
			
		}
		
	
		
		
		


		private String getMb(int start2) {
			return Math.round((double)start2 / (1000.0*1000.0))+"";
		}

		private boolean iscn(String string) {
			return string.length()!=2 || string.indexOf('_')>=0 || string.indexOf('X')>=0 || string.indexOf('Y')>=0|| string.indexOf('Z')>=0; 
		}
}
