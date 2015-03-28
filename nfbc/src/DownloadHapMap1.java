import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import conversion.ConvertIllumina;



public class DownloadHapMap1 {
	
	
	
	public static void main(String args[]){
		try{
		download("NA06994");
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public static void main1(String[] args){
	//	ConvertIllumina.delete = false;
			File output = new File("output");
			output.mkdir();
			File data1 = new File(output,"data1");
			data1.mkdir();
			
			
			
							List<String> chr1 = //Arrays.asList("22".split(":"));
								Arrays.asList("1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22".split(":"));
			Collections.reverse(chr1);
			final String[] chrom = chr1.toArray(new String[0]);

			final String[] pop ="yri"
//				"asw:chb:chd:gih:jpt:lwk:mex:mkk:tsi:yri"
				//"mkk:mex:tsi:asw:chd:yri:jpt+chb:gih"// "CEU:YRI:"
				//"tsi:lwk:mkk:asw:jpt+chb:mex:gih:chd"
			.split(":");
			   

				
				for(int ik=0; ik<pop.length; ik++){
					try{
			String[] type =new String[]{""};// new String[] {"unr.", "D.",""};
			String[] type1 =new String[] {""};// new String[] {"UNRELATED", "DUOS", "TRIOS"};
			File outd = new File(data1,pop[ik]);
			outd.mkdir();
		    for(int i=0; i<chrom.length; i++){
		    //	Set<String> outb_s = new H
		    	StringBuffer outb = new StringBuffer();
		    	File[] outf = new File[type.length];
		    	boolean first = true;
		    	for(int j=0; j<type.length; j++){
		    		 outf[j] = readPage(chrom[i], pop[ik], type[j], type1[j]);
		    		 if(outf[j]!=null){
		    			 if(!first){
		    				 outb.append(":");
		    			 }
		    			 else first = false;
		    			 outb.append(outf[j].getAbsolutePath());
		    		 }
		    		
		    		
		    	}
		    	File resf = new File(data1,chrom[i]+".zip");
	    		if(!resf.exists()){
	    		ConvertIllumina.main(
	    				("--build build36  --splStr \\s --mode compress  " +
	    				"--file  "+outb.toString()+" --numThreads 4 --output "+outd.getAbsolutePath()+" --format wide --default_type Genotype").split("\\s+")
	    				);
	    		File outf1 = new File(outd, chrom[i]+".zip");
	    		File targetDir = outf[0].getParentFile().getParentFile();
	    		System.err.println("target dir is "+targetDir);
	    		outf1.renameTo(new File(targetDir, chrom[i]+".zip"));
	    		  for(int j=0; j<outf.length; j++){
	    			   if(outf[j]!=null){
	    				   outf[j].delete();
	    			   }
	    		   }
	    		        //if(true) System.exit(0);
	    		}
	    		
	    		
		    }
				
		}catch(Exception exc){
			exc.printStackTrace();
		}
				}
	}
	
	
	
	public static File readPage(String chr, String pop, String type, String type1) throws Exception{
		
		String name =
			
			"genotypes_chr"+chr+"_"+pop.toUpperCase()+"_r27_nr.b36_fwd.txt.gz";
			//"hapmap3_r2_b36_fwd.consensus.qc.poly.chr"+chr+"_"+pop+"."+type+"phased.gz";
		String url = 
			"http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes/latest_phaseII+III_ncbi_b36/forward/non-redundant/"+name;
			//"http://hapmap.ncbi.nlm.nih.gov/downloads/genotypes//forward/non-redundant/";
//			"http://ftp.hapmap.org/genotypes//forward/non-redundant/"+name;
			
		//	"http://ftp.hapmap.org/phasing/2009-02_phaseIII/HapMap3_r2/"+pop.toUpperCase()+"/"+type1+"/"+name;
//		String url = "http://ftp.hapmap.org/phasing/2009-02_phaseIII/HapMap3_r2/YRI/UNRELATED/";
		File popf = new File(pop);
		popf.mkdir();
		File chrf = new File(popf, chr+"");
		chrf.mkdir();
		File outf = new File(chrf, name);
		if(outf.exists()) return outf;
		System.err.println(url);
		boolean success = download1(url, outf);
		if(!success){
			 name = "hapmap3_r2_b36_fwd.consensus.qc.poly.chr"+chr+"_"+pop+"."+type+"phased.gz";
			 url = "http://ftp.hapmap.org/phasing/2009-02_phaseIII/HapMap3_r2/"+pop.toUpperCase()+"/"+name;
			 success = download1(url, outf);
		}
		if(success) return outf;
		else return null;
	}

	/*public static void downloadAll() throws Exception {
		File f = new File("infile.txt");
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String st = br.readLine();
		
		final int[] len = new int[chrom.length];
		for(int i=0; (st = br.readLine())!=null; i++){
			String[] str = st.split("\\s+");
			len[i] = (Integer.parseInt(str[1]));
		}
		Thread th = new Thread(new Runnable(){
			
		public void run(){
		
		for(int i=0; i<chrom.length; i++){
			double max =(int)  Math.ceil((double)len[i]/(1000.0*1000.0));
			System.err.println(chrom[i]+" "+max);
			for(int j=0;  j<max; j+=5){
//				10_5000001_10000001
				String coord;
				if(j+5 >= max){
					coord =  chrom[i]+"_"+(j*1000*1000+1)+"_"+(len[i]);
				}
				else coord= chrom[i]+"_"+(j*1000*1000+1)+"_"+((j+5)*1000*1000+1);
				System.err.println(coord);
				try{
				download(coord);
				//if(true) return;
 				Thread.sleep(1000);
				}catch(Exception exc){
					exc.printStackTrace();
					System.exit(0);
				}
			}
		}
		}
		});
		th.start();
	}*/
	public static void download(String id) throws Exception{
		
		String name = id+".SLX.maq.SRP000031.2009_08.bam.bas";
		download1("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/data/"+id+"/alignment/"+name, new File(name));
		//download1("ftp://ftp.sanger.ac.uk/pub4/humgen/cnv/chr"+coord+".zip", new File(coord+".zip"));
	}
	
	public static boolean download1(String urlst, File res) {
		try{
			System.err.println(urlst);
	URL url = new URL(urlst);
	URLConnection con = url.openConnection();
	 
	 if(res.exists() && res.length()>0) return true;
     con.setRequestProperty("User-agent","java");
     con.connect();
     
     /*int response = con.getR
     if ((response != HttpURLConnection.HTTP_ACCEPTED) && (response != HttpURLConnection.HTTP_OK)) {
         //if something went wrong
         throw new IOException("Could not connect to update server.");
     }
     else {*/
   //  con.getInputStream()
    	 BufferedInputStream inputStream = new BufferedInputStream(con.getInputStream());
    	 OutputStream os = new BufferedOutputStream(new FileOutputStream(res));
         byte[] buf = new byte[200];
         int size = con.getContentLength();
         int read=1;
         while(read!=-1){
        	 
         if(size > 200) {
             int len1 = read = inputStream.read(buf,0,200);
            if(read!=-1) os.write(buf,0,len1);
         }
         else {
        	 int len1 = read = inputStream.read(buf,0,size);
            if(read!=-1) os.write(buf,0,len1);
         }
       
       
         }
         os.close();
         inputStream.close();
         if(con instanceof HttpURLConnection){
     ((HttpURLConnection)con).disconnect();
         }
     return res.exists();
		}catch(Exception exc){
			System.err.println("WARNING : "+exc.getMessage());
			exc.printStackTrace();
			return false;
		}
         
   //  }
   // pw.close();
}
	/* public boolean checkForUpdate() throws IOException{

	        try {
	            URL url = new URL("http://www.broad.mit.edu/mpg/haploview/uc/version.txt");
	            HttpURLConnection con = (HttpURLConnection)url.openConnection();
	            con.setRequestProperty("User-agent",Constants.USER_AGENT);
	            con.connect();

	            int response = con.getResponseCode();

	            if ((response != HttpURLConnection.HTTP_ACCEPTED) && (response != HttpURLConnection.HTTP_OK)) {
	                //if something went wrong
	                throw new IOException("Could not connect to update server.");
	            }
	            else {
	                //all is well
	               
	                String data = "";
	                if(read != 0)  {
	                    data = new String(buf);
	                    double newestVersion = Double.parseDouble(data);

	                    if(newestVersion > Constants.VERSION) {
	                        this.newVersion = newestVersion;
	                        this.newVersionAvailable = true;
	                    }
	                    else {
	                        this.newVersionAvailable = false;
	                        this.newVersion = Constants.VERSION;
	                    }

	                }
	            }
	            con.disconnect();

	        } catch(MalformedURLException mue) {
	            //System.err.println("the following url exception occured:" + mue);
	        }*/
}
