import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.OutputStream;
import java.net.URL;

import sun.net.www.protocol.ftp.FtpURLConnection;


public class Download {
	
	public static void main(String[] args){
		try{
		downloadAll();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}

	public static void downloadAll() throws Exception {
		File f = new File("infile.txt");
		final String[] chrom = "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:21:22:23".split(":");
		
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
	}
	public static void download(String coord) throws Exception{
	URL url = new URL("ftp://ftp.sanger.ac.uk/pub4/humgen/cnv/chr"+coord+".zip");
	 FtpURLConnection con = (FtpURLConnection)url.openConnection();
	 File res = new File(coord+".zip");
	 if(res.exists()) return;
    // con.setRequestProperty("User-agent",Constants.USER_AGENT);
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
         con.close();
         
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
