package nfbc;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import net.sf.jftp.config.Settings;
import net.sf.jftp.net.BasicConnection;
import net.sf.jftp.net.ConnectionHandler;
import net.sf.jftp.net.ConnectionListener;
import net.sf.jftp.net.DataConnection;
import net.sf.jftp.system.StringUtils;

import com.jcraft.jsch.Logger;
import com.sshtools.j2ssh.configuration.SshConnectionProperties;

// this class download a file via anonymous ftp and shows output.
//
// if you want to use the api in a more complex way, please do at least take a look at the
// FtpConnection, FtpTransfer, ConnectionHandler, DirPanel (blockedTransfer, transfer)
// and ConnectionListener sourcecode.
public class FTPDownload implements Logger, ConnectionListener
{

 // is the connection established?
 private boolean isThere = false;

 public static long time = 0;

 // connection pool, not necessary but you should take a look at this class
 // if you want to use multiple event based ftp transfers.
 private ConnectionHandler handler = new ConnectionHandler();
 SftpConnection con;
 String user;
//creates a FtpConnection and downloads a file
 public FTPDownload(String host, String user, String pw)
 {
     this.user = user;
    // the ftp client default is very small, you may want to increase this
    Settings.bufferSize = 16384; 

    long current = System.currentTimeMillis();
    //System.out.println("1) "+(System.currentTimeMillis()-current)+"ms.");

    // register app as Logger, debug() and debugRaw below are from now on called by
    // FtpConnection
  //  Log.setLogger(this);

    // create a FtpConnection - note that it does *not* connect instantly
    SshConnectionProperties prop = new SshConnectionProperties();
    prop.setHost(host);
   con = new SftpConnection(prop);
  //  ScpClient con = new ScpClient();

    //System.out.println("2) "+(System.currentTimeMillis()-current)+"ms.");

    // set updatelistener, interface methods are below
    con.addConnectionListener(this);
    con.login(user, pw);
   // con.chdir("/home/lcoin/bcos");
    // set handler
    //con.setConnectionHandler(handler);

    // connect and login. this is from where connectionFailed() may be called for example
   // con.login("cdemon","........");

    //System.out.println("3) "+(System.currentTimeMillis()-current)+"ms.");

    // login calls connectionInitialized() below which sets isThere to true
   /* while(!isThere)
    {
        try { Thread.sleep(10); }
        catch(Exception ex) { ex.printStackTrace(); }
    }*/

    //System.out.println("4) "+(System.currentTimeMillis()-current)+"ms.");

    // download the file - this method blocks until the download has finished
    // if you want non-blocking, multithreaded io, just use
    //
    // con.handleDownload(file);
    //
    // which spawns a new thread for the download
   
 }
 public void disconnect(){
     con.disconnect();
 }
 private int whichIndex(String st, List<String> str){
     for(int i=0; i<str.size(); i++){
         if(str.get(i).indexOf(st)>=0) return i;
     }
     return -1;
 }
 public   void getFiles1(String base, String experiment, String type) throws Exception{
   
     con.chdir(base);
     File out = new File(experiment);
     out.mkdir();
     String[] f = con.sortLs();
     for(int i=0; i<f.length; i++){
         if(f[i].startsWith(experiment)){
             con.chdir(base+"/"+f[i]);
             String[] f1 = con.sortLs();
             for(int j=0; j<f1.length; j++){
                 if(f1[j].endsWith("zip")){
//                     InputStream is =con.getDownloadInputStream(f1[j]);
                     ZipInputStream in = new ZipInputStream(con.getDownloadInputStream(f1[j]));
                    BufferedReader br = new BufferedReader(new InputStreamReader(in));
                     while(true){
                            ZipEntry ent = in.getNextEntry();
                            if(ent==null) break;
                            if(ent.isDirectory() || ent.getName().indexOf(type)<0) {
                             //   System.err.println(ent.getName());
                                in.closeEntry();
                                continue;
                            }
                            String st = "";
                            while((st!=null)){
                               st =  br.readLine();
                                System.err.println(st);
                            }
                           /* byte[] buf = new byte[smbBuffer];
                            int len = 0;
                            int reallen = 0;
        
                            //System.out.println(file+":"+getLocalPath()+outfile);
                            while(true)
                            {
                                len = in.read(buf, 0, 100);
                                if(len == StreamTokenizer.TT_EOF)
                                {
                                    break;
                                }
        
                               // out.write(buf, 0, len);
                                reallen += len;
        
                                //System.out.println(file + ":" + StringUtils.getFile(file));
                               
                                {
                                    con.fireProgressUpdate(StringUtils.getFile(f1[j]),
                                                       DataConnection.GET, reallen);
                                }
                            }*/
                     }
                  // con.download(f1[j],out.getAbsolutePath()+"/");    
                 }
             }
        
         }
     }             
 }
 public static int smbBuffer = 32000;
 
 public   void getFiles2(String base, String experiment){
	 con.chdir(base+"/"+experiment);
     File user = new File(System.getProperty("user.dir"));
     File out = new File(user, experiment);
     out.mkdir();
     String[] f = con.sortLs();
     String[] size = con.size;
     for(int i=0; i<f.length; i++){
    	 if(f[i].startsWith(".")) continue;
    	 File localDir = new File(out, f[i]);
    	 
    	 if(!localDir.exists() || (!size[i].equals(localDir.length()+""))){
    		 System.err.println("downloading "+f[i]+" "+size[i]+" "+localDir.length());
              con.download(f[i], localDir.getParentFile().getAbsolutePath()+"/");
              
    	 }
    	 else{
    		 System.err.println("already exists "+f[i]+" "+size[i]+" "+localDir.length());
    	 }
     }
 }
 public   void getFiles(String experiment){
     String base = "/home/"+user+"/bcos";
     con.chdir(base);
     File out = new File(experiment);
     out.mkdir();
     String[] f = con.sortLs();
     for(int i=0; i<f.length; i++){
         System.out.println("changing to "+base+"/"+f[i]);
        try{
         con.chdir(base+"/"+f[i]);
        }catch(Exception exc){
            continue;
        }
         List<String> st = new ArrayList<String>(Arrays.asList(con.sortLs()));
         st.remove("./");
         st.remove("../");
         int ind = whichIndex("cover", st);
         if(ind <0) continue;
         System.err.println(st+" "+ind);;
         String file = base+"/"+f[i]+"/"+st.remove(ind);
         
         String[] res = con.read(file, "Job description");
        
         if(res!=null){
             if(res[0].trim().equals(experiment) && res.length>1){
                 System.err.println(Arrays.asList(res));
                 File outFile = new File(out, res[1].trim());
                 outFile.mkdir();
               
                 for(int ik=0; ik<st.size(); ik++){
                     File outF = new File(outFile, st.get(ik));
                     if(!outF.exists()){
                         con.download(st.get(ik), out+"/"+res[1].trim()+"/");
                     }
                 }
                          } else{
                 System.err.println("no match "+res[0]+" "+experiment);
             }
         }
     }
 }
 public void run(){
     String[] f = con.sortLs();
     System.out.println(Arrays.asList(f));
     
 }

 // download welcome.msg from sourceforge or any other given file
/* public static void main(String argv[])
 {
    if(argv.length < 2)
    {
        //FTPDownload f = new FTPDownload("ftp.kernel.org", "/welcome.msg");

        long x = 0;

        for(int i=0;i<5;i++) {
            FTPDownload f = new FTPDownload("localhost", "dsm.pdf");
            x += f.time;
        }

        System.out.println("5 runs took "+x+" ms, "+(long) x/5+" ms average.");
    }
    else
    {
        FTPDownload f = new FTPDownload(argv[0], argv[1]);
    }
 }*/

// ------------------ needed by ConnectionListener interface -----------------

// called if the remote directory has changed
 public void updateRemoteDirectory(BasicConnection con)
 {
    System.out.println("new path is: " + con.getPWD());
 }

 // called if a connection has been established
 public void connectionInitialized(BasicConnection con)
 {
    isThere = true;
 }
 
 // called every few kb by DataConnection during the trnsfer (interval can be changed in Settings)
 public void updateProgress(String file, String type, long bytes) {}

 // called if connection fails
 public void connectionFailed(BasicConnection con, String why) {System.out.println("connection failed!");}

 // up- or download has finished
 public void actionFinished(BasicConnection con) {}


 // ------------ needed by Logger interface  --------------

    // main log method
    public void debug(String msg) {} //{System.out.println(msg);}

    // rarely used
    public void debugRaw(String msg) {}//System.out.print(msg);}

    // methods below are not used yet.

        public void debug(String msg, Throwable throwable) {}

        public void warn(String msg) {}

        public void warn(String msg, Throwable throwable) {}

        public void error(String msg) {}

        public void error(String msg, Throwable throwable) {}

     public void info(String msg) {}

        public void info(String msg, Throwable throwable) {}

        public void fatal(String msg) {}

        public void fatal(String msg, Throwable throwable) {}

        public boolean isEnabled(int arg0) {
            // TODO Auto-generated method stub
            return false;
        }

        public void log(int arg0, String arg1) {
            // TODO Auto-generated method stub
            
        }

}
