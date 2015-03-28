package exec;


	import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
public 	class StreamGobbler1 extends Thread
	{
	    OutputStream is;
	    String type;
	    BufferedReader br;
	  //  OutputStream os;
	    
	  public  StreamGobbler1(OutputStream is, String type, BufferedReader redirect)
	    {
	        this.is = is;
	        this.type = type;
	        this.br = redirect;
	       // this.os = redirect;
	    }
	    
	 int cnt=0;
	  public void run()
	    {
	        try
	        {
	           // PrintWriter pw = null;
	           // if (os != null)
	           //     pw = new PrintWriter(os);
	                
	            OutputStreamWriter isr = new OutputStreamWriter(is);
	          //  BufferedWriter br = new BufferedWriter(isr);
	            String line=null;
	            while ( (line = br.readLine()) != null)
	            {
	             //   if (pw != null)
	            	isr.write(line);
	            	isr.write("\n");
	            	isr.flush();
	              //      pw.println(line);
	             //   System.out.println(cnt+" "+type + ">" + line);   
	                cnt++;
	            }
	            isr.close();
	           // br.close();
	           // if (pw != null)
	            //    pw.flush();
	        } catch (IOException ioe)
	            {
	            ioe.printStackTrace();  
	            }
	    }
	}
	


