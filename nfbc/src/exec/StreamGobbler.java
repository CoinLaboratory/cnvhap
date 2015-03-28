package exec;


	import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Writer;
public	class StreamGobbler extends Thread
	{
	    InputStream is;
	    String type;
	   Writer os;
	    
	 public   StreamGobbler(InputStream is, String type)
	    {
	        this(is, type, null);
	    }
	   public  StreamGobbler(InputStream is, String type, Writer redirect)
	    {
	        this.is = is;
	        this.type = type;
	       this.os = redirect;
	    }
	    
	  int cnt=0;
	   public void run()
	    {
	        try
	        {
	            PrintWriter pw = null;
	            if (os != null)
	              pw = new PrintWriter(os);
	                
	            InputStreamReader isr = new InputStreamReader(is);
	            BufferedReader br = new BufferedReader(isr);
	            String line=null;
	            while ( (line = br.readLine()) != null)
	            {
	               if (pw != null)
	                  pw.println(line);
	            //    System.out.println(cnt+" "+type + ">" + line);    
	                cnt++;
	            }
	           if (pw != null)
	               pw.flush();
	        } catch (IOException ioe)
	            {
	            ioe.printStackTrace();  
	            }
	    }
	}
	


