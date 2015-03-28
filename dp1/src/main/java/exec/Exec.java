package exec;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;


	public class Exec
	{
		
		public static void main(String[] args){
			try{
				 Runtime rt = Runtime.getRuntime();
				 String[] exec = new String[] {"R","--vanilla"};
		         Process proc = rt.exec(exec);
		         Writer err = new OutputStreamWriter(System.out);
		        
		         PipedInputStream pi = new PipedInputStream();
		         PipedOutputStream po = new PipedOutputStream();
		         po.connect(pi);
		    	 final    OutputStreamWriter osw = new OutputStreamWriter(po);
		     
		       //  PrintWriter pw = new PrintWriter(osw);
		         BufferedReader in = new BufferedReader(new InputStreamReader(pi));
		         // any error message?
		         StreamGobbler errorGobbler = new 
		             StreamGobbler(proc.getErrorStream(), "ERROR", err);            
		         
		         // any output?
		         StreamGobbler outputGobbler = new 
		             StreamGobbler(proc.getInputStream(), "OUTPUT", err);
		        
		         StreamGobbler1  inputGobbler = new StreamGobbler1(proc.getOutputStream(), "INPUT",in);
		         // kick them off
		      
		     Thread th = new Thread(){
		    	 public void start(){
		    		 try{
		     for(int i=0;i<5; i++){
		    	System.err.println("here"+i);
		    	    osw.write("x = "+i+"\n"); 
		    	    osw.write("y = 2\n"); 
		    	    osw.write("print(x+y)\n"); 
		    	//    osw.close();
		    	    osw.flush();
		    	//    po.flush();
		    	    
		    	  //  System.err.println(err.toString());
		     }
		    		 }catch(Exception exc){
		    			 exc.printStackTrace();
		    		 }
		    	 }
		     };
		  
		     inputGobbler.start();         
	         errorGobbler.start();
	      outputGobbler.start();
	      th.start();
		     osw.close();
		     System.err.println(err.toString());
		         errorGobbler.join();
		         inputGobbler.join();
		         outputGobbler.join();
		         th.join();
		         err.close();
		        in.close();
		       // System.err.println(err.toString());
		    
		         // any error???
		         int exitVal = proc.waitFor();
		          System.out.println("Process exitValue: " + exitVal);
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
		
		public static void exec(String[] exec1, BufferedReader br1, Writer err, Writer err2) throws Exception{
			  Runtime rt = Runtime.getRuntime();
		         Process proc = rt.exec(exec1);
		         // any error message?
		         StreamGobbler errorGobbler = new 
		             StreamGobbler(proc.getErrorStream(), "ERROR", err);            
		         
		         // any output?
		         StreamGobbler outputGobbler = new 
		             StreamGobbler(proc.getInputStream(), "OUTPUT", err);
		        
		         StreamGobbler1  inputGobbler = new StreamGobbler1(proc.getOutputStream(), "INPUT",br1);
		         // kick them off
		      inputGobbler.start();         
		         errorGobbler.start();
		      outputGobbler.start();
		     
		         errorGobbler.join();
		         inputGobbler.join();
		         outputGobbler.join();
		         err.close();
		         br1.close();
		         // any error???
		         int exitVal = proc.waitFor();
		          System.out.println("Process exitValue: " + exitVal);
		  }
		
		
	  
	}

