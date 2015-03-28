package exec;

import java.io.BufferedReader;
import java.io.Writer;


	public class Exec
	{
		
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

