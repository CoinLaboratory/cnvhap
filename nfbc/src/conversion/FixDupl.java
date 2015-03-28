package conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;

public class FixDupl {
File in;
PrintWriter out;
public void FixDupl(File in) throws Exception{
	File outf = new File(in.getParentFile(), in.getName()+"_1");
	BufferedReader br = new BufferedReader(new FileReader(in));
	
	
}

}
