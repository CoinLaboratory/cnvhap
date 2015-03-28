package nfbc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class Analysis {
    File file;
    
    public Analysis(File f) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(f));
    }
}
