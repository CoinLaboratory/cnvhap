import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;


public class ExtractDataFiles {
	public static void main(String[] args){
		try{
			String d = System.getProperty("user.dir");
	        File dir1 = new File(d);
			ExtractDataFiles sif = new ExtractDataFiles(dir1, args[0], args[1]);
			sif.extract();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	private void extract() {
		for(int i=1; i<include.length; i++){
			if(include[i].getContents().equals("true")){
				pw.println(input[i].getContents());
			}
		}
		pw.close();
	 wb.close();
	}
	PrintWriter pw;
	Cell[] input ;
	Cell[] include;
	 Workbook wb ;
	ExtractDataFiles(File dir, String name, String experiment) throws Exception{
		String st = "";
		File ff = new File(dir, name);
		pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "dataFiles.txt"))));
		wb= Workbook.getWorkbook(ff);
    	
    	Sheet ws = wb.getSheet(experiment);
    	 for(int i=0; i<ws.getRows(); i++){
    		Cell[] c =  ws.getRow(i);
    		 if(c[0].getContents().startsWith("inputDir")){
    			 input = c;
    		 }
    		 else if(c[0].getContents().startsWith("include")){
    			 include = c;
    		 }
    	 }
    	
	}
	
}
