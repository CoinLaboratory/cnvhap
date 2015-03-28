package lc1.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;

public class ConvertXLSToCSV {
public static void main(String[] args){
	try{
		//if(true)return;
		String delim = "\t";
		File[] f = ( new File(System.getProperty("user.dir"))).listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && ((new File(pathname, "data.xls"))).exists();
			}
			
		});
		for(int k=0; k<f.length; k++){
			Workbook wb = Workbook.getWorkbook(new File(f[k], "data.xls"));
			String[] nme = wb.getSheetNames();
			File out = new File(f[k], "data");
			out.mkdir();
			inner: for(int j=0; j<nme.length; j++){
			// String[] exp = experiment;
				if(!nme[j].equals("paper")) continue inner;
				Sheet ws = wb.getSheet(nme[j]);
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(out, nme[j]+".csv"))));
				for(int i=0; i<ws.getColumns(); i++){
					Cell[] c = ws.getColumn(i);
					for(int l=0; l<c.length; l++){
						pw.print(c[l].getContents()+delim);
					}
					pw.println();
				}
				pw.close();
			}
		}
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
