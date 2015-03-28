package calc;

import java.io.File;

import jxl.Workbook;

public class TransferExcelToCsv {

	
	public static void main(String[] args){
		try{
			File ff = new File("CNV data.xlsx");
//			File ff = new File(args[0]);
			Workbook wb = 
					Workbook.getWorkbook(ff);
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
}
