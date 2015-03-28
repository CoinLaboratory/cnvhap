package lc1.util;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

public class PrintLatexTable {
static void print(File f, String[] col_head, String[] row_head, double[][] d){
	try{
	PrintWriter pw = new PrintWriter(new FileWriter(f));
	StringBuffer format = new StringBuffer("l");
	for(int k=0; k<col_head.length; k++){
		format.append("c");
	}
	pw.println("\\begin{tabular}{"+format+"}");
	
	for(int j=0; j<col_head.length; j++){
		pw.print("\t&"+col_head[j]);
	}
	pw.println("\\\\");
	for(int i=0; i<d.length; i++){
		pw.print(row_head[i]);
		for(int j=0; j<d[i].length; j++){
			pw.print("\t&"+String.format("%5.3g", d[i][j]));
		}
		pw.println("\\\\");
	}
	pw.println("\\end{tabular}");
	pw.close();
	
	}catch(Exception exc){
		exc.printStackTrace();
	}
}

public static void print(PrintWriter pw, String top, String[] col_head, String[] row_head, int[][] d){
	try{
	//PrintWriter pw = new PrintWriter(new FileWriter(f));
	StringBuffer format = new StringBuffer("|l|");
	for(int k=0; k<col_head.length; k++){
		format.append("c|");
	}
	
	pw.println("\\begin{tabular}{"+format+"}");
	pw.println("\\hline");
	pw.println("& \\multicolumn{"+col_head.length+"}{|c|}{"+top+"} \\\\");
	pw.println("\\hline");
	for(int j=0; j<col_head.length; j++){
		pw.print("\t&"+col_head[j]);
	}
	pw.println("\\\\");
	pw.println("\\hline");
	for(int i=0; i<d.length; i++){
		pw.print(row_head[i]);
		for(int j=0; j<d[i].length; j++){
			pw.print("\t&"+String.format("%5.3g", (double) d[i][j]));
		}
		pw.println("\\\\");
		
	}
	pw.println("\\hline");
	pw.println("\\end{tabular}");
//	pw.close();
	pw.println();
	
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
