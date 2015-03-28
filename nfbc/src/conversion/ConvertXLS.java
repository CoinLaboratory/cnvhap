package conversion;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;

public class ConvertXLS {
	
	public static void main(String[] args){
		try{
		File dir = new File(System.getProperty("user.dir"));
		File xls = new File(dir, args[0]);
	
	//	for(int i=1; i<args.length; i++){
			//String[] snpid = args[i].split(":");
		Workbook wb = Workbook.getWorkbook(xls);
		int len =  wb.getNumberOfSheets();
		for(int i=0; i<len; i++){
			Sheet sh = wb.getSheet(i);
			String names = sh.getName();
			if(names.startsWith("SNPS")){
				if(true){
		ConvertXLS cxls = new ConvertXLS(xls,names);
		//		}
		(new CompressDir1(cxls.dir)).run();
		}
		}
		}
		
	//	CompressDir.delete(cxls.dir);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
File dir;
Workbook workbook;
//Sheet[] sheets;
//List<String> indiv;
List<List<String>> name = new ArrayList<List<String>>();
List<String> samples;
public ConvertXLS(File in,  String snp_name) throws Exception{
	File par = in.getParentFile();
	String name1 = snp_name.split("_")[1];
	File dir1 = new File(par, in.getName().split("\\.")[0]);
	 dir1.mkdir();
	File   dir2 = new File(dir1, name1);
	 dir2.mkdir();
	Workbook wb = Workbook.getWorkbook(in);
	Sheet snp_sheet = wb.getSheet(snp_name);
	Cell[] snp_header = snp_sheet.getRow(0);
	int snp_id = indexOf(snp_header,"id");
	int chrom_id = indexOf(snp_header,"chr");

	//int row_id = snpid==null ? 1 : indexOf(snp_sheet.getColumn(snp_id),snpid.get(0));
	//Sheet rs_sheet = wb.getSheet(snp_sheet.getRow(row_id)[snp_id].getContents());
	String chrom = snp_sheet.getRow(1)[chrom_id].getContents();
	dir = new File(dir2,chrom.substring(3));
	dir.mkdir();
	//Cell[] row = rs_sheet.getRow(0);
	PrintWriter pw_name = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Name"))));
	PrintWriter pw_snp = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "SNPS"))));
	List<String> snps = new ArrayList<String>();
	List<List<Double>> bins = new ArrayList<List<Double>>();
//	List<Double> bins1 = new ArrayList<Double>();
//	int bin_id = -1;
	pw_name.println("Genotype\tOriginal");
	pw_name.println("chr\tstart\tend\tid\tbins");
		
	

	for(int i=1; i<snp_sheet.getRows(); i++){
		Cell[] cell = snp_sheet.getRow(i);
		if(cell[0].getContents().startsWith("#")) continue;
		if(cell.length<=1) continue;
		PrintWriter pw_snp1 = i==0 ? pw_name : pw_snp;
	//	if(snpid==null || snpid.contains(cell[snp_id].getContents()) || i==0){
	
		if(i>0){
		
				snps.add(cell[snp_id].getContents());
				pw_snp1.print(cell[0].getContents());
				for(int k=1; k<4; k++){
					pw_snp1.print("\t"+cell[k].getContents());
				}
				if(cell.length>4 && cell[4].getContents().length()>0){
					List<Double> bins_ = getBins(Arrays.asList(cell).subList(4, cell.length),wb.getSheet(cell[snp_id].getContents()));
					Collections.sort(bins_);
					bins.add(bins_);
					pw_snp1.print("\t"+bins_.get(0));
					for(int k=1; k<bins_.size(); k++){
						pw_snp1.print(":"+bins_.get(k));
					}
				}
				else{
					bins.add(null);
					
				}
				pw_snp1.println();
				
			//}
			
			
		}
		
	//}
	}
	Set<String> samples = new HashSet<String>();
	for(int i=0; i<snps.size(); i++){
		Sheet sheet = wb.getSheet(snps.get(i));
		if(sheet==null){
			sheet = wb.getSheet(snps.get(i).split("_")[0]);
		}
		Cell[] col = sheet.getColumn(0);
		for(int k=1; k<col.length; k++){
			if(col[k].getContents().length()>0)
				samples.add(col[k].getContents());
		}
	}
	this.samples = new ArrayList<String>(samples);
	PrintWriter pw_samp = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "Samples"))));
	for(int i=0; i<samples.size(); i++){
		pw_samp.println(this.samples.get(i));
	}
	pw_samp.close();
	
	pw_snp.close();
	pw_name.println("id");
	
	StringBuffer nulst =new StringBuffer("null");
	for(int i=1; i<this.name.size(); i++){
		nulst.append("\tnull");
	}
	null_string = nulst.toString();
	for(int i=0; i<snps.size(); i++){
		if(bins.get(i)==null || bins.get(i).size()<0){
		writeSheet(wb.getSheet(snps.get(i)), snps.get(i));
		}
		else{
			 
			Sheet sh = wb.getSheet(snps.get(i).split("_")[0]);
			
				writeSheet(sh,  bins.get(i) ,snps.get(i));
			pw_name.println(max);
		}
	}
	pw_name.close();
}
private List<Double> getBins(List<Cell> subl, Sheet sh) {
	List<Double> l = new ArrayList<Double>();
	if(subl.get(0).getContents().equals("all")){
		SortedSet<Double> d = new TreeSet<Double>();
		for(int k=1; k<sh.getRows(); k++){
			Cell[] rows = sh.getRow(k);
			for(int kk=1; kk<rows.length; kk++){
				String str = rows[kk].getContents();
				if(!str.equals("FAIL"))
				d.add(Double.parseDouble(str));
			}
		}
		
		return cluster(d);
	}
	else{
	for(int i=0; i<subl.size();i++){
		String c = subl.get(i).getContents();
		if(c.length()>0){
			l.add(Double.parseDouble(c));
		}
	}
	}
	return l;
}
private List<Double> cluster(SortedSet<Double> d) {
	if(d.size()>95) {
		double d1 = d.first();
		double d2 = d.last();
		double step = (d2 - d1)/95;
		List<Double> l = new ArrayList<Double>();
		for(int k=0; k<95; k++){
			l.add(d1+step*k);
		}
		return l;
	}
	
	
	return new ArrayList<Double>(d);
}
String null_string;
public void writeSheet(Sheet sh, String name) throws Exception{
	Map<String, String> m = readSheet(sh);
	PrintWriter pw_snp = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name))));
	for(int i=0; i<this.samples.size(); i++){
		String val = m.get(samples.get(i));
		if(val==null) val = null_string;
		pw_snp.println(val);
	}
	pw_snp.close();
	
}
public void writeSheet(Sheet sh,  List<Double> bins,String name) throws Exception{
	Map<String, String> m = readSheet(sh, bins);
	{
		PrintWriter pw_snp = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, name))));
		for(int i=0; i<this.samples.size(); i++){
			String val = m.get(samples.get(i));
			if(val==null) val = null_string;
			pw_snp.println(val);
		}
		pw_snp.close();
	}
	if(bins.size()>2){
		PrintWriter pw_snp = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "bins_"+name))));
		int i=0;
		for(; i<bins.size(); i++){
			if(i>0){
				pw_snp.print(bins.get(i-1)+"<"+"\t");
			}
			else{
				pw_snp.print("0<=\t");
			}
			pw_snp.println(conv.get(i+1)+"\t<="+bins.get(i));
		}
		pw_snp.println(bins.get(i-1)+"<\t"+conv.get(i+1));
		pw_snp.close();
	}
	
}
public Map<String, String> readSheet(Sheet sh){
	Map<String, String> m = new HashMap<String, String>();
	for(int k=1;k<sh.getRows(); k++){
		Cell[] ce = sh.getRow(k);
		String id = ce[0].getContents();
		StringBuffer sb = new StringBuffer(ce[1].getContents());
		for(int kk=2; kk<=name.size(); kk++){
			sb.append("\t"+ce[kk].getContents());
		}
		m.put(id, sb.toString());
		
	}
	return m;
}

public Map<String, String> readSheet(Sheet sh,List<Double>bins){
	Map<String, String> m = new HashMap<String, String>();
	Cell[] cell = sh.getRow(0);
	
	int len =0;
	for(;len<cell.length; len++){
		if(cell[len].getContents().length()==0) break;
	}
	for(int k=1;k<sh.getRows(); k++){
		
		Cell[] ce = sh.getRow(k);
		String id = ce[0].getContents();
		if(id.length()==0) break;
		String[] contents = transform(ce[1].getContents(),bins);
		StringBuffer sb1 = new StringBuffer(contents[0]);
		StringBuffer sb2 = new StringBuffer(contents[1]);
		for(int kk=2; kk<len; kk++){
			contents = transform(ce[kk].getContents(),bins);
			sb1.append(contents[0]);
			sb2.append(";"+contents[1]);
		}
		m.put(id,sb1.toString()+"\t"+sb2.toString());
		
	}
	return m;
}
int max=0;
private String[] transform(String contents, List<Double> l) {
	// TODO Auto-generated method stub
	//if(true) return contents+"\t";
	if(contents.toLowerCase().equals("fail")) return new String[] {"null", contents};
 double d = Double.parseDouble(contents.replaceAll(",", ""));
 if(Double.isNaN(d)) return null;
 int i=0;
 for(; i<l.size(); i++){
	 if(d <= l.get(i)) break;
 }
if(i+1>max) max = i+1;
 String str =  ""+conv.get(i+1);
 return new String[] {str, contents};
}

static Convert conv = new Convert();
static class Convert{
static int A = (int)'A';
static int Z = (int) 'Z';
static int _ = (int) '_';
public int get(char ch){
	return l.indexOf(ch);
}
public char get(int i){
	return l.get(i);
}
List<Character> l = new ArrayList<Character>();
	Convert(){
		l.add('_');
		for(int i = A; i<=Z; i++){
			l.add((char) i);
		}
	 for(int i=31; i<126; i++){
		 if(!l.contains((char) i)){
			 l.add((char)i);
		 }
	 }
	}

}
private int indexOf(Cell[] sh, String string) {
	for(int k=0; k<sh.length; k++){
		if(sh[k].getContents().equals(string)) return k;
	}
	return -1;
}
}
