package lc1.pca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import lc1.util.Compressor;
import Jama.Matrix;
import Jama.SingularValueDecomposition;
import cern.colt.matrix.linalg.Algebra;

public class CalcPCA {

	public static void main(String[] args){
		try{
		File  dir1 = new File(System.getProperty("user.dir"));
	
		CalcPCA reg = new CalcPCA(dir1, "Log R");
		reg.read();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	File[] zipFiles;
	
	File dir;
	
	Matrix X;
	
	
	List<String> indiv = new ArrayList<String>();
	ZipFile zf = null;
	CalcPCA(File dir, String colId) throws Exception{
		final Pattern p = Pattern.compile("[0-9]");
		zipFiles = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File arg0) {
				String name = arg0.getName();
				Matcher mat = p.matcher(name);
				
				mat.useAnchoringBounds(false);
				boolean m1 = mat.find();
			//	boolean m = p.matcher(name).matches();
				return name.endsWith("zip") && m1;// && name.startsWith("1.");
			}
			
		});
		int k=0;
		int numEnt = 0;
		zf = new ZipFile(zipFiles[0]);
		indiv = Compressor.getIndiv(zf,"Samples",0);//.subList(0, 100);
		//indiv = indiv.subList(0, toIndex)
		String[] str = Compressor.readZipFrom(zf, "Name").get(0).split("\t");
		int lrrCol =0;
		for(;lrrCol<str.length; lrrCol++){
			if(str[lrrCol].indexOf(colId)>=0) break;
		}
		this.lrrCol = lrrCol;
		row = new double[indiv.size()];
		numRows = indiv.size();
		
		for(int j=0; j<zipFiles.length; j++){
			numEnt+=zf.size() - 2;
//			if(zf.getEntry("Names")!=null) numEnt--;
			if(zf.getEntry("SNPS")!=null) numEnt--;		
	//		if(zf.getEntry("Samples")!=null) numEnt--;
			zf.close(); 
			if(j<zipFiles.length-1)
			zf = new ZipFile(zipFiles[j+1]);
		}
		X = transpose ?
				new Matrix(numEnt, indiv.size()):
			new Matrix(indiv.size(),numEnt);
	}
	boolean transpose = true;
	
	static List<String> extr = Arrays.asList(new String[] {"Name", "SNPS", "Samples"});
	int k=0;
	final double[] row;
	final int lrrCol;
	final int numRows;
	public void read() throws Exception{
		//int ind =0;
		zf = new ZipFile(zipFiles[0]);
		indiv = Compressor.getIndiv(zf,"Samples",0);
		int numEnt =0;
		
		for(int j=0; j<zipFiles.length; j++){
			Enumeration en = zf.entries();
		
			for(; en.hasMoreElements(); ){
				ZipEntry z = (ZipEntry) en.nextElement();
				String nme = z.getName();
				if(!extr.contains(nme)){
					readZip(zf, nme, row, lrrCol);
					for(int i=0; i<numRows; i++){
						if(transpose)
						X.set(k, i, row[i]);
						else
							X.set(i, k, row[i]);
					}
					k++;
				}
			}
			numEnt+=zf.size();
			if(zf.getEntry("Names")!=null) numEnt--;
			if(zf.getEntry("SNPS")!=null) numEnt--;		
			if(zf.getEntry("Samples")!=null) numEnt--;
			zf.close(); 
			if(j<zipFiles.length-1)
			zf = new ZipFile(zipFiles[j+1]);
		}
		//DoubleMatrix2D matrix1 = Statistic.covariance(X);
	      
		  SingularValueDecomposition svd =   X.svd();
		 double[] d =  svd.getSingularValues();
//			  new SingularValueDecomposition(X.transpose());//alg.transpose(X));
		  System.err.println("h "+svd.rank());
			
	}
	Algebra alg = new Algebra();
	
	public static void readZip(final ZipFile f, String st, double[] res, int col){
        try{
          BufferedReader br = new  BufferedReader(new InputStreamReader(f.getInputStream(f.getEntry(st))));
          String str = "";
          double sum=0;
          int cnt=0;
          for(int i=0; (str = br.readLine())!=null && i<res.length;i++){
              String[] stri = str.split("\t");
              double v = Double.parseDouble(stri[col]);
              if(!Double.isNaN(v)){
            	sum+=v;
            	cnt++;
              }
            //  else{
            //	  res[i] = Double.NaN;
            //  }
              res[i] =v;
          }
          double avg = sum/cnt;
          br.close();
          sum=0;
          cnt =0;
          for(int i=0; i<res.length; i++){
        	  double v = res[i];
        	  if(!Double.isNaN(v)){
              	sum+=Math.pow(v -avg,2);
              	cnt++;
              }
          }
          double stddev = Math.sqrt(sum/cnt);
          for(int i=0; i<res.length; i++){
        	  double v = res[i];
        	  if(!Double.isNaN(v)){
              	res[i] = (v - avg)/stddev;
              }
        	  else{
        		  res[i] = 0;
        	  }
          }
        }catch(Exception exc){
            exc.printStackTrace();
            //return null;
        }
    }
	
	
	
}
