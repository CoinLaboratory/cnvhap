package lc1.dp.appl;

import java.io.File;
import java.lang.reflect.Field;
import java.util.Arrays;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.LikelihoodDataCollection;
import lc1.possel.HLALD;
import lc1.util.Constants;

import org.apache.commons.cli.CommandLine;

public class CalcLD {
public static void main(String[] args){
	try{
		
		
		
           File f = new File(".");
       	String[] str1 = f.getAbsolutePath().split("/");
		//allnfbc_1_9_102428377_109966491/res/allnfbc
		String[]str = str1[str1.length-2].split("_");
		String chr =str[str.length-3];
		int start =  Integer.parseInt(str[str.length-2]);
		int end =  Integer.parseInt(str[str.length-1]);
	//	end = start +(end-start)/5;
		int len = args.length;
		  String[] args2 = new String[] {"--paramFile", "param.txt", "--column", "1", "--mid ", chr+":"+
					 start+":"+end};
		  String[]  args1 = new String[len+args2.length];
		  System.arraycopy(args, 0, args1, 0, len);
		  System.arraycopy(args2, 0, args1, len, args2.length);
     	   CommandLine par1 = Constants.parseInitial(args1);
     	   CommandLine para = Constants.parse(args1,1,1,  null);
     	  Field fie = Constants.class.getField("mid");
     	  //Constants.restrictKb(0);
      	fie.set(null, new String[][] {args1[args1.length-1].split(":")});
    //  	Constants.ldtype = new String[] {"2;0"};
      	Arrays.fill(Constants.loess,false);
    	Arrays.fill(Constants.gc,false);
     	   
	CalcLD ld = new CalcLD(f, str[0],chr,start, end );
	ld.run();
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
	DataCollection dc;
	File f;
	CalcLD(File f, String nme, String chr, int start, int end) throws Exception{
		this.f = f;
	
	
		File f1 = new File(f,"res/"+nme+"/geno/"+chr+".zip");
		int[] rest = Constants.restrictKbMax();
		int[] mid = new int[] {start-rest[0], end+rest[1] };
		this.dc =null;/* new  LikelihoodDataCollection(
				f1, (short) 0, chr.indexOf('X') >= 0
						|| chr.indexOf('Y') > 0 ? 1 : 2,
				new int[][] { mid },null);*/
	}
	
	public void run(){
		HLALD.calcLD(dc, f, false, false);
	}
	
}
