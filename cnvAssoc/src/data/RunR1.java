package data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipFile;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.Rengine;

import cnvtrans.LoopCall;


public class RunR1 {
	 Integer start;// =21856055;//0;//156000000;// 0;//	200758448;
 Integer end;
	static int max_num = Integer.MAX_VALUE;
	public final int maxIndiv;
	BufferedReader r;
	static Rengine re = new Rengine(new String[] { "--vanilla" }, true,
			new LoopCall());
			
	int cnt =0;
	
	int count=0;
	final boolean pleio,multiGen;
	
	
	
	
	
	
	
	public boolean nextLine() throws Exception{
		if(count>max_num) return false;
		String current = this.snps.readLine();
		
		if(current==null) return false;
		this.curr_snp =  current.split("\\s+");
		count++;
		return true;
	}
	 double[][] weights;
	 double rescale;
	
	 
	 public void processLine(String snp_id)  throws Exception{
		
		//re.assign("rsid",new String[] {curr_snp[0], curr_snp[1], curr_snp[3]});
		this.read(zf, snp_id, geno);
		re.assign("rsid", snp_id);
		transf(geno,weights,rescale);
		assign(geno,"genotmp",null,false,colsToInclude);
	
		eval("if(length(datanme$alleleCols)>0) genotmp = mergeAlleleCols(genotmp,datanme$alleleCols)");
		eval("if(length(datanme$gCols)>0) genotmp = fixGenoCols(genotmp, gCols)");
			
				
		eval("genotmp = fixNonNumeric(genotmp)");
		eval("datanme$geno[,1,] = t(genotmp)");
		
		
		if(this.expandData>1) {
			assign(weights, "datanme2$geno[,2,]",null,true,colsToInclude);
		}
		File snps = new File("snps");
		if(CHECK && !snps.exists()){
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(snps)));
			for(int i=0; i<geno[0].length; i++){
				for(int k=0; k<geno.length; k++){
					pw.print(geno[k][i]);
					if(k<geno.length-1) pw.print("\t");
				}
				pw.println();
			}
			pw.close();
		    REXP re1 = 	re.eval("datanme$geno");
			System.err.println();
		}
	}
	
	public void processGenoCov(String snp_id, int[] colsToInclude) throws Exception{
		this.read(zf, snp_id, geno);
		transf(geno,weights,rescale);
		assign(geno,"genotmp",null,false,colsToInclude);
	}

	
	private void transf(String[][] geno2, double[][] weights, double rescale) {
		//Boolean getAvg = false;
		int len = geno2[0].length;
		boolean expand = expandData>1;
		for(int k=0; k<geno2.length; k++){	
			if(expand){
				   Arrays.fill(weights[k],0);
				}
			for(int j=0; j<len; j++){
				String str = geno2[k][j];
				String repl=str;
				if(str.length()==0 && def_values_for_geno!=null){
					str = def_values_for_geno[k];
				}
				if(str.indexOf(',')>=0){
					repl	=""+getAvg(str);
					if(expand){
						String[] spl = str.split(",");
						int max_ind =-1;
						double sum=0;
						double max = 0;
						for(int kk=0; kk<spl.length; kk++){
							if(spl[kk].equals("-")) max_ind = kk;
							else{
								int v =(int) Math.floor(spl[kk].length()==0 ? 0 : Integer.parseInt(spl[kk])* (rescale/this.scaling));
								
								sum+=v;
								weights[k][kk*len+j] = v;
								if(max_ind<0 && v> max) max_ind = kk;
							}
						}
						if(sum>rescale) throw new RuntimeException("!!");
						if(sum<rescale){
							
							weights[k][max_ind*len+j] = rescale-sum;
						}
					}
				}else if(expand){
					double d = Double.parseDouble(repl);
					int floor = (int) Math.floor(d);
					int ceil = (int) Math.ceil(d);
					if(ceil==floor) weights[k][floor*len+j] = Math.round(1.0*rescale);
					else{
						double fl =Math.round((ceil -d)*rescale); 
						weights[k][floor*len+j] = fl;
						weights[k][ceil*len +j] = rescale-fl;
					}
				}
				geno2[k][j] = repl;
			}
		}
		
	}
	public double getAvg(String str){
		String[] str1 = str.split(",");
		int maxind = -1;
		double tot = 0;
		double cnt=0;
		for(int k=0; k<str1.length; k++){
			if(str1[k].length()>0){
				if(str1[k].equals("-")) maxind = k;
				else{
				   tot+= Double.parseDouble(str1[k])*(double)k;	
				   cnt+=Double.parseDouble(str1[k]);
				}
			}
		}
		if(maxind>=0){
			tot+=(scaling-cnt)*maxind;
		}
		double res = tot/scaling;
		if(inverse){
			res = Math.round(res);
		}
		return res;
	}
	
String chrom;
public void setChrom(String chrom){
	this.chrom = chrom;
}

	//int[] alias;
	final String[][] geno ;
	String[] geno_header = null;
	
	public static boolean CHECK  =false;

	public void initialise() throws Exception{
		String st="";
		while((st=this.r.readLine())!=null){

			if(st.indexOf("GETGENOCOV")>=0) break;
			if(st.indexOf("POSTMERGE")>=0) break;
			if(st.indexOf("REPEAT")>=0) break;
			if(st.length()>0 && ! st.startsWith("#")){
				eval(st);
			}
		}
		boolean cont = st.indexOf("GETGENOCOV")>=0;
		if(cont){
			List<String> torun = new ArrayList<String>();
			while((st=this.r.readLine())!=null){
				if(st.indexOf("ENDGETGENOCOV")>=0) break;
				else{
					if(!st.startsWith("#"))torun.add(st);
				}
			}
			 String[] genocovs = eval("genocov").asStringArray();
			if(genocovs!=null && genocovs.length>0){
			 double[] inds = eval("genoTrans").asDoubleArray();
		
			
			for(int k=0; k<genocovs.length; k++){
				 int[] colsToInclude1  = new int[]{(int)inds[k]-1};
				// Arrays.fill(colsToInclude1,false);
				// colsToInclude1[(int)inds[k]-1] = true;
				 this.processGenoCov(genocovs[k], colsToInclude1);
				 Iterator<String> it = torun.iterator();
				 for(;it.hasNext();){eval(it.next());}
			}
			}
			 
			while((st=this.r.readLine())!=null){
				
				if(st.indexOf("POSTMERGE")>=0) break;
				if(st.indexOf("REPEAT")>=0) break;
				if(st.length()>0 && ! st.startsWith("#")){
					eval(st);
				}
			}
		}
		if(CHECK){
		String[] fams = eval("families").asStringArray();
		double[][] phenos = eval("datanme$pheno").asDoubleMatrix();
		String[] nme = eval("dimnames(datanme$pheno)[[2]]").asStringArray();
		REXP inclM = eval("datanme$incl");
		REXP lrr_ind = eval("datanme$ind");
		REXP phenoCovna = eval("phenoCovna");
		REXP varna = eval("varna");
//		REXP re1 = eval("which(caseInd[,1])");
//		REXP re2 = eval("which(!caseInd[,1])");
//		REXP re3 = eval("which(datanme$pheno[,1]==0)");
//		REXP re4 = eval("which(datanme$pheno[,1]==1)");
		System.err.println("h");
		}
	}
	
	public void premerge() throws Exception{
		String st="";
		while((st=this.r.readLine())!=null){
			if(st.indexOf("REPEAT")>=0) break;
			if(st.length()>0 && ! st.startsWith("#")){
				eval(st);
			}
		}
		if(CHECK){
			double[][] offs = eval("datanme$offsets").asDoubleMatrix();
			System.err.println("h");
		}
	}
	
	
	public  String[] evalJoinString(String type, int kk1, int k1, int ncol, String form){
		REXP re =  eval( "joinRes("+type+"[,"+kk1+","+k1+",],"+ncol+",\""+form+"\")");
		if(re==null) return null;
		else return re.asStringArray();
	}
	public  REXP evalPasteJoinString(String type, int kk1, int k1, int ncol, String form){
		return eval( "paste(joinRes("+type+"[,"+kk1+","+k1+",],"+ncol+",\""+form+"\"),collapse=\"/\")");
		
	}
	public  REXP evalJoinString(String type, int kk1, int k1, String form){
		return eval( "joinResVec("+type+"[,"+kk1+","+k1+"],\""+form+"\")");
		
	}
	
	
	final Map<String, String> snpsToDo;
	public void run(File outp) throws Exception{
		String st = "";
		//String[] sn = re.eval("snpstodo").asStringArray();
		
		double[] extralevelp = re.eval("extralevelp").asDoubleArray();
		double maf_thresh = eval("maf_thresh").asDouble();
		List<String> rep = new ArrayList<String>();
		List<String> rep1 = new ArrayList<String>();
		List<List<String>> rep2 = new ArrayList();
		while((st=this.r.readLine())!=null){
			if(st.indexOf("CONDITIONONMAF")>=0) break;
			if(st.length()>0 && ! st.startsWith("#")){
				rep.add(st);
			}
		}
		while((st=this.r.readLine())!=null){
			if(st.indexOf("EXTRALEVEL")>=0) break;
			if(st.length()>0 && ! st.startsWith("#")){
				rep1.add(st);
			}
		}
		while(st!=null){
			List<String> l = new ArrayList<String>();
			rep2.add(l);
			inner: while((st=this.r.readLine())!=null){
				if(st.indexOf("EXTRALEVEL")>=0) break inner;
				if(st.length()>0 && ! st.startsWith("#")){
					l.add(st);
				}
			}
		}
		String[] pheno1 = eval("dimnames(datanme$pheno)[[1]]").asStringArray();
		String[] pheno2= pleio ? eval("c(dimnames(datanme$pheno)[[1]],\"all\")").asStringArray(): eval("dimnames(datanme$pheno)[[1]]").asStringArray();
		 double[] spl = eval("spl").asDoubleArray();
		 String splFormat = getFormat(spl);
		String[] strats = eval("dimnames(datanme$strat)[[2]]").asStringArray();
		String[] type1 = eval("geno_header").asStringArray();
		String[] type2 = multiGen ? eval("c(geno_header,\"combined\")").asStringArray() : eval("geno_header").asStringArray();
		String[] family = eval("families").asStringArray();
		REXP comb = re.eval("singleOutput");
		REXP zipOutp = re.eval("zipOutput");
		//REXP maxn = re.eval("maxnum");
		//if(maxn!=null) max_num = maxn.asInt();
		boolean combinePhens = comb!=null && comb.asBool().isTRUE();
		boolean zip = zipOutp!=null && zipOutp.asBool().isTRUE();
		//REXP ph = eval("pheno");
		PrintWriter[][] pw = new PrintWriter[pheno2.length][type2.length];
		PrintWriter[][] pw_extra = new PrintWriter[pheno2.length][type2.length];
		File[][] outf_extra = new File[pheno2.length][type2.length];
		String pref = this.chrom+"_"+this.start+"_"+this.end;
		for(int kk=0; kk<pw.length; kk++){
			if(kk==0 || !combinePhens){
		for(int k=0; k<pw[kk].length; k++){
		  //for(int kj=0; kj<pw[kk][k].length; kj++){
			File outf = new File(outp,"res"+"_"+(combinePhens && pheno2.length>1 ? "all": pheno2[kk] )+pref+"_"+type2[k]+".txt"+(zip ? ".gz" : ""));
			 outf_extra[kk][k] = new File(outp,"extra_"+(combinePhens && pheno2.length>1 ? "all": pheno2[kk] )+pref+"_"+type2[k]+".txt"+(zip ? ".gz" : ""));
			//System.err.println("out f is "+outf+".gz");
		    pw[kk][k]= 	
		    	 new PrintWriter(new BufferedWriter(
		    			 zip ? new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf))):
		    				 new FileWriter(outf)
		    		));
		    pw_extra[kk][k] = 	
		    	 new PrintWriter(new BufferedWriter(
		    			 zip ? new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outf_extra[kk][k]))):
		    				 new FileWriter(outf_extra[kk][k])
		    		));
		    	
		    	//new PrintWriter(new BufferedWriter(new FileWriter(outf)));
		    pw[kk][k].println(splFormat);
		    pw[kk][k].println(Arrays.asList(strats)+",meta");
			pw[kk][k].println("pheno\tchrom\tstart\tend\trsid\tnosnps\tmaf_control\tmaf_case\thwe_control\thwe_case\tbeta\tpvalue\tnocontrol\tnocase\ttype");	
		//}
		}}else{
			pw[kk] = pw[0];
			pw_extra[kk] = pw_extra[0];
			 outf_extra[kk] = outf_extra[0];
		}
			
		}
		   int ncol = eval("dim(res)[4]").asInt();
		//   eval("print(length(spl))");
		   int ncol1 = eval("length(spl)").asInt()-1;
		while(nextLine()){
			if(snpsToDo!=null){
				String phenToExcl = snpsToDo.get(curr_snp[3]);
				if(phenToExcl==null) continue;
				else if(phenToExcl.equals("none")){
					eval("phenSubInds = 1:length(phenNames1)");
				}
				else{
					eval("phenSubInds = which(phenNames1!=\""+phenToExcl+"\")");
				//	int[] phe_sub_inds = this.re.eval("phenSubInds").asIntArray();
				//	System.err.println("h");
				}
			}
			//else{
			//	eval("phenSubInds = 1:length(phenNames1)");
			//}
			if(snpsToDo!=null && ! snpsToDo.containsKey(curr_snp[3])) continue;
			if(end!=null && Integer.parseInt(curr_snp[1])>end) continue;
			if(start==null || Integer.parseInt(curr_snp[1])>=start){
			//	kk++;
				try{
				processLine(curr_snp[rs_ind]);
				}catch(Exception exc){
					System.err.println("problem with "+curr_snp[rs_ind]);
					exc.printStackTrace();
				}
				for(int k=0; k<rep.size(); k++){
					eval(rep.get(k));
				}
				double maf = eval("max(maf)").asDouble();
				if(maf<maf_thresh) continue;
				for(int k=0; k<rep1.size(); k++){
					eval(rep1.get(k));
				}
				
				double[] minp = eval("apply(res[,,,2,drop=F],2,min,na.rm=T)").asDoubleArray();
			//	double[] minp2 = eval("apply(res2[,,,2,drop=F],2,min,na.rm=T)").asDoubleArray();
				double mp = 1.0;
				for(int k=0; k<minp.length; k++){
					if(!Double.isNaN(minp[k]) && minp[k]<mp) mp = minp[k];
				}
				if(minp==null){
					//REXP re = eval("as.matrix(res["+(kk+1)+",,])");
				   System.err.println("RES WAS NULL");
				   continue;
				}
				for(int j=0; j<rep2.size(); j++){
					if(mp<=extralevelp[j]){
						List<String>l = rep2.get(j);
						for(int k=0; k<l.size(); k++){
							eval(l.get(k));
						}
					}
				}
				
				
			if(CHECK){
				String[] str = eval("names(data)").asStringArray();
				REXP geno = eval("data");
			    REXP gens = eval("as.matrix(datanme$geno[,datanme$ind])");
			    double[] incl = eval("as.numeric(datanme$incl[,1])").asDoubleArray();
			    double[][] phens = eval("datanme$pheno").asDoubleMatrix();
			//    double[] phen = eval("as.numeric(datanme$pheno[,1])").asDoubleArray();
			 //   double[] offset = eval("datanme$offsets[,1]").asDoubleArray();
			  //  REXP cas = eval("which(datanme$caseInd[,1])");
			   // REXP cont = eval("which(!datanme$caseInd[,1] & !is.na(datanme$pheno[,1]))");
			 //  
			    REXP strat = eval("datanme$strat");
			    System.err.println("h");
			}
		
			
		 //for(int kj=0; kj<strats.length;kj++){
		    //	int kj1 = kj+1;
		 
		
			
		 inner1: for(int kk=0; kk<pheno2.length; kk++ ){
			int kk1 = kk+1;
			for(int k=0; k<type2.length; k++){
				int k1 = k+1;
			//	REXP re1 = eval("res");
				String[] res  = evalJoinString("res",kk1,k1,ncol, "%5.3g");
				//String[] res2 = evalJoinString("res2",kk1,k1,ncol, "%5.3g");
				//String[] res = (expandData>1 ?   res2 : res1);
				
				//	eval("joinRes(res[,"+kk1+","+k1+",],"+ncol+")");
		    	if(res==null ){
		    		System.err.println("was null " +pheno2[kk]);
		    		continue inner1;
		    	}
		   
		       REXP countsCase = null,countsControl = null, hweCase = null,hweControl = null, maf_case = null,maf_control = null;
				if(k<type1.length && minp[kk]<extralevelp[0] && kk<pheno1.length){
				   countsCase = evalPasteJoinString("countsCase",kk1,k1,ncol1,"%5.0f");
				   countsControl = evalPasteJoinString("countsControl",kk1,k1,ncol1,"%5.0f");
				  
			  
				  hweCase = evalJoinString("hwe_case",kk1,k1, "%5.3g");
			   	 hweControl = evalJoinString("hwe_control",kk1,k1, "%5.3g");
			 	  maf_case = evalJoinString("maf_case",kk1,k1, "%5.3g");
				   	 maf_control = evalJoinString("maf_control",kk1,k1, "%5.3g");
			   	//	   eval("joinRes(hwe_case[,"+kk1+",])");
				 //  hweControl = eval("joinRes(hwe_control[,"+kk1+",])");
				  // maf_case = eval("joinRes(maf_case[,"+kk1+",])");
				   //maf_control = eval("joinRes(maf_control[,"+kk1+",])");
				}
		
				if(true){
				pw[kk][k].println(pheno2[kk]+"\t"+curr_snp[0]+"\t"+curr_snp[1]+"\t"+curr_snp[2]+"\t"+curr_snp[3]+"\t"+
						(nosnp_ind < 0 ? "1" : curr_snp[this.nosnp_ind])+"\t"+
						format1("%5.3g", maf_control,k)+"\t"+
						format1("%5.3g", maf_case,k)+"\t"+
						format1("%5.3g", hweControl,k)+"\t"+
						format1("%5.3g", hweCase,k)+"\t"+
						format("%5.3g",res[0])+"\t"+
						format("%5.3g",res[1])+"\t"+
						//format("%5.3g",res[2])+"\t"+
						//format("%5.3g",res[3])+"\t"+
						format1("%5.0f", countsControl,k)+"\t"+
						format1("%5.0f", countsCase,k)+"\t"+
						(inverse || kk>=family.length ? "inverse":		family[kk])
						);
				}
				pw[kk][k].flush();
				if(kk<pheno1.length && k<type1.length && rep2.size()>=2 && minp[kk]<extralevelp[1] && eval("idsCase")!=null){
				//	double[] spl = eval("spl").asDoubleArray();
					String[]  idsCase = eval("joinRes(idsCase[,"+(kk+1)+","+(k+1)+",])").asStringArray();
					String[]  idsControl= eval("joinRes(idsControl[,"+(kk+1)+","+(k+1)+",])").asStringArray();
				//	String[]  idsControl1= eval("as.matrix(idsControl["+(kk+1)+",,"+k+"])").asStringArray();
					boolean print =false;
					inner: for(int j=0; j<idsCase.length; j++){
						if(idsControl[j].length()>0 || idsCase[j].length()>0 ){
							print = true; break inner;
						}
					}
					if(print){
						pw_extra[kk][k].println("pheno is "+pheno2[kk]);
						for(int j=0; j<idsCase.length; j++){
							pw_extra[kk][k].println("control, count="+j+": "+idsControl[j]);
							pw_extra[kk][k].println("case, counnt="+j+": "+idsCase[j]);
						}
					}
				}
			}
		    }
			}
			//}
		}
		for(int kk=0; kk<pw.length; kk++){
			for(int k=0; k<pw[kk].length; k++){
			//	for(int kj=0; kj<pw[kk][k].length; kj++){
				pw[kk][k].close();
				pw_extra[kk][k].close();
				if(outf_extra[kk][k].length()==0) outf_extra[kk][k].delete();
			//}
		}
		}
		r.close();
		re.end();
	}
	


	private String getFormat(double[] spl) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<spl.length; i++){
			sb.append(String.format("%5.3g ",spl[i]));
		}
		return sb.toString();
	}

	private String format(String string, String string2) {
	return string2.replaceAll(" ","");
	}
	private String format1(String string, REXP maf_control, int k) {
		if(maf_control==null || maf_control.asString()==null) {
			return "NA";
		}
		else {
			String res =  maf_control.asString().replaceAll(" ", "");
			return res;
		}
	}
	private String format(String form,REXP re, int k) {
		if(re==null || re.asDoubleMatrix()==null) return "NA";
		double[] ds = re.asDoubleMatrix()[k];
		if(ds==null) return "NA";
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<ds.length; i++){
			sb.append(String.format(form, ds[i]).trim());
			if(i<ds.length-1) sb.append(",");
		}
		return sb.toString();
	}

	private String getFormat(int length, String st) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<length; i++){
			sb.append(st);
			if(i<length-1)sb.append(",");
		}
		return sb.toString();
	}

	String[] curr_snp = null;
	
	

	public static void main(String[] args){
		long tim = System.currentTimeMillis();
		String zips = args[0];
		File f = new File(zips);
		int ind = f.getName().indexOf('*');
		int ind1 = f.getName().indexOf('[');
		int ind2 = f.getName().indexOf(']');
		String[] str = null;
		File pf = f.getParentFile();
		if(ind>=0){
		
			final String bef = f.getName().substring(0,ind);
			final String aft = f.getName().substring(ind+1);
			str = pf.list(new FilenameFilter(){

				@Override
				public boolean accept(File dir, String name) {
					// TODO Auto-generated method stub
					return name.startsWith(bef) && name.endsWith(aft);
				}
				
			});
		}else if(ind1>=0 && ind2>=0){
			final String bef = f.getName().substring(0,ind1);
			final String[] bet = f.getName().substring(ind1+1,ind2).split("-");
			int st = Integer.parseInt(bet[0]);
			int end = Integer.parseInt(bet[1]);
			final String aft = f.getName().substring(ind2+1);
			str = new String[end-st+1];
			for(int i=st; i<=end; i++){
				str[i-st] = bef+i+aft;
			}
		}
		else{
			str = new String[] {f.getName()};
		}
		for(int i=0; i<str.length; i++){
			args[0] = (new File(pf,str[i])).getAbsolutePath();
			main1(args);
		}
		long tim1 = System.currentTimeMillis();
		System.err.println("time in millisecondes "+(tim1-tim));
		System.exit(0);
	}
	/* args[3] is true or false */
	public static void main1(String[] args){
		
		System.err.println(Arrays.asList(args));
		try{
		
			String[] dirs = args[0].split(":");
			String[] pdirs = args[1].split(":");
		
			File[] dir = new File[dirs.length];
			File[] work = new File[dirs.length];
			//File user1 = new File(System.getProperty("user.dir"));
			for(int k=0; k<dirs.length; k++){
				dir[k] = new File(dirs[k]);
				work[k] = new File(pdirs[k]);
			}
			File rscript = new File(args[2]);
			
				String[] m = args[3].split(":");
			int[] mid = new int[] {Integer.parseInt(m[0]), m[1].equals("max") ? Integer.MAX_VALUE-10 : Integer.parseInt(m[1])};
			
				main(dir,work, rscript, mid);
			
		}
	catch(Exception exc){
		exc.printStackTrace();
		//System.exit(0);
	}
	
}
	
	
	
	public static void main(File[] f, File[] phenoDir, File rscript, int[] midr) throws Exception{
	
		String[] nme = new String[f.length];
		RunR1[] rr = new RunR1[nme.length];
		File[] zip = new File[nme.length];
		StringBuffer combname =new StringBuffer();
		for(int k=0; k<f.length; k++){
			String nme1 = phenoDir[k].getName().replace('-', '_');
			if(f.length==1)nme1 = "datanme";
		   rr[k] = new RunR1(f[k], rscript, phenoDir[k],midr,nme1);
		  combname.append(phenoDir[k].getName());
		  if(k<f.length-1) combname.append("_");
		}
		RunR1 merge;
		if(rr.length>1){
		   merge = new RunRMerged(rr);
		}else{
			merge = rr[0];
		}
		
		   merge.setChrom(rr[0].chrom);
		   merge.premerge();
		
		File resDir = new File("resDir");
		resDir.mkdir();
		File resDir1 = new File(resDir,combname.toString());
		resDir1.mkdir();
		merge.run(resDir1);
	/*	for(int k=0; k<rr.length; k++){
			  rr[k].premerge();
				rr[k].run(fs[k]);
		}*/
		
	}
	ZipFile zf;
	//Enumeration enume;
	
BufferedReader snps;
int rs_ind=-1;
int nosnp_ind = -1;
String nme;



String[] def_values_for_geno=null;


public double scaling;

//new RunR1(nme1,zip[k], expt, rscript, null, f[k].getName());
/* f is path to zip file.  working dir is path to pheno limit files etc */
int[] colsToInclude=null;
private BufferedReader getReader(String entry){
	 BufferedReader snps = null;
	try{
		snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(entry))));
	    }catch(Exception exc){
	    	System.err.println("no entry"+entry+", trying "+chrom+".zip");
	    	try{
	    		snps = new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(chrom+"/"+entry))));
	    	}catch(Exception exc1){
	    		System.err.println("still no entry"+chrom+"/"+entry+",will fail");
	    		exc1.printStackTrace();
	    		System.exit(0);
	    	}
	    }
	    return snps;
}

public RunR1(File f,File rscript,  File workingDir, int[] startend,  String nme) throws Exception{
		//	File f = zip;
	 
	        chrom = f.getName().substring(0,f.getName().lastIndexOf('.'));
	        String wn = workingDir.getName();
	        if(!wn.equals("."))  chrom = chrom.replaceAll(workingDir.getName(), "").replace('-', '_');
			this.nme = nme;
		this.r = new BufferedReader(new FileReader(rscript));
		re.assign("workingDir", workingDir.getAbsolutePath());
		System.err.println("opening "+f.getAbsolutePath());
		 zf = new ZipFile(f);
		  start = startend[0];
		  end = startend[1];
 	    
		 snps = getReader("SNPS");
		 {
 	    File snpsF = new File(f.getParentFile(),"snps.txt");
 	    if(!snpsF.exists() || snpsF.length()==0) {
 	    	snpsToDo = null;
 	    }
 	    else{
 	    	snpsToDo = new HashMap<String, String>();
 	    	BufferedReader br = new BufferedReader(new FileReader(snpsF));
 	    	String st = "";
 	    	while((st = br.readLine())!=null){
 	    		if(st.trim().startsWith("#")) continue;
 	    		String[] str = st.split("\\s+");
 	    		snpsToDo.put(str[0], str.length>1 ? str[1] : "none");
 	    	}
 	    }
		 }
		
 	  
 	    
		List<String > l2 =  read(zf, "Name");
		List<String> l20 =Arrays.asList(l2.get(0).split("\t")); 
		List<String> l21 =Arrays.asList(l2.get(1).split("\t")); 
		List<String> l22 =Arrays.asList(l2.get(2).split("\t"));
		fix(l20);fix(l21);fix(l22);
		//
		double[] maxg = new double[l20.size()];
		double[] ming = new double[l20.size()];
		if(l2.size()>3 ){
			this.def_values_for_geno =   l2.get(3).split("\t") ;
    		for(int k=0; k<maxg.length; k++){
	    		String str = def_values_for_geno[k];
		    	maxg[k] = str.indexOf(',')>=0 ? str.split(",").length : Integer.parseInt(str);
		    	ming[k] = 0;
		    }
			String[] spl = def_values_for_geno[0].split(",");
			for(int j=0; j<spl.length; j++){
				
				if(spl[j].length()>0 && Double.parseDouble(spl[j])>0) this.scaling =  Double.parseDouble(spl[j]);
			}
		}else{
			Arrays.fill(maxg,4);
			for(int k=0; k<maxg.length; k++){
				ming[k] =0;
				if(l20.get(k).indexOf("count")>=0) maxg[k] = 4;
				else if(l20.get(k).indexOf("state")>=0) maxg[k] = 2;
				else if(l20.get(k).indexOf("geno")>=0) maxg[k] = 2;
				else if(l20.get(k).indexOf("GType")>=0) maxg[k] = 2;
				else if(l20.get(k).indexOf("Log")>=0){
					ming[k] = -5;
					maxg[k] = 5;
				}
				else if(l20.get(k).indexOf("Freq")>=0) maxg[k] = 1;
				
			}
		}
		re.assign("maxg", maxg);
		re.assign("ming", ming);
		rs_ind = l21.indexOf("rsid");
		nosnp_ind = l21.indexOf("nosnps");
		if(rs_ind<0) rs_ind = l21.indexOf("snpid");
		if(rs_ind<0) rs_ind = l21.indexOf("id");
		if(rs_ind <0) throw new RuntimeException("could not find rs index");
		int nosamp = this.countLines(zf, "Samples");
		String[][] samples = new String[l22.size()][nosamp];
	 
		this.read(zf,"Samples", samples);
		assign(samples,"samples", l22.toArray(new String[0]),false, null);
		Set<String> resNames = new HashSet<String>(Arrays.asList("Samples:Name:SNPS"));
		
		
			geno_header = l20.toArray(new String[0]);
			re.assign("geno_header", geno_header);
			//this.next(current);
			this.geno = new String[l20.size()][nosamp];
		this.initialise();
		double[] spl  = re.eval("spl").asDoubleArray();
		expandData = re.eval("expandData").asBool().isTRUE() ? spl.length-1 : 1;
		inverse = re.eval("inverseRegress").asBool().isTRUE();
		
		this.maxIndiv =nosamp+10;//(int) eval("max_indiv").asDouble();
	
		if(expandData>1) this.weights = new double[l20.size()][geno[0].length*expandData];
		this.rescale = eval("rescale").asDouble();
		this.pleio = eval("multiPhen").asBool().isTRUE();
		this.multiGen = eval("multiGen").asBool().isTRUE();
		double[] inds = eval("datanme$ind").asDoubleArray();
		if(inds==null) throw new RuntimeException("inds is null - datanme$ind undefined in R");
		   this.colsToInclude = new int[inds.length];
			
if(inds.length==0){
	 throw new RuntimeException("no geno indices selected");
}
		Arrays.fill(colsToInclude,-1);
		for(int k=0; k<inds.length; k++){
			colsToInclude[k]=(int)inds[k]-1;
			System.err.println("cols to include "+inds[k]);
		}
		System.err.println("h");
	}
final boolean inverse;
final int expandData;
	
	
	

	protected RunR1(BufferedReader r2, String string, Integer start2, Integer end2,
		String chrom2, int expand, boolean inverse, int maxIndiv, boolean pleio, boolean multiGen, Map<String, String> snpsToDo) {
this.r = r2;
this.pleio = pleio;
this.nme = string;
this.multiGen = multiGen;
this.maxIndiv = maxIndiv;
start = start2;
		end = end2;
		chrom = chrom2;
		this.geno = null;
		this.inverse = inverse;
		this.expandData = expand;
		this.snpsToDo = snpsToDo;
}

	private void assign(String[][] samples, String string, String[] header, boolean transpose, int[] col_ids1) {
		StringBuffer sb = new StringBuffer();
		int[] col_ids = col_ids1;
		if(col_ids1==null){
			col_ids = new int[samples.length];
			for(int k=0; k<col_ids.length; k++){
				col_ids[k] = k;
			}
		}
		for(int k=0; k<col_ids.length; k++){
			String nme = "x_"+k;
			
			re.assign(nme, samples[col_ids[k]]);
			sb.append(nme+",");
			
		}
		//REXP r = eval(string);
		String str  = 
			transpose ? string+"= as.matrix(rbind("+sb.substring(0,sb.length()-1)+"))":
			string+"= as.matrix(cbind("+sb.substring(0,sb.length()-1)+"))";
		
		eval(str);
		
		if(header!=null){
			sb = new StringBuffer();
			for(int k=0; k<header.length; k++){
				sb.append("\""+header[k]+"\",");
			}
			str = sb.substring(0,sb.length()-1);
			str = "dimnames("+string+")[[2]] = c("+str+")";
			eval(str);
		}
		//double[][] gen1 = re.eval(string).asDoubleMatrix();
	}
	private void assign(double[][] samples, String string, String[] header, boolean transpose, int[] cols) {
		StringBuffer sb = new StringBuffer();
		for(int k=0; k<cols.length; k++){
		
			String nme = "x_"+k;
			
			re.assign(nme, samples[cols[k]]);
			sb.append(nme+",");
		
		}
		String str  = 
			transpose ? string+"= as.matrix(rbind("+sb.substring(0,sb.length()-1)+"))":
			string+"= as.matrix(cbind("+sb.substring(0,sb.length()-1)+"))";
		//String str  = string+"= as.matrix(cbind("+sb.substring(0,sb.length()-1)+"))";
		eval(str);
		
		if(header!=null){
			sb = new StringBuffer();
			for(int k=0; k<header.length; k++){
				sb.append("\""+header[k]+"\",");
			}
			str = sb.substring(0,sb.length()-1);
			str = "dimnames("+string+")[[2]] = c("+str+")";
			eval(str);
		}
		//double[][] gen1 = re.eval(string).asDoubleMatrix();
	}



	
	StringBuffer sb = new StringBuffer();
	
	int numtoclose=0;
	protected REXP eval(String str1) {
		if(str1.startsWith("#")) return null;
		String str = str1.replaceAll("datanme", this.nme);
		if(str1.indexOf('}')>=0){
			sb.append(str1);
			if(str1.indexOf('{')<0) numtoclose--;
			if(numtoclose==0){
				String str2 = sb.toString();
			    if(!str.startsWith("if(debug")) {re.eval(str2);
					if(LoopCall.print) System.err.println(">"+str2);
			    }
			  sb = new StringBuffer();
			}
			
		}
		else if(str.indexOf('{')>=0){
			sb.append(str);
			if(!str.endsWith("{")) sb.append(";");
			numtoclose++;
		}
		else if(numtoclose>0){
			if(str.indexOf('#')<0){
				sb.append(str+";");
			}
		}
		else{
			if(numtoclose!=0) throw new RuntimeException("!!");
			if(LoopCall.print) 
			  if(!str.startsWith("if(debug")) {
				  System.err.println(">"+str);
			     return re.eval(str);
			  }
		}
		return null;
		
	}

	private String readNext() throws Exception{
		String st = r.readLine().trim();;
		while(st.startsWith("#") && st.length()==0){
			st = r.readLine().trim();;
		}
		return st;
	}
	//final int thin;
	
	
	private static void assign(BufferedReader br) throws Exception{
		String str = "";
		for(int k=0; (str = br.readLine())!=null; k++){
			
			re.assign("x_"+k, str);
		}
	}
	//final String st1,st2,st3,st4;
	//List<String> indiv;
	
	
	public static void fix(List<String> l21) {
		// TODO Auto-generated method stub
		for(int k=0; k<l21.size(); k++){
			l21.set(k, l21.get(k).trim());
		}
	}

	private List<String> read(ZipFile zf,  String probe, int i) throws Exception {
		List<String> res = new ArrayList<String>();
		
		BufferedReader br = getReader(probe);
			//new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st.split("\t")[0]);
		}
		br.close();
		return res;
		
		
	}
	public  List<String> read(ZipFile zf,  String probe) throws Exception {
		List<String> res = new ArrayList<String>();
		
		BufferedReader br =getReader(probe);
			//new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		while((st =br.readLine() )!=null){
			res.add(st);
		}
		br.close();
		return res;
		
		
	}
	private double[] read(ZipFile zf,  String probe, int i, double[] res, int[] alias) throws Exception {
		//List<String> res = new ArrayList<String>();
		BufferedReader br =getReader(probe);
			//new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		for(int k=0; (st =br.readLine() )!=null;k++){
			res[alias==null ? k : alias[k]]=Double.parseDouble(st.split("\t")[i]);
		}
		br.close();
		return res;
		
	}
	private int countLines(ZipFile zf, String probe) throws Exception{
		BufferedReader br = getReader(probe);
			//new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		int i=0;
		int max_cols =0;
		for(;(st=br.readLine())!=null; i++){
			int len = st.split("\t").length;
			if(len>max_cols) max_cols = len;
		}
		br.close();
		return i;
	}
	private void read(ZipFile zf,  String probe,  String[][] res) throws Exception {
		//List<String> res = new ArrayList<String>();
		BufferedReader br =getReader(probe); 
			//new BufferedReader(new InputStreamReader(zf.getInputStream(zf.getEntry(probe))));
		String st = "";
		int k=0;
		for( k=0; k<res[0].length &&  (st =br.readLine() )!=null;k++){
			String[] str = st.split("\\s+");
			for(int j=0; j<res.length; j++){
				if(j<str.length)res[j][k] =str[j];
				else res[j][k] ="";
			}
		}
		if(k<res[0].length){
			throw new RuntimeException("!!");
		}
		br.close();
		
		
	}
	
}
