package data;

import org.rosuda.JRI.REXP;

public class RunRMerged extends RunR1 {
	RunR1[] rr;
	public RunRMerged(RunR1[] rr) throws Exception{
		super(rr[0].r,"merge",rr[0].start,rr[0].end,rr[0].chrom,rr[0].expandData, 
				rr[0].inverse,rr[0].maxIndiv,rr[0].pleio, rr[0].multiGen, rr[0].snpsToDo);
       this.rr = rr;
       this.initialise();
       this.nosnp_ind = rr[0].nosnp_ind;
     //  this.current = rr[0].current;
    //   this.snps = rr[0].snps;
		// TODO Auto-generated constructor stub
	}
	
	public String mergeString(RunR1[] rr,String type, int len,String suff){
		 StringBuffer st = new StringBuffer(type);
		 if(len<2) throw new RuntimeException("!!");
		 if(len==2){
			 st.append("=rbind(");
		 }else{
			 st.append("=abind(");
		 }
		 for(int k=0; k<rr.length; k++){
			 st.append(rr[k].nme+suff+"$"+type);
			 if(k<rr.length-1) st.append(",");
			 
		 }
		 if(len==2){
			 st.append(")");
		 }else{
			 st.append("along=3)");
		 }
		st.append("; dimnames("+type+")[[2]] = dimnames("+rr[0].nme+suff+"$"+type+")[[2]]");
	     return st.toString();
	}
	public void initialise() throws Exception{
		initialise(""); initialise("2");
	}
	public void initialise(String suff) throws Exception{
		String[] str1 = rr[0].eval("names(datanme)").asStringArray();
		StringBuffer sb = new StringBuffer("datanme = list(");
		for(int k=0; k<str1.length; k++){
			//if(rr[0].eval("is.matrix(datanme$"+str1[k]+")").asBool().isTRUE()){
			REXP dim = rr[0].eval("dim(datanme$"+str1[k]+")");
			
			if(dim!=null && dim.asIntArray()!=null){
				 eval(mergeString(rr,str1[k],dim.asIntArray().length,suff));
			}else{
				 eval(str1[k]+"="+rr[0].nme+suff+"$"+str1[k]);
			}
				 sb.append("\""+str1[k]+"\"="+str1[k]);
			if(k<str1.length-1) sb.append(",");
		}
		sb.append(")");
		/* eval(mergeString(rr,"pheno"));
	        eval(mergeString(rr,"incl"));
	        eval(mergeString(rr,"pheno_cov"));
	        eval(mergeString(rr,"caseInd"));
	        eval("ind = "+rr[0].nme+"$ind");
			String str = "datanme = list(\"pheno\"=pheno,\"incl\" = incl,\"pheno_cov\" = pheno_cov," +
					"\"offsets\" =NULL, \"geno\" = NULL,\"ind\" = ind)";*/
			eval(sb.toString());

		if(CHECK){
		 String[] fams = eval("families").asStringArray();
		 double[][] phenos = eval("datanme$pheno").asDoubleMatrix();
		 REXP inclM = eval("datanme$incl");
		 System.err.println("h");
		}
	}
	
	@Override 
	public boolean  nextLine() throws Exception{
		boolean res = rr[0].nextLine();
		this.curr_snp =  rr[0].curr_snp;
		if(nosnp_ind >=0){
		StringBuffer no_snps = new StringBuffer(rr[0].nosnp_ind<0 ? "1" : rr[0].curr_snp[rr[0].nosnp_ind]);
		for(int k=1; k<rr.length; k++){
			boolean res1 = rr[k].nextLine();
			if(res1!=res) throw new RuntimeException("!!");
			no_snps.append(","+(rr[k].nosnp_ind<0 ? "1" : rr[k].curr_snp[rr[k].nosnp_ind]));
			if(!rr[k].curr_snp[0].equals(rr[0].curr_snp[0])) throw new RuntimeException("!!");
		}
		this.curr_snp[this.nosnp_ind] = no_snps.toString();
		}
		return res;
	}
	@Override
	public  void  processLine(String snp_id) throws Exception{
		for(int i=0; i<rr.length; i++){
		 rr[i].processLine(snp_id);
		}
		eval(mergeString(rr,"geno",2,""));
		eval("datanme$geno = geno");
		eval(mergeString(rr,"weights",2,"2"));
		eval("datanme2$weights = weights");
	}
}
