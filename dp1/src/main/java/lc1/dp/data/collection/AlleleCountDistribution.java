package lc1.dp.data.collection;

import cern.jet.stat.Probability;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.PseudoDistribution;
import JSci.maths.statistics.*;

public class AlleleCountDistribution extends PseudoDistribution {
	int a_count,b_count;
	int b_count_low,a_count_low;
	String alleles;
	String quality;
	double bfreq;

	public AlleleCountDistribution(String alleles,String quality, short data_index) {

		a_count = count(alleles,'A');
		b_count = count(alleles,'B');
		this.alleles = alleles;
		this.quality = quality;
		setDataIndex(data_index);


	}
	
	public AlleleCountDistribution(String alleles,String quality,short data_index, double bfreq) {

		a_count = quality.length()==alleles.length()?countHighQuality(alleles,'A',quality):0;
		b_count = quality.length()==alleles.length()?countHighQuality(alleles,'B',quality):0;
		b_count_low = count(alleles,'B');
		a_count_low = count(alleles,'A');
		b_count=b_count_low;
		a_count=a_count_low;
		setDataIndex(data_index);
		this.bfreq = bfreq;
		this.alleles = alleles;
		this.quality = quality;


	}


	protected static int count(String st, char character) {
		int cnt =0;
		for(int i=0; i<st.length(); i++){
			if( st.charAt(i)==character) cnt++;
		}
		return cnt;
	}
	
	protected static int countHighQuality(String st, char character, String quality) {
		int cnt =0;
		for(int i=0; i<st.length(); i++){
			if( st.charAt(i)==character && ((int)quality.charAt(i)-33)>=30) cnt++;
		}
		return cnt;
	}

	@Override
	public void addCount(int objIndex, double value) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public double[] calcDistribution(double[] distribution,
			EmissionStateSpace emStSp, int pos) {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public PseudoDistribution clone() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public PseudoDistribution clone(double swtch) {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public double[] counts() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public Integer fixedInteger() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public int getMax() {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public String getPrintString() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public void initialise() {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public boolean isMeasured() {
		if(true) throw new RuntimeException("!!");
		return false;
	}

	@Override
	public double logProb() {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public double[] probs() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public double probs(int objI) {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public double scoreB(int j, int i) {
		if(true) throw new RuntimeException("!!");
		return 0;
	}
	//final EmissionStateSpace emstsp;
//	@Override
	public double testScoreBR(EmissionStateSpace emStSp, int j, int i) {
		
		int bcount = emStSp.getBCount(j);
		int tot = emStSp.getCN(j);
		if(tot==0) return 0.0;
		int acount = tot - bcount;
		double binprob = (double)bcount/(double)tot;
//		double binprob=bfreq;
		double probBB = -1;
		double probAB = -1;
		double probAA = -1;
		
		
		if(tot==0 || a_count+b_count==0 || quality.length()==0) 
			return 0.0;
		else
		{
			double[] numQuality = phredToProb(quality);
			if (acount!=0 && bcount!=0) //AB
			{
				probBB = this.rightProb('B', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
				probAA = this.rightProb('A',quality,alleles);
			}
			else if(acount==0){ //BB
				probBB = this.rightProb('B', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
				probAA = this.wrongProb('B', quality, alleles);
			}
			else if(bcount==0){//AA
				probBB = this.wrongProb('A', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
				probAA = this.rightProb('A', quality, alleles);
			}
		}
		double[] normalizedProb=normalize(new double[]{probBB,probAB,probAA});
//		int bcount = emStSp.getBCount(j);
//		int tot = emStSp.getCN(j);
//		if(tot==0) return 0.0;
//		double binprob = (double)bcount/(double)tot;

		
//		double d =  Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
		double d =  normalizedProb[1];
		return d;
	}

//	@Override
	public double oldScoreBR(EmissionStateSpace emStSp, int j, int i) {
		
		int bcount = emStSp.getBCount(j);
		int tot = emStSp.getCN(j);
		if(tot==0) return 0.0;
		int acount = tot - bcount;
//		double binprob = (double)bcount/(double)tot;
		double binprob=bfreq;
		double probBB = -1;
		double probAB = -1;
		double probAA = -1;
		
		double[] normalizedProb = new double[] {0,0,0};
		if(a_count+b_count==0 || quality.length()==0) 
			return 1.0;
		else if (tot==0)
			return 0.0;
		else
		{
			double[] numQuality = phredToProb(quality);
			if (a_count!=0 && b_count!=0) //AB
			{
				probBB = this.rightProb('B', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob);
				probAA = this.rightProb('A',quality,alleles);
			}
			else if(a_count==0){ //BB
				probBB = this.rightProb('B', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob);
				probAA = this.wrongProb('B', quality, alleles);
			}
			else if(b_count==0){//AA
				probBB = this.wrongProb('A', quality, alleles);
				probAB = Probability.binomial(this.b_count, this.a_count+this.b_count, binprob);
				probAA = this.rightProb('A', quality, alleles);
			}
			normalizedProb=normalize(new double[]{probBB,probAB,probAA});
		}
		
//		int bcount = emStSp.getBCount(j);
//		int tot = emStSp.getCN(j);
//		if(tot==0) return 0.0;
//		double binprob = (double)bcount/(double)tot;

		
//		double d =  Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
		double d =  normalizedProb[1];
		return d;
	}
	
//	public double[] getGenoProb(String[] geno, int stspSize){
//		double[] probs=new double[stspSize];
//		double[] temp=new double[stspSize];
//		int[] countAB=this.countAB(geno[0]);
//		if(countAB[0]==0 && countAB[1]==0) throw new RuntimeException("no allele data");
//		if(countAB[0]!=0 && countAB[1]!=0) { //AB
//			temp[0]=this.homoProb('B', geno[1], geno[0]);
//			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5)-Probability.binomial(countAB[0]-1, countAB[0]+countAB[1], 0.5);
//			temp[2]=this.homoProb('A', geno[1], geno[0]);
//		}
//		else if(countAB[0]==0){ //BB
//			temp[0]=this.homoProb('B', geno[1], geno[0]);
//			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5);
//			temp[2]=this.aveProb(geno[1]);    			
//		}
//		else if(countAB[1]==0){//AA
//			temp[0]=this.aveProb(geno[1]); 
//			temp[1]=Probability.binomial(countAB[0], countAB[0]+countAB[1], 0.5)-Probability.binomial(countAB[0]-1, countAB[0]+countAB[1], 0.5);
//			temp[2]=this.homoProb('A', geno[1], geno[0]);    			
//		}
//		else throw new RuntimeException("!!");
//		probs=this.normalize(temp);
//
//		return probs;
//	}

//	@Override
	public double scoreBR(EmissionStateSpace emStSp, int j, int i) {
		if(a_count+b_count==0) 
			return 1.0;
		int bcount = emStSp.getBCount(j);
		int tot = emStSp.getCN(j);
		double bHapProb = (double)bcount/(double)tot;
		if(tot==0) 
		{
//			bHapProb = 1e-3;
			return 0.0;
		}
		
//		double binprob = (double)bcount/(double)tot;
//		double binprob = this.bfreq * this.rightProb('B', quality, alleles)* bHapProb;
//		int pseudocount = 0;
//		int multiplier = 15 -tot;
//		BetaDistribution currBeta = new BetaDistribution((bcount+1+pseudocount)*multiplier,(tot-bcount+1+pseudocount)*multiplier);
//		double normalize = 0;
//		for (float ii=0.01f;ii<1;ii+=0.01)
//		{
//			double prob = currBeta.probability(ii);
//			if (prob>normalize)
//				normalize=prob;
//		}
//				
//					
//		double binprob =  currBeta.probability(this.bfreq);///currBeta.probability((bcount+1+0.01-1)/(bcount+1+0.0000001+tot-bcount+1-2));// Probability.beta(bcount+1, tot-bcount+1,this.bfreq);//-Probability.beta(bcount+1, tot-bcount+1,this.bfreq-1e-4);//(tot*0.05+0.04)*(bcount+0.05);
////		double normalize = 1;
//		normalize*=1.1;
//		
//		if (bcount==0)
//		{
//			if (tot!=0 || (tot==0 && pseudocount!=0))
//				normalize = pseudocount==0?currBeta.probability(0.0001):(currBeta.probability(0.5)*1.5);
//		}else
//		{
//			if (bcount/tot==1 || pseudocount!=0)
//				normalize = pseudocount==0?currBeta.probability(0.9999):(currBeta.probability(0.5)*1.5);
//			else
//				normalize = currBeta.probability(((double)bcount+pseudocount)/((double) tot+pseudocount));
//		}
//		binprob/=normalize;
//		double binprob=1;
//		if (tot==0 || tot==1)
//			binprob = this.bfreq;
//		else if (tot==2)
//			binprob = 1-Math.pow(1-this.bfreq,2);
//		else
//			binprob = 1-Math.pow(1-this.bfreq, 3);
		double binprob=bHapProb;
		double d =  (new BinomialDistribution(this.a_count+this.b_count,binprob)).probability(this.b_count);// Probability.binomial(this.b_count, this.a_count+this.b_count, binprob)-(this.b_count>0?Probability.binomial(this.b_count-1, this.a_count+this.b_count, binprob):0);
		return d;
		//-Probability.binomial(countAB[0]-1, countAB[0]+countAB[1], binprob);

	}
	



	public int[] countAB(String ab){
		int[] count=new int[]{0,0};
		for(int i=0; i<ab.length(); i++){
			if(ab.charAt(i)=='A') count[0]++;
			else count[1]++;
		}
		if(count[0]==0 && count[1]==0) throw new RuntimeException("no alleles");
		return count;
	}

	public double homoProb(char errorAllele, String quality, String base){
		double prob=1.0;
		for(int i=0; i<base.length(); i++){
			if(base.charAt(i)==errorAllele) prob*=(Math.pow(10, (-((int)quality.charAt(i)-33)/10)));
		}
		return prob;
	}
	
	public double rightProb(char errorAllele, String quality, String base){
		double prob=1.0;
		for(int i=0; i<base.length(); i++){
			if(base.charAt(i)==errorAllele) prob*=(1-Math.pow(10, (-((int)quality.charAt(i)-33)/10)));
		}
		return prob;
	}
	
	public double wrongProb(char errorAllele, String quality, String base){
		double prob=0.0;
		for(int i=0; i<base.length(); i++){
			if(base.charAt(i)==errorAllele) prob+=(Math.pow(10, (-((int)quality.charAt(i)-33)/10)));
		}
		return prob/base.length();
	}
	
	
	
	public double[] phredToProb(String quality)
	{
		double[] prob = new double[quality.length()];
		for (int i=0;i<quality.length();i++)
			prob[i] = Math.pow(10, (-((int)quality.charAt(i)-33)/10));
		return prob;
	}

	public double aveProb(String quality){
		double prob=0.0;
		for(int i=0; i<quality.length(); i++){
			
			prob+=(1-Math.pow(10, (-((int)quality.charAt(i)-33)/10)));
		}
		return prob/(double)quality.length();
	}

	public double[] normalize(double[] val){
		double[] nor=new double[val.length];
		double sum=0.0;
		for(int i=0; i<val.length; i++){
			sum+=val[i];
		}
		for(int i=0; i<nor.length; i++){
			nor[i]=val[i]/sum;
		}
		return nor;
	}



	@Override
	public void setCounts(int i1, double cnt) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public void setFixedIndex(int k) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public void setProb(double[] prob) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public void setProbs(int to, double d) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public double sum() {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public PseudoDistribution swtchAlleles() {
		if(true) throw new RuntimeException("!!");
		return null;
	}

	@Override
	public double totalCount() {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public void transfer(double pseudoC1) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public void validate() {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public void addCounts(ProbabilityDistribution probabilityDistribution) {
		if(true) throw new RuntimeException("!!");

	}

	@Override
	public double sample() {
		if(true) throw new RuntimeException("!!");
		return 0;
	}

	@Override
	public void setParamsAsAverageOf(ProbabilityDistribution[] tmp) {
		if(true) throw new RuntimeException("!!");

	}

}
