package lc1.stats;

import java.util.Arrays;
import java.util.List;

import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.Constants;

public class HWECalc {

	EmissionStateSpace emstsp;
	
	List<String> l;
	DataCollection dc;
	int[] no_cop;
	int[] numA;
	//int[] numB;
	 
	double total;
	 
	public HWECalc(EmissionStateSpace emstsp){
		this.emstsp = emstsp;
		no_cop = new int[emstsp.defaultList.size()];
		numA = new int[emstsp.defaultList.size()];
		for(int k=0; k<no_cop.length; k++){
			ComparableArray compa = (ComparableArray) emstsp.get(k);
			no_cop[k] = compa.noCopies(true);
			numA[k] = no_cop[k] - (int) compa.noB();
		//	numB[k] = (int) compa.noB();
		}
	
	}
	
	public void setList(List<String> l, DataCollection dc){
		this.l = l;
		this.dc = dc;
		total = l.size();
	}
	
	double[] count = new double[5];  //counts of cn from 0 to 4
	double[][] count1 = new double[5][5]; //A allele count within each of 4 cn states
	
	double[] exp = new double[5];  //counts of cn from 0 to 4
	double[][] exp1 = new double[5][5]; //A allele count within each of 4 cn states
	
	double p,q,r;  //allele freq 0,1,2
	double[] pp = new double[5];
	
	//double[] qq = new double[4];
	public void count(int i){
		Arrays.fill(count,0);
		Arrays.fill(exp,0);
		for(int k=0; k<count1.length; k++){
			Arrays.fill(count1[k], 0);
			Arrays.fill(exp1[k], 0.0);
		}
		p=0;
		q =0;
		r=0;
		Arrays.fill(pp,0.0);
	//	Arrays.fill(qq,0.0);
		for(int j=0; j<l.size(); j++){
			PseudoDistribution dist =((HaplotypeEmissionState) dc.dataL.get(l.get(j))).emissions[i];
			Integer k=dist.fixedInteger();
			if(k!=null){
				addCount(k, 1.0);
			}
			else{
				double[] probs = dist.probs();
				for(int k1=0; k1<probs.length; k1++){
					addCount(k1, probs[k1]);
				}
			}
		}
		double tot = p+q+r;
		p = p/tot;
		q = q/tot;
		r = r/tot;
		for(int k=1; k<pp.length; k++){
			pp[k] = count[k]==0 ? 0: pp[k]/ (count[k]*(double)no_cop[k]);
		}
		
			exp[0] = Math.pow(p,2)*total;
			exp[1] = 2*p*q*total;
			exp[2] = Math.pow(q,2)*total+2*p*r*total;
			exp[3] = 2*q*r*total;
			exp[4] = Math.pow(r, 2)*total;
	/*	double sum = Constants.sum(exp);
		if(Math.abs(this.total - sum) >0.1) {
			throw new RuntimeException("!!");
		}*/
			exp1[0][0] = 1.0 *count[0];
			exp1[1][0] = (1-pp[1])*count[1];
			exp1[1][1] = pp[1]*count[1];
			exp1[2][0] = Math.pow(1-pp[2],2)*count[2];
			exp1[2][1] = 2*(1-pp[2])*pp[2]*count[2];
			exp1[2][2] = Math.pow(pp[2],2)*count[2];
			exp1[3][0] = Math.pow(1-pp[3],3)*count[3];
			exp1[3][1] = 3*pp[3]*Math.pow(1-pp[3],2)*count[3];
			exp1[3][2] = 3*(1-pp[3])*Math.pow(pp[3],2)*count[3];
			exp1[3][3] = Math.pow(pp[3],3)*count[3];
			exp1[4][0] = Math.pow(1-pp[4],4)*count[4];
			exp1[4][1] = 4*pp[4]*Math.pow(1-pp[4],3)*count[4];
			exp1[4][2] = 6*Math.pow(1-pp[4],2)*Math.pow(pp[4],2)*count[4];
			exp1[4][3] = 4*(1-pp[4])*Math.pow(pp[4],3)*count[4];
			exp1[4][4] = Math.pow(pp[4],4)*count[4];
		}
	
	private void addCount(int k,  double d ) {
		count[no_cop[k]]+=d;
		count1[no_cop[k]][numA[k]]+=d;
		if(no_cop[k]==0){
			p+=2*d;
		}
		else if(no_cop[k]==1){
			p+=d;
			q+=d;
			
		}
		else if(no_cop[k]==2){
			q+=2*d;
		}
		else if(no_cop[k]==3){
			q+=d;
			r+=d;
		}
		else if(no_cop[k]==4){
			r+=2*d;
		}
		pp[no_cop[k]]+=(numA[k])*d;
	}
ChiSq chisq = new ChiSq();

public double sig(){
	return sig(count, exp, 2);
}
public double sig(int i){
	return sig(count1[i], exp1[i],1);
}
	public double sig(double[] count, double[] exp, int deg){
		double sum=0;
		double diff =Math.abs(Constants.sum(count)-Constants.sum(exp)); 
		if(diff>0.01){
			throw new RuntimeException("!!");
		}
		for(int i=0; i<count.length; i++){
			if(exp[i]>0){
				sum+=Math.pow(count[i]-exp[i], 2) /exp[i];
			}
		}
		double p =  chisq.chi2prob(deg, sum);
		return p;
	}
	
}
