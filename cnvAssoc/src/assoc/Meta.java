package assoc;
import java.util.Arrays;

import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;



public class Meta {
	int len;
	ChiSquaredDistributionImpl chisq = new ChiSquaredDistributionImpl(1);

	Meta(int noStudy) {
		this.len = noStudy;
		this.variance = new double[len];
		this.wfixed = new double[len];
		this.wbeta = new double[len];
		this.include = new boolean[len];
		this.Lambda = new double[len];
		Arrays.fill(Lambda, 1);
		wrandom = new double[len];
		wbetarandom = new double[len];
	}
	double pFisher, tau2, zscorerandom,betaCombined, seCombined;
	double[] wrandom, wbetarandom;
	double pRandom, pCombined;
	double pFix, betarandom, sebetarandom;
	double pQ, betafixed, sebetafixed;
	double[] variance;
	double[] Lambda;
	double sumwbeta,sumw,sumwsquare, Q, countstudy,sumwbetarandom,sumwrandom, chisqsum;
	double[] wfixed, wbeta;
	boolean[] include;
	double zscorefix;
	public void calc(double[] beta, double[] se, double[] p) throws Exception{
	this.sumwbeta = 0;
	this.sumw = 0;
	this.sumwsquare = 0;
	this.Q = 0;
	this.countstudy = 0;
	this.sumwbetarandom = 0;
	this.sumwrandom = 0;
	//infomax = 0;
	this.chisqsum = 0;

	// holder=OR[1];
	// holder=lowerCI.length;
	// for(int i=0; i < OR.length; i++){

	// holder= this.OR[i];
	// }

	for (int i = 0; i < se.length; i++) {

		variance[i] = Math.pow(se[i], 2);

		if (Double.isNaN(se[i])
				|| se[i] < 0
			//	|| this.info[i] < this.info_thresh
				//|| (freq_b[i] != null && (freq_b[i] < maf_thresh || freq_b[i] > maf_thresh1))
				) {
			include[i] = false;
			wfixed[i] = 0;
			wbeta[i] = 0;
		} else {
		//	if (Double.isNaN(this.info[i])) {
		//		throw new RuntimeException("info should be defined here");
		//	}
			if (Double.isNaN(beta[i])) {
				throw new RuntimeException("beta should be defined here");
			}
//			if (this.info[i] > infomax)
//				infomax = this.info[i];
			chisqsum += (-2 * Math.log(p[i])) / Lambda[i];
			wfixed[i] = 1 / variance[i];
			wbeta[i] = wfixed[i] * beta[i];
			countstudy = countstudy + 1;
			include[i] = true;
		}
		sumwbeta = sumwbeta + wbeta[i];
		sumw = sumw + wfixed[i];
		sumwsquare = sumwsquare + Math.pow(wfixed[i], 2);
	}
	if (countstudy == 0) {
		
		pQ = Double.NaN;
		this.betafixed = Double.NaN;
		this.sebetafixed = Double.NaN;
		
		pFix = Double.NaN;
		betarandom = Double.NaN;
		this.sebetarandom = Double.NaN;
		
		pRandom = Double.NaN;
		this.pCombined = Double.NaN;
	//	return freq_b;
	}
	betafixed = sumwbeta / sumw;

	if (sumw != 0) {
		sebetafixed = Math.pow(1 / sumw, 0.5);
		double zscorefix;
		zscorefix = Math.abs(betafixed / sebetafixed);
		
			//FIX THIS!!!!
			chisq.setDegreesOfFreedom(1.0);
			pFix = 	chisq.cumulativeProbability(Math.pow(zscorefix, 2));
		// 2 * (1 - pal.statistics.NormalDistribution.cdf(zscorefix, 0,
		// 1));
	} else {
		pFix = Double.NaN;
		sebetafixed = Double.NaN;
		
		zscorefix = Double.NaN;
	}
	double[] Qi = new double[beta.length];
	for (int i = 0; i < beta.length; i++) {
		if (!Double.isNaN(se[i])) {
			Qi[i] = wfixed[i] * Math.pow(beta[i] - betafixed, 2);
			Q = Q + Qi[i];
		} else
			Qi[i] = 0;
	}
	
	if (sumw != 0) {
chisq.setDegreesOfFreedom(countstudy);
		pQ = chisq.cumulativeProbability(Q);
		chisq.setDegreesOfFreedom(2*countstudy);
		pFisher = chisq.cumulativeProbability(this.chisqsum
				);
	} else {
		pQ = Double.NaN;
		pFisher = Double.NaN;
	}
	if (Q <= countstudy - 1) {
		tau2 = 0;
	} else {
		tau2 = (Q - (countstudy - 1)) / (sumw - (sumwsquare / sumw));
	}

	for (int i = 0; i < se.length; i++) {
		if (Double.isNaN(se[i])) {
			wrandom[i] = 0;
			wbetarandom[i] = 0;
		} else {
			wrandom[i] = 1 / (variance[i] + tau2);
			wbetarandom[i] = wrandom[i] * beta[i];
		}
		sumwbetarandom = sumwbetarandom + wbetarandom[i];
		sumwrandom = sumwrandom + wrandom[i];
	}

	betarandom = sumwbetarandom / sumwrandom;
	sebetarandom = Math.pow(1 / sumwrandom, 0.5);
	zscorerandom = Math.abs(betarandom / sebetarandom);
	if (sumw != 0) {
		chisq.setDegreesOfFreedom(1.0);
		pRandom =chisq.cumulativeProbability(Math.pow(zscorerandom, 2));
			
		// 2 * (1 - pal.statistics.NormalDistribution.cdf(
		// zscorerandom, 0, 1));
	} else {
		pRandom = Double.NaN;
	}
	
	if (Double.isNaN(pQ)) {
		pCombined = Double.NaN;
		betaCombined = Double.NaN;
		seCombined = Double.NaN;
	} else if (pQ < 0.01) {
		pCombined = pRandom;
		betaCombined = this.betarandom;
		seCombined = this.sebetarandom;
	} else {
		pCombined = this.pFix;
		betaCombined = this.betafixed;
		seCombined = this.sebetarandom;
	}
	}
}
