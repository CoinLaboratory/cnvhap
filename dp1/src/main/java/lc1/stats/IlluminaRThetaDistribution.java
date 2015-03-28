package lc1.stats;

import lc1.util.Constants;

public class IlluminaRThetaDistribution extends IlluminaDistribution {

	public IlluminaRThetaDistribution(short data_index) {
		super(data_index);
		// TODO Auto-generated constructor stub
	}
	public IlluminaRThetaDistribution(IlluminaRThetaDistribution r){
		super(r);
	}
	
	@Override
    public PseudoDistribution swtchAlleles() {
       transformBack();
       super.swtchAlleles();
		transformToThetaR();
        return this;
    }
	public void transformToThetaR() {
		int index = this.data_index;
		double aCount = this.r;
		double bCount = this.b;
//		System.err.println("cnts "+aCount+" "+bCount);
		if(aCount+bCount<=0){
			throw new RuntimeException("!!");
		}
		this.r =aCount+bCount;
		if(Constants.transformRToLogR) r = Math.log((r)/2.0);
	
		this.b = Math.atan2(bCount,aCount)*(Constants.oneeighty/Math.PI);
		//if(b<-1.0) throw new RuntimeException("!!");
//		System.err.println(r+" "+b+" "+aCount+" "+bCount+" ");
		if(r>Constants.maxR(index)) Constants.maxR[index] = r;
		if(r<Constants.minR(index)) Constants.minR[index] = r;
		if(b>Constants.maxB(index)) Constants.maxB[index] = b;
		if(b<Constants.minB(index)){
			Constants.minB[index] = b;
		}
	}
	public void transformBack(){
		double x = Math.tan(b *(Math.PI/Constants.oneeighty));
		double r1 = 	Constants.transformRToLogR ? Math.exp(r)*2.0 : r;
		b = r1/(x+1);
		r =x*b;
	}
	
	
	
	
	
}
