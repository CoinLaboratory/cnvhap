import lc1.stats.ChiSquareDistribution;


public class Test {
public static void main(String[] args){
	
	try{
		double inp = 300;
		double d= Math.log10(1.0 - ChiSquareDistribution.cdf(inp,1));
		System.err.println(d);
		
		double d1= ChiSquareDistribution.cdfLog(inp,1);
		System.err.println(Math.pow(10,d1));
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
}
