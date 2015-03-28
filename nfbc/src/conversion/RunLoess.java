package conversion;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.descriptive.rank.Median;

public class RunLoess {
Median perc = new Median();

public static void main(String[] args){
	try{
		File dir = new File(".");
		RunLoess loess = new RunLoess(dir, 0.35);
		loess.run();
	}catch(Exception exc){
		exc.printStackTrace();
	}
}
static double lowBound = -0.3;
static double highBound = 0.3;
	File dir;
	double mult = 1000;
	double charLength = (700*1000) / mult;
	double varthresh;
	double[] loess;
	double[] y;
	double[] x,w;
	int cnt;
	File[] in;
	LoessInterpolator lo;
	//ZipFile zf;
	double span;
	int length;
	public RunLoess(File dir, double var_thresh) throws Exception{
		this.in = dir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
			return !pathname.getName().endsWith(".out") && !pathname.getName().startsWith("loess");
			}
		}
		);
		this.dir = dir;
		Arrays.sort(in);
		//System.err.println(dir.getAbsolutePath()+" "+Arrays.asList(in));
		var = new double[in.length];
		median = new double[in.length];
		median_loess = new double[in.length];
		readSNPS();
		this.var_thresh = var_thresh;
	}
	double var_thresh;
	
	//static class MedianCalc{
	//	List<Double> l = new ArrayList<Double>();
	//}
	
	double[] var;
	double[] median;
	double[] median_loess;
	double overallMedian;
	
	public void run() throws Exception{
		for(int k=0; k<in.length; k++){
			System.err.println("read "+k);
			read(k);
		}
		overallMedian = perc.evaluate(median);
		System.err.println("overall med "+overallMedian);
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "loess"))));
		PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, "summary.out"))));
		for(int k=0; k<loess.length; k++){
			pw.println(loess[k]/(double) numbercounted);
		}
		for(int k=0; k<in.length; k++){
			pw1.print( in[k].getName()+(k==in.length-1 ? "\n" : "\t"));
		}
		pw1.println("variance");
		for(int k=0; k<in.length; k++){
			pw1.print(String.format("%5.3g", var[k])+(k==in.length-1 ? "\n" : "\t"));
		}
		pw1.println("median");
		for(int k=0; k<in.length; k++){
			pw1.print(String.format("%5.3g", median[k])+(k==in.length-1 ? "\n" : "\t"));
		}
		pw1.println("overall median\t"+String.format("%5.3g",overallMedian));
		pw1.println("median_loess");
		for(int k=0; k<in.length; k++){
			pw1.print(String.format("%5.3g", median_loess[k])+(k==in.length-1 ? "\n" : "\t"));
		}
	
		pw.close();
		pw1.close();
	}
	
	
	int offset;
	
	public void readSNPS() throws Exception{
		BufferedReader br = Utils.getBufferedReader(in[0]);
		String nxt =br.readLine();
		double sum=0;
		double cnt_ =0;
		double x0 = Double.parseDouble(nxt.split("\\s+")[0])/mult;
		double xn=x0;
		int k=0;
		String str =nxt;
		for( k=0; nxt!=null; k++){
			 nxt = br.readLine();
			if(nxt==null){
				String[] str_ = str.split("\\s+");
				xn = Double.parseDouble(str_[0])/mult;
			}
			str = nxt;
			//sum+=y[k];
		}
			double width = xn-x0;
			
			double len = k;
			double minSpan =10.0/len;
			if(charLength >width) charLength = width;
			span =Math.max( charLength/(width),minSpan); //minimum of 10 probes
		
			double spacing = width/len;
			int noSnps = (int) Math.round(span*len);
			System.err.println(width+" "+len+" "+spacing+" "+noSnps);
			offset = noSnps;
			length = k;
			x = new double[k+2*offset];
			y = new double[k+2*offset];//7000.0/mult
			w = new double[k+2*offset];
			for(int j=0; j<offset; j++){
				x[j] = x0-spacing *(offset-j);
			}
			for(int j=0; j<offset; j++){
				x[offset+length+j] = xn+spacing *(j+1);
			}
			
			width = x[x.length-1] - x[0];
			span =Math.max( charLength/(width),minSpan);
			span = Math.min(span, 1.0);
			if(width==0){
				noSnps = 1;
				span = 1.0;
				spacing = 3000;
			}
			
			loess = new double[k];
			
			lo = new LoessInterpolator(span,4);
			
	}
	
	int numbercounted=0;
	
	public void read(int i) throws Exception{
		BufferedReader br = Utils.getBufferedReader(in[i]);
		String str ="";
		double sum=0;
		double cnt_ =0;
		Arrays.fill(w, 1.0);
		Arrays.fill(y, 0);
		for(int k=0; ((str = br.readLine())!=null); k++){
			int k1 = k+offset;
			String[] str_ = str.split("\\s+");
			double y_ = Double.parseDouble(str_[1]);
			
			if(i==0){
				x[k1] = Double.parseDouble(str_[0])/mult;
				if(k1>0 && x[k1]-x[k1-1]<1e-5){
					x[k1] =x[k1-1]+1e-5;
				}
			}
			if(!Double.isNaN(y_)){
				sum+=y_;
				cnt_++;
				y[k1] = y_;
				
			}
			else{
				y[k1]=0;
				w[k1]=0;
			}
			//sum+=y[k];
		}
		median[i] = perc.evaluate(y, offset, length);
		double med = median[i];
		
		
		double mean = sum/cnt_;
		sum=0;
		double sum1=0;
		double cnt1 =0;
		for(int k=0; k<length; k++){
			int k1 =k+offset;
			double y_ = y[k1];
			if(!Double.isNaN(y_)){
				double c = Math.pow(y_ - mean,2);
				sum+=c;
			//	y[k1] = y_-med;
				if(y_-med > highBound || y_-med<lowBound){
					w[k1] =0;
				}
				else{
					sum1+=c;
					cnt1++;
				}
			}
		}
		var[i] = sum/cnt_;
		double var1 = sum1/cnt1;
		NormalDistributionImpl norm = new NormalDistributionImpl(median[i],Math.sqrt(var1));
		try{
		for(int k1=0; k1<offset; k1++){
			y[k1] = norm.inverseCumulativeProbability(Math.random());
		}
		
		for(int k1=offset+length; k1<2*offset+length; k1++){
			y[k1] = norm.inverseCumulativeProbability(Math.random());
		}
		}catch(Exception exc){
			System.err.println("prob with "+var1);
		}
		double[] smooth  =  lo.smooth(x, y,w);
		median_loess[i] = perc.evaluate(smooth, offset, length);
		if(var[i] < var_thresh){
			double medloess = median_loess[i];
			for(int k=0; k<length; k++){
				int k1 =k+offset;
				loess[k]+=(smooth[k1] - medloess);
			}
			numbercounted++;
		}
		
		
		br.close();
	}
	
}
