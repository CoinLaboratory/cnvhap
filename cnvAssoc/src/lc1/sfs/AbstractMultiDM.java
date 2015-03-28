/**
 * 
 */
package lc1.sfs;


public abstract class AbstractMultiDM{
	
	public static class Entry{
		double count=0;  //no of times cell used;
		double prob; //prior prob cell used;
		double bf=0;  //P(D | p1..pn) * prob 
		
		public String toString(){
			return String.format("%5.3g %5.3g %5.3g", new Object[] {prob, count, bf}).replaceAll(" ", "_");
		}
		
		public void addCount(double total){
			count+=bf/total;
		}
		public double calcBF(double lh){ //calculate the bayes factor
			this.bf = lh * prob;
			return bf;
		}
		
		double[] vals; //evaluation of set of functions 
		public void initialise(){
			this.count=0;
		}
		public Entry(double d){
			this.prob = d;
		}
		public void transfer(double total){
			prob = count/total;
			//if(Double.isNaN(prob)){
			//	throw new RuntimeException("!!");
			//}
			count=0;
			
		}
		
	}
	
	 int dim;
	
	 
	 
	 public void validate(){
		 try{
			 double s = this.sum();
			 if(Math.abs(s-1.0)>1e-5){
				 throw new RuntimeException("sum not equal to one: "+s);
			 }
		 }catch(Exception exc){
			 exc.printStackTrace();
		 }
	 }
	 
	abstract Entry get(int[] inds);
	
	public double[] getV(int[] inds){
		return get(inds).vals;
	}
	
	public double getBF(int [] inds, double lh){
		return this.get(inds).calcBF(lh);
	}
	public void addCounts(int[] inds, double sum){
		this.get(inds).addCount(sum);
	}
	
	public void calcLHWeightedFuncts(int[] inds, double[] lh, ProbMap[] probMaps){
		Entry ent = this.get(inds);
		double pr = ent.bf;
		double[] vals = ent.vals;
		if(pr<-1e-5){
			throw new RuntimeException("!!");
		}
		for(int j=0; j<vals.length; j++){
			lh[j]+=vals[j]*pr;
			probMaps[j].add(vals[j], pr);
		}
	}
	
	public  abstract void mult(double p);
	public abstract double transfer();
	public abstract double sum();
	abstract void set(int[] inds, double[] v);
	
	static class OneDMultiDM extends AbstractMultiDM{
	    Entry[] v;
		OneDMultiDM(int[] len, double p){
			this.dim = 1;
			v = new Entry[len[0]];
			double p1 = p/(double)v.length;
			for(int k=0; k<v.length; k++){
				v[k] = new Entry(p1);
			}
		}
		@Override
		Entry get(int[] inds) {
			// TODO Auto-generated method stub
			return v[inds[0]];
		}

		@Override
		void set(int[] inds, double[] v) {
			// TODO Auto-generated method stub
			this.v[inds[0]].vals = v;
		}
		@Override
		public double transfer() {
			double sum = 0;
			for(int k=0; k<v.length; k++){
				sum+=v[k].count;
			}
			if(sum>0){
				for(int k=0; k<v.length; k++){
					v[k].transfer(sum);
				}
			}
			return sum;
		}
		@Override
		public void mult(double p) {
			for(int k=0; k<v.length; k++){
				v[k].prob = v[k].prob*p;
			}
			
		}
		@Override
		public double sum() {
			double sum = 0;
			for(int k=0; k<v.length; k++){
				sum+=v[k].prob;
			}
			return sum;
		}
	}
	static class TwoDMultiDM extends AbstractMultiDM{
	Entry[][] v;
	TwoDMultiDM(int[] len, double p){
		dim = 2;
		v = new Entry[len[0]][len[1]];
		double p1 = p/((double)len[0]*(double)len[1]);
		for(int k=0; k<len[0]; k++){
			for(int k1=0; k1<len[1]; k1++){
				v[k][k1] = new Entry(p1);
			}
		}
	}
		@Override
		Entry get(int[] inds) {
			return v[inds[0]][inds[1]];
		}

		@Override
		void set(int[] inds,double[] v) {
			// TODO Auto-generated method stub
			this.v[inds[0]][inds[1]].vals = v;
		}
		
		public double transfer(int j){
	  		double sum = 0;
			for(int k=0; k<v[j].length; k++){
				sum+=v[j][k].count;
			}
			if(sum>0){
				for(int k=0; k<v[j].length; k++){
					v[j][k].transfer(sum);
				}
			}
			return sum;
		}
		@Override
		public double sum() {
			double sum = 0;
			for(int k=0; k<v.length; k++){
				for(int k1=0; k1<v[k].length; k1++){
					sum+=v[k][k1].prob;
				}
			}
			return sum;
		}
		public void mult(int j, double p){
			for(int k=0; k<v[j].length; k++){
				v[j][k].prob = v[j][k].prob* p;
			}
		}
		@Override
		public void mult(double p){
			for(int j=0; j<v.length; j++){
				for(int k=0; k<v[j].length; k++){
					v[k][j].prob = v[k][j].prob* p;
				}
			}
		}
		@Override
		public double transfer() {
			double sum = 0;
			int len1 = this.v.length;
			double[] sums = new double[len1];
			for(int k=0; k<len1; k++){
				sums[k] = transfer(k);
				sum+=sums[k];
			}
			if(sum>0){
				for(int k=0; k<len1; k++){
					mult(k,sums[k]/sum);
				}
			}
			return sum;
		}
	}
	public static  AbstractMultiDM getMultiDM(int[] len2, int dim, double p) {
		if(dim==1) return new OneDMultiDM(len2,p);
	//	else if(dim==2) return new AbstractMultiDM.TwoDMultiDM(len2,p);
		else return new MultiDM(len2, dim,p);
	}
	static class MultiDM extends AbstractMultiDM{
		   AbstractMultiDM[] matr;
		   MultiDM(int[] len, int dim, double p){
			  this.dim = dim;
			  matr = new AbstractMultiDM[len[this.dim-1]];
			  for(int i=0; i<matr.length; i++){
				matr[i] = getMultiDM(len, dim-1,p/(double)matr.length);
			  }
			} 
		   
		

			@Override
			Entry get(int[] inds) {
				return matr[inds[dim-1]].get(inds);
			}

			@Override
			void set(int[] inds, double[] v) {
				matr[inds[dim-1]].set(inds, v);
			}



			@Override
			public void mult(double p) {
				for(int k=0; k<matr.length; k++){
					matr[k].mult(p);
				}
				
			}



			@Override
			public double transfer() {
				double[] sums = new double[matr.length];
				double sum=0;
				for(int k=0; k<matr.length; k++){
					sums[k] = matr[k].transfer();
					sum+=sums[k];
				}
				if(sum>0){
					for(int k=0; k<matr.length; k++){
						matr[k].mult(sums[k]/sum);
					}
				}
				return sum;
			}



			@Override
			public double sum() {
				double sum=0;
				for(int k=0; k<matr.length; k++){
					sum+= matr[k].sum();
				}
				return sum;
			}
		}
	
	
}

