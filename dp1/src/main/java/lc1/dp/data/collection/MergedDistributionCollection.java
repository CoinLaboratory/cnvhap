package lc1.dp.data.collection;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;

import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.illumina.IlluminaProbR;
import lc1.dp.illumina.MultipleDataRegression;
import lc1.dp.illumina.MultipleDataRegressionRB;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;

public class MergedDistributionCollection extends AbstractDistributionCollection {

	Info[][] map;
	
	public AbstractDistributionCollection[] dc;
	//List<DoubleMatrix1D> coeff;
	//List<Double>zeroR;
	int numTypes=0;
//	final RegressParams[] paramR, paramR_probeOnly, paramRVar, paramRVar_probeOnly;
//	final RegressParams[] paramBaf,  paramBafVar, paramCovar;
	//nme is used to decide if same regression model can be used
	final String[] nme;
	
	
	final double[] pseudoMod;
	public MergedDistributionCollection(AbstractDistributionCollection[] dc,
			Info[][] map, String[] nme) {
		this.dc = dc;
		this.map = map;
		this.nme = nme;
	this.probeOnly= new Boolean[dc.length];
		//	DoubleMatrix1D[] lrrCoeff = new DoubleMatrix1D[dc.length];
	//	 coeff = new ArrayList<DoubleMatrix1D>();
	//	zeroR = new ArrayList<Double>();
			   type_alias = new int[dc.length];	
			  numTypes =0;
			  Class[] clazz = new Class[dc.length]; 
			 Arrays.fill(type_alias, -1);
		inner: for(int i=0; i<dc.length; i++){
			if(dc[i]==null) continue;
			if(nonnull<0) nonnull =i;
			//lrrCoeff[i] = ((DistributionCollection)dc[i]).probRB.getAssumedLRRCoeff(
				//	Constants.r_mean(i), Constants.r_x(i), Constants.backgroundCount(i)).viewColumn(0);
			clazz[i] = ((DistributionCollection)dc[i]).probRB.getClass();
			for(int j=0; j<i; j++){
				if(dc[j]!=null){
				if( clazz[i].equals(clazz[j]) && nme[i].equals(nme[j]) &&
						((DistributionCollection)dc[i]).probRB.canMergeWith(((DistributionCollection)dc[j]).probRB)
					//	Constants.basisNme(i).equals(Constants.basisNme(j))
						){
					type_alias[i] = type_alias[j];
					continue inner;
				}
				}
			}
			type_alias[i] = numTypes;
		//	coeff.add(lrrCoeff[i]);
		//	zeroR.add(Constants.r_mean(i)[0]);
			numTypes++;
		}
			 nonNull=new int[numTypes];
		     pos1 = new int[numTypes][];
		  probRB = new IlluminaProbR[numTypes][];
		  this.pseudoMod = new double[numTypes];
		  for(int i=0; i<type_alias.length; i++){
			  if(type_alias[i]>=0)
			  pseudoMod[type_alias[i]] = Constants.pseudoMod(i);
		  }
	//	  probeOnly = new boolean[numTypes];
		 this.dataIndex =new List[numTypes];
		 for(int i=0; i<dataIndex.length; i++){
			 dataIndex[i] = new ArrayList<Integer>();
			 for(int k=0; k<type_alias.length; k++){
				 if(type_alias[k]==i){
					 dataIndex[i].add(k);
				 }
			 }
		 }
		// TODO Auto-generated constructor stub
	}
	
	
	private boolean equals(DoubleMatrix1D doubleMatrix1D,
			DoubleMatrix1D doubleMatrix1D2) {
		for(int i=0; i<doubleMatrix1D.size(); i++){
			if(Math.abs(doubleMatrix1D.get(i) - doubleMatrix1D2.get(i))>0.001) return false;
		}
		return true;
	}
	private boolean equals(double[] doubleMatrix1D,
			double[] doubleMatrix1D2) {
		for(int i=0; i<doubleMatrix1D.length; i++){
			if(Math.abs(doubleMatrix1D[i] - doubleMatrix1D2[i])>0.001) return false;
		}
		return true;
	}


	final int[] type_alias;

	public void initialise(){
		for(int i=0; i<dc.length; i++){
			if(dc[i]!=null) dc[i].initialise();
		}
	}
	
	public void print(File out){
		for(int i=0; i<dc.length; i++){
			
			if(dc[i]!=null){
				File f1 = new File(out, nme[i]);
				f1.mkdir();
				//pw.println("Illumina dists for index\t"+i+"\t"+nme[i]);
				dc[i].print(f1);
			}
		}
	}
	short zero = (short)0;
	@Override
	public void addBCount(short data_index, int j,
			double weight, double val, int i) {
		this.dc[data_index].addBCount(zero,  j, weight, val, 
		this.map[i][data_index].relative_position		
		);
		
	}

	@Override
	public void addRCount(short data_index, int backgroundCount, int no_cop,
			double weight, double val, int i) {
		this.dc[data_index].addRCount(zero, backgroundCount, no_cop, weight, val, 
				this.map[i][data_index].relative_position		
				);
		
	}

	
	
	@Override
	public String getBName(int ij) {
		// TODO Auto-generated method stub
		return this.dc[0].getBName(ij);
	}

	

	/*@Override
	public Color[] getColR() {
		// TODO Auto-generated method stub
		return dc[nonnull].getColR();
	}*/
	 int nonnull=-1;

	/*@Override
	public Color[] getColB() {
		// TODO Auto-generated method stub
		return dc[0].getColB();
	}*/
	
	/*@Override
	public Shape[] getShapeB() {
		// TODO Auto-generated method stub
		return dc[0].getShapeB();
	}*/

	@Override
	public ProbabilityDistribution getDistribution(short data_index,
			int cn_bottom, int cn_top, int pos) {
		Info inf = this.map[pos][data_index];
		if(inf==null) return null;
		return this.dc[data_index].getDistribution(zero, cn_bottom, cn_top,  
				inf.relative_position		
				);
	}
	
	public void print(int i){
		Info[] inf = this.map[i];
		for(int k=0; k<inf.length; k++){
			if(inf[k]!=null){
			System.err.println(k);
				this.dc[k].print(inf[k].relative_position);
			}
		}
		//probRB.print(i);
	}
	
	public ProbabilityDistribution getDistribution1(short data_index,
			int cn_bottom, int cn_top, int pos) {
		return this.dc[data_index].getDistribution1(zero, cn_bottom, cn_top,  
				this.map[pos][data_index].relative_position		
				);
	}

	@Override
	public ProbabilityDistribution getDistribution(short data_index, int j,
			int i) {
		Info inf = this.map[i][data_index];
		if(inf==null) return null;
		return this.dc[data_index].getDistribution(zero, j,  
				inf.relative_position		
				);
	}
	
	@Override
	public ProbabilityDistribution getDistributionBfrac(short data_index, int j,
			int i) {
		return this.dc[data_index].getDistributionBfrac(zero, j,  
				this.map[i][data_index].relative_position		
				);
	}

	@Override
	public void maximisationStep(final int i,final double[] pseudo, final double []pseudo_global,List tasks) {
		if(Constants.sepModels()){
			for (int di=0; di<this.dc.length; di++){
				dc[di].maximisationStep(i, pseudo, pseudo_global, tasks);
				
			}
		}
		else{
			if(pseudo[2]<1e5){
			for (int di=0; di<this.dc.length; di++){
				if(dc[di]!=null){
				((DistributionCollection)dc[di]).probRB.transferMix(pseudo[2] *Constants.pseudoMod(di));//.maximisationStep(i, pseudo,  pseudo_mix, tasks);
				}
			}
			}
			if(Constants.measureGlobal()){
				for (int di=0; di<this.dc.length; di++){
					if(dc[di]!=null){
					((DistributionCollection)dc[di]).regressGlobal(pseudo_global, i);
					}
				}
			}
			for(int k=0; k<this.map.length; k++){
				final int kk=k;
				tasks.add(new Callable(){

					@Override
					public Object call() throws Exception {
					//	maximiseMixtures(kk, pseudo_mix);
				//		System.err.println("max "+i);
						regress(kk, pseudo,i);
						/*for (int di=0; di<this.dc.length; di++){
							((DistributionCollection)dc[di]).probRB.rMaximistation(kk,
									pseudo, Constants.r_train(2), Constants.b_train(2), 
									Constants.rho_train(3), kk);
							
						}*/
						return null;
					}
					
				});
			}
		}
	}
	
	
	
	int[] nonNull;
	int[][]pos1;
	IlluminaProbR[][] probRB;
	
	final Boolean[] probeOnly;
	List<Integer>[] dataIndex;
	
	int[][] trainingGroups = null;
	
	public Double minQuality(int i){
		Info[] inf = this.map[i];
		double v = 0.0;
		for(int k=0; k<inf.length; k++){
			if(inf[k]!=null && dc[k]!=null){
				Double v1 = this.dc[k].minQuality(inf[k].relative_position);
				if(v1!=null && v1>v){
					v = v1;
				}
			}
		}
		return v;
	}
	/*pseudoR1 is for pseudo counts tying individual clusters to regression based averages */
    public void regress(int pos, double[] pseudo1, int it){
    	//pseudo, double[][] pseudoR, double[][] pseudoB, double[] pseudoRho
//     		){
    	//System.err.println("maximising "+pos)
    	
    	Info[] inf = this.map[pos];
    	Arrays.fill(nonNull,0);
    	Arrays.fill(probeOnly, true);
    	for(int i=0; i<inf.length; i++){
    		if(inf[i]!=null && dc[i]!=null) {
    		//	if(((DistributionCollection)this.dc[i]).probRB.probeOnly[inf[i].relative_position]!=null){
    			 nonNull[this.type_alias[i]]++;
    		//	}
    			 probeOnly[this.type_alias[i]] =  probeOnly[this.type_alias[i]] & ((DistributionCollection)dc[i]).probeOnly;
    		}
    	}
    
    	for(int i=0; i<pos1.length; i++){
    		pos1[i] = new int[nonNull[i]];
    		probRB[i] = new IlluminaProbR[nonNull[i]];
    		
    	}
       
         Arrays.fill(nonNull, 0);
       //  Arrays.fill(probeOnly, true);
         for(int i=0; i<inf.length; i++){
     		if(inf[i]!=null && dc[i]!=null){
     			//if(((DistributionCollection)this.dc[i]).probRB.probeOnly[inf[i].relative_position]!=null){
     				probRB[type_alias[i]][nonNull[type_alias[i]]] = ((DistributionCollection)this.dc[inf[i].data_index]).probRB;
     				pos1[type_alias[i]][nonNull[type_alias[i]]] = inf[i].relative_position;
     	//			probeOnly[type_alias[i]] =
     	//				probeOnly[type_alias[i]] && 
     		//			((DistributionCollection)this.dc[inf[i].data_index]).probRB.probeOnly[inf[i].relative_position];
     				nonNull[type_alias[i]]++;
     			//}
     		}
     	}
         for(int i=0; i<probRB.length; i++){
         if(probRB[i].length>0 && (i==0 || probRB[i]!=probRB[0])){
        	 double mod = this.pseudoMod[i];
        	 double[] pseudo; 
        	 if(mod!=1.0){
        		 pseudo = new double[pseudo1.length];
        	
        	 for(int k=0; k<pseudo.length; k++){
        		 pseudo[k] = pseudo1[k] *mod;
        	 }
        	 }
        	 else{
        		 pseudo = pseudo1;
        	 }
        	 try{
        		
         MultipleDataRegression reg =probeOnly[i] ?  new MultipleDataRegression(probRB[i], pos1[i],pseudo,0, probRB[i][0].prior ):
        	 new MultipleDataRegressionRB(probRB[i], pos1[i],pseudo,0 ,probRB[i][0].prior);
        	
 		    
        	boolean setToExclude = reg.makeRegress( true,false);
        	if(setToExclude){
        		for(int ik=0; ik<probRB[i].length; ik++){
        			probRB[i][ik].setToExclude(pos);
        		}
        	}
        	 }catch(Exception exc){
        		 System.err.println("prob at i "+DataCollection.datC.loc.get(i)+DataCollection.datC.snpid.get(i));
        		 exc.printStackTrace();
        		 System.exit(0);
        	 }
 	   
         }
        }
    }
	

	@Override
	public Double sampleB(int data_index, int obj_index, int pos) {
		return this.dc[data_index].sampleB(zero, obj_index, map[pos][data_index].relative_position);
	}

	@Override
	public Double sampleR(int data_index, int cn_bg, int cn_fg, int pos) {
		return this.dc[data_index].sampleR(zero, cn_bg, cn_fg, map[pos][data_index].relative_position);
		
	}
	
	@Override
	public double scoreB(short data_index, int j, double b, int i) {
		return this.dc[data_index].scoreB(zero, j, b, map[i][data_index].relative_position);

	}

	@Override
	public double scoreR(short data_index, int cn_bg, int cn_fg,
		double r, 	int i) {
		return this.dc[data_index].scoreR(zero, cn_bg, cn_fg, r,map[i][data_index].relative_position);
		
	}

	@Override
	public void addRBCount(short data_index, int j, double weight, double valR,
			double valB, int i) {
		// TODO Auto-generated method stub
		this.dc[data_index].addRBCount(data_index, j, weight, valR, valB, this.map[i][data_index].relative_position	);
	}

	@Override
	public ProbabilityDistribution2 getDistributionRB(short data_index, int noCop, int noB,
			int i) {
		Info inf =map[i][data_index]; 
		if(inf==null){
			return null;
		}
		if(dc[data_index]==null){
			return null;
		}
		// TODO Auto-generated method stub
		return this.dc[data_index].getDistributionRB(data_index, noCop,noB, 
				inf.relative_position);
	}
	@Override
	public  ProbabilityDistribution2 getDistributionRBGlob(short data_index, int noCop,int noB) {
		
		if(dc[data_index]==null) return null;
		// TODO Auto-generated method stub
		return this.dc[data_index].getDistributionRBGlob(data_index, noCop,noB);
	}
	@Override
	public void addRBCount(short data_index, int noCop, int noB, double val,
			Double r, Double b, int i) {
		Info inf =map[i][data_index]; 
		if(inf==null){
			return;
		}
		if(dc[data_index]==null) return ;
		// TODO Auto-generated method stub
		this.dc[data_index].addRBCount(data_index, noCop, noB, val, r, b, 
				inf.relative_position);
		//.getDistributionRB(data_index, noCop,noB, 
			//	inf.relative_position);
		
		//.getDistributionRB(this.data_index, noCop, noB,i).addCount(r,b, val);
		
	}

	@Override
	public double scoreRB(short data_index, int j, double r, double b, int i) {
		// TODO Auto-generated method stub
		return this.dc[data_index].scoreRB(data_index, j, r, b, map[i][data_index].relative_position);
	}
	
	public Double getFrac( int data_index,int i, boolean r){
		Info inf = map[i][data_index];
		if(inf==null || dc[data_index]==null) return null;
		else return this.dc[data_index].getFrac(0, inf.relative_position, r);
	}
	
	public String[] getForm(int data_index, int i) {
		Info inf = map[i][data_index];
		if(inf==null || dc[data_index]==null) return null;
		return this.dc[data_index].getForm(0, inf.relative_position);
	}
	public Double getFracGlob(int data_index, boolean pi, boolean r){
		
		if( dc[data_index]==null) return null;
		else return this.dc[data_index].getFracGlob(0, pi, r);
	}


	@Override
	public String[] getFormGlob(int data_index) {
		if( dc[data_index]==null) return null;
		else return this.dc[data_index].getFormGlob(0);
	}


	@Override
	public void addCollection(AbstractDistributionCollection dc) {
		MergedDistributionCollection mdc = (MergedDistributionCollection)dc;
		for(int i=0; i<this.dc.length; i++){
			this.dc[i].addCollection(mdc.dc[i]);
		}
		
	}


	@Override
	public void reverse() {
		for(int i=0; i<this.dc.length; i++){
			dc[i].reverse();
		}
		
	}

}
