package lc1.dp.illumina;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;
import java.util.concurrent.Callable;

import lc1.dp.data.collection.Info;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.stats.IlluminaRDistribution;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.PseudoDistribution;
import lc1.stats.PseudoMixture;
import lc1.stats.SimpleExtendedDistribution;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.util.Constants;

public class DistributionCollection extends AbstractDistributionCollection {
public final boolean probeOnly;
	public static AbstractDistributionCollection dc=null;
	
	//public IlluminaProbB probB;
//public	IlluminaProbR probR;
	public IlluminaProbR probRB;
	
	@Override
	public void addMixture(HaplotypeEmissionState st, Boolean[] probeOnly) {
		int ploidy = st.noCop();
		PseudoDistribution[] emissions = st.emissions;
		// = Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(null);
		for(int k=0; k<emissions.length; k++){
			
				PseudoDistribution nullDist =  Emiss.getSpaceForNoCopies(ploidy).getHWEDist1(probeOnly[k]==null || probeOnly[k] ?0.0: null);
		      
			SimpleExtendedDistribution1 dist = Constants.allowLocalDist() ?this.probRB.mixeR[k] : probRB.mixeglobalR;
			if(emissions[k]==null) {
				emissions[k] = nullDist;
			}
			else if(dist.probs[0]<1){
			if(emissions[k] instanceof IlluminaRDistribution) emissions[k] =  new PseudoMixture(
					 new PseudoDistribution [] {emissions[k], 	nullDist}, dist);
			}
			
			
		}
	}
	
public void setMinMax(double minR, double maxR, double minB, double maxB){
	//if(probR!=null) probR.setMinMax(minR, maxR);
	if(probRB!=null) probRB.setMinMax(minR, maxR, minB, maxB);
	//if(probR1!=null) probR1.setMinMax(minR, maxR);
}
	
//public DistributionCollection clone(){
	//return new DistributionCollection(probR.clone(), probB.clone());
	//return dc;
//}
public void initialise(){
	//if(probB!=null) this.probB.initialiseBCounts();
	//if(probR!=null) this.probR.initialiseRCounts();
	if(probRB!=null) this.probRB.initialiseRCounts();
	//if(probR1!=null) probR1.initialiseRCounts();
}
	public void print(File out){
		try{
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(out, "clusters.txt"))));
		//if(probB!=null)probB.print(pw);
		//if(probR!=null) probR.print(pw);
		if(probRB!=null) probRB.print(pw);
		pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	public void print(int i){
		probRB.print(i);
	}
public DistributionCollection(DistributionCollection dc, boolean sameMix){
	//if(dc.probB!=null) this.probB = new IlluminaProbB(dc.probB);
	//if(dc.probR!=null) this.probR = new IlluminaProbR(dc.probR);
	this.probeOnly = dc.probeOnly;
	if(dc.probRB!=null) this.probRB = new IlluminaProbR(dc.probRB, sameMix);
		 this.data_index = dc.data_index;
		//this.lrr = dc.lrr;
		//this.zeroR = dc.zeroR;
}

public  void setIndex(short i){
	this.data_index = i;
}






public DistributionCollection(
			List<Integer> backgroundCN, int index, int noSnps1,  File dir) {
	int noSnps = Constants.allowLocalDist() ? noSnps1 : 0;
	this.probeOnly = Constants.probeOnly(index);
	EmissionStateSpace stSp = Emiss.getSpaceForNoCopies(Constants.maxPloidy());
	
	    if(Constants.joint){
	    	File clust1 = Constants.readGlobalClusterFile() ? null : Constants.getClusterFile(index);
	    	
	    
		this.probRB = 
			
				            	 
				            	/* Constants.format()[index].toLowerCase().equals("illuminaxy") ?
				            				(IlluminaProbR)		 new IlluminaXY(stSp,
				            					                true,index,noSnps,  probeOnly, dir
				            					             ):*/
				            					            	 
				            	 Constants.format()[index].toLowerCase().equals("depth") ?
				            			 
				            			 new SequenceDepth(stSp, true, index, noSnps,  null):
				            				 Constants.format()[index].toLowerCase().equals("allelecount") ?
				            						 new IlluminaProbRBCounts(stSp, true,index,noSnps, 
				            		                            clust1==null ?  new File(dir, "clusters.txt"): clust1
				            		                         ) :
			probeOnly ? 
			new IlluminaProbR(stSp, true,index,noSnps, 
                clust1==null ?  new File(dir, "clusters.txt"): clust1
             ) : 
            		new IlluminaProbRB(stSp, true,index,noSnps, 
                            clust1==null ?  new File(dir, "clusters.txt"): clust1
                         ) ;
            	 
	    }
	    else{
		/*this.probB=
	    		 Constants.format()[index].toLowerCase().equals("ascn") ?
	    		  new IlluminaProbTheta(stSp, true, 1.0, index, noSnps):
	   		  new IlluminaProbB(stSp, true, 1.0, index, noSnps);*/
	      }
	if(Constants.plot>0){
		this.r_formulae = new String[noSnps][Constants.probeOnly[index] ? 2: 4];
	
	}
	//lrr = probRB.getAssumedLRRCoeff(
		//	Constants.r_mean(index), Constants.r_x(index), Constants.backgroundCount(this.data_index)).viewColumn(0); 
	//zeroR = Constants.r_mean(index)[0];
}



	public DistributionCollection( IlluminaProbR clone3) {
		// TODO Auto-generated constructor stub
		//this.probB = clone2;
		//this.probR = clone;
		this.probRB = clone3;
		this.probeOnly = !(clone3 instanceof IlluminaProbRB); 
	//	lrr = probRB.getAssumedLRRCoeff(
		//		Constants.r_mean(0), Constants.r_x(0), Constants.backgroundCount(this.data_index)).viewColumn(0);
	
		//this.probR1= clone3;
	}
	//final DoubleMatrix1D lrr;
	//double zeroR;

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#getDistribution(short, int, int, int)
	 */
	public ProbabilityDistribution getDistribution(short data_index, int cn_bottom,
			int cn_top, int pos) {
		return null;//this.probR.getDistribution(cn_bottom, cn_top,pos);
	}
	public Double getFrac( int data_index,int i, boolean r){
		if(r){
		if(probRB!=null && probRB.mixeR[i]!=null) 
		 return probRB.mixeR[i].probs[0];
		else return 1.0;//probR.mixe[i].probs[0];
		}
		else{
		//	if(probRB!=null && probRB.mixeB[i]!=null) 
			//	 return probRB.mixeB[i].probs[0];
				//else
					return 1.0;//probR.mixe[i].probs[0];
		}
	}
	
	public Double getFracGlob(int ind, boolean po, boolean r){
		if(r){
		if(probRB!=null) 
		 return 
			 probRB.mixeglobalR.probs[0];
		else return 1.0;
		}
		else{
			if(probRB!=null) 
				 return 
							 probRB.mixeglobalR.probs[0]
							                          ;
				else return 1.0;
		}
	}
	
	public static Double getFracS( int data_index,int i, boolean r){
		 if(dc==null) return 1.0;
			else{
				if(Constants.allowLocalDist)
				return dc.getFrac(data_index, i, r);
				else return dc.getFracGlob(data_index, true, r);
			}
	}
	public static Double getFracGlobS(int ind, boolean po, boolean r){
		 if(dc==null) return 1.0;
		else return dc.getFracGlob(ind, po, r);
	}
	
	

	public ProbabilityDistribution getDistribution1(short data_index, int cn_bottom,
			int cn_top, int pos) {
		return null;//this.probR1.getDistribution(cn_bottom, cn_top,pos);
	}

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#scoreB(short, int, int)
	 */
	public double  scoreB(short data_index, int j,double b, int i) {
		return 0.0;//this.probB.calcB(j, b, i);
		
	}
	
	public double  scoreRB(short data_index, int j,double r, double b, int i) {
		return this.probRB.calcR(j, r, b, i);
		
	}

	//short zero = (short)0;
	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#scoreR(short, int, int, int)
	 */
	public double scoreR(short data_index, int backgroundCount, int no_cop,
			double r, int i) {
		return 0.0;//this.probR.calcR(backgroundCount, no_cop, r, i);
	}

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#addBCount(short, int, int, double, int)
	 */
	public void addBCount(short data_index,  int j,
			double weight, double val, int i) {
		//this.probB.addBCount(j, weight, val, i);
		
	}
	
	public void addRBCount(short data_index,  int j,
			double weight, double valR, double valB, int i) {
	//	if(Constants.allowLocalDist())
		this.probRB.addRCount(j, weight, valR,valB,  i);
		
		
	}
	
	public void drop(List<Integer> toDrop) {
	this.probRB.drop(toDrop);
	}

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#addRCount(short, int, int, double, int)
	 */
	public void addRCount(short data_index, int backgroundCount, int no_cop,
			double weight, double val, int i) {
		//this.probR.addRCount(backgroundCount, no_cop, weight, val, i);
		
	}
/* (non-Javadoc)
 * @see lc1.dp.illumina.AbstractDistributionCollection#getDistribution(short, int, int)
 */
	public ProbabilityDistribution getDistribution(short data_index, int j, int i) {
		
		return null;//0.0;//this.probB.getDistribution(j,i);
	}
	
	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#getDistribution(short, int, int)
	 */
		public ProbabilityDistribution2 getDistributionRB(short data_index, int n,int noB, int i) {
			if(Constants.allowLocalDist())
			return this.probRB.getDistribution(n,noB, i);
			else {
				//if(probRB.probeOnly[i]==null) return 
				return this.getDistributionRBGlob(data_index, n, noB);
			}
		}
		
public ProbabilityDistribution2 getDistributionRBGlob(short data_index, int noCop,int noB) {
			 return this.probRB.getDistributionGlobal(noCop, noB);
		
		}
		
		
		@Override
		public void addRBCount(short data_index, int noCop, int noB, double val,
				Double r, Double b, int i) {
		//	if(Constants.allowLocalDist())
			this.probRB.addCount(noCop, noB, i, r, b, val);
			//.getDistributionRB(this.data_index, noCop, noB,i).addCount(r,b, val);
			
		}
	
	public ProbabilityDistribution getDistributionBfrac(short data_index, int j, int i) {
		return null;//this.probB.b[i][j];//getDistribution(j,i);
	}
	short  data_index;
	 

	
	
	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#maximisationStep(int, double, java.util.List)
	 */
    @Override
	public  void maximisationStep( final int i,final double[] pseudo, final double[] pseudoG, List tasks){
  
    	 this.regressGlobal(pseudoG, i);
      if(Constants.getMinValue(Constants.r_train(0))<1000 || Constants.getMinValue(Constants.b_train(0))<1000 || 
    		  Constants.getMinValue(Constants.r_train(0))<1000 || Constants.getMinValue(Constants.b_train(0))<1000 	  ){
    	  if(probRB!=null){
    		  double mp = pseudo[2]*Constants.pseudoMod(this.data_index);
				if(mp<1e4){
					probRB.transferMix(mp);
					
				}
    		  for(int k=0; k<this.probRB.r.length; k++){
  				final int kk=k;
  				tasks.add(new Callable(){

  					@Override
  					public Object call() throws Exception {
  						
  				//		System.err.println("counts mixe "+probRB.mixe[kk].counts[0]+" "+probRB.mixe[kk].counts[1]+" "+probRB.mixe[kk].probs[0]);
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
      
      
     //}
         
    }
    public void regress(int pos, double[] pseudo1, int it){
    	//pseudo, double[][] pseudoR, double[][] pseudoB, double[] pseudoRho
//     		){
    	double[] pseudo;
    	double mod = Constants.pseudoMod(this.data_index);
    	if(mod!=1){
    		pseudo = new double[pseudo1.length];
    		for(int i=0; i<pseudo.length; i++){
    			pseudo[i] = pseudo1[i];
    		}
    	}
    	else{
    		pseudo = pseudo1;
    	}
    	IlluminaProbR[]probRB = new IlluminaProbR[] {this.probRB};
    	int[] pos1 = new int[] {pos};
    //	if(this.probRB.probeOnly[pos]==null) return;
         MultipleDataRegression reg = 
        	 this.probeOnly ?  new MultipleDataRegression(probRB, pos1,pseudo,0, probRB[0].prior):
        	 new MultipleDataRegressionRB(probRB, pos1,pseudo,0, probRB[0].prior

 		    		 
         );
        	
       // boolean probeOnly = this.probRB.probeOnly[pos];	
        		
 	
       boolean setToExclude =  reg.makeRegress( true, false);
      if(this.r_formulae!=null){
    	  reg.rst(r_formulae[pos]);
         reg.bst(r_formulae[pos]);
      }
      
       if(setToExclude){
   		for(int ik=0; ik<probRB.length; ik++){
   			probRB[ik].setToExclude(pos);
   		}
   	}
        } 
    public String[] getFormG(int i2){
	  
    	return this.probRB.globalRst;
    }
    @Override
    public String[] getForm(int i2, int i) {
		// TODO Auto-generated method stub
    	
    		if(r_formulae!=null){
    			return this.r_formulae[i];
    		}
    		else return super.getForm(i2, i);
	}
    

	public  static String[] getFormGlobS(int ind){
		 if(dc==null) return new String[] {"",""};
		else return dc.getFormGlob(ind);
	}
	public String[] getFormGlob(int ind){
	 if(probRB==null) return null;
	 else {
		 return probRB.globalRst;
	 }
	}
    
    public static String[] getFormS(int i2, int i) {
		if(dc==null) return null;
		else{
			if(Constants.allowLocalDist)
			return dc.getForm(i2, i);
			else return null;
		}
	}
    
    public void regressGlobal(double[] pseudo, int it){
    	double[] pseudo1;
    	double mod = Constants.pseudoMod(this.data_index);
    	if(mod!=1){
    		pseudo1 = new double[pseudo.length];
    		for(int i=0; i<pseudo1.length;i ++){
    			pseudo1[i] = pseudo[i] *mod;
    		}
    	}
    	else{
    		pseudo1 = pseudo;
    	}
    	try{
    	
    	this.probRB.maximiseGlobal(pseudo1, it);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    
    }
    
    
   

	

	

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#getBName(int)
	 */
	public String getBName(int ij) {
		//String name =  this.probB.b[0][ij].name();
		return  "null";
	}

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#sampleR(int, int, int, int)
	 */
	public Double sampleR(int data_index, int cn_bg, int cn_fg, int pos) {
		return null;//this.probR.sampleR(cn_bg, cn_fg, pos);
	}

	/* (non-Javadoc)
	 * @see lc1.dp.illumina.AbstractDistributionCollection#sampleB(int, int, int)
	 */
	public Double sampleB(int data_index, int obj_index, int pos) {
		return null;// this.probB.sampleB(obj_index, pos);
	}

	/*@Override
	public Color[] getColB() {
		if(probRB!=null) return probRB.color;
		return null;//this.probB.color;
	}*/
	
	/*@Override
	public Shape[] getShapeB() {
		if(probRB!=null) return probRB.shape;
		return null;//this.probB.shape;
	}*/

	/*@Override
	public Color[] getColR() {
		// TODO Auto-generated method stub
		 return probRB.color;
	}*/

	@Override
	public Double minQuality(int relative_position) {
	  SimpleExtendedDistribution dist = 
		  Constants.allowLocalDist() ? 
		  this.probRB.mixeR[relative_position] : this.probRB.mixeglobalR;
	  if(dist==null) return null;
	  else return dist.probs(1);
	}

	public void append(DistributionCollection dc2, Info[][] info)  throws Exception{
		if(info!=null){
		IlluminaProbR[] dists = new IlluminaProbR[] {this.probRB, dc2.probRB};
		IlluminaProbR.getMerged(dists, info);
		}
		else	this.probRB.append(dc2.probRB);
		
	}

	

	

	

	
	
}
