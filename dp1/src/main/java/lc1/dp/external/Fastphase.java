  package lc1.dp.external;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.zip.ZipFile;

import javax.swing.JOptionPane;

import com.sun.jmx.snmp.tasks.ThreadService;

import lc1.dp.appl.CNVHap;
import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.Sampler;
import lc1.dp.data.collection.AssociationCalculator;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.HWECalculator;
import lc1.dp.data.collection.IlluminaRDataCollection;
import lc1.dp.data.collection.LinearRegressionCalculator;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.illumina.AbstractDistributionCollection;
import lc1.dp.illumina.DistributionCollection;
import lc1.dp.model.CachedHMM;
import lc1.dp.model.CollapsedHMM;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.model.CompressHMM;
import lc1.dp.model.ExpandedHMM;
import lc1.dp.model.FreeHaplotypeHMM;
import lc1.dp.model.HaplotypeHMM;
import lc1.dp.model.MarkovModel;
import lc1.dp.model.PairMarkovModel;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.swing.HMMPanel;
import lc1.dp.swing.HWEPlot;
import lc1.dp.swing.HaplotypePanel1;
import lc1.dp.swing.IndividualPlot;
import lc1.dp.swing.RateMatrixPanel;
import lc1.dp.swing.SignificancePlot;
import lc1.dp.transition.ExponentialTransitionProbs;
import lc1.stats.Dirichlet;
import lc1.stats.SimpleExtendedDistribution;
import lc1.util.Constants;
import lc1.util.PseudoIterator;
/** This is the primary class for running fastphase */
public class Fastphase {
    final public static Logger logger = Logger.getAnonymousLogger();
    static long time = System.currentTimeMillis();
    List<CompoundMarkovModel> best = null; 
    Double bestSc = null;
    public static MarkovModel model;  
    
    private static void change(Comparable[] comp){
           Comparable tmp = comp[0];
           comp[0] = comp[1];
           comp[1] = tmp;
       }
   
final int noSites;
    
   private static Double[] geom(Double[] d){
        double[] prod = new double[] {1.0, 0.0, 1.0};
        
        for(int i=0; i<d.length; i++){
            prod[0]+=Math.log(d[i]);
            prod[1]+=d[i];
            prod[2] = Math.min(prod[2], d[i]);
            
        }
        double len = (double)d.length;
        return new Double[] {Math.exp(prod[0]/len),prod[1]/len ,prod[2]};
    }
    
    
  final DataCollection dataC;
    
    static boolean isSame(Comparable[] emiss){
        return emiss[0]==emiss[1];
    }
   
   
   /*private double getMax(HaplotypeHMM hapHMM, int no){
       if(hapHMM instanceof ExponentialHaplotypeHMM) return 1;
       Double[][] prob = hapHMM.getHittingProb(hapHMM.noSnps);
       double max = Double.NEGATIVE_INFINITY;
       for(int i=0; i<no; i++){
           double max_i = getMax(prob[prob.length-1-i]);
           if(max_i > max) max = max_i;
       }
       return  max;
   }
 /*  private boolean converged(HaplotypeHMM hapHMM, double thres, int no){
       double max = getMax(hapHMM, no);
       logger.info("MAX is "+max);
       if(max <thres) return true;
       else{
          
           return false;
       }
   }*/
   private double getMax(Double[] doubles) {
       double max = Double.NEGATIVE_INFINITY;
       for(int i=0; i<doubles.length; i++){
           if(doubles[i]>max) max = doubles[i];
       }
       return max;
}
/** List of all models trained */
   final int modelLength;
   public void trainCases(int numIt1, PairMarkovModel hmm_union,PairMarkovModel hmm_cases, int rep, boolean emiss, boolean trans) throws Exception{
       logger.info("training cases");
       //TrainingElement train = new TrainingElement((HaplotypeHMM) ((HaplotypeHMM)hmm_union.m1).clone(hmm_cases.m1, 0.1), null);
      // for(int k1=0; k1<this.modelLength; k1++){
     //      train.train(numIt1, 1,true);
     //  }
   }
   
  public TrainingElement getTrainingElement(FreeHaplotypeHMM hmm){
      return new TrainingElement(hmm);
  }
   
   public Iterator<TrainingElement> getTrainingElements(final Iterator<MarkovModel> hmm) throws Exception{
       final TrainingElement first = new TrainingElement(hmm.next());
       return new Iterator<TrainingElement>(){
           boolean isfirst = true;
           
            public boolean hasNext() {
                return isfirst ||  hmm.hasNext();
            }
            public TrainingElement next() {
                if(isfirst){
                    isfirst =false;
                    return first;
                }
                else return new TrainingElement(hmm.next());
            }
            public void remove() {}
           
       };
   }
   public void prime(HaplotypeHMM hmm, int prime){
       List<EmissionState> l = new ArrayList<EmissionState>();
       for(Iterator<EmissionState> it = this.dataC.dataLvalues(); it.hasNext();){
           l.add(it.next());
       }
       EmissionStateSpace emStSp2 = Emiss.getEmissionStateSpace(1);
       EmissionStateSpace emStSp1 = ((CompoundEmissionStateSpace)emStSp2).getMembers()[0];
       int[][] map = new int[emStSp2.size()][];
       for(int k=0; k<map.length; k++){
           ComparableArray comp = (ComparableArray)emStSp2.get(k);
           map[k] = new int[comp.size()];
           for(int i=0; i<map[k].length; i++){
               map[k][i] = emStSp1.get(comp.get(i));
           }
       }
       int start = Constants.nextInt(hmm.noSnps-prime);
       for(int j=1; j<2;
      // hmm.modelLength(); 
       j++){
           EmissionState ld = l.get(j);
           HaplotypeEmissionState hes = (HaplotypeEmissionState)hmm.getState(j);
          
          // EmissionState lde = ld.state_indices;
           for(int i=start; i<ld.length() && i<start +prime; i++){
              double[] probs2 =  new double[emStSp1.size()];//hes.emissions[i].probs;
              Arrays.fill(probs2, 0);
              double[] probs1 = ld.getEmiss(i);
              double sum=0;
              for(int k=0; k<probs1.length; k++){
                
                 if(probs1[k]!=0){
                     int[] mapping = map[k];
                     for(int k1=0; k1<mapping.length; k1++){
                         probs2[mapping[k1]]+=probs1[k];
                         sum+=probs1[k];
                     }
                 }
              }
              for(int k=0; k<probs2.length; k++){
                  probs2[k] = probs2[k]/sum;
              }
              hes.emissions[i] = new SimpleExtendedDistribution(new Dirichlet(probs2, 10000));
            //  hes.emissions[i].transferProbToPseudo();
           }
          
       }
   }
   int overallIndex =0;
   public void train(TrainingElement train, boolean swtch) throws Exception{
       if(swtch){
           dataC.extractFromTrioData();
           this.updateData();
       }
         
          if(Constants.prime()>0){
              if(train.hapHMM instanceof FreeHaplotypeHMM) prime((FreeHaplotypeHMM)train.hapHMM, Constants.prime());
              else if(train.hapHMM instanceof CompoundMarkovModel){
                  MarkovModel mm = ((CompoundMarkovModel)train.hapHMM).getMemberModels()[1];
                  if(!mm.getName().startsWith("allele")) throw new RuntimeException("!! "+mm.getName());
                  prime((FreeHaplotypeHMM)mm, Constants.prime());
              }
          }
          
         
          logger.info("before simple emissions");
  boolean expanded = false;
    int[] numIt = Constants.numIt();
    boolean do_update = true;
    if(Constants.sum(numIt)==0 && Constants.plot>0){
    	numIt[0] = 1;
    	do_update = false;
    }
    if( Constants.writeHMM()>=1 && false)  writeHMM(train.hmm,
            new File(hmmFile.getAbsolutePath()+train.hapHMM.getName()+"_bef_"+overallIndex), dataC );
    int i_overall=0;
   
    for(int i=0; i<numIt.length; i++){
           train.initialise(numIt[i], Constants.precision(),expanded);
           for(int i1=0; i1<train.maxIt;i1++){
               if(Constants.printPlots() &&  i1==train.maxIt-1 && Constants.plot>=1 && i==numIt.length-1){
                   CNVHap.cnvf.pack();
                   CNVHap.cnvf.setVisible(true);
                   train.firePropertyChange("setToPlot", null, 2);
              //   DickFormat.cnvf.setToPlot(2);
               }
               else if(Constants.printPlots() &&  i1==train.maxIt-2 && Constants.plot>=1&& i==numIt.length-1){
                   CNVHap.cnvf.pack();
                   CNVHap.cnvf.setVisible(true);
                   train.firePropertyChange("setToPlot", null, 1);
//                   DickFormat.cnvf.setToPlot(1);
               }
              
             //  boolean isLast =  
                   boolean updateEmissions = i_overall==Constants.updateEmissionsIt();
            	   train.train(i_overall, updateEmissions, do_update);
            	  
            		   if(updateEmissions){
                 		  DistributionCollection.dc = null;
                 	   }
            	   i_overall++;
              // Logger.global.info("writing");
            //   if(isLast)DickFormat.cnvf.write(i1);
              // if(true){
               //    System.exit(0);
              // }
             //  if(Constants.plot()>1 || (isLast && Constants.plot()>0) )train.doPlots(isLast);
             //  if(Constants.plot()>0) DickFormat.cnvf.update();
              // if(isLast) break;
               
          }
           train.firePropertyChange("done", null, 1);
        if(i<numIt.length-1){
            logger.info("switching hmm ");;
         if(Constants.modifyFrac3!=null)   Emiss.switchHWE();
           if(Constants.expand().length>0 && Constants.expand()[0]>0){
        	//   Constants.m
        	  expanded = true;
        	   train.expandHMM(Constants.expand().length>1);
        	 
        	   //note this will change inversion code behaviour!
//        	   train.switchHMM();
           }
           else  train.switchHMM();
          
        }
    }
   
   
                
        
         // if(swtch){
         //     logger.info("arranging according to pedigree");;
          //    dataC.arrangeDataAccordingToPedigree();
//              this.updateData();
        //  }
          
   }
 int index_count =0;
   public void trainGraphical(final TrainingElement train, boolean swtch, boolean expanded) throws Exception{
      
          logger.info("before simple emissions");
  
   // int[] numIt = Constants.numIt();
          train.initialise(100, Constants.precision(),expanded);
   /*CNVHap.cnvf.jbutton_it.addActionListener(new ActionListener(){

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		try {
			train.train(index_count);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		index_count++;
	}
	   
   });
   CNVHap.cnvf.jbutton_plot1.addActionListener(new ActionListener(){

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			CNVHap.cnvf.pack();
            CNVHap.cnvf.setVisible(true);
            train.firePropertyChange("setToPlot", null, 1);
		}
		   
	   });
   CNVHap.cnvf.jbutton_plot2.addActionListener(new ActionListener(){

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			CNVHap.cnvf.pack();
           CNVHap.cnvf.setVisible(true);
           train.firePropertyChange("setToPlot", null, 2);
		}
		   
	   });
   CNVHap.cnvf.jbutton_sw.addActionListener(new ActionListener(){

		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			if(Constants.expand().length>0 && Constants.expand()[0]>0){
	        	//   Constants.m
	        	
	        	   train.expandHMM(Constants.expand().length>1);
	        	  
	        	   //note this will change inversion code behaviour!
//	        	   train.switchHMM();
	           }
	           else  train.switchHMM();
			 try {
				train.initialise(100, Constants.precision());
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		   
	   });*/
   
  
   
   
                
        
         // if(swtch){
         //     logger.info("arranging according to pedigree");;
          //    dataC.arrangeDataAccordingToPedigree();
//              this.updateData();
        //  }
          
   }
   
   boolean tofinish = false;
  
   public Callable train(final Iterator<TrainingElement> hmm) throws Exception{
      return new Callable(){
           public Object call(){
               try{
           boolean swtch = (Constants.sampleWithPedigree() && !Constants.trainWithPedigree());
           for(int j=0; hmm.hasNext(); j++){
        	   tofinish = false;
               TrainingElement train = hmm.next();
               train.modifyModelWithData();
            /*  if(Constants.plot>1) CNVHap.cnvf.jbutton_finish.addActionListener(new ActionListener(){

					@Override
					public void actionPerformed(ActionEvent arg0) {
						tofinish = true;
						
					}
       			   
       		   });*/
               if(Constants.plot>=2 && false){
            	   trainGraphical(train,swtch, false);
            	  // inner: while(!tofinish){
            	//	 Thread.sleep(100);
            		  
            	 //  }
               }
               else train(train, swtch);
               if( Constants.writeHMM()>=1)  writeHMM(train.hmm,
                       new File(hmmFile.getAbsolutePath()), dataC );
          
               if(Constants.sample() && ! Constants.keepBest() ){
                   logger.info("sampling from HMM");
                   monitor.sampleFromHMM(train.hmm, Constants.noSamples(),j, train.data1);
               }
               if(Constants.sample() && Constants.keepBest()){
                   if(bestSc==null || train.logProb > bestSc){
                   best = train.hmm;
                   bestSc = train.logProb;
                   }
               }
               if(Constants.saveModel()){
            	   try{
            	   CompressHMM hmmC = new CompressHMM(new File(dir, "hmm"), (HaplotypeHMM)train.hapHMM, (DataCollection)dataC);
            	   hmmC.run();
            	   {File out = new File(dir, "hmm.out");
            	   ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(out));
            	   os.writeObject(train.hapHMM);
            	   os.close();
            	   }
            	   File out1 = new File(dir, "dc.out");
            	   ObjectOutputStream os1 = new ObjectOutputStream(new FileOutputStream(out1));
            	   os1.writeObject(dataC.dc);
            	   os1.close();
            	   }catch(Exception exc){
            		   exc.printStackTrace();
            	   }
               }
           }
           if(Constants.keepBest() && Constants.sample()){
               logger.info("sampling from HMM");
               monitor.sampleFromHMM(best, Constants.noSamples(), 0, null);
           }
           }catch(Exception exc){
               exc.printStackTrace();
           }
           return null;
           }
       };
   }
   
   
   
  




  /*  public void trainIncr( int numIt1, int numIt2, int numFounders, double u) throws Exception{
           int numToAdd = 1;
         
           for(int i=0; i<this.modelLength; i++){
                  TrainingElement train = new TrainingElement(i);
                 // train.hapHMM.setTrain(true, false)
                  double res = Double.NEGATIVE_INFINITY;
                  while(true){
                           train.train(numIt1, 1,true);
                        //  if(numRepeats==1) Fastphase.this.writeHMM(train.hmm);
                          if(//converged(train.hapHMM, 0.01, numToAdd) || 
                                  train.hapHMM.modelLength()-1 >= numFounders) break;
                          else{
                              if(train.hapHMM instanceof FreeHaplotypeHMM) train.hapHMM.fix();
                              train.addState(numToAdd, init);
                             // if(train.hapHMM.modelLength()>4) train.hapHMM.setTrain(true, true);
                          }
                  }
                  
                  logger.info("mid point");
                 // monitor.monitor(this);
                  //res = train.train(numIt1, 0.1, 1, true, true);
                  if( train.switchHMM()) train.train(numIt1, 1,  true);
               //   models[i] = train.hmm;
               //   score[i] = res;
           }
          // logger.info("finished training "+Arrays.asList(score));
        //  outp.close();
       }*/
    static int pseudoCount = 0;
public static Set<Integer> marks =null;
       /** class for training a single hmm */
      public class TrainingElement{
           public MarkovModel hapHMM ;
           public CompoundMarkovModel swHMM;
           public CompoundMarkovModel swHMM2;
           public HaplotypeHMM countHMM;
           public List<CompoundMarkovModel>  hmm;
          //final public SimpleExtendedDistribution[] probHomoIsHemi;
           List<EmissionState>[] data1;
          // double[][] weights;
           PseudoIterator[] pseudoIterator;
           int num, trainingRounds, maxIt;
           double prec;
           BaumWelchTrainer[] bwt;
          // IndividualPlot[][] plot;
          // HaplotypePanel[] panel;
           public void makeBWT(){
        	
               String[] dc = dataC.getUnderlyingDataSets();
              // plot = new IndividualPlot[bwt.length][dc.length];
               //panel = new HaplotypePanel[bwt.length];
             if(Constants.plot()>=1) {
            	 CNVHap.cnvf.clearTabs();
            	 for(Iterator<PropertyChangeListener> it = 
            		 Arrays.asList(this.pcs.getPropertyChangeListeners()).iterator(); it.hasNext();){
            	 this.pcs.removePropertyChangeListener(it.next());
            	 }
            	 for(int i=0; i<bwt.length; i++){
            		 if(bwt[i]!=null){
            		 for(Iterator<PropertyChangeListener> it = 
                		 Arrays.asList(bwt[i].pcs.getPropertyChangeListeners()).iterator(); it.hasNext();){
            			 bwt[i].pcs.removePropertyChangeListener(it.next());
                	 }
            		 }
            	 }
             }
             
               for(int i=0; i<bwt.length; i++){
                   if(modelLengths[i] && data1[i]!=null && data1[i].size()>0){
                       bwt[i] = new BaumWelchTrainer(hmm.get(i),  data1[i].toArray(new EmissionState[0]), Constants.weights__, 
                    		   (dataC).getUnderlyingDataSets(), ((DataCollection)dataC).weights() 
                    		 );
                       if(Constants.plot()>0){
                           List<Integer> loc = dataC.loc;
                           List<String> snpid = dataC.snpid;
                          
                          /* Set<Set<Short>> indices =  bwt[i].fillDataIndices();
                           for(Iterator<Set<Short>> it = indices.iterator(); it.hasNext(); ){
                        	  Set<Short>data_index =  it.next();
                               StringBuffer name = new StringBuffer();
                               for(Iterator<Short> it1 = data_index.iterator(); it1.hasNext();){
                            	   name.append(dc[it1.next()]+" "+(i+1)+"_");
                               }
                             
                                   //  CNVHap.cnvf.addTab(panel.getName()+"original",panel1);
                                     // else DickFormat.cnvl.add(panel);
                               
                           }*/
                           if(true){
                        	   if(bwt[i].trainDists || Constants.plot(i)){
                        		  
                            	   Set<Short> s = bwt[i].getDataIndices(true);
                            	   
                            	   Iterator<Short> it1 = s.iterator();
                            	   Integer[] di = new Integer[s.size()];
                            	   for(int ii=0; ii<di.length; ii++){
                            		   di[ii] = (int) it1.next();
                            	   }
                            	   List<AssociationCalculator>[][]  ac = 
                            		   Constants.calcAssoc? 
                            		   ((DataCollection)dataC).getArmitageCalculator() : null;//new ArmitageCalculator((MergedDataCollection)dataC);
                              	 
                            	   IndividualPlot plot=null;
                            	  
                            	   if(true || DistributionCollection.dc!=null){
                                   plot =
                                	   new IndividualPlot(bwt[i], dc, loc, snpid, dataC.intensityProbes(),Constants.experiment(), dir, i+1, (short)-1, 
                                			   false);
                                   this.bwt[i].addPropertyChangeListener(plot);
                                   this.addPropertyChangeListener(plot);
                                  
                                  
                                  if(CNVHap.cnvf!=null) 
                                       CNVHap.cnvf.addTab("Data",plot);
                                       
                            	  // }
                            	   }
                            	   
                                 //  else DickFormat.cnvl.add(plot);
                            	   List<SignificancePlot> l = new ArrayList<SignificancePlot>();
                            	   
                            	  // SignificancePlot.splots = l;
                                  if(ac!=null){
                                	//  for(int kk=0; kk<ac.length; kk++){
                                		  if(ac[i]!=null && ac[i].length>0){
                                			
                                			  for(int k=0; k<ac[i].length; k++){
                                				  //if(ac[i][k]!=null)
                                				  for(int j=0; j<ac[i][k].size(); j++){
                                					  AssociationCalculator acl = ac[i][k].get(j);
                                					 acl.setModel(bwt[i].hmm);
                                				  
                                				  
                                	// ac[i].get(k).
                                	  SignificancePlot splot= new SignificancePlot(bwt[i], acl,loc,snpid,acl.dc1.name,dir, l.size()==0,k);
                                	
                                	  
                                	  if(l.size()>0){
                                		l.get(0).addPropertyChangeListener(splot);  
                                	  }
                                	  l.add(splot);
                                				  }
                                	  
                                			  }
                                			  plot.setSignficancePlot(l);
                                	if(Constants.plotDistributions()
                                        			  //&& ((AssociationCalculator)ac[i][k].get(j) instanceof LinearRegressionCalculator)
                                        			  ){
                                		 for(int j=0; j<ac[i][0].size(); j++){
                                		  AssociationCalculator acl = ac[i][0].get(j);
                                        		  for(int kj=0; kj<((LinearRegressionCalculator)acl).width; kj++){
                                       		   IndividualPlot plot_dist =
                                           	   new IndividualPlot(bwt[i], dc, loc, snpid, dataC.intensityProbes(), "allDistr_"+acl.getPheno(kj)
                                           			   , dir, i+1,(short) kj,  Constants.plotMerged());
                                       		    plot_dist.setSignficancePlot(l);
                                              		this.bwt[i].addPropertyChangeListener(plot_dist);
                                              		this.addPropertyChangeListener(plot_dist);
                                              		CNVHap.cnvf.addTab("DataDist_"+acl.types_all.get(kj),plot_dist);
                                        		  }
                                		 }
                                     }
                                	
                                	
                                	  if(CNVHap.cnvf!=null) {
                                		  for(int ii=0; ii<l.size(); ii++){
                                          CNVHap.cnvf.addTab("Association_"+(i+1)+"_"+l.get(ii).getName(),l.get(ii));
                                          this.bwt[i].addPropertyChangeListener(l.get(ii));
                                          this.addPropertyChangeListener(l.get(ii));
                                         
                                		  }
                               	   		}
                                	  }
                                	  }
                                 // }
                                  HaplotypePanel1 panel = null;
                                  if(Constants.showHaps()){ 
                                	//  short[] di = new short[data_index.size()];
                                	//  Iterator<Short> it1 = data_index.iterator();
                                	//  for(int kk=0; kk<di.length; kk++){
                                	//	  di[kk] = it1.next();
                                	//  }
                               	   try{
                                        panel = new HaplotypePanel1(bwt[i], dc, loc, dataC.snpid, dataC.alleleA, dataC.alleleB, "haplotypes", i+1,  dir, this.hapHMM, i+1);
                               	   }
                                     //  HaplotypePanel1 panel1 = new HaplotypePanel1((DataCollection) dataC);
                                       catch(Exception exc){
                                       	exc.printStackTrace();
                                       }
                                       bwt[i].addPropertyChangeListener(panel);
                                       this.addPropertyChangeListener(panel);
                                       //if(Constants.plot()>1)
                                         if(CNVHap.cnvf!=null) 
                                              CNVHap.cnvf.addTab(panel.getName(),panel);
                                         if(l!=null) {
//                                    		 panel.addPropertyChangeListener(splot);
                                    		 for(int ii=0; ii<l.size(); ii++){
                                    		 if(ii==0) {
                                    			 l.get(ii).addPropertyChangeListener(panel);
                                    		//	 if(plot!=null) l.get(ii).addPropertyChangeListener(plot);
                                    		 }
                                    		
                                    		 }
                                    		if(plot!=null && l.size()>0) plot.setSignficancePlot(l.get(0));
                                    		
                                    	 }
                                  }
                                 if(Constants.calcHWE) {
                                	  HWECalculator[] hwe = ((DataCollection)dataC).getHWECalculator();
                                	 // for(int kk=0; kk<hwe.length; kk++){
                                		  if(hwe[i]!=null){
                                	  HWEPlot splot= new HWEPlot(bwt[i].hmm, hwe[i],loc,snpid,"all",dir);
		                                 	 if(panel!=null) {
		//                                 		 panel.addPropertyChangeListener(splot);
		                                 		 splot.addPropertyChangeListener(panel);
		                                 		if(plot!=null) splot.addPropertyChangeListener(plot);
		                                 	 }
		                                 	  if(CNVHap.cnvf!=null) {
		                                           CNVHap.cnvf.addTab("HWE "+(i+1),splot);
		                                           this.bwt[i].addPropertyChangeListener(splot);
		                                          // plot.addPropertyChangeListener(splot);
		                                           this.addPropertyChangeListener(splot);
		                                	   		}
                                		  }
                                	 // }
                                  }
                               }
                           }
                           /*else{
                           for(short data_index=0; data_index<dc.length;data_index++ ){
                        	   if(Constants.plot(data_index)){
                               String name = dc[data_index]+" "+(i+1);
                               Integer[] di = new Integer[]{(int) data_index};
                               if(bwt[i].trainDists){
                            	   Set<Short> s = bwt[i].getDataIndices(true);
                            	   if(s.contains(data_index)){
                                   IndividualPlot plot = new IndividualPlot(bwt[i], dataC.getUnderlyingDataSets(),loc, snpid, Arrays.asList(di), name, dir);
                                 
                                  
                                 
                                  
                                   this.bwt[i].addPropertyChangeListener(plot);
                                   this.addPropertyChangeListener(plot);
                                  if(CNVHap.cnvf!=null) 
                                       CNVHap.cnvf.addTab(name,plot);
                            	   }
                                 //  else DickFormat.cnvl.add(plot);
                                  
                               }
                        	   }
                               
                           }
                           }*/
                         
                         
                           //plot[i].plot();
                       }
                       if(data1[i].size()==0) throw new RuntimeException("!!!");
                   }
                   }
               Logger.global.info(" mem is after BWT  "+Runtime.getRuntime().freeMemory());
               if(Constants.showHMM && Constants.plot>=1)
               {
                   HMMPanel pan = CNVHap.cnvf.addHMMTab((FreeHaplotypeHMM)hapHMM, dataC.getLocations(), dataC.alleleA, dataC.alleleB(), dataC.size(), dataC.pheno(),dir);
                   this.addPropertyChangeListener(pan);
                   for(int i=0; i<bwt.length; i++){
                	   if(bwt[i]!=null){
                		   this.bwt[i].addPropertyChangeListener(pan);
                	   }
                   }
                   if(Constants.globalTrans()){
                   RateMatrixPanel pn = CNVHap.cnvf.addRateMatrixTab((FreeHaplotypeHMM)hapHMM, dir);
                   this.addPropertyChangeListener(pn);
                   for(int i=0; i<bwt.length; i++){
                	   if(bwt[i]!=null){
                		   this.bwt[i].addPropertyChangeListener(pn);
                	   }
                   }
                   }
               }
               if(Constants.plot()>=1){
                   CNVHap.cnvf.pack();
                   CNVHap.cnvf.setVisible(true);
               }
               
           }
           TrainingElement(MarkovModel haphmm){
               this.pcs = new PropertyChangeSupport(this);
                   Class[] clazz = new Class[] { ExponentialTransitionProbs.class};
                  HaplotypeHMM hmm_i = null;
                  if(modelSwitch && false){
                   this.swHMM = new PairMarkovModel(new MarkovModel[] { 
                           hmm_i},
                               new int[] {0}, PairEmissionState.class, true,false);
                   this.swHMM2 = new PairMarkovModel(new MarkovModel[] { 
                           hmm_i},
                               new int[] {0,0}, PairEmissionState.class, true, false);
                  }
                setHMM(haphmm, 0);
               bwt = new BaumWelchTrainer[modelLengths.length];
               makeBWT();
           }
       
       public void modifyModelWithData(){
    	    if(Constants.useDataAsModel()!=null){
     	  List<String> datamodel = Constants.useDataAsModel() ==null ? null : Arrays.asList(Constants.useDataAsModel());
     		((HaplotypeHMM) hapHMM).modify(
     			
     				DataCollection.datC, datamodel);
     		  this.makeGenotypeHMM(1);
              this.makeBWT();
     	  }
    	    }
           public boolean expandHMM(boolean swtch) {
        	   int[] expand = Constants.expand();
        	   if(swtch && !Constants.newTrans())   ((FreeHaplotypeHMM)this.hapHMM).swtch();
        	  this.hapHMM = new ExpandedHMM((FreeHaplotypeHMM)this.hapHMM, expand, dataC);
        	//modifyModelWithData();
        	  if(Constants.useHMMFile(1)!=null){
        			 try{
        				 if(Constants.useHMMFile(1).endsWith(".zip")){
        			 ZipFile hmmFile = new ZipFile(new File(Constants.useHMMFile(1)));
        			((HaplotypeHMM) hapHMM).modify(hmmFile, Constants.transferHMM(), dataC.loc);
        			 hmmFile.close();
        				 }
        				 else{
        					 ObjectInputStream ois = new ObjectInputStream(new FileInputStream(new File(Constants.useHMMFile(1)))); 
        						hapHMM = (FreeHaplotypeHMM) ois.readObject();
        						hapHMM.initialiseCounts();
        						System.err.println("done reading");
        				 }
        			 }catch(Exception exc){
        				 exc.printStackTrace();
        			 }
        		 }
        	 //  ((FreeHaplotypeHMM)this.hapHMM).expand(swtch);
        	   try{
        	  if(Constants.CHECK) hapHMM.validate(hapHMM.noSnps);
               //    hapHMM.updateEmissionStateSpaceDist(1);
               //    hmmFile =  new File(dir, "results_"+hapHMM.getClass().toString()+"_"+(hapHMM.modelLength()-1)+".txt");
                   this.makeGenotypeHMM(1);
                   this.makeBWT();
        	   }catch(Exception exc){
        		   exc.printStackTrace();
        		   return false;
        	   }
        	   return true;
           }
       
           public boolean switchHMM() {
               try{
            
            
                  logger.info("SWITCHING!!!!");
                  if( hapHMM instanceof FreeHaplotypeHMM){
                      this.hapHMM =  (MarkovModel) hapHMM.clone(true);
//                          new FreeHaplotypeHMM((FreeHaplotypeHMM)hapHMM) ;
                  }
                  else if( hapHMM instanceof PairMarkovModel){
                      hapHMM.validate(hapHMM.noSnps);
                   /*   boolean max = true;
                      
                     MarkovModel[] mm = ((PairMarkovModel)hapHMM).getMemberModels();
                     for(int ik=0; ik<1; ik++){
                         max= max && ( (FreeHaplotypeHMM)mm[ik]).trans.transProbs[2] instanceof FreeTransitionProbs1;
                      } 
                     if(!max){
                     for(int ik=0; ik<1; ik++){
                        mm[ik] =(MarkovModel)( (FreeHaplotypeHMM)mm[ik]).clone(true);
                     }
                     this.hapHMM =   new PairMarkovModel(mm, new int[] {0,1}, AlleleCopyPairEmissionState.class,
                             Emiss.stateSpace() , false       
                             //)
                             );
                     }
                     else{*/
//                      CachedHMM hmmC = new CachedHMM((CompoundMarkovModel)hapHMM);
                         hapHMM = new FreeHaplotypeHMM(hapHMM);
                     //}
                  }
                  hapHMM.validate(hapHMM.noSnps);
              //    hapHMM.updateEmissionStateSpaceDist(1);
              //    hmmFile =  new File(dir, "results_"+hapHMM.getClass().toString()+"_"+(hapHMM.modelLength()-1)+".txt");
                  this.makeGenotypeHMM(1);
                  this.makeBWT();
                 // if(Constants.writeHMM()) writeHMM(hmm, new File(hmmFile.getAbsolutePath()+"_before"),  dataC );
              }catch(Exception exc){
                  exc.printStackTrace();
              }
              return true;
           }
           private void  setHMM(MarkovModel hmm, int phase) {
               this.hapHMM = hmm;
              // if(((EmissionState)hmm.getState(1)).isFixed()){
              //     modify(hmm, dataC);
             //  }
/*               if(hmm instanceof FreeHaplotypeHMM && Constants.modifyWithData()>0){
                     try{
                   modify((FreeHaplotypeHMM)hmm, dataC, null);
                     }catch(Exception exc){
                         exc.printStackTrace();
                         System.exit(0);
                     }
               }
               
                   else if(hmm instanceof CompoundMarkovModel){
                       MarkovModel mm = ((CompoundMarkovModel)hmm).getMemberModels()[1];
                     //  if(!mm.getName().startsWith("allele")) throw new RuntimeException("!! "+mm.getName());
                       modify((FreeHaplotypeHMM)mm, dataC, 1);
                   }*/
              this.makeGenotypeHMM(phase);
             //  if(bwt!this.bwt.hmm!=null) this.bwt.hmm = this.hmm;
            }
           public void modify(FreeHaplotypeHMM hmm, DataCollection data, Integer index){
              // if(!hmm.modifyWithData()) return;
               Logger.global.info("modifying!!");
             EmissionState maf = null;;//data.calculateMaf1();
             EmissionStateSpace emStSp = maf.getEmissionStateSpace();
             if(index!=null){
                 maf = ((PairEmissionState)maf).getMemberStates(true)[index];
             }
             if(Constants.modifyWithData()==2){
            	 double numF = 0;
            	 inner: for(int i=0; i<hmm.noSnps; i++){
            		 Dirichlet dir = new Dirichlet(maf.getEmiss(i), Constants.u_global[0]);
                   for(int j=1; j<hmm.modelLength(); j++){
                       HaplotypeEmissionState st = (HaplotypeEmissionState) hmm.getState(j);
                      
                    	   st.emissions[i].change(dir);
                   }
            	 }
                  
             }
          //   CompoundEmissionStateSpace ces = ((CompoundEmissionStateSpace)maf.getEmissionStateSpace());
          outer:  for(int kj=0; kj < emStSp.copyNumber.size(); kj++){
            double numF = 0;
             List<EmissionState> emStates = new ArrayList<EmissionState>();
             for(int j=1; j<hmm.modelLength(); j++){
                 EmissionState st = (EmissionState) hmm.getState(j);
                 Integer fixed = st.getFixedInteger(0);
                 if(fixed==null)  fixed = st.getBestIndex(0);
                 Integer cn = emStSp.getCN(fixed);
                 if(cn==emStSp.copyNumber.get(kj)){
                     numF++;
                     emStates.add(st);
                 }
             }
             int[] set = emStSp.getGenoForCopyNo(emStSp.copyNumber.get(kj));
             if(set.length==1 || emStSp.copyNumber.get(kj)!=1) continue outer;
             //int[] set = new int[] {0,1};
             double avg = 1.0 / (double)set.length;//hmm.getEmissionStateSpace().size();
            
             inner: for(int i=0; i<hmm.noSnps; i++){
                 Integer fixed = maf.getFixedInteger(i);
                 if(fixed!=null){
                     for(int j=1; j<hmm.modelLength(); j++){
                         EmissionState st = (EmissionState) hmm.getState(j);
                         st.setFixedIndex(i, fixed);
                     }
                 }
                 else{
                  double[] probs = maf.getEmiss(i);
                  if(Constants.CHECK && Math.abs(Constants.sum(probs)-1.0) >0.01) throw new RuntimeException("!!");
                  int[] count = new int[probs.length];
                  int sum=0;
                  SortedMap<Double, Integer> lessThan = new TreeMap<Double, Integer>();
                  double summ = 0;
                  for(int kk=0; kk<set.length; kk++){
                      summ+=probs[set[kk]];
                  }
                  if(summ==0) continue inner;
                  for(int kk=0; kk<set.length; kk++){
                      int k = set[kk];
                      double prob = (probs[k]/summ);
                     double raw = prob*numF;
                    
                     int cnt = (int)
                         Math.max(Constants.modifyWithData(), 
                         prob < avg ?  Math.ceil(prob*numF) :  Math.floor(prob*numF));
                     if(cnt < raw){
                         lessThan.put(raw-cnt, k);
                     }
                     count[k] = cnt;
                     sum+=cnt;
                 }
                 while(sum < numF){
                     try{
                     int k = lessThan.remove(lessThan.lastKey());
                     count[k]++;
                     sum++;
                     }catch(Exception exc){
                         System.err.println(index);
                         exc.printStackTrace();
                         System.exit(0);
                     }
                 }
                 while(sum>numF){
                     int maxIndex = Constants.getMax(count);
                     count[maxIndex]--;
                     sum--;
                 }
                 int j=0;
                 for(int kk=0; kk<set.length; kk++){
                     int k  = set[kk];
                     for(int kk1=0; kk1<count[k]; kk1++){
                        EmissionState st = (EmissionState) emStates.get(j);
                         st.setFixedIndex(i, k);
                         j++;
                     }
                 }
                 }
           //      System.err.println("done");
             }
            }
           //  if(Constants.numRep()>1) ((FreeHaplotypeHMM)hmm).reorderStates(true, true);
           }
         private CompoundMarkovModel makeBasicModel(int no_cop, MarkovModel toCopy){
             int[] no_copies = new int[no_cop];
             MarkovModel[] m = new MarkovModel[] {toCopy};
             for(int i=0; i<no_copies.length; i++){
               //  m[i] = toCopy;
                 no_copies[i] = 0;
             }
               return  new PairMarkovModel(m, no_copies, PairEmissionState.class,   true, false);
             
            
         }
         
        
          
         private MarkovModel makeHMM(int no_cop, MarkovModel hapHMM){
             CompoundMarkovModel result = null;
             if(no_cop==1 || no_cop==2 || true){
                 result =  makeBasicModel(no_cop, hapHMM);
             }
             else {
                 /*if(no_cop==3){
                     CompoundMarkovModel mm = hmm.get(1);
                     MarkovModel[] m = new MarkovModel[] {mm.unwrapModel(), swHMM};
                     result = new PairMarkovModel(m, new int[] {0,1}, HalfTrioEmissionState.class,  Emiss.getEmissionStateSpace(2), true);
                 }
                 else if(no_cop==4){
                     throw new RuntimeException("should not ask for this one as we can split into two models");
                 }*/
               /*  else if(no_cop==5){
                     CompoundMarkovModel inner = new PairMarkovModel(
                             new MarkovModel[] {hmm.get(0).unwrapModel(), swHMM}, new int[] {0,1}, HalfTrioEmissionState.class, Emiss.getEmissionStateSpace(2), true);
                     MarkovModel[] m = new MarkovModel[] {hmm.get(2).unwrapModel(),inner};
                     result = new PairMarkovModel(m, new int[] {0,1}, TrioEmissionState2.class, Emiss.getEmissionStateSpace(4), true);
                 }
                 else if(no_cop==6){
                     CompoundMarkovModel mm = hmm.get(1);
                     if(mm==null) throw new RuntimeException("is null!!");
                   //  CompoundMarkovModel swHMM2 = makeBasicModel(2, swHMM, hemiHomo);
                     MarkovModel[] m = new MarkovModel[] {mm.unwrapModel(), swHMM2};
                    // if(emissionStateSpace[5]==null) throw new RuntimeException("!!");
                     result =// new PairMarkovModel(m, new int[] {0,0,1}, TrioEmissionState2.class,  emissionStateSpace[5]);
                         new PairMarkovModel(m, new int[] {0,0,1}, TrioEmissionState2.class, Emiss.getEmissionStateSpace(5), true);
                 }*/
                throw new RuntimeException("!! "+no_cop);
               
             }
             if(Constants.fast() && result.noCopies()>1 && hapHMM.modelLength()>2 && !result.allOneLength()) result = new CollapsedHMM(result);
            // return result;
             Logger.global.info(" mem is "+Runtime.getRuntime().freeMemory());
             if(Constants.cache()){
                 
                 MarkovModel res1 =  new CachedHMM(result);
                 Logger.global.info(" mem is "+Runtime.getRuntime().freeMemory());
                 return res1;
             }
             else return result;
         }
           
           private void makeGenotypeHMM(int phase){
        	   Fastphase.model = hapHMM;
               this.hmm =  new ArrayList<CompoundMarkovModel>(){
                   {
                       for(int i=0; i<modelLengths.length; i++){
                           add(null);
                       }
                   }
                   public int size(){
                       return modelLengths.length;
                   }
                   public CompoundMarkovModel get(int i){
                       if(i>=modelLengths.length) return null;
                       if(!modelLengths[i]) return null;
                       if(super.get(i)==null){ 
                           set(i, (CompoundMarkovModel) makeHMM(i+1, hapHMM));
                       }
                       return super.get(i);
                   }
               };
               this.data1 = new List[modelLengths.length];
            //   this.weights = new double[modelLengths.length][];
               for(int i=0; i<hmm.size(); i++){
                   if(modelLengths[i]){
                       //hmm[i] = (CachedHMM) makeHMM(i+1);
                       if(data[i]!=null){
                           data1[i] = new ArrayList<EmissionState>();
                          // weights[i] = new double[data[i].size()];
                          for(int ik=0; ik < data[i].size(); ik++){
                               String key = data[i].get(ik);
                               data1[i].add( dataC.dataL(key));
                               //weights[i][ik] = dataC.weight(key);
                           }
                          
                         
                       }
                   }
               }
           }
           
           public void setPseudoWeights(double[][] pseudo){
        	   if( hapHMM instanceof PairMarkovModel){
                   boolean max = true;
                   
                  MarkovModel[] mm = ((PairMarkovModel)hapHMM).getMemberModels();
                  mm[0].setPseudoCountWeights(pseudo);
                  mm[1].setPseudoCountWeights(new double[][] {pseudo[1], pseudo[1], pseudo[2], pseudo[2]});
               }
               else{
                   this.hapHMM.setPseudoCountWeights(pseudo);
               }
        	   if(swHMM!=null){
                   double[] ps = new double[] {1e10, 0.001,  Constants.u_exp(), 1e10};
                       this.swHMM.setPseudoCountWeights(new double[][] {ps, ps,ps,ps});
                      this.swHMM.initialiseCounts();
               }
           }
           
           public void initialise(
        		  int round){
              
               this.hapHMM.initialiseCounts();
               for(int ik=0; ik<hmm.size(); ik++){
                   if(hmm.get(ik)!=null){
                       hmm.get(ik).initialiseCounts();
                   }
               }
               if(swHMM!=null){
                  
                      this.swHMM.initialiseCounts();
               }
              ((DataCollection) dataC).initialisationStep();
              if(DistributionCollection.dc!=null) DistributionCollection.dc.initialise();
              // if(dataC instanceof LikelihoodDataCollection){
                //   ((LikelihoodDataCollection)dataC).initialisationStep();
               //}
           
           }
          public double logProb = Double.NEGATIVE_INFINITY;
         /* public void doPlots(boolean last){
              for(int ik=bwt.length-1; ik>=0 ;ik--){
                  for(int ij=0; ij<plot[ik].length; ij++){
                  if(this.plot[ik][ij]!=null){
                      plot[ik][ij].update();
                      if(last && false){
                          plot[ik][ij].writeToFile();
                      }
                  }
                  }
              }
              Logger.global.info("here");
              // DickFormat.cnvf.pack();
          
          }*/
          
          public boolean  train(int i, boolean updateEmissions, boolean train) throws Exception{
        	 /* if(i==10){
        		
        		  MultipleDataRegression.switchSet();
        	  }*/
        	
              firePropertyChange("pre_exp", null,null);
             if(Constants.stop(i) && CNVHap.cnvf!=null){
            	 try{
            	 String[] poss = CNVHap.cnvf.getPoss();
            	 String sel = CNVHap.cnvf.currentTab();
            	 if(sel!=null){
String s = (String) JOptionPane.showInputDialog(null,"Show pane","", JOptionPane.PLAIN_MESSAGE, null,
		poss, sel);
CNVHap.cnvf.setTab(s);
             }
            	// int res = JOptionPane.(CNVHap.cnvf.jframe, "Continue?");
             }catch(Exception exc){
            	 exc.printStackTrace();
             }
             }
              double[][] pseudo_ = new double[pseudoIterator.length][];
          
              for(int k=0; k<pseudo_.length; k++){
            	  pseudo_[k] = pseudoIterator[k].next();
              }
              double[] pseudo = pseudo_[0];
              double[] pseudoG = pseudo_[3];
              double lp =   0;
              setPseudoWeights(pseudo_);
             initialise(i);
              List tasks  = new ArrayList();
              boolean isLast = i>maxIt ;
        //  for(int ik=bwt.length-1; ik>=0 ;ik--){
            	  for(int ik=0; ik<bwt.length ;ik++){
                  if(bwt[ik]!=null){
                     if(bwt[ik].hmm instanceof CachedHMM) ((CachedHMM)bwt[ik].hmm).refresh();
                   
                      lp+=bwt[ik].expectationStep(pseudo, isLast, updateEmissions);
                      
                      if(Constants.trainDists() && train){
                    	  if(DistributionCollection.dc!=null){
                    		  if(Constants.dcFile==null){
                    	  DistributionCollection.dc.maximisationStep(i, 
                    			  new double[] {pseudo[2], pseudo[3], pseudo[4]},
                    			  new double[] {pseudoG[2], pseudoG[3], pseudoG[4]},
                    		 tasks);
                    		  }
                    	  }
                    	  if(DataCollection.datC.bg()!=null) DataCollection.datC.bg().transferCountsToProbs(pseudo[4]);
                  //        bwt[ik].maximisationStep(i, pseudo[3], tasks); //trains the SkewNormal distributions
                        
                      }
                  }
              }
            //  isLast = isLast || (num>=1 && Math.abs(lp-logProb) < prec);
            if(train){
            try{
                  BaumWelchTrainer.involeTasks(tasks, true);
                  firePropertyChange("dist_maximisation", null,null);
              }catch(Exception exc){
                  exc.printStackTrace();
              }
            }
              tasks.clear();
             if(train) BaumWelchTrainer.maximisationStep(hapHMM, i, tasks);
          
           
              Logger.global.info("START ########################################");
              if(swHMM!=null && i >Constants.indexToTrainSWHMM() && train){
                      BaumWelchTrainer.maximisationStep(swHMM, i, tasks);
              }
              try{
                  if(train) BaumWelchTrainer.involeTasks(tasks, false);
                  firePropertyChange("hmm_maximisation", null, null);
//                  es.invokeAll(tasks);
              }catch(Exception exc){
                  exc.printStackTrace();
              }
             
              if(Constants.writeHMM()>1)  writeHMM(hmm,
                      new File(hmmFile.getAbsolutePath()+"."+overallIndex), dataC );
              overallIndex++;
              for(int ik=0; ik<bwt.length; ik++){
                  if(bwt[ik]!=null){
                      hmm.get(ik).refresh();
                  }
              }
              long t1 = System.currentTimeMillis();
          
              logger.info("log prob is "+lp+" at "+i+" with "+(this.hapHMM.modelLength()-1)+" states "+(t1-time)+" at "+Constants.print(pseudo));
              logger.info("log prob is "+lp+" at "+i+" with "+(this.hapHMM.modelLength()-1)+" states "+(t1-time)+" at "+Constants.print(pseudoG));
           //   this.hapHMM.
              if(marks!=null){
              for(Iterator<Integer> it = marks.iterator(); it.hasNext();){
            	  double rate = (((FreeHaplotypeHMM)this.hapHMM).trans.getRate(it.next()));
            	   logger.info("rate "+"\t "+String.format("%5.3g\t%5.3g mb", new Object[] {rate,(-1*Math.log(1-rate)/Constants.probCrossOverBetweenBP)/1e6}));
              }
              }
            //  time = t1;
            if(lp+1 < logProb && Constants.CHECK) {
                logger.info("warning: log prob should always be increasing ".toUpperCase()+lp +" should be > "+logProb);
            }
           // if( (swHMM==null || this.swHMM.converged()) && Math.abs(lp-logProb) < prec  )num++;
           // else num=0;
            logProb = lp;
            //initialise(i);
            if(i>=maxIt  ){// && max_next >=max){
              //  if(Constants.plot())doPlots(true, i, tasks);
               return true;
            }
            trainingRounds++;
            return false;
          }
          final PropertyChangeSupport pcs;
        //List<HMMPanel> hmmplots = new ArrayList<HMMPanel>();

      void addPropertyChangeListener(PropertyChangeListener arg0){
            this.pcs.addPropertyChangeListener(arg0);
        }

       void firePropertyChange(String propertyName, Object oldValue, Object newValue){
            pcs.firePropertyChange(propertyName, oldValue, newValue);
        }
          
          
        public void initialise( int maxIt, double prec, boolean expanded) throws Exception{
             pseudoIterator = new PseudoIterator[] {
            	
            	 Constants.pseudo__(expanded), Constants.pseudo1__, Constants.pseudo2__, Constants.pseudo3__
             };
             for(int i=0; i<pseudoIterator.length; i++){
            	 try{
            	 pseudoIterator[i].refresh();
            	 }catch(Exception exc){
            		 System.err.println("problem with "+i);
            		 exc.printStackTrace();
            		 System.exit(0);
            	 }
             }
             this.maxIt = maxIt;
             this.prec = prec;
               if( Constants.writeHMM()>1)  writeHMM(hmm,
                       new File(hmmFile.getAbsolutePath()+"_-1_"+overallIndex), dataC );
             // if(do_train) logger.info("training for  "+maxIt+" iterations");;
               if(maxIt==0) return;
               logProb = Double.NEGATIVE_INFINITY;
            trainingRounds =0;
              num=0;
    
           //    doPlots(true);
             
           }
           
       }
   
 
    
   public void updateData(){
       for(int i=0; i<data.length; i++){
           if(data[i]!=null){
               data[i].clear();
           }
       }
       for(Iterator<String> it = dataC.getKeys().iterator();it.hasNext(); ){
           String key = it.next();
           EmissionState dat =  dataC.dataL(key);
         //  System.err.println(dat.noCopies()-1);
           this.data[dat.noCop()-1].add(key);
       }
   }
    
    private List<String>[] data; //index is the no of copies
   final File hmmFile;
    File hmmOutputFile, phasedFile,  initialHMMFile;
   final boolean[] modelLengths;  //whether there is a model with this no copies;
  
   final File dir;
   boolean modelSwitch = false;
    int noIndiv=0;
    /** trainNum is the index (exclusive) of those used for training.  If this is null all are used in training and phasing,
     *  if this is not null, the first up to trainNum are used for training, and those after for testing */
    public Fastphase( DataCollection data, Sampler monitor, int numRep, File dir) throws Exception{
//Double[] phenV = data.dataL.get("21603").phenValue();
    //	HaplotypeEmissionState nxt = (HaplotypeEmissionState) data.dataL.get("NA12892");
    	if(data.dataL.size()==0) throw new RuntimeException("no elements in data collection "+data.name);
    	this.modelLength = numRep;
    String dcFile = Constants.dcFile();
    if(dcFile!=null){
    	ObjectInputStream ois = new ObjectInputStream(new FileInputStream(new File(Constants.dcFile()))); 
		AbstractDistributionCollection res = (AbstractDistributionCollection) ois.readObject();
    	data.setDC(res);
    }
            try{
            ConsoleHandler handler = new ConsoleHandler();
            FileHandler handlerF = new FileHandler(dir.getAbsolutePath()+"/stderr_"+Constants.fast(), false);
            logger.addHandler(handlerF);
            Formatter formatter = 
            new Formatter(){
                public String format(LogRecord record){
                    return record.getSourceClassName()+":\n"+record.getMessage()+"\n";
                }
            };
            handler.setFormatter(formatter);
            handlerF.setFormatter(formatter);
            }catch(Exception exc){
                exc.printStackTrace();
                System.exit(0);
            }
       
   // 	HaplotypeEmissionState hes = (HaplotypeEmissionState) data.dataL.get("82027");
        this.dir = dir;
        this.dataC = data;
        hmmFile =  new File(dir, "results_hmm"+//haphmm.getClass().toString()+"_" +
               ".txt");
        SortedSet<Integer> len = new TreeSet<Integer>();
        for(Iterator<EmissionState> it = dataC.dataLvalues(); it.hasNext();){
            EmissionState dat  = it.next();
            int noCopies = dat.noCop();
            len.add(noCopies);
            noIndiv+=noCopies;
        }
        System.err.println("size "+((DataCollection)dataC).size());
        modelLengths = new boolean[len.last()];
        this.data = new List[len.last()];
      //  this.emissionStateSpace = new List[len.last()];
        for(int i=0; i<modelLengths.length; i++){
            modelLengths[i] = len.contains(i+1);
            if(modelLengths[i] && i>=2) modelSwitch = true;
            if(modelLengths[i] && (i==2 || i==3 || i==4)){
                modelLengths[0] =true; 
            }
            if(modelLengths[i] && (i==2 || i==3 || i==4 || i==5)){
                modelLengths[1]=true;
            }
          /* if(modelLengths[i] && i==5){
               modelLengths[2]=true;
               modelLengths[0] = true;
               modelLengths[1] = true;
           }*/
        }
        for(int i=0; i<modelLengths.length; i++){
            if(modelLengths[i]){
                this.data[i] = new ArrayList<String>();
                
               
            }
        }
        this.updateData();
        StringBuffer sb_is = new StringBuffer();  StringBuffer sb_ds = new StringBuffer();
        StringBuffer sb_is1 = new StringBuffer(); 
        for(int i=0; i<data.length(); i++){
            sb_is.append("%8i ");
            sb_ds.append("%8.2g ");
        }
        this.sb_i = sb_is.toString();
        this.sb_d = sb_ds.toString();
        this.hmmOutputFile = new File(dir, "hmm_output");
        this.noSites = data.length();
        this.phasedFile = new File(dir, "phased.txt");
        this.monitor = monitor;
        this.initialHMMFile = new File(dir, "results_initial.txt");
    }
    
   
   
   
  
  


    private static boolean noNull(Comparable[] c){
        for(int i=0; i<c.length; i++){
            if(c[i]==null) return false;
        }
        return true;
    }
    
 
    
    final String sb_i;
    final String sb_d;
    
  
   
    
   
    public boolean equal(Comparable[] o1, Comparable[] o2){
           return Arrays.asList(o1).equals(Arrays.asList(o2));
       }
   
   
    final Sampler monitor;
    static Calendar cal = new GregorianCalendar();
  
   
public  static void writeHMM(List phmms, File textHMMOutput, DataCollection data) throws Exception{
   CompoundMarkovModel phmm =null;
   for(int i=phmms.size()-1; i>=0; i--){
     
       if(phmms.get(i)!=null){
          
          phmm  = ((CompoundMarkovModel)phmms.get(i));
          break;
       }
   }
   
    logger.info("WRITING HMM "+textHMMOutput.getAbsolutePath());
    textHMMOutput.mkdir();
   // PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(textHMMOutput)));
    //pw.println(cal.getTime());
    if(DistributionCollection.dc!=null){
    DistributionCollection.dc.print(textHMMOutput);
    }
    if(true){
    //	phmm.printProbR(pw);
    	HaplotypeHMM mm = (HaplotypeHMM)(phmm.getMemberModels()[0]);
    	if(mm instanceof FreeHaplotypeHMM)
    	((FreeHaplotypeHMM)mm).trans.printR(textHMMOutput);
    	else
    		((ExpandedHMM)mm).trans.printR(textHMMOutput);
    	/*if(Constants.writeHMM()<=1){
    		pw.close();
    		return;
    	}*/
    }
     if(data instanceof IlluminaRDataCollection){
        ((IlluminaRDataCollection)data).printDist(textHMMOutput);
    }
   
   /* if(data instanceof MergedDataCollection){
       DataC[] da = ( (MergedDataCollection)data).ldl;
       for(int i=0; i<da.length; i++){
           if(da[i] instanceof CCIlluminaDataCollection){
               ((CCIlluminaDataCollection)da[i]).printDist(pw);
             
           }
           else if(da[i] instanceof IlluminaRDataCollection){
               ((IlluminaRDataCollection)da[i]).printDist(pw);
             
           }
         
           pw.println(da[i].getClass());
         //  ((RatioDataCollection)da[i]).printDist(pw);
           //}
       }
    }*/
    //boolean restrict = false;
    StringBuffer sb = new StringBuffer();
    List<Integer> posi = new ArrayList<Integer>();
    List<Integer> cols = null;
    for(int i=0; i<data.loc.size(); i++){
        posi.add(i);
    }
    /*if(restrict){
        cols = new ArrayList<Integer>();
        for(int i=0; i<maf1.noSnps(); i++){
            if(maf_null1.get(i)!=0){
                cols.add(i);
            }
        }
    }*/
    Object[] num =subList(data.getLocations(), cols).toArray();
    Object[] name = data.snpid!=null ? subList(data.snpid, cols).toArray(): null;
    
   // EmissionState maf_i = data.maf;
   //if(maf_i instanceof PairEmissionState) maf_i = ((PairEmissionState)maf_i).getMemberStates(true)[1];
   //((HaplotypeEmissionState)maf_i).print(pw, "MAF            ");
   
   // Object[] mafs =   HaplotypeEmissionState.emissions(maf_i, cols);//subList(Arrays.asList(data.maf.emissions), cols).toArray();
   //Double[] maf = (Double[])mafs[0];
   //String[] maf_st = (String[])mafs[1];
    Object[] posi_i = subList(posi, cols).toArray();
   // Double[] maf = FindLines.getMaf(this.data, this.noSites);
    StringBuffer sb1= new StringBuffer();
    StringBuffer sb2= new StringBuffer();
   for(int i=0; i<data.loc.size(); i++){
       sb1.append("%8d ");
       sb.append("%8s ");
       sb2.append("%8.2g ");
    }
  // SimpleExtendedDistribution[] phh = ((PairMarkovModel)((WrappedModel)phmm).unwrapModel()).probHomoIsHemi;
   PrintWriter hmmpw = new PrintWriter(new BufferedWriter(new FileWriter(new File(textHMMOutput, "hmm.txt"))));
   //if(name!=null && name.length>0) hmmpw.println("name     "+String.format(sb.toString(), name));
   //if(posi_i.length>0)hmmpw.println("posi     "+String.format(sb.toString(), posi_i));
   //if(num.length>0)hmmpw.println("position "+String.format(sb1.toString(), num));
   //if(maf.length>0)pw.println("m_allele "+String.format(sb.toString(), maf_st));  
   //if(maf.length>0)pw.println("m_allele "+String.format(sb2.toString(), maf));  
  // if(maf1.size()>0)pw.println("m_allele1 "+String.format(sb1.toString(), subList(maf1, cols).toArray()));  
   //hmmpw.println(phmm.info());
  
   phmm.print(hmmpw, cols, data.noAllelles());
  /* if(phh!=null){
   SimpleExtendedDistribution[] homoHemi  =(SimpleExtendedDistribution[]) 
   subList(Arrays.asList(phh), cols).toArray(new SimpleExtendedDistribution[0]);
    Double[] homoHemiL = new Double[homoHemi.length];
    for(int i=0; i<homoHemi.length; i++){
    homoHemiL[i] = homoHemi[i].probs[0];
    if(homoHemiL[i]< 1e-3) homoHemiL[i] = 0.0;
}
   pw.println("homohemi  "+String.format(sb2.toString(), homoHemiL));
   }*/
   hmmpw.close();
   logger.info("finishing wrinting hmm");
   
   
}

public static List subList(List l, List<Integer>cols){
    if(cols==null) return l;
    List res = new ArrayList();
    for(int i=0; i<cols.size();i++){
        res.add(l.get(cols.get(i)));
    }
    return res;
}
public void sample() {
    // TODO Auto-generated method stub
    
}


    
}



