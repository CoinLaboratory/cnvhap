package lc1.dp.swing;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;

import lc1.dp.core.ProbMultivariate;
import lc1.dp.core.Sampler;
import lc1.dp.data.collection.HWECalculator;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.stats.Mixture;
import lc1.stats.ProbabilityDistribution;
import lc1.stats.SkewNormal;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.svg.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.AxisChangeEvent;
import org.jfree.chart.event.AxisChangeListener;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYZDataset;


public class HWEPlot extends JTabbedPane implements PropertyChangeListener{
    
    final List<Integer> index;
    final String[] names1;
    
   static  int width = 1200;
   static int height = 200;
   static Dimension dim2 = new Dimension(width, height);
   static final  Dimension dim = new Dimension(width, height);
   static final  Dimension dim_ = new Dimension(width, height*2);
   static final  Dimension dimG = new Dimension(800, 800);
  //static EmissionStateSpace[] stateEm=  new EmissionStateSpace[2];//null;
   CompoundEmissionStateSpace emStSp;
 //EmissionStateSpace emStSp1;
// final MarkovModel hmm;
   
  /* static{
       if(Constants.transMode(1)==null && Constants.saveStates()){
           stateEm[0] =   Emiss.getStateEmissionStateSpace(new int[] {Constants.numF(0)});
           stateEm[1] =   Emiss.getStateEmissionStateSpace(new int[] {Constants.numF(0), Constants.numF(0)});
       }
   }*/
 //  final double[] cn_count;
  // final int[] cn_alias;  //for distribution
    JComponent jpB;
    JComponent[] jB, jR;
   // JComponent jG; 
  
    JComponent[] chartB, chartR;
   // JComponent[][]chartDist;
    //JComponent[][]chartG;
    ChartPanel[] chartBNew, chartRNew;
    //ChartPanel[][]chartGNew;
    JComponent jpBS, jpRS;
   // JScrollPane jGS;
    public static  Font font16 = new Font("SansSerif", Font.PLAIN, 16);
    public static  Font font10 = new Font("SansSerif", Font.PLAIN, 10);
    public static  Font font4 = new Font("SansSerif", Font.PLAIN, 6);
    final int noSnps;
   // BaumWelchTrainer bwt;
    List<Integer> location;
    
  final int firstHalf;
    public XYSeriesCollection getRSeriesCollection(int k){
      XYSeriesCollection r = rdc.get(k);
      if(r==null){
          
        	   r = new XYSeriesCollection();
          
           //if(bwt.probDists[index].mvf!=null)
           //for(int ij=0; ij<bwt.probDists[index.get(i)].mvf.pairs.size(); ij++){
        	   if(k<firstHalf){
        		   r.addSeries( new XYSeries("HWE overall -log10(pval)"));
        	   }
        	   else{
        		   r.addSeries( new XYSeries("HWE between CN states"));
        	   for(int i=2; i<ac.type_len; i++){
               r.addSeries( new XYSeries("cn = "+(i-1)));
        	   }
        	   }
             // if(k==0) addRBSeries(bwt.probDists[index.get(i)].mvf.getCompoundName(ij));
          // }
           
            rdc.set(k, r);
      }
           return r;
      
   }
  
/*public DefaultXYZDataset1 getBSeriesCollection(int k){
           
           DefaultXYZDataset1 b =  bdc.get(k);
           if(b==null){
               
                	b = new DefaultXYZDataset1();
                	
                	bdc.set(k, b);
             
           }
           return b;
   }*/
    
   final List<XYSeriesCollection> rdc;
  // final Map<Integer, DefaultXYZDataset1>bdc;  // rdc ordered by distribution, bdc ordered by emissionstatespace order
   
  // XYSeriesCollection current_b, current_r;
   public File chartDistF, chartBF, chartRF;
   public static File makeOrDelete(File outdir, String st){
       File fi = new File(outdir, st);
       if(!fi.exists()) fi.mkdir();
       else{
       File[] f = fi.listFiles();
       for(int i=0; i<f.length; i++){
           f[i].delete();
       }
       }
       return fi;
   }
 //List<Integer> indicesToInclude = new ArrayList<Integer>();
 //String[][] probesToPlot;
 //int[][] snp_alias;
 static int divsize = 1;
 HWECalculator ac;
private String[] names2;
//private ColorAdapter[] ca_rstates;
    public HWEPlot(MarkovModel hmm, HWECalculator ac, List<Integer> loc,List<String> snpid,   String name, File outdir){
        super();
        this.hmm = hmm;
        this.ac = ac;
        this.emStSp =Emiss.getSpaceForNoCopies(ac.ploidy);
       // ac.dc1.getEmStSpace();
       
        List<Integer> index = new ArrayList<Integer>();
       // DataCollection[]ldl = ((MergedDataCollection)ac.dc1).ldl;
        for(int i=0; i<ac.types.size(); i++){
        	//for(int j=i+1; j<ac.types.size(); j++){
            	index.add(i);//new int[] {i,j});
           // }
        }
        this.firstHalf = index.size();
        if(ac.type_len>1){
        	for(int i=0; i<firstHalf; i++){
        		index.add(index.get(i));
        	}
        }
        this.names1 = new String[ index.size()];
      
		this.names2 = new String[ index.size()];
        for(int i=0; i<names1.length; i++){
        	int ind = index.get(i);
        	names1[i] = ac.types.get(ind);//+" vs "+ac.types.get(ind[1]);
        	names2[i] = names1[i];//"Log odds of being "+ac.types.get(ind[0])+" against "+
        //	ac.types.get(ind[1]);
        	if(i<this.firstHalf){
        		names1[i]=("HWE overall: "+names1[i]);
        		names2[i]=names2[i]+" stratified by Copy Number";
        	}
        	else{
        		names1[i]+=("HWE within copy number state "+names1[i]);
        		names2[i]=names1[i];
        		
        	}
        }
        this.setName(name);
       // this.hmm = bwt.hmm;
        this.index=  index;
       // this.bwt = bwt;
        this.location = loc;
        this.ca_b = new ColorAdapter[index.size()];
        this.ca_r = new ColorAdapter[index.size()];
    //   this.ca_rstates = new ColorAdapter[index.size()];
       // int maxi =0;
       ColorAdapter hmmca = ColorAdapter.get((((CompoundMarkovModel)hmm).getMarkovModel(0)));
        for(int i=0; i<index.size(); i++){
        	ca_r[i] = i>=firstHalf ? hmmca : new ColorAdapter();
        	ca_b[i] = new ColorAdapter();
        	 
        }
       // this.data_index_alias = new int[maxi+1];
       // Arrays.fill(data_index_alias, -1);
       // for(int i=0; i<index.size(); i++){
       // 	data_index_alias[index.get(i)] = i;
       // }
        this.noSnps = loc.size()-1;
        //chartDist = new JPanel[index.size()][3];
      //  JComponent[] jpDist = getBottomPane(chartDist);
        chartBF = makeOrDelete(outdir, "plot_"+this.getName()+"_predictionsB");
        chartRF = makeOrDelete(outdir, "plot_"+this.getName()+"_predictionsR");
        chartDistF = makeOrDelete(outdir, "plot_"+this.getName()+"_distributions");
       
        //this.emStSp1 = emStSp.getMembers()[0];
        rdc = new ArrayList<XYSeriesCollection>();
        for(int i=0; i<2*firstHalf; i++){
        	rdc.add(null);
        }
        //bdc = new HashMap<Integer, DefaultXYZDataset1>();
    /*  probesToPlot = Constants.snpsToPlot1(index);
    //   snp_alias = new int[index.size()][];
       this.rb = new Map[index.size()];
       for(int i=0; i<rb.length; i++){
    	   rb[i] =  new HashMap<Integer, XYSeriesCollection>();
       }
       if(probesToPlot.length>0 && probesToPlot[0].length >0 && probesToPlot[0][0].equals("all")){
    	   this.probesToPlot = new String[probesToPlot.length][loc.size()];
    	   snp_alias = new int[probesToPlot.length][loc.size()];
    	   for(int i1=0; i1<probesToPlot.length; i1++){
    	  for(int i=0; i<loc.size(); i++){
            	snp_alias[i1][i] =i;
            	probesToPlot[i1][i] = snpid.get(i);
            	this.rb[i1].put(snp_alias[i1][i], new XYSeriesCollection());
           }
    	   }
       }
       else{
        for(int i=0; i<probesToPlot.length; i++){
        	snp_alias[i] = new int[probesToPlot[i].length];
        	for(int ii=0; ii<probesToPlot[i].length; ii++){
        	snp_alias[i][ii] = snpid.indexOf(probesToPlot[i][ii]);
        	this.rb[i].put(snp_alias[i][ii], new XYSeriesCollection());
        	}
        }
       }*/
       // this.include = new Boolean[]{true};
        jB = new JPanel[this.index.size()];
        jR = new JPanel[this.index.size()];
        if(Constants.plot()>=1){
          
             chartB = new JPanel[this.index.size()];
             chartR= new JPanel[this.index.size()];
            // chartG = new JPanel[probesToPlot.length][];
             chartBNew = new ChartPanel[this.index.size()];
             chartRNew= new ChartPanel[this.index.size()];
             
           /*  chartGNew = new ChartPanel[probesToPlot.length][];
             for(int ii=0; ii<chartGNew.length; ii++){
            	 chartGNew[ii] = new ChartPanel[probesToPlot[ii].length];
            	 chartG[ii] = new JPanel[probesToPlot[ii].length];
             }*/
         
             //jG = new JPanel[this.index.size()];
            jpBS= getJFrame(chartB,  true);
            jpRS= getJFrame(chartR, true);
           //  jG = getJFrame(chartG,  true);
             /*jpBS = //new JSplitPane(JSplitPane.VERTICAL_SPLIT,  
            		 getSplitPane(jB,this.names1,JSplitPane.VERTICAL_SPLIT);
            		// jpDist[1]);
           
             jpRS =// new JSplitPane(JSplitPane.VERTICAL_SPLIT, 
            		 getSplitPane(jR,names1,JSplitPane.VERTICAL_SPLIT);//,jpDist[0]);*/
          //   jB.setMinimumSize(dim2);
           //  jGS = new JScrollPane(jG, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

             //  jR.setMinimumSize(dim2);
             jpBS.setMinimumSize(dim2);
           //  jGS.setMinimumSize(dim2);
          //   jpBS.setDividerLocation(400);
          //   jpRS.setDividerLocation(400);
             jpB = jpBS;//,"b allele", 2*width);
           
             this.addTab("HWE values", new JScrollPane(jpRS, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
           //  this.addTab("CN Category log odds",new JScrollPane( jpBS, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED));
          //   if(probesToPlot.length>0){
          //  	 this.addTab("Scatter", jGS);
          //   }
             this.setSelectedIndex(0);
            // setLayout(new BorderLayout());
          //   add(BorderLayout.CENTER, jpBS);
          //   add(BorderLayout.AFTER_LAST_LINE, jpDist);
        }
        else{
        	for(int i=0; i<jB.length; i++){
        	jB[i] = new JPanel();
        	jR [i]= new JPanel();
        	 jB[i].setLayout(new BorderLayout());
       	  jR[i].setLayout(new BorderLayout());
        	}
        //	jG = new JPanel();
        	 
        	//  jG.setLayout(new BorderLayout());
        }
       // for(int i=0; i<include.length; i++){
        	//if(include(i)) this.indicesToInclude.add(i);
        //	if(include(i))numToInclude++;
      //  }
       
    }
    
   

    public static  JComponent getSplitPane(JComponent[] jg2, String[] names, int split) {
    	if(jg2.length==1) {
    		JPanel jres = new JPanel();
    		jres.setLayout(new BorderLayout());
    		jres.setToolTipText(names[0]);
    		jres.setName(names[0]);
    		JComponent jsp = split==JSplitPane.VERTICAL_SPLIT ? new JScrollPane(
    				jg2[0],JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED) : jg2[0];
        				
        				
    		jres.add(BorderLayout.CENTER,jsp
    				
    		);
    		JPanel jtp = new NamePanel(names[0]);
    		jtp.setSize(jg2[0].size());
    	//	jtp.setText(names[0]);
    		jres.add(BorderLayout.BEFORE_FIRST_LINE,jtp);
    		return jres;
    	 	}
    	else{
    		int mid =(int) Math.ceil(((double)jg2.length)/2.0);
    		JComponent[] left = new JComponent[mid];
    		JComponent[] right = new JComponent[jg2.length-mid];
    		System.arraycopy(jg2, 0, left, 0, mid);
    		System.arraycopy(jg2, mid, right, 0, right.length);
    		String[] leftnames = new String[left.length];
    		String[] rightnames = new String[right.length];
    		System.arraycopy(names, 0, leftnames, 0, mid);
    		System.arraycopy(names, mid, rightnames, 0, right.length);
    		JSplitPane jpBS =    new JSplitPane1(split,  
    				getSplitPane(left, leftnames,split), getSplitPane(right,rightnames,split));
    		
    		jpBS.setOneTouchExpandable(true);
    		  jpBS.setDividerSize(divsize);
    		  jpBS.setDividerLocation((double)left.length / (double)names.length);
    		  return jpBS;
    	}
    }
    
    /*public JComponent[] getBottomPane(JComponent[][] chartDist){
    	  JTabbedPane pane = new JTabbedPane();
    	  JTabbedPane pane1 = new JTabbedPane();
        for(int i=0; i<chartDist.length; i++){
           pane.add( Constants.format()[index.get(i)],getBottomPane(chartDist[i]));
           pane1.add( Constants.format()[index.get(i)],chartDist[i][2]);
        }
    
        return new JComponent[]{pane, pane1};
    }*/
    public JComponent getBottomPane(JComponent[] chartDist){
        for(int i=0; i<chartDist.length; i++){
            chartDist[i] = new JPanel();
           chartDist[i].setMinimumSize(dim);
        }
        
        JComponent tabs =  getTabbedPane(
                new JScrollPane(chartDist[0], JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED),
                new JScrollPane(chartDist[1], JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED)
                          );
     
        tabs.setMinimumSize(dim);
        return tabs;
    }
    public JComponent getTabbedPane(JComponent r_pair, JComponent r_single){
        JTabbedPane pane = new JTabbedPane();
        pane.addTab("r_pair", r_pair);
        pane.addTab("r_single", r_single);
        return pane;
    }
   
    
    
  
    int plot = 0;
    public void setToPlot(int i){
       plot= i;
    }
   public void writeCurrentCharts(int i) {
       try{
           if(!chartBF.exists()) chartBF.mkdir();
           if(!chartBF.exists()) chartRF.mkdir();
           JComponent  currentRPanel = this.getCurrentRPanel(i);
           JComponent  currentBPanel = this.getCurrentBPanel(i);
           if(currentRPanel==null) return ;
          this.writeToZipFile(currentBPanel, chartBF, "odds");
          this.writeToZipFile(currentRPanel, chartRF, "sig");
       }catch(Exception exc){
           exc.printStackTrace();
       }
   }
   
   public void writeToZipFile(Component charts, File dir, String id) throws Exception{
        File out = new File(dir, (id)+".png");
        SVGGraphics2D g = new SVGGraphics2D(out,charts);//ImageConstants.PNG); 
        g.setDeviceIndependent(Constants.plot()==1);
        g.startExport();
        charts.print(g);
        g.endExport();
   }
    
    public void writeToZipFile(Component[] charts, File dir, String[] id) throws Exception{
       for(int i=0; i<charts.length; i++){
           File out = new File(dir, (id[i])+".png");
           charts[i].setSize(500,500);
           ImageGraphics2D g  = new ImageGraphics2D(out,charts[i], ImageConstants.PNG); 
           g.setDeviceIndependent(Constants.plot()==1);
           //GIFGraphics2D g = new GIFGraphics2D(out,charts[i]); 
           g.startExport();
           charts[i].print(g);
           g.endExport();
       }
    }
   
  
  /* static JFrame getJFrame(JComponent jpBS, String name, int width){
       JFrame jf = new JFrame(name);
       jf.getContentPane().add(jpBS);
       jf.setSize(width,800);
       jf.setDefaultCloseOperation(jf.EXIT_ON_CLOSE);
       return jf;
   }
    public boolean include(int k){
    	
        if(include!=null){
            if(include[k]==null){
            	include[k] = true;
            
            }
            return include[k];
        }
        else{
        	
        	return true;
        }
    	
        
    }*/
   // Boolean[] include;
    JComponent[] jp;
    JComponent getJFrame( JComponent[][] cp, boolean all){
    	 jp= new JComponent[cp.length];
    	for(int i=0; i<cp.length; i++){
    		jp[i] =getJFrame(cp[i], all);
    	}
    	return this.getSplitPane(jp,names1, JSplitPane.HORIZONTAL_SPLIT);
    }
    
    
    JComponent getJFrame( JComponent[] cp, boolean all){
        JPanel jpBS = new JPanel();
    
     	
     	  jpBS.setLayout(new BoxLayout(jpBS, BoxLayout.Y_AXIS));
       
        
    //    int height =  dim.height;
         for(int i=0; i<cp.length; i++){
             if(all ){
             cp[i] = new JPanel();
             jpBS.add(cp[i]);
           // cp[i].setMinimumSize(new Dimension(dim.width,height*cp.length ));
                // jpBS.add(cp[i]);
             }
         }
    //    this.pack();
         
         return jpBS;
     }
    
  
    public double[] noSamples;
    
    
   /* public ChartPanel[] getDistCharts(int index){
    	if(noSamples==null){
    		noSamples = new double[bwt.probDists[index].mvf.single.size()];
    	}
        ChartPanel[] cp = new ChartPanel[3];
        bwt.probDists[index].mvf.updateSampleSize(noSamples);
          cp[0] = 
        	  Constants.format()[index].startsWith("geno") ? null : 
        	  plotR(bwt.probDists[index].s1, bwt.probDists[index].mvf, null,index,  0, this.lu_rpair, "Conditional distributions over LRR", "LRR", bwt.probDists[index].ca_r);
        
          cp[1] = plotR(bwt.probDists[index].mvf.single, null, 
        		 noSamples,index, 1, lu_r, "Conditional distributions over LRR", "LRR",  bwt.probDists[index].ca_r);
          cp[2] = plotR(bwt.probDists[index].s2, null, null,index,2, lu_b, "Conditional distributions over BAF","BAF", bwt.probDists[index].ca_b);
          
          return cp;
    }*/
    
  
    
    public void reinitialise() {
      //  this.currentIndividual = l;
    	
    		this.ac.initialise();
    		/*for(int i=0; i<rb.length; i++){
    		for(Iterator<XYSeriesCollection> it = this.rb[i].values().iterator(); it.hasNext();){
    			XYSeriesCollection xys = it.next();
    			 for(int k=0; k<xys.getSeriesCount(); k++){
    		            xys.getSeries(k).clear();
    		        }
    		}
    		}*/
    	
      /* DefaultXYZDataset1 current_b = this.getBSeriesCollection(l);
       XYSeriesCollection current_r  = this.getRSeriesCollection(l);
       //for(int j=0; j<current_b.length; j++){
        for(int k=0; k<current_r.getSeriesCount(); k++){
            current_r.getSeries(k).clear();
        }
       //}
      // for(int j=0; j<current_b.length; j++){
       // for(int k=0; k<current_b.getSeriesCount(); k++){
        //    current_b.getS
       // }
      // }*/
        
    }
    
  /*  private void addRBSeries(String compoundName) {
    	for(int i=0; i<this.rb.length; i++){
    	for(Iterator<XYSeriesCollection> it = rb[i].values().iterator(); it.hasNext();){
    		XYSeriesCollection nxt = it.next();
    		nxt.addSeries(new XYSeries(compoundName));
    		// TODO Auto-generated method stub
    	}
    	
        }
    }*/
    final MarkovModel hmm;
    //final int[] data_index_alias; //converts from index to position in index
 //  final Map<Integer, XYSeriesCollection>[] rb;; 
 // short[] ind_ = new short[1];
   /** i is index of individual
     * stateDistribution has the distribution over genotypes
     *  */
   
    public synchronized void addedInformation(StateDistribution dist, int ll, int i,String name){
    	HaplotypeEmissionState state = (HaplotypeEmissionState)this.ac.dc1.dataL.get(name);
    	   double[] prob =  PairEmissionState.pool.getObj(Emiss.getSpaceForNoCopies(state.noCop()).genoListSize());
    	   double sumNonMix =   Sampler.getProbOverStates(dist, this.hmm, state, i,prob);
    	 double sum =   Sampler.getProbOverStates(dist, this.hmm, state, i,prob);
    	 if(sumNonMix/sum> Constants.imputedThresh(1)){
    		 this.ac.scoreChi1(prob, state.dataIndex(),i);
    	 }
       
        PairEmissionState.pool.returnObj(prob);
    }
    
   
    private int get(List<Integer> loc, int i) {
    	for(int j=0; j<loc.size(); j++){
    		if(loc.get(j)>=i) return j;
    	}
    	return loc.size();
    }
    public void updateData(){
    	 double lowP = 1.0;
    	 int pos=0;
    	 int lowInd =0;
    	//Map<Integer> rs_pos = new HashMap<Integer, Double>();
    	for(int i=0; i<location.size(); i++){
    	 double x = location.get(i);
    	for(int ll=0; ll<firstHalf; ll++){
    		 int ind = this.index.get(ll);
    		   Double[]sign_ =  this.ac.getSignificance(i,  ind);
    		 if(ac.type_len>1){ //plots state significance 
    			// DefaultXYZDataset1 current_b = this.getBSeriesCollection(ll+firstHalf);
       		     XYSeriesCollection current_r  = this.getRSeriesCollection(ll+firstHalf);
       		     int jj=0;
       		  
    		   for(int j=1; j<ac.type_len; j++){
    	        double sign =sign_[j];
    	        double p = Math.max( 1e-200 ,sign);
     	       
     	     
    	      
    	        //int[] b_alias =  ((CompoundMarkovModel)bwt.hmm).probB(ind).b_alias;
    	       current_r.getSeries(j-1).addOrUpdate(x,  Math.min(Constants.maxLogP(),-Math.log10(p)));
    	      /*      Double[][] odds = this.ac.oddsRatio(i,ind[0], ind[1], j);
    	       for(int i1=0; i1<odds[0].length;i1++){
    	    	   //if(odds[i1]!=null){
    	    		   
    	    	   current_b.updateSeries(jj, i,this.location.get(i),odds[0][i1],Math.max(0, Math.log(odds[1][i1]))*1000);//.addOrUpdate(x, odds[i1]);
    	    	   jj++;
    	    	   }
    	    	  */
    	       }
    		 }
    		   double sign = sign_[0];
   	        if(sign < lowP){
   	        	lowP = sign;
   	        	pos = i;
   	        	lowInd =ll;
   	        	
   	        
   	        }
   	        XYSeriesCollection current_r  = this.getRSeriesCollection(ll);
   	     double p = Math.max( 1e-200 ,sign);
   	        //int[] b_alias =  ((CompoundMarkovModel)bwt.hmm).probB(ind).b_alias;
   	       current_r.getSeries(0).addOrUpdate(x, Math.min(Constants.maxLogP(),-Math.log10(p)));
    	       /* Double[][] odds = this.ac.oddsRatio(i,ind[0], ind[1], 0);
    	        DefaultXYZDataset1 current_b = this.getBSeriesCollection(ll);
    	       for(int i1=0; i1<odds[0].length;i1++){
    	    	   current_b.updateSeries(i1, i,this.location.get(i),odds[0][i1],Math.max(0, Math.log(odds[1][i1]))*1000);//.addOrUpdate(x, odds[i1]);
    	    	 
    	       }
    	       // }*/
    	}
    	}
    	if(false){
    	System.err.println("lowest pvalue is "+lowP+" "+this.location.get(pos));
    	int max = ac.dc1.snpid.size();
    	int[] snp_a = new int[Constants.expandScatter*2+1];
    	Arrays.fill(snp_a, -1);
    	String[] rs = new String[Constants.expandScatter*2+1];
    	double[] pvalue = new double[Constants.expandScatter*2+1];
    	snp_a[Constants.expandScatter] = pos;
    	rs[Constants.expandScatter] = this.ac.dc1.snpid.get(pos);
    	
    	 XYSeries current_r  = this.getRSeriesCollection(lowInd).getSeries(0);
    	 pvalue[Constants.expandScatter] = current_r.getY(current_r.indexOf(this.location.get(pos))).doubleValue();
    //		 -Math.log10(lowP);
    	for(int i=1; i<Constants.expandScatter+1; i++){
    		int l =Constants.expandScatter-i; 
    		snp_a[l] = pos -i;
    		rs[l] = snp_a[l]<0  ? "" : ac.dc1.snpid.get(snp_a[l]);
    		
    		pvalue[l] =snp_a[l]<0  ? 1.0 : current_r.getY(current_r.indexOf(this.location.get(snp_a[l]))).doubleValue();
    		
    		int r=  Constants.expandScatter+i;
    		snp_a[r] = pos +i;
    	
    		rs[r] = snp_a[r]>=max ? "" :ac.dc1.snpid.get(snp_a[r]);
    		pvalue[r] = snp_a[r]>=max ? 1.0 : current_r.getY(current_r.indexOf(this.location.get(snp_a[r]))).doubleValue();
    		
    		if(snp_a[r]>=max) snp_a[r] = -1;
    	}
    	//if(IndividualPlot.updateProbes)
    	this.firePropertyChange("pval", null,new Object[] {snp_a, rs, pvalue});
    	}
    }
    private void sum(double[] emiss, double[] prior_st) {
		Arrays.fill(prior_st, 0.0);
		for(int i=0; i<emiss.length; i++){
			prior_st[emStSp.getCN(i)]+=emiss[i];
		}
		
	}

	double[] res;
	Color zeroCol = new Color(0,0,0);
   /* private synchronized double[] getProbOverStates(StateDistribution emissionC,
            MarkovModel hmm, HaplotypeEmissionState obj, int i) {
    	
    	 EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
         if(res==null)
       	  res = new double[emstsp.defaultList.size()];
      Arrays.fill(res,  0.0);
      double sum = 0;
      for(int j=1; j<emissionC.dist.length; j++){
          double p = emissionC.dist[j];
          EmissionState state_j = (EmissionState) hmm.getState(j);
        
          if(p>0){
        	  sum+=obj.calcDistribution((EmissionState) state_j, i,  p);
        	  EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
        	  double[] dist = state_j.distribution;
              for(int k=0; k<dist.length; k++){
            	  res[emstsp.get(emstsp1.get(k))]+=dist[k];
              }
          }
      }
      for(int k=0; k<res.length; k++){
          res[k] = res[k] / sum;
      }
        return res;
    }*/
    
    public JComponent getCurrentBPanel(int l){
        if(this.chartB!=null) return chartB[l];
        else return new JPanel();
    }
  
  /*  public JComponent getCurrentGPanel(int l, int k){
        if(this.chartG!=null) return chartG[l][k];
        else return new JPanel();
    }*/
  
    public JComponent getCurrentRPanel(int l){
        if(this.chartR!=null) return chartR[l];
        else return new JPanel();
    }
    
    public static  Font font= new Font("Tahoma", Font.PLAIN, 10);
    
    public void update(){
        for(int i=0; i<chartR.length; i++){
            if(chartRNew[i]==null) continue;
            JComponent currentRPanel = this.getCurrentRPanel(i);
            JComponent currentBPanel = this.getCurrentBPanel(i);
         
            if(currentRPanel==null) continue;
            currentBPanel.setMinimumSize(dim);
            currentRPanel.setMinimumSize(dim);
           currentBPanel.setSize(dim);
            currentRPanel.setSize(dim);
              currentRPanel.removeAll(); 
             currentBPanel.removeAll();
             
             currentRPanel.add(chartRNew[i]);
           //  currentBPanel.add(chartBNew[i]);
        }
        /*for(int i=0; i<chartG.length; i++){
        	for(int ik=0; ik<chartG[i].length; ik++){
	        	if(chartGNew[i][ik]==null) continue;
	        	JComponent currentGPanel = this.getCurrentGPanel(i,ik);
	        	if(currentGPanel==null) continue;
	        	currentGPanel.setMinimumSize(dimG);
	        	currentGPanel.setSize(dimG);
	        	currentGPanel.removeAll();
	        	currentGPanel.add(chartGNew[i][ik]);
        	}
        }*/
    }
    public static Stroke dashed =
        new BasicStroke(0.2f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, new float[] {10.0f, 1.0f, 10.0f}, 1.0f);
    public static Stroke dotted =
        new BasicStroke(0.8f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f, new float[] {2.0f, 20.0f}, 0.0f);
 
   /* public void check(XYSeriesCollection datas){
        Map<Double, Double> s = new HashMap<Double, Double>();
        for(int i=0; i<datas.getSeriesCount(); i++){
            XYSeries series_i = datas.getSeries(i);
            for(int k=0; k<series_i.getItemCount(); k++){
                double x = series_i.getX(index).doubleValue();
                double y = series_i.getY(index).doubleValue();
                Double y_1 = s.get(x);
                if(y_1!=null) throw new RuntimeException("already in !"+x+" "+y+" "+y_1);
               s.put(x, y);
               
            }
        }
    }*/
    
    public Shape getShape(char ch){
      return new Rectangle(2,2);
       // Graphics2D gr = (Graphics2D)this.getGraphics();
       // return this.font10.createGlyphVector(gr.getFontRenderContext(), new char[] {ch}).getOutline();
    }
   static char[] shapes = new char[] {'x', 'o', '^', '|','m','#','+','=',  '%','*','a','b','c','d','e','f'};
    public static  Font font8 = new Font("SansSerif", Font.PLAIN, 10);
    public static  Font font6 = new Font("SansSerif", Font.PLAIN, 8);
static Color transPc = Color.getHSBColor(0f, 1f, 1f);
Color c = new Color(255, 255, 255, 0);// - this is for outline only!
    public  JFreeChart graphB(XYZDataset datas_, String title, Color[] ca , ColorAdapter ca1) {
    	XYZDataset[] datas = new XYZDataset[]{datas_};
  
        final JFreeChart chart = ChartFactory.createBubbleChart(//title, xAxisLabel, yAxisLabel, dataset, orientation, legend, tooltips, urls)
                title+"_bubble",
                "Position on chromosome "+Constants.chrom0(), // domain axis label
                "Value ", // range axis label
                datas[0], // data
                PlotOrientation.VERTICAL, true, // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
      
        XYPlot plot = (XYPlot) chart.getPlot();
        for(int i=0; i<datas.length; i++){
        	  plot.setDataset(i,datas[i]);
        }
      
      
       
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderPaint(Color.white);
        plot.setBackgroundPaint(Color.white);
     //   whiteIndex.clear();
            NumberAxis yAxis = (NumberAxis)plot.getRangeAxis();
            NumberAxis xAxis = (NumberAxis)plot.getDomainAxis();
          //  yAxis.setRange(new Range(-2.0,2.0))
            chart.getTitle().setFont(font8);
            yAxis.setTickLabelFont(font6);
            yAxis.setLabelFont(font6);
            xAxis.setTickLabelFont(font6);
            xAxis.setLabelFont(font6);
            if(location.size()>0){
	            xAxis.setAutoRange(false);
	            xAxis.setLowerBound(location.get(0)-1000);
	            xAxis.setUpperBound(location.get(location.size()-1)+1000);
            }
         
            
            
         
            plot.setDomainGridlinePaint(Color.WHITE);
           plot.setRangeGridlinePaint(Color.WHITE);
          
   		//XYBubbleRenderer rend1 = new XYBubbleRenderer(
   			//	XYBubbleRenderer.SCALE_ON_DOMAIN_AXIS);
   		//chart.getXYPlot().setRenderer(0, rend);
   		//chart.getXYPlot().setRenderer(1, rend1);
            for(int kk=0; kk<datas.length; kk++){
           // int ik=0;
            	 XYBubbleRenderer rend_ = new XYBubbleRenderer(
            				XYBubbleRenderer.SCALE_ON_DOMAIN_AXIS);
            	 
            	
         	//	Color c1 = new Color(255, 255, 255, 0);// - this is for outline only!
         		rend_.setAutoPopulateSeriesFillPaint(false);
         	
         	//	rend_.setSeriesFillPaint(0, c);
         		
         	//	rend_.setPaint(c);
         		

         	//	rend_.setBaseFillPaint(c);
         		
         		//rend_.setBaseOutlinePaint(Color.GREEN);
         		
         	
            //	XYBubbleRenderer rend_ = kk==0 ? rend:rend1;
            plot.setRenderer(kk, rend_);
         

            for(int i=0; i<datas[kk].getSeriesCount(); i++){
            	String[] series_name = datas[kk].getSeriesKey(i).toString().split("_");
            
                Color colors = ca1==null ? ca[i] :ca1.getColor(series_name[1]) ;//.getColor(datas[kk].getSeriesKey(i).toString());
             //   if(whiteIndex.contains(i))     colors = Color.WHITE;
              //  rend_.setSeries
               
               // rend_.setSeriesFillPaint(i, zeroCol);
            //   rend_.setSeriesFillPaint(i,transPc);
        		//rend_.setSeriesPaint(0, null);// setSeriesFillPaint(series, paint)
             	
                rend_.setSeriesOutlinePaint(i, colors);
               if(ca==null) rend_.setSeriesPaint(i, HMMPanel.modify(colors, Double.parseDouble(series_name[0])/2.0));
               else rend_.setSeriesPaint(i, colors);
                		//new Color(colors.getRed(), colors.getGreen(), colors.getBlue(),Integer.parseInt(series_name[0])*254));
                //if(ca==null) rend_.setSeriesShape(i, IndividualPlot.shapes[Integer.parseInt(series_name[0])]);
             //   renderer[kk].setSeriesShape(i, shape[kk]);
              
              //  ik++;
            }
            }
          /*  renderer1.setShapesVisible(false);
            renderer2.setShapesVisible(false);
            renderer2.setLinesVisible(true);
            renderer2.setLinesVisible(true);
            //renderer1.setStroke(new BasicStroke(BasicStroke.JOIN_MITER));
           // renderer1.setStroke(dashed);
            renderer2.setStroke(dotted);*/
            
           // plot.setRenderer(0,renderer3);
           // plot.setRenderer(1, renderer);
          /*  for(int i=0; i<8; i++){
            	plot.setRenderer(i, renderer);
            }
            if(lu[0]!=null && lu[1]!=null && false){
                plot.setRenderer(1, renderer1);
                plot.setRenderer(2, renderer2);
            }*/
            chart.getXYPlot().mapDatasetToDomainAxis(1, 0);
    		chart.getXYPlot().mapDatasetToRangeAxis(1, 0);
    		chart.getXYPlot().mapDatasetToDomainAxis(0, 0);
    		chart.getXYPlot().mapDatasetToRangeAxis(0, 0);

              chart.setBackgroundPaint(Color.white);
          //    plot.setBackgroundPaint(Color.black);
             return chart;
}
    
  /*  public  JFreeChart graph(XYSeriesCollection datas,  String title, ColorAdapter ca ) {
        
        final JFreeChart chart = ChartFactory.createXYLineChart(
                title+" chromosome "+Constants.chrom0(),
                "LRR ", // domain axis label
                "BAF ", // range axis label
                datas, // data
                PlotOrientation.VERTICAL, false, // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
      
        XYPlot plot = (XYPlot) chart.getPlot();
        
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderPaint(Color.white);
        plot.setBackgroundPaint(Color.white);
     //   whiteIndex.clear();
            NumberAxis yAxis = (NumberAxis)plot.getRangeAxis();
            NumberAxis xAxis = (NumberAxis)plot.getDomainAxis();
            chart.getTitle().setFont(font8);
            yAxis.setTickLabelFont(font6);
            yAxis.setLabelFont(font6);
            xAxis.setTickLabelFont(font6);
            xAxis.setLabelFont(font6);
           xAxis.setAutoRange(false);
           xAxis.setRange(new Range(Constants.minR(), Constants.maxR()));
        //  if(Constants.ascn()) yAxis.setRange(new Range(Constants.minR(), Constants.maxR()));
            plot.setDomainGridlinePaint(Color.WHITE);
           plot.setRangeGridlinePaint(Color.WHITE);
            final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
            renderer.setShapesFilled(Boolean.TRUE);
            renderer.setLinesVisible(false);
            final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
            final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
            
            
           if(dotted!=null) renderer1.setBaseStroke(dotted);
            
            int ik=0;
            for(int i=0; i<datas.getSeriesCount(); i++){
                Color colors = ca.getColor(datas.getSeriesKey(ik).toString());
             
                renderer.setSeriesPaint(i, colors);
                renderer.setSeriesShape(i, this.getShape(shapes[0]));
                renderer1.setSeriesPaint(i, colors);
                renderer2.setSeriesPaint(i, colors);
                ik++;
            }
          
            renderer1.setShapesVisible(false);
            renderer2.setShapesVisible(false);
            renderer2.setLinesVisible(true);
            renderer2.setLinesVisible(true);
            //renderer1.setStroke(new BasicStroke(BasicStroke.JOIN_MITER));
            renderer1.setStroke(dashed);
          if(dotted!=null)  renderer1.setStroke(dotted);
            plot.setRenderer(0,renderer);
          
  
              chart.setBackgroundPaint(Color.white);
             return chart;
}*/
    static  JFreeChart getChart(XYSeriesCollection exp, 
             XYSeriesCollection obs, String name, String yaxis, ColorAdapter ca,
             Map<Double[], String[]> sn_params
            ) {
     
         final JFreeChart chart = ChartFactory.createXYLineChart(
                 name,
                 yaxis, // domain axis label
                 "Frequency expected", // range axis label
                 obs, // data
                 PlotOrientation.HORIZONTAL, false, // include legend
                 true, // tooltips?
                 false // URL generator? Not required...
                 );
         
         ((NumberAxis) ((XYPlot) chart.getPlot()).getRangeAxis())
                 .setAutoRangeIncludesZero(false);
         chart.setBackgroundPaint(Color.white);
//         chart.getLegend().setAnchor(org.jfree.chart.Legend.SOUTH);
         final XYPlot plot = chart.getXYPlot();
         // plot.setRenderer(new XYLineAndShapeRenderer());
         // ( (XYLineAndShapeRenderer)plot.getRenderer()).set.s
         plot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
         // plot.getDomainAxis().setFixedAutoRange()
         plot.setDataset(1, exp);
         plot.mapDatasetToRangeAxis(1, 1);
       //  domainAxis.setRange(min[0], min[1]);
         final ValueAxis axis1 =Constants.logplot() ?  new LogarithmicAxis("Conditional probability distribution") : new NumberAxis("Conditional probability distribution");
         final ValueAxis axis2 =Constants.logplot() ? new LogarithmicAxis("Cumulative observed count") : new NumberAxis("Cumulative observed count");
         axis1.setTickLabelsVisible(true);
         axis2.setTickLabelsVisible(true);
         axis1.setTickLabelFont(font4);
         axis2.setTickLabelFont(font4);
         if(Constants.logplot()){
        	 axis1.setAutoTickUnitSelection(false);
        	 axis2.setAutoTickUnitSelection(false);
         }
       
        // axis1.setStandardTickUnits(new NumberTickUnit());
         plot.setRangeAxis(0,axis1);
         plot.setRangeAxis(1, axis2);
         final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
         final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
         
         int ik=0;
         for(int i=0; i<obs.getSeriesCount(); i++){
            // if(ik==color.length) ik=0;
          Color colors = ca.getColor(obs.getSeries(i).getKey().toString());
          renderer1.setSeriesPaint(i, colors);
          renderer2.setSeriesPaint(i, colors);
         // renderer.setSeriesShape(i font16.createGlyphVector(arg0, arg1));
          ik++;
         }
         // renderer2.setToolTipGenerator(new
         // StandardCategoryToolTipGenerator());
         if(dotted!=null)renderer2.setBaseStroke(dotted);
         renderer1.setShapesVisible(false);
         renderer2.setShapesVisible(false);
         plot.setRenderer(0, renderer1);
         plot.setRenderer(1, renderer2);
         plot.setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
      //  LegendTitle legend =(LegendTitle) chart.getLegend();
        
     // legend.setItemFont(font8);
         if(sn_params!=null){
         for(Iterator<Map.Entry<Double[], String[]>> para = sn_params.entrySet().iterator(); para.hasNext();){
        	 Map.Entry<Double[], String[]> ent = para.next();
        	 String[] val = ent.getValue();
        	 Double[] x = ent.getKey();
        	 if(x!=null && Constants.logplot() && Constants.addAnnotationToDistGraphs()){
	     	XYTextAnnotation annot = new XYTextAnnotation(val[0]
					+ "", x[0], x[1]);
			annot.setFont(font10);
			annot.setPaint(ca.getColor(val[1]));
			plot.addAnnotation(annot);
        	 }
         }
         }
         return chart;
     }
    
    
   // XYSeriesCollection[] lu_rpair =new XYSeriesCollection[2];
    XYSeriesCollection[] lu_r= new XYSeriesCollection[2];
    XYSeriesCollection[] lu_b = new XYSeriesCollection[2];
    public ChartPanel plotR(List<ProbabilityDistribution> ss1,
            ProbMultivariate mvf, double[] sum, int data_index, int index1, XYSeriesCollection[] lowerupper,
            String name, String namey, ColorAdapter ca) {
        List<Double[]> cis = new ArrayList<Double[]>();
        int len = ss1==null ? mvf.pairs.size() : ss1.size();
        for(int i=0; i<len; i++){
            cis.add(new Double[] {null, null});
        }
        double[] range = new double[] {0.05, 0.95};
        List<String> names = new ArrayList<String>();
        if(mvf!=null){
            for(int i=0; i<mvf.pairs.size(); i++){
                names.add(mvf.getCompoundName(i));
            }
        }
        else{
            for(int i=0; i<ss1.size(); i++){
                names.add(ss1.get(i).name());
            }
        }
        Map<Double[], String[]> det = new HashMap<Double[], String[]>();
        XYSeriesCollection[] obs = mvf==null ?  plot(ss1, sum, cis, range, names, det) : plot(mvf.pairs, sum, cis, range, names, det);
         XYSeriesCollection theor = obs[1];
       
          //  Range xrange = theor.getDomainBounds(false);
           // double y_max = getYMax(theor);
            lowerupper[0] = new XYSeriesCollection();
            lowerupper[1] = new XYSeriesCollection();
          //  double  min = xrange.getLowerBound();
        //   double max = xrange.getUpperBound();
          
            for(int i=0; i<theor.getSeriesCount(); i++){
                XYSeries lowers = new XYSeries(theor.getSeries(i).getKey());
                XYSeries uppers = new XYSeries(theor.getSeries(i).getKey());
               
                Double[] ci =cis.get(i);
                lowers.add(this.location.get(0), ci[0]); lowers.add(this.location.get(noSnps), ci[0]);
                uppers.add(this.location.get(0), ci[1]); uppers.add(this.location.get(noSnps),  ci[1]);
                lowerupper[0].addSeries(lowers);
                lowerupper[1].addSeries(uppers);
            }
       
        JFreeChart chart = getChart(obs[0], obs[1],name,namey, ca, det);
     
          //  chartDist[data_index][index1].removeAll();
            final ChartPanel cp  = new ChartPanel(chart,
            		 1200, //width
            		 400, //height
            		 1200, //mindrawWidth
                     200, //mindrawHeight
                     1200, //maxDrawWith
                     400,//maxDrawHeight
                     ChartPanel.DEFAULT_BUFFER_USED,
                     true,  // properties
                     true,  // save
                     true,  // print
                     true,  // zoom
                     true   // tooltips		
            
            );

          //  chartDist[index].setMinimumSize(cp.getMinimumSize());
          //  chartDist[data_index][index1].add(cp);
            cp.setMinimumSize(dim);
            cp.setSize(dim);
          //  chartDist[data_index][index1].setMinimumSize(dim);
         //   chartDist[data_index][index1].setSize(dim);
         return cp;
//            Rectangle r = chartDist[index].getBounds();
          //  chartDist[index].repaint(1000);//(jpDist.getGraphics());
    }




    private double getYMax(XYSeriesCollection theor) {
        double max = Double.NEGATIVE_INFINITY;
        for(int i=0; i<theor.getSeriesCount(); i++){
            XYSeries series = theor.getSeries(i);
            for(int j=0; j<series.getItemCount(); j++){
                double y = series.getY(j).doubleValue();
                if(y>max){
                    max = y;
                }
            }
        }
        return max;
    }
    private double[] getCI(XYSeries ser, double pl, double pm) {
        double max = ser.getY(ser.getItemCount()-1).doubleValue();
        double[] res = new double[2];
        for(int i=0; i<ser.getItemCount(); i++){
            if(ser.getY(i+1).doubleValue() / max> pl ){
                res[0] = i;
                break;
            }
        }
        for(int i=ser.getItemCount()-1; i>=0;i--){
            System.err.println(ser.getY(i-1).doubleValue()/max);
            if(ser.getY(i-1).doubleValue() / max< pm ){
                res[1] = i;
                break;
            }
        }
        return res;
    }


    private static double sum(Collection<ProbabilityDistribution> s1) {
        double sum=0;
        for(Iterator<ProbabilityDistribution> it = s1.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
                sum+=pdist.sum();
        }
        return sum;
    }




    public static final List<ProbabilityDistribution> extract(Map<String, List<ProbabilityDistribution>> m){
        List<ProbabilityDistribution> res = new ArrayList<ProbabilityDistribution>();
        for(Iterator<List<ProbabilityDistribution>> it = m.values().iterator(); it.hasNext();){
            res.add(it.next().get(0));
        }
        return res;
    }
    public static final Map<String, List<ProbabilityDistribution>> transform(Collection<ProbabilityDistribution> s2){
        Map<String, List<ProbabilityDistribution>> m = new HashMap<String, List<ProbabilityDistribution>>();
        
        for(Iterator<ProbabilityDistribution> it = s2.iterator(); it.hasNext();){
            ProbabilityDistribution pdist = it.next();
            String st = pdist.toString();
            List<ProbabilityDistribution > l = m.get(st);
            if(l==null){
                m.put(pdist.toString(), l = new ArrayList<ProbabilityDistribution>());
            }
            l.add(pdist);
        }
        return m;
    }
    
   static boolean rescale = false;
public static double plotThresh = 1e-3;
    public static XYSeriesCollection[] plot(List<ProbabilityDistribution> coll, double[] sums, List<Double[]> cis,double[] mm, List<String> names,
         Map<Double[], String[]> det){
        XYSeriesCollection datas1 = new XYSeriesCollection();
        XYSeriesCollection datas2 = new XYSeriesCollection();
    //    double tot =0;
        for(int ij=0; ij<coll.size(); ij++){
        	ProbabilityDistribution pdist1 =  coll.get(ij);
        	SkewNormal pdist = pdist1 instanceof SkewNormal? (SkewNormal) pdist1 :  (SkewNormal) ((Mixture)pdist1).dist[0];
            if(pdist !=null ){
               double sum = 
                   sums==null ? 
                  pdist.sum() : sums[ij];
         //      tot+=sum;
               //   if(sum==0) sum=1;
                  XYSeries obs = new XYSeries(names.get(ij));
                  String key = names.get(ij).toString();
                 boolean swtch =   key.startsWith("1.0");
                  double cum = pdist.plotObservations(names.get(ij), true, obs, swtch);
                  
                   XYSeries theor = new XYSeries(names.get(ij));
                   pdist.plotTheoretical(names.get(ij), false, theor);
                  
                    double max = 0;
                    double theoSum=0;
                  for(int ik=0; ik<theor.getItemCount(); ik++){
                      theoSum+=theor.getY(ik).doubleValue();
                      if(theor.getY(ik).doubleValue()>max) max = theor.getY(ik).doubleValue();
                  }
                  double cumSum =0;
                  Double[] ci = cis.get(ij);
                  if(Math.abs(pdist.shape())<0.001){
                      pdist.ci(ci, mm);
                  }
                  else{
                  for(int ik=0; ik<theor.getItemCount(); ik++){
                      cumSum+=theor.getY(ik).doubleValue();
                      if(ci[0]==null && cumSum/theoSum>mm[0]){
                          ci[0] = theor.getX(ik).doubleValue();
                      }
                      if(ci[1]==null && cumSum/theoSum>=mm[1]){
                          ci[1] = theor.getX(ik).doubleValue();
                          break;
                      }
                  }
                  }
                  Double[] mode = new Double[] {pdist.location,0.5*sum};//( !Constants.logplot() ? 1.0 : sum)};//0.1*(sum/max) };
                  det.put(mode, new String[] {pdist.toString(), names.get(ij)});
                
               //   int obsCount =0;
                  System.err.println("key "+key);
                  /*if(swtch){
                	  for(int ik=obs.getItemCount()-1; ik>=0; ik--){
                		  double v = obs.getY(ik).doubleValue();
                		  double newV =cum - v; 
                		  if(newV > plotThresh){
                		  obs.update(ik, new Double( newV));
                		//  obsCount++;
                		  }
                		  else{
                			  obs.remove(ik);
                		  }
                	  }
                  }*/
                  for(int ik=theor.getItemCount()-1; ik>=0; ik--){
                	  double frac = (theor.getY(ik).doubleValue()/max);
                		/*  if(swtch){
                			  frac = 1.0 - frac;
                          }*/
                	  if(cum<plotThresh && ! Constants.logplot())   {
                      	theor.remove(ik);
                      }
                	  else if(rescale){
                        double the =  frac* sum;//( !Constants.logplot() ? 1.0 : sum);
                        if(the>plotThresh ||  Constants.logplot()){
                          theor.update(ik, new Double(the));
                        }
                	  }
                    
                    }
                /*  if(normalise){
                  for(int ik=0; ik<obs.getItemCount(); ik++){
                     obs.update(ik, new Double((obs.getY(ik).doubleValue()/sum)));
                  }
                  }*/
                 // if(theor.getItemCount()>0)
                    datas2.addSeries(theor);
                  //if(obs.getItemCount()>0)  
                  datas1.addSeries(obs);
            }
        }
        
        return new XYSeriesCollection [] {datas1, datas2};
//      
        /*if(true){
            try{
                Thread.currentThread().wait(10000);
            }catch(Exception exc){
                exc.printStackTrace();
            }
        }*/
    }
    
    long lastupdate = -1;
    
 public void changeAxis(Range range, boolean dox,  int k1){
		
		/*for(int kk=0; kk<chartRNew.length; kk++){
			if(chartRNew[kk]==null) continue;
			//if(kk!=k1){
				ValueAxis axisr, axisb;
				if(dox){
				
					axisr = chartRNew[kk].getChart().getXYPlot().getDomainAxis();
					axisb = chartBNew[kk].getChart().getXYPlot().getDomainAxis();
					
				}
				else{
					axisr = chartRNew[kk].getChart().getXYPlot().getRangeAxis();
					axisb = chartBNew[kk].getChart().getXYPlot().getRangeAxis();
				}
				if(!axisr.getRange().equals(range)){
					axisr.setRange(range, true, false);
				}
				if(!axisb.getRange().equals(range)){
					axisb.setRange(range, true, false);
				}
				chartRNew[kk].updateUI();
				chartBNew[kk].updateUI();
			//}
		}*/
    }
 final ColorAdapter[] ca_r, ca_b;
 final Color[] ca_b1 = EmissionStateSpace.col;
  //int cnt1 =0;int cnt=0;
    public synchronized void propertyChange(PropertyChangeEvent arg0) {
    	try{
        String nme = arg0.getPropertyName();
     //   Logger.global.info("prop change "+nme);
        if(nme.equals("setToPlot")){
            int level = (Integer) arg0.getNewValue();
            setToPlot(level);
            return;
        }
        if(Constants.plot()<=1 && plot==0) return;
     
        if(nme.equals("pre_exp")){
            
           
                reinitialise();
              
           
         
        }
         
        else if(nme.equals("pre_exp")){
          /*if(plot==2 && Constants.printPlots()){
        	setup(this.cnt);
          }*/
        }
        else if(nme.equals("emiss")){
            if(plot<2 && Constants.plot()<2) return;
            Object[] obj = (Object[]) arg0.getNewValue();
            Integer l = (Integer)obj[1];
            StateDistribution dist = (StateDistribution) obj[0];
            Integer i = (Integer)obj[2];
            Integer index = (Integer)obj[3];
            String name = (String) obj[4];
            Integer ploi = (Integer) obj[5];
            if(ploi.intValue()!=this.ac.ploidy) return;
           addedInformation( dist, l, i, name);
        }
        else if(nme.equals("finished")){
          /*  if(plot<2 && Constants.plot()<2) return;
      
            Object[] obj = (Object[]) arg0.getNewValue();
            
            final  Integer l = (Integer)obj[1];*/
            
          
            
             //if(!include(l)) return ;
            // int ind1 = indicesToInclude.indexOf(l);
          
        }
        else if(nme.equals("expec_i")){
        
            
        }
        else if(nme.equals("expectation1")){
        	   this.updateData();
        	   Dimension dim_4 = new Dimension(dim.width, dim.height*this.rdc.size());
        		 Properties p = new Properties();
  	  			p.setProperty("PageSize", "A4");	
  	  			AffineTransform at = new AffineTransform();
           			   VectorGraphics g  = new ImageGraphics2D(new File(chartRF, "HWE_.png"),dim_4, ImageConstants.PNG);
           			   g.setProperties(p);
   	                   g.setDeviceIndependent(Constants.plot()==1);
   	                   g.startExport();
        	   for(int l=0; l<rdc.size(); l++){
        	   final int ind1 = l;
               String st = this.names1[ind1];
             
            String datname = "";//(datind>=0) ?  this.names1[datind] : "";
          	
          
               String r_name =  
              	datname+"_"+
              	 st;
               String b_name = datname+"_"+this.names2[ind1];
               try{
              XYSeriesCollection current_r = this.getRSeriesCollection(l) ;
              if(current_r==null) throw new NullPointerException("!!");
             
              final ChartPanel cr = new ChartPanel(graph(current_r ,  r_name, ca_r[l]),
              		 1200, //width
              		 200, //height
              		 1200, //mindrawWidth
                       100, //mindrawHeight
                       1200, //maxDrawWith
                       400,//maxDrawHeight
                       ChartPanel.DEFAULT_BUFFER_USED,
                       true,  // properties
                       true,  // save
                       true,  // print
                       true,  // zoom
                       true   // tooltips		
              
              );
              cr.getChart().getXYPlot().getRangeAxis().setAutoRange(true);
              cr.getChart().getXYPlot().getRangeAxis().setAutoRangeMinimumSize(3.0);
              
         
             cr.setMinimumSize(dim);
              cr.setSize(dim);
         
              if(Constants.plot()>1){
                  chartRNew[l] = cr;
                  if(Constants.plot==2){
                  	cr.getChart().getXYPlot().getDomainAxis().addChangeListener(new AxisChangeListener(){
      				

      				public void axisChanged(AxisChangeEvent event) {
      					
      					ValueAxis axis = ((ValueAxis)event.getAxis());
      					boolean dox = axis==cr.getChart().getXYPlot().getDomainAxis();
      					Range range = axis.getRange();
      					if(dox){
      						//System.err.println("fired property change in "+getName());
      						HWEPlot.this.firePropertyChange("axis", null, range);
      					}
      					changeAxis(range, dox,  ind1);
      				}
                 	
                 });
                  }
                 
              }
              if(plot==2){
                      try{
                      	if(Constants.printPlots())
                      	{
                      	   at.setToTranslation(0, l*dim.getHeight());
      	                   g.setTransform(at);
      	                  cr.print(g);
                    
                   
                      	}
                                               }catch(Exception exc){
                                                  
                                                  
                             exc.printStackTrace();
                         }
                   
                    
                
                  
           
                 
            }
               }catch(Exception exc){
                   System.err.println(this.index+" "+r_name);
                   exc.printStackTrace();
                   
               }
        	   }
        	   g.endExport();
            if(plot<1 && Constants.plot()<2) return;
         
            /*for(int i=0; i<this.probesToPlot.length; i++){
            	  
            	for(int ii=0; ii<probesToPlot[i].length; ii++){
            		if(snp_alias[i][ii]>=0){
            	 XYSeriesCollection current_r = this.rb[i].get(snp_alias[i][ii]);
            	this.chartGNew[i][ii] =  new ChartPanel(graph(current_r, 
            			// this.names1[index.get(i)]+"_"+
            			probesToPlot[i][ii]+"_"+this.location.get(this.snp_alias[i][ii]), ca_r[0]),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            			Constants.scatterWidth(),
            	          false,
            	            true,  // properties
            	            true,  // save
            	            true,  // print
            	            true,  // zoom
            	            true   // tooltips		
            	);
            		}
            //	chartGNew[i][ii].setM
            	}
            	
              
            }*/
            if(Constants.printPlots() && plot==2){
            	/* File out = new File(chartDistF, ("scatter")+".png");
             /   ImageGraphics2D g  = new ImageGraphics2D(out,this.jG, ImageConstants.PNG); 
                g.setDeviceIndependent(Constants.plot()==1);
                g.startExport();
              
                jG.print(g);
                g.endExport();*/
            }
          //  ChartPanel[][] cp = this.getDistCharts();
           
            if(Constants.plot()>1){
                update();
                this.updateUI();
            }
            if(plot==2){
            	 
                try{
               
                   /*     if(gB_!=null){
                      for(int i=0; i<gB_.length; i++){
                    	  if(gB_[i]!=null){
                    	  gB_[i].endExport();
                    	  gR_[i].endExport();
                    	  }
                      }
                        }*/
                }catch(Exception exc){
                    exc.printStackTrace();
                }
                this.plot = 0;
            }
            
        }
        else if(nme.equals("dist_maximisation")){
           
        }
        else if(nme.equals("hmm_maximisation")){
            
        }
    //    else throw new RuntimeException("!! "+nme);
     //   Logger.global.info("exiting "+nme);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
    	
    }
   /*private ChartPanel[][] getDistCharts() {
	ChartPanel[][] res = new ChartPanel[index.size()][];
	for(int i=0; i<res.length; i++){
		res[i] = getDistCharts(i);
	}
	return res
}*/



public  JFreeChart graph(XYSeriesCollection datas_, String title, ColorAdapter ca ) {
    	XYSeriesCollection[] datas = new XYSeriesCollection[]{datas_};
    	AbstractXYItemRenderer[] renderer = new AbstractXYItemRenderer[datas.length];
    	//Shape[] shape = new Shape[datas.length];
    	for(int i=0; i<datas.length; i++){
    		String[] str = Constants.plotType(i);
    		try{
    		renderer[i] = new XYLineAndShapeRenderer();
    			//(AbstractXYItemRenderer)Class.forName("org.jfree.chart.renderer.xy."+str[0]).getConstructor(new Class[0]).newInstance(new Object[0]);
    		}catch(Exception exc){
    			exc.printStackTrace();
    			
    		}
    		if(renderer[i] instanceof XYLineAndShapeRenderer){
    			XYLineAndShapeRenderer rend = (XYLineAndShapeRenderer)renderer[i];
	    	//	if(str[1].equals("nodot")){
	    		//	rend.setBaseShapesVisible(false);
	    		//}
	    		//else{
	    			rend.setBaseShapesVisible(true);
	    		/*	if(str[1].length()==1){
	    				shape[i] = getShape(str[1].charAt(0));
	    			}
	    			else {
	    			//	ShapeUtilities.
	    				shape[i] = new Rectangle(2,2);
	    			}*/
	    			
	    	//	}
	    		//if(str[2].equals("noline")){
	    			rend.setBaseLinesVisible(false);
	    		//}
	    		//else{
	    			//rend.setBaseLinesVisible(true);
	    		//}
    		}
	else{
    			XYBarRenderer rend = (XYBarRenderer) renderer[i];
    			rend.setDrawBarOutline(true);
    			rend.setBase(str.length>1 ? Double.parseDouble(str[1]): 0);
    		}
    	}
        final JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                "Position on chromosome "+Constants.chrom0(), // domain axis label
                "Value ", // range axis label
                datas[0], // data
                PlotOrientation.VERTICAL, true, // include legend
                true, // tooltips?
                false // URL generator? Not required...
                );
      
        XYPlot plot = (XYPlot) chart.getPlot();
        for(int i=0; i<datas.length; i++){
        	  plot.setDataset(i,datas[i]);
        }
      
      
       /* Set<Integer> whiteIndex = new HashSet<Integer>();
        for(int i=datas.getSeriesCount()-1; i>=0 ;i--){
            if(datas.getSeries(i).getItemCount()==0){
               whiteIndex.add(i);
            }
        }*/
        chart.setBorderVisible(false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderPaint(Color.white);
        plot.setBackgroundPaint(Color.white);
     //   whiteIndex.clear();
            NumberAxis yAxis = (NumberAxis)plot.getRangeAxis();
            NumberAxis xAxis = (NumberAxis)plot.getDomainAxis();
          //  yAxis.setRange(new Range(-2.0,2.0))
            chart.getTitle().setFont(font8);
            yAxis.setTickLabelFont(font6);
            yAxis.setLabelFont(font6);
            xAxis.setTickLabelFont(font6);
            xAxis.setLabelFont(font6);
            if(location.size()>0){
	            xAxis.setAutoRange(false);
	            xAxis.setLowerBound(location.get(0)-1000);
	            xAxis.setUpperBound(location.get(location.size()-1)+1000);
            }
         
            
            
           /* if(lu[0]!=null && lu[1]!=null && false){
                plot.mapDatasetToRangeAxis(0, 0);
                plot.setDataset(1,lu[0]);
                plot.setDataset(2, lu[1]);
                plot.mapDatasetToRangeAxis(1, 0);
                plot.mapDatasetToRangeAxis(2, 0);
            }*/
         
            plot.setDomainGridlinePaint(Color.WHITE);
           plot.setRangeGridlinePaint(Color.WHITE);
           // final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
            /*renderer.setShapesFilled(Boolean.TRUE);
            renderer.setLinesVisible(false);*/
     /*       final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
            final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
            final XYBarRenderer renderer3= new XYBarRenderer();
            renderer3.setDrawBarOutline(true);*/
           // renderer3.setBase(1.0);
         //   renderer3.setB
            for(int kk=0; kk<datas.length; kk++){
           // int ik=0;
            plot.setRenderer(kk, renderer[kk]);
            for(int i=0; i<datas[kk].getSeriesCount(); i++){
            	String seriesKey= datas[kk].getSeriesKey(i).toString();
                Color colors = ca.getColor(seriesKey);
             //   if(whiteIndex.contains(i))     colors = Color.WHITE;
                
                renderer[kk].setSeriesPaint(i, colors);
                renderer[kk].setSeriesFillPaint(i,colors);
               
        //        renderer[kk].setSeries
              //[kk].setSeriesShape(i, shape[kk]);
              
              //  ik++;
            }
            }
          /*  renderer1.setShapesVisible(false);
            renderer2.setShapesVisible(false);
            renderer2.setLinesVisible(true);
            renderer2.setLinesVisible(true);
            //renderer1.setStroke(new BasicStroke(BasicStroke.JOIN_MITER));
           // renderer1.setStroke(dashed);
            renderer2.setStroke(dotted);*/
            
           // plot.setRenderer(0,renderer3);
           // plot.setRenderer(1, renderer);
          /*  for(int i=0; i<8; i++){
            	plot.setRenderer(i, renderer);
            }
            if(lu[0]!=null && lu[1]!=null && false){
                plot.setRenderer(1, renderer1);
                plot.setRenderer(2, renderer2);
            }*/
  
              chart.setBackgroundPaint(Color.white);
          //    plot.setBackgroundPaint(Color.black);
             return chart;
}

	// AbstractVectorGraphicsIO     gB,gR;
    final static int numPerPage = 10;
    
   // AbstractVectorGraphicsIO[] gB_; //= new AbstractVectorGraphicsIO[];
 //   AbstractVectorGraphicsIO[] gR_;// = new HashMap<Integer, AbstractVectorGraphicsIO>();

    
    
}
