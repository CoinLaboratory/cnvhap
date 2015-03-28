package lc1.dp.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Shape;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextField;

import lc1.dp.core.BaumWelchTrainer;
import lc1.dp.core.DP;
import lc1.dp.core.Sampler;
import lc1.dp.core.StatePath;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.model.CompoundMarkovModel;
import lc1.dp.model.MarkovModel;
import lc1.dp.states.CompoundState;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.dp.states.PairEmissionState;
import lc1.dp.states.PhasedDataState;
import lc1.dp.states.State;
import lc1.stats.StateDistribution;
import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageGraphics2D;
import org.freehep.graphicsio.emf.EMFGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;

public class HaplotypePanel1 extends JPanel implements PropertyChangeListener{
   
	public final String[] names;
	boolean anonymiseSamples = false;
	final Set<Integer> toAnnotate;
	final boolean[] toAnnotateSnps;
	boolean drawNames = Constants.haplotypeHeight()>=10;
	
	static String imageType = Constants.haploImageType();
	/*public class JPanel1 extends JPanel {
    	
       	 @Override
       	 public void paint(Graphics g){
       		 super.paint(g);
       		   VectorGraphics vg = VectorGraphics.create(g);
       		 vg.setColor(Color.orange);
       		 vg.setLineWidth(1.0);
       		// double height_start = locP.getHeight()+(genP==null ? genP.getHeight() : 0);
       		 for(int i=0; i<rs_pos.length; i++){
       			 if(rs_pos[i]>=0){
       			//double loc_st = locP.getX(rs_pos[i], true);
       			double loc_end = locP.getX(rs_pos[i], false);
       			//if(genP!=null){
       			//	vg.drawLine(loc_st, 0, loc_st, genP.getHeight());
       				
       			//}
       			// vg.drawLine
       			vg.drawLine(loc_end, 0, loc_end, this.getHeight());
       			 }
       		 }
       		 
       	 }
    }*/



	//int no_copies;
    final MarkovModel hmm;
    static boolean replace = false;

 
   static int loc_height1 = 80;
   static int loc_height2 = 20;
  BaumWelchTrainer bwt;
  Dimension dim4, dim4_name;
  int noSnps;
  public static  Font font16 = new Font("SansSerif", Font.PLAIN, 5);
  public static  Font font12 = new Font("SansSerif", Font.PLAIN, 12);
  public static  Font font6 = new Font("SansSerif", Font.PLAIN, 6);
 
 List<Integer> location;

final  List<Character> maj;
final  List<Character> min;
 boolean isImputed(int k, int i){
	 if(k<0) return true;
     return this.bwt.data[k].getFixedInteger(i)==null;
 }

 
 
 final double heig, width;
 final List<Integer> toincl;

 double mult = Constants.haplotypeHeight();//35;
 double mult1 = Constants.haplotypeWidth();
private int[] offset;





final TopleftPanel tlp;
final TopRightPanel trp;
class TopleftPanel extends JPanel{
	private JTextField locP_name;
	private JTextField indexP_name;
	private JTextField genP_name;
	
	
	TopleftPanel(TopRightPanel hapP){
		 setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		 dim4_name = new Dimension(50, (int)Math.round(heig));
		 locP_name = getTextField("location");
	     locP_name.setSize(new Dimension(dim4_name.width, hapP.locP.getSize().height));
	     if(Constants.showEnsThresh()>0){
	    	       genP_name = getTextField("Genes");
	    	     genP_name.setSize(new Dimension(dim4_name.width, hapP.genP.getSize().height));
	     }
	     indexP_name = getTextField("type");
	     indexP_name.setSize(dim4_name);
	     if(genP_name!=null){
		    	
	    	 add(genP_name);
	     }
    	    add(locP_name);
	     add(indexP_name);
	}
	final Font font = new Font("SansSerif", Font.PLAIN, 10);
	public void writeToVg(VectorGraphics gT) {
		
		if(drawNames){
			gT.setFont(font);
			int offset =0;
			if(genP_name!=null){
				genP_name.print(gT);
				offset+=genP_name.getSize().height;  
			}
			at.setToTranslation(0, offset);
			gT.setTransform(at);
		    locP_name.print(gT);
		    offset+=locP_name.getHeight();
		    at.setToTranslation(0, offset);
			gT.setTransform(at);
		    indexP_name.print(gT);
		}
		
	}
}

class TopRightPanel extends JPanel{
	 LocPanel locP ;
	 final LogoPanel indexP;
	 GenePanel genP;
	 int height;
	
	 
	 public TopRightPanel(  List<Integer> loc)
	  throws Exception{
		  setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		   dim4 = new Dimension((int)width,(int) heig);
	      
	      String[] datast1 = bwt.getIndexString();
	     
	   
	    
	     locP= new LocPanel(loc, (int)width, (int)loc_height2,0, "haplotype");
	      if(Constants.showEnsThresh()>0){
	    	     genP=  new GenePanel(Constants.chrom0(),  loc.get(0), loc.get(loc.size()-1),   font10i,null,(int)width);
	    	     Dimension dim5 = new Dimension((int)width, (int)loc_height1);
	    	     genP.setSize(dim5);
	    	     genP.setPreferredSize(dim5);
	    	     add(genP);
	       }
	     add(locP);
	     indexP =  new LogoPanel(-1, 1,false);
	  
	    // Arrays.fill(indexP.state1[0], (short) -1);
	     Arrays.fill(indexP.cert, 1.0);
	    
	     indexP.hap1 = new String[][] {datast1};
	     indexP.setMinimumSize(dim4);
	     indexP.setPreferredSize(dim4);
	     indexP.setSize(dim4);
	  
	     
	     if(genP!=null){
	    	
	    	 height+=genP.getSize().height;
	     }
	   
	  
	    add(locP);
	    add(indexP);
	
	     height+=locP.getSize().height;
	     height+=indexP.getSize().height;
	     
	      
	  }

	public void writeToVg(VectorGraphics gt, int shift) {
		int offset_x = drawNames ? dim4_name.width : 0;
		offset_x +=shift;
			int offset =0;
			if(genP!=null){
				at.setToTranslation(dim4_name.width, offset);
				gt.setTransform(at);
				genP.print(gt);
				offset+=genP.getSize().height;
			}
			at.setToTranslation(offset_x, offset);
			gt.setTransform(at);
		    locP.print(gt);
		    offset+=locP.getHeight();
		    at.setToTranslation(offset_x, offset);
			gt.setTransform(at);
		    indexP.print(gt);
	}
}

class HapPanel extends JPanel{
	
	 JPanel[] jpBS;
	 Map<Integer, LogoPanel> chartG;
	 int[] height;
	 
	 public synchronized LogoPanel removePanel(int l){
	      return chartG.remove(l);
	  }
	 public synchronized LogoPanel getCurrentGPanel(int l){
	      LogoPanel cg = chartG.get(l);
	     boolean hash   = toAnnotate.contains(l);
	      if(cg==null) 
	    {
	          chartG.put(l, cg= new  LogoPanel(l,  bwt.data[l].noCop(), hash));
	      }
	     return cg;
	  }
	 
	 JComponent[] getJFrame( Map<Integer, LogoPanel> cp, 
				JComponent[] jpBS ,
				
				boolean add, int[] height){
		  DataCollection dc1 = DataCollection.datC;
		 //    int height =  dim.height;
		      for(int i=0; i<bwt.data.length; i++){
		         
		        // cp[i].setMinimumSize(new Dimension(dim.width,height*cp.length ));
		          if(toincl.contains(i)){
		        	 
		        	  int index = -1;
		        	
		        		index = ((HaplotypeEmissionState)bwt.data[i]).dataIndex();
		        	
		        	  if(add){
		        		  {
				              LogoPanel cg = new LogoPanel(i,  bwt.data[i].noCop(), toAnnotate.contains(i));
				              cg.setMinimumSize(dim4);
				              cg.setPreferredSize(dim4);
				              cg.setSize(dim4);
				              cp.put(i, cg);
				              jpBS[index].add(cg);
		        		  }
			              
		        	  }
		        	  height[index]+=dim4.height;
		          }
		      }
		 //    this.pack();
		      return jpBS;
		  }
	 
	List<Integer> breaks = new ArrayList<Integer>();
	 public HapPanel(  List<Integer> loc)
	  throws Exception{
		  
		for(int k=1; k<loc.size(); k++){
			if(loc.get(k)-loc.get(k-1)>10000){
				breaks.add(k);
			}
		}
	      
	  String[] datast1 = bwt.getIndexString();
	      setLayout(new BorderLayout());
	      chartG =new HashMap<Integer, LogoPanel>();
	      height = new int[names.length];
	      if(Constants.plot>=2){
	     jpBS = new JPanel[names.length];
	   
	    
	    
	     offset = new int[names.length];
	     Arrays.fill(offset, 0);
	     for(int i=0; i<jpBS.length; i++){
	    	 jpBS[i] = new JPanel();
	    	 jpBS[i].setLayout(new BoxLayout(jpBS[i], BoxLayout.Y_AXIS));
	     }
	    
	      getJFrame(chartG,  jpBS, Constants.plot()>1, height);
	     }
	      else{
	    	  for(int i=0; i<bwt.data.length; i++){
			         if(toincl.contains(i)){
			        	  int index = ((HaplotypeEmissionState)bwt.data[i]).dataIndex();
			        	  height[index]+=dim4.height;
			          }
			      }  
	      }
	  }
	 
	 public void paint(Graphics g){
   	  super.paint(g);
   	  double height = this.getHeight();
	  VectorGraphics vg = VectorGraphics.create(g);
   	  for(int i=0; i<this.breaks.size(); i++){
	   	  double x_start = getXpos(i);
	   
	   	  vg.setColor(Color.RED);
	   	  vg.setLineWidth(10);
	   	  vg.drawLine(x_start, 0, x_start, height);
   	  }
     }


	 
}


class NamePanel1 extends JPanel{
	
	private HashMap<Integer, JTextField> chartG_name;
	private JPanel[] jpBS_name;
	
	public synchronized JTextField getCurrentGPanelName(int l){
	      JTextField cg = chartG_name.get(l);
	      if(cg==null) 
	    {
	          chartG_name.put(l, cg= new  JTextField( bwt.data[l].getName()));
	      }
	     return cg;
	  }
	JComponent[] getJFrame( 
			Map<Integer, JTextField> cp_name,
			 JComponent[] jpBS_name, 
			
			boolean add){
	     
	      for(int i=0; i<bwt.data.length; i++){
	          if(toincl.contains(i)){
	        	  int index = ((HaplotypeEmissionState)bwt.data[i]).dataIndex();
	        	  if(add){
	        		
	        		  {
			              JTextField cg_name = getTextField(bwt.data[i].name);
			              cg_name.setMinimumSize(dim4_name);
			              cg_name.setPreferredSize(dim4_name);
			              cg_name.setSize(dim4_name);
			              cp_name.put(i, cg_name);
			              jpBS_name[index].add(cg_name);
	        		  }
		              
	        	  }
	        	 // height[index+1]+=dim4.height;
	          }
	      }
	 //    this.pack();
	      return jpBS_name;
	  }
	
	public NamePanel1(HapPanel hapP)
	  throws Exception{
	  
		
		   
	      chartG_name = new HashMap<Integer, JTextField>();
	  
	    
	    // Arrays.fill(indexP.state1[0], (short) -1);
	  
	     jpBS_name = new JPanel[names.length];
	
	   
	     offset = new int[names.length];
	     Arrays.fill(offset, 0);
	     for(int i=0; i<jpBS_name.length; i++){
	    	 jpBS_name[i] = new JPanel();
	    	 jpBS_name[i].setLayout(new BoxLayout(jpBS_name[i], BoxLayout.Y_AXIS));
	     }
	     
	     
	    
	   
	      getJFrame( chartG_name, jpBS_name,  Constants.plot()>1);
	    //  getJFrame(chartG_name, jpBS_name, Constants.plot()>1, height);
	    
	     
	     
	   
	     
	  }


	
}
 

 public static  JComponent getSplitPane(JComponent[] jg2, String[] names) {
	if(jg2.length==1) {
		JPanel jres = new JPanel();
		jres.setLayout(new BorderLayout());
		jres.setToolTipText(names[0]);
		jres.setName(names[0]);
		jres.add(BorderLayout.CENTER,
				//new JScrollPane(
						jg2[0]
					//	JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				//JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED)
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
		JSplitPane sjp =   new JSplitPane1(JSplitPane.VERTICAL_SPLIT,  
				getSplitPane(left, leftnames), getSplitPane(right, rightnames));
		sjp.setToolTipText(leftnames[0]);
		 sjp.setDividerLocation((double)left.length / (double)names.length);
		return sjp;
	}
}

  
 
 
 public HaplotypePanel1(BaumWelchTrainer bwt, String[] names1, List<Integer> loc,List<String> snpid,  List<Character> maj,List<Character> min, String name, int no_copies,  File outdir, MarkovModel hmm, int ploidy)
  throws Exception{
    this.setName(name);
  
    toAnnotateSnps = new boolean[snpid.size()];
    if( Constants.rsToAnnotate()!=null ){
    Collection<String> toann = Arrays.asList(Constants.rsToAnnotate());
    for(int i=0; i<snpid.size(); i++){
    	if(toann.contains(snpid.get(i))){
    		toAnnotateSnps[i] = true;
    	}
    	else toAnnotateSnps[i] = false;
    }
    }
    this.toAnnotate = new HashSet<Integer>();
    if(Constants.annotateSamples!=null){
    Collection<String> s1 =Arrays.asList(Constants.annotateSamples);
    for(int i=0; i<bwt.data.length; i++){
    	if(s1.contains(bwt.data[i].getName())){
    		toAnnotate.add(i);
    	}
    }
    }
    this.names = names1;
  
    this.ca = ColorAdapter.get(hmm);
    this.maj = maj;
    this.noSnps = loc.size();
    this.logoF = new File(outdir, "plot_haplotype_"+name+"_"+ploidy);
   
    if(!logoF.exists()) logoF.mkdir();
    else{
       File[] f = logoF.listFiles();
       for(int i=0; i<f.length; i++){
           f[i].delete();
       }
    }
    this.min = min;
    this.reg_length = loc.get(loc.size()-1) - loc.get(0);
      this.hmm = bwt.hmm;
   
      heig = (int) Math.round(mult*Constants.maxPloidy());
      width = noSnps*mult1;
      toincl = new ArrayList<Integer>();
      toincl.clear();
      for(int l=0; l<bwt.data.length; l++){
    	  if(bwt.data[l].noCop()==ploidy)
        	  toincl.add(l);
      }
   
 
      this.bwt = bwt;
      int[] core = getIndices(Constants.core(),loc);
      this.start = core[0];
      this.end = core[1];
      this.location = loc.subList(start, end);
   
    String[] datast =  DataCollection.datC.getUnderlyingDataSets();
  String[] datast1 = bwt.getIndexString();
      setLayout(new BorderLayout());
      this.trp = new TopRightPanel(loc);
      this.tlp = new TopleftPanel(trp);
       hapP = new HapPanel(loc);
      nameP = new NamePanel1(hapP);
    
     offset = new int[datast.length];
     Arrays.fill(offset, 0);
  if(Constants.plot>=2){
     top =   new JSplitPane1(JSplitPane.HORIZONTAL_SPLIT,  
   		tlp, trp);
    bottom =   new JSplitPane1(JSplitPane.HORIZONTAL_SPLIT,  
       		HaplotypePanel1.getSplitPane(nameP.jpBS_name,names), HaplotypePanel1.getSplitPane(hapP.jpBS,names));
    top.setDividerLocation(0.2);
    bottom.setDividerLocation(0.2);
    JSplitPane sjp =   new JSplitPane1(JSplitPane.VERTICAL_SPLIT,  
       		top,
       	new JScrollPane(
       		bottom,
			JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
	JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED
			));
       						
  
      add(BorderLayout.CENTER, sjp);
 
     // this.setSize(dim6);
    bottom.addPropertyChangeListener(this); top.addPropertyChangeListener(this);
  }
  }

 private int[] getIndices(int[] core, List<Integer> loc) {
	int st = core[0];
	int end = core[1];
	int st_i =0;
	int end_i=0;
	
	for(end_i=0; end_i<loc.size(); end_i++){
		int pos = loc.get(end_i);
		if(pos < st) st_i++;
		if(pos > end ){
		
			break;
		}
			
	}
	end_i = end_i-1;
	return new int[] {st_i, end_i};
}



JSplitPane bottom;
 JSplitPane top;
 

  
  private JTextField getTextField(String string) {
	 JTextField tf =  new JTextField(string);
	// tf.setSize(dim4_name);
	 tf.setBackground(Color.white);
	 tf.setFont(font8);
	 tf.setBorder(null);
	 return tf;
}
  
  public static  Font font8 = new Font("SansSerif", Font.PLAIN,(int) (Constants.haplotypeHeight*1.2));
VectorGraphics[] g_;//,g1_;
VectorGraphics[] g_1;
  
  private int get(List<Integer> loc, int i) {
	for(int j=0; j<loc.size(); j++){
		if(loc.get(j)>=i) return j;
	}
	return loc.size();
}

  
 /* JComponent[] getJFrame( Map<Integer, LogoPanel> cp, JComponent[] jpBS, DataCollection dc){
	     
	  //    int height =  dim.height;
	       for(int i=0; i<dc.indiv().size(); i++){
	           HaplotypeEmissionState hes = (HaplotypeEmissionState) dc.data.get(dc.indiv().get(i));
	           int index = ((HaplotypeEmissionState)hes).dataIndex();
	         // cp[i].setMinimumSize(new Dimension(dim.width,height*cp.length ));
	           if(toincl.contains(i)){
	               LogoPanel cg = new LogoPanel(i,  hes.noCop());
	               cg.set(hes);
	               cg.setMinimumSize(dim4);
	               cg.setPreferredSize(dim4);
	               cg.setSize(dim4);
	               cp.put(i, cg);
	               jpBS[index].add(cg);
	           }
	       }
	  //    this.pack();
	       return jpBS;
	   }*/
  
  public void addedInformation1(PhasedDataState[] samples, int ll){
	 // if(usefb) return;
	  LogoPanel l_i = hapP.getCurrentGPanel(ll);
	  JTextField lp_name = nameP.getCurrentGPanelName(ll);
	  EmissionStateSpace[] spaces = ((CompoundEmissionStateSpace)samples[1].getEmissionStateSpace()).getMembers();
      for(int k=0; k<l_i.no_copies; k++){
        
          String[] hap_1 = l_i.hap1[k];
          short[] state_1 = l_i.state1[k];
          for(int i=0; i<samples[0].noSnps(); i++){
              ComparableArray compa =(ComparableArray) samples[0].getElement(i);
              hap_1[i] =((lc1.dp.data.representation.AbstractEmiss) compa.get(k)).toStringShort();
              ComparableArray compaSt =(ComparableArray) samples[1].getElement(i);
             if(!useforwardBackward){
              int sta = spaces[k].getGenotype( compaSt.get(k)).shortValue();
              state_1[i] =(short)(sta+1);
             }
           
          }
      }
     // Logger.global.info("done");
  }

  public static  Font font10= new Font(Font.MONOSPACED, Font.PLAIN, 4);
  public static  Font font10i= new Font("Comic", Font.ITALIC, 10);
 
  
  
  
  final double reg_length;
public final ColorAdapter ca;
  static Color background_color = Color.white;
   
 

 
  public double getXpos(int i){
	  if(Constants.useLocInHap()){
		  return this.trp.locP.getX(i, true);
	  }
	  else
	  return i*mult1;
  }
  
 final int start, end;
  class LogoPanel extends JPanel{
      private String[][] hap1; //first index is by individual; second is first or second haplotype
      double[] cert;
      short[][] state1;
      
     // String name;
   //   double y_offset = 
    //	  Constants.haplotypeHeight < 10 ? 0 :  
    //	  (int) Math.floor((((double)Constants.haplotypeHeight)/3.0));
      int index;
      int no_copies;
      int data_index;
      public void set(HaplotypeEmissionState hes){
    	  CompoundEmissionStateSpace emStSp = Emiss.getSpaceForNoCopies(no_copies);
    	  EmissionStateSpace emStSp1 = emStSp.getMembers()[0];
    	  for(int i=0; i<hes.noSnps(); i++){
    		  int[] indices =  emStSp.getMemberIndices(hes.getBestIndex(i));
    		 
    		  for(int j=0; j<no_copies; j++){
    			  hap1[j][i] = ((Emiss)emStSp1.get(indices[j])).toStringShort();
    			  cert[i] =1.0;
    			  int cn = emStSp1.getCN(indices[j]);
    			  state1[j][i] =(short) (cn==0 ? 1 : (cn==1 ? 3 : 2));
    		  }
    	  }
      }
      final boolean hash;
      LogoPanel(int i,int no_copies, boolean hash){
          this.index = i;
          this.hash = hash;
          data_index =i<0 ? -1 : ((HaplotypeEmissionState)bwt.data[i]).dataIndex();
          this.no_copies = no_copies;
          //this.name = i>=0 ? names[data_index]+"_"+name : name;
          this.setMinimumSize(dim4);
          this.setSize(dim4);
          this.cert = new double[noSnps];
          this.hap1 = new String[no_copies][noSnps];
          this.state1 = new short[no_copies][noSnps];
          this.setSize(dim4);
         // if(this.hap1[0].length<5){
        //	  throw new RuntimeException("!!");
         // }
      }
      Character uncertain = '-';
     
      
      @Override
      public void paint( Graphics g ) {
          g.setColor(Color.white);
          g.fillRect(0,0, this.getWidth(), this.getHeight());
          VectorGraphics vg = VectorGraphics.create(g);
       //   int noSnps = hap1[0].length;
      //   double  wid = mult1;//((double)this.getWidth()  ) / ((double) noSnps);
         double height  = ((double)this.getHeight()  ) / (double) no_copies;
         vg.setColor(Color.black);
         vg.setFont(font16);
        double base = this.getHeight();
        Font font = 
       			  font10 ;//: font10i;
     
         FontRenderContext frc = vg.getFontRenderContext();
      
         for(int i=0; i<noSnps; i++){
        	 
        	double x_start = getXpos(i);
        	
        	double widx;
        	if(i==noSnps-1){
        		
           	 widx = i==0 ? 5 : x_start -getXpos(i-1) ;  //assume width is same as previous for last
        	}
        	else{
        		double x_end = i==noSnps-1 ? 1.0:getXpos(i+1);
        	 widx = x_end -x_start;
        	}
             for(int j=0; j<no_copies; j++)
             {
            	 try{
            	 if(Constants.showHapAllele()){
            		 
                 String  ch = hap1[j][i];
                 if(ch==null) continue;
               //  if(ch.equals("_")) continue;
                 Character majo, mino;
                 if(maj==null || maj.size()==0){
                     majo = mino = null;
                 }
                 
                 else{
                     mino = min.get(i);
                     majo = maj.get(i);
                 }
                 if(replace && majo!=null && !majo.equals(uncertain) && !mino.equals(uncertain)){
                     ch = ch.replace('A', majo).replace('B', mino);
                 }
                 ch = ch.replace("_", "-");
                 double ht1 = height / ((double)ch.length());
                for(int k=0; k< ch.length(); k++){
                 if(cert[i]>Constants.hapl_cert_thresh()){
	                 GlyphVector gv = font.createGlyphVector(frc, new char[] {ch.charAt(k)});
	                 Shape outline = gv.getOutline();
	                 Rectangle2D obounds = outline.getBounds2D();
	                 AffineTransform at = new AffineTransform();
	                 double extra =  ((double)k)*ht1;
	                 at.setToTranslation(x_start, base - (j+1)*height + extra ); 
	                 //1
	                 double sc = cert[i];
	                 //sc = 1;
	                 at.scale((widx/obounds.getWidth()), (ht1/obounds.getHeight())*sc);
	                 at.translate(-obounds.getMinX(), -obounds.getMinY());
	                 vg.setColor(ca.getColor3((state1[j][i])));
	                 outline  = at.createTransformedShape(outline);
	                 vg.fill(outline);
	                 vg.draw(outline);
                 }
                }
             }else{
            	 vg.setColor(ca.getColor((state1[j][i]-1)));
            	// vg.f
            	 vg.setLineWidth(0);
            	 vg.fillRect(x_start, base -(j+1)*height, widx, height);
        	 }
            	 }catch(Exception exc){
            		 exc.printStackTrace();
            	 }
        	 }
             vg.setColor(Color.LIGHT_GRAY);
             vg.setLineWidth(Constants.haplotypeHeight/10.0);
             int y1 =0;
             int y2 = y1+this.getHeight();
             double wid = Math.min(this.getWidth(),getXpos(noSnps-1)+1 );
             (vg).drawLine(0, y1,wid , y1);
             (vg).drawLine(0, y2,  wid, y2);
             if(hash || toAnnotateSnps[i] ){
         		vg.setColor(Color.white);
         		 vg.setLineWidth(1);
         		vg.drawLine(x_start,0, x_start+widx, this.getHeight());
         		
         	}
         }
      }
  
  }
  
 
 
 // final LocPanel locP;
 /* public void update(int i){
      JComponent currentGPanel = this.getCurrentGPanel(i);
      if(currentGPanel==null) return;
       currentGPanel.setMinimumSize(dim4);
        currentGPanel.removeAll();
        if(chartGNew[i]!=null){
       currentGPanel.add(this.chartGNew[i]);
        }
  }*/
  /*public void update(){
      for(int i=0; i<chartG.length; i++){
       
      }
  }*/
 /* public LogoPanel getCurrentLogoPanel(int i){
      if(chartGNew!=null){
          return chartGNew[i];
      }
      else{
          return new LogoPanel(i, bwt.data[i].name);
      }
  }
  public void addChart(final int i){
      if(((HaplotypeEmissionState)bwt.data[i]).dataIndex()==this.data_index){
          Dimension dim4 = new Dimension(this.noSnps*20, 120);
          
       LogoPanel cg = new LogoPanel(i, bwt.data[i].name);
           cg.setMinimumSize(dim4);
           cg.setPreferredSize(dim4);
             
            
              chartGNew[i] = cg;
      }
     // }
  }*/
  final File logoF;
 /* public void writeToFile(){
      try{
          
          this.writeToZipFile(this.jpBS, logoF, this.getName());
         
      }catch(Exception exc){
          exc.printStackTrace();
      }
  }*/
  VectorGraphics g = null;
  
  public boolean plot = false;
 
 
  
  

  double[] res;
  /*private synchronized double[] getProbOverStates(StateDistribution emissionC,
          MarkovModel hmm, HaplotypeEmissionState obj, int i) {
	  EmissionStateSpace emstsp =Emiss.getSpaceForNoCopies(obj.noCop());
      if(res==null)
    	  res = new double[emstsp.defaultList.size()];
      Arrays.fill(res,  0.0);
      double sum = 0;
      for(int j=1; j<emissionC.dist.length; j++){
          double p = emissionC.dist[j];
          EmissionState state_j = (EmissionState) hmm.getState(j);
          EmissionStateSpace emstsp1 = state_j.getEmissionStateSpace();
          if(p>0){
              sum+=obj.calcDistribution((EmissionState) state_j, i,  p);
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
  final HapPanel hapP;
  final NamePanel1 nameP;
  public void addedInformation(StateDistribution emissionC, int ll, int i){
      LogoPanel lp =  this.hapP.getCurrentGPanel(ll);   
      JTextField lp_name = nameP.getCurrentGPanelName(ll);
	  double[] cert_i =  lp.cert;
	
       //   Arrays.fill(cert_i, 1.0);
	  HaplotypeEmissionState sta = (HaplotypeEmissionState) bwt.data[ll];
      
          if(Constants.modelbg() && ll==0){
        	  double[] prbs = (sta.emissions[i]).probs();
        	  int maxind1 = Constants.getMax(prbs);
        	  State st = this.hmm.getState(Constants.getMax(emissionC.dist));
        	 EmissionState[] states =  ((CompoundState)st).getMemberStates(true);
        	 double cert = prbs[maxind1];///Constants.sum(emissionC.dist);
        	 cert_i[i] = cert;
        	 for(int k=0; k<states.length; k++){
        		 lp.state1[k][i] = (short)states[k].getIndex();
        	  lp.hap1[k][i] = states[k].noCop()+"";
        	  
        	 }
          }
          else{
        	  if(useforwardBackward){
	        	  int max_ind = Constants.getMax(emissionC.dist);
	        	  State st = this.hmm.getState(Constants.getMax(emissionC.dist));
	         	 EmissionState[] states =  ((CompoundState)st).getMemberStates(true);
	         	
	         	 for(int k=0; k<states.length; k++){
	         		 lp.state1[k][i] = (short)states[k].getIndex();
	         	   //  lp.hap1[k][i] = states[k].noCop()+"";
	         	    cert_i[i] = emissionC.dist[max_ind];
	         	 }
        	  }
        	  else{
        	  double[] prob =  PairEmissionState.pool.getObj(Emiss.getSpaceForNoCopies(sta.noCop()).genoListSize());
           	 Sampler.getProbOverStates(emissionC, bwt.hmm, sta, i,prob, Constants.isLogProbs());
              int maxind = Constants.getMax(prob);
              cert_i[i] = prob[maxind];
              PairEmissionState.pool.returnObj(prob);
        	  }
          }
  }
  static boolean useforwardBackward = Constants.fbHapPanel();
 int cnt =0;
 int cnt1 = 0;
 
 
 //boolean printAll = true;
 AffineTransform at = new AffineTransform();
 
 
 public synchronized void propertyChange(PropertyChangeEvent arg0) {
	  try{
      String nme = arg0.getPropertyName();
      if(nme.equals( JSplitPane.DIVIDER_LOCATION_PROPERTY)){
    	Object src =  arg0.getSource();
    	int location = ((Integer) arg0.getNewValue());
    	if(src==top ){
    		if(bottom.getDividerLocation()!=location) bottom.setDividerLocation(location);
    	}
    	else{
    		if(top.getDividerLocation()!=location) top.setDividerLocation(location);
    	}
    	return;
      }
     
      else  if(nme.equals("setToPlot")){
          int level = (Integer) arg0.getNewValue();
       
          if(level==2){
        	  System.err.println("SETTING TO PLOT");
        	  setToPlot();
        	  if(Constants.printPlots()){
		        	   Properties p = new Properties();
		        	  g_ = new VectorGraphics[names.length];
		        
		      	    for(int i=0; i<g_.length; i++){
		      	    	int height = hapP.height[i];
		      	    	if(height>0){
		      	    	g_[i] = getVectorGraphics(logoF, "hapl_"+names[i],new Dimension(dim4.width+(this.drawNames ? dim4_name.width : 0), height));
		      	    			
		      	    		
		      	    	g_[i].setProperties(p);
		      			g_[i].startExport();
		      	    	}
		      	    }
		      	  
        	  }
          }
         
//          if(plot){
           
  //        }
          return;
      }
      if(Constants.plot()<=1 && ! plot) return;
      
     if(nme.equals("emiss")){
    	  if(Constants.plot()<=1 && ! plot) return;
          Object[] obj = (Object[]) arg0.getNewValue();
          Integer l = (Integer)obj[1];
          StateDistribution dist = (StateDistribution) obj[0];
          Integer i = (Integer)obj[2];
            addedInformation(dist, l, i);
           // update(i);
      }
 
     else if(nme.equals("pre_exp")){
      //   this.setup();
     }
      else if(nme.equals("finished")){
          Object[] obj = (Object[]) arg0.getNewValue();
        
          Integer l = (Integer)obj[1];
         
          if(plot ){
        	 /* if(Constants.printPlots() && l==0 && SignificancePlot.snp_alias.length<this.noSnps && Constants.showScatter()){
        	    	 Properties p = new Properties();
        	    	int start = SignificancePlot.snp_alias[SignificancePlot.start_];
              	int end= SignificancePlot.snp_alias[SignificancePlot.end_];
              	int width1 =(int) Math.round( mult*(end-start));
        	    //	g1_ = new VectorGraphics[names.length];
     	      	   // for(int i=0; i<g1_.length; i++){
     	      	   // 	g1_[i] = new ImageGraphics2D(new File(logoF, "hapl_"+names[i]+"_sig.png"),new Dimension(width1+(this.drawNames ? dim4_name.width : 0), hapP.height[i]),ImageConstants.PNG);
     	      	   // 	g1_[i].setProperties(p);
     	      		//	g1_[i].startExport();
     	      	   // }
        	    }*/
        	  int index = ((HaplotypeEmissionState)bwt.data[l]).dataIndex();
           //   if(!this.bwt.dataIndex(l).equals(this.data_index)) return;
        	  {
	              LogoPanel panel ;
	              JTextField field;
	              if(Constants.plot()<2){
	                  panel = hapP.removePanel(l);
	                  field = nameP.chartG_name.remove(l);
	              }
	              else{
	                  panel = hapP.getCurrentGPanel(l);
	                  field = nameP.getCurrentGPanelName(l);
	              }
	              panel.setSize(dim4);
	              field.setSize(dim4_name);
	              {
		              VectorGraphics vg = this.g_[index];
		             
		              if(vg!=null){
		             if(drawNames){ 
		            	 vg.setFont(this.font8);
		            	 at.setToTranslation(0, offset[index]);
		            	 vg.setTransform(at);
		            	 vg.setColor(Color.gray);
		            	 vg.fillRect(0, 0, field.getWidth(), field.getHeight());
		            	 vg.setColor(Color.white);
		            	 
		            	 vg.drawString(anonymiseSamples ? (l+1)+"" : field.getText(), 0,field.getHeight()*0.8);
		             }
		              at.setToTranslation(drawNames ? dim4_name.width : 0, offset[index]);
		              vg.setTransform(at);
		              panel.print(vg);
	              }
	              }
	            /*  if(g1_!=null && g1_[index]!=null){
	            	    int start = SignificancePlot.snp_alias[SignificancePlot.start_];
	            	    VectorGraphics vg = this.g1_[index];
	            	   
			              at.setToTranslation((drawNames ? dim4_name.width : 0)-start, offset[index]);
			              vg.setTransform(at);
			              panel.print(vg);
			              if(drawNames){
				              at.setToTranslation(0, offset[index]);
				              vg.setTransform(at);
				              field.print(vg);
		            	    }
	              }*/
	              
        	  }
              offset[index]+=dim4.height;
            
              }
            
        
       //   addChart(l);
         
      }
      else if(nme.equals("expec_i")){
    	  if(plot || Constants.plot()>=2){
	          Object[] obj = (Object[]) arg0.getNewValue();
	          Integer ll = (Integer)obj[0];
	          String name = bwt.data[ll].getName();
	          if(!Constants.modelbg() || ll>0){
	        	  CompoundEmissionStateSpace emStSp = Emiss.getSpaceForNoCopies(bwt.data[ll].noCop());
	        	
	        	  PhasedDataState[] sam = new PhasedDataState[] {
	                  SimpleScorableObject.make(name, noSnps, emStSp,(short)-1),
	                  SimpleScorableObject.make(name, noSnps, 
	                		  this.hmm.getStateSpace()
	                		,(short)-1)
	                 
	          };
	         
	        
	          DP dp = (DP)obj[1];
	          StatePath sp = dp.getStatePath(false);
	          Sampler.sample(sp, (CompoundMarkovModel)hmm, false, sam);
	          this.addedInformation1(sam, ll);
	          }
    	  }
      }
      else if(nme.equals("expectation1")){ //end of expectation calculation
          //update();
          if(plot){
              try{
            	  for(int i=0; i<g_.length; i++){
            		  if(g_[i]!=null)
            		  g_[i].endExport();
            	  }
            	 /* if(g1_!=null){
            	  for(int i=0; i<g1_.length; i++){
            		  if(g1_[i]!=null) g1_[i].endExport();
            	  }
            	  }*/
              }catch(Exception exc){
                  exc.printStackTrace();
              }
              
              if(Constants.printPlots()){
            	
         	    Properties p = new Properties();
         	    AffineTransform at = new AffineTransform();
        		p.setProperty("PageSize", "A4");	
        		{
               VectorGraphics gT = this.getVectorGraphics(logoF, "genes", 
        			new Dimension(dim4.width+(this.drawNames ? dim4_name.width : 0), this.trp.height));
        		gT.setProperties(p);
        		gT.startExport();
        			if(this.drawNames) this.tlp.writeToVg(gT);
        			trp.writeToVg(gT,0);
        	
        		
        	    gT.endExport();
        		}
        		try{
        			int start=0;
        			int end;
        			if(SignificancePlot.snp_alias!=null){
        				start = SignificancePlot.snp_alias[SignificancePlot.start_];
                    	 end= SignificancePlot.snp_alias[SignificancePlot.end_];
        			}
        			else end = DataCollection.datC.size();
        		
                	int width1 =(int) Math.round( mult*(end-start));
                	Dimension d = new Dimension(width1+(this.drawNames ? dim4_name.width : 0), this.trp.height);
                	if(d.width>0 && d.height>0){
	                	 VectorGraphics gT = this.getVectorGraphics(logoF, "genes1", d);
//	             			new ImageGraphics2D(new File(logoF, "genes1.png"),d,ImageConstants.PNG);
	             		gT.setProperties(p);
	             		gT.startExport();
	             		
	             			trp.writeToVg(gT,-start);
	             			if(this.drawNames) this.tlp.writeToVg(gT);
	             		
	             	    gT.endExport();
                	}
        		}catch(Exception exc){
        			exc.printStackTrace();
        		}
        	   
        	    
        	    
              }
              
            this.plot = false;
          }
          this.updateUI();
      }
      else if(nme.equals("done")){
    	//  this.setup1();
      }
	  }catch(Exception exc){
		exc.printStackTrace();  
	  }
      
  }
private VectorGraphics getVectorGraphics(File logoF2, String string,
		Dimension dimension) throws Exception {
	File outF = new File(logoF2, string+"."+this.imageType);
	
	
	 if(imageType.equals("pdf")){
		return new PDFGraphics2D(outF, dimension);
	}
	else if(imageType.equals("emf")){
		return new EMFGraphics2D(outF, dimension);
	}
	else{
		return new ImageGraphics2D(outF,dimension,this.imageType.toUpperCase());
	}
}




public void setToPlot() {
   this.plot = true;
    
}
  
}

