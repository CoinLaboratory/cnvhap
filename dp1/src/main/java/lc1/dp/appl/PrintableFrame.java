package lc1.dp.appl;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.io.File;
import java.util.Properties;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuItem;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.PageConstants;
import org.freehep.graphicsio.emf.EMFGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
import org.jfree.chart.ChartPanel;

public class PrintableFrame {
	 public class MyMenuItem extends JMenuItem implements ActionListener{

			public MyMenuItem(String string) {
				super(string);
				this.addActionListener(this);
			}

			public void actionPerformed(ActionEvent e) {
				filechooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
				 filechooser.setMultiSelectionEnabled(false);
				 filechooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
				 filechooser.showOpenDialog(null);
				 File outFile =  filechooser.getSelectedFile();
				 try{
				writeToFile(outFile);
				 }catch(Exception exc){
					 exc.printStackTrace();
				 }
			}

			}
	 JFrame jframe;
	 ChartPanel[] cp;
	 
	 PrintableFrame(ChartPanel[] cp){
		 
	 }
	 static JFileChooser filechooser = new JFileChooser();
	 public void writeToFile(File out) throws Exception{
			boolean emf=false;
	    	if(out.getName().endsWith(".emf")) emf  = true;
	    	else if (!out.getName().endsWith(".pdf")){
	    		out = new File(out.getParent(), out.getName()+".pdf");
	    	}
	    	//(Iterator<String> phen = pheno.iterator(); phen.hasNext();i++){
	    		//String phent = phen.next();
	    	
	    	
	        	
	        		  Properties p = new Properties();
	                  p.setProperty("PageSize","A4");
	              //    ExportDialog ep;
	                  p.setProperty(PageConstants.ORIENTATION, PageConstants.LANDSCAPE);
	                  VectorGraphics g = emf ? new EMFGraphics2D(out, cp[0].getSize()) : new PDFGraphics2D(out, 
	                		  jframe.getContentPane().getSize()); 
	                  
	                  g.setProperties(p); 
	                  g.startExport(); 
	                 /* if(top.gp!=null){top.gp.print(g); 
	                  AffineTransform at = new AffineTransform();
	                  at.setToTranslation(0,top.gp.getHeight());
	                  g.setTransform(at);
	                  }*/
	                  
	                 // gp.print(g); 
	                 
	                  int height =cp[0].getHeight();
	                  for(int ii=0; ii<cp.length; ii++){
	                	  AffineTransform at = new AffineTransform();
		                  at.setToTranslation(0,height);
		                  g.setTransform(at);
		                  cp[ii].print(g);
		                  height+=cp[ii].getHeight();
	                  }
	                  g.endExport();
	/*            EMFGraphics2D     g = new EMFGraphics2D(new File(out, top.getName()+".emf"),top);//ImageConstants.PNG); 
	            g.startExport(); 
	            app_l.get(i).armitage[k].print(g); 
	            g.endExport();*/
	        	
	    	//}
	    }
}
