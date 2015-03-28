package lc1.dp.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Properties;
import java.util.logging.Logger;

import javax.swing.JFrame;
import javax.swing.JPanel;

import lc1.util.Constants;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.emf.EMFGraphics2D;
import org.freehep.graphicsio.pdf.PDFGraphics2D;
/*@Author Lachlan Coin*/
public class DomainPanel extends JPanel{
    
    
    static WindowAdapter APP_CLOSER = new WindowAdapter() {
        public void windowClosing(WindowEvent e) {
            System.exit(0);
        }
    };
    static double[][] getMath(int len){
        double[][] res = new double[len][len];
        for(int i=0; i<len; i++){
            for(int j=0; j<len; j++){
                res[i][j] =Constants.rand.nextDouble();
            } 
        }
        return res;
    }
    final int start; final int end;
    static class ConfigFrame extends JFrame {
        static final long serialVersionUID = 34525;

       

        ConfigFrame(DomainPanel domP, String title) {
         super(title);
           // JTabbedPane config = new JTabbedPane();
            getContentPane().setLayout(new BorderLayout());
            getContentPane().add(domP, BorderLayout.CENTER);
         //   this.setMinimumSize(dim);
           // this.setSize(dim);
           
         
          //  config.addTab("LD", domP );
           // config.setSize(500,500);
           
            addWindowListener(APP_CLOSER);
            
        }
    }
    public static void plot(double[][] mat, int st, int end, String title, boolean show) {
        try{
            
            DomainPanel domP = new DomainPanel(mat, st, end);
            domP.setPreferredSize(dim);
            if(show){
                ConfigFrame frame = new ConfigFrame(domP, title);
                frame.setVisible(show);
                frame.pack();
            }
          // frame.setPreferredSize(dim);
          
            Properties p = new Properties();
            p.setProperty("PageSize","A4");
          boolean pdf = false;
            VectorGraphics g1 = pdf ?new PDFGraphics2D(new File(title+".pdf"), dim): new EMFGraphics2D(new File(title+".emf"), dim); 
            g1.setProperties(p); 
            g1.startExport(); 
            domP.print(g1); 
            g1.endExport();
            //frame.
            }catch(Exception exc){
                exc.printStackTrace();
            }
    }
  
    public static void main(String[] args){
      plot(getMath(40), 20, 30, "test", true);
    }

    Color  background_color = Color.WHITE,
            line_color = Color.BLACK,
            font_color = Color.BLACK;
    int x_start =0;
    int x_end;

    Font small_font,
        large_font,
        small_italic_font,
        large_italic_font;
    FontMetrics fm_small,
                fm_large,
                fm_small_italic,
                fm_large_italic;

    Logger logger = Logger.global.getAnonymousLogger();
   // int[] y_loc;
  //  int domain_thickness = 30;
   // double ratio;
    Color[] colors = new Color[] {Color.RED, Color.GREEN, Color.CYAN, Color.BLUE, Color.PINK, Color.YELLOW, Color.MAGENTA, Color.ORANGE};
    void mediumFonts() {
        small_font        = new Font( "Helvetica", Font.PLAIN,  10 );
        large_font        = new Font( "Helvetica", Font.PLAIN,  11 );
        small_italic_font = new Font( "Helvetica", Font.ITALIC, 10 );
        large_italic_font = new Font( "Helvetica", Font.ITALIC, 11 );
        fm_small        = getFontMetrics( small_font );
        fm_large        = getFontMetrics( large_font );
        fm_small_italic = getFontMetrics( small_italic_font );
        fm_large_italic = getFontMetrics( large_italic_font );
    }
    double[][] mat;
    
    
    
    DomainPanel(double[][] d, int start ,int end){
        
  //      this.set
        this.start = start;
        this.end = end;
        this.mat = d;
        mediumFonts();
    }
    int index=0;
    
   int[] y_loc;
  //  double offset = 0;
   // int domain_thickness = 10;
 //  double ratio;
    static int len = 500;
    static int width = 500;
   double wid ;
    static Dimension dim = new Dimension(len,width);
   
    public void paint( Graphics g ) {
        g.setColor(this.background_color);
        g.fillRect(0,0, getWidth(), getHeight());
        VectorGraphics vg = VectorGraphics.create(g);
        g.setColor(Color.BLACK);
        g.drawLine(0,0,len, width);
        double wid =  ((double)width ) / ((double) mat.length);
        for(int i=0; i<mat.length; i++){
            double x_start =  ((double)i) *wid;
            for(int j=0; j< mat[i].length;
            j++){
                double y_start = ((double)j) *wid;  
                Color c = //i==pos || j==pos  ? Color.GREEN: 
                    Color.RED;
                Color newC = new Color(c.getRed(),  c.getGreen(), c.getBlue(),(int)Math.round(mat[i][j]*255));
                vg.setColor(newC);
                vg.fillRect(x_start, y_start, wid, wid);
                /*if(i==pos || j==pos){
                    Color c1 =Color.GREEN;
                    //Color newC1 = new Color(c.getRed(),  c.getGreen(), c.getBlue(),(int)Math.round(mat[i][j]*255));
                    vg.setColor(c1);
                    vg.drawRect(x_start, y_start, wid, wid);
                }*/
            }
        }
        Color c1 =Color.GREEN;
        //Color newC1 = new Color(c.getRed(),  c.getGreen(), c.getBlue(),(int)Math.round(mat[i][j]*255));
        vg.setColor(c1);
        
        vg.drawRect(wid*start, 0, wid*(end-start+1), width);
        vg.drawRect(0, wid*start, width, wid*(end-start+1));
        
              //vg.drawLine(100.0,100.0,200.0,200.0);
    }
 
}
