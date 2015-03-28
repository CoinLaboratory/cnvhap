package lc1.ensj;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

public class SummaryFileReader {
    public static void main(String[] args){
        SummaryFileReader sfr = new SummaryFileReader(new File("."));
    }
    public SummaryFileReader(File dir){
        File[] f = dir.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
               return pathname.isDirectory();
            }
            
        });
        for(int i=0; i<f.length; i++){
            try{
            read(new File(f[i], "summary4.txt"));
            }catch(Exception exc){
                exc.printStackTrace();
            }
        }
        siteRes1.normalise();
        siteRes2.normalise();
        System.err.println(siteRes1.toString());
        System.err.println(siteRes2.toString());
            double tot  = this.m.get(0.5).res1[1];
        for(Iterator<OverallResults> it = m.values().iterator(); it.hasNext();){
            OverallResults oe = it.next();
            oe.res1[1] = tot;
            System.err.println(oe);
        }
    }
    
    Map<Double, OverallResults> m = new TreeMap<Double, OverallResults>();
    SiteResults siteRes1 = new SiteResults(1);
    SiteResults siteRes2 = new SiteResults(2);
    int cnt =0;
        public void read(File f) throws Exception{
        
            BufferedReader br = new BufferedReader(new FileReader(f));
            String st = "";
            double tot =0;
            while((st = br.readLine())!=null){
                if(st.startsWith("thresh")){
                    double thresh = Double.parseDouble(st.split("\\s+")[2]);
                  
                    for(int i=0; i<6; i++){
                        st = br.readLine();
                    }
                    OverallResults oe = m.get(thresh);
                    if(oe==null){
                        m.put(thresh, oe = new OverallResults(thresh));
                    }
                    oe.append(st.split("\\s+"));
                }
                else if(st.startsWith("cnt thresh 1")){
                    siteRes1.append(br);
                    br.readLine();
                    siteRes2.append(br);
                    cnt++;
                }
            }
            br.close();
        }
        
        
       
        class OverallResults{
           final double threshold;
           OverallResults(double th){
               this.threshold = th;
           }
            double[] res1 = new double[] {0,0}; //power
            double[] res2 = new double[] {0,0}; //fp
            double[] res3 = new double[] {0,0}; //misclass
            public void append(String[] st){
                for(int i=0; i<st.length; i++){
                    st[i] = st[i].replace(':', ' ');
                }
                double d1 =  Double.parseDouble(st[2].trim());
                res1[0] += (d1-Double.parseDouble(st[1].trim()));
                res1[1] +=d1;
                res2[0] += Double.parseDouble(st[4].trim());
                res2[1] += Double.parseDouble(st[5].trim());
                res3[0] += Double.parseDouble(st[7].trim());
                res3[1] += Double.parseDouble(st[8].trim());
               // res1[0] = res1[1] -res1[0];
            }
            public String toString(){
                StringBuffer sb = new StringBuffer();
                sb.append(this.threshold+": "+res1[0]/res1[1]+" // "+res2[0]/res2[1]+" //  "+res3[0]/res3[1]);
                return sb.toString();
            }
        }
        class SiteResults{
            int cntThresh;
            double[] thresh = new double[] {0.5, 0.6, 0.7, 0.8, 0.9, 0.95};
            double[][] res = new double[6][2];
            int[] count = new int[6];
            SiteResults(int cnt){
                for(int i=0; i<res.length; i++){
                    Arrays.fill(res[i], 0.0);
                }
                Arrays.fill(count, 0);
                this.cntThresh = cnt;
            }
            public void normalise() {
                for(int i=0; i<res.length; i++){
                    res[i][0] = res[i][0]/count[i];
                    res[i][1] = res[i][1]/count[i];
                }
                
            }
            public void append(BufferedReader br) throws Exception{
                boolean failed = false;
                for(int i=0; i<res.length; i++){
                    String[] str = br.readLine().split("\\s+");
                    try{
                    double d1 = Double.parseDouble(str[11]);
                    double d2 = Double.parseDouble(str[12]);
                   
                    
                    if(Double.isNaN(d1) || Double.isNaN(d2)){
                       Logger.global.warning("!! "+Arrays.asList(str));
                    }
                    else{
                    res[i][0] += d1;
                        res[i][1] += d2;
                        count[i]++;
                    }
                    }catch(Exception exc){
                        exc.printStackTrace();
                    }
                        
                }
                
            }
            public String toString(){
                StringBuffer sb = new StringBuffer("cnt thresh "+cntThresh+"\n");
                for(int i=0; i<res.length; i++){
                   sb.append(thresh[i]+": "+res[i][0]+" "+res[i][1]+"\n");
                }
                return sb.toString();
            }
        }
}
