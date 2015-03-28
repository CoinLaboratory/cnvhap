package lc1.CGH;




public class AgilentProbeReader extends AbstractAberatiionReader {
    
public AgilentProbeReader(long lengthLim) {
        super(lengthLim, "");
        // TODO Auto-generated constructor stub
    }
String[] cols =  new String[] {"ChromNum", "GeneName"};
    
 //   ProbeUID    Name    ChromNum    GeneName    GeneSymbol  Sequence    SystematicName  Start   IsExonic    InTCAG  IntervalType    Is185kProbe
    
    @Override
    public String[] getCols() {
       return cols;
    }
       
    public String getName(String[] str){
        return  "";
     }
     public String getChr(String[] str){
         return str[col[0]];
     }
     public String getStart(String[] str){
         return str[col[1]].split(":")[1].split("-")[0];
     }
     public String getEnd(String[] str){
         return str[col[1]].split(":")[1].split("-")[1];
     }
     public int getNoProbes(String[] str){
         return -1;
     }  
     public  String getProbeId(String[] str){
         return str[1];
     }
    
     @Override
     public boolean exclude(String[] str){
         return false;
     }
     @Override
     public double getNoCopy(String[] str){
         return 0;
                 
      }
    
   /* AgilentProbeReader (File dir, String chromosome, Location region)throws Exception{
        File f1 = new File(dir, "IC_Custom244K_238459_ProbesTable.txt.gz" );
        BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f1))));
        String st = br.readLine();
        String chr = "Chr";
        String[] names = st.split("\\s+");
        while((st = br.readLine())!=null){
            String[] str = st.split("\\t");
            String[] loc = str[loc_int].split(":");
            String chrom = loc[0].substring(3);
            if(chromosome!=null && !chrom.equals(chromosome)) continue;
           String[] pos = loc[1].split("-");
            long min =Long.parseLong(pos[0]);
           long max = Long.parseLong(pos[1]);
            Location location = max > min ? new Location(chrom, min,max) :new Location(chrom, max,min) ;
            if(region==null || region.overlaps(location)>0){
                add(location);
            }
        }
        br.close();
        this.sort();
    }*/
}
