package lc1.CGH;

import java.io.File;
import java.util.Set;


public class MultiAberationReader extends Locreader {
    public MultiAberationReader( File dir, String chromosome, Location region, Set<String> indiv)throws Exception{
        this(dir, chromosome, region, ExcelAberationReader.class, indiv);
    }
    
    public MultiAberationReader( File dir, String chromosome, Location region, Class clazz, Set<String> indiv)throws Exception{
       super(Long.MAX_VALUE,"");
     //   TarInputStream dir1;
        boolean exists = dir.exists();
       /* if(dir.exists()){
            dir1 = new TarInputStream(new FileInputStream(dir));
        }
        else{
            dir1 = new TarInputStream(new GZIPInputStream(new FileInputStream(new File(dir.getAbsolutePath()+".gz"))));
        }
      
      //   TarEntry nxt =  dir1.getNextEntry();
         BufferedReader br=null;
        while(nxt!=null){
            if(nxt.getName().endsWith(".xls")){
                br = new BufferedReader(new InputStreamReader(dir1));
                AbstractAberatiionReader lr = (AbstractAberatiionReader)
                    clazz.getConstructor(new Class[]{long.class}).newInstance(new Object[]{Long.MAX_VALUE});
                lr.initialise(br,  chromosome, region, 4,1,nxt.getName(), indiv);
                this.addAll(lr);
               // br.close();
            }
            nxt = dir1.getNextEntry();
        }
        this.sort();
       br.close();*/
    }
        
    
}
