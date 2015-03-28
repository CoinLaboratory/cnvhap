package conversion;

import java.io.File;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class AccessZip {
public static void main(String[] args){
   
        String[] string = args[0].split(":");
        for(int j=0; j<string.length; j++){
            try{
            File f = new File(string[j]);
            ZipFile zf = new ZipFile(f);
            System.err.println(f.exists()+" "+string[j]);
            if(!f.exists()) throw new RuntimeException("does not exist "+f);
            int index = which( Compressor.getIndiv(zf, "Samples"), args[2]);
    
            Set<String> entries = getEntries(zf, new HashSet<String>(Arrays.asList(new String[] {"SNPS", "Samples", "Name"})));
            Set<String> names = getIndiv(zf);
               Set<String> entriesO = new HashSet<String>(entries);
               System.err.println(entriesO.size()+" "+names.size());
               entriesO.removeAll(names);
               names.removeAll(entries);
               if(entriesO.size()>0){
                   System.err.println("extra entries "+entriesO);
               }
               if(names.size()>0){
                   System.err.println("extra names "+names);
               }
             
               //for(Enumeration en = zf.entries(); en.hasMoreElements();){
               //    System.err.println(en.nextElement()+",");
              // }
               System.err.println(Compressor.getIndiv(zf, args[1]+"").get(index));
        }catch(Exception exc){
            exc.printStackTrace();
        }
        }
    
}
private static int which(List<String> indiv, String string) {
   for(int i=0; i<indiv.size(); i++){
       if(indiv.get(i).startsWith(string)) return i;
       
   }
   return -1;
}
public static  Set<String> getEntries(ZipFile zf, Set<String>excl){
    Set<String> res = new HashSet<String>();
    for(Enumeration en = zf.entries(); en.hasMoreElements();){
        res.add(((ZipEntry)en.nextElement()).getName());
    }
    res.removeAll(excl);
    return res;
}
public  static Set<String> getIndiv(ZipFile zf) throws Exception{
    Set<String> l  = new HashSet<String>();
    for(Iterator<String> it = Compressor.getIndiv(zf, "SNPS").iterator(); it.hasNext();){
        String st = it.next().split("\t")[0];
      boolean a =   l.add(st);
      if(!a) throw new RuntimeException("!! "+st);
    }
   
    return l;
    
}
}
