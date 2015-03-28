import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import lc1.dp.data.collection.DataCollection;
import lc1.ensj.GenesForRegion;


public class SplitInputFile {
	public static String splitM = "genome";
	public static int step = 300000;
	public static int overlap = 100;
	public static int halfoverlap = (int) Math.round((double) overlap/2.0);
	 static GenesForRegion grf = new GenesForRegion();
	public static void main(String[] args){
		try{
			String d = System.getProperty("user.dir");
	        File dir1 = new File(d);
			SplitInputFile sif = new SplitInputFile(dir1, args);
			if(splitM.equals( "local")) sif.splitLocal();
			else sif.splitGlobal();
			sif.pw.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	private void splitGlobal() throws Exception{
		String currentChr=null;
		String currentString;
		SortedSet<Integer> positions = new TreeSet<Integer>();
		if(buildF.length>1) throw new RuntimeException("!!");
		BufferedReader br = new BufferedReader(new FileReader(buildF[0]));
		currentString = br.readLine();
		while(currentString!=null){
			if(currentChr!=null) writePositions(currentChr, positions);
			String[] str = currentString.split("\\s+");
			currentChr = str[0];
			
			positions.clear();
		
			inner: while(str[0].equals(currentChr)){
				positions.add(Integer.parseInt(str[1]));
				currentString = br.readLine();
				if(currentString==null) break inner;
				str = currentString.split("\\s+");
			}
		}
		writePositions(currentChr, positions);
		pw.close();
		
	}
	
	Map<String, Integer> centromere = new HashMap<String, Integer>();
	
	public int getCentromere(String st) throws Exception{
		Integer res = centromere.get(st);
		if(res==null){
			centromere.put(st, res = grf.getCentromere(st));
		}
		return res.intValue();
	}
	private void writePositions(String chr, SortedSet<Integer> positions) throws Exception {
		if(chr.equals("M") || chr.equals("XY")) return;
		int centr = getCentromere(chr.substring(3));
		int last = positions.last();
		int mv = Integer.MAX_VALUE;
		System.err.println("centromere is "+centr+" chr:"+chr);
		List<Integer>[] pos = new List[] {new ArrayList<Integer>(positions.headSet(centr)), 
							new ArrayList<Integer>(positions.tailSet(centr))};
		for(int j=0; j<pos.length; j++){
			double num = Math.ceil((double)pos[j].size()/(double)step); //number of groups
			int step1 = (int) Math.ceil((double) pos[j].size()/num);
			int end_core = j==0 ? 0 : centr;
			for(int i=0; i<pos[j].size(); i+=step1){
				int start_core = end_core;
				int start =i==0 && j==0 ? 0 :(i==0 && j==1 ? centr+1 :   pos[j].get(i));
				int end;
				if(j==0 && i+step1>=pos[j].size()){
					end = centr;
					end_core = centr;
				}
				else if(j==1 && i+step1>=pos[j].size()){
					int ist = i+step1;
					int sze = pos[j].size();
					int len = sze -i;
					end  = pos[j].get(pos[j].size()-1)+1;
					end_core = end;
				}
				else{
					end = pos[j].get(Math.min(pos[j].size()-1, i+step1+overlap))+1;
					end_core = pos[j].get(Math.min(pos[j].size()-1, i+step1+halfoverlap))+1;
				}
				pw.println("--mid\t"+chr.substring(3)+":"+start+":"+end+"\t--core\t"+start_core+":"+end_core);
			}
		}
		pw.flush();
		
	}
	final File[] buildF;
	public  void splitLocal() throws Exception{
		SplitLocal sif = new SplitLocal();
		for(Iterator<String> it = sif.loc.keySet().iterator(); it.hasNext();){
			sif.split(it.next());
			
		}
	}
	
	class SplitLocal{
		double thresh = 3000000;
		Map<String, SortedSet<Integer>> loc = new TreeMap<String, SortedSet<Integer>>();
		
		SplitLocal() throws Exception{
			for(int i=0; i<buildF.length; i++){
				BufferedReader br = DataCollection.getBufferedReader(buildF[i]);
				String st = "";
				while((st =br.readLine())!=null){
					String[] str = st.split("\t");
					SortedSet<Integer> loc1 = loc.get(str[0]);
					if(loc1==null){
						loc.put(str[0], loc1 = new TreeSet<Integer>());
					}
					loc1.add(Integer.parseInt(str[1]));
				}
				br.close();
			}
		}
		
		void split(String chr){
			//end = new StringBuffer();
			Iterator<Integer> it = loc.get(chr).iterator();
			int last=-1;
			int start_=0;
		//	start.add(start_);
			int end_;
			boolean first = true;
			inner: while( it.hasNext()){
				Integer i = it.next();
				if(first){
					start_ = i-1;
					first = false;
				}
				if(last<0 || (i-last)< thresh){
					last = i;
				}
				else{
					end_ = last+1;
					//end.append((last+1)+":");
					last =i;
					pw.println("--mid  "+chr+":"+start_+":"+end_);
					start_ = last-1;
					//start.add(last-1);
				}
			}
		//	end.append(last+1);
			
			//System.err.println(chr+":"+end.toString());
		}
		
	}
	
	
	
	SplitInputFile(File dir, String[] buildF) throws Exception{
		String st = "";
		String outFile = "split"+step+"_"+overlap+".txt";
		this.buildF = new File[buildF.length];
		for(int i=0; i<buildF.length; i++){
			this.buildF[i] = new File(dir, buildF[i]);
		}
		pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(dir, outFile))));
		
	}
	
	
	final PrintWriter pw;
	
	
	
	
   
}
/*
  List<Integer>loc = new ArrayList<Integer>();
    for(int ik=0; ik<files.length; ik++){
        loc.addAll(readPositions(chr, new File(baseDir+"/"+files[ik]+"/"+build)));
    }
    Collections.sort(loc);
//  System.err.println(step);
    inner: for(int j=0; ; j++){
        int k = step*j;
        if(k>=loc.size()) break inner;
        int start = loc.get(k);
        int endIndex = k+step+overlap;
        int endKb =  endIndex>=loc.size() ? loc.get(loc.size()-1)+1000 : loc.get(endIndex);
        String identString = st+"-"+Math.abs(opt_string.hashCode());
        File logFile = CreateSubmissions.print1(user, outputFile, identString, chr, start, endKb+1, cols, files, build,dataP,fixedP);
      print(logFile,chr,start-1, endKb+1,0, logIn, false
             );
       cnt++;
    }
       static Set<String> done ;
            static String dF ;
            done = new HashSet<String>();//getDone(doneF);

 private Collection<? extends Integer> readPositions(String chr, File file) throws Exception{
       if(!file.exists() && !(new File(file.getParentFile(), file.getName()+".gz").exists())){
    	   return new HashSet<Integer>();
       }
       BufferedReader br = DataCollection.getBufferedReader(file);
       System.err.println(file.getName());
       String st  = "";
       String ch = "chr"+chr;
       List<Integer> res = new ArrayList<Integer>();
       while((st = br.readLine())!=null){
           String[] str = st.split("\\s");
           if(str[0].equals(ch)){
               res.add(Integer.parseInt(str[1]));
           }
       }
       return res;
    } 
    
     public static Set<String> getDone(File f){
        Set<String> done = new HashSet<String>();
        File[] fl = f.listFiles(new FileFilter(){

            public boolean accept(File pathname) {
                return pathname.isDirectory() || pathname.getName().endsWith(".tar.gz") || pathname.getName().endsWith(".tar");
            }
            
        });
        for(int i=0; i<fl.length; i++){
            String name = fl[i].getName();
            try{
        //   BufferedReader br =  AberationFinder.getBufferedReader(fl[i], "phased.txt");
        //   br.readLine();
            }catch(Exception exc){
                continue;
            }
            int gzIndex = name.indexOf(".gz");
            if(gzIndex>=0) name = name.substring(0,gzIndex);
            gzIndex = name.indexOf(".tar");
            if(gzIndex>=0) name = name.substring(0,gzIndex);
            done.add(name);
        }
        return done;
    }
*/