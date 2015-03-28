package lc1.possel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import lc1.CGH.AberationFinder;
import lc1.dp.data.collection.DataCollection;
import lc1.dp.data.collection.SimpleDataCollection;
import lc1.dp.data.representation.Emiss;
import lc1.dp.states.EmissionState;
import lc1.util.Constants;
import edu.mit.wi.haploview.HaploView;

public class PlotLD {
	 public static int value = 4;


	public static void main(String[] args){
	        try{
	        	   File user = new File(System.getProperty("user.dir"));
	        	 
	        	SimpleDataCollection sdt = 
	                SimpleDataCollection.readFastPhaseOutput(AberationFinder.getBufferedReader(user, "phased_states.txt"),Emiss.class, 
	                		
	                		Emiss.getStateEmissionStateSpace(new int[] {4,4}));
	            sdt.removeKeyIfStartsWith("NA");
	        //    sdt.calculateMaf(true);
	            int len = sdt.length();
	           
	            File snpFile  = new File(user, "snp.txt");
	            if(snpFile.exists()){
	                sdt.snpid= new ArrayList<String>();
	                sdt.readPosInfo(snpFile, new int[] {1,0}, false, new List[] {sdt.snpid, sdt.loc}, new Class[] {String.class, Integer.class});
	            }
	            int[] test = new int[] {1};
	            for(int i=0; i<test.length; i++){
	            	value = test[i];
	            writeEHH(sdt, new File(user, "sweep"),Constants.chrom0()+"", null, "test");
	            }
	                                                                                               
	        }catch(Exception exc){
	        	exc.printStackTrace();
	        }
}
	 
	 
	 public static int writeEHH(
		        
		        DataCollection sdt2, File f, String chr, Set<String> indiv, String name) throws Exception{
		        if(!f.exists()) {
		            f.mkdir();
		        }
		      
		        EmissionState maf1 = null;//sdt2.calculateMaf1();
		       // sdt2.maf = maf1;
		        Collection<Integer> toD = maf1.getConstantPos();
		     // toD.addAll(maf1.getMultPos());
		        if(toD.size()>0) {
		            Logger.global.info
		           // throw new RuntimeException
		            ("multi deletions "+chr+" "+name);
		        }
		        File emFile  = new File(f, name+".emphase");
		        File snpFile = new File(f, name+".snp");
		        int max = sdt2.loc.get(sdt2.loc.size()-1);
		        int min = sdt2.loc.get(0);
		        //int lenR = (int) Math.round( (double)(max - sdt2.loc.get(mid[0]))/1000.0);
		        //int lenL = (int) Math.round( (double)(sdt2.loc.get(mid[1])-min)/1000.0);
		        File pngFile = new File(f, name+".emphase.LD.PNG");//sdt2.loc.get(0)+""+sdt2.loc.get(sdt2.loc.size()-1));
		        //File newPng = new File(f, name+
		          //      "-"+lenL+"-"+lenR+"kb"+
		            //    ".LD.PNG");//
		        
		        File track = new File(f, name+".track");
		        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(track)));
		       // pw.print(mid[0]+1);
		      //  for(int i=mid[0]+2; i<=mid[1]+1; i++){
		       //     pw.print(" "+i);
		       // }
		        pw.println();
		        pw.close();
		       // sdt2.snpid.set(mid[0], "deletion");
		        //sdt2.snpid.set(mid[1], "deletion");
		       //for(int i=mid[0]; i<mid[1]; i++){
		       //    toD.add(i);
		       //}
		     //   if(mid[0] < mid[1]){
		      //         toDrop.add(mid[0]);
		     //        mid[1]--;
		     //     }
		       String blocktype = "GAB";
		       if(true){ 
		    	   sdt2.writeDickFormat(emFile, false, indiv, toD);
		           sdt2.writeSNPFile(snpFile, chr, false, toD);
		      // Logger.global.info("hwe "+f+" "+chr+"_"+mid[0]+"_"+mid[1]);
		       HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
		                "-info", snpFile.getAbsolutePath(),
		    //            "-dprime",
		               // "-blocks", track.getAbsolutePath(),
		             //   "-blockoutput", blocktype,
		           //    "-blockCutHighCI", "1.0",
		            //   "-blockCutLowCI","1.0",
		              //  "-ldcolorscheme", "RSQ",
		            //  "-hwcutoff", "0",
		     // "-png" ,"-n"
		               });
		        //File newPng = new File(f, name+
		        //      "-"+lenL+"-"+lenR+"kb"+
		          //    ".LD.PNG");//
		     /* HaploView.main(new String[] {"-haps", emFile.getAbsolutePath(), 
		                "-info", snpFile.getAbsolutePath(),
		                "-dprime",
		               // "-blocks", track.getAbsolutePath(),
		                "-blockoutput", blocktype,
		                //"-blockCutHighCI", "1.0",
		                //"-blockCutLowCI","1.0"
		                "-ldcolorscheme", "RSQ",
		             //   "-hwcutoff", "0.1"
		             "-png" ,"-n"
		               });*/
		       }
		        
		        
		        sdt2.writeDickFormat(emFile, false, indiv, null);
		        sdt2.writeSNPFile(snpFile, chr, true, null);
		      //  Thread.sleep(1000);
		        String blockExt = blocktype.equals("GAB") ? ".GABRIELblocks" :".4GAMblocks";
		        File blockFile = new File(f, name+".emphase"+blockExt);
		        int cnt =0;
		        while(!blockFile.exists() && cnt < 3){
		            Thread.sleep(1000);
		            blockFile = new File(f, name+".emphase"+blockExt);
		            cnt++;
		        }
		        if(blockFile.exists()){
		        BufferedReader br = new BufferedReader(new FileReader(blockFile));
		        String st = "";
		       /* while((st = br.readLine())!=null){
		            if(st.startsWith("BLOCK")){
		                String[] str = st.split(":")[1].trim().split("\\s+");
		                int[] res = new int[str.length];
		                boolean contains = false;
		                for(int jk=0; jk<res.length; jk++){
		                    res[jk] = Integer.parseInt(str[jk]);
		                    if(res[jk]==mid[0] || res[jk]==mid[1]) contains = true;
		                }
		                if(contains) return sdt2.loc.get(res[res.length-1]) - sdt2.loc.get(res[0]);
		            }
		        }*/
		        }
		        return 0;
		       // if(true) Thread.currentThread().wait(1000000);
		       // String[] cmd = new String[] {"mv", f.getName()+"/"+pngFile.getName(), f.getName()+"/"+newPng.getName()};
		       // StringWriter err = new StringWriter();
		       // ProcessTools.exec(cmd, null, err, err);
		      //  System.err.println(err.toString());
		         // TODO Auto-generated method stub
		         
		        // TODO Auto-generated method stub
		    }
}
