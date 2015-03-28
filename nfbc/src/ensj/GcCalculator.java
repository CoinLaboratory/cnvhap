package ensj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Sequence;
import org.ensembl.driver.CoreDriver;
import org.ensembl.driver.SequenceAdaptor;
import org.ensembl.registry.Registry;
import org.ensembl.variation.driver.VariationDriver;

import conversion.OptionBuild;

public class GcCalculator {

	private CoreDriver coreDriver;
	private VariationDriver variationDriver;

	// String build1 = conversion.OptionBuild.build;
	// if(!build.startsWith(build1.split("ild")[1].substring(0,2))) {
	// throw new RuntimeException(build+" !! "+build1);
	// }

	public void getChromosomeCoordinates() {

	}

	public void calculateGCcontent(String dir, String name, int windowSize)
			throws Exception {

		Registry dr = Registry.createDefaultRegistry();
		coreDriver = dr.getGroup("human").getCoreDriver();
		variationDriver = dr.getGroup("human").getVariationDriver();
		String[] database = coreDriver.getConfiguration().getProperty(
				"database").split("_");
		String build = database[database.length - 1];

		// File in = new File(dir, name+".txt");
		// BufferedReader br = conversion.Utils.getBufferedReader(in);
		File outf = new File(dir, name + ".gc.txt");
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(
				outf, true)));
		
		BufferedReader br1 = conversion.Utils.getBufferedReader(outf);
		String st1 = br1 == null ? null : br1.readLine();

//		CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
		int[] chrIncluded = new int[] { 1 };
		 int[] chrLengths = new int[]{247199719};
//		int[] chrLengths = new int[] { 10000 };

		for (int chr : chrIncluded)
			for (int pos = (0 + windowSize / 2); pos <= chrLengths[chr-1]
					- windowSize / 2; pos += windowSize)

			{
				try {

					// String chr = str[0].substring(3);
					// int pos = Integer.parseInt(str[1]);
					Location chromosomeLoc = new Location("chromosome:" + chr
							+ ":" + (pos - Math.round(windowSize / 2)) + "-"
							+ (pos + Math.round(windowSize / 2)));
//					System.out.println("chromosome:" + chr + ":" + (pos - Math.round(windowSize / 2)) + "-"
//							+ (pos + Math.round(windowSize / 2)));
					SequenceAdaptor se = coreDriver.getSequenceAdaptor();
					Sequence seq = se.fetch(chromosomeLoc);
					String str1 = seq.getString();

					se.closeAllConnections();
					int gc = 0;
					int at = 0;
					for (int i = 0; i < str1.length(); i++) {
						char ch = str1.charAt(i);
						if (ch == 'G' || ch == 'C')
							gc++;
						else if (ch == 'A' || ch == 'T')
							at++;
					}
					double gcc = ((double) gc) / ((double) gc + (double) at);
					out.println(pos+"\t" + String.format("%5.3g", gcc));
					out.flush();
				} catch (Exception exc) {
					exc.printStackTrace();
					out.flush();
					System.out.println("chromosome:" + chr
							+ ":" + (pos - Math.round(windowSize / 2)) + "-"
							+ (pos + Math.round(windowSize / 2)));
					out.println(pos+"\t" + Double.NaN);
				}
			}

		out.close();
		coreDriver.closeAllConnections();

	}

	public static void main(String[] args) throws Exception {
		GcCalculator gcc = new GcCalculator();
		gcc.calculateGCcontent(
				"C:\\Users\\Evan\\Desktop\\chr1 GC content", "chr1",1000);
	}

}
