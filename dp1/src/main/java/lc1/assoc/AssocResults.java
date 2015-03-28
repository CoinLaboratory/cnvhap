package lc1.assoc;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import javax.swing.table.TableModel;

public interface AssocResults {

    Set<String> getPhenoTypes(String type, Collection<String> allowed);

    List<Integer> loc();

    String[] getPheno(String pheno, String output, int i);
    public int[] getHeaderIndices(String pheno, String output);
    public String[] getPheno(int[] ind, int i);
    TableModel[] getRS(int parseInt, String pheno) throws Exception;

    String chrom(int i);

	String getName();

	

}
