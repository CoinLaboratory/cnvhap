package lc1.dp.data.collection;


import java.util.List;
import java.util.logging.Logger;

import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;

public class ListTableModel extends AbstractTableModel implements TableModelListener{
    List<Object[] > tab;
    final boolean invert;
    int no_cols;
    int no_rows;
    public ListTableModel(List<Object[]> res, boolean invert) {
       this.tab = res;
       this.invert = invert;
       no_rows = res.size();
       no_cols =0;
       for(int i=0; i<res.size(); i++){
    	   int a = res.get(i).length;
    	   if(a> no_cols){
    		   no_cols=  a;
    	   }
       }
       this.addTableModelListener(this);
    }
  /*  @Override
    public void fireTableCellUpdated(int row, int column){
    	v
    }*/
    public int getRowCount() {
       return invert ? no_cols : no_rows;
    }

    public int getColumnCount() {
       return invert ? no_rows : no_cols;//tab.get(0).length;
    }

    public Object getValueAt(int rowIndex, int columnIndex) {
    	if(invert){
    		Object[] obj = tab.get( columnIndex );
    		return rowIndex >=obj.length ? "" : obj[rowIndex];
    	}
    	else{
    		Object[] obj = tab.get( rowIndex );
    		return columnIndex >=obj.length ? "" : obj[columnIndex];
    	}
    
    	
     
    }
@Override
   public boolean isCellEditable(int row, int col){
	   return true;
   }

public void tableChanged(TableModelEvent e) {
	//tab.get(e.getFirstRow())[e.getColumn()] = e.
	Logger.global.info("h");
}

}
