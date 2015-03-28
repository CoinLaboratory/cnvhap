package lc1.util;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Executor {
static List<ExecutorService	> l = new ArrayList<ExecutorService>();
public static ExecutorService getEs(Class clazz, int numT){
	   ExecutorService es =  Executors.newFixedThreadPool(numT);
	   l.add(es);
	   return es;
}

public static ExecutorService getCachedEs(Class clazz){
	   ExecutorService es = Executors.newCachedThreadPool();
	   l.add(es);
	   return es;
}
public static void shutdown(){
	for(int i=0; i<l.size(); i++){
		l.get(i).shutdown();
	}
}

}
