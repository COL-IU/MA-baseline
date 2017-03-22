import java.io.*;
import java.util.*;

public class CountRuns{

    public static void main(String[] args){
	Seq ecoliSeq = new FastaReader().parseFasta(args[0]).get("ecoli");
	ArrayList<Gene> geneList = new PttParser().parse(args[1], "ecoli");
	ArrayList<RunLengthGeneCounter> runLengthGeneCounterList = ecoliSeq.getRuns(geneList);
	for(int i=0; i<runLengthGeneCounterList.size();i++){
	    runLengthGeneCounterList.get(i).print();
	}
	
    }

}
