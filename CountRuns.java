import java.io.*;
import java.util.*;

public class CountRuns{

// args[0] -> fasta file
// args[1] -> ptt file
// args[2] -> fasta header
    public static void main(String[] args){
	Seq seq = new FastaReader().parseFasta(args[0]).get(args[2]); 
	ArrayList<Gene> geneList = new PttParser().parse(args[1], args[2]);
	ArrayList<RunLengthGeneCounter> runLengthGeneCounterList = seq.getRuns(geneList);
	for(int i=0; i<runLengthGeneCounterList.size();i++){
	    runLengthGeneCounterList.get(i).print();
	}
	
    }

}
