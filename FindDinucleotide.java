import java.io.*;

public class FindDinucleotide{
    
    private Seq seq;
  
    public static void main(String[] args){
	new FindDinucleotide(args[0],args[1]);
    }
  
    public FindDinucleotide(String file, String header){
	this.seq = new FastaReader().parseFasta(file).get(header);
	this.process();
    }
    
    public void process(){
    
	int[] ApN = this.seq.countDinucleotide('A');
	int[] CpN = this.seq.countDinucleotide('C');
	int[] GpN = this.seq.countDinucleotide('G');
	int[] TpN = this.seq.countDinucleotide('T');
	
	printCounts('A', ApN);
	printCounts('C', CpN);
	printCounts('G', GpN);
	printCounts('T',TpN);
    
    }

    public void printCounts(char firstBase, int[] counts){
	System.out.println("FROM:\t" + firstBase + " ->");
	for(int i=0; i<counts.length;i++){
	    System.out.println(this.getBase(i) + " :\t" + counts[i]);
	}
	System.out.println();
    }

    public char getBase(int i){
	if(i==0)
	    return 'A';
	else if(i==1)
	    return 'C';
	else if(i==2)
	    return 'G';
	else if(i==3)
	    return 'T';
	return 'X';
    }
}
