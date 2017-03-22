import java.util.*;

/**
 * Seq is the base class for holding dna sequence datas.
 * 
 * @author Heewook Lee
 * @date Jan 25 2007
 * @updataed July 03, 2007
 * @ver 1.1
 */
public class Seq{
    
    private StringBuffer seqBuffer; // global variable for holding sequence data in StringBuffer Obj
    private String seqName; // golbal variable for holding sequence name in Sting Obj
    
    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     * @param seq sequence in String obj
     * 
     */
    public Seq(String fastaName, String seq){
        seqName = fastaName;
        seqBuffer = new StringBuffer(seq);
    }
    
    
    public int[] countDinucleotide(char firstBase){
	int[] count = new int[4];//ACGT
	for(int i=0;i<(seqBuffer.length()-1);i++){
	    if(seqBuffer.charAt(i)==firstBase){
		char nb = seqBuffer.charAt(i+1);
		if(nb=='A')
		    count[0]++;
		else if(nb == 'C')
		    count[1]++;
		else if(nb =='G')
		    count[2]++;
		else if(nb =='T')
		    count[3]++;
	    }
	}
	return count;
    }
    
    /*
     * A: 0   C: 1   G: 2   T: 3
     *
     */
    public int[] countAllDinucleotide(){
	int[] dimerCount = new int[16];
	for(int i=0; i<(seqBuffer.length()-1); i++){
	    dimerCount[this.baseToNumber(seqBuffer.charAt(i)) * 4 + this.baseToNumber(seqBuffer.charAt(i+1))]++;
	}
	dimerCount[this.baseToNumber(seqBuffer.charAt(seqBuffer.length()-1)) * 4 + this.baseToNumber(seqBuffer.charAt(0))]++;
	return dimerCount;
    }
    
    public int[] countAllTrinucleotide(){
	int[] trimerCount = new int[64];
	for(int i=0; i<(seqBuffer.length()-2); i++){
	    trimerCount[this.baseToNumber(seqBuffer.charAt(i)) * 16 
			+ this.baseToNumber(seqBuffer.charAt(i+1)) * 4 
			+ this.baseToNumber(seqBuffer.charAt(i+2)) ]++;
	}
	trimerCount[this.baseToNumber(seqBuffer.charAt(seqBuffer.length()-2)) * 16 
		    + this.baseToNumber(seqBuffer.charAt(seqBuffer.length()-1)) * 4
		    + this.baseToNumber(seqBuffer.charAt(0)) ]++;
	trimerCount[this.baseToNumber(seqBuffer.charAt(seqBuffer.length()-1)) * 16 
		    + this.baseToNumber(seqBuffer.charAt(0)) * 4
		    + this.baseToNumber(seqBuffer.charAt(1)) ]++;
	return trimerCount;
    }
    
    public int baseToNumber(char base){
	if(base == 'A' || base == 'a')
	    return 0;
	else if(base == 'C' || base == 'c')
	    return 1;
	else if(base == 'G' || base == 'g')
	    return 2;
	else if(base == 'T' || base == 't')
	    return 3;
	else
	    return -1;
    }

    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     * @param seq sequence in String obj
     * 
     */
    public Seq(String fastaName, StringBuffer seq){
        seqName = fastaName;
        seqBuffer = seq;
    }
    
    /*
     * Constructor
     *
     * @param fastaName name of the fasta seq
     *
     */
    public Seq(String fastaName){
	seqName = fastaName;
	seqBuffer = new StringBuffer();
    }

    public Seq(){

    }


    public Gene toGene(){
	return new Gene(seqName, 1, seqBuffer.length(), true, seqName);
    }

    /*
     * getSeqName() method simple returns the seqName
     *
     * @return seqName String obj seqName is returned
     */
    public String getSeqName(){
        return seqName;
    }
    
    /*
     * getSeq() method simple returns the sequence buffer
     *
     * @return seqBuffer StringBuffer obj sequence is returned
     */
    public StringBuffer getSeq(){
        return new StringBuffer(this.seqBuffer.toString());
    }

    public int getSeqLength(){
	return this.seqBuffer.length();
    }
    
    public Seq deepCopy(){
	Seq newer = new Seq();
	newer.seqBuffer = this.getSeq();
	newer.seqName = this.seqName;
	return newer;
    }

    public Seq getRevSeq(){
	return new Seq(this.seqName, new StringBuffer(this.revSeqInternal().toString()));
    }
    /*
     * formatSeq() method returns the sequence in fasta format without the seq name tag 
     * each line containts 50 neucleotide bases.
     *   
     * @return outBuffer StringBuffer obj sequnce is returned with newline character every 50 bases
     */
    public StringBuffer formatSeq(){
	StringBuffer outBuffer = new StringBuffer();
	int lowInd = 0;
	for(int i=0; i<seqBuffer.length();i++){
            if((i % 50) == 49){
                outBuffer.append(seqBuffer.substring(lowInd,i+1));
                outBuffer.append("\n");
                lowInd += 50;
            }
        }
	if( (seqBuffer.length()%50) != 0){
            outBuffer.append(seqBuffer.substring(lowInd));
            outBuffer.append("\n");
        }
        return outBuffer;
    }
    
    public ArrayList<IntegerCounter> countRuns(){
	//System.err.println("TOTAL LEN : " + this.seqBuffer.length());
	ArrayList<IntegerCounter> runCountArr = new ArrayList<IntegerCounter>();
	int curRunLen = 0;
	char curRunChar = '0';
	for(int i=0;i<this.seqBuffer.length();i++){
	    char curChar = Character.toUpperCase(seqBuffer.charAt(i));
	    if(curChar == curRunChar)//if it's a currentRun
		curRunLen++;
	    else{/* need to update now*/
		if(curRunLen > 0){
		    if(runCountArr.size() < curRunLen)
			updateRunArr(runCountArr,curRunLen);
		    runCountArr.get(curRunLen-1).add1(curRunChar);
		}
		//if(curRunLen == 10)
		//    System.err.println(i + "\t" + curRunChar);
		curRunLen = 1;
		curRunChar = curChar;
	    }
	}
	if(curRunLen > 0 )
	    runCountArr.get(curRunLen-1).add1(curRunChar);
	
	return runCountArr;
    }
    
    public void updateRunLengthArr(ArrayList<RunLengthGeneCounter> arr, int runLen ){
	int l = runLen + 1;
	while(arr.size() < runLen){
	    arr.add(new RunLengthGeneCounter(l));
	    l++;
	}
    }
    
    private void initRunLengthArr(ArrayList<RunLengthGeneCounter> arr){
	for(int i=0; i<10; i++){
	    arr.add(new RunLengthGeneCounter(i+1));
	}
    }

    public ArrayList<RunLengthGeneCounter> getRuns(ArrayList<Gene> geneList){
	ArrayList<RunLengthGeneCounter> runLengthCountArr = new ArrayList<RunLengthGeneCounter>();
	this.initRunLengthArr(runLengthCountArr);
	int curRunLen = 0;
	char curRunChar = '0';
	for(int i=0;i<this.seqBuffer.length();i++){
	    char curChar = Character.toUpperCase(seqBuffer.charAt(i));
	    if(curChar == curRunChar)//if it's a currentRun
		curRunLen++;
	    else{/* need to update now*/
		
		if(runLengthCountArr.size() < curRunLen)
		    this.updateRunLengthArr(runLengthCountArr, curRunLen);
		if(curRunLen > 3)
		    runLengthCountArr.get(curRunLen-1).addPosition( ( i - curRunLen + 1), curRunChar, geneList);
		
		curRunLen = 1;
		curRunChar = curChar;
	    }
	}
	
	if(runLengthCountArr.size() < curRunLen)
	    this.updateRunLengthArr(runLengthCountArr, curRunLen);
	if(curRunLen > 1 )
	    runLengthCountArr.get(curRunLen-1).addPosition( ( this.seqBuffer.length() - curRunLen + 1), this.seqBuffer.charAt(this.seqBuffer.length()-1) , geneList);
	
	return runLengthCountArr;
    }




    /* ACGT, row<ACGT>, col<ACGT>
     *
     * This method outputs a 3 dimensional array: basically 4 2-dimensional array
     * where each array is 4 by 4 <ACGT> vs <ACGT> where row is the 5' neighbor base and col is the 3' neighbor base
     * 
     */
    public int[][][] getNeighboringBaseCountMatrix(){
	int[][][] nbcMatrix = new int[4][4][4];
	
	char startBase = Character.toUpperCase(this.seqBuffer.charAt(0));
	char endBase = Character.toUpperCase(this.seqBuffer.charAt(this.seqBuffer.length()-1));
	
	char prevBase = endBase;
	char curBase = startBase;//Character.toUpperCase(this.seqBuffer.charAt(0));
	char nextBase = '0';
	for(int i=1;i<this.seqBuffer.length();i++){
	    nextBase = Character.toUpperCase(this.seqBuffer.charAt(i));
	    nbcMatrix[this.baseToIndex(curBase)][this.baseToIndex(prevBase)][this.baseToIndex(nextBase)]++;
	    
	    prevBase = curBase;
	    curBase = nextBase;
	}
	nextBase = startBase;
	nbcMatrix[this.baseToIndex(curBase)][this.baseToIndex(prevBase)][this.baseToIndex(nextBase)]++;
	
	return nbcMatrix;
    }

    private int baseToIndex(char b){
	if(b == 'A')
	    return 0;
	else if(b == 'C')
	    return 1;
	else if(b == 'G')
	    return 2;
	else if(b == 'T')
	    return 3;
	else
	    return -1;
    }
    

    private void updateRunArr(ArrayList<IntegerCounter> arr, int curRunLen){
	while(arr.size() < curRunLen){
	    arr.add(new IntegerCounter());
	}
    }

    public ArrayList<Integer>[] getMotifLoci(String motif, boolean bothStrands){
	int fc = 0;
	int rc = 0;
	ArrayList<Integer>[] lociArr = null;
	if(bothStrands){
	    lociArr = new ArrayList[2];
	    lociArr[1] = new ArrayList<Integer>(); /* for rev strand */
	}else
	    lociArr = new ArrayList[1];
	lociArr[0] = new ArrayList<Integer>(); /* for fwd strand */

	int i = 0;
	while(i<this.seqBuffer.length()){
	    int tmp =this.seqBuffer.indexOf(motif, i);
	    if(tmp < 0)
		i = this.seqBuffer.length();
	    else{
		lociArr[0].add(new Integer(tmp+1));
		fc++;
		i = tmp + 1;
	    }
	}
	if(bothStrands){
	    StringBuffer rev = revSeqInternal();
	    i = 0;
	    while(i<rev.length()){
		int tmp = rev.indexOf(motif, i);
		if(tmp < 0)
		    i = this.seqBuffer.length();
		else{
		    lociArr[1].add(new Integer(rev.length() - tmp));
		    rc++;
		    i = tmp + 1;
		}
	    }
	}
	
	System.err.println("Fwd:\t" + fc);
	System.err.println("Rev:\t" + rc);

	return lociArr;
    }


    //thi is to support extended nucleotide bases (IUPAC codes)
    public ArrayList<Integer>[] getMotifLociExtendedAlphabets(String motif, boolean bothStrands){
	//System.out.println("SEQ:\t" + this.seqBuffer + "\tMOTIF:\t" + motif);
	int fc = 0;
	int rc = 0;
	ArrayList<Integer>[] lociArr = null;
	if(bothStrands){
	    lociArr = new ArrayList[2];
	    lociArr[1] = new ArrayList<Integer>(); /* for rev strand */
	}else
	    lociArr = new ArrayList[1];
	lociArr[0] = new ArrayList<Integer>(); /* for fwd strand */

	int i = 0;
	while(i<this.seqBuffer.length()){
	    //System.err.println("fromIndex:\t" + i);
	    int tmp =this.indexOfNextMotif(motif, i);
	    if(tmp < 0){
		//System.err.println("HIT\t" + i);
		i = this.seqBuffer.length();
	    }else{
		lociArr[0].add(new Integer(tmp+1));
		fc++;
		i = tmp + 1;
	    }
	}
	if(bothStrands){
	    //StringBuffer rev = revSeqInternal();
	    Seq rev = this.getRevSeq();
	    i = 0;
	    while(i<rev.getSeqLength()){
		int tmp = rev.indexOfNextMotif(motif, i);
		if(tmp < 0){
		    //System.err.println("REVHIT\t" + i);
		    i = this.seqBuffer.length();
		}else{
		    lociArr[1].add(new Integer(rev.getSeqLength() - tmp));
		    rc++;
		    i = tmp + 1;
		}
	    }
	}
	
	System.err.println("Fwd:\t" + fc);
	System.err.println("Rev:\t" + rc);

	return lociArr;
    }
    
    private int indexOfNextMotif(String motif, int fromIndex){
	for(int i = fromIndex; i<(this.seqBuffer.length()-motif.length()+1); i++){
	    if(isMotifMatchAt(motif, i))
		return i;
	}
	return -1;//if not match is found.
    }

    private boolean isMotifMatchAt(String motif, int index){
	int i;
	int j = 0;
	for(i=index; i<this.seqBuffer.length() && j<motif.length();i++ ){
	    if(!this.matchBaseAt(i,motif.charAt(j)))
		return false;
	    j++;
	}
	//System.err.println("Do we ever match Motif??");
	return true;
    }

    //0-base index 
    private boolean matchBaseAt(int position, char base){
	if(Seq.matchBase(this.seqBuffer.charAt(position), base))
	    return true;
	return false;
    }

    public static boolean matchBase(char refBase, char base){
	char ru = Character.toUpperCase(refBase);
	char bu = Character.toUpperCase(base);
	boolean returnFlag = false;
	if(bu == 'R'){
	    if(ru == 'A' || ru == 'G')
		return true;
	}else if(bu == 'Y'){
	    if(ru == 'C' || ru == 'T')
		return true;
	}else if(bu == 'S'){
	    if(ru == 'G' || ru == 'C')
		return true;
	}else if(bu == 'W'){
	    if(ru == 'A' || ru == 'T')
		return true;
	}else if(bu == 'K'){
	    if(ru == 'G' || ru == 'T')
		return true;
	}else if(bu == 'M'){
	    if(ru == 'A' || ru == 'C')
		return true;
	}else if(bu == 'B'){
	    if(ru == 'C' || ru == 'G' || ru == 'T')
		return true;
	}else if(bu == 'D'){
	    if(ru == 'A' || ru == 'G' || ru == 'T')
		return true;
	}else if(bu == 'H'){
	    if(ru == 'A' || ru == 'C' || ru == 'T')
		return true;
	}else if(bu == 'V'){
	    if(ru == 'A' || ru == 'C' || ru == 'G')
		return true;
	}else if(bu == 'N')
	    return true;
	else if(bu == 'A' || bu == 'C' || bu == 'G' || bu == 'T'){
	    if(ru == bu)
		return true;
	}
	return false;
    }

    public void printACGT(Gene g){
	int[] acgtCounts = this.countACGT(g);
	System.out.println(g.getHeader() + "\t" + acgtCounts[0] + "\t" + acgtCounts[1] + "\t" + acgtCounts[2] + "\t" + acgtCounts[3] + "\t" + acgtCounts[4]);
    }

    public int[] countACGT(Gene g){
	int[] acgtCounts = new int[5];
	for(int i=g.getStart()-1;i<g.getEnd();i++){
	    char tmp = seqBuffer.charAt(i);
	    if(g.fwd()){
		if(tmp=='A' || tmp=='a')
		    acgtCounts[0]++;
		else if(tmp == 'C' || tmp == 'c')
		    acgtCounts[1]++;
		else if(tmp == 'G' || tmp == 'g')
		    acgtCounts[2]++;
		else if(tmp == 'T' || tmp == 't')
		    acgtCounts[3]++;
		else
		    acgtCounts[4]++;
	    }else{//revcomp
		if(tmp=='A' || tmp=='a')
		    acgtCounts[3]++;
		else if(tmp == 'C' || tmp == 'c')
		    acgtCounts[2]++;
		else if(tmp == 'G' || tmp == 'g')
		    acgtCounts[1]++;
		else if(tmp == 'T' || tmp == 't')
		    acgtCounts[0]++;
		else
		    acgtCounts[4]++;
	    }
	}
	//System.out.println(g.getHeader() + "\t" + acgtCounts[0] + "\t" + acgtCounts[1] + "\t" + acgtCounts[2] + "\t" + acgtCounts[3] + "\t" + acgtCounts[4]);
	return acgtCounts;
    }

    public int[] countMotif(String motif, boolean bothStrands){
	int[] counts;
	if(bothStrands)
	    counts = new int[2];
	else
	    counts = new int[1];
	
	for(int i=0; i<counts.length;i++){
	    counts[i] = 0;
	}

	int i=0;
	while(i<this.seqBuffer.length()){
	    int tmp = this.seqBuffer.indexOf(motif, i);
	    if(tmp < 0)
		i = this.seqBuffer.length();
	    else{
		counts[0]++;
		i = tmp + 1;
	    }
	}
	
	if(bothStrands){
	    StringBuffer rev = revSeqInternal();
	    i = 0;
	    while(i<rev.length()){
		int tmp = rev.indexOf(motif, i);
		if(tmp < 0)
		    i = this.seqBuffer.length();
		else{
		    counts[1]++;
		    i = tmp + 1;
		}
	    }
	
	}
	return counts;
    }

    /*
     * appendSeq( String ) method appends the sequence 
     * at the end of the sequence held by Seq obj.
     *
     */
    public void appendSeq(String tempSeq){
        seqBuffer.append(tempSeq);
    }

    public void insertSeq(int offset, String str){
	seqBuffer.insert(offset, str);
    }

    public void insertSeq(int offset, char c){
	seqBuffer.insert(offset, c);
    }

    public void deleteSeqAt(int offset){
	seqBuffer.deleteCharAt(offset);
    }

    public void replaceSeqAt(int offset, char c){
	seqBuffer.replace(offset, offset+1, ""+c);
    }
    
    public void replaceSeqAt(int offset, String str){
	seqBuffer.replace(offset, offset+1, str);
    }
    
    public char seqAt(int index){
	return seqBuffer.charAt(index-1);
    }


    public String getPutationStringForNthBase(int mutTypeIndex, int chosenNthBase, int numLines){
	char[] targets = new char[2];
	if(mutTypeIndex%2 == 0){ //if an indices is even, it starts from AT, otherwise, it starts from GC
	    targets[0] = 'A';
	    targets[1] = 'T';
	}else{
	    targets[0] = 'G';
	    targets[1] = 'C';
	}
	int count = 0;
	for(int i=0; i<this.seqBuffer.length(); i++){
	    char curChar = Character.toUpperCase(seqBuffer.charAt(i));
	    if(curChar == targets[0] || curChar == targets[1]){
		//count++;
		if(count == chosenNthBase){
		    return "ecoli " + i + " " + curChar + " " + makePutationArray(curChar, getMutChar(curChar, mutTypeIndex),numLines) + curChar + "\n";
		    //System.out.println(mutTypeIndex + "\t" + chosenNthBase + "\t" + targets[0]+"/"+targets[1] + "\t" + i);
		    //return i;
		}
		count++;
	    }
	}
	return null;//return -1;
    }    


    /*
     * mutTypeIndex 
     * 0 , 1 : Transition : AT -> GC,  GC -> AT
     * 2, 3, 4, 5 : Transversion : AT -> TA, GC -> TA, AT -> CG, GC -> CG
     *
     */
    public int getIndexForNthBase(int mutTypeIndex, int chosenNthBase, int numLines){
	char[] targets = new char[2];
	if(mutTypeIndex%2 == 0){ //if an indices is even, it starts from AT, otherwise, it starts from GC
	    targets[0] = 'A';
	    targets[1] = 'T';
	}else{
	    targets[0] = 'G';
	    targets[1] = 'C';
	}
	int count = 0;
	for(int i=0; i<this.seqBuffer.length(); i++){
	    char curChar = Character.toUpperCase(seqBuffer.charAt(i));
	    if(curChar == targets[0] || curChar == targets[1]){
		//count++;
		if(count == chosenNthBase){
		    System.out.println("ecoli " + i + " " + curChar + " " + makePutationArray(curChar, getMutChar(curChar, mutTypeIndex),numLines) + curChar );
		    //System.out.println(mutTypeIndex + "\t" + chosenNthBase + "\t" + targets[0]+"/"+targets[1] + "\t" + i);
		    return i;
		}
		count++;
	    }
	}
	return -1;
    }

    public int getIndexForNthBase(int mutTypeIndex, int chosenNthBase){
	return this.getIndexForNthBase(mutTypeIndex, chosenNthBase, 34);
    }


    private char getMutChar(char curChar, int mutTypeIndex){
	if(mutTypeIndex == 0){
	    if(curChar == 'A')
		return 'G';
	    else
		return 'C';
	}else if(mutTypeIndex == 1){
	    if(curChar == 'G')
		return 'A';
	    else
		return 'T';
	}else if(mutTypeIndex == 2){
	    if(curChar == 'A')
		return 'T';
	    else
		return 'A';
	}else if(mutTypeIndex == 3){
	    if(curChar == 'G')
		return 'T';
	    else
		return 'A';
	}else if(mutTypeIndex == 4){
	    if(curChar == 'A')
		return 'C';
	    else
		return 'G';
	}else if(mutTypeIndex == 5){
	    if(curChar == 'G')
		return 'C';
	    else
		return 'G';
	}else
	    return 'X';
	    
    }

    private String makePutationArray(char consChar, char mutChar, int numLines){
	StringBuffer bf = new StringBuffer();
	int mutLine = (int)(Math.random()*numLines);
	for(int i=0; i<numLines;i++){
	    if(i== mutLine)
		bf.append(mutChar + " ");
	    else
		bf.append(consChar + " ");
	}
	return bf.toString();
    }

    
    /*
     *
     * getSubSeq(int, int) method gets the substring of current seq
     * beginIndex is inclusive and end index is exclusive.
     */
    public String getSubSeq(int begin, int end){
	return seqBuffer.substring(begin-1, end-1);
    }

    public String getSubSeqWithDelim(int begin, int end, String delim){
	String subseq = this.getSubSeq(begin,end);
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<subseq.length();i++){
	    bf.append(subseq.charAt(i));
	    if(i<(subseq.length()-1))
		bf.append("\t");
	}
	return bf.toString();
    }	
    public Seq getSubSeqObj(int begin, int end){
	return new Seq(this.seqName, this.getSubSeq(begin, end));
    }
    
    public int countGC(){
	int count = 0;
	for(int i=0; i<this.seqBuffer.length(); i++){
	    char curChar = seqBuffer.charAt(i);
	    if(curChar == 'G' || curChar == 'C' 
	       || curChar == 'g' || curChar == 'c')
		count++;
	}
	return count;
    }

    /*
     * revSeq() method returns the reverse complement sequence of the sequens held by the Seq Obj
     *
     * @return outBuffer StringBuffer holding the reverse complement seq of the seq held by the Seq Obj
     */
    public StringBuffer revSeq(){
        StringBuffer tempBuffer = revSeqInternal();
        StringBuffer outBuffer = new StringBuffer();
        int lowInd = 0;
        for(int i=0; i<tempBuffer.length();i++){
            if((i % 50) == 49){
                outBuffer.append(tempBuffer.substring(lowInd,i+1));
                outBuffer.append("\n");
                lowInd += 50;
            }
        }
	if( (tempBuffer.length()%50) != 0){
	    outBuffer.append(tempBuffer.substring(lowInd));
	    outBuffer.append("\n");
	}
        return outBuffer;
    }
    
    public String getGeneSequenceAndAppendXBases(Gene g, int x){
	if(x == 0)
	    return this.getGeneSequence(g);
	int s = g.getStart()-x;
	int soff = 0;
	int e = g.getEnd()+x;
	int eoff = 0;
	if(s < 1){
	    soff = -1*(s) + 1;
	    s = 1;
	}
	if(e > this.getSeqLength()){
	    eoff = e - this.getSeqLength();
	    e = this.getSeqLength();
	}
	
	if(g.getContig().equals(this.seqName)){
	    String tmp = (soff>0 ? this.getSubSeq(this.getSeqLength()-soff+1,this.getSeqLength()+1) : "")
		+ this.getSubSeq(s, e+1)
		+ (eoff>0 ? this.getSubSeq(1,1+eoff) : "");
	    if(g.fwd())
		return tmp;
	    else
		return revComp(new StringBuffer(tmp));
	}else
	    return null;
    }

    public String getGeneSequence(Gene g){
	if(g.getContig().equals(this.seqName)){
	    if(g.fwd())
		return this.getSubSeq(g.getStart(), g.getEnd()+1);
	    else
		return revComp(new StringBuffer(this.getSubSeq(g.getStart(), g.getEnd()+1)));
	}else
	    return null;
    }
    
    public Seq getGeneSequenceObj(Gene g){
	return new Seq(g.getGName(), new StringBuffer(this.getGeneSequence(g)));
    }

    public static String revComp(StringBuffer sb){
	StringBuffer tempBuffer = sb.reverse();
	for(int i=0; i<tempBuffer.length();i++){
            char curChar = tempBuffer.charAt(i);
	    tempBuffer.setCharAt(i, rev(curChar));
        }
	return tempBuffer.toString();
    }
    
    private StringBuffer revSeqInternal(){
        StringBuffer tempBuffer = seqBuffer.reverse();
        for(int i=0; i<tempBuffer.length();i++){
            char curChar = tempBuffer.charAt(i);
	    tempBuffer.setCharAt(i, rev(curChar));
        }
        return tempBuffer;
    }

    /*public SynType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	
      }*/

    public MutType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons){
	return this.isSyn(start,fwd,position,mutbase,codons,Character.toUpperCase(seqBuffer.charAt(position-1)));
    }
    
    
    
    /*
     * Changed to carry mutType information
     *
     */
    public MutType isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	try{
	    //System.out.print(mutbase + " ");
	    mutbase = Character.toUpperCase(mutbase);
	    //if(mutbase == consensus)
	    //		return 1;
	    //else{
	    int codonPos = this.getStartOfGene(start,fwd,position);
	    StringBuffer triple = new StringBuffer();
	    StringBuffer refTriple = new StringBuffer();
	    if(fwd){
		if(codonPos == 0){
		    refTriple.append(consensus);
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		    triple.append(mutbase);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		}else if(codonPos == 1){
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    refTriple.append(consensus);
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    triple.append(mutbase);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		}else if(codonPos == 2){
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
		    refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    refTriple.append(consensus);
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
		    triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
		    triple.append(mutbase);
		}
	    }else{//reverse
		if(codonPos == 0){
		    refTriple.append(rev(consensus));
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		    triple.append(rev(mutbase));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		}else if(codonPos == 1){
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    refTriple.append(rev(consensus));
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    triple.append(rev(mutbase));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		}else if(codonPos == 2){
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
		    refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    refTriple.append(rev(consensus));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
		    triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
		    triple.append(rev(mutbase));
		}
	    }
	    
	    return codons.isSyn(refTriple.toString(), triple.toString());
	    
		//}
	}catch(Exception e){
	    System.err.println("|" + start + "\t" + fwd + "\t" + position + "\t" + mutbase +"\t" +codons);//t start, boolean fwd, int position, char mutbase, Codons codons);
	    //System.err.println(ct.toStringBuffer().toString());
	    e.printStackTrace();
	    System.exit(1);
	}
	return null;
    }




    /*
     * start: in case of (-) stranded gene, start here is the actual start --> so 2nd column of fraggenescan output
     * in other words, start is always the 5' end of the gene(start position of the gene)
     * 
     * RETURN 1 if syn, 0 if non-synonymous, -1 if codon contains bases other than a t g c
     */
    /*public int isSyn(int start, boolean fwd, int position, char mutbase, Codons codons, char consensus){
	try{
	    //System.out.print(mutbase + " ");
	    mutbase = Character.toUpperCase(mutbase);
	    if(mutbase == consensus)
		return 1;
	    else{
		int codonPos = this.getStartOfGene(start,fwd,position);
		StringBuffer triple = new StringBuffer();
		StringBuffer refTriple = new StringBuffer();
		if(fwd){
		    if(codonPos == 0){
			refTriple.append(consensus);
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
			triple.append(mutbase);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position+1)));
		    }else if(codonPos == 1){
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			refTriple.append(consensus);
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			triple.append(mutbase);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position)));
		    }else if(codonPos == 2){
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
			refTriple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			refTriple.append(consensus);
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-3)));
			triple.append(Character.toUpperCase(seqBuffer.charAt(position-2)));
			triple.append(mutbase);
		    }
		}else{//reverse
		    if(codonPos == 0){
			refTriple.append(rev(consensus));
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
			triple.append(rev(mutbase));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-3))));
		    }else if(codonPos == 1){
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
			refTriple.append(rev(consensus));
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
			triple.append(rev(mutbase));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position-2))));
		    }else if(codonPos == 2){
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
			refTriple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
			refTriple.append(rev(consensus));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position+1))));
			triple.append(rev(Character.toUpperCase(seqBuffer.charAt(position))));
			triple.append(rev(mutbase));
		    }
		}
		
		return codons.isSyn(refTriple.toString(), triple.toString());
		
	    }
	}catch(Exception e){
	    System.err.println("|" + start + "\t" + fwd + "\t" + position + "\t" + mutbase +"\t" +codons);//t start, boolean fwd, int position, char mutbase, Codons codons);
	    //System.err.println(ct.toStringBuffer().toString());
	    e.printStackTrace();
	    System.exit(1);
	}
	return 0;
    }
    */

    private static char rev(char curChar){
	if(curChar == 'T')
	    return 'A';
	else if(curChar == 't')
	    return 'a';
	else if(curChar == 'A')
	    return 'T';
	else if(curChar == 'a')
	    return 't';
	else if(curChar == 'G')
	    return 'C';
	else if(curChar == 'g')
	    return 'c';
	else if(curChar == 'C')
	    return 'G';
	else if(curChar == 'c')
	    return 'g';
	else
	    return curChar;
    }

    private String getCodon(boolean fwd, int codonPos, int position, char mutbase){
	if(codonPos == 0)
	    ;
	return new String();
    }
    
    /*
     *
     * @param
     * start    --> start position of target gene
     * fwd      --> true if fwd, rev otherwise
     * position --> position of polymorphism
     *
     * @returns
     * integer value of codon position for the given polymorphic pos.
     * 0 --> 1st pos of triplet
     * 1 --> 2nd pos of triplet
     * 2 --> 3rd pos of triplet
     *
     * start = 5
     * target position = 12
     * 
     * (fwd)       * (start pos)        ** (target pos)
     * (rev)       * (target pos)       ** (start pos)
     * 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18   --> position on genome
     * 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17   --> index 
     *
     * fwd strand
     * 12th position is the 2nd position of a codon.
     * --> (targetPos-start)%3 = (12-5)%3 = 1 
     *
     * rev strand
     * 5th pos is the 2nd pos of a codon 
     * --> (start - targetPos)%3 = (12-5)%3 = 1
     *
     */
    private int getStartOfGene(int start, boolean fwd, int position){
	if(fwd){//if the gene is on the fwd strand
	    return (position-start)%3;
	}else{//rev strand
	    return (start-position)%3;
	}
    }

    private int mapToRevIndex(int pos){
	return this.seqBuffer.length()-pos-1;
    }


    public void GC_AT_ratio(int binsize){

	for(int i=0;i<this.seqBuffer.length();i=i+binsize){
	    String segment = "";
	    if(i+binsize > this.seqBuffer.length())
		segment = this.seqBuffer.substring(i);
	    else
		segment = this.seqBuffer.substring(i,i+binsize);
	    int count = 0;
	    for(int j=0; j<segment.length();j++){
		char curChar = segment.charAt(j);
		if(curChar == 'G' || curChar == 'C'|| curChar == 'g' || curChar == 'c')
		    count++;
	    }
	    if(i+binsize > this.seqBuffer.length())
		System.out.println((i+1) + "\t" + (i+binsize) + "\t" + count + "\t" + (this.seqBuffer.length()-i-count));
	    else
		System.out.println((i+1) + "\t" + (i+binsize) + "\t" + count + "\t" + (binsize-count));
	}
    }

    public static void main(String[] args){
	Hashtable<String, Seq> hash = new FastaReader().parseFasta(args[0]);
	Seq s = hash.get("ecoli");
	s.GC_AT_ratio(Integer.parseInt(args[1]));
    }

    public int[] getATCGXArray(int start, int end, int[] atcgxArray){
	//int[] atcgxArray = new int[5];
	String tmp = this.getSubSeq(start,end+1);
	int A = 0;
	int T = 0;
	int C = 0;
	int G = 0;
	int X = 0;
	for(int i=0;i<tmp.length();i++){
	    if(tmp.charAt(i) == 'A')
		A++;
	    else if(tmp.charAt(i) == 'T')
		T++;
	    else if(tmp.charAt(i) == 'C')
		C++;
	    else if(tmp.charAt(i) == 'G')
		G++;
	    else
		X++;
	}
	
	atcgxArray[0] += A;
	atcgxArray[1] += T;
	atcgxArray[2] += C;
	atcgxArray[3] += G;
	atcgxArray[4] += X;
	
	return atcgxArray;
    }

    public static int[] getATCGXArray(String tmp, int[] atcgxArray){
	int A = 0;
	int T = 0;
	int C = 0;
	int G = 0;
	int X = 0;
	for(int i=0;i<tmp.length();i++){
	    if(tmp.charAt(i) == 'A')
		A++;
	    else if(tmp.charAt(i) == 'T')
		T++;
	    else if(tmp.charAt(i) == 'C')
		C++;
	    else if(tmp.charAt(i) == 'G')
		G++;
	    else
		X++;
	}
	
	atcgxArray[0] += A;
	atcgxArray[1] += T;
	atcgxArray[2] += C;
	atcgxArray[3] += G;
	atcgxArray[4] += X;
	
	return atcgxArray;
    }


    public void printATCGCounts(){
	this.printATCGCounts(1,this.seqBuffer.length());
    }

    public void printATCGBothStrands(int start, int end){
	System.out.println("#A#T#C#G(Referance strand)[ " + start + " , "+end+" ]");
	this.printATCGCounts(start, end);
	
	System.out.println("\t#A#T#C#G(opposite strand)[ " + start + " , " + end + " ]");
	this.printATCGCountsRevCompBasedOnForwardCoordinates(start, end);
    }

    public int[] countATCGCounts(int start, int end){
	String tmp = this.getSubSeq(start,end+1);
	int[] atcgx = new int[5];
	//int A = 0;
	//int T = 0;
	//int C = 0;
	//int G = 0;
	//int X = 0;
	for(int i=0;i<tmp.length();i++){
	    if(tmp.charAt(i) == 'A')
		atcgx[0]++;//A++;
	    else if(tmp.charAt(i) == 'T')
		atcgx[1]++;//T++;
	    else if(tmp.charAt(i) == 'C')
		atcgx[2]++;//C++;
	    else if(tmp.charAt(i) == 'G')
		atcgx[3]++;//G++;
	    else
		atcgx[4]++;//X++;
	}
	return atcgx;
    }

    public void printATCGCounts(int start, int end){
	int[] atcgx = new int[5];
	if(end >= start){
	    atcgx = countATCGCounts(start, end);
	}else{
	    int[] tmp = countATCGCounts(start, this.seqBuffer.length());
	    int[] tmp2 = countATCGCounts(1, end);
	    for(int i=0;i<atcgx.length;i++)
		atcgx[i] = tmp[i] + tmp2[i];
	}
	
	System.out.println("\t" + atcgx[0] + "\t" + atcgx[1] + "\t" + atcgx[2] + "\t" + atcgx[3] +"\t" + atcgx[4] );
	
    }

    public int[] countATCGCountsRevCompBasedOnForwardCoordinates(int start, int end){
	String tmp = this.getSubSeq(start,end+1);
	int[] atcgx = new int[5];
	//int A = 0;
	//int T = 0;
	//int C = 0;
	//int G = 0;
	//int X = 0;
	//we are counting the reverse so T for A, A for T, C for G, G for C.
	for(int i=0;i<tmp.length();i++){
	    if(tmp.charAt(i) == 'A')
		atcgx[1]++;//T++;
	    else if(tmp.charAt(i) == 'T')
		atcgx[0]++;//A++;
	    else if(tmp.charAt(i) == 'C')
		atcgx[3]++;//G++;
	    else if(tmp.charAt(i) == 'G')
		atcgx[2]++;//C++;
	    else
		atcgx[4]++;//X++;
	}
	return atcgx;	
    }

    public void printATCGCountsRevCompBasedOnForwardCoordinates(int start, int end){
	int[] atcgx = new int[5];
	if(end >= start){
	    atcgx = countATCGCountsRevCompBasedOnForwardCoordinates(start, end);
	}else{
	    int[] tmp = countATCGCountsRevCompBasedOnForwardCoordinates(start, this.seqBuffer.length());
	    int[] tmp2 = countATCGCountsRevCompBasedOnForwardCoordinates(1, end);
	    for(int i=0;i<atcgx.length;i++)
		atcgx[i] = tmp[i] + tmp2[i];
	}
	
	System.out.println("\t" + atcgx[0] + "\t" + atcgx[1] + "\t" + atcgx[2] + "\t" + atcgx[3] +"\t" + atcgx[4] );
	
    }

    /*
    //start and end are both inclusive and 1-base cordinate.
    public void countTriplets(int start, int end, boolean fwd){
	String tmpSeq = null;
	if(start > end){
	    if(fwd){
		tmpSeq = getSubSeq(start,this.seqBuffer.length()+1) 
		    + getSubSeq(1, end+1); 
	    }else{
		tmpSeq = Seq.revComp(new StringBuffer(getSubSeq(start,this.seqBuffer.length()+1) + getSubSeq(1,end+1)));
	    }
	}//wraparound
	else{
	    if(fwd)
		tmpSeq = getSubSeq(start,end+1);
	    else
		tmpSeq = Seq.revComp(getSubSeqObj(start,end+1).getSeq());
	}
	
	int[] tripletCounts = new int[64]; 
	
	//for(int i=0; i<tmpSeq.length()-2; i++){
	//    tripletCounts[Triplets.triplet2Index(tmpSeq.substring(i,i+3))]++;
	//}
	
	int skippedTriplets = 0;
	for(int i=0; i<tmpSeq.length()-2; i++){
	    if(isNFree(tmpSeq.substring(i,i+3)))
		tripletCounts[Triplets.triplet2Index(tmpSeq.substring(i,i+3))]++;
	    else
		skippedTriplets++;
	}
	
	
	System.out.println("SKIPPED\t" + skippedTriplets);
	Triplets.printCounts(tripletCounts);
	
	}*/

        //start and end are both inclusive and 1-base cordinate.
    public int[] countTriplets(int start, int end, boolean fwd){
	String tmpSeq = null;
	if(start > end){
	    if(fwd){
		tmpSeq = getSubSeq(start,this.seqBuffer.length()+1) 
		    + getSubSeq(1, end+1); 
	    }else{
		tmpSeq = Seq.revComp(new StringBuffer(getSubSeq(start,this.seqBuffer.length()+1) + getSubSeq(1,end+1)));
	    }
	}//wraparound
	else{
	    if(fwd)
		tmpSeq = getSubSeq(start,end+1);
	    else
		tmpSeq = Seq.revComp(getSubSeqObj(start,end+1).getSeq());
	}
	
	int[] tripletCounts = new int[65]; //last index contains : #skipped triplet
	
	//for(int i=0; i<tmpSeq.length()-2; i++){
	//    tripletCounts[Triplets.triplet2Index(tmpSeq.substring(i,i+3))]++;
	//}
	
	//int skippedTriplets = 0;
	for(int i=0; i<tmpSeq.length()-2; i++){
	    if(isNFree(tmpSeq.substring(i,i+3)))
		tripletCounts[Triplets.triplet2Index(tmpSeq.substring(i,i+3))]++;
	    else
		tripletCounts[64]++;//skippedTriplets++;
	}
	
	
	//System.out.println("SKIPPED\t" + skippedTriplets);
	//Triplets.printCounts(tripletCounts);
	return tripletCounts;
    }

    public static boolean isNFree(String triplet){
	if( isNucleotide(triplet.charAt(0)) 
	   && isNucleotide(triplet.charAt(1)) 
	   && isNucleotide(triplet.charAt(2)) )
	    return true;
	return false;
    }
    
    public static boolean isNucleotide(char b){
	if(Triplets.char2digit(b) < 0)
	    return false;
	return true;
    }
	
}

class IntegerCounter{
    
    private int count;
    private int GCCount;
    private int ATCount;
    
    public IntegerCounter(){
	this(0);
    }
    
    public IntegerCounter(int initial){
	this.count = initial;
    }
    
    public int getCount(){
	return this.count;
    }
    
    public int getGCCount(){
	return this.GCCount;
    }

    public int getATCount(){
	return this.ATCount;
    }
    
    public void add1(){
	this.count++;
    }
    
    public void add1(char c){
	this.count++;
	if(c == 'G' || c == 'g' || c == 'C' || c == 'c')
	    this.GCCount++;
	else
	    this.ATCount++;
    }
    
    public void sub1(){
	this.count--;
    }
    
    public void addN(int n){
	this.count += n;
    }
    
    
}


class RunLengthGeneCounter{

    private int runLength;
    //these two List should be synced for indicies.
    private ArrayList<Gene> geneList;
    private ArrayList<Integer> positionList;
    private ArrayList<Character> charList;
	
    public RunLengthGeneCounter(int len){
	this.runLength = len;
	this.geneList = new ArrayList<Gene>();
	this.positionList = new ArrayList<Integer>();
	this.charList = new ArrayList<Character>();
    }
    
    public void addPosition(int pos, char base, ArrayList<Gene> genes){
	boolean found = false;
	for(int i=0; i<genes.size();i++){
	    if(genes.get(i).contains(pos)){
		this.geneList.add(genes.get(i));
		this.positionList.add(new Integer(pos));
		this.charList.add(new Character(base));
		found = true;
		break;
	    }
	}
	
	/*if(!found){
	   this.geneList.add(null);
	    this.positionList.add(null);
	    }*/
    }
    
    public void print(){
	for(int i=0; i<this.geneList.size();i++){
	    if(this.geneList.get(i) != null){
		System.out.println(this.runLength + "\t" + this.charList.get(i) + "\t" + this.positionList.get(i).intValue() + "\t" + this.geneList.get(i).toStringWithGeneName());
		
	    }
	}
    }

}    
