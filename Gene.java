public class Gene implements Comparable<Gene>{

    private String contig;
    private int start;//smaller number
    private int end;//larger number
    private boolean fwd; //true if fwd, false otherwise.
    private String gName;

    public int compareTo(Gene other){
	
	if(this.start < other.start)
	    return -1;
	else if(this.start == other.start)
	    return 0;
	else
	    return 1;
	    
    }

    private String pttline;
    /*
    public Gene(String ct, String line){
	String[] tokens = line.split("\\t");
	}*/

    public void merge(int e, String geneName){
	this.end = e;
	this.gName = this.gName + "_" + geneName;
    }

    public void setEnd(int e){
	this.end = e;
    }

    public Gene(String ct, String coreGeneLine){
	String[] coreGeneToks = coreGeneLine.split("\\t");
	this.contig =ct;
	this.start = Integer.parseInt(coreGeneToks[5]);
	this.end = Integer.parseInt(coreGeneToks[6]);
	this.fwd = (coreGeneToks[4].equals("Clockwise") ? true : false);
	this.gName = coreGeneToks[1];
    }

    public Gene(String ct, int st, int end, boolean fwd, String geneName){
	this.contig = ct;
	this.start = st;
	this.end = end;
	this.fwd = fwd;
	this.gName = geneName;
    }

    public Gene(String ct, int st, int end, boolean fwd, String geneName, String pl){
	this.contig = ct;
	this.start = st;
	this.end = end;
	this.fwd = fwd;
	this.gName = geneName;
	this.pttline = pl;
    }

    public String getPttline(){
	return this.pttline;
    }

    public String getGeneSequenceAndAppendXBases(Seq seq, int x){
	return seq.getGeneSequenceAndAppendXBases(this, x);
    }

    public String getGeneSequence(Seq seq){
	return seq.getGeneSequence(this);
    }

    public String isLeading(int oriS, int oriE, int terS, int terE, int gLen){
	if( (this.end > terS && this.end <terE) || (this.start >terS && this.start <terE) )
	    return "TRAP";

	if(terS < oriE){ //zero is in between oriE and terS ------oriE----0-----terS : oriE > terS
	    if(this.fwd){
		//----oriS---oriE--====>--0--------terS  OR ----oriS---oriE--====0===>--------terS
		if(this.start > oriE && (this.end <= gLen || this.end < terS) ) 
		    return "Leading";
		//----oriS---oriE----0--=====>-----terS
		else if(this.start >= 0 && this.end <terS)
		    return "Leading";
		else
		    return "Lagging";
	    }else{
		//---terS----terE--<======------oriS----oriE----0----
		if(this.start >terE && this.end <oriS)
		    return "Leading";
		else
		    return "Lagging";
		
	    }
	}else{ // -------0---oriS---oriE-------------------terS---terE--------
	    if(this.fwd){
		// ------0----oriS--oriE---=====>----terS---terE-------
		if(this.start > oriE && this.end < terS)
		    return "Leading";
		else
		    return "Lagging";
	    }else{
		//---terS----terE----<=====----0---oriS----oriE--- OR ---terS----terE----<===0===---oriS----oriE---
		if(this.start > terE && (this.end <=gLen || this.end < oriS) )
		    return "Leading";
		//---terS----terE----0-<====--oriS----oriE---
		else if(this.start > 0 && this.end <oriS)
		    return "Leading";
		else
		    return "Lagging";
	    }
	}
    }

    public Gene(String ct, int st, int end, boolean fwd, int frame){
	if(frame == 3)
	    frame = 0;
	if( (end-st+1)%3 != 0 ){
	    if(fwd){
		if( st%3 == frame)
		    end = end - ((end-st+1)%3);
		else if((end + 1) %3 == frame)
		    st = st + ((end-st+1)%3);
		else{
		    int oldst = st;
		    st = ( (st+1)%3 == frame ? st+1 : st+2);
		    end = ( (oldst+1)%3 == frame ? (end - ((end-oldst+1)%3) + 1) : (end- ((end-oldst+1)%3) + 2));
		}
	    }else{
		if( (end+1)%3 == frame)
		    st = st + ((end-st+1)%3);
		else if(st%3 == frame)
		    end = end - ((end-st+1)%3);
		else{
		    int oldst = st;
		    st = ( (st+1)%3 == frame ? st+1 : st+2);
		    end = ( (oldst+1)%3 == frame ? (end - ((end-oldst+1)%3) + 1) : (end- ((end-oldst+1)%3) + 2));
		}
	    }
	    //System.err.println(ct + "\t" + st + "\t" + end + "\t" + fwd + "\t" + frame); 
	}
	
	if( (end-st+1)%3 != 0)
	    System.err.println("Are you  Kidding me?");
	this.contig = ct;
	this.start = st;
	this.end = end;
	this.fwd = fwd;
	
	
    }
    
    public Gene(String ct, String[] fragGeneScanTokens){
	this(ct, 
	     Integer.parseInt(fragGeneScanTokens[0]), 
	     Integer.parseInt(fragGeneScanTokens[1]), 
	     ( fragGeneScanTokens[2].equals("+") ? true : false )
	     , Integer.parseInt(fragGeneScanTokens[3]));
    }

    public String getGName(){
	return this.gName;
    }
    public String getContig(){
	return this.contig;
    }

    public String getUniqueString(){
	return this.contig + "_" + this.gName + "_" + this.start + "_" + this.end + "_" + (this.end-this.start+1);
    }
    
    public int getStart(){
	return this.start;
    }
    
    public int getEnd(){
	return this.end;
    }
    
    public int length(){
	return this.end - this.start + 1;
    }
    
    public boolean fwd(){
	return this.fwd;
    }
    public String directionStr(){
	if(this.fwd)
	    return "fwd";
	else
	    return "rev";
    }

    public String toString(){
	return this.contig + "\t" + this.start + "\t" + this.end + "\t" + this.fwd;
    }
    
    public String toStringWithGeneName(){
	return this.contig + "\t" + this.gName + "\t" + this.start + "\t" + this.end + "\t" + this.fwd;
    }

    public String getHeader(){
	return this.start + "\t" + this.end + "\t" + (this.fwd ? "+" : "-") + "\t" + (this.end - this.start + 1) + "\t" + this.gName;
    }

    public String toPttString(){
	return this.start + ".." + this.end + "\t" + (this.fwd ? "+" : "-") + "\t-\t-\t" + this.gName + "\t-\t-\t-\t-"; 
    }

    public boolean contains(int pos){
	if( (pos >= this.start) &&
	    (pos <= this.end) )
	    return true;
	return false;
    }

}
