package gov.ca.water.gisutil;
import java.io.*;
import java.util.Vector;

/**
 * Usage:
 * int i=0;
 * while(true){
 *   String line = rmaOutReader.getNextLine();
 *   if(line==null) break;
 *   String parts[] = line.split("\\t+");
 *	    if(startTime==null) startTime = parts[1];
 * // 	    String printString = "parts:";
 * // 	    for(int i=0; i<=parts.length-1; i++){
 * // 		printString += parts[i]+",";
 * // 	    }
 * 	    double stage, xvel, yvel;
 * 	    try{
 * 		stage = Double.parseDouble(parts[2]);
 *	    }catch (java.lang.NumberFormatException e){
 *		stage = Constants.MISSING;
 *	    }
 *	    try{
 *		xvel = Double.parseDouble(parts[3]);
 *	    }catch (java.lang.NumberFormatException e){
 *		xvel = Constants.MISSING;
 *	    }
 *	    try{
 *		yvel = Double.parseDouble(parts[4]);
 *	    }catch (java.lang.NumberFormatException e){
 *		yvel = Constants.MISSING;
 *	    }
 *
 *	    stageRDA.put(i, stage);
 *	    xvelRDA.put(i, xvel);
 *	    yvelRDA.put(i, yvel);
 *	    //System.out.println(printString);
 *	    i++;
 *	}
 *	rmaOutReader.close();
 */
public class AsciiFileReader{
    private LineNumberReader _asciiIn;
    private FileReader _aInFile;
    
    public AsciiFileReader(String path){
	try{
	    _aInFile = new FileReader(path);
	    _asciiIn = new LineNumberReader(_aInFile);
	}catch(IOException e){
	    System.out.println("IOException caught in AsciiFileReader constructor: "+e);
	}
    }

//     /*
//      * Won't work...
//      */
//     public void rewind(){
//  	_asciiIn.setLineNumber(0);
//  	try{
//  	    _asciiIn.reset();
//  	}catch(java.io.IOException e){
//  	    System.out.println("error in AsciiFileReader.rewind: IOException: "+e);
//  	}
//     }

    /*
     * copied from http://stackoverflow.com/questions/326390/how-to-create-a-java-string-from-the-contents-of-a-file
     */
    public static String readFileIntoString( String file ) throws IOException {
	BufferedReader reader = new BufferedReader( new FileReader (file));
	String line  = null;
	StringBuilder stringBuilder = new StringBuilder();
	String ls = System.getProperty("line.separator");
	while( ( line = reader.readLine() ) != null ) {
	    stringBuilder.append( line );
	    stringBuilder.append( ls );
	}
	return stringBuilder.toString();
    }//readFileIntoString

    public String getNextLine(){
	String returnString = null;
	try{
	    returnString =  _asciiIn.readLine();
	}catch(Exception e){
	    System.out.println("Error while reading file. "+e);
	}
	return returnString;
    }

    public String[] getLines(String[] lines, String commentCharacter){
	Vector<String> linesVector = new Vector<String>();

	int numLines=0;
	String line=" ";
	while(line!=null){
	    line=getNextLine();
	    if(line != null){
		if(line.indexOf(commentCharacter)!=0){
		    numLines++;
		    linesVector.add(line);
		}
	    }else{
		break;
	    }
	}

	//System.out.println("numLines="+numLines);
	lines = new String[numLines];
	//System.out.println("length of string array="+lines.length);
	for(int i=0; i<=numLines-1; i++){
	    lines[i]=linesVector.get(i);
	    //System.out.println("added, line="+lines[i]);
	}
	close();
	return lines;
    }//getLines

    public void close(){
	try{
	    _asciiIn.close();
	}catch(Exception e){
	    System.out.println("Error while trying to close file");
	}
    }

}
