package gov.ca.water.gisutil;
import java.io.*;
import java.util.Vector;

import org.picocontainer.defaults.ObjectReference;

/**
 *
 */
public class AsciiFileWriter{
	FileWriter _aOutFile = null;
	BufferedWriter _asciiOut = null;
	String _filename;
	Vector _allLines = null;

	public AsciiFileWriter(String filename){
		_allLines = new Vector();
		_filename=filename;
		open();
	}//constructor

	public AsciiFileWriter(String filename, boolean append){
		_allLines = new Vector();
		_filename=filename;
		if(append){
			openAppend();
		}else{
			open();
		}
	}

	protected void openAppend(){
		try{
			_aOutFile = new FileWriter(_filename,true);
			_asciiOut = new BufferedWriter(_aOutFile);
		}catch(IOException e){
			System.out.println("error occurred while opening file "+_filename + e.getMessage());
		}

	}

	protected void open(){
		try{
			_aOutFile = new FileWriter(_filename);
			_asciiOut = new BufferedWriter(_aOutFile);
		}catch(IOException e){
			System.out.println("error occurred while opening file "+_filename + e.getMessage());
		}
	}//open


	public void writeLine(String line){
		try{
			_asciiOut.write(line);
			_asciiOut.newLine();
		}catch(Exception e){
			System.out.println("exception caught in AsciiFileWriter.write while reading ");
			System.out.println("file "+_filename+e.getMessage());
		}
	}//writeLine

	public void close(){
		try{
			_asciiOut.close();
		}catch(java.io.IOException e){
			System.out.println("exception caught while trying to close file: "+
					_filename+e.getMessage());
		}
	}//close

}//AsciiFileWriter
