package gov.ca.water.gisutil;

import java.awt.Polygon;

import DWR.CSDP.CsdpFunctions;
import DWR.CSDP.ResizableIntArray;

/**
 * NOT USED.
 * 
 * Reads a network file that was created to outline data to extract.
 * Then read a bathymetry data containing the data to extract.
 * write the data that exists inside the polygon to another bathymetry (.prn) file.
 * @author btom
 *
 */
public class ExtractShipChannelLeveesFromYoloBypassDEM {
	public static void main(String[] artgs) {
		new ExtractShipChannelLeveesFromYoloBypassDEM();
	}
	
	public ExtractShipChannelLeveesFromYoloBypassDEM() {
		Polygon boundaryPolygon = getBoundaryPolygon();
		extractShipChannelData(boundaryPolygon);
	}

	/*
	 * Read the network file that outlines the data and create a polygon object using the centerline coordinates as vertices.
	 */
	private Polygon getBoundaryPolygon() {
		String directory = "d:/dsm2gisreference/csdpFiles/final/";
//		String filename = "sacDeepShipLeveesFromyoloDEM.cdn";
		String filename = "tomPaineExtract.cdn";
		AsciiFileReader asciiFileReader = new AsciiFileReader(directory+filename);
		int i=0;
		ResizableIntArray xValues = new ResizableIntArray();
		ResizableIntArray yValues = new ResizableIntArray();
		while(true) {
			String line = asciiFileReader.getNextLine();
			if(line==null) break;
			if(line.indexOf(";") < 0 && line.indexOf("\"")<0) {
				String[] parts = line.trim().split(",");
				if(parts.length==2) {
					xValues.put(i, (int)(CsdpFunctions.feetToMeters(Double.parseDouble(parts[0]))));
					yValues.put(i, (int)(CsdpFunctions.feetToMeters(Double.parseDouble(parts[1]))));
					i++;
				}
			}
		}
		System.out.println("Done reading yolo bathymetry file");
		int numVertices = i;
		int[] x = new int[numVertices];
		int[] y = new int[numVertices];
		for(int j=0; j<numVertices; j++) {
			x[j] = xValues.get(j);
			y[j] = yValues.get(j);
		}
		
		
		asciiFileReader.close();
		Polygon polygon = new Polygon(x,y,numVertices);
		System.out.println("Polygon="+polygon);
		return polygon;
	}//getBoundaryPolygon

	/*
	 * Read the bathymetry file, and write the points that are contained within the polygon to another bathymetry file.
	 */
	private void extractShipChannelData(Polygon boundaryPolygon) {
		String directory = "d:/dsm2gisreference/csdpFiles/final/";
//		String filename = "dem_yolo_merge3_mod_20161202.prn";
//		String outfilename = "sacShipChanLeveeData.prn";
		String filename = "withDEMSouthDelta.prn";
		String outfilename = "allTomPaine.prn";
		AsciiFileReader asciiFileReader = new AsciiFileReader(directory+filename);
		AsciiFileWriter asciiFileWriter = new AsciiFileWriter(directory+outfilename);
//		String headerLines = "";
//		String linesToWrite = "";
		System.out.println("extractShipChannelData");
		int i=0;
		while(true) {
			if(i % 10000 == 0) System.out.println(i);
			String line = asciiFileReader.getNextLine();
			if(line==null) break;
			if(line.indexOf(";")>=0) {
//				headerLines += line+"\n";
			}else {
				String[] parts = line.split(",| ");
				double x = Double.parseDouble(parts[0]);
				double y = Double.parseDouble(parts[1]);
				if(boundaryPolygon.contains(x, y)) {
//					linesToWrite+=line+"\n";
					asciiFileWriter.writeLine(line);
				}
			}
			i++;
		}
//		asciiFileWriter.writeLine(headerLines);
//		asciiFileWriter.writeLine(linesToWrite);
		asciiFileReader.close();
		asciiFileWriter.close();
	}

}//class ExtractShipChannelLeveesFromSuisunMarshDEM
