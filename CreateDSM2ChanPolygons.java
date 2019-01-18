package gov.ca.water.gisutil;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.Vector;

import org.geotools.data.Transaction;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FeatureSource;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureImpl;
import org.geotools.map.FeatureLayer;
import org.geotools.map.Layer;
import org.geotools.map.MapContent;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.swing.data.JFileDataStoreChooser;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;

import DWR.CSDP.App;
import DWR.CSDP.Centerline;
import DWR.CSDP.CenterlinePoint;
import DWR.CSDP.CsdpFileMetadata;
import DWR.CSDP.CsdpFrame;
import DWR.CSDP.CsdpFunctions;
import DWR.CSDP.Network;
import DWR.CSDP.NetworkOutput;


/**
 * Read a shapefile containing waterbody polygons--currently using hydrography.shp
 * 
 * Polygons will have to be widened after creating them, because we don't know what the water levels were at the time of the lidar survey.
 * 
 * Read CSDP network file
 * Create polygon layer containing data that belongs to each DSM2 channel
 * 
 * Polygons may need to be enlarged, because the stage when the surveys were taken to create the hydrography layer is unknown.
 * Possible procedure: 
 * 1. Calculate volume using GIS between each polygon and the 2m DEM layer.
 * 2. Enlarge polygons by moving points away from CSDP centerline.
 * 3. Calculate volume again and compare. 
 * 		a. If no change, then our points are inside the levee for the given elevation, and we're done.
 * 		b. If the volume increases, keep moving points until it stops increasing.
 * 		c. If volume never stops increasing, they our points may be on the other side of the levee. 
 * 
 * @author btom
 *
 */
public class CreateDSM2ChanPolygons {

	/*
	 * Half of the length of the perpendicular line segments that are created at the upstream and downstream
	 * ends of each CSDP centerline
	 */
	private static final double PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH=100000.0;
	//	private Vector<MultiLine2D> _allLeveeCenterlines = new Vector<MultiLine2D>();
	//Polygons read from shapefile often have duplicate names
	private Vector<Long> _allWaterbodyPolygonObjectIDs = new Vector<Long> ();
	private Vector<String> _allWaterbodyPolygonNames = new Vector<String>();
	private Vector<Polygon> _allWaterbodyPolygons = new Vector<Polygon>();
	private Vector<String> _allCenterlineNames = new Vector<String>(); 
	private Hashtable<String, CenterlineGeometry> _allCenterlineGeometries = new Hashtable<String, CenterlineGeometry>();
	private GeometryFactory _geometryFactory = new GeometryFactory();
	private SimpleFeatureType _polygonShapefileFeatureType;
	Hashtable<String, LineString> perpendicularLineHashtable = new Hashtable<String, LineString>();

	private static final boolean DEBUG = false;


	/**
	 * GeoTools Quickstart demo application. Prompts the user for a shapefile and displays its
	 * contents on the screen in a map frame
	//        	Attribute 0: MULTILINESTRING ((-121.56761661599995 38.19186989900004, -121.56759192199996 38.19184188100007, -121.56757179399995 38.191821135000055, -121.56756256999995 38.19180980400006, -121.56754775599995 38.19179034600006, -121.56754139199995 38.191782303000025, -121.56753561399995 38.19177534700003, -121.56752532599995 38.191763672000036, -121.56749985699997 38.19173546600007, -121.56748264599997 38.191716719000055, -121.56747654199995 38.19170957100005, -121.56747372299998 38.191705942000056, -121.56747110999999 38.19170229200006, -121.56746660999994 38.19169503000006, -121.56746299199995 38.191687924000064, -121.56745795899997 38.19167718500006, -121.56745368199995 38.19166941700007, -121.56745105399995 38.191665244000035, -121.56744483599999 38.19165636500003, -121.56742774999998 38.19163360600004, -121.56742234599994 38.19162526100007, -121.56741789399996 38.19161697200008, -121.56741099199996 38.19160266600005, -121.56740808599994 38.19159702600007, -121.56740470499994 38.191590903000076, -121.56740079299999 38.19158432800003, -121.56739633199999 38.19157736400007, -121.56739132599995 38.19157008800005, -121.56738580899997 38.19156258100003, -121.56737982499999 38.191554917000076, -121.56737343099996 38.19154715600007, -121.56736669699995 38.191539341000066, -121.56735968899994 38.19153150300008, -121.56734512899999 38.19151580500005, -121.56728570199999 38.19145285800005, -121.56724085199994 38.191405700000075, -121.56720387299998 38.19136651400004, -121.56716083099997 38.19)1321893000065, -121.56714126599996 38.191301114000055, -121.56711688099995 38.19127420700005, -121.566912463 38.19104552300007, -121.566891667 38.19102245800008, -121.56686182599998 38.19098996100007, -121.56663084399997 38.190740591000065, -121.56652618499999 38.19062833000004, -121.56644377699996 38.190540250000026, -121.56617664199996 38.19025263800006, -121.56612501899997 38.19019728600006, -121.56595036199997 38.19001102300007, -121.56530724699996 38.18932677600003, -121.56503415099996 38.18903834100007, -121.56500465799996 38.18900677600004, -121.56488374799994 38.18887644400007, -121.56483867699995) 38.18882827200008, -121.56448387399996 38.188453587000026, -121.56444522699996 38.188413368000056, -121.56431173699997 38.188275770000075, -121.56427330999998 38.188235614000064, -121.56424310899996 38.18820362500003, -121.56398453799994 38.18792578500006, -121.56390585499997 38.187841598000034, -121.56363867299996 38.18755496500006, -121.56358520699996 38.18749718500004, -121.56332303099998 38.18721260800004, -121.56328505199997 38.187171867000075, -121.56319968899999 38.187082272000055, -121.56316856199999 38.18705003100007, -121.56313737199997 38.18701832800008, -121.56295793499999 38.18683724600004, -121.562934714 38.18681401300006, -121.56291201299996 38.186791725000035, -121.56289043499999 38.18677117200008, -121.56284175699994 38.18672588900006, -121.56283125799996 38.18671571500005, -121.56282122899995 38.18670550200005, -121.56281639399998 38.18670033300003, -121.56280707699995 38.18668983100008, -121.56279814999999 38.18667912400008, -121.56278047799998 38.18665739900007, -121.56277108199998 38.186646533000044, -121.56276612399995 38.18664111000004, -121.56268562199995 38.186543061000066, -121.56263499899995 38.186464915000045, -121.56261856399999 38.18643874600008))
	//        	Attribute 1: 611
	//        	Attribute 2: Brannan Andrus
	//        	Attribute 3: BALMD
	//        	Attribute 4: Non-urban
	//        	Attribute 5: Project
	//        	Attribute 6: dryland
	//        	Attribute 7: 
	//        	Attribute 8: structural, derived from 2007-08 LIDAR
	//        	Attribute 9: 
	//        	Attribute 10: Project
	//        	Attribute 11: NO
	//        	Attribute 12: 0.00738785711777        	
	 */
	public static void main(String[] args){
		new CreateDSM2ChanPolygons();
	}

	/*
	 * Constructor
	 */
	public CreateDSM2ChanPolygons() {
		//get input file path
		File inputShapefile = JFileDataStoreChooser.showOpenFile("shp", null);
		//get output file path
		String directory = "D:/dsm2GisReference/csdpFiles/final/gisVolTesting/";
		String outputShapefileFilename = "testPolygonShapefileOutput.shp";
		String shapefileDirectory = "D:/dsm2GisReference/data_source_dem_v4";
		File outputShapefile = getNewShapeFile(new File(shapefileDirectory+File.separator+outputShapefileFilename));
		createWaterbodyPolygonVector(inputShapefile);

		App csdpApp = new App();
		CsdpFileMetadata cfm = new CsdpFileMetadata();
		cfm.setHDatum(CsdpFileMetadata.UTMNAD83);
		cfm.setHUnits(CsdpFileMetadata.METERS);
		cfm.setHZone(10);
		cfm.setVDatum(CsdpFileMetadata.NAVD1988);
		cfm.setVUnits(CsdpFileMetadata.USSURVEYFEET);

		CsdpFunctions.setBathymetryMetadata(cfm);
		Network network = csdpApp.justReadNetwork(directory, "COMBINED_NETWORK_FILE.cdn");
		network.sortCenterlineNames();
		for(int i=0; i<network.getNumCenterlines(); i++) {
			String centerlineName = network.getCenterlineName(i);
			Centerline centerline = network.getCenterline(centerlineName);
			_allCenterlineNames.addElement(centerlineName);
			_allCenterlineGeometries.put(centerlineName, new CenterlineGeometry(centerline));
		}
		//geometry object will be Polygon
		Hashtable<String, Geometry> modelGridChannelPolygons = createModelGridChannelPolygons(network);
		//		String directory = "d:/dsm2GisReference/csdpFiles/final/";
		//write features to network files, so they can be verified with the csdp.
		Hashtable<String, Geometry> wbpHashtable = new Hashtable<String, Geometry>();
		for(int i=0; i<_allWaterbodyPolygonNames.size(); i++) {
			String name  = _allWaterbodyPolygonNames.get(i);
			Geometry geometry = _allWaterbodyPolygons.get(i);
			wbpHashtable.put(name, geometry);
		}

		writeFeaturesToNetworkFile(wbpHashtable, directory, "originalWaterbodyPolygons", null);
		System.out.println("wrote original waterbody polygons");
		writeFeaturesToShapefile(outputShapefile, modelGridChannelPolygons);
		System.out.println("wrote modelgridchannelpolygons");
//		Vector<String> centerlineNamesToWrite = new Vector<String>();
//		centerlineNamesToWrite.add("309");
//		centerlineNamesToWrite.add("310");
//		centerlineNamesToWrite.add("434");
		Vector<String> centerlineNamesToWrite = null; 
		writeFeaturesToNetworkFile(modelGridChannelPolygons, directory, "modelGridChannelPolygons", centerlineNamesToWrite);
		System.out.println("wrote diagnosticPolygonNetwork");
		writeNetworkWithPerpendicularLines(network, perpendicularLineHashtable, directory, "perpendicularLineNetwork");
		System.out.println("wrote perpendicular lines");
		System.out.println("Done");
		System.exit(0);
	}//constructor

	/*
	 * Read waterbody polygons from shapefile
	 */
	private void createWaterbodyPolygonVector(File file) {
		if (file == null) {
			return;
		}
		//use 8 for hydrography layer, 7 for dsm2 layer.
		int featureNameIndex = 8;

		FileDataStore store = null;
		try {
			store = FileDataStoreFinder.getDataStore(file);
			SimpleFeatureSource featureSource = store.getFeatureSource();

			// Create a map content and add our shapefile to it
			MapContent map = new MapContent();
			map.setTitle("Quickstart");

			Style style = SLD.createSimpleStyle(featureSource.getSchema());
			Layer layer = new FeatureLayer(featureSource, style);
			FeatureSource layerFeatureSource = layer.getFeatureSource();
			FeatureCollection layerFeatureCollection = layerFeatureSource.getFeatures();
			FeatureIterator<SimpleFeatureImpl> layerFeatureIterator = layerFeatureCollection.features();
			int featureIndex = 0;
			while(layerFeatureIterator.hasNext()) {
				SimpleFeatureImpl simpleFeatureImpl = layerFeatureIterator.next();
				//	        	FeatureId id = simpleFeatureImpl.getIdentifier();
				//	        	FeatureType ft = simpleFeatureImpl.getType();
				_polygonShapefileFeatureType = simpleFeatureImpl.getType();
				//	        	System.out.println("feature id, type="+id+","+ft);
				try {
					if(simpleFeatureImpl.getAttribute(0) instanceof Polygon) {
						System.out.println("It's a polygon, not a multipolygon.");
						System.exit(0);
					}
					MultiPolygon multiPolygon = (MultiPolygon) simpleFeatureImpl.getAttribute(0);
					String featureName = (String) simpleFeatureImpl.getAttribute(featureNameIndex);
					featureName+="_"+featureIndex;
					System.out.println("FeatureName = "+featureName);
					long objectID = (long) simpleFeatureImpl.getAttribute(1);
					//add each polygon in the multipolygon, one at a time (break the MultiPolygon into Polygon objects)
					int numPolygons = multiPolygon.getNumGeometries();
					//	        		if(numPolygons>1) {
					//	        			System.out.println("polygon " + objectID + ":" + featureName + " has " + numPolygons + "polygons");
					//	        		}
					for(int i=0; i<numPolygons; i++) {
						Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
						_allWaterbodyPolygonNames.addElement(featureName);
						_allWaterbodyPolygons.addElement(polygon);
						_allWaterbodyPolygonObjectIDs.addElement(objectID);

					}
				}catch(ClassCastException e1) {
					e1.printStackTrace();
					System.out.println("createWaterbodyPolygonVector: Class cast exception caught:" + e1.toString());
					System.out.println("Exiting");
					System.exit(0);
				}
				featureIndex++;
			}
			//	        System.out.println("num polygons found:"+_allWaterbodyPolygons.size());
		}catch(IOException e) {
			System.out.println("IOException caught!!!");
		}finally {
			store.dispose();
		}
		//
		//        System.out.println("added "+_allWaterbodyPolygonNames.size()+" values to vectors");
		//        System.exit(0);
	}//createWaterbodyPolygonVector

	/**
	 * Stores centerline vertex coordinates. Used to determine intersections with waterbody polygon objects.
	 * @author btom
	 *
	 */
	private class CenterlineGeometry{
		Centerline _centerline;
		LineString _centerlineLineString;
		public CenterlineGeometry(Centerline centerline) {
			_centerline = centerline;
			Coordinate[] coordinates = new Coordinate[centerline.getNumCenterlinePoints()];
			for(int i=0; i<centerline.getNumCenterlinePoints(); i++) {
				CenterlinePoint centerlinePoint = centerline.getCenterlinePoint(i);
				double centerlinePointX = CsdpFunctions.feetToMeters(centerlinePoint.getXFeet());
				double centerlinePointY = CsdpFunctions.feetToMeters(centerlinePoint.getYFeet());
				coordinates[i]= new Coordinate(centerlinePointX, centerlinePointY);
			}
			_centerlineLineString = _geometryFactory.createLineString(coordinates);
		}
		//    	public boolean intersectsPolygon(MultiPolygon polygon) {
		//    		return polygon.intersects(_centerlineLineString);
		//    	}
		//    	public boolean intersectsPolygon(Polygon polygon) {
		//    		return polygon.intersects(_centerlineLineString);
		//    	}
		public boolean intersectsPolygon(Geometry polygon) {
			return polygon.intersects(_centerlineLineString);
		}

	}//class CenterlineGeometry

	/*
	 * for each centerline in the network, find polygons that intersect, then cut and merge them
	 */
	private Hashtable<String, Geometry> createModelGridChannelPolygons(Network network) {
		String directory = "D:/dsm2GisReference/csdpFiles/final/";
		Hashtable<String, Geometry> modelGridChannelPolygons = new Hashtable<String, Geometry>(); 
		for(int i=0; i<_allCenterlineNames.size(); i++) {
			String centerlineName = _allCenterlineNames.get(i);
			CenterlineGeometry centerlineGeometry = _allCenterlineGeometries.get(centerlineName);
			//    		Vector<Polygon> allIntersectingWaterbodyPolygons = new Vector<Polygon>();
			//    		// channelPolygon will be a MultiPolygon initially, then will change to a Polygon
			Geometry channelPolygon = null;
			//find all polygons that intersect the centerline
			for(int j=0; j<_allWaterbodyPolygons.size(); j++) {
				Polygon waterbodyPolygon = _allWaterbodyPolygons.get(j);
				String waterbodyPolygonName = _allWaterbodyPolygonNames.get(j);
				long waterbodyObjectID = _allWaterbodyPolygonObjectIDs.get(j);
				if(centerlineGeometry.intersectsPolygon(waterbodyPolygon)) {
					//    				System.out.println("Intersection found: centerline, waterbodyObjectID, polygonName="+centerlineName+":"+ waterbodyObjectID+":"+ waterbodyPolygonName);
					if(channelPolygon==null) {
						channelPolygon = waterbodyPolygon;
					}else {
						//    					System.out.println("channelPolygon type:"+channelPolygon.getClass().getName());
						//    					System.out.println("waterbodyPolygon type: "+waterbodyPolygon.getClass().getName());
						//    					System.out.println("do they intersect? "+channelPolygon.intersects(waterbodyPolygon));
						//if the two polygons don't intersect, the result will be a MultiPolygon rather than a Polygon!
						channelPolygon = channelPolygon.union(waterbodyPolygon);
					}
				}
			}

			//now cut the polygon into pieces using the centerline endpoints.

			//what is the line: maybe perpendicular line to last centerline segment
			//If two channel node: make the angle the average of the connecting line segments angles
			//otherwise, make it perpendicular to the line segment.

			//for now, just to see if it will work, make it perpendicular to last line segment, and passing through endpoint.
			Centerline centerline = network.getCenterline(centerlineName);
			int numCenterlinePoints = centerline.getNumCenterlinePoints();
			int secondPointIndex = 1;
			int secondToLastPointIndex = numCenterlinePoints-2;
			CenterlinePoint upstreamPoint = centerline.getCenterlinePoint(0);
			CenterlinePoint secondPoint = centerline.getCenterlinePoint(secondPointIndex);

			CenterlinePoint downstreamPoint = centerline.getCenterlinePoint(numCenterlinePoints-1);
			CenterlinePoint secondToLastPoint = null;
			if(numCenterlinePoints < 2) {
				System.out.println("error: number of centerline points < 2");
			}else {
				secondToLastPoint = centerline.getCenterlinePoint(secondToLastPointIndex);
			}

			//create perpendicular line segment that is very long and centered on centerline endpoint.
			//the polygon should extend a little beyond the dsm2 channel (distance is determined by 
			//outsideChannelOffset) in the longitudinal direction at the upstream and downstream ends. 
			//After splitting the polygon, this will prevent external polygons from intersecting with the 
			//centerlineGeometry object.
			double outsideChannelOffset = .001;

			double upstreamX = CsdpFunctions.feetToMeters(upstreamPoint.getXFeet());
			double upstreamY = CsdpFunctions.feetToMeters(upstreamPoint.getYFeet());
			double secondPointX = CsdpFunctions.feetToMeters(secondPoint.getXFeet());
			double secondPointY = CsdpFunctions.feetToMeters(secondPoint.getYFeet());
			double secondToLastPointX = CsdpFunctions.feetToMeters(secondToLastPoint.getXFeet());
			double secondToLastPointY = CsdpFunctions.feetToMeters(secondToLastPoint.getYFeet());
			double downstreamX = CsdpFunctions.feetToMeters(downstreamPoint.getXFeet());
			double downstreamY = CsdpFunctions.feetToMeters(downstreamPoint.getYFeet());

			//it is possible for a centerline to have duplicate points at the upstream or downstream ends. Correct this
			//by changing the secondPoint and secondToLastPointValues.
			while(upstreamX==secondPointX && upstreamY==secondPointY) {
				if(numCenterlinePoints < secondPointIndex+2) {
					System.out.println("Error processing centerline "+centerlineName+": trying to increment secondPointIndex,"
							+ "but centerline does not have enough points. exiting.");
					System.exit(0);
				}else {
					secondPointIndex++;
					secondPoint = centerline.getCenterlinePoint(secondPointIndex);
					secondPointX = CsdpFunctions.feetToMeters(secondPoint.getXFeet());
					secondPointY = CsdpFunctions.feetToMeters(secondPoint.getYFeet());
				}
			}
			while(downstreamX==secondToLastPointX && downstreamY==secondToLastPointY) {
				if(secondToLastPointIndex==0) {
					System.out.println("Error processing centerline "+centerlineName+": trying to decrement "
							+ "secondToLastPointIndex, but I'm already at the upstream end. exiting");
					System.exit(0);
				}else {
					secondToLastPointIndex--;
					secondToLastPoint = centerline.getCenterlinePoint(secondToLastPointIndex);
					secondToLastPointX = CsdpFunctions.feetToMeters(secondToLastPoint.getXFeet());
					secondToLastPointY = CsdpFunctions.feetToMeters(secondToLastPoint.getYFeet());
				}
			}
			
			//theta2 is theta + 90 degrees
			//theta1Upstream is the angle (radians) of the upstream line segment
			//theta2Upstream is the angle (radians) of the very long perpendicular line segment that sits just outside
			//of the upstream end of the centerline, and is used for splitting polygons.
			double theta1Upstream = CsdpFunctions.getTheta(secondPointX, upstreamX, secondPointY, upstreamY);
			double theta1Downstream = CsdpFunctions.getTheta(secondToLastPointX, downstreamX, secondToLastPointY, downstreamY);
			double theta2Upstream = Math.PI/2.0 + theta1Upstream;
			double theta2Downstream = Math.PI/2.0 + theta1Downstream;

			//adjust the upstream and downstream centerline endpoints by moving them slightly outside the the centerline 
			//in the longitudinal direction
			double upstreamXOutside = upstreamX+outsideChannelOffset*Math.cos(theta1Upstream);
			double upstreamYOutside = upstreamY+outsideChannelOffset*Math.sin(theta1Upstream);
			double downstreamXOutside = downstreamX+outsideChannelOffset*Math.cos(theta1Downstream);
			double downstreamYOutside = downstreamY+outsideChannelOffset*Math.sin(theta1Downstream);

			//now find the endpoints of the line segments that will be used for splitting 
			double upstreamEndpoint1X = upstreamXOutside + PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.cos(theta2Upstream);
			double upstreamEndpoint1Y = upstreamYOutside + PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.sin(theta2Upstream);
			double upstreamEndpoint2X = upstreamXOutside - PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.cos(theta2Upstream);
			double upstreamEndpoint2Y = upstreamYOutside - PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.sin(theta2Upstream);

			double downstreamEndpoint1X = downstreamXOutside + PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.cos(theta2Downstream);
			double downstreamEndpoint1Y = downstreamYOutside + PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.sin(theta2Downstream);
			double downstreamEndpoint2X = downstreamXOutside - PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.cos(theta2Downstream);
			double downstreamEndpoint2Y = downstreamYOutside - PERPENDICULAR_LINE_SEGMENT_HALF_LENGTH * Math.sin(theta2Downstream);

			Coordinate[] upstreamLineSegmentCoordinates = new Coordinate[2];
			Coordinate[] downstreamLineSegmentCoordinates = new Coordinate[2];
			upstreamLineSegmentCoordinates[0] = new Coordinate(upstreamEndpoint1X, upstreamEndpoint1Y);
			upstreamLineSegmentCoordinates[1] = new Coordinate(upstreamEndpoint2X, upstreamEndpoint2Y);
			downstreamLineSegmentCoordinates[0] = new Coordinate(downstreamEndpoint1X, downstreamEndpoint1Y);
			downstreamLineSegmentCoordinates[1] = new Coordinate(downstreamEndpoint2X, downstreamEndpoint2Y);

			//These are the lines that are
			//  1. perpendicular to the centerline line segments at the upstream and downstream ends
			//  2. pass through a point that is slightly outside the centerline in the longitudinal direction
			//  3. used to split the waterbody polygons
			LineString upstreamBoundaryLine =  _geometryFactory.createLineString(upstreamLineSegmentCoordinates);
			LineString downstreamBoundaryLine =  _geometryFactory.createLineString(downstreamLineSegmentCoordinates);

			//now split the channel polygon on the upstream line, keep the part that intersects the centerline,
			//then repeat for the downstream line
			channelPolygon = splitAndMergePolygons(centerlineName, channelPolygon, centerlineGeometry, upstreamBoundaryLine);
			channelPolygon = splitAndMergePolygons(centerlineName, channelPolygon, centerlineGeometry, downstreamBoundaryLine);
			//now expand the polygon using the centerline as a reference
			System.out.println("DON'T FORGET TO EXPAND POLYGONS...");

			perpendicularLineHashtable.put(centerlineName+"_uPerp", upstreamBoundaryLine);
			perpendicularLineHashtable.put(centerlineName+"_dPerp", downstreamBoundaryLine);

			if(channelPolygon!=null) {
				modelGridChannelPolygons.put(centerlineName, channelPolygon);
			}else {
				System.out.println("no waterbody polygons are available for centerline "+centerlineName);
			}
		}//for each centerline

		return modelGridChannelPolygons;
	}//createModelGridChannelPolygons

	/*
	 * Given
	 * channelPolygon: the union of all polygons that intersect the CenterlineGeometry object.
	 * centerlineGeometry: the CenterlineGeometry object that contains the centerline information
	 * boundaryLine: A line which is perpendicular to a centerline endpoint and passes through the endpoint. Used for splitting.
	 * 
	 * Split the channelPolygon using the boundary line. Return a Polygon which is the union of all polygons that intersect 
	 * the centerlineGeometry object
	 * 
	 * Will return a Polygon or a MultiPolygon. When you union two polygons that don't intersect, you get a MultiPolygon
	 */
	private Geometry splitAndMergePolygons(String centerlineName, Geometry channelPolygon, CenterlineGeometry centerlineGeometry, Geometry boundaryLine) {
		Geometry newChannelPolygon = null;
		if(channelPolygon!=null) {
			//	    	System.out.println("splitAndMergePolygons: initial polygon area="+channelPolygon.getArea());
			GeometryCollection splitPolygonGeometryCollection = PolygonTools.splitPolygon(channelPolygon, boundaryLine);
			newChannelPolygon = null;

			//			this isn't looping through all geometriesl. numGeometries is 1
			//			splitPolygonGeometryCollection.getFactory().toGeometryArray(splitPolygonGeometryCollection);
			//			System.out.println("numGeometries="+splitPolygonGeometryCollection.getNumGeometries());
			//			System.out.println("geometry 0="+splitPolygonGeometryCollection.getGeometryN(0));
			//			System.out.println("geometry 1="+splitPolygonGeometryCollection.getGeometryN(1));

			Hashtable<String, Geometry> featuresToWrite = new Hashtable<String, Geometry>();
			for(int j=0; j<splitPolygonGeometryCollection.getNumGeometries(); j++) {
				Polygon polygon = (Polygon) splitPolygonGeometryCollection.getGeometryN(j);
				if(centerlineGeometry.intersectsPolygon(polygon)) {
					featuresToWrite.put(centerlineName+"_"+j+"_int", polygon);
					//    				System.out.println("splitAndMergePolygons: intersection found. newChannelPolygon="+newChannelPolygon);

					if(newChannelPolygon==null) {
						newChannelPolygon = (Geometry) polygon.clone();
					}else {
						newChannelPolygon = newChannelPolygon.union(polygon);
					}
				}else {
					featuresToWrite.put(centerlineName+"_"+j+"_noint", polygon);
				}
				//	    	System.out.println("splitAndMergePolygons: initial, return polygon areas="+channelPolygon.getArea()+","+newChannelPolygon.getArea());
			}
			if(DEBUG ) {
				if(centerlineName.equals("309") || centerlineName.equals("310") || centerlineName.equals("434") || 
						centerlineName.equals("24") || centerlineName.equals("6")) {
					writeFeaturesToNetworkFile(featuresToWrite, "d:/dsm2gisreference/csdpfiles/final/gisVolTesting", centerlineName+"_resultsOfPolygonSplits", null);
				}
			}
		}
		return newChannelPolygon;
		}//splitAndMergePolygons

		private void writeNetworkWithPerpendicularLines(Network network, Hashtable<String, LineString> perpendicularLineHashtable, 
				String directory, String filename) {
			Enumeration<String> pEnumeration = perpendicularLineHashtable.keys();
			while(pEnumeration.hasMoreElements()) {
				String oldCenterlineName = new String(pEnumeration.nextElement());
				String newCenterlineName = oldCenterlineName.replaceAll("\\s+", "_");
				LineString perpendicularLine = perpendicularLineHashtable.get(oldCenterlineName);
				Coordinate[] coordinates = perpendicularLine.getCoordinates();
				network.addCenterline(newCenterlineName);
				Centerline centerline = network.getCenterline(newCenterlineName);
				for(int i=0; i<coordinates.length; i++) {
					Coordinate coordinate = coordinates[i];
					centerline.addCenterlinePointFeet(CsdpFunctions.metersToFeet(coordinate.x), CsdpFunctions.metersToFeet(coordinate.y));
				}
			}
			network.sortCenterlineNames();
			NetworkOutput noutput = NetworkOutput.getInstance(directory, filename,
					CsdpFunctions.getNetworkFiletype(), network, null);
			noutput.writeData();

		}

		/*
		 * Write polygon features to a CSDP network file, for diagnostic purposes only.
		 * centerlineNames to export: if not null, should be a list of centerline names which are the only names to be exported.
		 */
		private void writeFeaturesToNetworkFile(Hashtable<String, Geometry> modelGridFeatures, String directory, String filename, 
				Vector<String> centerlineNamesToExport) {
			CsdpFrame csdpFrame = null;
			Network network = new Network(filename, csdpFrame);
			Enumeration<String> keys = modelGridFeatures.keys(); 
			while(keys.hasMoreElements()) {
				String oldCenterlineName = new String(keys.nextElement());
				String newCenterlineName = oldCenterlineName.replaceAll("\\s+", "_");

				Geometry channelPolygon = modelGridFeatures.get(oldCenterlineName);
				//channelPolygon could be a Polygon or a MultiPolygon
				if(channelPolygon instanceof Polygon && channelPolygon instanceof MultiPolygon) {
					System.out.println("instance of both!!!");
					System.exit(0);
				}
				if(centerlineNamesToExport==null || centerlineNamesToExport.contains(newCenterlineName)) {
					if(channelPolygon instanceof Polygon) {
						network.addCenterline(newCenterlineName);
						Centerline centerline = network.getCenterline(newCenterlineName);
						Coordinate[] coordinates = channelPolygon.getCoordinates();
						if(coordinates.length>0) {
							for(int j=0; j<coordinates.length; j++) {
								Coordinate coordinate = coordinates[j];
								centerline.addCenterlinePointFeet(CsdpFunctions.metersToFeet(coordinate.x), CsdpFunctions.metersToFeet(coordinate.y));
							}
						}else {
							network.removeCenterline(newCenterlineName);
						}
					}else if(channelPolygon instanceof MultiPolygon) {
						for(int j=0; j<channelPolygon.getNumGeometries(); j++) {
							String newCenterlineNameWithIndex = newCenterlineName+"_"+j;
							network.addCenterline(newCenterlineNameWithIndex);
							Centerline newCenterline = network.getCenterline(newCenterlineNameWithIndex);
							Polygon polygonN = (Polygon) channelPolygon.getGeometryN(j);
							Coordinate[] coordinates =  polygonN.getCoordinates();
							if(coordinates.length>0) {
								for(int k=0; k<coordinates.length; k++) {
									newCenterline.addCenterlinePointFeet(CsdpFunctions.metersToFeet(coordinates[k].x), CsdpFunctions.metersToFeet(coordinates[k].y));
								}
							}else {
								network.removeCenterline(newCenterlineNameWithIndex);
							}
						}
					}
				}
			}

			network.sortCenterlineNames();
			NetworkOutput noutput = NetworkOutput.getInstance(directory, filename,
					CsdpFunctions.getNetworkFiletype(), network, null);
			noutput.writeData();
		}//writeFeaturesToNetworkFile

		/*
		 * Most of this code copied from http://docs.geotools.org/stable/tutorials/feature/csv2shp.html
		 */
		private void writeFeaturesToShapefile(File shapefile, Hashtable<String, Geometry> modelGridFeatures) {

			System.out.println("CreateDSM2ChanPolygons.writeFeaturesToShapefile");
			/*
			 * We use the DataUtilities class to create a FeatureType that will describe the data in our
			 * shapefile.
			 * 
			 * See also the createFeatureType method below for another, more flexible approach.
			 */
			//        SimpleFeatureType TYPE = null;
			//		try {
			//
			////			Change this to polygon type
			//			
			////			TYPE = DataUtilities.createType("Location",
			////			        "location:Point:srid=4326," + // <- the geometry attribute: Point type
			////			                "name:String," + // <- a String attribute
			////			                "number:Integer" // a number attribute
			////			);
			//			TYPE = DataUtilities.createType("Location",
			//			        "location:Point:srid=4326," + // <- the geometry attribute: Point type
			//			                "name:String," + // <- a String attribute
			//			                "number:Integer" // a number attribute
			//			);
			//		} catch (SchemaException e) {
			//			System.out.println("CreateDSM2ChanPolygons.writeFeaturesToShapefile: ShemaException caught");
			//			e.printStackTrace();
			//		}

			//		SimpleFeatureCollection collection = FeatureCollections.newCollection();
			DefaultFeatureCollection defaultFeatureCollection = new DefaultFeatureCollection();
			SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(_polygonShapefileFeatureType);
			System.out.println("allCenterlineNames.size="+_allCenterlineNames.size());
			for(int i=0; i<_allCenterlineNames.size(); i++) {
				String centerlineName = _allCenterlineNames.get(i);
				Geometry channelPolygon = modelGridFeatures.get(centerlineName);
				featureBuilder.add(channelPolygon);
				//    		System.out.println("writeFeaturesToShapefile: channelPolygon"+channelPolygon);
				//    		SimpleFeature feature = featureBuilder.buildFeature(null);
				SimpleFeature feature = featureBuilder.buildFeature(centerlineName);
				//    		System.out.println("Adding Polygon object to feature builder for centerlineName="+centerlineName);
				defaultFeatureCollection.add(feature);
			}


			ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();
			Map<String, Serializable> params = new HashMap<String, Serializable>();
			try {
				params.put("url", shapefile.toURI().toURL());
			} catch (MalformedURLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			params.put("create spatial index", Boolean.TRUE);

			ShapefileDataStore newDataStore = null;
			SimpleFeatureSource featureSource = null;
			Transaction transaction = new DefaultTransaction("create");
			String typeName = null;
			try {
				newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);
				newDataStore.createSchema(_polygonShapefileFeatureType);
				/*
				 * You can comment out this line if you are using the createFeatureType method (at end of
				 * class file) rather than DataUtilities.createType
				 */
				newDataStore.forceSchemaCRS(DefaultGeographicCRS.WGS84);
				typeName = newDataStore.getTypeNames()[0];
				featureSource = newDataStore.getFeatureSource(typeName);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			if (featureSource instanceof SimpleFeatureStore) {
				SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;

				featureStore.setTransaction(transaction);
				try {
					featureStore.addFeatures(defaultFeatureCollection);
					transaction.commit();

				} catch (Exception problem) {
					problem.printStackTrace();
					try {
						transaction.rollback();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}

				} finally {
					try {
						transaction.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				//            System.exit(0); // success!
			} else {
				System.out.println(typeName + " does not support read/write access");
				//            System.exit(1);
			}
		}//writeFeaturesToShapefile

		/**
		 * Copied from http://docs.geotools.org/stable/tutorials/feature/csv2shp.html
		 * Prompt the user for the name and path to use for the output shapefile
		 * 
		 * @param csvFile
		 *            the input csv file used to create a default shapefile name
		 * 
		 * @return name and path for the shapefile as a new File object
		 */
		private static File getNewShapeFile(File csvFile) {
			String path = csvFile.getAbsolutePath();
			String newPath = path.substring(0, path.length() - 4) + ".shp";

			JFileDataStoreChooser chooser = new JFileDataStoreChooser("shp");
			chooser.setDialogTitle("Save shapefile");
			chooser.setSelectedFile(new File(newPath));

			int returnVal = chooser.showSaveDialog(null);

			if (returnVal != JFileDataStoreChooser.APPROVE_OPTION) {
				// the user cancelled the dialog
				System.exit(0);
			}

			File newFile = chooser.getSelectedFile();
			if (newFile.equals(csvFile)) {
				System.out.println("Error: cannot replace " + csvFile);
				System.exit(0);
			}

			return newFile;
		}

	}//class CreateDSM2ChanPolygons


	//probably not needed any more. uses a different approach, which involves using levee centerlines. This approach was abandoned
	//because it seemed easier to work with waterbodies, because a waterbody in the channel will intersect CSDP centerline objects,
	//whereas levee centerlines will not. 


	//    /**
	//     * Stores a Vector of Line2D objects. Used to determine closest multilinestring to a given point
	//     * @author btom
	//     *
	//     */
	//    class MultiLine2D{
	//    	Vector<Line2D.Double> _allLine2D;
	//    	public MultiLine2D(MultiLineString multiLineString) {
	//        	_allLine2D = new Vector<Line2D.Double>();
	//        	Coordinate[] coordinates = null;
	//        	try {
	//        		coordinates = multiLineString.getCoordinates();
	//        	}catch(NullPointerException npe) {
	//        		System.out.println("NullPointerException caught while calling getCoordinates method, but this doesn't seem to be a problem");
	//        	}
	//        	for(int j=1; j<coordinates.length; j++) {
	//        		Coordinate c1 = coordinates[j-1];
	//        		Coordinate c2 = coordinates[j];
	//        		Line2D.Double line2D = new Line2D.Double(c1.x, c1.y, c2.x, c2.y);
	//        		_allLine2D.addElement(line2D);
	//        	}
	//    	}//constructor
	//
	//    	public double ptLineDist(double x, double y) {
	//    		double returnValue = Double.MAX_VALUE;
	//    		for(int i=0; i<_allLine2D.size(); i++) {
	//    			returnValue = Math.min(returnValue, _allLine2D.get(i).ptLineDist(x, y));
	//    		}
	//    		return returnValue;
	//    	}
	//    }//class MultiLine2D
	//    
	//    private MultiLine2D[] findClosestTwoMultiLine2DObjects(double x, double y) {
	//    	double[] minDistances = new double[]{Double.MAX_VALUE, Double.MAX_VALUE};
	//    	MultiLine2D[] closestTwoMultiLine2DObjects = null;
	//    	for(int i=0; i<_allLeveeCenterlines.size(); i++) {
	//    		MultiLine2D multiLine2D = _allLeveeCenterlines.get(i);
	//    		double lineDist = multiLine2D.ptLineDist(x, y);
	//    		double diff1 = lineDist-minDistances[0];
	//    		double diff2 = lineDist-minDistances[1];
	//    		if(diff1 < 0 || diff2 < 0) {
	//    			if(diff1<=diff2) {
	//    				minDistances[1]=lineDist;
	//    				closestTwoMultiLine2DObjects[1]=multiLine2D;
	//    			}else {
	//    				minDistances[0]=lineDist;
	//    				closestTwoMultiLine2DObjects[0]=multiLine2D;
	//    			}
	//    		}
	//    	}
	//    	return closestTwoMultiLine2DObjects;
	//    }//findClosestTwoMultiLine2DObjects
	//    
	//    private void createLeveeCenterlineVector() {
	//        File file = JFileDataStoreChooser.showOpenFile("shp", null);
	//        if (file == null) {
	//            return;
	//        }
	//        try {
	//	        FileDataStore store = FileDataStoreFinder.getDataStore(file);
	//	        SimpleFeatureSource featureSource = store.getFeatureSource();
	//	
	//	        // Create a map content and add our shapefile to it
	//	        MapContent map = new MapContent();
	//	        map.setTitle("Quickstart");
	//	
	//	        Style style = SLD.createSimpleStyle(featureSource.getSchema());
	//	        Layer layer = new FeatureLayer(featureSource, style);
	//	        FeatureSource layerFeatureSource = layer.getFeatureSource();
	//	        FeatureCollection layerFeatureCollection = layerFeatureSource.getFeatures();
	//	        FeatureIterator<SimpleFeatureImpl> layerFeatureIterator = layerFeatureCollection.features();
	//	        while(layerFeatureIterator.hasNext()) {
	//	        	SimpleFeatureImpl simpleFeatureImpl = layerFeatureIterator.next();
	////	        	FeatureId id = simpleFeatureImpl.getIdentifier();
	//	//        	FeatureType ft = simpleFeatureImpl.getType();
	//	//        	System.out.println("feature id, type="+id+","+ft);
	//	        	
	//	        	MultiLineString multiLineString = (MultiLineString) simpleFeatureImpl.getAttribute(0);
	//	        	_allLeveeCenterlines.addElement(new MultiLine2D(multiLineString));
	//	        }
	//        }catch(IOException e) {
	//        	System.out.println("IOException caught!!!");
	//        }
	//    }
	//
	//    
	//}//class CreateDSM2ChanPolygons

	//	public static void main(String[] args) {
	//		File file = new File("d:/dsm2GisReference/gisChannelVolumes/deltaLeveesCenterlines.shp");
	//		Map map = new HashMap();
	//		AsciiFileWriter afw = new AsciiFileWriter("d:/dsm2GisReference/gisChannelVolumes/data.txt");
	//		try {
	//			map.put( "url", file.toURL() );
	//			DataStore dataStore;
	//			dataStore = DataStoreFinder.getDataStore(map);
	//	
	//			String[] typeNames = dataStore.getTypeNames();
	//			afw.writeLine("typeNames.length="+typeNames.length);
	//			for(int i=0; i<typeNames.length; i++) {
	//				String typeName = typeNames[i];
	//				SimpleFeatureSource featureSource = dataStore.getFeatureSource( typeName );
	//				SimpleFeatureCollection collection = featureSource.getFeatures();
	//		
	//				ReferencedEnvelope env = collection.getBounds();
	//				double left = env.getMinX();
	//				double right = env.getMaxX();
	//				double top = env.getMaxY();
	//				double bottom = env.getMinY();
	//				afw.writeLine("typeName, featureSource, collection="+typeName+","+featureSource+","+collection);
	//				SimpleFeatureIterator featureSourceCollectionIterator = collection.features();
	//				int featureIndex = 0;
	//				while(featureSourceCollectionIterator.hasNext()) {
	//					SimpleFeature simpleFeature = featureSourceCollectionIterator.next();
	//					Collection values = simpleFeature.getValue();
	//					Iterator valuesIterator = values.iterator();
	//					while(valuesIterator.hasNext()) {
	//						Object valueElement = values.iterator().next();
	//						System.out.println("valueElement="+valueElement);
	//					}
	//					afw.writeLine("feature "+ featureIndex + "=" + simpleFeature);
	//					System.out.println("featureType="+simpleFeature.getFeatureType());
	//					GeometryAttribute geometryAttribute = simpleFeature.getDefaultGeometryProperty();
	//					System.out.println("geometryAttribute descriptor="+geometryAttribute.getDescriptor());
	//					//					List featureSourceFeature = simpleFeature.getAttributes();
	////					Iterator featureSourceFeatureIterator = featureSourceFeature.iterator();
	////					afw.writeLine("=====================================BEGIN=================================");
	////					int feature2Index = 0;
	////					while(featureSourceFeatureIterator.hasNext()) {
	////						Object o = featureSourceFeatureIterator.next();
	////						if(o==null) {
	////							afw.writeLine("null");
	////						}else {
	////							String s = o.toString();
	////							afw.writeLine("("+feature2Index+"):"+","+s);
	////						}
	////						feature2Index++;
	////					}
	//					afw.writeLine("=====================================END===================================");
	//					featureIndex++;
	//				}
	//				
	//			}
	//		} catch (Exception e) {
	//			System.out.println("Error");
	//			// TODO Auto-generated catch block
	//			e.printStackTrace();
	//			
	//		} finally {
	//			afw.close();
	//		}
	//		System.out.println("Done");
	//	}
	//}


