package gov.ca.water.gisutil;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.util.LineStringExtracter;
import com.vividsolutions.jts.operation.polygonize.Polygonizer;


/**
 * copied from https://gis.stackexchange.com/questions/189976/jts-split-arbitrary-polygon-by-a-line
 * @author btom
 *
 */
public class PolygonTools {

	public static Geometry polygonize(Geometry geometry) {
		List lines = LineStringExtracter.getLines(geometry);
		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(lines);
		Collection polys = polygonizer.getPolygons();
		Polygon[] polyArray = GeometryFactory.toPolygonArray(polys);
		return geometry.getFactory().createGeometryCollection(polyArray);
	}

	public static GeometryCollection splitPolygon(Geometry poly, Geometry line) {
		Geometry nodedLinework = poly.getBoundary().union(line);
		Geometry polys = polygonize(nodedLinework);

		
//		this seems to be the problem. only returning one polygon 
		
		// Only keep polygons which are inside the input

		List output = new ArrayList();
		for (int i = 0; i < polys.getNumGeometries(); i++) {
			Polygon candpoly = (Polygon) polys.getGeometryN(i);
			if (poly.contains(candpoly.getInteriorPoint())) {
				output.add(candpoly);
			}
		}
		return poly.getFactory().createGeometryCollection(GeometryFactory.toGeometryArray(output));
	}
} 