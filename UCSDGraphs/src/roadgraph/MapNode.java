package roadgraph;

import java.util.ArrayList;
import java.util.List;

import geography.GeographicPoint;

/**
 * 
 * @author Alena Ryzhkova
 *
 */
public class MapNode implements Comparable<MapNode>{
	private GeographicPoint vertex;
	
	private Double distance = Double.MAX_VALUE;
	public Double getHeuristicDistance() {
		return heuristicDistance;
	}

	public void setHeuristicDistance(Double heuristicDistance) {
		this.heuristicDistance = heuristicDistance;
	}

	private Double heuristicDistance = Double.MAX_VALUE;
	
	private List<MapEdge> edges; 	// List of outgoing neighbors;
	
	public MapNode(GeographicPoint vertex) {
		this.vertex = vertex;
		edges = new ArrayList<>();
	}
	
	public Double getDistance() {
		return distance;
	}

	public void setDistance(Double distance) {
		this.distance = distance;
	}

	public boolean addEdge(MapEdge edge) {
		if(edge.getStart().equals(vertex)) {
			edges.add(edge);
			return true;
		}
		return false;
	}

	public GeographicPoint getVertex() {
		return vertex;
	}

	public List<MapEdge> getEdges() {
		return edges;
	}

	@Override
	public int compareTo(MapNode another) {
		return distance.compareTo(another.distance);
	}

}
