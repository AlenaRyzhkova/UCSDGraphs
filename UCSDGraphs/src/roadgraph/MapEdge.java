package roadgraph;

import geography.GeographicPoint;

/**
 * 
 * @author Alena Ryzhkova
 *
 */
public class MapEdge {
	private MapNode start;
	private MapNode end;
	
	private String roadName;
	private String roadType;
	private double lenght;
	
	public MapEdge(MapNode from, MapNode to, String roadName,
			String roadType, double length) {
		this.start = from;
		this.end = to;
		this.roadName = roadName;
		this.roadType = roadType;
		this.lenght = length;
	}

	public String getRoadName() {
		return roadName;
	}

	public void setRoadName(String roadName) {
		this.roadName = roadName;
	}

	public String getRoadType() {
		return roadType;
	}

	public void setRoadType(String roadType) {
		this.roadType = roadType;
	}

	public GeographicPoint getStart() {
		return start.getVertex();
	}

	public GeographicPoint getEnd() {
		return end.getVertex();
	}
	
	public MapNode getEndNode() {
		return end;
	}
	
	public MapNode getStartNode() {
		return start;
	}

	public double getLenght() {
		return lenght;
	}
	
	
}
