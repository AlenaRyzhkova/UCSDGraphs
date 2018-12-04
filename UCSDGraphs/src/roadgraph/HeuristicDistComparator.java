package roadgraph;

import java.util.Comparator;

public class HeuristicDistComparator implements Comparator<MapNode> {

	@Override
	public int compare(MapNode node0, MapNode node1) {
		Double tdis0 = node0.getDistance() + node0.getHeuristicDistance();
		Double tdis1 = node1.getDistance() + node1.getHeuristicDistance();
		return tdis0.compareTo(tdis1);
	}

}
