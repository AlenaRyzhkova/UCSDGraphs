package roadgraph;


import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;

/**
 * @author UCSD MOOC development team
 * @author Alena Ryzhkova
 * 
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
public class MapGraph {
	// AdjacencyList implementation of graph
	HashMap<GeographicPoint, MapNode> adjList;
	
	
	/** 
	 * Create a new empty MapGraph 
	 */
	public MapGraph()
	{
		adjList = new HashMap<>();
	}
	
	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
		return adjList.size();
	}
	
	/**
	 * Return location of the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
		return adjList.keySet();
	}
	
	private Collection<MapNode> getNodeVertices(){
		return adjList.values();
	}
	
	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
		int numEdges = 0;
		Set<GeographicPoint> vertecies = adjList.keySet();
		for(GeographicPoint v: vertecies) {
			numEdges += adjList.get(v).getEdges().size();
		}
		return numEdges;
	}

	
	
	/** Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does 
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location)
	{
		if(location == null || isVertex(location))
			return false;
		adjList.put(location, new MapNode(location));
		return true;
	}
	
	/**
	 * Return MapNode associated with given location (geographic point)
	 * @param point for searching among vertices
	 * @return MapNode for this Geographic point
	 * @throws IllegalArgumentException If the point have not been added as node in the MapGraph
	 */
	private MapNode getVertex(GeographicPoint loc) throws IllegalArgumentException{
		if(!isVertex(loc))
			throw new IllegalArgumentException("Location " + loc + " is not in the MapGraph");
		return adjList.get(loc);
	}
	
	private boolean isVertex(GeographicPoint loc) {
		if(getVertices().contains(loc))
			return true;
		else
			return false;
	}
	
	/**
	 * Adds a directed edge to the graph from pt1 to pt2.  
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * @param from The starting point of the edge
	 * @param to The ending point of the edge
	 * @param roadName The name of the road
	 * @param roadType The type of the road
	 * @param length The length of the road, in km
	 * @throws IllegalArgumentException If the points have not already been
	 *   added as nodes to the graph, if any of the arguments is null,
	 *   or if the length is less than 0.
	 */
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName,
			String roadType, double length) throws IllegalArgumentException {
		
		// check that the input is correct and complete
		if(roadName==null)
			throw new IllegalArgumentException("Not specified roadName");
		if(roadType==null)
			throw new IllegalArgumentException("Not specified roadType");
		if(length<=0)
			throw new IllegalArgumentException("Invalid length of edge");
		
		MapNode start = getVertex(from);
		MapNode end = getVertex(to);
		
		start.addEdge(new MapEdge(start, end, roadName, roadType, length));
		
	}
	

	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return bfs(start, goal, temp);
	}
	
	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, 
			 					     GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		//Check that start and end point in the MapGraph
		if(!isVertex(start)||!isVertex(goal)) {
			System.out.println("Start or/and goal vertex isn't in the MapGraph");
			return null;
		}
		
		// Initialize everything
		MapNode startNode = getVertex(start);
		MapNode goalNode = getVertex(goal);
		HashMap<MapNode, MapNode> parent = new HashMap<>();
		
		boolean found = bfsSearch(startNode, goalNode, parent, nodeSearched);
		
		if(!found) {
			System.out.println("No path exist");
			return null;
		}
		
		return constructPath(startNode, goalNode, parent);
	}
	
	private static boolean bfsSearch(MapNode startNode, MapNode goalNode, HashMap<MapNode, MapNode> parent,
			Consumer<GeographicPoint> nodeSearched) {
		
		Queue<MapNode> queue = new LinkedList<>();
		HashSet<MapNode> visited = new HashSet<>();
		boolean found = false;
		
		queue.add(startNode);
		while(!queue.isEmpty()) {
			MapNode curr = queue.remove();
			if(curr.equals(goalNode)) {
				found = true;
				break;
			}
			List<MapEdge> neighbors = curr.getEdges();
			for(MapEdge n: neighbors) {
				MapNode neighborNode = n.getEndNode();
				if(!visited.contains(neighborNode)) {
					visited.add(neighborNode);
					parent.put(neighborNode, curr);
					queue.add(neighborNode);
					nodeSearched.accept(neighborNode.getVertex());
				}
			}
		}
		return found;
	}
	
	private static List<GeographicPoint> constructPath(MapNode startNode, MapNode goalNode, HashMap<MapNode, MapNode> parent){
		LinkedList<GeographicPoint> path = new LinkedList<>();
		MapNode curr = goalNode;
		while(curr!=startNode) {
			path.addFirst(curr.getVertex());
			curr=parent.get(curr);
		}
		path.addFirst(startNode.getVertex());
		return path;
	}

	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		// You do not need to change this method.
        Consumer<GeographicPoint> temp = (x) -> {};
        return dijkstra(start, goal, temp);
	}
	
	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, 
										  GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		//Check that start and end point in the MapGraph
		if(!isVertex(start)||!isVertex(goal)) {
			System.out.println("Start or/and goal vertex isn't in the MapGraph");
			return null;
		}
				
		// Initialize everything
		MapNode startNode = getVertex(start);
		MapNode goalNode = getVertex(goal);
		HashMap<MapNode, MapNode> parent = new HashMap<>();
		
		boolean found =  dijkstraSearch(startNode, goalNode, parent, nodeSearched);
		if(!found) {
			System.out.println("No path exist");
			return null;
		}
		return constructPath(startNode, goalNode, parent);
	} 
	
	private boolean dijkstraSearch(MapNode startNode, MapNode goalNode, HashMap<MapNode, MapNode> parent,
			Consumer<GeographicPoint> nodeSearched) {
		
		// Initialize priority queue and visited set
		PriorityQueue<MapNode> queue = new PriorityQueue<>();
		HashSet<MapNode> visited = new HashSet<>();
		boolean found = false;
		// Initialize distance to infinity
		for(MapNode node: getNodeVertices()) {
			node.setDistance(Double.MAX_VALUE);
		}
		
		startNode.setDistance(0.0);
		queue.add(startNode);
		
		while(!queue.isEmpty()) {
			MapNode curr = queue.remove();
			if(!visited.contains(curr)) {
				visited.add(curr);
				if(curr.equals(goalNode)) {
					found = true;
					break;
				}
				List<MapEdge> neighbors = curr.getEdges();
				for(MapEdge n: neighbors){
					double newPath = n.getStartNode().getDistance() + n.getLenght();
					if(newPath<n.getEndNode().getDistance()) {
						n.getEndNode().setDistance(newPath);
						parent.put(n.getEndNode(), curr);
						queue.add(n.getEndNode());
						nodeSearched.accept(n.getEndNode().getVertex());
					}
				}
			}
		}
		System.out.println("Dijkstra: Visited.size: " + visited.size());
		return found;
	}

	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return aStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, 
											 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		//Check that start and end point in the MapGraph
		if(!isVertex(start)||!isVertex(goal)) {
			System.out.println("Start or/and goal vertex isn't in the MapGraph");
			return null;
		}
						
		// Initialize everything
		MapNode startNode = getVertex(start);
		MapNode goalNode = getVertex(goal);
		HashMap<MapNode, MapNode> parent = new HashMap<>();
				
		boolean found =  aStarSearchEngine(startNode, goalNode, parent, nodeSearched);
		if(!found) {
			System.out.println("No path exist");
			return null;
		}
		return constructPath(startNode, goalNode, parent);
		
	}

	
	private boolean aStarSearchEngine(MapNode startNode, MapNode goalNode, HashMap<MapNode, MapNode> parent,
			Consumer<GeographicPoint> nodeSearched) {
		// Initialize priority queue and visited set
		Comparator<MapNode> comparator = new HeuristicDistComparator();
		PriorityQueue<MapNode> queue = new PriorityQueue<>(comparator);
		HashSet<MapNode> visited = new HashSet<>();
		boolean found = false;
		// Initialize distance to infinity
		for(MapNode node: getNodeVertices()) {
			node.setDistance(Double.MAX_VALUE);
			node.setHeuristicDistance(Double.MAX_VALUE);
		}
				
		startNode.setDistance(0.0);
		startNode.setHeuristicDistance(0.0);
		queue.add(startNode);
		while(!queue.isEmpty()) {
			MapNode curr = queue.remove();
			if(!visited.contains(curr)) {
				visited.add(curr);
				if(curr.equals(goalNode)) {
					found = true;
					break;
				}
				List<MapEdge> neighbors = curr.getEdges();
				for(MapEdge n: neighbors){
					double newPath = n.getStartNode().getDistance() + n.getLenght();
					double hPath =  n.getEnd().distance(goalNode.getVertex());
					if(newPath<n.getEndNode().getDistance()) {
						n.getEndNode().setDistance(newPath);
						n.getEndNode().setHeuristicDistance(hPath);
						parent.put(n.getEndNode(), curr);
						queue.add(n.getEndNode());
						nodeSearched.accept(n.getEndNode().getVertex());
					}
				}
			}
		}
		System.out.println(goalNode.getVertex() + ": " + goalNode.getDistance());
		//System.out.println("AStar Visited.size: " + visited.size());
		return found;
	}
	
	/**
	 * Calculate distance between start and goal locations, using A* algorithm
	 * 
	 * @param start location
	 * @param goal location
	 * @return distance between start location and goal location
	 * @throws IllegalArgumentException
	 */
	private double calculateDist(MapNode startNode, MapNode goalNode){
		HashMap<MapNode, MapNode> parent = new HashMap<>();
		Consumer<GeographicPoint> temp = (x) -> {};
		boolean found =  aStarSearchEngine(startNode, goalNode, parent,temp);
		//If we can't find path we return infinity distance (equal infinity cost travel between these vertices)
		if(!found)
			return Double.MAX_VALUE;
		return goalNode.getDistance();
		
	}
	
	private static void printTSPAdjMatrix(double[][] tspAddmatrix) {
		for(int i=0; i<tspAddmatrix.length; i++) {
			for(int j=0; j<tspAddmatrix[i].length; j++) {
				String e;
				if(tspAddmatrix[i][j]==Double.MAX_VALUE)
					System.out.print("inf ");
				else {
					System.out.printf("%.2f", tspAddmatrix[i][j]);
					System.out.print(" ");
				}
			}
			System.out.println();
		}
	}
	
	/**
	 * Solve Traveling Salesperson Problem using greedy algorithm.
	 * @param list of vertices
	 * @return route which visit each vertex and return back to start
	 */
	public List<GeographicPoint> tspGreedyAlg (List<GeographicPoint> vertices) throws IllegalArgumentException{
		
		// Initialization
		
		List<GeographicPoint> tspRoute = new ArrayList<>();
		
		// Create MapNodes from locations
		ArrayList<MapNode> verticesNodes = new ArrayList<>();
		for(GeographicPoint loc: vertices) {
			if(!isVertex(loc))
				throw new IllegalArgumentException();
			verticesNodes.add(getVertex(loc));
		}
		
		//Create adjacency matrix for traveling salesperson problem
		int numVert = vertices.size();
		double[][] tspAddmatrix = new double[numVert][numVert];
		for(int i=0;i<numVert; i++) {
			for(int j=0; j<numVert; j++)
				if(i==j)
					// Initialize the diagonal of the matrix by infinity
					// so ... when we look for a minimum, we don't choose the same vertex
					tspAddmatrix[i][j]=Double.MAX_VALUE;	
				else
					tspAddmatrix[i][j]=calculateDist(verticesNodes.get(i), verticesNodes.get(j));
		}
		
		System.out.println("Print tsp adjacency matrix");
		printTSPAdjMatrix(tspAddmatrix);
		
		// Engine of TSP Greedy algorithm
		
		// Zero step
		int i = 0;	// Choose starting node
		int next = 0;
		double minDist = Double.MAX_VALUE;
		
		// we start from zero node, so we mustn't choose it up to the end of the route
		// that's why we set infinity value to the zero column
		for(int k=0; k<numVert; k++)
			tspAddmatrix[k][i] = Double.MAX_VALUE;
		
		do {
			i = next;
			tspRoute.add(verticesNodes.get(i).getVertex());
			next = 0;
			minDist = Double.MAX_VALUE;
			// Look for node with minimum distance from current
			for(int j=0; j<numVert;j++) {
				if(minDist>tspAddmatrix[i][j]) {
					minDist=tspAddmatrix[i][j];
					next = j;
				}
			}
			
			//prepare Matrix for next iteration
			// we go to the "next" done, so we mustn't choose it in future, so we set infinity in "next"'s column
			for(int k=0; k<numVert; k++)
				tspAddmatrix[k][next] = Double.MAX_VALUE;
			System.out.println("\nCurrent vertex idx = " + i + ", next = " + next + ", minDist = " + minDist);
			printTSPAdjMatrix(tspAddmatrix);
			
		}while(minDist<Double.MAX_VALUE);
		// back to home
		tspRoute.add(verticesNodes.get(0).getVertex());
		
		return tspRoute;
	}
	
	public static void printRoute(List<GeographicPoint> route) {
		for(GeographicPoint loc: route)
			System.out.println("[" + loc + "]");
	}
	
	public static void main(String[] args)
	{
		/*System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("DONE.");
		
		// You can use this method for testing.  
		
		*/
		/* Here are some test cases you should try before you attempt 
		 * the Week 3 End of Week Quiz, EVEN IF you score 100% on the 
		 * programming assignment.
		 */
		
		/*MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);
		
		GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);
		
		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroutBFS = simpleTestMap.bfs(testStart, testEnd);
		System.out.println(testroutBFS);
		List<GeographicPoint> testroute = simpleTestMap.dijkstra(testStart,testEnd);
		List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart,testEnd);
		
		
		MapGraph testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		
		// A very simple test using real data
		testStart = new GeographicPoint(32.869423, -117.220917);
		testEnd = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		System.out.println(testroute);
		
		
		// A slightly more complex test using real data
		testStart = new GeographicPoint(32.8674388, -117.2190213);
		testEnd = new GeographicPoint(32.8697828, -117.2244506);
		System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		testroute = testMap.dijkstra(testStart,testEnd);
		testroute2 = testMap.aStarSearch(testStart,testEnd);
		*/
		
		
		/* Use this code in Week 3 End of Week Quiz */
		/*
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);
		
		
		List<GeographicPoint> route = theMap.dijkstra(start,end);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		*/
		
		/*Testing TSP greedy algorithm*/
		
		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		
		GraphLoader.loadRoadMap("data/maps/moscow_northwest.map", firstMap);
		System.out.println("DONE.\nCreate list of points for TSP algorithm");
		
		List<GeographicPoint> pointsForTSP = new ArrayList<>();
		pointsForTSP.add(new GeographicPoint(55.8508307, 37.3651173));
		pointsForTSP.add(new GeographicPoint(55.8378601, 37.3523172));
		pointsForTSP.add(new GeographicPoint(55.8335013, 37.3961817));
		pointsForTSP.add(new GeographicPoint(55.8271157, 37.4376596));
		pointsForTSP.add(new GeographicPoint(55.8454759, 37.3918207));
		
		printRoute(pointsForTSP);
		
		System.out.println("\nStart TSP Greedy algorithm\n");
		List<GeographicPoint> tspResult = firstMap.tspGreedyAlg(pointsForTSP);
		printRoute(tspResult);
	}
	
}
