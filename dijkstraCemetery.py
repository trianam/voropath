def _trijkstra(self, startA, endA):
    '''
    Alternative Dijkstra algorithm for calculating
    shortest path in a graph with some subpath
    not availables. In this case triplets that
    hit an obstacle.
    '''
    start = tuple(startA)
    end = tuple(endA)
    endTriplet = (end,end,end) #special triplet for termination
    inf = float("inf")
    path = []
    Q = []
    dist = {}
    prev = {}

    #create triplets
    for node0 in self._graph.nodes():
        for node1 in self._graph.neighbors(node0):
            for node2 in filter(lambda node: node!=node0, self._graph.neighbors(node1)):
                if not self._triangleIntersectPolyhedrons(np.array(node0), np.array(node1), np.array(node2)):
                    triplet = (node0,node1,node2)
                    if node0 != start:
                        dist[triplet] = inf
                    else:
                        dist[triplet] = 0
                    Q[:0] = [triplet]

    Q[:0] = [endTriplet]
    dist[endTriplet] = inf

    while Q:
        Q = sorted(Q, key=lambda el: dist[el])
        u = Q[0]
        Q = Q[1:]

        if u == endTriplet or dist[u] == inf:
            break

        for v in filter(lambda tri: u[1] == tri[0] and u[2] == tri[1], Q):
            tmpDist = dist[u] + self._graph[u[0]][u[1]]['weight']
            if tmpDist < dist[v]:
                dist[v] = tmpDist
                prev[v] = u

        if u[2] == end:
            tmpDist = dist[u] + self._graph[u[0]][u[1]]['weight'] + self._graph[u[1]][u[2]]['weight']
            if tmpDist < dist[endTriplet]:
                dist[endTriplet] = tmpDist
                prev[endTriplet] = u

    u = endTriplet
    while u in prev:
        u = prev[u]
        path[:0] = [u[1]]

    if path:
        path[len(path):] = [end]
        path[:0] = [start]

    return np.array(path)


def _dijkstraPlus(self, startA, endA):
    '''
    Not working, first tentative of alternative
    Dijkstra algorithm
    '''
    
    start = tuple(startA)
    end = tuple(endA)
    inf = float("inf")
    path = []
    Q = []
    dist = {}
    prev = {}

    dist[start] = 0
    for node in self._graph.nodes():
        Q[:0] = [node]
        if node != start:
            dist[node] = inf

    while Q:
        Q = sorted(Q, key=lambda el: dist[el])
        u = Q[0]
        Q = Q[1:]

        if u == end or dist[u] == inf:
            break

        for v in self._graph.neighbors(u):
            if v in Q:
                if (u == start) or (not self._triangleIntersectPolyhedrons(np.array(prev[u]), np.array(u), np.array(v))):
                    alt = dist[u] + self._graph.edge[u][v]['weight']
                    if alt < dist[v]:
                        dist[v] = alt
                        prev[v] = u

    u = end
    while u in prev:
        path[:0] = [u]
        u = prev[u]

    if path:
        path[:0] = [u]

    return np.array(path)
