#include <vector>
#include <set>
#include <map>
#include <iostream>
using std::vector;
using std::set;
using std::map;
using std::cout; using std::endl;

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/biconnected_components.hpp>

namespace boost {
  struct edge_component_t {
    typedef edge_property_tag kind;
  } edge_component;
}

using namespace boost;

template<typename Gr, typename Vert>
class myBFSvisitor : public default_bfs_visitor {
public:
    myBFSvisitor(vector<Vert> & mv, vector<Vert> & ut) { myverts = &mv; untouched = &ut; }
    void finish_vertex(Vert u, const Gr & g) { myverts->push_back(u); untouched->erase( find(untouched->begin(), untouched->end(), u) ); }
private :
    vector<Vert> * myverts, *untouched;
};

typedef adjacency_list < vecS, vecS, undirectedS, no_property, property<edge_component_t, std::size_t> > Graph;

void findBiconn(Graph & g, vector<vector<int> > & retcomps, vector<int> & retCollapseInto) {
    ///////////////////// find biconnected components
    property_map < Graph, edge_component_t >::type component = get(edge_component, g);
    std::size_t num_comps = biconnected_components(g, component);
    std::vector<int> art_points;
    articulation_points(g, std::back_inserter(art_points));
    cout << "found num-comps " << num_comps << " num-arts " << art_points.size() << endl;

    map<int, set<int> > biconns;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        int a = source(*ei,g), b = target(*ei,g), ci = component[*ei];
        if(biconns.find(ci) == biconns.end()) biconns[ci] = set<int> ();
        biconns[ci].insert(a); biconns[ci].insert(b);
    }
    vector< vector<int> > bcomps;
    for( map<int, set<int> >::iterator it = biconns.begin(); it != biconns.end(); it++ ) bcomps.push_back( vector<int>(it->second.begin(), it->second.end()) );
    for(int i=0; i < bcomps.size(); i++) {
        cout << i << " Biconnected component of size " << bcomps[i].size() << " : ";
        for(vector<int>::iterator sit = bcomps[i].begin(); sit != bcomps[i].end(); sit++) cout << *sit << " ";
        cout << endl;
    }

    ///////////////////// find order in which components will be collapsed
    map<int, set<int> > edges;
    for(int i=0; i < bcomps.size(); i++)
        for(int k=i+1; k < bcomps.size(); k++) {
            vector<int> intxn;
            std::back_insert_iterator<vector<int> > intit (intxn);
            set_intersection(
                    bcomps[i].begin(), bcomps[i].end(),
                    bcomps[k].begin(), bcomps[k].end(), intit);
            if(intxn.size() > 1) { cout << "not expecting intersection > 1" << endl; exit(0); }
            if(intxn.size() != 1) continue;
            if(edges.find(i) == edges.end()) edges[i] = set<int>();
            if(edges.find(k) == edges.end()) edges[k] = set<int>();
            edges[i].insert(k);
            edges[k].insert(i);
            cout << "anedge " << i <<" -- "<< k << endl;
        }
    set<int> untouched; for(int i=0; i < bcomps.size(); i++) untouched.insert(i);
    map<int,vector<int> > children;
    vector<int> order; //order.push_back(0);
    while(untouched.size() > 0) {
        int oi = -999;
        for(int i=0; i < order.size(); i++) {
            if(untouched.find(order[i]) == untouched.end()) continue;
            oi = order[i]; break;
        }
        if(oi == -999) { oi = *(untouched.begin()); order.insert(order.begin(), oi); }
        else {
            if(children.find(oi) != children.end()) { cout <<"children problem"<< endl; exit(0); }
            children[oi] = vector<int>();
            for(set<int>::iterator it=edges[oi].begin(); it != edges[oi].end(); ++it)
                if(untouched.find(*it) != untouched.end() && find(order.begin(),order.end(),*it) == order.end()) {
                    order.insert(find(order.begin(), order.end(), oi), *it);
                    children[oi].push_back(*it);
                }
            untouched.erase(oi);
        }
        for(int i=0; i < order.size(); i++) cout << order[i] <<" "; cout << endl;
    }

    map<int,int> parents; // children -> parent
    for(map<int,vector<int> >::iterator it=children.begin(); it != children.end(); ++it)
        for(int i=0; i < it->second.size(); i++)
            parents[ it->second [i] ] = it->first;

    vector<int> b_arts; // find articulation points per component based upon parent component
    for(int i=0; i < bcomps.size(); i++) {
        vector<int> intxn(0);
        if(parents.find(i) != parents.end()) {
            std::back_insert_iterator<vector<int> > intit (intxn);
            set_intersection(
                bcomps[i].begin(), bcomps[i].end(),
                bcomps[ parents[i] ].begin(), bcomps[ parents[i] ].end(), intit);
        }
        if(intxn.size() > 1) { cout << "not expecting intersection > 1 " << intxn[0] << " " << intxn[1] << endl; cout.flush(); exit(0); }
        if(intxn.size() == 1) b_arts.push_back(intxn[0]);
        else { // try finding something that is not already a collapse-point
            for(int k=0; k < bcomps[i].size(); k++)
                if( find(b_arts.begin(),b_arts.end(),bcomps[i][k]) == b_arts.end() || k == bcomps[i].size()-1 )
                    { b_arts.push_back(bcomps[i][k]); break; }
        }
    }

    retcomps.clear(); retCollapseInto.clear(); // reorder components and artipts according to the order calcd
    for(int i=0; i < order.size(); i++) {
        int ci = order[i];
        retcomps.push_back(vector<int>( bcomps[ci].begin(),bcomps[ci].end() ));
        retCollapseInto.push_back( b_arts[ci] );
    }
}

void findBC(vector<vector<int> > & edges, vector<vector<int> > & components, vector<int> & collapseInto) {
    int numnodes = -1;
    for(int ei=0; ei < edges.size(); ei++) {
        if(numnodes < edges[ei][0]) numnodes = edges[ei][0];
        if(numnodes < edges[ei][1]) numnodes = edges[ei][1];
    }
    cout << "NUMNODES " << numnodes <<" "<< edges.size() << endl;
    Graph g(numnodes+1);
    for(int ei=0; ei < edges.size(); ei++) { add_edge(edges[ei][0], edges[ei][1], g); }

    components.clear(); collapseInto.clear();
    findBiconn(g, components, collapseInto);

    //for(int i=0; i < components.size(); i++) {
    //    cout << "[";
    //    for(int k=0; k < components[i].size(); k++) cout << components[i][k] << ", ";
    //    cout << "] " << collapseInto[i] << endl;
    //}
}
