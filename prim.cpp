#include <iostream>
#include <list>
#include <set>
#include <utility>
#include <vector>
#include<limits>
#include <random>
#include <cassert>
#include <math.h>



using namespace std;


#define INF numeric_limits<double>::infinity()
#define NINF numeric_limits<double>::infinity()*-1
#define MAXSIZE 2<<19
//typedef pair<set<unsigned int>, double> edge;



struct vertex
{
    size_t name;
    double x, y, z, w;
    list < pair<vertex*, double> > adj; // Adjacency list
    // For use in the heap data structure:
    //vertex *left = nullptr;
    //vertex *right = nullptr;
    double key; 

    vertex(unsigned int name, double x, double y, double z, double w, double key) : name(name), x(x), y(y), z(z), w(w), key(key) {} //Key = distance
};


void exchange(vertex **a, vertex **b)
{
    vertex *temp = *a;
    *a = *b;
    *b = temp;
}




class Heap
{
    
public:
    vertex** vertices; 
    size_t max_size; 
    size_t size; 
    // Constructor
    Heap(size_t capacity);
    // Destructor
    ~Heap();
  
    // Reorganizes heap to maintain priority order
    void heapify(size_t i);

     // Inserts a new vertex
    bool insert(vertex* v);

    // Extract and return the vertex w minimum key
    vertex* extract_min();

    // Expels vertex at i from queue and reorganizes heap
    void del(size_t i);

    // Returns vertex with lowest key
    vertex* peek() { return vertices[0]; }
    
    
  
    //TODO: void update_key(int i, int val);
  
    // USEFUL FUNCTIONS
  
    size_t parent(size_t i) { return (i-1)/2; }
  
    size_t left(size_t i) { return (2*i + 1); }
  
    size_t right(size_t i) { return (2*i + 2); }

    void print_tree() {
        unsigned int level = 0;
        unsigned int n_printed = 0;

        cout << "\n\n--------------PRINTING TREE------------------" << endl;
        
        cout << "SIZE IS: " << size<< endl;
        cout << "MAX_SIZE IS: " << max_size << endl;
        cout << "---------------------------------------------" << endl;

        while (n_printed < size )
        {
            //cout << "n_printed IS: " << n_printed << endl;
            //cout << "level IS: " << level << endl;

            size_t j = (1<<level) - 1;

            char margin[3*level + 1];
            for (size_t i = 0; i < 3*level; i++)
            {
                margin[i] = '-';
            }
            margin[3*level] = '\0';
            
            while (j < (1<<(level+1)) -1 && n_printed < size)
            {
                cout << margin << "LEVEL " << level << ": " << vertices[j]->name << endl;
                n_printed++;
                j++;
            }
            ++level;
        }
        cout << "\n\n" << endl;

        
    }
};
  
Heap::Heap(size_t capacity)
{
    size = 0;
    max_size = capacity;
    vertices = new vertex*[capacity];
}

Heap::~Heap()
{
    delete[] vertices;
}

// Inserts new vertex. Returns true iff successful.
bool Heap::insert(vertex* v)
{
    //cout << "INSERT() EXECUTED FOR VERTEX " << v->name << endl;
    if (size == max_size)
    {
        cerr << "\n Cannot insert more. Heap is full" << endl;
        return false;
    }
    if (size == 0)
    {
        vertices[0] = v;
        size++;
        return true;
    }
    
    // Insert at the end
    size_t it = size;
    //cout << "it IS: " << it << endl;
    vertices[it] = v;
    //vertices[it]->name = it;
    ++size;

    // Restructure if necessary. 
    while (it != 0 && vertices[parent(it)]->key > vertices[it]->key)
    {
        exchange(&vertices[it], &vertices[parent(it)]);
        it = parent(it);
    }

    

    return true;
    
}


void Heap::heapify(size_t i) 
{
    size_t le = left(i);
    size_t ri = right(i);
    if (le >= size && ri >= size) {
        return;
    }
    size_t min = i;
    if (vertices[le]->key < vertices[i]->key)
        min = le;
    if (vertices[ri]->key < vertices[min]->key)
        min = ri;
    
    if (min != i)
    {
        exchange(&vertices[i], &vertices[min]);
        heapify(min);
    }

}

vertex* Heap::extract_min()
{
    assert(size > 0);

    if (size == 1)
    {
        --size;
        return vertices[0];
    }

    vertex* res = vertices[0];
    vertices[0] = vertices[size - 1];
    --size;
    heapify(0);

    return res;
    
}


class UndirectedGraph 
{	
    unsigned int N;
    unsigned int M;
    int dim;


public:
    list < vertex > vertices;

	UndirectedGraph(int dim);

	void draw_edge(struct vertex* a, struct vertex* b, double weight);  

    vertex* add_vertex(unsigned int name, double x, double y, double z, double w, double key);

	double prim(); 

	void print_edges(); //FOR debugging

    vertex* source;
    
};


UndirectedGraph::UndirectedGraph(int dim) 
{
	this->N = 0;
    this->M = 0;
	this->dim = dim;
}

vertex* UndirectedGraph::add_vertex(unsigned int name, double x, double y, double z, double w, double key)
{
    vertex v = vertex(name, x, y, z, w, key);
    this->vertices.push_back(v);
    ++this->N;
    return &this->vertices.back();
    
}


void UndirectedGraph::draw_edge(struct vertex* a, struct vertex* b, double weight)
{
	assert(a);
    assert(b);
    assert(a != b);

    a->adj.push_back(make_pair(b, weight));
    b->adj.push_back(make_pair(a, weight));
    ++this->M;

}

void UndirectedGraph::print_edges()
{
	for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        cout << "Length of adj list for vertex #" << it->name << " is: " << it->adj.size()<< endl;
        cout << "Printing adj list for vertex #" << it->name << endl;
        for (auto adj_v = it->adj.begin(); adj_v != it->adj.end(); adj_v++)
        {
            cout << "Next: " << adj_v->first->name << endl;
        }

    }
}

double UndirectedGraph::prim()
{
    double length = 0; 
    set<vertex*> S;
    Heap priority_q(MAXSIZE);
    assert(this->source->key == 0);
    while (priority_q.size > 0)
    {
        vertex* v = priority_q.extract_min();
        S.insert(v);
        for (auto it = v->adj.begin(); it != v->adj.end(); ++it)
        {
            vertex* w = it->first;
            double l = it->second;
            if (S.count(w)) 
                continue;
            if (w->key > l) 
            {
                w->key = l;
                // Update prev if necessary
                length += l;
                priority_q.insert(w);
            }
        }
    }
    return length;   
}





int main() 
{

    const int N = 128;
    const int dim = 2;
    assert(dim > 0 && dim <= 4);
    UndirectedGraph myGraph(dim);
    vertex* s;
    
    //Add vertices
    for (size_t i = 0; i < N; ++i) {
        std::mt19937 gen{std::random_device{}()};
        std::uniform_real_distribution<> dist{0.0, 1.0};
        double loc[4];
        for (int i = 0; i < dim; ++i) {
            loc[i] = dist(gen);
        }
        for (int i = dim; i < 4-dim; i++) {
            loc[i] = 0;
        }

        if (i == 0) 
        {
            //SOURCE
            s = myGraph.add_vertex(1, loc[0], loc[1], loc[2], loc[3], 0);
            continue;
        }
        myGraph.add_vertex(i+1, loc[0], loc[1], loc[2], loc[3], INF);
    }
    myGraph.source = s;

  // Add edges
    for (auto a= myGraph.vertices.begin(); a != myGraph.vertices.end(); a++)
    {
        //cout << "Coordinates of Source are: (" << a->x << "," << a->y << "," << a->z << "," << a->w << ") " << endl;

        for(auto b = myGraph.vertices.begin(); b != myGraph.vertices.end(); b++)
        {
            if (b->name <= a->name)
            {
                continue;
            }
            myGraph.draw_edge(&(*a), &(*b), sqrt(pow((b->x - a->x), 2) + pow((b->y - a->y), 2) + pow((b->z - a->z), 2) + pow((b->w - a->w), 2)));
        }
    }


    // Now execute prim's algorithm
    cout << myGraph.prim() << endl;
  
    
    
    //myGraph.print_edges();
}




//TODO: delete()