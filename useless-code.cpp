// EXTRA CODE I MIGHT USE LATER!!



    //EXAMPLE USAGE
	UndirectedGraph myGraph(2);
    
    vertex *v1 = myGraph.add_vertex(1, 3, 4, 5);
    vertex *v2 = myGraph.add_vertex(2, 6, 2, 7);
    vertex *v3 = myGraph.add_vertex(3, -1, 0, 1);
    vertex *v4 = myGraph.add_vertex(4, -2, 0, 10);
    vertex *v5 = myGraph.add_vertex(5, -2, 1, 10);

	myGraph.draw_edge(v1, v2, 100);
    myGraph.draw_edge(v1, v3, 200);
    myGraph.draw_edge(v4, v2, 500);
    myGraph.draw_edge(v4, v3, 500);
    myGraph.draw_edge(v5, v4, 8298);

    //Heap myHeap = Heap();
    /*
    myHeap.insert(v2);
    myHeap.insert(v3);
    myHeap.insert(v1);
    myHeap.insert(v4);
    myHeap.insert(v5); */

   /*  cout <<  "Minimum is #" << myHeap.return_min()->name <<endl;
    cout <<  "Minimum is #" << myHeap.return_min()->name <<endl;
    cout <<  "Minimum is #" << myHeap.return_min()->name <<endl;
 */


    //cout << "Number of vertices: " << myGraph.N << endl;
    //cout << "Number of edges: " << myGraph.M << endl;
    //cout << "Printing all edges: " << endl;
    myGraph.print_edges();
	return 0;









/* 
class Heap
{
    public:
        Heap();

        void heapify(size_t i);

        vertex* return_min();

        vertex* insert(vertex* new_vertex);

        vertex* del(int i);
        
        void update(size_t i, double new_key);

        bool heapified = false;

        size_t heap_size = 0;

        vertex* arr[MAXSIZE];

        //HELPERS THAT PROVIDE O(1) PERFORMANCE !

        int parent(int i) { return (i-1) / 2;}

        int left (int i) { return (2*i +1);}

        int right (int i) { return (2*i + 2);}

        void interchange(vertex **v, vertex **w)
        {
            vertex* temp = *v;
            *v = *w;
            *w = temp;
        }

};

Heap::Heap()
{
    this->heap_size = 0;  
} 


vertex* Heap::insert(vertex* new_vertex)
{
    if (heap_size > MAXSIZE) {
        throw std::runtime_error("Can't add more vertices since max capacity is reached");
    }
    
    ++heap_size;
    size_t i = heap_size - 1;
    this->arr[i] = new_vertex;

    while( arr[parent(i)]->key > arr[i]->key &&  i != 0 )// Repeat until you reach the parent or order respected
    {
        interchange(&arr[i], &arr[parent(i)]);
        i = parent(i);
    }
    this->heapified = true;
    
}


void Heap::heapify(size_t i)
{
    size_t le = left(i);
    size_t ri = right(i);
    size_t smallest = i;

    if (arr[le]->key < arr[i]->key && le < heap_size)
    {
        smallest = le;
    } 
    if (arr[ri]->key < arr[smallest]->key && ri < heap_size)
    {
        smallest = ri;
    }
    if (smallest != i)
    {
        interchange(&arr[i], &arr[smallest]);
        heapify(smallest); // Recursive call
    }
    
    this->heapified = true;
}

void Heap::update(size_t i, double new_key)
{
    assert(new_key <= arr[i]->key);
    arr[i]->key = new_key;
    while (i != 0 && arr[i]->key < arr[parent(i)]->key)
    {
        interchange(&arr[i], &arr[parent(i)]);
        i = parent(i);
    }
    this->heapified = false;
}

vertex* Heap::del(int i)
{
    update(i, NINF);
    this->heapified = false;
    return return_min();
}

vertex* Heap::return_min() 
{
    assert(this->heapified == true);
    if (heap_size <= 0) {
        cerr << "Invalid heap size" << endl;
    } 
    if (heap_size == 1)
    {
        --heap_size;
        return arr[0];
    } 

    vertex* root = arr[0];
    arr[0] = arr[heap_size - 1];
    --heap_size;
    heapify(0);

    return root;
}
 */
















/*
void UndirectedGraph::prim(unsigned int source)
{
	double dist[N] = {INF}
    int prev[N] = {0}
    list<vertex> S; 
    
}
*/










// MO

    Heap myHeap(100);
    UndirectedGraph myGraph(2);
    

    // TODO: CReate a bunch of vertices making sure to set key = 
    vertex *v1 = myGraph.add_vertex(1, 3, 4, INF);
    vertex *v2 = myGraph.add_vertex(2, 6, 2, 1);
    vertex *v3 = myGraph.add_vertex(3, -1, 0, INF);
    vertex *v4 = myGraph.add_vertex(4, -2, 0, 2);
    vertex *v5 = myGraph.add_vertex(5, -2, 1, INF);

	myGraph.draw_edge(v1, v2, 100);
    myGraph.draw_edge(v1, v3, 400);
    myGraph.draw_edge(v4, v2, 500);
    myGraph.draw_edge(v4, v3, 40);
    myGraph.draw_edge(v5, v4, 8298);

    myHeap.insert(v1);
    myHeap.insert(v2);
    myHeap.insert(v3);
    myHeap.insert(v4);

    //cout << "VERTICES[0]: " << myHeap.vertices[0]->name << endl;
    //cout << "VERTICES[1]: " << myHeap.vertices[1]->name << endl;

    myHeap.print_tree();
    
	cout << "EXECUTING PEEK(): " << myHeap.peek()->name << endl;
    cout << "EXECUTING EXTRACTMIN():" << endl; 
    cout << "OUTPUT: " << myHeap.extract_min()->name << endl;
    myHeap.print_tree();









    (para el primero:)
    generar weight para el primero aleatoriamente
    vertex* source = myGraph.add_vertex(name=1, x=3, y=4, key=0, weight_de_arriba);
    UndirectedGraph myGraph(dim, source);


    for (N-1)
    {
        WEIGHT = RANDOM
        if(WEIGHT < 2/dth root (N) ){
            vertex *v1 = myGraph.add_vertex(name=1, x=3, y=4, key=infinity, weight);
            agregar
        }
    }