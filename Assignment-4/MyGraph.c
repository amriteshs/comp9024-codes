#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

// macro definitions
#define SIZE 9999
#define TRUE 1
#define FALSE 0

// a vertex is a 2D point
typedef struct Vertex {
    int x;  // x-coordinate
    int y;  // y-coordinate
} Vertex;

// each edge is a pair of vertices (end-points)
typedef struct Edge {
    Vertex *p1;     // first end point
    Vertex *p2;     // second end point
} Edge;

typedef struct VertexNode {
    Vertex *v;                  // vertex
    struct VertexNode *next;    // pointer to next node
    struct VertexNode *prev;    // pointer to previous node
    struct AdjNode *headadj;    // first element of adjacency list
    struct AdjNode *tailadj;    // last element of adjacency list
    struct HeapNode *loc;       // location of vertex in heap
    int visited;                // to check if node has been visited
    double minDist;             // minimum distance for Djikstra's algorithm
} VertexNode;

typedef struct AdjNode {
	VertexNode *vnode;      // vertex node
    struct AdjNode *next;   // next node in adjacency list
} AdjNode;

typedef struct GraphRep *Graph;

typedef struct GraphRep {   // graph header
	VertexNode *vertices;   // an array of vertices or a linked list of vertices
	VertexNode *tailv;      // last vertex in graph
	int nV;                 // #vertices
	int nE;                 // #edges
} GraphRep;

typedef struct QueueNode {
    VertexNode *vnode;      // vertex node
    struct QueueNode *next; // next node in queue
} QueueNode;

typedef struct Queue {
    QueueNode *head;    // queue head
    QueueNode *tail;    // queue tail
} Queue;

typedef struct HeapNode {
	double key;                 // key of task
	VertexNode *vnode;          // vertex node
	struct HeapNode *parent;    // parent of node
	struct HeapNode *left;      // left child of node
	struct HeapNode *right;     // right child of node
} HeapNode;

// data type for a priority queue (heap)
typedef struct BinaryHeap {    // this is heap header
	int size;           // count of items in the heap
    HeapNode *root;     // root node of heap
    HeapNode *tailLeaf; // last leaf node in heap
} BinaryHeap;

// create a new vertex node
VertexNode *newVertexNode(Vertex *vx)
{
    VertexNode *new = malloc(sizeof(VertexNode));
	assert(new != NULL);

	new -> v = vx;
	new -> next = NULL;
	new -> prev = NULL;
	new -> headadj = NULL;
	new -> tailadj = NULL;
	new -> loc = NULL;
	new -> visited = FALSE;
	new -> minDist = 0.0;

	return new;
}

// create a new node for adjacency list
AdjNode *newAdjNode(VertexNode *vn)
{
    AdjNode *new = malloc(sizeof(AdjNode));
	assert(new != NULL);

	new -> vnode = vn;
	new -> next = NULL;

	return new;
}

// create a new queue node
QueueNode *newQueueNode(VertexNode *vn)
{
    QueueNode *new = malloc(sizeof(QueueNode));
	assert(new != NULL);

	new -> vnode = vn;
	new -> next = NULL;

	return new;
}

// create a new queue
Queue *newQueue(QueueNode *qnode)
{
    Queue *new = malloc(sizeof(Queue));
	assert(new != NULL);

	new -> head = qnode;
	new -> tail = qnode;

	return new;
}

// create a new heap node to store an item
HeapNode *newHeapNode(double k, VertexNode *vn)
{
	HeapNode *new = malloc(sizeof(HeapNode));
	assert(new != NULL);

	new -> key = k;
	new -> vnode = vn;
    new -> parent = NULL;
    new -> left = NULL;
    new -> right = NULL;

	return new;
}

// create a new empty heap-based priority queue
BinaryHeap *newHeap()
{
    // this function creates an empty binomial heap-based priority queue
	BinaryHeap *T = malloc(sizeof(BinaryHeap));
	assert (T != NULL);

	T -> size = 0;
    T -> root = NULL;
    T -> tailLeaf = NULL;

	return T;
}

// restore heap property by swapping node along an upward path from point of insertion
// since heap has height O(log n), time complexity is O(log n)
HeapNode *UpHeap(BinaryHeap *T, HeapNode *hnode)
{
    // swapping at root node
    if (hnode -> parent == NULL) {
        T -> root = hnode;
        hnode -> vnode -> loc = hnode;

        return hnode;
    }

    // if key of node >= key of its parent, no swapping takes place
    if (hnode -> key >= hnode -> parent -> key) {
        hnode -> vnode -> loc = hnode;

        return hnode;
    }

    // key of node < key of its parent
    // swapping takes place
    HeapNode *temp;

    // node is left child of its parent
    if (hnode -> parent -> left == hnode) {
        // swap node and its parent
        temp = newHeapNode(hnode -> key, hnode -> vnode);
        temp -> left = hnode -> parent;
        temp -> right = hnode -> parent -> right;

        temp -> left -> left = hnode -> left;

        if (temp -> left -> left != NULL) {
            temp -> left -> left -> parent = temp -> left;
        }

        temp -> left -> right = hnode -> right;

        if (temp -> left -> right != NULL) {
            temp -> left -> right -> parent = temp -> left;
        }

        if (temp -> right != NULL) {
            temp -> right -> parent = temp;
        }

        if (hnode -> parent -> parent != NULL) {
            temp -> parent = hnode -> parent -> parent;

            if (hnode -> parent -> parent -> left == hnode -> parent) {
                hnode -> parent -> parent -> left = temp;
            } else {
                hnode -> parent -> parent -> right = temp;
            }
        }

        if (temp -> left != NULL) {
            temp -> left -> parent = temp;
        }

        hnode = temp;
        hnode -> vnode -> loc = hnode;
        UpHeap(T, hnode);   // recursive call to parent
        // thus, swapping proceeds in upward fashion

        return hnode -> left;
    }

    // node is right child of its parent
    // swap node and its parent
    temp = newHeapNode(hnode -> key, hnode -> vnode);
    temp -> left = hnode -> parent -> left;
    temp -> right = hnode -> parent;

    temp -> right -> left = hnode -> right;

    if (temp -> right -> left != NULL) {
        temp -> right -> left -> parent = temp -> right;
    }

    temp -> right -> right = hnode -> right;

    if (temp -> right -> right != NULL) {
        temp -> right -> right -> parent = temp -> right;
    }

    if (temp -> left != NULL) {
        temp -> left -> parent = temp;
    }

    if (hnode -> parent -> parent != NULL) {
        temp -> parent = hnode -> parent -> parent;

        if (hnode -> parent -> parent -> left == hnode -> parent) {
            hnode -> parent -> parent -> left = temp;
        } else {
            hnode -> parent -> parent -> right = temp;
        }
    }

    if (temp -> right != NULL) {
        temp -> right -> parent = temp;
    }

    hnode = temp;
    hnode -> vnode -> loc = hnode;
    UpHeap(T, hnode);   // recursive call to parent
    // thus, swapping proceeds in upward fashion

    return hnode -> right;
}

// restore heap property by swapping node along a downward path from root
// since heap has height O(log n), time complexity is O(log n)
HeapNode *DownHeap(BinaryHeap *T, HeapNode *hnode)
{
    // swapping at root
    if (hnode -> left == NULL && hnode -> right == NULL) {
        hnode -> vnode -> loc = hnode;

        return hnode;
    }

    if (hnode -> key < hnode -> left -> key) {
        // swapping takes place between node and its left child
        HeapNode *temp = newHeapNode(hnode -> left -> key, hnode -> left -> vnode);
        temp -> left = hnode;
        temp -> right = hnode -> right;

        hnode -> left = hnode -> left -> left;
        hnode -> right = hnode -> left -> right;

        if (hnode -> parent == NULL) {
            T -> root = temp;
        } else {
            if (hnode -> parent -> left == hnode) {
                hnode -> parent -> left = temp;
            } else {
                hnode -> parent -> right = temp;
            }
        }

        hnode -> parent = temp;
        temp -> vnode -> loc = temp;
        DownHeap(T, temp);  // recursive call to child
        // thus, swapping proceeds in downward fashion

        return hnode;
    }

    // swapping takes place between node and its right child
    HeapNode *temp = newHeapNode(hnode -> right -> key, hnode -> right -> vnode);
    temp -> left = hnode -> left;
    temp -> right = hnode;

    hnode -> left = hnode -> right -> left;
    hnode -> right = hnode -> right -> right;

    if (hnode -> parent == NULL) {
        T -> root = temp;
    } else {
        if (hnode -> parent -> left == hnode) {
            hnode -> parent -> left = temp;
        } else {
            hnode -> parent -> right = temp;
        }
    }

    hnode -> parent = temp;
    temp -> vnode -> loc = temp;
    DownHeap(T, temp);  // recursive call to child
    // thus, swapping proceeds in downward fashion

    return hnode;
}

// insert new heap node to binary heap
// UpHeap() is called which takes O(log n) time, since upward traversal through heap takes O(log n) time
// therefore, time complexity - O(log n), where n is number of vertices in graph
HeapNode *Insert(BinaryHeap *T, double k, VertexNode *vn)
{
    // create new heap node
    HeapNode *hnode = newHeapNode(k, vn);

    // insert at root if heap is empty
	if (T -> size == 0) {
		T -> root = hnode;
		T -> tailLeaf = UpHeap(T, hnode);   // last leaf node in heap
		T -> size += 1;

		return hnode;
	}

	// heap only contains root node
	// insertion takes place at left child of root node
	if (T -> size == 1) {
		hnode -> parent = T -> root;
        T -> root -> left = hnode;
        T -> tailLeaf = UpHeap(T, hnode);   // last leaf node in heap
        T -> size += 1;

        return hnode;
	}

	// last leaf node is left child of its parent
	// insertion takes place at right child of its parent
    if (T -> tailLeaf -> parent -> right == NULL) {
        T -> tailLeaf -> parent -> right = hnode;
        hnode -> parent = T -> tailLeaf -> parent;
        T -> tailLeaf = UpHeap(T, hnode);   // last leaf node in heap
        T -> size += 1;

        return hnode;
    }

    // last leaf node is right child of its parent
    // traverse upward till root or leaf node is found
    // move up until find root or leaf node
    while (TRUE) {
        // last leaf node is root node
        if (T -> tailLeaf == T -> root) {
            // traverse through left child of last leaf node
            // insertion takes place here
            while (T -> tailLeaf -> left != NULL) {
                T -> tailLeaf = T -> tailLeaf -> left;
            }

            T -> tailLeaf -> left = hnode;
            hnode -> parent = T -> tailLeaf;
            break;
        }

        // otherwise, if last leaf node is left child of its parent
        // update it as right child of its parent
        if (T -> tailLeaf -> parent -> left == T -> tailLeaf) {
            T -> tailLeaf = T -> tailLeaf -> parent -> right;

            // traverse downward till left child becomes null
            // insertion takes place at left child
            while (T -> tailLeaf -> left != NULL) {
                T -> tailLeaf = T -> tailLeaf -> left;
            }

            T -> tailLeaf -> left = hnode;
            hnode -> parent = T -> tailLeaf;
            break;
        }

        T -> tailLeaf = T -> tailLeaf -> parent;
    }

	T -> tailLeaf = UpHeap(T, hnode);   // last leaf node in heap
	T -> size += 1;

	return hnode;
}

// move the last leaf node to current node
// update the last leaf node
// last leaf node found by recursively traversing down the heap (therefore, will take O(log n) time
// thus, time complexity will be O(log n)
HeapNode *RemoveMinRecursive(BinaryHeap *T, HeapNode *hnode)
{
    // heap node has no child nodes
	if (hnode -> left == NULL && hnode -> right == NULL) {
        // root node
        // make heap root null
		if (hnode -> parent == NULL){
			T -> size -= 1;
			T -> root = NULL;

			return hnode;
		}

		// leaf node
		// compute second last leaf node in heap
        HeapNode *newTailLeaf;

        if (T -> size <= 1) {
            // heap has 0 or 1 nodes; last leaf node becomes null after deletion
            newTailLeaf = NULL;
        } else if (T -> size == 2) {
            // heap has 2 nodes; last leaf node set to heap root after deletion
            newTailLeaf = T -> root;
        } else {
            if (T -> tailLeaf -> parent -> right == T -> tailLeaf) {
                // last leaf node is right child of its parent
                // set it to left child of its parent
                newTailLeaf = T -> tailLeaf -> parent -> left;
            } else {
                // last leaf node is left child of its parent
                newTailLeaf = T -> tailLeaf -> parent;
                int flag = FALSE;

                while (newTailLeaf -> parent != NULL) {
                    // traverse upward till node with right child is encountered
                    if (newTailLeaf -> parent -> right == newTailLeaf) {
                        newTailLeaf = newTailLeaf -> parent -> left;
                        flag = TRUE;
                    }

                    newTailLeaf = newTailLeaf -> parent;
                }

                while (newTailLeaf -> right != NULL) {
                    // traverse downward through the node till node with no right child is encountered
                    newTailLeaf = newTailLeaf -> right;
                }
            }
        }

        // last leaf node is current node
        if (T -> tailLeaf == hnode) {
            T -> tailLeaf = newTailLeaf;
            T -> size -= 1;

            return NULL;
        }

        // otherwise, create new heap node from last leaf node
        HeapNode *newhn = newHeapNode(T -> tailLeaf -> key, T -> tailLeaf -> vnode);

        if (T -> tailLeaf -> parent -> left == T -> tailLeaf) {
            T -> tailLeaf -> parent -> left = NULL;
        } else {
            T -> tailLeaf -> parent -> right = NULL;
        }

        if (newTailLeaf == hnode) {
            // update last leaf node to new node if second last leaf node is current node
            T -> tailLeaf = newhn;
        } else {
            // otherwise, update last leaf node to second last leaf node
            T -> tailLeaf = newTailLeaf;
        }

        T -> size -= 1;

        return newhn;
	}

	// node has only left child
	if (hnode -> right == NULL) {
        // create new heap node
		HeapNode *newhn = newHeapNode(hnode -> left -> key, hnode -> left -> vnode);
		newhn -> parent = hnode -> parent;

		// compute second last leaf node
        HeapNode *newTailLeaf;

        if (T -> size <= 1) {
            // heap has 0 or 1 nodes; last leaf node becomes null after deletion
            newTailLeaf = NULL;
        } else if (T -> size == 2) {
            // heap has 2 nodes; last leaf node set to heap root after deletion
            newTailLeaf = T -> root;
        } else {
            if (T -> tailLeaf -> parent -> right == T -> tailLeaf) {
                // last leaf node is right child of its parent
                // set it to left child of its parent
                newTailLeaf = T -> tailLeaf -> parent -> left;
            } else {
                // last leaf node is left child of its parent
                newTailLeaf = T -> tailLeaf -> parent;
                int flag = FALSE;

                while (newTailLeaf -> parent != NULL) {
                    // traverse upward till node with right child is encountered
                    if (newTailLeaf -> parent -> right == newTailLeaf) {
                        newTailLeaf = newTailLeaf -> parent -> left;
                        flag = TRUE;
                    }

                    newTailLeaf = newTailLeaf -> parent;
                }

                while (newTailLeaf -> right != NULL) {
                    // traverse downward through the node till node with no right child is encountered
                    newTailLeaf = newTailLeaf -> right;
                }
            }
        }

        // set last leaf node to second last leaf node
        T -> tailLeaf = newTailLeaf;

        // if new node is at root, update root
		if (newhn -> parent == NULL) {
            T -> root = newhn;
        }

        T -> size -= 1;

		return newhn;
	}

	// current node has both left and right child nodes
    if (hnode -> right != NULL) {
        // key of left child <= key of right child
        if (hnode -> right -> key <= hnode -> left -> key) {
            // create new node from current node's right child
            HeapNode *newhn = newHeapNode(hnode -> right -> key, hnode -> right -> vnode);
            newhn -> left = hnode -> left;
            newhn -> right = RemoveMinRecursive(T, hnode -> right); // recursive call to right child of current node
            newhn -> parent = hnode -> parent;
            newhn -> vnode -> loc = newhn;

            if (newhn -> left != NULL) {
                newhn -> left -> parent = newhn;
            }

            if (newhn -> right != NULL) {
                newhn -> right -> parent = newhn;
            }

            if (newhn -> parent == NULL) {
                T -> root = newhn;
            }

            return newhn;
        }
    }

    // key of left child > key of right child
    // create new node from current node's left child
    HeapNode *newhn = newHeapNode(hnode -> left -> key, hnode -> left -> vnode);
    newhn -> left = RemoveMinRecursive(T, hnode -> left);   // recursive call to left child of current node

    if (T -> tailLeaf != hnode -> right) {
        newhn -> right = hnode -> right;
    }

    newhn -> parent = hnode -> parent;

    if (newhn -> left != NULL) {
        newhn -> left -> parent = newhn;
    }

    if (newhn -> right != NULL) {
        newhn -> right -> parent = newhn;
    }

    if (newhn -> parent == NULL) {
        T -> root = newhn;
    }

    return newhn;
}

// RemoveMinRecursive() is called, which takes O(log n) time
// this is because it traverses downward through the tree
// therefore, time complexity - O(log n), where n is number of vertices in graph
HeapNode *RemoveMin(BinaryHeap *T)
{
	HeapNode *hnode = T -> root;
	RemoveMinRecursive(T, hnode);

 	return hnode;
}

// A vertex node stores a vertex and other information, and you need to expand this type

// The above types serve as a starting point only. You need to expand them and add more types.
// Watch the lecture video between 7:50pm-8:20 or so for my discussion about this assignment

// create a new empty graph; takes constant time
// therefore, time complexity - O(1)
Graph CreateEmptyGraph()
{
    Graph g = malloc(sizeof(GraphRep));
    assert (g != NULL);

    g -> vertices = NULL;
    g -> tailv = NULL;
    g -> nV = 0;
    g -> nE = 0;

    return g;
}

// insert an edge into the graph
// takes O(m + n) time in worst case, because we need to traverse through adjacency list to check if edge already exists in graph
// therefore, time complexity - O(m + n), where m and n are number of edges and vertices respectively
int InsertEdge(Graph g, Edge *e)
{
    VertexNode *vn1 = NULL;
    VertexNode *vn2 = NULL;
    int flag1 = FALSE;
    int flag2 = FALSE;

    // graph already has vertices
    if (g -> vertices != NULL) {
        for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
            if (node -> v -> x == e -> p1 -> x && node -> v -> y == e -> p1 -> y) {
                // first vertex of edge
                vn1 = node;
                flag1 = TRUE;
            } else if (node -> v -> x == e -> p2 -> x && node -> v -> y == e -> p2 -> y) {
                // second vertex of edge
                vn2 = node;
                flag2 = TRUE;
            }

            if (flag1 && flag2) {
                break;
            }
        }
    }

    AdjNode *adj1;
    AdjNode *adj2;

    // if both vertices exist in graph
    if (vn1 != NULL && vn2 != NULL) {
        AdjNode *adj = NULL;

        // check if edge already exists by checking existence of second vertex in adjacency list of first vertex
        if (vn1 -> headadj != NULL) {
            for (AdjNode *node = vn1 -> headadj; node != NULL; node = node -> next) {
                if (node -> vnode == vn2) {
                    adj = node;
                    break;
                }
            }
        }

        // edge already exists
        if (adj != NULL) {
            return 0;
        }

        // add both vertices to each other's adjacency lists
        adj1 = newAdjNode(vn2);
        vn1 -> tailadj -> next = adj1;
        vn1 -> tailadj = adj1;

        adj2 = newAdjNode(vn1);
        vn2 -> tailadj -> next = adj2;
        vn2 -> tailadj = adj2;
    } else if (vn1 != NULL) {
        // second vertex does not exist
        // create second vertex
        // add both vertices to each other's adjacency lists
        vn2 = newVertexNode(e -> p2);
        g -> tailv -> next = vn2;
        g -> tailv = vn2;

        adj1 = newAdjNode(vn2);
        vn1 -> tailadj -> next = adj1;
        vn1 -> tailadj = adj1;

        adj2 = newAdjNode(vn1);
        vn2 -> headadj = adj2;
        vn2 -> tailadj = adj2;

        g -> nV++;
    } else if (vn2 != NULL) {
        // first vertex does not exist
        // create first vertex
        // add both vertices to each other's adjacency lists
        vn1 = newVertexNode(e -> p1);
        g -> tailv -> next = vn1;
        g -> tailv = vn1;

        adj1 = newAdjNode(vn2);
        vn1 -> headadj = adj1;
        vn1 -> tailadj = adj1;

        adj2 = newAdjNode(vn1);
        vn2 -> tailadj -> next = adj2;
        vn2 -> tailadj = adj2;

        g -> nV++;
    } else {
        // both vertices do not exist
        // create new vertex nodes for both
        vn1 = newVertexNode(e -> p1);
        vn2 = newVertexNode(e -> p2);

        if (g -> vertices != NULL) {
            g -> tailv -> next = vn1;
        } else {
            g -> vertices = vn1;
        }

        vn1 -> next = vn2;
        g -> tailv = vn2;

        // add both vertices to each other's adjacency lists
        adj1 = newAdjNode(vn2);
        vn1 -> headadj = adj1;
        vn1 -> tailadj = adj1;

        adj2 = newAdjNode(vn1);
        vn2 -> headadj = adj2;
        vn2 -> tailadj = adj2;

        g -> nV += 2;
    }

    g -> nE++;

    return 1;
}

// delete an edge from the graph
// takes O(m + n) time in worst case, because we need to traverse through adjacency list to check if edge exists in graph
// therefore, time complexity - O(m + n), where m and n are number of edges and vertices respectively
void DeleteEdge(Graph g, Edge *e)
{
    VertexNode *vn1 = NULL;
    VertexNode *vn2 = NULL;
    int flag1 = FALSE;
    int flag2 = FALSE;

    if (g -> vertices != NULL) {
        for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
            if (node -> v -> x == e -> p1 -> x && node -> v -> y == e -> p1 -> y) {
                // first vertex of edge
                vn1 = node;
                flag1 = TRUE;
            } else if (node -> v -> x == e -> p2 -> x && node -> v -> y == e -> p2 -> y) {
                // second vertex of edge
                vn2 = node;
                flag2 = TRUE;
            }

            if (flag1 && flag2) {
                break;
            }
        }
    }

    // edge does not exist; do nothing
    if (vn1 == NULL || vn2 == NULL) {
        return;
    }

    // free edge
    //free(e);

    int flag = FALSE;
    AdjNode *tmpA;
    VertexNode *tmpV;

    // delete second vertex from adjacency list of first vertex
    // and update connections in the adjacency list
    if (vn1 -> headadj -> vnode == vn2) {
        flag = TRUE;
        tmpA = vn1 -> headadj;

        if (vn1 -> headadj == vn1 -> tailadj) {
            vn1 -> headadj = NULL;
            vn1 -> tailadj = NULL;
        } else {
            vn1 -> headadj = vn1 -> headadj -> next;
        }

        // free vertex node
        free(tmpA);
    } else {
        for (AdjNode *node = vn1 -> headadj; node -> next != NULL; node = node -> next) {
            if (node -> next -> vnode == vn2) {
                if (node -> next == vn1 -> tailadj) {
                    vn1 -> tailadj = node;
                }

                tmpA = node -> next;
                node -> next = node -> next -> next;
                free(tmpA);     // free vertex node
                flag = TRUE;
                break;
            }
        }
    }

    // edge does not exist since second vertex not found in adjacency list of first vertex
    // do nothing
	if (!flag) {
        return;
	}

	// delete first vertex from adjacency list of second vertex
    // and update connections in the adjacency list
	if (vn2 -> headadj -> vnode == vn1) {
        tmpA = vn2 -> headadj;

        if (vn2 -> headadj == vn2 -> tailadj) {
            vn2 -> headadj = NULL;
            vn2 -> tailadj = NULL;
        } else {
            vn2 -> headadj = vn2 -> headadj -> next;
        }

        // free vertex node
        free(tmpA);
    } else {
        for (AdjNode *node = vn2 -> headadj; node -> next != NULL; node = node -> next) {
            if (node -> next -> vnode == vn1) {
                if (node -> next == vn2 -> tailadj) {
                    vn2 -> tailadj = node;
                }

                tmpA = node -> next;
                node -> next = node -> next -> next;
                free(tmpA);     // free vertex node
                break;
            }
        }
    }

    // remove first vertex from graph if its adjacency list is empty, i.e., it is an isolated vertex
	if (vn1 -> headadj == NULL) {
        if (vn1 == g -> vertices) {
            if (vn1 == g -> tailv) {
                g -> vertices = NULL;
                g -> tailv = NULL;
            } else {
                g -> vertices = g -> vertices -> next;
            }
        } else {
            for (VertexNode *node = g -> vertices; node -> next != NULL; node = node -> next) {
                if (node -> next == vn1) {
                    if (node -> next == g -> tailv) {
                        g -> tailv = node;
                    }

                    tmpV = node -> next;
                    node -> next = node -> next -> next;
                    free(tmpV);     // free vertex node
                    break;
                }
            }
        }

        g -> nV--;
        free(vn1);
	}

	// remove second vertex from graph if its adjacency list is empty, i.e., it is an isolated vertex
	if (vn2 -> headadj == NULL) {
        if (vn2 == g -> vertices) {
            if (vn2 == g -> tailv) {
                g -> vertices = NULL;
                g -> tailv = NULL;
            } else {
                g -> vertices = g -> vertices -> next;
            }
        } else {
            for (VertexNode *node = g -> vertices; node -> next != NULL; node = node -> next) {
                if (node -> next == vn2) {
                    if (node -> next == g -> tailv) {
                        g -> tailv = node;
                    }

                    tmpV = node -> next;
                    node -> next = node -> next -> next;
                    free(tmpV);
                    break;
                }
            }
        }

        g -> nV--;
        free(vn2);
	}

	g -> nE--;
}

// BFS is performed to find reachable vertices, which takes O(m + n) time in worst case
// therefore, time complexity - O(m + n), where m and n are number of edges and vertices respectively
void ReachableVertices(Graph g, Vertex *v)
{
    // return if graph is empty
    if (g -> vertices == NULL) {
        return;
    }

    // initialize vertex node to null
    VertexNode *vn = NULL;

    for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
        if (node -> v -> x == v -> x && node -> v -> y == v -> y) {
            // set vertex node to current node if vertex is found
            vn = node;
            break;
        }
    }

    // return if vertex does not exist in graph
    if (vn == NULL) {
        return;
    }

    // create a new queue with the vertex node as its head
    Queue *q = newQueue(newQueueNode(vn));
    int flag1 = FALSE;
    int flag2 = FALSE;

    // traverse queue till it becomes empty
    while (q -> head != NULL) {
        // set current node to head of queue
        QueueNode *curr = q -> head;

        // check if current node has been visited
        if (!(curr -> vnode -> visited)) {
            // set current node as visited
            curr -> vnode -> visited = TRUE;

            // print reachable vertices
            if (!flag1) {
                flag1 = TRUE;
            } else {
                if (!flag2) {
                    printf("(%d,%d)", curr -> vnode -> v -> x, curr -> vnode -> v -> y);
                    flag2 = TRUE;
                } else {
                    printf(",(%d,%d)", curr -> vnode -> v -> x, curr -> vnode -> v -> y);
                }
            }

            // check if adjacency lisy of current list is empty
            if (curr -> vnode -> headadj != NULL) {
                // add all (unvisited) vertices in adjacency list of current node to queue
                for (AdjNode *node = curr -> vnode -> headadj; node != NULL; node = node -> next) {
                    if (!(node -> vnode -> visited)) {
                        QueueNode *newqn = newQueueNode(node -> vnode);
                        q -> tail -> next = newqn;
                        q -> tail = newqn;
                    }
                }
            }
        }

        // move to next node in queue by updating head
        q -> head = q -> head -> next;
    }

    // reset all vertices as unvisited
    for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
        node -> visited = FALSE;
    }

    printf("\n");
}

// print shortest path found
void PrintPath(VertexNode *vn)
{
    if (vn == NULL) {
        return;
    }

    PrintPath(vn -> prev);

    if (vn -> prev == NULL) {
        printf("(%d,%d)", vn -> v -> x, vn -> v -> y);
    } else {
        printf(",(%d,%d)", vn -> v -> x, vn -> v -> y);
    }
}

// calculate shortest path from vertex u to vertex v
// all vertices are visited in O(m + n) time
// for each (unvisited) vertex, either Insert() or UpHeap() / DownHeap() are called, which take O(log n) time
// therefore, time complxity - O((m + n)log n), where m and n are number of edges and vertices in graph respectively
void ShortestPath(Graph g, Vertex *u, Vertex *v)
{
    // return if graph is empty
    if (g -> vertices == NULL) {
        return;
    }

    VertexNode *srcNode = NULL;
    VertexNode *goalNode = NULL;
    int flag1 = FALSE;
    int flag2 = FALSE;

    // traverse graph to find source and goal vertex nodes
    for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
        if (node -> v -> x == u -> x && node -> v -> y == u -> y) {
            // source vertex
            srcNode = node;
            flag1 = TRUE;
        }

        if (node -> v -> x == v -> x && node -> v -> y == v -> y) {
            // goal vertex
            goalNode = node;
            flag2 = TRUE;
        }

        if (flag1 && flag2) {
            break;
        }
    }

    // return if either source or goal node are not present in graph
    if (srcNode == NULL || goalNode == NULL) {
        return;
    }

    // print node and return if source and goal nodes are same
    if (srcNode == goalNode) {
        printf("(%d,%d)\n", goalNode -> v -> x, goalNode -> v -> y);
        return;
    }

    // create new binary heap and insert source node with initial distance 0
    BinaryHeap *Q = newHeap();
    Insert(Q, 0.0, srcNode);

    while (Q -> root != NULL) {
        // remove node with minimum distance value from heap and set it as current node
        HeapNode *curr = RemoveMin(Q);

        // check if current node has been visited
        if (!(curr -> vnode -> visited)) {
            // set current node as visited
            curr -> vnode -> visited = TRUE;

            // traverse though adjacency list of current node
            for (AdjNode *node = curr -> vnode -> headadj; node != NULL; node = node -> next) {
                // continue if adjacent node is same as source node
                if (node -> vnode == srcNode) {
                    continue;
                }

                // calculate Euclidean distance between current heap node and its current adjacent node
                double distance = sqrt(pow(((curr -> vnode -> v -> x) - (node -> vnode -> v -> x)), 2) + pow(((curr -> vnode -> v -> y) - (node -> vnode -> v -> y)), 2));

                if (node -> vnode -> minDist == 0.0) {
                    // insert adjacent node with updated distance to binary heap if visited for first time
                    node -> vnode -> minDist = distance + curr -> vnode -> minDist;
                    node -> vnode -> prev = curr -> vnode;
                    Insert(Q, node -> vnode -> minDist, node -> vnode);
                } else if (node -> vnode -> minDist > (distance + curr -> vnode -> minDist)) {
                    if (node -> vnode -> loc != NULL) {
                        // adjacent node already present in heap
                        // perform upheap and downheap operations to update key
                        // and add adjacent node to current node's path
                        node -> vnode -> loc -> key = distance + node -> vnode -> minDist;
                        UpHeap(Q, node -> vnode -> loc);
                        DownHeap(Q, node -> vnode -> loc);
                        node -> vnode -> minDist = distance + curr -> vnode -> minDist;
                        node -> vnode -> prev = curr -> vnode;
                    } else {
                        // insert adjacent node with updated distance to binary heap if it has been previously removed from heap
                        Insert(Q, distance + node -> vnode -> minDist, node -> vnode);
                    }
                }
            }
        }

        // set current node's heap location as null (since it has been removed)
        curr -> vnode -> loc = NULL;
    }

    // if goal node has no previous vertex, then no path exists
    // otherwise, print path
    if (goalNode -> prev != NULL) {
        PrintPath(goalNode);
        printf("\n");
    }

    // reset all attributes of vertex node
    for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
        node -> visited = FALSE;
        node -> minDist = 0.0;
        node -> loc = NULL;
        node -> prev = NULL;
    }
}

// adjacency list is traversed to free all nodes, which takes O(m + n) time
// therefore, time complexity - O(m + n), where m is number of edges in graph
void FreeGraph(Graph g)
{
    // start from first vertex of graph
    VertexNode *vn = g -> vertices;

    while (vn != NULL) {
        // start from first element of adjacency list
        AdjNode *adj = vn -> headadj;

        // free each element in adjacency list of vertex node
        while (adj != NULL) {
            AdjNode *tmpA = adj;
            adj = adj -> next;
            free(tmpA);
        }

        // free vertex node
        VertexNode *tmpV = vn;
        vn = vn -> next;
        free(tmpV);
    }

    free(g);
}

// BFS is performed to traverse though adjacency list, which takes O(m + n) time in worst case
// therefore, time complexity - O(m + n), where m and n are number of edges and vertices respectively
void ShowGraph(Graph g)
{
    // return if graph is empty
    if (g -> vertices == NULL) {
        return;
    }

    // create new queue with first vertex of graph as its head
    Queue *q = newQueue(newQueueNode(g -> vertices));
    int flag = FALSE;

    // traverse through graph's vertices
    for (VertexNode *vn = g -> vertices; vn != NULL; vn = vn -> next) {
        // continue if vertex has already been visited
        if (vn -> visited) {
            continue;
        }

        // if queue becomes empty, set next unvisited vertex as head
        // this happens when unconnected components are present in graph
        if (q -> head == NULL) {
            q -> head = newQueueNode(vn);
        }

        // traverse queue till it becomes empty
        while (q -> head != NULL) {
            // set queue head as current node
            QueueNode *curr = q -> head;

            // check if current node has been visited
            if (!(curr -> vnode -> visited)) {
                // set current vertex as visited
                curr -> vnode -> visited = TRUE;

                // print vertices
                if (!flag) {
                    printf("(%d,%d)", curr -> vnode -> v -> x, curr -> vnode -> v -> y);
                    flag = TRUE;
                } else {
                    printf(",(%d,%d)", curr -> vnode -> v -> x, curr -> vnode -> v -> y);
                }

                // check if adjacency list of vertex exists
                if (curr -> vnode -> headadj != NULL) {
                    // add all (unvisited) vertices of adjacency list to queue
                    for (AdjNode *node = curr -> vnode -> headadj; node != NULL; node = node -> next) {
                        if (!(node -> vnode -> visited)) {
                            QueueNode *newqn = newQueueNode(node -> vnode);
                            q -> tail -> next = newqn;
                            q -> tail = newqn;
                        }
                    }
                }
            }

            // move to next node in queue by updating head
            q -> head = q -> head -> next;
        }
    }

    // reset all vertices as unvisited
    for (VertexNode *node = g -> vertices; node != NULL; node = node -> next) {
        node -> visited = FALSE;
    }

    printf("\n");
}

//sample main for testing
int main()
{
    Graph g1;
    Edge *e_ptr;
    Vertex *v1, *v2;

    // Create an empty graph g1;
    g1 = CreateEmptyGraph();

    // Create first connected component
    // Insert edge (0,0)-(0,10)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 0;
    v2 -> y = 10;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (0,0)-(5,6)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 5;
    v2 -> y = 6;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (0, 10)-(10, 10)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 10;

    v2 -> x = 10;
    v2 -> y = 10;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (0,10)-(5,6)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 10;

    v2 -> x = 5;
    v2 -> y = 6;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (0,0)-(5,4)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 5;
    v2 -> y = 4;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (5, 4)-(10, 4)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 5;
    v1 -> y = 4;

    v2 -> x = 10;
    v2 -> y = 4;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (5,6)-(10,6)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 5;
    v1 -> y = 6;

    v2 -> x = 10;
    v2 -> y = 6;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (10,10)-(10,6)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 10;
    v1 -> y = 10;

    v2 -> x = 10;
    v2 -> y = 6;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (10, 6)-(10, 4)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 10;
    v1 -> y = 6;

    v2 -> x = 10;
    v2 -> y = 4;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Create second connected component
    // Insert edge (20,4)-(20,10)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 20;
    v1 -> y = 4;

    v2 -> x = 20;
    v2 -> y = 10;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (20,10)-(30,10)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 20;
    v1 -> y = 10;

    v2 -> x = 30;
    v2 -> y = 10;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    // Insert edge (25,5)-(30,10)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 25;
    v1 -> y = 5;

    v2 -> x = 30;
    v2 -> y = 10;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    if (InsertEdge(g1, e_ptr) == 0) {
        printf("edge exists\n");
    }

    //Display graph g1
    ShowGraph(g1);

    // Find the shortest path between (0,0) and (10,6)
    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 10;
    v2 -> y = 6;

    ShortestPath(g1, v1, v2);

    free(v1);
    free(v2);

    // Delete edge (0,0)-(5, 6)
    e_ptr = (Edge *) malloc(sizeof(Edge));
    assert(e_ptr != NULL);

    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 5;
    v2 -> y = 6;

    e_ptr -> p1 = v1;
    e_ptr -> p2 = v2;

    DeleteEdge(g1, e_ptr);

    free(e_ptr);
    free(v1);
    free(v2);

    // Display graph g1
    ShowGraph(g1);

    // Find the shortest path between (0,0) and (10,6)
    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 10;
    v2 -> y = 6;

    ShortestPath(g1, v1, v2);

    free(v1);
    free(v2);

    // Find the shortest path between (0,0) and (25,5)
    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v2 = (Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    v2 -> x = 25;
    v2 -> y = 5;

    ShortestPath(g1, v1, v2);

    free(v1);
    free(v2);

    // Find reachable vertices of (0,0)
    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v1 -> x = 0;
    v1 -> y = 0;

    ReachableVertices(g1, v1);

    free(v1);

    // Find reachable vertices of (20,4)
    v1 = (Vertex *) malloc(sizeof(Vertex));
    assert(v1 != NULL);

    v1 -> x = 20;
    v1 -> y = 4;

    ReachableVertices(g1, v1);

    free(v1);

    // Free graph g1
    FreeGraph(g1);

    return 0;
}
