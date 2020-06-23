#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define SIZE 9999
#define TRUE 1
#define FALSE 0

// This template is only a guide
// You need to include additional fields, data structures and functions

// data type for heap nodes
typedef struct HeapNode {
	// each node stores the priority (key), name, execution time,
	// release time and deadline of one task
	int key;                    // key of this task
	int TaskName;               // task name
	int Etime;                  // execution time of this task
	int Rtime;                  // release time of this task
	int Dline;                  // deadline of this task
	int order;                  // order of node
	struct HeapNode *child;     // child of node
	struct HeapNode *sibling;   // sibling of node
	struct HeapNode *parent;    // parent of node
} HeapNode;

// data type for a priority queue (heap)
typedef struct BinomialHeap {    // this is heap header
	int  size;              // count of items in the heap
    struct HeapNode *root;  // root node of heap
} BinomialHeap;

// create a new heap node to store an item (task)
HeapNode *newHeapNode(int k, int n, int c, int r, int d)
{
    // k:key, n:task name, c: execution time, r:release time, d:deadline
	HeapNode *new;
	new = malloc(sizeof(HeapNode));

	assert(new != NULL);

	new -> key = k;
	new -> TaskName = n;
	new -> Etime = c;
	new -> Rtime = r;
	new -> Dline = d;

    new -> order = 0;
    new -> child = NULL;
    new -> sibling = NULL;
    new -> parent = NULL;

	return new;
}

// create a new empty heap-based priority queue
BinomialHeap *newHeap()
{
    // this function creates an empty binomial heap-based priority queue
	BinomialHeap *T;
	T = malloc(sizeof(BinomialHeap));

	assert (T != NULL);

	T -> size = 0;
    T -> root = NULL;

	return T;
}

// merge two binomial heaps in non-decreasing order
// traversal takes place through the tree root nodes of heap
// one of the heaps has size 1 and the other has size n
// therefore, time complexity - O(log n), where n is number of tasks in binomial heap
HeapNode *Merge(BinomialHeap *T1, BinomialHeap *T2)
{
    // if one of the heaps is empty, return the other's root
    if (T1 -> root == NULL) {
        return T2 -> root;
    }

    if (T2 -> root == NULL) {
        return T1 -> root;
    }

    HeapNode *newR;
    HeapNode *N1;
    HeapNode *N2;

    // set heap root with lower order as root of result heap
    if (T1 -> root -> order <= T2 -> root -> order) {
        newR = T1 -> root;
        N1 = T1 -> root -> sibling;
        N2 = T2 -> root;
    } else {
        newR = T2 -> root;
        N1 = T1 -> root;
        N2 = T2 -> root -> sibling;
    }

    HeapNode *curr = newR;

    // append heap node with lower order as sibling of last tree root node in result heap
    // iterate to sibling of that heap node and repeat process
    while (N1 != NULL && N2 != NULL) {
        if (N1 -> order <= N2 -> order) {
            curr -> sibling = N1;
            N1 = N1 -> sibling;
        } else {
            curr -> sibling = N2;
            N2 = N2 -> sibling;
        }

        curr = curr -> sibling;
    }

    // append remaining heap as sibling of last tree root node in result heap
    if (N1 != NULL) {
        curr -> sibling = N1;
    }

    if (N2 != NULL) {
        curr -> sibling = N2;
    }

    return newR;
}

// find union of two binomial heaps
// traversal takes place through the tree root nodes of heap in Merge()
// one of the heaps has size 1 and the other has size n
// therefore, time complexity - O(log n), where n is number of tasks in binomial heap
HeapNode *Union(BinomialHeap *T1, BinomialHeap *T2)
{
    BinomialHeap *T = newHeap();
    T -> root = Merge(T1, T2);

    // return if heap is empty
    if (T -> root == NULL) {
        return NULL;
    }

    HeapNode *prev = NULL;              // initially, set previous node to NULL
    HeapNode *curr = T -> root;         // initially, set current node to root node of heap
    HeapNode *next = curr -> sibling;   // initially, set next node to sibling of root node of heap

    // iterate over tree root nodes of heap till next node becomes NULL
    while (next != NULL) {
        if (curr -> order != next -> order) {
            // move ahead if current node and next node have different orders
            prev = curr;
            curr = next;
        } else {
            if (next -> sibling != NULL && curr -> order == next -> sibling -> order) {
                // move ahead if current node and sibling of next node have different orders
                prev = curr;
                curr = next;
            } else {
                // order of current node and next node is same
                if (curr -> key <= next -> key) {
                    // set next node as child of current node if key of current node <= key of next node
                    curr -> sibling = next -> sibling;
                    next -> parent = curr;
                    next -> sibling = curr -> child;
                    curr -> child = next;
                } else {
                    // set root of heap as next node if current node is heap root when key of current node > key of next node
                    if (prev != NULL) {
                        prev -> sibling = next;
                    } else {
                        T -> root = next;
                    }

                    // set next node as child of current node if key of current node > key of next node
                    curr -> parent = next;
                    curr -> sibling = next -> child;
                    next -> child = curr;
                    curr = next;
                }

                // increment order of current node
                curr -> order += 1;
            }
        }

        next = curr -> sibling;
    }

    return T -> root;
}

// insert task in binomial heap
// traversal takes place through the tree root nodes of heap in Merge() called inside Union()
// the new heap created has size 1 and the heap with which its union is taken has size n
// therefore, time complexity - O(log n), where n is number of tasks in binomial heap
void Insert(BinomialHeap *T, int k, int n, int c, int r, int d)
{
    // k: key, n: task name, c: execution time, r: release time, d: deadline
    // no need to check if this task already exists in T

    // create a heap of order 0 with the node containing task details
    BinomialHeap *newT = newHeap();
    newT -> root = newHeapNode(k, n, c, r, d);

    // perform union operation on existing heap and the 0-order heap
    T -> root = Union(T, newT);
    T -> size += 1;

    free(newT);
}

// remove node with minimum key from binomial heap
// traversal takes place through the tree root nodes of heap in Union()
// therefore, time complexity - O(log n), where n is number of tasks in binomial heap
HeapNode *RemoveMin(BinomialHeap *T)
{
    if (T -> root == NULL) {
        return NULL;
    }

    // traverse tree root nodes of heap to find node with minimum key
    HeapNode *minNode = T -> root;

    for (HeapNode *curr = T -> root -> sibling; curr != NULL; curr = curr -> sibling) {
        if (curr -> key < minNode -> key) {
            minNode = curr;
        }
    }

    // traverse to node with minimum key and unlink it from the node that has it as its sibling
    if (minNode == T -> root) {
        T -> root = minNode -> sibling;
    } else {
        for (HeapNode *curr = T -> root; curr != NULL; curr = curr -> sibling) {
            if (curr -> sibling == minNode) {
                curr -> sibling = minNode -> sibling;
                break;
            }
        }
    }

    // unlink all nodes that have minimum node as parent
    for (HeapNode *curr = minNode -> child; curr != NULL; curr = curr -> sibling) {
        curr -> parent = NULL;
    }

    // create new heap of order 0 that contains minimum node
    BinomialHeap *newT = newHeap();

    newT -> root = minNode -> child;
    minNode -> child = NULL;
    minNode -> sibling = NULL;

    T -> root = Union(T, newT);
    T -> size -= 1;

    // set order of minimum node to 0
    minNode -> order = 0;

    return minNode;
}

// find minimum key value in binomial heap
// traversal takes place through the tree root nodes of heap
// therefore, time complexity - O(log n), where n is number of tasks in binomial heap
int Min(BinomialHeap *T)
{
    // traverse tree root nodes of heap to find minimum key value
    int min = T -> root -> key;

    for (HeapNode *curr = T -> root -> sibling; curr != NULL; curr = curr -> sibling) {
        if (curr -> key < min) {
            min = curr -> key;
        }
    }

    return min;
}

// task scheduling operation using binomial heap-based priority queues
// Insert() takes O(log n) time and takes place n times, where n is number of tasks
// Min() takes O(log n) time, where n is number of tasks
// RemoveMin() takes O(log n) time and takes place n times for each queue, where n is number of tasks
// only 2 priority queues used; the 3rd one not required for this approach; execution will still be completed in O(n log n) time
// this is because we are only running the while loop until the deadline priority queue becomes empty, which contains n tasks and calls RemoveMin() which takes O(log n) time
// therefore, time complexity - O(n log n), where n is number of tasks in binomial heap
int TaskScheduler(char *f1, char *f2, int m)
{
    int x;
    int len = 0;
    int lenR = 0;
    int list[SIZE];
    int result[SIZE];

    FILE *fp1 = fopen(f1, "r"); // open source file in read mode

    // exit if file does not exist
    if (!fp1) {
        printf("file1 does not exist\n");

        return -1;
    }

    // scan the integers present in file and store them in an array
    while (!feof(fp1)) {
        if (fscanf(fp1, "%d", &x) == 1) {   // number encountered in file
            list[len++] = x;
        } else {                            // non-numeric characters encountered in file
            char c = fgetc(fp1);

            if ((c >= 0 && c < '0') || (c > '9')) {
                printf("input error when reading the attribute of the task %d\n", list[(int)(len / 4)] + 1);

                return -1;
            }
        }
    }

    // close file1
    fclose(fp1);

    BinomialHeap *Q1 = newHeap();   // priority queue with release times as keys
    BinomialHeap *Q2 = newHeap();   // priority queue with deadlines as keys

    // insert all tasks in Q1 (earlier release time has higher priority)
    for (int i = 0; i < len; i += 4) {
        if (list[i] >= 0 && list[i + 2] >= 0 && list[i + 1] > 0 && list[i + 1] > 0) {
            Insert(Q1, list[i + 2], list[i], list[i + 1], list[i + 2], list[i + 3]);
        } else {
            printf("input error when reading the attribute of the task %d\n", list[i]);

            return -1;
        }
    }

    int flag = TRUE;
	int curr_time = 0;                      // set current time to 0 initially
	int core;                               // core being used
	int next_schedule_point[m + 1];         // array to store next possible scheduling time on each core

	// initialize scheduling points on all cores to 0
	for (int i = 1; i <= m; i++) {
        next_schedule_point[i] = Min(Q1);
	}

	// loop runs until all tasks have been executed
	while (flag) {
	    // remove all tasks from Q1 that have current time as release time and insert in Q2 (earlier deadline has higher priority)
        while (Q1 -> root != NULL && curr_time == Min(Q1)) {
            HeapNode *minNodeR = RemoveMin(Q1);
            Insert(Q2, minNodeR -> Dline, minNodeR -> TaskName, minNodeR -> Etime, minNodeR -> Rtime, minNodeR -> Dline);
        }

        // move to core 1 to schedule a task on it
        core = 1;

        // schedule tasks till Q2 becomes empty or all cores are in use
        while (Q2 -> root != NULL && core <= m) {
            // if current core has a task running on it, move to next core
            if (curr_time < next_schedule_point[core]) {
                core++;
                continue;
            }

            // remove task with earliest deadline from Q2
            HeapNode *minNodeD = RemoveMin(Q2);

            if (curr_time + minNodeD -> Etime > minNodeD -> key) {
                // exit when deadline is missed by a task; no feasible schedule
                return 0;
            } else {
                // store task name, core on which task is scheduled, and start time of task to result
                result[lenR++] = minNodeD -> TaskName;
                result[lenR++] = core;
                result[lenR++] = curr_time;

                // set core's next scheduling point to the time when the task on it completes its execution
                next_schedule_point[core] = curr_time + minNodeD -> Etime;
            }

            // move to next core to schedule a task on it
            core++;
        }

        if (Q1 -> root == NULL && Q2 -> root == NULL) {
            flag = FALSE;
        }

        // increment current time by 1
        curr_time++;
    }

    FILE *fp2 = fopen(f2, "w");     // open destination file in write mode

    // print result to file
    for (int i = 0; i < lenR; i += 3) {
        if (i == lenR - 3) {
            fprintf(fp2, "%d Core%d %d", result[i], result[i + 1], result[i + 2]);
        } else {
            fprintf(fp2, "%d Core%d %d ", result[i], result[i + 1], result[i + 2]);
        }
    }

    // close file2
    fclose(fp2);

    return 1;
}

int main() //sample main for testing
{
    int i;

    i = TaskScheduler("samplefile1.txt", "feasibleschedule1.txt", 4);

    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is a feasible schedule on 4 cores */

    i = TaskScheduler("samplefile1.txt", "feasibleschedule2.txt", 3);
    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is no feasible schedule on 3 cores */

    i = TaskScheduler("samplefile2.txt", "feasibleschedule3.txt", 5);
    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is a feasible schedule on 5 cores */

    i = TaskScheduler("samplefile2.txt", "feasibleschedule4.txt", 4);
    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is no feasible schedule on 4 cores */

    i = TaskScheduler("samplefile3.txt", "feasibleschedule5.txt", 2);
    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is no feasible schedule on 2 cores */

    i = TaskScheduler("samplefile4.txt", "feasibleschedule6.txt", 2);
    if (i == 0) {
        printf("No feasible schedule!\n\n");
    }
    /* There is a feasible schedule on 2 cores */

    return 0;
}
