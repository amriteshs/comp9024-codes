#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// macro definitions
#define SIZE 9999
#define TRUE 1
#define FALSE 0

// all the basic data structures and functions are included in this template
// you can add your own auxiliary functions as you like

// data structures representing DLList

// data type for nodes
typedef struct DLListNode {
    int value;  // value (int) of this list item
    struct DLListNode *prev;
    // pointer to previous node in list
	struct DLListNode *next;
    // pointer to next node in list
} DLListNode;

//data type for doubly linked lists
typedef struct DLList {
	int size;      // count of items in list
	DLListNode *first; // first node in list
	DLListNode *last;  // last node in list
} DLList;

// create a new DLListNode
DLListNode *newDLListNode(int it)
{
	DLListNode *new;

	new = malloc(sizeof(DLListNode));
	assert(new != NULL);

	new -> value = it;
	new -> prev = new -> next = NULL;

	return new;
}

// create a new empty DLList
DLList *newDLList()
{
	DLList *L;

	L = malloc(sizeof(struct DLList));
	assert (L != NULL);

	L -> size = 0;
	L -> first = NULL;
	L -> last = NULL;

	return L;
}

// create a DLList from a text file
// time complexity - O(n), where n is size of DLList created, since all elements of list are traversed through once
DLList *CreateDLListFromFileDlist(const char *filename)
{
    int list[SIZE];
    int len = 0;

    if (filename == "stdin") {  // input is given by user
        char *num = malloc(sizeof(char *));     // allocate memory to user input
        char *end;

        // take user input until empty line is given as input
        while (*fgets(num, 16, stdin) != '\n') {
            list[len++] = strtol(num, &end, 10);    // convert numeric string to integer and store in array
        }

        free(num);      // free memory space allocated to user input
        fflush(stdin);  // clear standard input buffer
    } else {                    // input is taken from file
        // open file in read mode
        FILE *fp = fopen(filename, "r");

        // assertion fails if file does not exist
        assert(fp != NULL);

        // scan first occurrence of integer in file and store as first element of array
        fscanf(fp, "%d", &list[len]);

        // scan each occurrence of integer in file and store it in array until the end of file is encountered
        while (!feof(fp)) {
            len++;
            fscanf(fp, "%d", &list[len]);
        }

        // close file
        fclose(fp);
    }

    DLList *L = newDLList();
    DLListNode *N;
    DLListNode *tmp;

    for (int i = 0; i < len; i++) {
        // create new list node setting its value to current array element
        N = newDLListNode(list[i]);

        if (i == 0) {
            L -> first = N;     // set node created in first loop iteration as the first (or head) node of list being created
        } else {
            if (i == len - 1) {
                L -> last = N;  // set node created in last loop iteration as the last (or tail) node of list being created
            }

            N -> prev = tmp;    // make previous node of current node point to temp node (currently pointing to current node)
            tmp -> next = N;    // make next node of temp node point to current node
        }

        tmp = N;    // update temp node and make it point to current node in each loop iteration
    }

    // set size of array as the size of list being created
    L -> size = len;

    return L;
}

// clone a DLList
// time complexity - O(n), where n is size of DLList u, since list is being linearly traversed once
DLList *cloneList(DLList *u)
{
    int len = 0;
    DLList *L = newDLList();
    DLListNode *N;
    DLListNode *tmp;

    for (DLListNode *node = u -> first; node != NULL; node = node -> next) {
        len++;
        N = newDLListNode(node -> value);   // create new list node setting its value to value of current list element

        if (L -> first == NULL) {
            L -> first = N;     // set current node as first node of new list
        } else {
            if (node -> next == NULL) {
                L -> last = N;  // when next node of current node is null, i.e., last list element is encountered, set it as last node of new list
            }

            N -> prev = tmp;    // make previous node of current node point to temp node (currently pointing to current node)
            tmp -> next = N;    // make next node of temp node point to current node
        }

        tmp = N;    // update temp node and make it point to current node in each loop iteration
    }

    // set size of array as the size of list being created
    L -> size = len;

    return L;
}

// compute the union of two DLLists u and v
// time complexity - O(m * n), where m and n are the sizes of DLLists u and v respectively
// this is because DLList v is traversed through, and in each iteration, DLList u is traversed through
DLList *setUnion(DLList *u, DLList *v)
{
    int ctr = 0;
    DLList *L = cloneList(u);       // create new list L by cloning list u
    DLListNode *N;
    DLListNode *tmp = L -> last;    // make temp node point to last node of list L initially

    // traverse through list v to find list elements with values that do not occur in list L (or list u)
    for (DLListNode *i = v -> first; i != NULL; i = i -> next) {
        int flag = TRUE;    // set flag initially to TRUE for each loop iteration

        for (DLListNode *j = u -> first; j != NULL; j = j -> next) {
            // set flag to FALSE and break loop if list u has node with same value as that of current node of v
            if (i -> value == j -> value) {
                flag = FALSE;
                break;
            }
        }

        // if flag is TRUE, i.e., value of current node of v does not occur in list u, update list L by appending the node to the tail of L
        if (flag) {
            ctr++;  // increment counter

            N = newDLListNode(i -> value);  // create new node setting its value to that of current node of v
            N -> prev = tmp;                // make previous node of current node point to temp node
            tmp -> next = N;                // make next node of temp node point to current node
            tmp = tmp -> next;              // move pointer to temp node to its next node
        }
    }

    // increment size of L (initially set to size of list u) by the number of elements of list v appended to it
    L -> size += ctr;

    return L;
}

// compute the intersection of two DLLists u and v
// time complexity - O(m * n), where m and n are the sizes of DLLists u and v respectively
// this is because DLList u is traversed through, and in each iteration, DLList v is traversed through
DLList *setIntersection(DLList *u, DLList *v)
{
    int ctr = 0;
    DLList *L = newDLList();
    DLListNode *N;
    DLListNode *tmp;

    // find nodes with values common to both lists u and v
    for (DLListNode *i = u -> first; i != NULL; i = i -> next) {
        for (DLListNode *j = v -> first; j != NULL; j = j -> next) {
            // if nodes with values common to both lists u and v found, create new node with same value, append it to list L and break loop
            if (i -> value == j -> value) {
                ctr++;  // increment counter

                N = newDLListNode(i -> value);  // create a new node setting its value to that of current node of list u

                if (L -> first == NULL) {
                    L -> first = N;     // set current node as first node of L if it has not been set yet
                } else {
                    N -> prev = tmp;    // make previous node of current node point to temp node
                    tmp -> next = N;    // make next node of temp node point to current node
                }

                tmp = N;    // make temp node point to current node

                break;
            }
        }
    }

    // set number of common elements found as the size of list L
    L -> size = ctr;

    return L;
}

// free up all space associated with list
// time complexity - O(n), where n is size of DLList L, since list is being linearly traversed once
void freeDLList(DLList *L)
{
    DLListNode *node = L -> first;  // make current node point to first list node initially
    DLListNode *tmp;

    // memory deallocation is carried out in the order opposite to the order in which it was allocated
    // first, memory allocated to each node is freed up
    while (node != NULL) {
        tmp = node -> next; // temp node points to next node of current node
        free(node);         // space allocated to current node is freed up
        node = tmp;         // current node is updated to point to temp node, i.e., next node of current node
    }

    // next, memory allocated to the list itself is freed up
    free(L);
}

// display items of a DLList
// time complexity - O(n), where n is size of DLList u, since list is being linearly traversed once
void printDLList(DLList *u)
{
    // starting at first list node, traverse through to the last node
    for (DLListNode *node = u -> first; node != NULL; node = node -> next) {
        printf("%d\n", node -> value);  // print value of current node
    }

    printf("\n");
}

int main()
{
    DLList *list1, *list2, *list3, *list4;

    list1 = CreateDLListFromFileDlist("File1.txt");
    printDLList(list1);

    list2 = CreateDLListFromFileDlist("File2.txt");
    printDLList(list2);

    list3 = setUnion(list1, list2);
    printDLList(list3);

    list4 = setIntersection(list1, list2);
    printDLList(list4);

    freeDLList(list1);
    freeDLList(list2);
    freeDLList(list3);
    freeDLList(list4);

    printf("please type all the integers of list1\n");
    list1 = CreateDLListFromFileDlist("stdin");

    printf("please type all the integers of list2\n");
    list2 = CreateDLListFromFileDlist("stdin");

    list3 = cloneList(list1);
    printDLList(list3);
    list4 = cloneList(list2);
    printDLList(list4);

    freeDLList(list1);
    freeDLList(list2);
    freeDLList(list3);
    freeDLList(list4);

    return 0;
}
