#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// macro definitions
#define SIZE 9999
#define TRUE 1
#define FALSE 0
#define MAX(x, y) ((x > y) ? x : y)             // returns maximum of x and y
#define MIN(x, y) ((x < y) ? x : y)             // returns minimum of x and y
#define HEIGHT(x) ((x) ? (x -> height) : -1)    // if node x is NULL, returns height of x as -1, else returns x -> height

// all the basic data structures and functions are included in this template
// you can add your own auxiliary functions as you like

// data type for AVL tree nodes
typedef struct AVLTreeNode {
	int key;    // key of this item
	int value;  // value (int) of this item
	int height; // height of the subtree rooted at this node
	struct AVLTreeNode *parent; // pointer to parent
	struct AVLTreeNode *left;   // pointer to left child
	struct AVLTreeNode *right;  // pointer to right child
	int rank;
} AVLTreeNode;

// data type for AVL trees
typedef struct AVLTree {
	int size;           // count of items in AVL tree
	AVLTreeNode *root;  // root
} AVLTree;

// create a new AVLTreeNode
AVLTreeNode *newAVLTreeNode(int k, int v )
{
	AVLTreeNode *new;
	new = malloc(sizeof(AVLTreeNode));

	assert(new != NULL);

	new -> key = k;
	new -> value = v;
	new -> height = 0;      // height of this new node is set to 0
	new -> left = NULL;     // this node has no child
	new -> right = NULL;
	new -> parent = NULL;   // no parent

	return new;
}

// create a new empty AVL tree
AVLTree *newAVLTree()
{
	AVLTree *T;
	T = malloc(sizeof(AVLTree));

	assert (T != NULL);

	T -> size = 0;
	T -> root = NULL;

	return T;
}

// major function declarations
AVLTree *CreateAVLTree(const char *filename);
AVLTree *CloneAVLTree(AVLTree *T);
AVLTree *AVLTreesUnion(AVLTree *T1, AVLTree *T2);
AVLTree *AVLTreesIntersection(AVLTree *T1, AVLTree *T2);
int InsertNode(AVLTree *T, int k, int v);
int DeleteNode(AVLTree *T, int k, int v);
AVLTreeNode *Search(AVLTree *T, int k, int v);
void FreeAVLTree(AVLTree *T);
void PrintAVLTree(AVLTree *T);

// helper function declarations
void CloneAVLTreeNode(AVLTreeNode *cloneN, AVLTreeNode *N);
AVLTreeNode *AVLTreeFromArrNode(AVLTree *T, AVLTreeNode *P, int K[], int V[], int first, int last);
AVLTree *AVLTreeFromArr(int K[], int V[], int len);
void InorderTraversal(AVLTreeNode *N, int K[], int V[], int *idx);
AVLTreeNode *TriNodeRestructuring(AVLTree *T, AVLTreeNode *N);
int InsertNodeRec(AVLTree *T, AVLTreeNode *N, int k, int v);
int DeleteNodeRec(AVLTree *T, AVLTreeNode *N, int k, int v);
AVLTreeNode *SearchNode(AVLTreeNode *N, int k, int v);
void FreeAVLTreeNode(AVLTreeNode *N);
void PrintAVLTreeNode(AVLTreeNode *N);

// create a new AVL tree by taking input from either user or file
// all elements are scanned in linear O(n) time and stored in array
// subsequently during linear traversal of array (having size n), InsertNode() is called for each key-value pair (taking O(log(n)) time)
// therefore, time complexity - O(n * log(n)), where n is size of resulting AVL tree
AVLTree *CreateAVLTree(const char *filename)
{
    AVLTree *T = newAVLTree();
    int list[SIZE];
    int len = 0;

    if (filename == "stdin") {      // input given by user
        int x = 0;
        int sign = 1;
        int ctr = 0;
        int flagOB = FALSE;         // flag for opening brackets
        int flagCB = FALSE;         // flag for closing brackets
        int flagCM = FALSE;         // flag for comma
        int flagNM = FALSE;         // flag for numbers
        char *buf = (char *)malloc(256 * sizeof(char *));   // allocate memory to user input

        // take input until new line is encountered
        while(scanf("%[^\n]%*c", buf) == 1){
            // parse each line character by character
            for (int i = 0; buf[i] != '\0'; i++) {
                if (buf[i] == '(') {
                    // return error if "(" has "(", ",", or numbers preceding it, or ")" not preceding it (except first "(")
                    if (flagOB || ((ctr) ? !flagCB : flagCB) || flagCM || flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagOB = TRUE;
                    flagCB = FALSE;
                } else if (buf[i] == ',') {
                    // return error if "," has ")" or another "," preceding it, or numbers or "(" not preceding it
                    if (!flagOB || flagCB || flagCM || !flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagCM = TRUE;
                } else if (buf[i] == ')') {
                    // return error if ")" has another ")" preceding it, or "(", "," or numbers not preceding it
                    if (!flagOB || flagCB || !flagCM || flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    ctr++;
                    flagOB = FALSE;
                    flagCB = TRUE;
                    flagCM = FALSE;
                    flagNM = FALSE;
                } else if (buf[i] == '-' || (buf[i] >= '0' && buf[i] <= '9')) {
                    if (buf[i + 1] >= '0' && buf[i + 1] <= '9') {
                        if (buf[i] == '-') {
                            sign = -1;
                        } else {
                            x = x * 10 + (buf[i] - '0');
                        }
                    } else {
                        // store the number in an array and increment its index
                        list[len++] = sign * (x * 10 + (buf[i] - '0'));
                        x = 0;
                        sign = 1;

                        if (flagNM) {
                            // return error if number corresponding to value has ")" preceding it, or "(" or "," not preceding it
                            if (!flagOB || flagCB || !flagCM) {
                                printf("Error: Invalid expression syntax\n");
                                exit(EXIT_FAILURE);
                            }

                            flagNM = FALSE;
                        } else {
                            // return error if number corresponding to key has ")" or "," preceding it, or "(" not preceding it
                            if (!flagOB || flagCB || flagCM) {
                                printf("Error: Invalid expression syntax\n");
                                exit(EXIT_FAILURE);
                            }

                            flagNM = TRUE;
                        }
                    }
                }
            }
        }

        free(buf);      // free memory allocated to user input
        fflush(stdin);  // clear standard input buffer
    } else {
        int x;
        char c;
        int ctr = 0;
        int flagOB = FALSE;     // flag for opening brackets
        int flagCB = FALSE;     // flag for closing brackets
        int flagCM = FALSE;     // flag for commas
        int flagNM = FALSE;     // flag for numbers
        FILE *fp = fopen(filename, "r");    // open file in read mode

        // assertion fails if file does not exist
        assert(fp != NULL);

        // scan each occurrence of integer in file and store it in array until the end of file is encountered
        while (!feof(fp)) {
            if (fscanf(fp, "%d", &x) == 1) {    // if number is encountered in file
                if (flagNM) {
                    // return error if number corresponding to value has ")" preceding it, or "(" or "," not preceding it
                    if (!flagOB || flagCB || !flagCM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagNM = FALSE;
                } else {
                    // return error if number corresponding to key has ")" or "," preceding it, or "(" not preceding it
                    if (!flagOB || flagCB || flagCM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagNM = TRUE;
                }

                list[len++] = x;
            } else {        // if characters other than numbers are encountered in file
                c = fgetc(fp);

                if (c == '(') {
                    // return error if "(" has "(", ",", or numbers preceding it, or ")" not preceding it (except first "(")
                    if (flagOB || ((ctr) ? !flagCB : flagCB) || flagCM || flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagOB = TRUE;
                    flagCB = FALSE;
                } else if (c == ',') {
                    // return error if "," has ")" or another "," preceding it, or numbers or "(" not preceding it
                    if (!flagOB || flagCB || flagCM || !flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    flagCM = TRUE;
                } else if (c == ')') {
                    // return error if ")" has another ")" preceding it, or "(", "," or numbers not preceding it
                    if (!flagOB || flagCB || !flagCM || flagNM) {
                        printf("Error: Invalid expression syntax\n");
                        exit(EXIT_FAILURE);
                    }

                    ctr++;
                    flagOB = FALSE;
                    flagCB = TRUE;
                    flagCM = FALSE;
                    flagNM = FALSE;
                }
            }
        }

        // return error if ")" is not the last character encountered in file
        if (!flagCB) {
            printf("Error: Invalid expression syntax\n");
            exit(EXIT_FAILURE);
        }

        // close file
        fclose(fp);
    }

    // insert each key-value pair into AVL tree as a node
    for (int i = 0; i < len - 1; i += 2) {
        int n = InsertNode(T, list[i], list[i + 1]);
    }

    return T;
}

// recursively clone each node of original tree during preorder traversal
void CloneAVLTreeNode(AVLTreeNode *cloneN, AVLTreeNode *N)
{
    // base case; return if node is NULL
    if (N == NULL) {
        return;
    }

    // set height of cloned node to height of original node
    cloneN -> height = N -> height;

    // if original node has left child, create a new node having its key and value, and set it as left child of cloned node
    if (N -> left) {
        cloneN -> left = newAVLTreeNode(N -> left -> key, N -> left -> value);
        cloneN -> left -> parent = cloneN;
    }

    // if original node has right child, create a new node having its key and value, and set it as right child of cloned node
    if (N -> right) {
        cloneN -> right = newAVLTreeNode(N -> right -> key, N -> right -> value);
        cloneN -> right -> parent = cloneN;
    }

    CloneAVLTreeNode(cloneN -> left, N -> left);    // recursive call to left child of current node
    CloneAVLTreeNode(cloneN -> right, N -> right);  // recursive call to right child of current node
}

// create clone of an AVL tree
// CloneAVLTreeNode() is called here which takes O(n) time for preorder traversal of tree
// therefore, time complexity - O(n), where n is size of T
AVLTree *CloneAVLTree(AVLTree *T)
{
    // create new AVL tree
    AVLTree *cloneT = newAVLTree();

    // if T is empty, return the cloned tree
    if (T -> root == NULL) {
        return cloneT;
    }

    // otherwise, set root of cloned tree to a new node having key and value of the root of the original tree
    cloneT -> root = newAVLTreeNode(T -> root -> key, T -> root -> value);
    // set size of clone tree to size of original tree
    cloneT -> size = T -> size;

    // call recursive function to clone each node of original tree starting from root
    CloneAVLTreeNode(cloneT -> root, T -> root);

    return cloneT;
}

// recursively create an AVL tree from a sorted array (K stores keys and V stores corresponding values)
AVLTreeNode *AVLTreeFromArrNode(AVLTree *T, AVLTreeNode *P, int K[], int V[], int first, int last)
{
    // base case; return if start index becomes greater than end index
    if (first > last) {
        return NULL;
    }

    int mid = (first + last) / 2;   // find floor of mid-point of start and end indices
    AVLTreeNode *N = newAVLTreeNode(K[mid], V[mid]);    // create new node setting its key and value to that of the array indices corresponding to mid

    // if parent is NULL, set current node as root of tree
    if (P == NULL) {
        T -> root = N;
    }

    N -> left = AVLTreeFromArrNode(T, N, K, V, first, mid - 1); // recursive call to set left child of current node as mid-point of start and middle indices of current array slice
    N -> right = AVLTreeFromArrNode(T, N, K, V, mid + 1, last); // recursive call to set right child of current node as mid-point of middle and last indices of current array slice

    N -> parent = P;    // set P as parent of current node
    N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;   // set height of current node as max of the heights of its left and right children, incremented by 1

    // return current node
    return N;
}

// create AVL tree from sorted key and value arrays
AVLTree *AVLTreeFromArr(int K[], int V[], int len)
{
    AVLTree *T = newAVLTree();  // create new AVL tree

    // if AVL tree is not empty, otherwise empty tree is returned
    if (len > 0) {
        T -> root = AVLTreeFromArrNode(T, NULL, K, V, 0, len - 1);  // call recursive function setting initial start and end indices to 0 and (array length - 1) and parent as NULL
        T -> size = len;    // set size of tree to size of sorted array
    }

    return T;
}

// store keys and values of each node in separate arrays
// arrays will be automatically sorted since inorder traversal is being performed
void InorderTraversal(AVLTreeNode *N, int K[], int V[], int *idx)
{
    // base case; return if node is NULL
    if (N == NULL) {
        return;
    }

    InorderTraversal(N -> left, K, V, idx);     // recursive call to left child of current node

    K[*idx] = N -> key;     // set value at key array's current index to key of current node
    V[*idx] = N -> value;   // set value at value array's current index to value of current node
    (*idx)++;               // increment index

    InorderTraversal(N -> right, K, V, idx);    // recursive call to right child of current node
}

// the inorder traversals take O(m) and O(n) time respectively
// the merging of the two arrays takes O(m + n) time, since both arrays are being linearly traversed simultaneously
// finally, AVLTreeFromArr() is called, which takes O(m + n) time, since all elements of each combined array (key and value, respectively) are visited once
// therefore, time complexity - O(m + n), where m is size of T1 and n is size of T2
AVLTree *AVLTreesUnion(AVLTree *T1, AVLTree *T2)
{
    int idx;

    idx = 0;                                    // set index to 0 initially
    int *K1 = malloc(T1 -> size * sizeof(int)); // create array to store keys for tree T1
    int *V1 = malloc(T1 -> size * sizeof(int)); // create array to store values for tree T1
    InorderTraversal(T1 -> root, K1, V1, &idx); // perform inorder traversal of T1 so that final arrays are sorted

    idx = 0;                                    // reset index to 0 initially
    int *K2 = malloc(T2 -> size * sizeof(int)); // create array to store keys for tree T2
    int *V2 = malloc(T2 -> size * sizeof(int)); // create array to store values for tree T2
    InorderTraversal(T2 -> root, K2, V2, &idx); // perform inorder traversal of T2 so that final arrays are sorted

    idx = 0;                                    // reset index to 0 initially
    int *Ku = malloc((T1 -> size + T2 -> size) * sizeof(int));  // create array to store resultant keys after union of T1 and T2
    int *Vu = malloc((T1 -> size + T2 -> size) * sizeof(int));  // create array to store resultant values after union of T1 and T2

    int i1 = 0;     // initial index for K1 and V1
    int i2 = 0;     // initial index for K2 and V2
    int len = 0;    // initial index for Ku and Vu

    // both (K1, V1) and (K2, V2) arrays have elements remaining
    while (i1 < T1 -> size && i2 < T2 -> size) {
        if (K1[i1] < K2[i2] || (K1[i1] == K2[i2] && V1[i1] < V2[i2])) {
            // if (K1, V1) < (K2, V2) at current index, append (K1, V1) to (Ku, Vu) and move (K1, V1) to its next index
            Ku[len] = K1[i1];
            Vu[len] = V1[i1];
            i1++;
        } else if (K1[i1] > K2[i2] || (K1[i1] == K2[i2] && V1[i1] > V2[i2])) {
            // if (K1, V1) > (K2, V2) at current index, append (K2, V2) to (Ku, Vu) and move (K2, V2) to its next index
            Ku[len] = K2[i2];
            Vu[len] = V2[i2];
            i2++;
        } else {
            // if (K1, V1) == (K2, V2) at current index, append (K2, V2) to (Ku, Vu) and move both (K1, V1) and (K2, V2) to their next indices
            Ku[len] = K2[i1];
            Vu[len] = V2[i1];
            i1++;
            i2++;
        }

        // increment index of (Ku, Vu)
        len++;
    }

    // only (K1, V1) arrays have elements remaining
    while (i1 < T1 -> size) {
        // append remaining (K1, V1) to (Ku, Vu)
        Ku[len] = K1[i1];
        Vu[len] = V1[i1];
        i1++;
        len++;  // increment index of (Ku, Vu)
    }

    // only (K2, V2) arrays have elements // append remaining (K1, V1) to (Ku, Vu)remaining
    while (i2 < T2 -> size) {
        // append remaining (K2, V2) to (Ku, Vu)
        Ku[len] = K2[i2];
        Vu[len] = V2[i2];
        i2++;
        len++;  // increment index of (Ku, Vu)
    }

    return AVLTreeFromArr(Ku, Vu, len);
}

// the inorder traversals take O(m) and O(n) time respectively
// the merging of the two arrays takes O(m + n) time, since both arrays are being linearly traversed simultaneously
// finally, AVLTreeFromArr() is called, which takes O(m + n) time, since all elements of each combined array (key and value, respectively) are visited once
// therefore, time complexity - O(m + n), where m is size of T1 and n is size of T2
AVLTree *AVLTreesIntersection(AVLTree *T1, AVLTree *T2)
{
    int idx;

    idx = 0;                                    // set index to 0 initially
    int *K1 = malloc(T1 -> size * sizeof(int)); // create array to store keys for tree T1
    int *V1 = malloc(T1 -> size * sizeof(int)); // create array to store values for tree T1
    InorderTraversal(T1 -> root, K1, V1, &idx); // perform inorder traversal of T1 so that final arrays are sorted

    idx = 0;                                    // reset index to 0 initially
    int *K2 = malloc(T2 -> size * sizeof(int)); // create array to store keys for tree T2
    int *V2 = malloc(T2 -> size * sizeof(int)); // create array to store values for tree T2
    InorderTraversal(T2 -> root, K2, V2, &idx); // perform inorder traversal of T2 so that final arrays are sorted

    idx = 0;                                    // reset index to 0 initially
    int *Ki = malloc((T1 -> size + T2 -> size) * sizeof(int));  // create array to store resultant keys after intersection of T1 and T2
    int *Vi = malloc((T1 -> size + T2 -> size) * sizeof(int));  // create array to store resultant values after intersection of T1 and T2

    int i1 = 0;     // initial index for K1 and V1
    int i2 = 0;     // initial index for K2 and V2
    int len = 0;    // initial index for Ki and Vi

    // both (K1, V1) and (K2, V2) arrays have elements remaining
    while (i1 < T1 -> size && i2 < T2 -> size) {
        if (K1[i1] < K2[i2] || (K1[i1] == K2[i2] && V1[i1] < V2[i2])) {
            // if (K1, V1) < (K2, V2) at current index, move (K1, V1) to its next index
            i1++;
        } else if (K1[i1] > K2[i2] || (K1[i1] == K2[i2] && V1[i1] > V2[i2])) {
            // if (K1, V1) > (K2, V2) at current index, move (K2, V2) to its next index
            i2++;
        } else {
            // if (K1, V1) == (K2, V2) at current index, append (K2, V2) to (Ki, Vi) and move both (K1, V1) and (K2, V2) to their next indices
            Ki[len] = K2[i1];
            Vi[len] = V2[i1];
            i1++;
            i2++;
            len++;      // increment index of (Ki, Vi)
        }
    }

    return AVLTreeFromArr(Ki, Vi, len);
}

// traverse to first ancestor node after insertion or deletion operation has been performed, whose left and right subtrees' heights differ by more than 1
// restructure the subtree rooted at that node in order to balance the tree
// at a particular node, this operation takes O(1) time
AVLTreeNode *TriNodeRestructuring(AVLTree *T, AVLTreeNode *N)
{
    AVLTreeNode *z;     // first ancestor node with unbalanced subtrees
    AVLTreeNode *y;     // child node of z with greater height
    AVLTreeNode *x;     // child node of y with greater height

    if (HEIGHT(N -> left) > HEIGHT(N -> right)) {
        // left-left single rotation case
        // y is left child of z and x is left child of y
        if (HEIGHT(N -> left -> left) > HEIGHT(N -> left -> right)) {
            z = N;
            y = N -> left;

            z -> left = y -> right;     // update left child of z to point to right child of y
            y -> right = z;             // update right child of y to point to z

            // update parents of z and y
            if (z -> parent) {
                if (z -> parent -> left == z) {
                    z -> parent -> left = y;
                } else {
                    z -> parent -> right = y;
                }
            }

            y -> parent = z -> parent;
            z -> parent = y;

            if (z -> left) {
                z -> left -> parent = z;
            }

            // update heights of z and y
            z -> height = MAX(HEIGHT(z -> left), HEIGHT(z -> right)) + 1;
            y -> height = MAX(HEIGHT(y -> left), HEIGHT(y -> right)) + 1;

            // if z was root of tree initially, reset root of tree to y
            if (z == T -> root) {
                T -> root = y;
            }

            N = y;
            N -> right = z;

            return N;
        }

        // left-right double rotation case
        // y is left child of z and x is right child of y
        z = N;
        y = N -> left;
        x = N -> left -> right;

        z -> left = x -> right;     // update left child of z to point to right child of x
        y -> right = x -> left;     // update right child of y to point to left child of x
        x -> left = y;              // update left child of x to point to y
        x -> right = z;             // update right child of x to point to z

        // update parents of z, y and x
        if (z -> parent) {
            if (z -> parent -> left == z) {
                z -> parent -> left = x;
            } else {
                z -> parent -> right = x;
            }
        }

        x -> parent = z -> parent;
        z -> parent = x;
        y -> parent = x;

        if (z -> left) {
            z -> left -> parent = z;
        }

        if (y -> right) {
            y -> right -> parent = y;
        }

        // update heights of z, y and x
        z -> height = MAX(HEIGHT(z -> left), HEIGHT(z -> right)) + 1;
        y -> height = MAX(HEIGHT(y -> left), HEIGHT(y -> right)) + 1;
        x -> height = MAX(HEIGHT(x -> left), HEIGHT(x -> right)) + 1;

        // if z was root of tree initially, reset root of tree to x
        if (z == T -> root) {
            T -> root = x;
        }

        N = x;
        N -> left = y;
        N -> right = z;

        return N;
    }

    // right-left double rotation case
    // y is right child of z and x is left child of y
    if (HEIGHT(N -> right -> left) > HEIGHT(N -> right -> right)) {
        z = N;
        y = N -> right;
        x = N -> right -> left;

        z -> right = x -> left;     // update right child of z to point to left child of x
        y -> left = x -> right;     // update left child of y to point to right child of x
        x -> left = z;              // update left child of x to point to z
        x -> right = y;             // update right child of x to point to y

        // update parents of z, y and x
        if (z -> parent) {
            if (z -> parent -> left == z) {
                z -> parent -> left = x;
            } else {
                z -> parent -> right = x;
            }
        }

        x -> parent = z -> parent;
        z -> parent = x;
        y -> parent = x;

        if (z -> right) {
            z -> right -> parent = z;
        }

        if (y -> left) {
            y -> left -> parent = y;
        }

        // update heights of z, y and x
        z -> height = MAX(HEIGHT(z -> left), HEIGHT(z -> right)) + 1;
        y -> height = MAX(HEIGHT(y -> left), HEIGHT(y -> right)) + 1;
        x -> height = MAX(HEIGHT(x -> left), HEIGHT(x -> right)) + 1;

        // if z was root of tree initially, reset root of tree to x
        if (z == T -> root) {
            T -> root = x;
        }

        N = x;
        N -> left = z;
        N -> right = y;

        return N;
    }

    // right-right single rotation case
    // y is right child of z and x is right child of y
    z = N;
    y = N -> right;

    z -> right = y -> left;     // update right child of z to point to left child of y
    y -> left = z;              // update left child of y to point to z

    // update parents of z and y
    if (z -> parent) {
        if (z -> parent -> left == z) {
            z -> parent -> left = y;
        } else {
            z -> parent -> right = y;
        }
    }

    y -> parent = z -> parent;
    z -> parent = y;

    if (z -> right) {
        z -> right -> parent = z;
    }

    // update heights of z and y
    z -> height = MAX(HEIGHT(z -> left), HEIGHT(z -> right)) + 1;
    y -> height = MAX(HEIGHT(y -> left), HEIGHT(y -> right)) + 1;

    // if z was root of tree initially, reset root of tree to y
    if (z == T -> root) {
        T -> root = y;
    }

    N = y;
    N -> left = z;

    return N;
}

// recursively search for position to insert node and then perform insertion and tri-node restructuring (if required)
int InsertNodeRec(AVLTree *T, AVLTreeNode *N, int k, int v)
{
    // return 0 if desired key-value pair is equal to key-value pair of current node, i.e., node is already present in tree
    if (k == N -> key && v == N -> value) {
        return 0;
    }

    int flag;

    if (k < N -> key || (k == N -> key && v < N -> value)) {
        // desired key-value pair is less than key-value pair of current node
        if (N -> left) {
            // current node has a left child
            flag = InsertNodeRec(T, N -> left, k, v);   // recursive call in left child of current node

            // if node was found in subtree rooted at N, check if tri-node restructuring is required at N
            if (flag) {
                // tri-node restructuring is performed if difference between heights of left and right subtrees of N is greater than 1
                if (!(MAX(HEIGHT(N -> left), HEIGHT(N -> right)) - MIN(HEIGHT(N -> left), HEIGHT(N -> right)) <= 1)) {
                    N = TriNodeRestructuring(T, N);
                }

                // update height of N after restructuring (if performed)
                N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;
            }

            return flag;
        } else {
            // current node has no left child, therefore insertion can be performed at that position
            N -> left = newAVLTreeNode(k, v);                               // set left child of N as new node with key-value pair (k, v)
            N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;   // update height of N
            N -> left -> parent = N;                                        // set current node as parent of inserted node
            T -> size += 1;                                                 // increment size of tree by 1
        }
    } else if (k > N -> key || (k == N -> key && v > N -> value)) {
        // desired key-value pair is greater than key-value pair of current node
        if (N -> right) {
            // current node has a right child
            flag = InsertNodeRec(T, N -> right, k, v);      // recursive call in right child of current node

            // if node was found in subtree rooted at N, check if tri-node restructuring is required at N
            if (flag) {
                // tri-node restructuring is performed if difference between heights of left and right subtrees of N is greater than 1
                if (!(MAX(HEIGHT(N -> left), HEIGHT(N -> right)) - MIN(HEIGHT(N -> left), HEIGHT(N -> right)) <= 1)) {
                    N = TriNodeRestructuring(T, N);
                }

                // update height of N after restructuring (if performed)
                N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;
            }

            return flag;
        } else {
            // current node has no right child, therefore insertion can be performed at that position
            N -> right = newAVLTreeNode(k, v);                              // set right child of N as new node with key-value pair (k, v)
            N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;   // update height of N
            N -> right -> parent = N;                                       // set current node as parent of inserted node
            T -> size += 1;                                                 // increment size of tree by 1
        }
    }

    // return 1 if insertion is performed successfully, i.e., node is not already present in tree
    return 1;
}

// InsertNodeRec() takes O(log(n)) time, since in each recursion, it either calls its left child or right child while searching until the position to insert the node is found, eliminating the other subtree from its search
// TriNodeRestructuring() is called inside InsertNodeRec(), which takes O(1) time
// therefore, time complexity - O(log(n)), where n is size of T
int InsertNode(AVLTree *T, int k, int v)
{
    // if tree is empty, create new node with (k, v) as its key-value pair and set it as root of tree
    if (T -> root == NULL) {
        T -> root = newAVLTreeNode(k, v);
        T -> size += 1;     // increment size of tree by 1

        // return 1 since insertion is successfully performed
        return 1;
    }

    // call recursive function to search for position to insert node
    return InsertNodeRec(T, T -> root, k, v);
}

// recursively search for and delete the node with the corresponding key-value pair and perform tri-node restructuring at each ancestor node (if required)
int DeleteNodeRec(AVLTree *T, AVLTreeNode *N, int k, int v)
{
    // return 0 if the node is not present in tree
    if (N == NULL) {
        return 0;
    }

    int flag;

    if (k < N -> key || (k == N -> key && v < N -> value)) {
        // desired key-value pair is less than key-value pair of current node
        flag = DeleteNodeRec(T, N -> left, k, v);   // recursive call to search in left child of current node

        // if node was found in subtree rooted at N, check if tri-node restructuring is required at N
        if (flag) {
            // tri-node restructuring is performed if difference between heights of left and right subtrees of N is greater than 1
            if (!(MAX(HEIGHT(N -> left), HEIGHT(N -> right)) - MIN(HEIGHT(N -> left), HEIGHT(N -> right)) <= 1)) {
                N = TriNodeRestructuring(T, N);
            }

            // update height of N after restructuring (if performed)
            N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;
        }

        return flag;
    } else if (k > N -> key || (k == N -> key && v > N -> value)) {
        // desired key-value pair is greater than key-value pair of current node
        flag = DeleteNodeRec(T, N -> right, k, v);  // recursive call to search in right child of current node

        // if node was found in subtree rooted at N, check if tri-node restructuring is required at N
        if (flag) {
            // tri-node restructuring is performed if difference between heights of left and right subtrees of N is greater than 1
            if (!(MAX(HEIGHT(N -> left), HEIGHT(N -> right)) - MIN(HEIGHT(N -> left), HEIGHT(N -> right)) <= 1)) {
                N = TriNodeRestructuring(T, N);
            }

            // update height of N after restructuring (if performed)
            N -> height = MAX(HEIGHT(N -> left), HEIGHT(N -> right)) + 1;
        }

        return flag;
    } else {
        // desired key-value pair is equal to key-value pair of current node, i.e., node is present in tree
        if (N -> left == NULL && N -> right == NULL) {
            // N has no child nodes
            if (N -> parent) {
                // if N is not root node, set its parent's pointer to it as NULL
                if (N -> parent -> left == N) {
                    N -> parent -> left = NULL;
                } else {
                    N -> parent -> right = NULL;
                }
            } else {
                // if N is root node, set root of tree to NULL
                T -> root = NULL;
            }

            // free space allocated to N and decrement size of tree by 1
            free(N);
            T -> size -= 1;
        } else if (N -> left != NULL && N -> right == NULL) {
            // N only has left child node
            if (N -> parent) {
                // if N is not root node, update its parent and left child as pointing to each other
                N -> left -> parent = N -> parent;

                if (N -> parent -> left == N) {
                    N -> parent -> left = N -> left;
                } else {
                    N -> parent -> right = N -> left;
                }
            } else {
                // if N is root node, set root of tree to left child of N and set its parent as NULL
                N -> left -> parent = NULL;
                T -> root = N -> left;
            }

            // free space allocated to N and decrement size of tree by 1
            free(N);
            T -> size -= 1;
        } else if (N -> left == NULL && N -> right != NULL) {
            // N only has right child node
            if (N -> parent) {
                // if N is not root node, update its parent and right child as pointing to each other
                N -> right -> parent = N -> parent;

                if (N -> parent -> left == N) {
                    N -> parent -> left = N -> right;
                } else {
                    N -> parent -> right = N -> right;
                }
            } else {
                // if N is root node, set root of tree to right child of N and set its parent as NULL
                N -> right -> parent = NULL;
                T -> root = N -> right;
            }

            // free space allocated to N and decrement size of tree by 1
            free(N);
            T -> size -= 1;
        } else {
            // N has both left and right child nodes
            // in this case, inorder successor of N is found by traversing to the leftmost node in its right subtree
            AVLTreeNode *IS = N -> right;

            while (IS -> left) {
                IS = IS -> left;
            }

            // update key and value of N to key and value of inorder successor of N
            N -> key = IS -> key;
            N -> value = IS -> value;

            // call recursive function to delete inorder successor node from tree
            return DeleteNodeRec(T, IS, IS -> key, IS -> value);
        }
    }

    // return 1 if node to be deleted is present in tree
    return 1;
}

// DeleteNodeRec() takes O(log(n)) time, since in each recursion, it either calls its left child or right child while searching until the node to be deleted is found, eliminating the other subtree from its search
// TriNodeRestructuring() is called inside DeleteNodeRec(), which takes O(1) time
// therefore, time complexity - O(log(n)), where n is size of T
int DeleteNode(AVLTree *T, int k, int v)
{
    // call recursive function to search for node that needs to be deleted
    return DeleteNodeRec(T, T -> root, k, v);
}

// recursively search for desired node
AVLTreeNode *SearchNode(AVLTreeNode *N, int k, int v)
{
    // return NULL if node is not found
    if (N == NULL) {
        return NULL;
    }

    // search in left child of current node if desired key-value pair is less than current node's key-value pair
    if (k < N -> key || (k == N -> key && v < N -> value)) {
        return SearchNode(N -> left, k, v);
    }

    // search in right child of current node if desired key-value pair is greater than current node's key-value pair
    if (k > N -> key || (k == N -> key && v > N -> value)) {
        return SearchNode(N -> right, k, v);
    }

    // return current node if desired key-value pair is equal to current node's key-value pair
    return N;
}

// free memory space allocated to AVL tree
// SearchNode() is called which takes O(log(n)) time, since in each recursion, it either calls its left child or right child until the desired search is found, eliminating the other subtree from its search
// therefore, time complexity - O(log(n)), where n is size of T
AVLTreeNode *Search(AVLTree *T, int k, int v)
{
    // call recursive function to search for node with desired key-value pair starting from root; return NULL if no such node exists
    return SearchNode(T -> root, k, v);
}

// recursively free space allocated to each node during postorder traversal (since child nodes should be freed before the parent node)
void FreeAVLTreeNode(AVLTreeNode *N)
{
    // base case; return if node is NULL
    if (N == NULL) {
        return;
    }

    FreeAVLTreeNode(N -> left);     // recursive call to left child of current node
    FreeAVLTreeNode(N -> right);    // recursive call to right child of current node
    free(N);                        // free space allocated to current node structure
}

// free memory space allocated to AVL tree
// FreeAVLTreeNode() is called which takes O(n) time for postorder traversal of tree
// therefore, time complexity - O(n), where n is size of T
void FreeAVLTree(AVLTree *T)
{
    // call recursive function to free space allocated to each node
    FreeAVLTreeNode(T -> root);

    // finally, free memory allocated to tree structure itself
    free(T);
}

// recursively print each node's items during inorder traversal
void PrintAVLTreeNode(AVLTreeNode *N)
{
    // base case; return if node is NULL
    if (N == NULL) {
        return;
    }

    PrintAVLTreeNode(N -> left);    // recursive call to left child of current node
    printf("(%d, %d), %d\n", N -> key, N -> value, N -> height);    // print key, value, and height of current node
    PrintAVLTreeNode(N -> right);   // recursive call to left child of current node
}

// print key, value, and height for each node in AVL tree
// PrintAVLTreeNode() is called which takes O(n) time for inorder traversal of tree
// therefore, time complexity - O(n), where n is size of T
void PrintAVLTree(AVLTree *T)
{
    // call recursive function to print each node
    PrintAVLTreeNode(T -> root);

    printf("\n");
}

int main() //sample main for testing
{
    int i, j;
    AVLTree *tree1, *tree2, *tree3, *tree4, *tree5, *tree6, *tree7, *tree8;
    AVLTreeNode *node1;

    tree1 = CreateAVLTree("stdin");
    PrintAVLTree(tree1);

    FreeAVLTree(tree1);

    // you need to create the text file file1.txt
    // to store a set of items without duplicate items
    tree2 = CreateAVLTree("file1.txt");
    PrintAVLTree(tree2);

    tree3 = CloneAVLTree(tree2);
    PrintAVLTree(tree3);

    FreeAVLTree(tree2);
    FreeAVLTree(tree3);

    //Create tree4
    tree4 = newAVLTree();

    j = InsertNode(tree4, 10, 10);
    for (i = 0; i < 15; i++) {
        j = InsertNode(tree4, i, i);

        if (j == 0) {
            printf("(%d, %d) already exists\n", i, i);
        }
    }

    PrintAVLTree(tree4);

    node1 = Search(tree4, 20, 20);

    if (node1 != NULL) {
        printf("key= %d value= %d\n", node1 -> key, node1 -> value);
    } else {
        printf("Key 20 does not exist\n");
    }

    for (i = 17; i > 0; i--) {
        j = DeleteNode(tree4, i, i);

        if (j == 0) {
            printf("Key %d does not exist\n", i);
        }

        PrintAVLTree(tree4);
    }

    FreeAVLTree(tree4);

    //Create tree5
    tree5 = newAVLTree();

    j = InsertNode(tree5, 6, 25);
    j = InsertNode(tree5, 6, 10);
    j = InsertNode(tree5, 6, 12);
    j = InsertNode(tree5, 6, 20);
    j = InsertNode(tree5, 9, 25);
    j = InsertNode(tree5, 10, 25);

    PrintAVLTree(tree5);

    //Create tree6
    tree6 = newAVLTree();

    j = InsertNode(tree6, 6, 25);
    j = InsertNode(tree6, 5, 10);
    j = InsertNode(tree6, 6, 12);
    j = InsertNode(tree6, 6, 20);
    j = InsertNode(tree6, 8, 35);
    j = InsertNode(tree6, 10, 25);

    PrintAVLTree(tree6);

    tree7 = AVLTreesIntersection(tree5, tree6);
    tree8 = AVLTreesUnion(tree5, tree6);

    PrintAVLTree(tree7);
    PrintAVLTree(tree8);

    return 0;
}
