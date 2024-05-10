quick sort
#include <stdio.h>

void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

int partition(int array[], int left, int right) {
 
  int pivot = array[right];
 
  int i = (left - 1);

  for (int j = left; j < right; j++) {
    if (array[j] <= pivot) {
       
      i++;
     
      swap(&array[i], &array[j]);
    }
  }

  swap(&array[i + 1], &array[right]);
 
  return (i + 1);
}

void quickSort(int array[], int left, int right) {
  if (left < right) {
   
    int pi = partition(array, left, right);
   
    quickSort(array, left, pi - 1);
   
    quickSort(array, pi + 1, right);
  }
}

void printArray(int array[], int size) {
  for (int i = 0; i < size; ++i) {
    printf("%d  ", array[i]);
  }
  printf("\n");
}

int main() {
  int data[] = {3,1,2,5,4,7};
 
  int n = sizeof(data) / sizeof(data[0]);
 
  printf("Unsorted Array\n");
  printArray(data, n);
 
  quickSort(data, 0, n - 1);
 
  printf("Sorted array in ascending order: \n");
  printArray(data, n);
}

merge sort
#include<stdio.h>
void merge(int arr[], int p, int q, int r) {

    int n1 = q - p + 1;
    int n2 = r - q;

    int L[n1], M[n2];

    for (int i = 0; i < n1; i++)
        L[i] = arr[p + i];
    for (int j = 0; j < n2; j++)
        M[j] = arr[q + 1 + j];

    int i, j, k;
    i = 0;
    j = 0;
    k = p;

    while (i < n1 && j < n2) {
        if (L[i] <= M[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = M[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = M[j];
        j++;
        k++;
    }
}
void main()
{
   
}

topologocal sort
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

struct Stack {
int data;
struct Stack* next;
};

struct Graph {
int V;
struct List* adj;
};

struct List {
int data;
struct List* next;
};

struct Stack* createStackNode(int data)
{
struct Stack* newNode
= (struct Stack*)malloc(sizeof(struct Stack));
newNode->data = data;
newNode->next = NULL;
return newNode;
}

struct List* createListNode(int data)
{
struct List* newNode
= (struct List*)malloc(sizeof(struct List));
newNode->data = data;
newNode->next = NULL;
return newNode;
}

struct Graph* createGraph(int V)
{
struct Graph* graph
= (struct Graph*)malloc(sizeof(struct Graph));
graph->V = V;
graph->adj
= (struct List*)malloc(V * sizeof(struct List));
for (int i = 0; i < V; ++i) {
graph->adj[i].next = NULL;
}
return graph;
}

void addEdge(struct Graph* graph, int v, int w)
{
struct List* newNode = createListNode(w);
newNode->next = graph->adj[v].next;
graph->adj[v].next = newNode;
}

void topologicalSortUtil(struct Graph* graph, int v,
bool visited[],
struct Stack** stack)
{
visited[v] = true;

struct List* current = graph->adj[v].next;
while (current != NULL) {
int adjacentVertex = current->data;
if (!visited[adjacentVertex]) {
topologicalSortUtil(graph, adjacentVertex,
visited, stack);
}
current = current->next;
}

struct Stack* newNode = createStackNode(v);
newNode->next = *stack;
*stack = newNode;
}

void topologicalSort(struct Graph* graph)
{
struct Stack* stack = NULL;

bool* visited = (bool*)malloc(graph->V * sizeof(bool));
for (int i = 0; i < graph->V; ++i) {
visited[i] = false;
}

for (int i = 0; i < graph->V; ++i) {
if (!visited[i]) {
topologicalSortUtil(graph, i, visited, &stack);
}
}

while (stack != NULL) {
printf("%d ", stack->data);
struct Stack* temp = stack;
stack = stack->next;
free(temp);
}

free(visited);
free(graph->adj);
free(graph);
}

int main()
{
struct Graph* g = createGraph(6);
addEdge(g, 1, 2);
addEdge(g, 1, 3);
addEdge(g, 1, 5);
addEdge(g, 3, 4);
addEdge(g, 3, 5);
addEdge(g, 4, 5);

printf("Topological Sorting Order: ");
topologicalSort(g);

return 0;
}

##knapsack 
#include<stdio.h>
int max(int a, int b) {
   if(a>b){
      return a;
   } else {
      return b;
   }
}
int knapsack(int W, int wt[], int val[], int n) {
   int i, w;
   int knap[n+1][W+1];
   for (i = 0; i <= n; i++) {
      for (w = 0; w <= W; w++) {
         if (i==0 || w==0)
            knap[i][w] = 0;
         else if (wt[i-1] <= w)
            knap[i][w] = max(val[i-1] + knap[i-1][w-wt[i-1]], knap[i-1][w]);
         else
            knap[i][w] = knap[i-1][w];
      }
 
   }
   return knap[n][W];
}
int main() {
   int val[] = {20, 25, 40};
   int wt[] = {25, 20, 30};
   int W = 50;
   int n = sizeof(val)/sizeof(val[0]);
   printf("The solution is : %d", knapsack(W, wt, val, n));
   return 0;
}

bfs 
#include <stdio.h>
#include <stdlib.h>
#define SIZE 40

struct queue {
  int items[SIZE];
  int front;
  int rear;
};

struct queue* createQueue();
void enqueue(struct queue* q, int);
int dequeue(struct queue* q);
void display(struct queue* q);
int isEmpty(struct queue* q);
void printQueue(struct queue* q);

struct node {
  int vertex;
  struct node* next;
};

struct node* createNode(int);

struct Graph {
  int numVertices;
  struct node** adjLists;
  int* visited;
};

void bfs(struct Graph* graph, int startVertex) {
  struct queue* q = createQueue();

  graph->visited[startVertex] = 1;
  enqueue(q, startVertex);

  while (!isEmpty(q)) {
    printQueue(q);
    int currentVertex = dequeue(q);
    printf("Visited %d\n", currentVertex);

    struct node* temp = graph->adjLists[currentVertex];

    while (temp) {
      int adjVertex = temp->vertex;

      if (graph->visited[adjVertex] == 0) {
        graph->visited[adjVertex] = 1;
        enqueue(q, adjVertex);
      }
      temp = temp->next;
    }
  }
}

struct node* createNode(int v) {
  struct node* newNode = malloc(sizeof(struct node));
  newNode->vertex = v;
  newNode->next = NULL;
  return newNode;
}

struct Graph* createGraph(int vertices) {
  struct Graph* graph = malloc(sizeof(struct Graph));
  graph->numVertices = vertices;

  graph->adjLists = malloc(vertices * sizeof(struct node*));
  graph->visited = malloc(vertices * sizeof(int));

  int i;
  for (i = 0; i < vertices; i++) {
    graph->adjLists[i] = NULL;
    graph->visited[i] = 0;
  }

  return graph;
}

void addEdge(struct Graph* graph, int src, int dest) {
  struct node* newNode = createNode(dest);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;

  newNode = createNode(src);
  newNode->next = graph->adjLists[dest];
  graph->adjLists[dest] = newNode;
}

struct queue* createQueue() {
  struct queue* q = malloc(sizeof(struct queue));
  q->front = -1;
  q->rear = -1;
  return q;
}

int isEmpty(struct queue* q) {
  if (q->rear == -1)
    return 1;
  else
    return 0;
}

void enqueue(struct queue* q, int value) {
  if (q->rear == SIZE - 1)
    printf("\nQueue is Full!!");
  else {
    if (q->front == -1)
      q->front = 0;
    q->rear++;
    q->items[q->rear] = value;
  }
}

int dequeue(struct queue* q) {
  int item;
  if (isEmpty(q)) {
    printf("Queue is empty");
    item = -1;
  } else {
    item = q->items[q->front];
    q->front++;
    if (q->front > q->rear) {
      printf("Resetting queue ");
      q->front = q->rear = -1;
    }
  }
  return item;
}

void printQueue(struct queue* q) {
  int i = q->front;

  if (isEmpty(q)) {
    printf("Queue is empty");
  } else {
    printf("\nQueue contains \n");
    for (i = q->front; i < q->rear + 1; i++) {
      printf("%d ", q->items[i]);
    }
  }
}

int main() {
  struct Graph* graph = createGraph(6);
  addEdge(graph, 0, 1);
  addEdge(graph, 0, 2);
  addEdge(graph, 1, 2);
  addEdge(graph, 1, 4);
  addEdge(graph, 1, 3);
  addEdge(graph, 2, 4);
  addEdge(graph, 3, 4);

  bfs(graph, 0);

  return 0;
}

dfs 
#include <stdio.h>
#include <stdlib.h>

struct node {
  int vertex;
  struct node* next;
};

struct node* createNode(int v);

struct Graph {
  int numVertices;
  int* visited;

  struct node** adjLists;
};

void DFS(struct Graph* graph, int vertex) {
  struct node* adjList = graph->adjLists[vertex];
  struct node* temp = adjList;

  graph->visited[vertex] = 1;
  printf("Visited %d \n", vertex);

  while (temp != NULL) {
    int connectedVertex = temp->vertex;

    if (graph->visited[connectedVertex] == 0) {
      DFS(graph, connectedVertex);
    }
    temp = temp->next;
  }
}

struct node* createNode(int v) {
  struct node* newNode = malloc(sizeof(struct node));
  newNode->vertex = v;
  newNode->next = NULL;
  return newNode;
}

struct Graph* createGraph(int vertices) {
  struct Graph* graph = malloc(sizeof(struct Graph));
  graph->numVertices = vertices;

  graph->adjLists = malloc(vertices * sizeof(struct node*));

  graph->visited = malloc(vertices * sizeof(int));

  int i;
  for (i = 0; i < vertices; i++) {
    graph->adjLists[i] = NULL;
    graph->visited[i] = 0;
  }
  return graph;
}

void addEdge(struct Graph* graph, int src, int dest) {
  struct node* newNode = createNode(dest);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;

  newNode = createNode(src);
  newNode->next = graph->adjLists[dest];
  graph->adjLists[dest] = newNode;
}

void printGraph(struct Graph* graph) {
  int v;
  for (v = 0; v < graph->numVertices; v++) {
    struct node* temp = graph->adjLists[v];
    printf("\n Adjacency list of vertex %d\n ", v);
    while (temp) {
      printf("%d -> ", temp->vertex);
      temp = temp->next;
    }
    printf("\n");
  }
}

int main() {
  struct Graph* graph = createGraph(4);
  addEdge(graph, 0, 1);
  addEdge(graph, 0, 2);
  addEdge(graph, 1, 2);
  addEdge(graph, 2, 3);

  printGraph(graph);

  DFS(graph, 2);

  return 0;
}

prims 
#include <stdio.h>
#include <limits.h>

#define max_vertices 100

int minkey(int key[], int mstSet[], int vertices) {
    int min = int_max, minindex;

    for (int v = 0; v < vertices; v++) {
        if (!mstset[v] && key[v] < min) {
            min = key[v];
            minIndex = v;
        }
    }

    return minindex;
}

void printmst(int parent[], int graph[max_vertices][max_vertices], int vertices) {
    printf("edge \tweight\n");
    for (int i = 1; i < vertices; i++) {
        printf("%d - %d \t%d\n", parent[i], i, graph[i][parent[i]]);
    }
}

void primmst(int graph[max_vertices][max_vertices], int vertices) {
    int parent[max_vertices];
    int key[max_vertices];    
    int mstSet[max_vertices];

    for (int i = 0; i < vertices; i++) {
        key[i] = int_max;
        mstset[i] = 0;
    }

    key[0] = 0;
    parent[0] = -1;
    for (int count = 0; count < vertices - 1; count++) {
        int u = minkey(key, mstset, vertices);

        mstset[u] = 1;

        for (int v = 0; v < vertices; v++) {
            if (graph[u][v] && !mstset[v] && graph[u][v] < key[v]) {
                parent[v] = u;
                key[v] = graph[u][v];
            }
        }
    }

    printmst(parent, graph, vertices);
}

int main() {
    int vertices;

    printf("input the number of vertices: ");
    scanf("%d", &vertices);

    if (vertices <= 0 || vertices > max_vertices) {
        printf("invalid number of vertices. exiting.\n");
        return 1;
    }

    int graph[max_vertices][max_vertices];

    printf("input the adjacency matrix for the graph:\n");
    for (int i = 0; i < vertices; i++) {
        for (int j = 0; j < vertices; j++) {
            scanf("%d", &graph[i][j]);
        }
    }

    primmst(graph, vertices);

    return 0;
}

kruskal 
#include <stdio.h>
#include <stdlib.h>

int comparator(const void* p1, const void* p2)
{
const int(*x)[3] = p1;
const int(*y)[3] = p2;

return (*x)[2] - (*y)[2];
}

void makeSet(int parent[], int rank[], int n)
{
for (int i = 0; i < n; i++) {
parent[i] = i;
rank[i] = 0;
}
}

int findParent(int parent[], int component)
{
if (parent[component] == component)
return component;

return parent[component]
= findParent(parent, parent[component]);
}

void unionSet(int u, int v, int parent[], int rank[], int n)
{
u = findParent(parent, u);
v = findParent(parent, v);

if (rank[u] < rank[v]) {
parent[u] = v;
}
else if (rank[u] > rank[v]) {
parent[v] = u;
}
else {
parent[v] = u;

rank[u]++;
}
}

void kruskalAlgo(int n, int edge[n][3])
{
qsort(edge, n, sizeof(edge[0]), comparator);

int parent[n];
int rank[n];

makeSet(parent, rank, n);

int minCost = 0;

printf(
"Following are the edges in the constructed MST\n");
for (int i = 0; i < n; i++) {
int v1 = findParent(parent, edge[i][0]);
int v2 = findParent(parent, edge[i][1]);
int wt = edge[i][2];

if (v1 != v2) {
unionSet(v1, v2, parent, rank, n);
minCost += wt;
printf("%d -- %d == %d\n", edge[i][0],
edge[i][1], wt);
}
}

printf("Minimum Cost Spanning Tree: %d\n", minCost);
}

int main()
{
int edge[5][3] = { { 0, 1, 10 },
{ 0, 2, 6 },
{ 0, 3, 5 },
{ 1, 3, 15 },
{ 2, 3, 4 } };

kruskalAlgo(5, edge);

return 0;
}

travelling selasman
#include <stdio.h>
int tsp_g[10][10] = {
   {12, 30, 33, 10, 45},
   {56, 22, 9, 15, 18},
   {29, 13, 8, 5, 12},
   {33, 28, 16, 10, 3},
   {1, 4, 30, 24, 20}
};
int visited[10], n, cost = 0;

void travellingsalesman(int c){
   int k, adj_vertex = 999;
   int min = 999;
   
   visited[c] = 1;
   
   printf("%d ", c + 1);
   
   for(k = 0; k < n; k++) {
      if((tsp_g[c][k] != 0) && (visited[k] == 0)) {
         if(tsp_g[c][k] < min) {
            min = tsp_g[c][k];
         }
         adj_vertex = k;
      }
   }
   if(min != 999) {
      cost = cost + min;
   }
   if(adj_vertex == 999) {
      adj_vertex = 0;
      printf("%d", adj_vertex + 1);
      cost = cost + tsp_g[c][adj_vertex];
      return;
   }
   travellingsalesman(adj_vertex);
}

int main(){
   int i, j;
   n = 5;
   for(i = 0; i < n; i++) {
      visited[i] = 0;
   }
   printf("Shortest Path: ");
   travellingsalesman(0);
   printf("\nMinimum Cost: ");
   printf("%d\n", cost);
   return 0;
}

matrix multipication
#include<stdio.h>
#include<pthread.h>
#include<unistd.h>
#include<stdlib.h>
#define MAX 4
 
 
void *mult(void* arg)
{
    int *data = (int *)arg;
    int k = 0, i = 0;
     
    int x = data[0];
    for (i = 1; i <= x; i++)
           k += data[i]*data[i+x];
     
    int *p = (int*)malloc(sizeof(int));
         *p = k;
     
    pthread_exit(p);
}
 
int main()
{
 
    int matA[MAX][MAX];
    int matB[MAX][MAX];
     
     
    int r1=MAX,c1=MAX,r2=MAX,c2=MAX,i,j,k;
 
 
    for (i = 0; i < r1; i++)
            for (j = 0; j < c1; j++)
                   matA[i][j] = rand() % 10;
           
    for (i = 0; i < r1; i++)
            for (j = 0; j < c1; j++)
                   matB[i][j] = rand() % 10;
   
    for (i = 0; i < r1; i++){
        for(j = 0; j < c1; j++)
            printf("%d ",matA[i][j]);
        printf("\n");
    }
             
    for (i = 0; i < r2; i++){
        for(j = 0; j < c2; j++)
            printf("%d ",matB[i][j]);
        printf("\n");    
    }
     
     
    int max = r1*c2;
     
     
    pthread_t *threads;
    threads = (pthread_t*)malloc(max*sizeof(pthread_t));
     
    int count = 0;
    int* data = NULL;
    for (i = 0; i < r1; i++)
        for (j = 0; j < c2; j++)
               {
               
            data = (int *)malloc((20)*sizeof(int));
            data[0] = c1;
     
            for (k = 0; k < c1; k++)
                data[k+1] = matA[i][k];
     
            for (k = 0; k < r2; k++)
                data[k+c1+1] = matB[k][j];
             
                pthread_create(&threads[count++], NULL,
                               mult, (void*)(data));
                 
                    }
     
    printf("RESULTANT MATRIX IS :- \n");
    for (i = 0; i < max; i++)
    {
      void *k;
       
      pthread_join(threads[i], &k);
           
           
          int *p = (int *)k;
      printf("%d ",*p);
      if ((i + 1) % c2 == 0)
          printf("\n");
    }
 
     
 
  return 0;
}

#robin karp
#include <stdio.h>
#include <string.h>

#define d 256

void search(char pat[], char txt[], int q)
{
int M = strlen(pat);
int N = strlen(txt);
int i, j;
int p = 0;
int t = 0;
int h = 1;

for (i = 0; i < M - 1; i++)
h = (h * d) % q;

for (i = 0; i < M; i++) {
p = (d * p + pat[i]) % q;
t = (d * t + txt[i]) % q;
}

for (i = 0; i <= N - M; i++) {

if (p == t) {
for (j = 0; j < M; j++) {
if (txt[i + j] != pat[j])
break;
}

if (j == M)
printf("Pattern found at index %d \n", i);
}

if (i < N - M) {
t = (d * (t - txt[i] * h) + txt[i + M]) % q;

if (t < 0)
t = (t + q);
}
}
}

int main()
{
char txt[] = "HELLO WORLD";
char pat[] = "HELLO";

int q = 101;

search(pat, txt, q);
return 0;
}

#naive algorithm 
#include <stdio.h>
#include <string.h>

void search(char* pat, char* txt) {
    int M = strlen(pat);
    int N = strlen(txt);

    for (int i = 0; i <= N - M; i++) {
        int j;

        for (j = 0; j < M; j++) {
            if (txt[i + j] != pat[j]) {
                break;
            }
        }

        if (j == M) {
            printf("Pattern found at index %d\n", i);
        }
    }
}

int main() {
    char txt1[] = "AABAACAADAABAABA";
    char pat1[] = "AABA";
    printf("Example 1:\n");
    search(pat1, txt1);
   
    char txt2[] = "agd";
    char pat2[] = "g";
    printf("\nExample 2:\n");
    search(pat2, txt2);

    return 0;
}
