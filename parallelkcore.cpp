#include <iostream>
#include <fstream>
#include <algorithm>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>
#include <unordered_map>
#include <vector>
#define DEBUG 1

// A struct to represent an edge in the edge list
struct edge {
    unsigned int src;
    unsigned int tgt;
};

struct graph {
    unsigned int NODENUM;
    unsigned int EDGENUM;
    unsigned int *start_indices;
    unsigned int *end_indices;
    edge *edgeList;
	std::vector<std::vector<unsigned int>> adjList;
}g;

// DEFINITIONS
void buildAdjacencyList();
	// Given a set of labeled edges, find all virtual connected-components in an adjacency list
void extractComponents(unsigned int* edgeLabels, std::unordered_map<unsigned int,std::vector<unsigned int>> &components);
	// Recursively traverses an adjacency list
void dfsComponent(unsigned int src,std::vector<unsigned int> &nodes,bool *visited, unsigned int *edgeLabels);
	// Determines if a given edge has been already labeled
bool isEdgeLabeled(unsigned int src, unsigned int tgt, unsigned int *edgeLabels);
	// Prints a map of vectors which represented connected components
void printComponents(std::unordered_map<unsigned int,std::vector<unsigned int>> &components);

// Utility function to print a given array
template <class T>
void printArray(T *arr, unsigned int n) {
    for(unsigned int i = 0; i < n; i++) {
        std::cout<<arr[i]<<" ";
    }
    std::cout<<"\n";
}

void printEdgeLabels(unsigned int *edgeLabels) {
    std::cout<<"[ ";
    for(int i = 0; i < g.EDGENUM; i++) {
	std::cout<<"("<<(g.edgeList+i)->src<<","<<(g.edgeList+i)->tgt<<"): ";
        if(i + 1 == g.EDGENUM)
	    std::cout<<(int)edgeLabels[i]<<" ";
	else
	    std::cout<<(int)edgeLabels[i]<<", ";
    }
    std::cout<<" ]\n";
}

void printComponents(std::unordered_map<unsigned int,std::vector<unsigned int>> &components) {
	std::cout<<"NUMBER OF COMPONENTS: "<<components.size()<<"\n";
	for(std::unordered_map<unsigned int,std::vector<unsigned int>>::iterator it = components.begin();
					it != components.end(); it++) {
		std::cout<<"["<<it->first<<"] -> "; 
		std::vector<unsigned int> v = (std::vector<unsigned int>) it->second;
		for(unsigned int n = 0; n < v.size(); n++) {
			std::cout<<v[n]<<" ";
		}
		std::cout<<";\n";
	}
}

long long currentTimeMilliS = 0;

long long currentTimeStamp() {
    struct timeval te;
    gettimeofday(&te, NULL); // get current time
    long long milliseconds = te.tv_sec*1000LL + te.tv_usec/1000; // calculate milliseconds
    return milliseconds;
}

void reset() {
    currentTimeMilliS = currentTimeStamp();
}

long long getTimeElapsed() {
    long long newTime = currentTimeStamp();
    long long timeElapsed = newTime - currentTimeMilliS;
    currentTimeMilliS = newTime;
    return timeElapsed;
}

// Memory maps input file
void createMemoryMap(char *fileName) {
    unsigned int binFile = open(fileName, O_RDWR);
    long fileSizeInByte;
    struct stat sizeResults;
    assert(stat(fileName, &sizeResults) == 0);
    fileSizeInByte = sizeResults.st_size;
    // map edge list into memory    
    g.edgeList = (edge *)mmap(NULL, fileSizeInByte, PROT_READ | PROT_WRITE, MAP_SHARED, binFile, 0);
    close(binFile);
}

// In-memory edge list instead of memory mapped file
void createInMemoryEdgeList(char *fileName) {
        g.edgeList = new edge[g.EDGENUM];
        std::ifstream is;
        is.open(fileName, std::ios::in | std::ios::binary);
        unsigned int src, tgt;
        unsigned int updatedEdgeNum = g.EDGENUM;
        for(unsigned int i = 0; i < g.EDGENUM; i++) {
                is.read((char *)(&src), sizeof(unsigned int));
                is.read((char *)(&tgt), sizeof(unsigned int));
                assert(src >= 0 && src <= g.NODENUM);
                assert(tgt >= 0 && tgt <= g.NODENUM);
                (g.edgeList + i)->src = src;
                (g.edgeList + i)->tgt = tgt;
        }
        is.close();
}

// Reads and input file, removes self-loops, creates a new file with an undirected graph (src->tgt | tgt->src)
void doubleAndReverseGraph(char *inputFile, char *outputFile) {
    std::ifstream is;
    is.open(inputFile, std::ios::in | std::ios::binary);
    std::ofstream os;
    os.open(outputFile, std::ios::out | std::ios::binary | std::ios::app);
    unsigned int src, tgt;
    unsigned int updatedEdgeNum = g.EDGENUM;
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        // read the source -> target value and arrange bytes
	is.read((char *)(&src), sizeof(unsigned int));
        is.read((char *)(&tgt), sizeof(unsigned int));
        src = htonl(src);
        tgt = htonl(tgt);
	// check values are positive
        assert(src >= 0 && src <= g.NODENUM);
        assert(tgt >= 0 && tgt <= g.NODENUM);
        // Removes self loops
        if(src != tgt) {
            os.write((char *)&src, sizeof(unsigned int));
            os.write((char *)&tgt, sizeof(unsigned int));
        }
        else {
	    // discard edge if it is a self loop
            updatedEdgeNum--;
        }
    }
    // reset pointer position in file
    is.seekg(0, std::ios::beg);
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        is.read((char *)(&src), sizeof(unsigned int));
        is.read((char *)(&tgt), sizeof(unsigned int));
        src = htonl(src);
        tgt = htonl(tgt);
        assert(src >= 0 && src <= g.NODENUM);
        assert(tgt >= 0 && tgt <= g.NODENUM);
        // Removes self loops
        if(src != tgt) {
            os.write((char *)(&tgt), sizeof(unsigned int));
            os.write((char *)(&src), sizeof(unsigned int));
        }
    }
    g.EDGENUM = updatedEdgeNum;
    g.EDGENUM *= 2;
	is.close();
    os.close();
}

bool compareEdges(const edge& a, const edge& b) {
    if(a.src < b.src)
        return true;
    if(a.src == b.src) {
        if(a.tgt < b.tgt)
            return true;
    }
    return false;
}

// Compares indices according to their corresponding edges
int compareByEdges(const void * a, const void * b) {
    if ((g.edgeList + *(unsigned int *)a)->src < (g.edgeList + *(unsigned int *)b)->src)
        return -1;
    if ((g.edgeList + *(unsigned int *)a)->src == (g.edgeList + *(unsigned int *)b)->src){
        if ((g.edgeList + *(unsigned int *)a)->tgt < (g.edgeList + *(unsigned int *)b)->tgt)
            return -1;
        if ((g.edgeList + *(unsigned int *)a)->tgt == (g.edgeList + *(unsigned int *)b)->tgt)
            return 0;
        if ((g.edgeList + *(unsigned int *)a)->tgt > (g.edgeList + *(unsigned int *)b)->tgt)
            return 1;
    }
    if ((g.edgeList + *(unsigned int *)a)->src > (g.edgeList + *(unsigned int *)b)->src)
        return 1;
}

// Formats the graph by sorting it and tracing original indices in the graph
void formatGraph(unsigned int *originalIndices) {
    unsigned int *indices = new unsigned int[g.EDGENUM];
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        indices[i] = i;
    }
    qsort(indices, g.EDGENUM, sizeof(unsigned int), compareByEdges);
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        originalIndices[indices[i]] = i;
    }
    std::sort(g.edgeList, g.edgeList + g.EDGENUM, compareEdges);
    //qsort(g.edgeList, g.EDGENUM, sizeof(edge), compareByEdges);
    delete [] indices;
}

// Finds the start and end indices  of each node in the graph
void findStartAndEndIndices() {
    g.start_indices = new unsigned int[g.NODENUM + 1];
    g.end_indices = new unsigned int[g.NODENUM + 1];
    std::fill_n(g.start_indices, g.NODENUM + 1, 0);
    std::fill_n(g.end_indices, g.NODENUM + 1, 0);
    unsigned int i;
    unsigned int old = g.edgeList->src;
    g.start_indices[old] = 0;
    for(i = 0; i < g.EDGENUM; i++) {
            if((g.edgeList + i)->src != old) {
                    g.end_indices[old] = i - 1;
                    old = (g.edgeList + i)->src;
                    g.start_indices[old] = i;
            }
    }
    g.end_indices[old] = i - 1;
}

bool isGraphEmpty(unsigned int *edgeLabels) {
        for(unsigned int i = 0; i < g.EDGENUM; i++) {
                if(edgeLabels[i] == -1)
                        return false;
        }
        return true;
}

// Computes the degree of each node in the graph
void findDegree(unsigned int *edgeLabels, float *degree) {
    std::fill_n(degree, g.NODENUM + 1, 0);
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        // If edge hasn't been deleted yet. An edge is considered deleted
        // when it has been labeled.
        if(edgeLabels[i] == -1) {
                degree[(g.edgeList + i)->src]++;
                degree[(g.edgeList + i)->tgt]++;
        }
    }
    for(unsigned int i = 0; i < g.NODENUM + 1; i++) {
        degree[i] /= 2;
    }
}

// Given a list of labeled edgess, it identifies disjoint components
void extractComponents(unsigned int *edgeLabels, std::unordered_map<unsigned int,std::vector<unsigned int>> &components) {
	// std::unordered_map<unsigned int,std::vector<unsigned int>> components;
	bool *visited = new bool[g.adjList.size()];
	std::fill_n(visited,g.adjList.size(),false);
	for(unsigned int i = 0; i < g.EDGENUM; ++i) {
		// only pick root nodes from labels which hasn't been labeled
		if(edgeLabels[i] == -1) {
			edge* e = (g.edgeList+i);
			int currSrc = e->src;
			if(!visited[currSrc]) {
				components[currSrc] = std::vector<unsigned int>(0);
				dfsComponent(currSrc,components[currSrc],visited,edgeLabels);
			}
		}
	}
	if(DEBUG)
		printComponents(components);
	delete [] visited;
}

void dfsComponent(unsigned int src, std::vector<unsigned int> &nodes,bool *visited, unsigned int *edgeLabels) {
	visited[src] = true;
	for(unsigned int i =0; i < g.adjList[src].size(); ++i) {
		unsigned int currChild = g.adjList[src][i];
		if(!(visited[currChild] || isEdgeLabeled(src,currChild,edgeLabels))) {
			nodes.push_back(currChild);
			dfsComponent(currChild,nodes,visited,edgeLabels);
		}
	}
}

bool isEdgeLabeled(unsigned int src, unsigned int tgt, unsigned int *edgeLabels) {
	for(unsigned int i = 0; i < g.EDGENUM; ++i){
		edge *e = (g.edgeList+i);
		if(e->src > src) break;
		if(e->src == src && e->tgt == tgt) {
			return (edgeLabels[i] != -1);
		}
	}
	return false;
}

// Builds an adjacency list representation, only once
void buildAdjacencyList() {	
	for(unsigned int n = 0; n < g.NODENUM+1; n++) {
		g.adjList.push_back(std::vector<unsigned int>(0));
	}
	for(unsigned int i = 0; i < g.EDGENUM; ++i) {
		unsigned int src = (g.edgeList+i)->src;
		unsigned int tgt = (g.edgeList+i)->tgt;
		g.adjList[src].push_back(tgt);
	}
	if(DEBUG) {
		std::cout<<"ADJANCECY LIST: \n";
		for(unsigned int n = 0; n < g.NODENUM + 1;n++) {
			if(g.adjList[n].size() > 0) {
				std::cout<<"HEAD: "<<n<<": [ ";
				for(unsigned int a = 0; a < g.adjList[n].size(); ++a) {
						std::cout<<g.adjList[n][a]<<" ";
				}
				std::cout<<" ]\n";
			}
		}
	}
}

void scan(unsigned int *deg, unsigned int level, unsigned int* curr, long *currTailPtr) {

    // Size of cache line
    const long BUFFER_SIZE_BYTES = 2048;
    const long BUFFER_SIZE = BUFFER_SIZE_BYTES/sizeof(unsigned int);

    unsigned int buff[BUFFER_SIZE];
    long index = 0;

#pragma omp for schedule(static)
    for (unsigned int i = 0; i < g.NODENUM + 1; i++) {

        if (deg[i] == level) {

            buff[index] = i;
            index++;

            if (index >= BUFFER_SIZE) {
                long tempIdx = __sync_fetch_and_add(currTailPtr, BUFFER_SIZE);

                for (long j = 0; j < BUFFER_SIZE; j++) {
                    curr[tempIdx+j] = buff[j];
                }
                index = 0;
            }
        }

    }

    if (index > 0) {
        long tempIdx = __sync_fetch_and_add(currTailPtr, index);

        for (long j = 0; j < index; j++)
            curr[tempIdx+j] = buff[j];
    }

#pragma omp barrier

}


void processSubLevel(unsigned int *curr, long currTail,
        unsigned int *deg, unsigned int *edgeLabels, unsigned int level, unsigned int *next, long *nextTailPtr) {

    // Size of cache line
    const long BUFFER_SIZE_BYTES = 2048;
    const long BUFFER_SIZE = BUFFER_SIZE_BYTES/sizeof(unsigned int);

    unsigned int buff[BUFFER_SIZE];
    long index = 0;

#pragma omp for schedule(static)
    for (long i = 0; i < currTail; i++) {
        unsigned int v = curr[i];
        // Do nothing if node doesn't exist in the graph
        if(g.start_indices[v] == 0 && g.end_indices[v] == 0) {
            ;
        }
        else {
            //For all neighbors of vertex v
            for (unsigned int j = g.start_indices[v]; j <= g.end_indices[v]; j++) {
                if(edgeLabels[j] == -1) {
                    unsigned int u = (g.edgeList + j)->tgt;
                    if (deg[u] > level) {
                        unsigned int du =  __sync_fetch_and_sub(&deg[u], 1);
                        if (du == (level+1)) {
                            buff[index] = u;
                            index++;
                            if (index >= BUFFER_SIZE) {
                                long tempIdx = __sync_fetch_and_add(nextTailPtr, BUFFER_SIZE);
                                for(long bufIdx = 0; bufIdx < BUFFER_SIZE; bufIdx++)
                                    next [tempIdx + bufIdx] = buff[bufIdx];
                                index = 0;
                            }
                        }
                    }
                }
            }
        }
    }

    if (index > 0) {
        long tempIdx =  __sync_fetch_and_add(nextTailPtr, index);;
        for (long bufIdx = 0; bufIdx < index; bufIdx++)
            next [tempIdx + bufIdx] = buff[bufIdx];
    }

#pragma omp barrier

#pragma omp for schedule(static)
    for (long i = 0; i < *nextTailPtr; i++) {
        unsigned int u = next[i];
        if (deg[u] != level)
            deg[u] = level;
    }


#pragma omp barrier

}

//ParK to compute k core decomposition in parallel
void parKCore(unsigned int *deg, unsigned int *edgeLabels) {
    unsigned int *curr = new unsigned int[g.NODENUM];
    unsigned int *next = new unsigned int[g.NODENUM];

    long currTail = 0;
    long nextTail = 0;

#pragma omp parallel
{
    unsigned int tid = omp_get_thread_num();
    long todo = g.NODENUM;
    unsigned int level = 0;

    while (todo > 0) {

        scan(deg, level, curr, &currTail);

        while (currTail > 0) {
            todo = todo - currTail;

            processSubLevel(curr, currTail, deg, edgeLabels, level, next, &nextTail);

            if (tid == 0) {
                unsigned int *tempCurr = curr;
                curr = next;
                next = tempCurr;

                currTail = nextTail;
                nextTail = 0;
            }

#pragma omp barrier

        }

        level = level + 1;
#pragma omp barrier

    }
}
    delete [] curr;
    delete [] next;

}

void labelEdgesAndUpdateDegree(unsigned int peel, bool *isFinalNode, float *degree, unsigned int *edgeLabels) {
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        unsigned int src = (g.edgeList + i)->src;
        unsigned int tgt = (g.edgeList + i)->tgt;
        if(isFinalNode[src] && isFinalNode[tgt] && edgeLabels[i] == -1) {
            edgeLabels[i] = peel;
            degree[src] -= 0.5;
            degree[tgt] -= 0.5;
        }
    }
}

void labelAndDeletePeelOneEdges(float *degree, unsigned int *edgeLabels) {
    float *tmp = new float[g.NODENUM + 1];
    std::copy(degree, degree + g.NODENUM + 1, tmp);
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        unsigned int src = (g.edgeList + i)->src;
        unsigned int tgt = (g.edgeList + i)->tgt;
        if(edgeLabels[i] == -1) {
                if(tmp[src] == 1 || tmp[tgt] == 1) {
                        edgeLabels[i] = 1;
                        degree[src] -= 0.5;
                        degree[tgt] -= 0.5;
                }
        }
    }
    delete [] tmp;
}

void writeToFile(unsigned int *edgeIndices, unsigned int *edgeLabels) {
        std::ofstream outputFile;
        outputFile.open("graph-decomposition.csv");
        for(unsigned int i = 0; i < g.EDGENUM; i++) {
                outputFile<<(g.edgeList + edgeIndices[i])->src<<","<<(g.edgeList + edgeIndices[i])->tgt<<","<<edgeLabels[i]<<"\n";
        }
        outputFile.close();
}

void writeMetaData(unsigned int NODENUM, unsigned int EDGENUM, long long preprocessingTime, long long algorithmTime) {
        std::ofstream outputFile;
        outputFile.open("graph-decomposition-info.json");
        outputFile<<"{\n";
        outputFile<<"\"vertices\":"<<NODENUM<<",\n";
        outputFile<<"\"edges\":"<<EDGENUM<<",\n";
        outputFile<<"\"preprocessing-time\":"<<preprocessingTime<<",\n";
        outputFile<<"\"algorithm-time\":"<<algorithmTime<<"\n}";
        outputFile.close();
}

int main(int argc, char *argv[]) {
    char *tmpFile = (char*) "tmp.bin";
    remove(tmpFile);
    g.EDGENUM = atoi(argv[2]);
    g.NODENUM = atoi(argv[3]);
	std::unordered_map<unsigned int,std::vector<unsigned int>> components;
	reset();
    doubleAndReverseGraph(argv[1], tmpFile);
    if(DEBUG)
        std::cout<<"DOUBLED AND REVERSED GRAPH\n";
    unsigned int *originalIndices = new unsigned int[g.EDGENUM];
    unsigned int *edgeLabels = new unsigned int[g.EDGENUM];
    std::fill_n(edgeLabels, g.EDGENUM, -1);
    createMemoryMap(tmpFile);
	if(DEBUG)
        std::cout<<"CREATED MEMORY MAP\n";
    formatGraph(originalIndices);
    buildAdjacencyList();
	extractComponents(edgeLabels,components); // TODO - remove, test only
	long long preprocessingTime = getTimeElapsed();
    reset();
    if(DEBUG)
        std::cout<<"FORMATTED GRAPH\n";
    findStartAndEndIndices();
    if(DEBUG)
        std::cout<<"START AND END INDICES COMPUTED\n";
    float *degree = new float[g.NODENUM + 1];
    findDegree(edgeLabels, degree);
    unsigned int *core = new unsigned int[g.NODENUM + 1];
    while(!isGraphEmpty(edgeLabels)) {
        std::copy(degree, degree + g.NODENUM + 1, core);
		if(DEBUG) { 
			std::cout<<"CURRENT NODES DEGREES: ";
			printArray(degree,g.NODENUM+1);
		}    
        parKCore(core, edgeLabels);
        unsigned int mc = *std::max_element(core, core + g.NODENUM + 1);
        if(DEBUG)
 			std::cout<<"CURRENT MAXIMUM CORE : "<<mc<<"\n";
        bool *isFinalNode = new bool[g.NODENUM + 1];
        std::fill_n(isFinalNode, g.NODENUM + 1, false);
        for(unsigned int i = 0; i <= g.NODENUM; i++) {
            if(core[i] == mc) {
                isFinalNode[i] = true;
            }
        }
        labelEdgesAndUpdateDegree(mc, isFinalNode, degree, edgeLabels);
		if(DEBUG) {
			std::cout<<"CURRENT EDGES CORE VALUES: ";
			printArray(core, g.NODENUM+1);
			std::cout<<"CURRENT EDGE LABELS: ";
			printEdgeLabels(edgeLabels);
			components.clear();
			extractComponents(edgeLabels,components);
		}
        delete [] isFinalNode;
    }
    g.EDGENUM /= 2;
    unsigned int *originalLabels = new unsigned int[g.EDGENUM];
    if(DEBUG)
        std::cout<<"RECONSTRUCTING ORIGINAL LABELS\n";
    for(unsigned int i = 0; i < g.EDGENUM; i++) {
        originalLabels[i] = edgeLabels[originalIndices[i]];
    }
    long long algorithmTime = getTimeElapsed();
    writeToFile(originalIndices, originalLabels);
    writeMetaData(atoi(argv[3]), atoi(argv[2]), preprocessingTime, algorithmTime);
    remove(tmpFile);
    delete [] core;
    delete [] degree;
    delete [] g.start_indices;
    delete [] g.end_indices;
    delete [] edgeLabels;
    delete [] originalLabels;
    delete [] originalIndices;
    return 0;
}
