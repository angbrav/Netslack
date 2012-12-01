#include <stdio.h>  
#include <stdlib.h>
#include <math.h> 
#include <omp.h>

#define NDELAY	600
#define _GNU_SOURCE

#include <string.h>  

struct edge{
    	int source;
    	int sink;
    	double de;
};

struct node{
    
	double dn;
	double dno;
	double dnc;
	double m;
	double at;
	double rat;
	double slack;
	int position;
    	int f,t,t_tmp,f_tmp;
	char name[10];
	int from[20];
	int to[20];
};

double delay(double iat) { 
    int i; 
    double sc = 0.0; 
    
    for(i = 0; i < NDELAY; i++) { 
        sc += 1.0 / (pow(tan(i * M_PI/NDELAY), 2.0) + 1.0) + pow(sin(i * M_PI/NDELAY), 2.0); 
    } 
    sc /= NDELAY; 
    
    return(iat*sc); 
} 

double CalculateAT(int to[20], int nodes_to, struct edge edges[], struct node nodes[]){

        double aux=0,r=0;
        int i;

        for (i=0;i<nodes_to;i++){
                aux=nodes[edges[to[i]].source].at + nodes[edges[to[i]].source].dn + edges[to[i]].de;
                //printf("Parcial at , from node: %d, to node: %d, at and dn from node: %f - %f, de edge: %f\n",edges[to[i]].source, edges[to[i]].sink, nodes[edges[to[i]].source].at, nodes[edges[to[i]].source].dn,edges[to[i]].de);
                if (aux>r)
                        r=aux;
        }
        return delay(r);

}

double CalculateRAT(int from[20], int nodes_from, struct edge edges[], struct node nodes[]){

        double auxRAT=0,r=2000000000;
        int i;

        for (i=0;i<nodes_from;i++){
                auxRAT= nodes[edges[from[i]].sink].rat - nodes[edges[from[i]].source].dn - edges[from[i]].de;
                //printf("Parcial at , from node: %d, to node: %d, rat and dn from node: %f - %f, de edge: %f\n",edges[from[i]].source, edges[from[i]].sink, nodes[edges[from[i]].sink].rat, nodes[edges[from[i]].source].dn,edges[from[i]].de);
                if (auxRAT<r)
                        r=auxRAT;
        }
        return delay(r);
}

// MAIN FUNCTION

int main(int argc, char *argv[]){

	//DECLARATION OF VARIABLES

		

	char buffer[128], cha;
	char *token;
	FILE * myfile;
	int numNodes,numInputs,numOutputs,i,from,to,f, aux, cInputs=0, cOutputs=0, cNodes=0;
	struct node *N;
	int *inputs, *outputs;
	struct edge *E;
	int edgesIndex=0, numEdges=0;
	struct edge *e;
	double dei,deo,mAux,start, end;

	start = omp_get_wtime( );

	myfile = fopen(argv[1],"r");            //open file to read the network
        fgets(buffer,128,myfile);
    	while (!feof(myfile)){
		switch(buffer[0]){
			case 'e':
				numEdges++;
				break;
			default:
				break;
		}
        	fgets(buffer,128,myfile);
	}
    	fclose(myfile);

	myfile = fopen(argv[1],"r");            //open file to read the network
    	fgets(buffer,128,myfile);           //bufer por línea
	numNodes = atoi(strtok(buffer, " "));   //leer nro total de nodos
	numInputs= atoi(strtok(NULL, " "));     //leer nro de entradas
	numOutputs= atoi(strtok(NULL, " "));    //leer nro de salidas
    
	//CREATING ARRAY FOR NODES, INPUTS AND OUTPUTS

	N = (struct node *) malloc(numNodes * sizeof(struct node));  
	E = (struct edge *) malloc(numEdges * sizeof(struct edge));  
	inputs = (int *) malloc(numInputs * sizeof(int));
	outputs = (int *) malloc(numOutputs * sizeof(int));
    	int *L = (int *)malloc(numNodes*sizeof(int));
    	int *LI = (int *)malloc(numNodes*sizeof(int));
    	int *S = (int *)malloc(numNodes*sizeof(int));

	//OBTAINING  INPUT NODES AND SETTING AT FOR THEM FROM FILE
	aux = 0;
        fgets(buffer,128,myfile);
    	while (!feof(myfile)){
		cha = *strtok(buffer," ");
		switch(cha){
			case 'i':
				token = strtok(NULL, " ");
				inputs[cInputs]=atoi(token+1);
				cInputs++;
				break;
			case 'a':
       				token = strtok(NULL, " ");
				N[atoi(token+1)].at=atof(strtok(NULL, " "));
				break;
			case 'o':
				token = strtok(NULL, " ");
				outputs[cOutputs]=atoi(token+1);
				cOutputs++;
				break;
			case 'r':
        			token=strtok(NULL, " ");
				N[atoi(token+1)].rat=atof(strtok(NULL, " "));
				break;
			case 'n':
				token = strtok(NULL, " ");
				aux = atoi(token+1);
				strcpy(N[aux].name,token);
				N[aux].position = aux;
				N[aux].dno = atof(strtok(NULL, " "));
				N[aux].dnc = atof(strtok(NULL, " "));
				N[aux].m = atof(strtok(NULL, " "));
				N[aux].dn=0;
        			/*if (aux>=numInputs)
        				N[aux].at=0;
        			if (aux < (numNodes-numOutputs))
        				N[aux].rat=2000000000;*/
				N[aux].slack=0;
        			N[aux].f=0;
        			N[aux].t=0;
				break;
			case 'e':
				token = strtok(NULL, " ");
            			from = atoi(token+1);        //read in line to get origin node
            			token = strtok(NULL, " ");
				to = atoi(token+1);          //read in line to get destination node
				mAux = N[to].m;                                     //obtain m from destination node
            			dei=atof(strtok(NULL, " "));                        //obtain intrinsic edge delay
				deo=atof(strtok(NULL, " "));                        //obtain characteristic delay
				E[edgesIndex].de =  dei + mAux*deo;                           //calculate edge's delay: de= dei + m*deo
				E[edgesIndex].source = from;                                  //add destination node to the origin node's list of sink nodes
				E[edgesIndex].sink = to;

				N[from].from[N[from].f] = edgesIndex;
				N[from].f++;                                        //increase f for source node
				N[to].to[N[to].t] = edgesIndex;
            			N[to].t++;                                         //increase t for destination node

				edgesIndex++;
				
				break;
			default:
				break;
		}
        	fgets(buffer,128,myfile);
	}

    	fclose(myfile);
	end = omp_get_wtime( );
	printf("Time reading: %f\n",end - start);
	
	//printf("FILE READ SUCCESSFUL\n");
    	/*
	printf("\nEdges:\n");
	displayEdgesList(edgesList);
	for (i=0; i<numNodes; i++){
		printf("Node %d\n",N[i].name);
		printf("Source list: ");
		displayEdgesList(N[i].source);
		printf("Sink list: ");
		displayEdgesList(N[i].sink);
		printf("Temp sink list: ");
		displayEdgesList(N[i].sink_tmp);
	}*/

    	//INPUT FILE CLOSED, NO FURTHER READINGS
	//COMPUTE NODE DELAY 
    
	start = omp_get_wtime( );
	for (i=0; i<numNodes; i++){
		N[i].dn = N[i].dno + N[i].f*N[i].dnc;
		N[i].t_tmp = N[i].t;
		N[i].f_tmp = N[i].f;
	}
	
	/*for (i=0; i<numEdges; i++){
		printf("Node %d: %f\n",i,E[i].de);
	}*/	

	int countL=0, countS=0, insertS=0, nfrom, nto;
	
	for (i=0; i<numInputs; i++){
		S[insertS]=inputs[i];
		insertS++;
	}
	
	
	while (!(countS==numNodes)){
		nfrom = S[countS];
		countS++;
		L[countL]=nfrom;		
		countL++;
		for (i=0;i<N[nfrom].f;i++){
			nto=E[N[nfrom].from[i]].sink;
			N[nto].t_tmp--;
			if (N[nto].t_tmp == 0){
				S[insertS]=nto;
				insertS++;
			}
		}
	}
	insertS=0; countL=0; countS=0;

	for (i=0; i<numOutputs; i++){
        	S[insertS]=outputs[i];
                insertS++;
        }

        while (!(countS==numNodes)){
	        nto = S[countS];
                countS++;
                LI[countL]=nto;
                countL++;
                for (i=0;i<N[nto].t;i++){
               		nfrom=E[N[nto].to[i]].source;
                        N[nfrom].f_tmp--;
                        if (N[nfrom].f_tmp == 0){
                       		S[insertS]=nfrom;
                                insertS++;
                        }
                }
	}  

	free(S);
	/*printf("Sorted list: ");
	for (i=0; i<numNodes; i++){
		printf("%d -> ",LI[i]);
	}
	printf("\n");
	*/
	for (i=numInputs; i<numNodes; i++)
        	N[L[i]].at = CalculateAT(N[L[i]].to, N[L[i]].t, E, N);
    

	//COMPUTE RAT
    
   	for (i=numOutputs; i<numNodes; i++)
        	N[LI[i]].rat = CalculateRAT(N[LI[i]].from, N[LI[i]].f, E, N);
	
	//COMPUTE SLACK
    
    	for (i=0; i<numNodes; i++)
        	N[i].slack=N[i].rat-N[i].at;    

	end = omp_get_wtime( );
	printf("Time computing: %f\n",end - start);
	//      WRITE OUTPUT FILE
	char *fout=(char*) malloc(sizeof(char)*(strlen(argv[1])+2));

        for (i=0;i<strlen(argv[1])-3;i++){
            fout[i] = argv[1][i];
        }
        fout[strlen(argv[1])-3] = argv[1][strlen(argv[1])];
        strcat(fout,".out");
        myfile=fopen(fout,"w");

	
    	for (i=0; i<numNodes; i++)        
        	fprintf(myfile,"%s %.2f %.2f %.2f\n", N[i].name, N[i].rat, N[i].at, N[i].slack);    
    	fclose ( myfile );
    	//printf("END OF PROGRAM\n");   
	return(0);
}
