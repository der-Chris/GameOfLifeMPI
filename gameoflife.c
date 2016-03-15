#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

void show(unsigned* currentfield, int w, int h) {
    printf("\033[H");
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] ? "\033[07m  \033[m" : "  ");
        printf("\033[E");
    }
    fflush(stdout);
}


float convert2BigEndian( const float inFloat )
{
    float retVal;
    char *floatToConvert = ( char* ) & inFloat;
    char *returnFloat    = ( char* ) & retVal;

    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];

    return retVal;
}
//myCurrentfield, w, h, t, "output", newComm, rank, size, beginPos, endPos, myCoords, myW, myH
void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix, MPI_Comm newComm, int rank, int size, int *beginPos, int *endPos, int *myCoords, int myW, int myH) {
    char name[1024] = "\0";
    sprintf(name, "out/%s_%d_%d.vtk", prefix, rank, t);
    FILE* outfile = fopen(name, "w");

    beginPos[0] = beginPos[0] == 0 ? 0 : beginPos[0] - 1 * myCoords[0];
    beginPos[1] = beginPos[1] == 0 ? 0 : beginPos[1] - 1 * myCoords[1];


    /*Write vtk header */
    fprintf(outfile,"# vtk DataFile Version 3.0\n");
    fprintf(outfile,"frame %d\n", t);
    fprintf(outfile,"BINARY\n");
    fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
    fprintf(outfile,"DIMENSIONS %d %d %d \n", myW, myH, 1); //w h
    fprintf(outfile,"SPACING 1.0 1.0 1.0\n");//or ASPECT_RATIO
    fprintf(outfile,"ORIGIN %d %d 0\n", beginPos[0], beginPos[1]); //x y z
    fprintf(outfile,"POINT_DATA %d\n", myW * myH);
    fprintf(outfile,"SCALARS data float 1\n");
    fprintf(outfile,"LOOKUP_TABLE default\n");

    for (int y = 0; y < myH; y++) {
        for (int x = 0; x < myW; x++) {
            float value = currentfield[calcIndex(w, x,y)] != 0.0 ? 1.0:0.0;
            value = convert2BigEndian(value);
            fwrite(&value, 1, sizeof(float), outfile);
        }
    }
    fclose(outfile);
}

int getNumberOfAliveCells(unsigned* currentfield, int w, int h, int x, int y) {
    int aliveFields = 0;
    for (int j = y - 1; j <= y + 1; j = j + 1) {
        for (int i = x - 1; i <= x + 1; i = i + 1) {
            if (currentfield[calcIndex(w, i,j)]) {
                aliveFields += 1;
            }
        }
    }
    return aliveFields;
}

void calculateNeighbourProcesses(int *myCoords, MPI_Comm comm, int *upperProcess, int *lowerProcess, int *rightProcess, int *leftProcess, int *upperRightProcess, int *upperLeftProcess, int *lowerRightProcess, int *lowerLeftProcess)
{
    int *queryCoords = calloc(2, sizeof(unsigned));
    queryCoords[0] = myCoords[0] + 1;
    queryCoords[1] = myCoords[1];
    MPI_Cart_rank(comm, queryCoords, upperProcess);

    queryCoords[0] = myCoords[0] - 1;
    queryCoords[1] = myCoords[1];
    MPI_Cart_rank(comm, queryCoords, lowerProcess);

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] + 1;
    MPI_Cart_rank(comm, queryCoords, rightProcess);

    queryCoords[0] = myCoords[0];
    queryCoords[1] = myCoords[1] - 1;
    MPI_Cart_rank(comm, queryCoords, leftProcess);

    queryCoords[0] = myCoords[0] + 1;
    queryCoords[1] = myCoords[1] + 1;
    MPI_Cart_rank(comm, queryCoords, upperRightProcess);

    queryCoords[0] = myCoords[0] + 1;
    queryCoords[1] = myCoords[1] - 1;
    MPI_Cart_rank(comm, queryCoords, upperLeftProcess);

    queryCoords[0] = myCoords[0] - 1;
    queryCoords[1] = myCoords[1] + 1;
    MPI_Cart_rank(comm, queryCoords, lowerRightProcess);

    queryCoords[0] = myCoords[0] - 1;
    queryCoords[1] = myCoords[1] - 1;
    MPI_Cart_rank(comm, queryCoords, lowerLeftProcess);

    free(queryCoords);
}

void calculateBorders(unsigned *currentfield, int myW, int myH, int *upperBorder, int *lowerBorder, int *rightBorder, int *leftBorder, int *upperRightBorder, int *upperLeftBorder, int *lowerRightBorder, int *lowerLeftBorder) {

}


int evolve(unsigned* currentfield, unsigned* newfield, int w, int h, MPI_Status status, MPI_Comm comm, int rank, int size, int *beginPos, int *endPos, int *myCoords, int myW, int myH) {
    int changes = 0;

    int *lowerBorder = calloc(myW, sizeof(unsigned));
    int *upperBorder = calloc(myW, sizeof(unsigned));
    int *rightBoder = calloc(myH, sizeof(unsigned));
    int *leftBorder = calloc(myH, sizeof(unsigned));
    int *upperRightBorder = calloc(1, sizeof(unsigned));
    int *upperLeftBorder = calloc(1, sizeof(unsigned));
    int *lowerRightBorder = calloc(1, sizeof(unsigned));
    int *lowerLeftBorder = calloc(1, sizeof(unsigned));

    int *lowerBorder = calloc(myW, sizeof(unsigned));
    int *upperBorder = calloc(myW, sizeof(unsigned));
    int *rightBoder = calloc(myH, sizeof(unsigned));
    int *leftBorder = calloc(myH, sizeof(unsigned));
    int *upperRightBorder = calloc(1, sizeof(unsigned));
    int *upperLeftBorder = calloc(1, sizeof(unsigned));
    int *lowerRightBorder = calloc(1, sizeof(unsigned));
    int *lowerLeftBorder = calloc(1, sizeof(unsigned));

    int *upperProcess= calloc (1, sizeof(int));
    int *lowerProcess= calloc (1, sizeof(int));
    int *rightProcess= calloc (1, sizeof(int));
    int *leftProcess= calloc (1, sizeof(int));
    int *upperRightProcess= calloc (1, sizeof(int));
    int *upperLeftProcess= calloc (1, sizeof(int));
    int *lowerRightProcess= calloc (1, sizeof(int));
    int *lowerLeftProcess= calloc (1, sizeof(int));

    calculateNeighbourProcesses(myCoords, comm, upperProcess, lowerProcess, rightProcess, leftProcess, upperRightProcess, upperLeftProcess, lowerRightProcess, lowerLeftProcess);

    calculateBorders(currentfield, myW, myH, upperBorder, lowerBorder, rightBorder, leftBorder, upperRightBorder, upperLeftBorder, lowerRightBorder, lowerLeftBorder);



    //printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d] & [%d] and: [%d]\n",rank, size, sourceProcessH[0], destProcessH[0], sourceProcessV[0], destProcessV[0]);

    /*
     for (int y = 0; y < h; y++) {
     for (int x = 0; x < w; x++) {
     */
    //TODO add game of life rules

    //Send Upper Border -> receive lower Border
    int sendcount = myW;
    int *sendbuf = calloc(sendcount, sizeof(unsigned));;
    MPI_Datatype sendtype = MPI_UNSIGNED;
    int dest = (rank + 1) % size;
    int sendtag = 1;
    int recvcount = myW;
    MPI_Datatype recvtype = MPI_UNSIGNED;
    int source = (rank + size - 1) % size;
    int recvtag = 1;
    for (int i = 0; i < myW; i = i + 1)
    {
        sendbuf[i] = currentfield[calcIndex(myW, i, myH - 1)];
    }
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, lowerBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Lower Border -> receive Upper Border
    dest = (rank + size - 1) % size;
    source = (rank + 1) % size;
    for (int i = 0; i < myW; i = i + 1)
    {
        sendbuf[i] = currentfield[calcIndex(myW, i, 0)];
    }
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, upperBorder, recvcount, recvtype, source, recvtag, comm, &status);
    printf("DEBUG: Thread No: [%d] sending to [%d] and send [%d][%d][%d][%d]\n", rank, dest, sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3]);
    printf("DEBUG: Thread No: [%d] received from [%d] and send [%d][%d][%d][%d]\n", rank, source, upperBorder[0], upperBorder[1], upperBorder[2], upperBorder[3]);

    //Send Right Border -> receive Left Border
    sendcount = myH;
    dest = (rank + 1) % size;
    recvcount = myW;
    source = (rank + size - 1) % size;
    for (int i = 0; i < myH; i = i + 1)
    {
        sendbuf[i] = currentfield[calcIndex(myW, myW - 1, i)];
    }
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, leftBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Left Border -> receive Right Border
    dest = (rank + size - 1) % size;
    source = (rank + 1) % size;
    for (int i = 0; i < myW; i = i + 1)
    {
        sendbuf[i] = currentfield[calcIndex(myW, 0, i)];
    }
    MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, rightBoder, recvcount, recvtype, source, recvtag, comm, &status);
    printf("DEBUG: Thread No: [%d] sending to [%d] and send [%d][%d][%d][%d]\n", rank, dest, sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3]);
    printf("DEBUG: Thread No: [%d] received from [%d] and send [%d][%d][%d][%d]\n", rank, source, rightBoder[0], rightBoder[1], rightBoder[2], rightBoder[3]);

    //TODO if changes == 0, the time loop will not run!
    return changes;
}

void filling(unsigned* currentfield, int w, int h) {
    for (int i = 0; i < h*w; i++) {
        currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
}

void game(int w, int h, int timesteps, MPI_Status status, MPI_Comm comm, int rank, int size) {

    //Aufgabe B
    MPI_Comm newComm;
    int ndims = 2;
    int *dims = calloc(ndims, sizeof(unsigned));
    dims[0] = 2;
    dims[1] = size / 2;
    int *periods = calloc(ndims, sizeof(int));
    periods[0] = 1;
    periods[1] = 1;
    int reorder = 0;

    MPI_Cart_create(comm, ndims, dims, periods, reorder, &newComm);

    int *myCoords = calloc(ndims, sizeof(int));
    MPI_Cart_coords(newComm, rank, ndims, myCoords);

    int *beginPos = calloc(ndims, sizeof(int));
    int *endPos = calloc(ndims, sizeof(int));
    beginPos[0] = (h / 2) * myCoords[0];
    beginPos[1] = (w / (size / 2)) * myCoords[1];
    endPos[0] = ((h / 2) * (myCoords[0] + 1)) - 1;
    endPos[1] = ((w / (size / 2)) * (myCoords[1] + 1)) - 1;
    int myW = w / (size / 2);
    int myH = h / 2;
    printf("DEBUG: Thread No: [%d] My Cart Coords are [%d] [%d] My starting Coords are [%d] [%d] and My ending Coords are [%d] [%d]\n",rank, myCoords[0], myCoords[1], beginPos[0], beginPos[1], endPos[0], endPos[1]);

    //int beginPos = myCoords[0] * w;
    //int endPos = (myCoords[0] + 1) * w - 1;

    //printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]. Width is [%d] and Height is [%d]\n",rank, size, sourceProcess[0], destProcess[0], w, h);
    //printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]. My first Coordinate is [%d] and my last Coordinate is [%d]\n",rank, size, sourceProcess[0], destProcess[0], beginPos, endPos);


    unsigned *currentfield = calloc(myW * myH, sizeof(unsigned));
    unsigned *newfield = calloc(myW * myH, sizeof(unsigned));


    filling(currentfield, myW, myH);
    int t = 0;
    //for (int t = 0; t < timesteps; t++) {
    // TODO consol output
    //     show(currentfield, w, h);
    writeVTK(currentfield, w, h, t, "output", newComm, rank, size, beginPos, endPos, myCoords, myW, myH);

     int changes = evolve(currentfield, newfield, w, h, status, newComm, rank, size, beginPos, endPos, myCoords, myW, myH);
    /*
     if (changes == 0) {
     sleep(3);
     break;
     }

     // usleep(200000);

     //SWAP
     unsigned *temp = currentfield;
     currentfield = newfield;
     newfield = temp;
     }
     */

    free(currentfield);
    free(newfield);
}

int main(int c, char **v) {
    int rank, size;
    MPI_Status status;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&c, &v);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int w = 0, h = 0, timesteps = 10;
    if (c > 1) w = atoi(v[1]); ///< read width
    if (c > 2) h = atoi(v[2]); ///< read height
    if (c > 3) timesteps = atoi(v[3]);
    if (w <= 0) w = 40; ///< default width
    if (h <= 0) h = 40; ///< default height

    /*
     int length = w * h;
     w = (length / size) * (rank + 1);
     
     
     int posNeighbour = (rank + 1) % size;
     int negNeighbour = (rank - 1 + size) % size;
     printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]\n",rank, size, negNeighbour, posNeighbour);#
     */
    game(w, h, timesteps, status, comm, rank, size);
    
    MPI_Finalize();
    return 0;
}
