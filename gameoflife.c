#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <endian.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

#define calcIndex(width, x,y)  ((y + 1)*(width + 2) + (x + 1))

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

void writeVTK(unsigned* currentfield, int w, int h, int t, char* prefix, MPI_Comm newComm, int rank, int size, int *beginPos, int *endPos, int *myCoords, int myW, int myH) {
    MPI_Status status;
    MPI_File fhw;
    char name[1024] = "\0";

    char header[200] = "# vtk DataFile Version 3.0\n";
    char tempString[100];
    sprintf(tempString, "frame %d\n", t);
    strcat(header, tempString);
    strcat(header, "BINARY\n");
    strcat(header, "DATASET STRUCTURED_POINTS\n");
    sprintf(tempString, "DIMENSIONS %d %d %d \n", h, w, 1); //w h
    strcat(header, tempString);
    strcat(header, "SPACING 1.0 1.0 1.0\n");
    sprintf(tempString, "ORIGIN %d %d 0\n", 0, 0); //x y z
    strcat(header, tempString);
    sprintf(tempString, "POINT_DATA %d\n", w * h);
    strcat(header, tempString);
    strcat(header, "SCALARS data float 1\n");
    strcat(header, "LOOKUP_TABLE default\n");

    int headerSize = strlen(header);


    sprintf(name, "out/%s_%d.vtk", prefix, t);
    MPI_File_open(MPI_COMM_WORLD, name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);

    float *values = calloc(myH*myW, sizeof(float));

    int offset = 0;
    if (rank == 0)
    {
        MPI_File_write_at(fhw, offset, &header, headerSize, MPI_BYTE, &status);
        offset = headerSize;

        MPI_File_write_at(fhw, offset, values, myH * myW * sizeof(float), MPI_BYTE, &status);
    } else
    {
        offset = headerSize + (rank * (myW * myH * sizeof(float)));
    }
    for (int y = 0; y < myH; y++) {
        for (int x = 0; x < myW; x++) {
            float value = currentfield[calcIndex(myW, x,y)]; //!= 0.0 ? 1.0:0.0;
            values[(y * myW) + x] = convert2BigEndian(value);
        }
    }

    MPI_File_write_at(fhw, offset, values, myH * myW * sizeof(float), MPI_BYTE, &status);
    MPI_File_close(&fhw);
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

void calculateOwnBorders(unsigned *currentfield, int myW, int myH, unsigned *sendUpperBorder, unsigned *sendLowerBorder, unsigned *sendRightBorder, unsigned *sendLeftBorder, unsigned *sendUpperRightBorder, unsigned *sendUpperLeftBorder, unsigned *sendLowerRightBorder, unsigned *sendLowerLeftBorder) {
    for (int i = 0; i < myW; i = i + 1)
    {
        sendUpperBorder = &currentfield[calcIndex(myW, i, myH - 1)];
        sendLowerBorder = &currentfield[calcIndex(myW, i, 0)];
    }

    for (int i = 0; i < myH; i = i + 1)
    {
        sendRightBorder[i] = currentfield[calcIndex(myW, myW - 1, i)];
        sendLeftBorder[i] = currentfield[calcIndex(myW, 0, i)];
    }

    sendUpperRightBorder = &currentfield[calcIndex(myW, myW - 1, myH -1)];
    sendUpperLeftBorder = &currentfield[calcIndex(myW, 0, myH -1)];
    sendLowerRightBorder = &currentfield[calcIndex(myW, myW - 1, 0)];
    sendLowerLeftBorder = &currentfield[calcIndex(myW, 0, 0)];
}

void calculateNeighbourBorders(unsigned *currentfield, int myW, int myH, MPI_Comm comm, unsigned *sendUpperBorder, unsigned *sendLowerBorder, unsigned *sendRightBorder, unsigned *sendLeftBorder, unsigned *sendUpperRightBorder, unsigned *sendUpperLeftBorder, unsigned *sendLowerRightBorder, unsigned *sendLowerLeftBorder, unsigned *recvUpperBorder, unsigned *recvLowerBorder, unsigned *recvRightBorder, unsigned *recvLeftBorder, unsigned *recvUpperRightBorder, unsigned *recvUpperLeftBorder, unsigned *recvLowerRightBorder, unsigned *recvLowerLeftBorder, int *upperProcess, int *lowerProcess, int *rightProcess, int *leftProcess, int *upperRightProcess, int *upperLeftProcess, int *lowerRightProcess, int *lowerLeftProcess) {
    MPI_Status status;

    //Send Upper Border -> receive lower Border
    int sendcount = myW;
    MPI_Datatype sendtype = MPI_UNSIGNED;
    int dest = upperProcess[0];
    int sendtag = 1;
    int recvcount = myW;
    MPI_Datatype recvtype = MPI_UNSIGNED;
    int source = lowerProcess[0];
    int recvtag = 1;
    MPI_Sendrecv(sendUpperBorder, sendcount, sendtype, dest, sendtag, recvLowerBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Lower Border -> receive Upper Border
    dest = lowerProcess[0];
    source = upperProcess[0];
    MPI_Sendrecv(sendLowerBorder, sendcount, sendtype, dest, sendtag, recvUpperBorder, recvcount, recvtype, source, recvtag, comm, &status);

    //Send Right Border -> receive Left Border
    sendcount = myH;
    dest = rightProcess[0];
    recvcount = myH;
    source = leftProcess[0];
    MPI_Sendrecv(sendRightBorder, sendcount, sendtype, dest, sendtag, recvLeftBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Left Border -> receive Right Border
    dest = leftProcess[0];
    source = rightProcess[0];
    MPI_Sendrecv(sendLeftBorder, sendcount, sendtype, dest, sendtag, recvRightBorder, recvcount, recvtype, source, recvtag, comm, &status);

    //Send Upper Right Border -> receive Upper Left Border
    sendcount = 1;
    dest = upperRightProcess[0];
    recvcount = 1;
    source = upperLeftProcess[0];
    MPI_Sendrecv(sendUpperRightBorder, sendcount, sendtype, dest, sendtag, recvUpperLeftBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Upper Left Border -> receive Upper Right Border
    dest = upperLeftProcess[0];
    source = upperRightProcess[0];
    MPI_Sendrecv(recvUpperLeftBorder, sendcount, sendtype, dest, sendtag, sendUpperRightBorder, recvcount, recvtype, source, recvtag, comm, &status);

    //Send Lower Right Border -> receive Lower Left Border
    sendcount = 1;
    dest = lowerRightProcess[0];
    recvcount = 1;
    source = lowerLeftProcess[0];
    MPI_Sendrecv(sendLowerRightBorder, sendcount, sendtype, dest, sendtag, sendLowerLeftBorder, recvcount, recvtype, source, recvtag, comm, &status);
    //Send Lower Left Border -> receive Lower Right Border
    dest = lowerLeftProcess[0];
    source = lowerRightProcess[0];
    MPI_Sendrecv(sendLowerLeftBorder, sendcount, sendtype, dest, sendtag, sendLowerRightBorder, recvcount, recvtype, source, recvtag, comm, &status);
}

int calcAlive(unsigned *currentfield, unsigned *newfield, int myW, int myH, unsigned *recvUpperBorder, unsigned *recvLowerBorder, unsigned *recvRightBorder, unsigned *recvLeftBorder, unsigned *recvUpperRightBorder, unsigned *recvUpperLeftBorder, unsigned *recvLowerRightBorder, unsigned *recvLowerLeftBorder)
{
    int changes = 0;
    for (int y = 0; y < myH; y++) {
        for (int x = 0; x < myW; x++) {
            int alive = 0;
            //for schleifen ersetzen durch 8 statische berechnungen?
            for(int y1 = y-1; y1 <= y+1; y1++) {
                for(int x1 = x-1; x1 <= x+1; x1++) {
                    if(currentfield[calcIndex(myW, x1, y1)]) {
                        alive++;
                    }
                }
            }

            //remove inner cell value
            alive -= currentfield[calcIndex(myW, x,y)];
            newfield[calcIndex(myW, x,y)] = (alive == 3 || (alive == 2 && currentfield[calcIndex(myW, x,y)]));
            if(currentfield[calcIndex(myW, x,y)] != newfield[calcIndex(myW, x,y)])
            {
                changes = 1;
            }
        }
    }
    return changes;
}

void mergeBorders(unsigned *currentfield, int myW, int myH, unsigned *recvUpperBorder, unsigned *recvLowerBorder, unsigned *recvRightBorder, unsigned *recvLeftBorder, unsigned *recvUpperRightBorder, unsigned *recvUpperLeftBorder, unsigned *recvLowerRightBorder, unsigned *recvLowerLeftBorder) {

    for (int i = 0; i < myW; i = i + 1)
    {
        currentfield[calcIndex(myW, i, myH)] = recvUpperBorder[i];
        currentfield[calcIndex(myW, i, -1)] = recvLowerBorder[i];
    }

    for (int i = 0; i < myH; i = i + 1)
    {
        currentfield[calcIndex(myW, myW, i)] = recvRightBorder[i];
        currentfield[calcIndex(myW, -1, i)] = recvLeftBorder[i];
    }

    currentfield[calcIndex(myW, myW, myH)] = recvUpperRightBorder[0];
    currentfield[calcIndex(myW, -1, myH)] = recvUpperLeftBorder[0];
    currentfield[calcIndex(myW, myW, -1)] = recvLowerRightBorder[0];
    currentfield[calcIndex(myW, -1, -1)] = recvLowerLeftBorder[0];
}


int evolve(unsigned *currentfield, unsigned *newfield, int w, int h, MPI_Status status, MPI_Comm comm, int rank, int size, int *beginPos, int *endPos, int *myCoords, int myW, int myH) {
    int changes = 0;

    unsigned *sendUpperBorder = calloc(myW, sizeof(unsigned));
    unsigned *sendLowerBorder = calloc(myW, sizeof(unsigned));
    unsigned *sendRightBorder = calloc(myH, sizeof(unsigned));
    unsigned *sendLeftBorder = calloc(myH, sizeof(unsigned));
    unsigned *sendUpperRightBorder = calloc(1, sizeof(unsigned));
    unsigned *sendUpperLeftBorder = calloc(1, sizeof(unsigned));
    unsigned *sendLowerRightBorder = calloc(1, sizeof(unsigned));
    unsigned *sendLowerLeftBorder = calloc(1, sizeof(unsigned));

    unsigned *recvLowerBorder = calloc(myW, sizeof(unsigned));
    unsigned *recvUpperBorder = calloc(myW, sizeof(unsigned));
    unsigned *recvRightBorder = calloc(myH, sizeof(unsigned));
    unsigned *recvLeftBorder = calloc(myH, sizeof(unsigned));
    unsigned *recvUpperRightBorder = calloc(1, sizeof(unsigned));
    unsigned *recvUpperLeftBorder = calloc(1, sizeof(unsigned));
    unsigned *recvLowerRightBorder = calloc(1, sizeof(unsigned));
    unsigned *recvLowerLeftBorder = calloc(1, sizeof(unsigned));

    int *upperProcess = calloc (1, sizeof(int));
    int *lowerProcess = calloc (1, sizeof(int));
    int *rightProcess = calloc (1, sizeof(int));
    int *leftProcess = calloc (1, sizeof(int));
    int *upperRightProcess = calloc (1, sizeof(int));
    int *upperLeftProcess = calloc (1, sizeof(int));
    int *lowerRightProcess = calloc (1, sizeof(int));
    int *lowerLeftProcess = calloc (1, sizeof(int));

    calculateNeighbourProcesses(myCoords, comm, upperProcess, lowerProcess, rightProcess, leftProcess, upperRightProcess, upperLeftProcess, lowerRightProcess, lowerLeftProcess);

    calculateOwnBorders(currentfield, myW, myH, sendUpperBorder, sendLowerBorder, sendRightBorder, sendLeftBorder, sendUpperRightBorder, sendUpperLeftBorder, sendLowerRightBorder, sendLowerLeftBorder);

    calculateNeighbourBorders(currentfield, myW, myH, comm, sendUpperBorder, sendLowerBorder, sendRightBorder, sendLeftBorder, sendUpperRightBorder, sendUpperLeftBorder, sendLowerRightBorder, sendLowerLeftBorder, recvUpperBorder, recvLowerBorder, recvRightBorder, recvLeftBorder, recvUpperRightBorder, recvUpperLeftBorder, recvLowerRightBorder, recvLowerLeftBorder, upperProcess, lowerProcess, rightProcess, leftProcess, upperRightProcess, upperLeftProcess, lowerRightProcess, lowerLeftProcess);

    mergeBorders(currentfield, myW, myH, recvUpperBorder, recvLowerBorder, recvRightBorder, recvLeftBorder, recvUpperRightBorder, recvUpperLeftBorder, recvLowerRightBorder, recvLowerLeftBorder);

    changes = calcAlive(currentfield, newfield, myW, myH, recvUpperBorder, recvLowerBorder, recvRightBorder, recvLeftBorder, recvUpperRightBorder, recvUpperLeftBorder, recvLowerRightBorder, recvLowerLeftBorder);

    free(sendUpperBorder);
    free(sendLowerBorder);
    free(sendRightBorder);
    free(sendLeftBorder);
    free(sendUpperRightBorder);
    free(sendUpperLeftBorder);
    free(sendLowerRightBorder);
    free(sendLowerLeftBorder);

    free(recvLowerBorder);
    free(recvUpperBorder);
    free(recvRightBorder);
    free(recvLeftBorder);
    free(recvUpperRightBorder);
    free(recvUpperLeftBorder);
    free(recvLowerRightBorder);
    free(recvLowerLeftBorder);

    free(upperProcess);
    free(lowerProcess);
    free(rightProcess);
    free(leftProcess);
    free(upperRightProcess);
    free(upperLeftProcess);
    free(lowerRightProcess);
    free(lowerLeftProcess);

    //TODO if changes == 0, the time loop will not run!
    return changes;
}

void filling(unsigned* currentfield, int w, int h, int rank) {
    clock_t t = time(NULL);
    srand((int) t * rank);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            currentfield[calcIndex(w, x, y)] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
        }
    }
}

void game(int w, int h, int timesteps, MPI_Status status, MPI_Comm comm, int rank, int size) {

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
    //printf("DEBUG: Thread No: [%d] My Cart Coords are [%d] [%d] My starting Coords are [%d] [%d] and My ending Coords are [%d] [%d]\n",rank, myCoords[0], myCoords[1], beginPos[0], beginPos[1], endPos[0], endPos[1]);

    unsigned *currentfield = calloc((myW + 2) * (myH + 2), sizeof(unsigned));
    unsigned *newfield = calloc((myW + 2) * (myH + 2), sizeof(unsigned));

    filling(currentfield, myW, myH, rank);

    for (int t = 0; t < timesteps; t++) {
        writeVTK(currentfield, w, h, t, "output", newComm, rank, size, beginPos, endPos, myCoords, myW, myH);

        int changes = evolve(currentfield, newfield, w, h, status, newComm, rank, size, beginPos, endPos, myCoords, myW, myH);
        int allchanges = changes;
        MPI_Allreduce(&changes, &allchanges, 1, MPI_INT, MPI_SUM, comm);
        if (allchanges == 0) {
            sleep(3);
            break;
        }
        //SWAP
        unsigned *temp = currentfield;
        currentfield = newfield;
        newfield = temp;
    }

    free(dims);
    free(periods);
    free(myCoords);
    free(beginPos);
    free(endPos);
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

    int w = 0, h = 0, timesteps = 1;
    if (c > 1) w = atoi(v[1]); ///< read width
    if (c > 2) h = atoi(v[2]); ///< read height
    if (c > 3) timesteps = atoi(v[3]);
    if (w <= 0) w = 50; ///< default width
    if (h <= 0) h = 50; ///< default height
    
    //    printf("DEBUG: Thread No: [%d] Thread Size: [%d] My Neighbours are: [%d] and: [%d]\n",rank, size, negNeighbour, posNeighbour);#
    game(w, h, timesteps, status, comm, rank, size);
    
    MPI_Finalize();
    return 0;
}
